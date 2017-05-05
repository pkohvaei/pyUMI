#!/usr/bin/env python

import logging
import numpy as np
from collections import Counter
from operator import itemgetter
from matplotlib import pyplot as plt

import cPickle as pickle
import pysam

from MSolverLib import build_multimapping_hashtable, build_uniques_associations, find_barcode_matches\
    , gene_tagged_multi_maps
from util import prepare_file_name, initialize_iterator, plot_bar


def format_columns(text,width): return width - len(text)


logger = logging.getLogger(__name__)

samples = ['GACCGC', 'AAAACT', 'GGCGTC', 'AAAGTT', 'GTTCGA', 'ATATAG',
           'TAAAGT', 'ATCAAA', 'TCTGCA' , 'CCCTGG', 'TTAATC', 'CCGGAC']
organism = 'Mouse'


def multimap_multigene_generator(multimaps, variance):

    for multimap in multimaps:
        genes = np.unique(np.array([x[0] for x in multimaps[multimap][2] if x[0]]))
        if len(genes) > variance:
            xc = multimaps[multimap][0]
            xm = multimaps[multimap][1]
            yield multimap, xc, xm, genes


def find_gene_hits(uniques, matches, genes):

    hits = []
    resolvers = []
    for match in matches:
        if uniques[match][0] in genes:
            hits.append((uniques[match]))
            resolvers.append(match)
    return hits, resolvers


def rank1_multimap_solver(mm_generator=None, uniques_hashtable=None, uniques=None, report=False):

    rank1_solvers = []
    resolved_mm = []
    total = 0
    r_total = 0

    for mm in mm_generator:
        total += 1
        qname, xc, xm, genes = mm
        uniq_matches = find_barcode_matches(uniques_hashtable, xc, xm)
        if len(uniq_matches) > 0:
            hits, resolvers = find_gene_hits(uniques, uniq_matches, genes)
            if hits:
                unique_hit = np.unique(np.array([hit[0] for hit in hits]))
                if len(unique_hit) == 1:
                    rank1_solvers.append((resolvers[0], unique_hit[0]))
                    resolved_mm.append((qname, unique_hit[0]))
                    r_total += len(genes)

    if report:
        b = "Number of multi reads with unique gene association : %s"
        c = "Total number of mappings of multi reads with unique gene association : %s"

        logger.info("\t" + b %format(len(resolved_mm),",").rjust(format_columns(b, 70)))
        logger.info("\t" + c % format(r_total, ",").rjust(format_columns(c, 70)))

    return resolved_mm, rank1_solvers, r_total


def majority_fraction_generator(mm_generator=None, uniques_hashtable=None, uniques=None):

    for mm in mm_generator:
        mj=0
        qname, xc, xm, genes = mm
        uniq_matches = find_barcode_matches(uniques_hashtable, xc, xm)
        if len(uniq_matches) > 0:
            hits = find_gene_hits(uniques, uniq_matches, genes)
            if hits:
                h = [x for x, _, _ in hits]
                votes = sorted(Counter(h).items())
                if len(votes) == 1:
                    mj = 1
                else:
                    mj = majority_fraction(votes)

        yield qname,mj


def majority_fraction(votes):

    population = [frequency for (gene,frequency) in votes]
    return float(population[-1])/sum(population)


def majority_fraction_report(mf_generator=None, draw_bar=False, color='r', font_size=4):

    f_dist = [fraction for _, fraction in mf_generator]
    f_dist_rounded = [float(format(f, '.2f')) for f in f_dist]
    dist = Counter(f_dist_rounded)

    x_largest = max(dist.iteritems(), key=itemgetter(0))[0]
    y_largest = max(dist.iteritems(), key=itemgetter(1))[1]

    if draw_bar:
        fig = plt.figure(figsize=(15, 7.5))
        ax = fig.add_subplot(1, 1, 1)

        x_major_ticks = np.arange(0, x_largest, .1)
        x_minor_ticks = np.arange(0, x_largest, .05)

        y_major_ticks = np.arange(0, y_largest + 100, 200000)
        y_minor_ticks = np.arange(0, y_largest + 100, 50000)

        ax.set_xticks(x_major_ticks)
        ax.set_xticks(x_minor_ticks, minor=True)
        ax.set_yticks(y_major_ticks)
        ax.set_yticks(y_minor_ticks, minor=True)

        fig.suptitle('Majority voting gene probability distribution for multimapped reads')
        plt.xlabel('Dominant vote probability')
        plt.ylabel('Read frequency')
        plt.bar(dist.keys(), dist.values(), .02, color=color, align='center');

    return dist


def unique_savers_dist(solved_multimaps, multimaps_hashtable, uniques_hashtable, draw_bar=False, **kwargs):

    x_scale = kwargs.pop('xscale',10)
    bar_color = kwargs.pop('color','b')
    bar_title = kwargs.pop('title','Distribution of multi_mapped alignments supported by unique gene pointer')

    dist_list = []

    for s in solved_multimaps:
        qname, _ = s
        xc,xm,_ = multimaps_hashtable[qname]
        dist_list.append(len(find_barcode_matches(uniques_hashtable, xc, xm)))

    dist = Counter(dist_list)

    largest = max(dist.iteritems(), key=itemgetter(1))[1]

    if draw_bar:
        sorted_dist = sorted(dist.items())
        x = [i for (i, j) in sorted_dist]
        y = [j for (i, j) in sorted_dist]

        plot_bar(x, y, 'Distribution of number of single-gene associations in unique reads / multi read', 'Number of supporting votes', 'Read frequency', (10, 4))

        '''
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(1, 1, 1)

        x_major_ticks = np.arange(0, len(dist), 10)
        x_minor_ticks = np.arange(0, len(dist), 2)

        y_major_ticks = np.arange(0, largest + 100, 1000)
        y_minor_ticks = np.arange(0, largest + 100, 500)

        ax.set_xticks(x_major_ticks)
        ax.set_xticks(x_minor_ticks, minor=True)
        ax.set_yticks(y_major_ticks)
        ax.set_yticks(y_minor_ticks, minor=True)

        fig.suptitle(bar_title)
        plt.xlabel('Number of supporting votes')
        plt.ylabel('Read frequency')
        plt.bar(range(len(dist)), dist.values(), .2, color='b', align='center');
        '''
    return dist


'''
def resolve_rank1_multimaps(alignmentFile, cbs=None, variance=0):

    mm_hashtable = build_multimapping_hashtable(alignmentFile)
    uniq_hashtable, uniq_map = build_uniques_associations(alignmentFile, cell_barcodes=cbs)
    r = multimap_multigene_generator(mm_hashtable, variance)
    resolved, savers = rank1_multimap_solver(r, uniq_hashtable, uniq_map, report=True)
    unique_savers_dist(savers, mm_hashtable, uniq_hashtable, draw_bar=True)
'''

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to rank1 multimap solver package.')

    path = '/data/UMI/data/MUS/R1/'
    pt = '/data/UMI/data/MUS/'

    file_names = []
    samples = ['GACCGC']
    for sample in samples:

        f1 = prepare_file_name(sample, pt, 'sample_', '.bam')
        f2 = prepare_file_name(sample, path, 'u_tab_', '.pkl')
        f3 = prepare_file_name(sample, path, 'u_map_', '.pkl')
        f4 = prepare_file_name(sample, path, 'm_tab_', '.pkl')

        file_names.append((sample, f1, f2, f3, f4))

    for sample, f1, f2, f3, f4 in file_names:

        st = pysam.AlignmentFile(f1,"rb")

        u_tab = pickle.load(open(f2, 'r'))
        u_map = pickle.load(open(f3, 'r'))
        m_tab = pickle.load(open(f4, 'r'))

        '''
        m_tab = build_multimapping_hashtable(st, sample)
        u_tab, u_map = build_uniques_associations(st, [sample])
                pickle.dump(u_tab, open(f2, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
        pickle.dump(u_map, open(f3, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
        pickle.dump(m_tab, open(f4, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

        logger.info('Generated views for sample %s.' % sample)
        '''

        multi_reads, alignments, gene_tagged, all_gene_tagged, only_ref = gene_tagged_multi_maps(m_tab)

        variance = 1

        mm_gen = multimap_multigene_generator(m_tab, variance)
        resolved, savers = rank1_multimap_solver(mm_gen, u_tab, u_map, report=True)
        #dist = unique_savers_dist(resolved, m_tab, u_tab, draw_bar=True)

        #number of multi reads
        #number of mappings
        #number of unique reads
        #logger.info("\t" + "Generator for multimapped reads with at least %d distinct gene annotation(s)." %(variance+1))


