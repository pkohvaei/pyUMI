#!/usr/bin/env python

import logging
import time
from collections import Counter
import numpy as np
import pickle

import pysam

from util import initialize_iterator, prepare_file_name


logger = logging.getLogger(__name__)

samples = ['GACCGC', 'AAAACT', 'GGCGTC', 'AAAGTT', 'GTTCGA', 'ATATAG',
           'TAAAGT', 'ATCAAA', 'TCTGCA' , 'CCCTGG', 'TTAATC', 'CCGGAC']
organism = 'Mouse'
path = '/data/UMI/data/MUS/RP/'
pt = '/data/UMI/data/MUS/'


def repeat_families(repeats):

    families = []

    for gene in repeats:
        prefix = gene[:gene.rfind('dup') -1]
        families.append(prefix)

    return Counter(families).keys(), Counter(families).values()


def categorize_genes(genes):

    gene_list = [item for item in genes if item is not None]

    if not gene_list:
        return None

    g = []
    r = []
    for gene in gene_list:
        if 'ENSMUS' in gene:
            g.append(gene)
        else:
            r.append(gene)
    return r, g


def build_r2g_table(pysam_iterator):

    reads = initialize_iterator(pysam_iterator)

    r2g = {}
    for r in reads:
        if r.has_tag('GE'):
            ge = r.get_tag('GE')
            qname = r.query_name
            if qname in r2g:
                r2g[qname].append(ge)
            else:
                r2g.update({qname: [ge]})

    return r2g


def build_umi_based_r2g_table(st):

    st.reset()
    reads = st.fetch(until_eof=True)

    r2g = {}
    for r in reads:

        if not r.is_unmapped and r.get_tag('NH') > 1:
            xm = r.get_tag('XM')
            ge = None
            if r.has_tag('GE'):
                ge = r.get_tag('GE')
            qname = r.query_name

            if not xm in r2g:
                r2g.update({xm:{qname:[ge]}})
            else:
                if qname in r2g[xm]:
                    r2g[xm][qname].append(ge)
                else:
                    r2g[xm].update({qname: [ge]})

    return r2g


def get_non_isolated_umis(pysam_iterator):

    reads = initialize_iterator(pysam_iterator)

    umis = {}
    u1 = []
    for r in reads:
        if not r.is_unmapped:
            xm = r.get_tag('XM')
            nh = r.get_tag('NH')
            if nh > 1:
                qn = r.query_name
                if xm in umis:
                    umis[xm].append(qn)
                else:
                    umis.update({xm: [qn]})
            elif nh == 1:
                u1.append(xm)

    uniqs_umis = Counter(u1).keys()

    non_isolated = []
    pending = []
    for umi in umis:
        p = len(np.unique(umis[umi]))
        if p > 1:
            non_isolated.append(umi)
        else:
            pending.append(umi)

    intersect = list(set(pending).intersection(uniqs_umis))

    non_isolated += intersect
    isolated = list(set(pending) - set(intersect))

    return non_isolated, isolated


'''
def gene_filtered_read_generator(alignmentIterable, exist=True, pattern=''):
    """Generator of type pysam.alignedSegment.
       Takes an iterable of type pysam.alignedSegment as input.
       Filters the input file as follows based on 'exist' and 'pattern'.
       Yields filtered reads."""

    for alignedSegment in alignmentIterable:
        if alignedSegment.has_tag("GE"):
            tag = alignedSegment.get_tag("GE")
            if pattern in tag:
                if exist:
                    yield alignedSegment
            else:
                if not(exist):
                    yield alignedSegment
'''


'''
def repeat_stats(alignmentFile, domain='reads', report=True, draw_pie=False):

    domain_str = ''
    if domain == 'read':
        domain_str = ' bam file '
        gen = gene_filtered_read_generator(alignmentFile, exist=False, pattern='ENSMUS')
        alignmentFile.reset()
        total = len([i for i in ngene_annotated_generator(alignmentFile, annotate=True)])
    elif domain == 'unique':
        domain_str = ' uniquely mapped '
        gen = gene_filtered_read_generator(stats.uniques_generator(alignmentFile), exist=False, pattern='ENSMUS')
        alignmentFile.reset()
        _,total = stats.count_um(alignmentFile)
    elif domain == 'multimap':
        domain_str = ' multimapped '
        gen = gene_filtered_read_generator(stats.multimap_generator(alignmentFile), exist=False, pattern='ENSMUS')
        alignmentFile.reset()
        total,_ = stats.count_um(alignmentFile)
    else:
        logger.info("Domain should be one of: 'read' , 'unique' , 'multimap'.")
        return

    alignmentFile.reset()
    repeats_cnt = len([i for i in gen])

    if report:
        logger.info("\t" + "Total number of" + domain_str +  "reads:" + "%s" % format(total, ",").rjust(25))
        logger.info("\t" + "Number of repeats in" + domain_str + "reads:" + "%s" % format(repeats_cnt, ",").rjust(20))

    if draw_pie:
        labels = domain_str , domain_str + "repeats"
        plt.pie([total-repeats_cnt,repeats_cnt], labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)

    return total, repeats_cnt
'''

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    fh = logging.FileHandler('/data/UMI/data/MUS/RP/RepeatDistReport.log')
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)

    logger.info('Call to Repeat Solver Library.')

    logger.info('Repeat distribution report for %d samples.' %len(samples))
    logger.info('Organism : %s\n' %organism)

    file_names = []
    for sample in samples:

        f1 = prepare_file_name(sample, pt, 'sample_', '.bam')
        f2 = prepare_file_name(sample, path, 'r2r_', '.pkl')
        file_names.append((f1, f2))

    counter = 1

    for in_file, out_file in file_names:
        logger.info('\n')
        logger.info('-' * 100)
        logger.info('Sample no. %d : %s' % (counter, in_file))
        logger.info('-' * 100)
        logger.info('\n')

        start_time = time.time()

        st = pysam.AlignmentFile(in_file, "rb")

        non_iso, iso = get_non_isolated_umis(st)

        r2g_table = build_umi_based_r2g_table(st)

        mm = sum([len(item.keys()) for item in r2g_table.values()])

        r2r_table = {}
        total = 0

        for umi in non_iso:
            if umi in r2g_table:
                for qname in r2g_table[umi]:
                    tup = categorize_genes(r2g_table[umi][qname])
                    if tup:
                        r, _ = tup
                        r2r_table.update({qname: r})
                    total += 1

        pickle.dump(r2r_table, open(out_file, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

        repeat_tagged = []
        family_tagged = []

        for qname in r2r_table:
            repeats = r2r_table[qname]
            if repeats:
                repeat_tagged.append(qname)
                _, v = repeat_families(repeats)
                if min(v) > 1:
                    family_tagged.append(qname)

        con1 = len(repeat_tagged)
        per_con1 = str(round(float(con1)/total, 3))
        con2 = len(family_tagged)
        per_con2 = str(round(float(con2)/total, 3))

        end_time = time.time()

        logger.info('Total number of multi reads: %s\n' % (format(mm,',')))
        logger.info('-'*80)
        logger.info('Number of multi reads with non-isolated UMIs: %s\n' %(format(total,',')))
        logger.info('Non-isolated multi reads with at least 1 repeat copy: %s (%s%%)\n' % (format(con1,','), per_con1))
        logger.info('Non-isolated multi reads with at least 2 copies of the same repeat family: %s (%s%%)\n' % (format(con2,','), per_con2))

        counter += 1












