#!/usr/bin/env python

import logging
import numpy as np
from collections import Counter
from operator import itemgetter
from matplotlib import pyplot as plt

import cPickle as pickle
import pysam

from MSolverLib import gene_tagged_multi_maps
from util import prepare_file_name

import R1Solver as r1s
import R2Solver as r2s


def format_columns(text,width): return width - len(text)


logger = logging.getLogger(__name__)

samples = ['GACCGC', 'AAAACT', 'GGCGTC', 'AAAGTT', 'GTTCGA', 'ATATAG',
           'TAAAGT', 'ATCAAA', 'TCTGCA' , 'CCCTGG', 'TTAATC', 'CCGGAC']
organism = 'Mouse'


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    fh = logging.FileHandler('/data/UMI/data/MUS/multiSolverReport.log')
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)

    logger.info('Call to Multimap Solver Report module.\n')

    logger.info('Mapping solver report for %d samples.' %len(samples))
    logger.info('Organism : %s' % organism)

    path = '/data/UMI/data/MUS/R1/'
    pt = '/data/UMI/data/MUS/'

    file_names = []
    sample_counter = 1

    for sample in samples:

        f1 = prepare_file_name(sample, pt, 'sample_', '.bam')
        f2 = prepare_file_name(sample, path, 'u_tab_', '.pkl')
        f3 = prepare_file_name(sample, path, 'u_map_', '.pkl')
        f4 = prepare_file_name(sample, path, 'm_tab_', '.pkl')

        file_names.append((sample, f1, f2, f3, f4))

    for sample, f1, f2, f3, f4 in file_names:

        logger.info('\n')
        logger.info('-' * 100)
        logger.info('Sample no. %d - cellular barcode : %s\n' % (sample_counter, sample))
        logger.info('-' * 100)
        logger.info('\n')

        st = pysam.AlignmentFile(f1,"rb")

        u_tab = pickle.load(open(f2, 'r'))
        u_map = pickle.load(open(f3, 'r'))
        m_tab = pickle.load(open(f4, 'r'))

        multi_reads, alignments, gene_tagged, all_gene_tagged, only_ref = gene_tagged_multi_maps(m_tab)

        logger.info('Multi read / mapping statistics :')
        logger.info('-' * 40)
        logger.info('\n')

        txt1 = "Number of multi reads: %s\n"
        txt2 = "Number of associated mappings: %s\n"
        txt3 = "Total number of mappings of reads with at least one gene association: %s\n"
        txt4 = "Total number of all-gene-associated mappings among reads: %s\n"
        txt5 = "Total number of mappings among reads with no gene association: %s\n\n\n"

        l = 90
        logger.info(txt1 % format(multi_reads,",").rjust(format_columns(txt1, l)))
        logger.info(txt2 % format(alignments,",").rjust(format_columns(txt2, l)))
        logger.info(txt3 % format(gene_tagged,",").rjust(format_columns(txt3, l)))
        logger.info(txt4 % format(all_gene_tagged,",").rjust(format_columns(txt4, l)))
        logger.info(txt5 % format(only_ref,",").rjust(format_columns(txt5, l) + 2))

        variance = 1
        window_size = 1500

        mm_gen = r1s.multimap_multigene_generator(m_tab, variance)
        resolved1, savers1, r_t1 = r1s.rank1_multimap_solver(mm_gen, u_tab, u_map)
        f1 = str(round(float(len(resolved1))/multi_reads, 3)*100) + '%'
        f2 = str(round(float(r_t1)/alignments, 3)*100) + '%'

        txt1 = "Number of multi reads with unique gene association: %s\n"
        tf1 = "Fraction of all multi reads: %s\n"
        txt2 = "Total number of mappings of multi reads with unique gene association: %s\n"
        tf2 = "Fraction of all mappings of multi reads: %s\n\n\n"

        logger.info('Mapping solver1 report: multi mappings with a unique strong gene association :')
        logger.info('-' * 40)
        logger.info('\n')

        logger.info(txt1 % format(len(resolved1),",").rjust(format_columns(txt1, l)))
        logger.info(tf1 % f1.rjust(format_columns(tf1, l)))
        logger.info(txt2 % format(r_t1, ",").rjust(format_columns(txt2, l)))
        logger.info(tf2 % f2.rjust(format_columns(tf2, l) + 2))

        mm_reg = r2s.multimap_multiregion_generator(m_tab, variance)
        resolved2, savers2, r_t2 = r2s.rank2_multimap_solver(mm_reg, u_tab, u_map, window_size)
        f3 = str(round(float(len(resolved2))/multi_reads, 3)*100) + '%'
        f4 = str(round(float(r_t2)/alignments, 3)*100) + '%'

        logger.info('Mapping solver2 report: multi mappings with a unique strong reference association :\n')
        logger.info("( Reference assignment vicinity threshold : " + str(window_size) + " bps. )")
        logger.info('-' * 40)
        logger.info('\n')

        txt1 = "Number of acceptable multi read assignments to a unique genomic reference: %s\n"
        txt2 = "Total number of alignments of assigned multi reads: %s\n"

        logger.info(txt1 %format(len(resolved2),",").rjust(format_columns(txt1, l)))
        logger.info(tf1 % f3.rjust(format_columns(tf1, l)))
        logger.info(txt2 % format(r_t2, ",").rjust(format_columns(txt2, l)))
        logger.info(tf2 % f4.rjust(format_columns(tf2, l) + 2))

        total_r = r_t1 + r_t2
        total_m = len(resolved1) + len(resolved2)
        tm_f = str(round(float(total_m)/multi_reads, 3)*100) + '%'
        tr_f = str(round(float(total_r)/alignments, 3)*100) + '%'

        logger.info('Summary :')
        logger.info('-' * 40)
        logger.info('\n')

        txt1 = "Total number of multi reads with strong associations: %s\n"
        txt2 = "Total number of mappings of multi reads with strong associations: %s\n"

        logger.info(txt1 % format(total_m,",").rjust(format_columns(txt1, l)))
        logger.info(tf1 % tm_f.rjust(format_columns(tf1, l)))
        logger.info(txt2 % format(total_r, ",").rjust(format_columns(txt2, l)))
        logger.info(tf2 % tr_f.rjust(format_columns(tf2, l) + 2))

        logger.info('\n\n\n')
        # dist = unique_savers_dist(resolved, m_tab, u_tab, draw_bar=True)
        # logger.info("Generator for multimapped reads with at least %d distinct gene annotation(s)." %(variance+1))

        sample_counter += 1
