#!/usr/bin/env python

import logging
import pickle
from numpy import mean

from util import prepare_file_name

samples = ['GACCGC', 'AAAACT', 'GGCGTC', 'AAAGTT', 'GTTCGA', 'ATATAG',
           'TAAAGT', 'ATCAAA', 'TCTGCA' , 'CCCTGG', 'TTAATC', 'CCGGAC']
organism = 'Mouse'

path = '/data/UMI/data/MUS/DT/'

logger = logging.getLogger(__name__)


def distribution_report(dist_table):

    v = dist_table.values()
    I = [i for i, j, k, w in v]
    J = [j for i, j, k, w in v]
    K = [k for i, j, k, w in v]
    W = [w for i, j, k, w in v]

    total = sum(I)
    logger.info('Total number of reads: %s\n' % format(total,","))
    logger.info('Number of distinct UMIs: %s\n' % format(len(I), ","))
    logger.info('Smallest mapping population among UMIs: %d' % min(I))
    logger.info('Largest mapping population among UMIs: %d\n' % max(I))
    logger.info('Average number of references per UMI: %d' % mean(J))
    logger.info('Minimum number of references per UMI: %d' % min(J))
    logger.info('Maximum number of references per UMI: %d\n' % max(J))
    logger.info('Minimum mapping frequency among all references and all UMIs: %d' % min(K))
    logger.info('Maximum mapping frequency among all references and all UMIs: %d\n' % max(W))


def multi_read_stats_report(cnt_table):

    v = cnt_table.values()
    total = sum([i for j, i in v])
    logger.info('Total number of mappings: %s' % format(total,","))
    s = min([j for j, i in v])
    b = max([j for j, i in v])

    logger.info('Minimum frequency of mappings of a read: %d' % s)
    logger.info('Maximum frequency of mappings of a read: %d\n' % b)


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    fh = logging.FileHandler('/data/UMI/data/MUS/DT/umiDistReport.log')
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)

    logger.info('Call to UMI distribution report module. \n')

    logger.info('UMI distribution report for %d samples.' %len(samples))
    logger.info('Organism : %s\n' %organism)

    sample_counter = 1

    for sample in samples:
        logger.info('\n')
        logger.info('-' * 100)
        logger.info('Sample no. %d, cellular barcode : %s' % (sample_counter, sample))
        logger.info('-' * 100)
        logger.info('\n')

        f1 = prepare_file_name(sample, path, 'summ_ut_', '.pkl')
        f2 = prepare_file_name(sample, path, 'summ_mt_', '.pkl')
        f3 = prepare_file_name(sample, path, 'summ_cnt_', '.pkl')

        uniq_summary = pickle.load(open(f1, 'r'))
        multi_summary = pickle.load(open(f2, 'r'))
        multi_stat_summary = pickle.load(open(f3, 'r'))

        logger.info('UMI distribution among unique reads:\n')
        logger.info('-' * 40)
        distribution_report(uniq_summary)
        logger.info('\n')
        logger.info('UMI distribution among multi reads:\n')
        logger.info('-' * 40)
        multi_read_stats_report(multi_stat_summary)
        distribution_report(multi_summary)

        sample_counter += 1


