#!/usr/bin/env python

import logging

import pysam

from umiViews import tag_based_generator
from umiDistAnalyzer import isolated_umis, get_umi_list
from util import prepare_file_name, initialize_iterator
from MSolverLib import get_reads

samples = ['GACCGC', 'AAAACT', 'GGCGTC', 'AAAGTT', 'GTTCGA', 'ATATAG',
           'TAAAGT', 'ATCAAA', 'TCTGCA' , 'CCCTGG', 'TTAATC', 'CCGGAC']
organism = 'Mouse'

path = '/data/UMI/data/MUS/'

logger = logging.getLogger(__name__)


def scoring_umi_table(pysam_iter, values=[2]):

    pysam_iter.reset()
    g = tag_based_generator(pysam_iter, tag='NH', values=values)

    htab = {}
    for r in g:
        score = r.get_tag('AS')
        nm = r.get_tag('NM')
        xm = r.get_tag('XM')
        prop = (score, nm)
        if xm in htab:
            htab[xm].append(prop)
        else:
            htab.update({xm: [prop]})

    return htab


def scoring_hash_table(pysam_iter):

    pysam_iter.reset()
    g = tag_based_generator(pysam_iter, tag='NH', values=[2])

    htab = {}
    for r in g:
        qn = r.query_name
        score = r.get_tag('AS')
        nm = r.get_tag('NM')
        prop = (score, nm)
        if qn in htab:
            htab[qn].append(prop)
        else:
            htab.update({qn: [prop]})

    return htab


def iso_scores(scores_dict, iso_list):

    x = {your_key: scores_dict[your_key] for your_key in iso_list}

    return x


def nequal(scores_dict):

    non_equals = [item for item in scores_dict.values() if len(item) > 1 if not item[0] == item[1]]

    return non_equals


def best_alignment(as_nm_list):
    max_score = max([i for i,j in as_nm_list])
    l1 = [(i,j) for i,j in as_nm_list if i == max_score]
    best = []
    if len(l1) ==1:
        best = l1
    else:
        min_mismatch = min([j for i,j in l1])
        l2 = [(i,j) for i,j in l1 if j == min_mismatch]
        best = l2
    return best


def iso_report(pysam_iter):

    g = tag_based_generator(pysam_iter, tag='NH', values=[2])
    m_umis = get_umi_list(g)
    #iso_all = isolated_umis(pysam_iter)

    #iso = list(set(m_umis).intersection(iso_all))

    sht = scoring_hash_table(pysam_iter)
    #scores = iso_scores(sht, iso)
    nequals = nequal(sht)
    #nequals_iso = nequal(scores)

    a = len(sht)
    #b = len(nequals_iso)
    c = len(nequals)

    pysam_iter.reset()
    g = tag_based_generator(pysam_iter, tag='NH', values=[1], include=False)
    d = len(get_reads(g))
    #prct1 = str(round(float(b)/d * 100, 4)) + '%'
    prct2 = str(round(float(c)/d * 100, 1)) + '%'


    logger.info('Number of multi reads with two mappings: %s\n' % format(a,","))
    logger.info('Number of multi reads with two mappings and a possible best mapping: %s\n' % format(c, ","))
    logger.info('Fraction of multi reads: %s\n' % prct2)


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    fh = logging.FileHandler('/data/UMI/data/MUS/DT/isoDmappedReport.log')
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)

    logger.info('Call to double_mapped solver module. \n')

    logger.info('Double_mapped multi reads report for %d samples.' %len(samples))
    logger.info('Organism : %s\n' %organism)

    samples = ['CCGGAC']
    sample_counter = 1

    for sample in samples:
        logger.info('\n')
        logger.info('-' * 100)
        logger.info('Sample no. %d, cellular barcode : %s' % (sample_counter, sample))
        logger.info('-' * 100)
        logger.info('\n')

        in_file = prepare_file_name(sample, path, 'sample_', '.bam')
        pysam_iter = pysam.AlignmentFile(in_file, "rb")
        iso_report(pysam_iter)

        sample_counter += 1


