#!/usr/bin/env python

import logging
from collections import Counter
import pickle
import time
import multiprocessing as mp

import pysam


from util import prepare_file_name, initialize_iterator, pickle_it


samples = ['GACCGC', 'AAAACT', 'GGCGTC', 'AAAGTT', 'GTTCGA', 'ATATAG',
           'TAAAGT', 'ATCAAA', 'TCTGCA' , 'CCCTGG', 'TTAATC', 'CCGGAC']

path = '/data/UMI/data/MUS/DT/'
pt = '/data/UMI/data/MUS/'

umi_length = 10


logger = logging.getLogger(__name__)


''' Functions for extracting umi information form the bam file '''


def tag_reads(pysam_iter):
    all_r = []
    for r in pysam_iter:
        if not r.is_unmapped:
            all_r.append(r.query_name + '_' + r.get_tag('XM'))
    return all_r


def isolated_umis(pysam_iter):

    reads = initialize_iterator(pysam_iter)
    x = Counter(tag_reads(reads)).keys()
    groups = {}

    for item in x:
        qname, xm = item.split('_')
        if xm in groups:
            groups[xm] += 1
        else:
            groups.update({xm: 1})

    isolated = [key for key, value in groups.iteritems() if value == 1]

    return isolated


def get_umi_list(pysam_iterator, include='mapped'):

    reads = initialize_iterator(pysam_iterator)
    umis = None

    if include == 'unmapped':
        umis = [r.get_tag('XM') for r in reads]
    elif include == 'mapped':
        umis = [r.get_tag('XM') for r in reads if not r.is_unmapped]
    elif include == 'unique':
        umis = [r.get_tag('XM') for r in reads if not r.is_unmapped if r.get_tag('NH') == 1]
    elif include == 'multi':
        umis = [r.get_tag('XM') for r in reads if not r.is_unmapped if not r.get_tag('NH') == 1]

    umis_c = Counter(umis).keys()

    return umis_c


''' Functions for building hash tables for multi_reads and unique reads '''


def update_dt(dt, xm, prop):

    if xm in dt:
        dt[xm].append(prop)


def update_counter(counter, xm, qname):

    if xm in counter:
        if qname in counter[xm]:
            counter[xm][qname] += 1
        else:
            counter[xm].update({qname: 1})


def build_umi_dist_tables(pysam_iterator, unique_umi_index=None, multi_umi_index=None):

    """Function for constructing umi distribution hash tables.
    Takes as input a Pysam iterator.
    Builds separate tables for unique and multi reads.
    Builds an auxiliary table for multi reads.
    Outputs:
    Two dictionaries in the form of {str:[(str, int, int)]},
    with the mapping: {umi: (reference, reference start, number of mismatches)},
    and one dictionary in the form of {str : {str: int}},
    with the mapping: {umi: {query name: frequency}}
    """

    reads = initialize_iterator(pysam_iterator)

    if not unique_umi_index:
        unique_umi_index = get_umi_list(reads, include='unique')

    reads = initialize_iterator(pysam_iterator)

    if not multi_umi_index:
        multi_umi_index = get_umi_list(reads, include='multi')

    u_dt = {}
    for k in unique_umi_index: u_dt.update({k: []})
    m_dt = {}
    for k in multi_umi_index: m_dt.update({k: []})
    cnt_m_dt = {}
    for k in multi_umi_index: cnt_m_dt.update({k: {}})

    reads = initialize_iterator(pysam_iterator)

    for read in reads:

        if not read.is_unmapped:
            xm = read.get_tag('XM')
            qname = read.query_name
            ref = read.reference_name
            pos = read.reference_start
            nm = read.get_tag('NM')
            prop = (ref, pos, nm)
            nh = read.get_tag("NH")

            if nh == 1:
                update_dt(u_dt, xm, prop)
            else:
                update_dt(m_dt, xm, prop)
                update_counter(cnt_m_dt, xm, qname)

    return u_dt, m_dt, cnt_m_dt


def summarize_dist_chapter((chapter_name, dist_chapter)):

    a = [i for (i, j, k) in dist_chapter]
    b = Counter(a)

    return chapter_name, (len(a), len(b), min(b.values()), max(b.values()))


def summarize_dist_table(dist_table):

    """Function for summarizing umi distribution information hash table.
    Accepts as input a dictionary in the form of {str: [(str, int, int)]}.
    Outputs a dictionary in the form of {str: (int, int, int, int)}.
    Each entry in the new dictionary corresponds to the following:
    umi: (number of reads, number of references, smallest read population among regions, biggest read population among regions)
    """

    pool = mp.Pool()

    chapter_names = dist_table.keys()
    chapters = []
    for chapter in chapter_names:
        dct = dist_table[chapter]
        chapters.append(dct)

    summary = dict(pool.map(summarize_dist_chapter, zip(chapter_names, chapters)))

    pool.close()

    return summary


def summarize_cnt_chapter((chapter_name, chapter)):

    return chapter_name, (len(chapter), sum(chapter.values()))


def summarize_cnt_table(cnt_table):

    """Function for summarizing multi read umi statistics information hash table.
    Accepts as input a dictionary in the form of {str: {str: int}}.
    Outputs a dictionary in the form of {str: (int, int)}.
    Each entry in the new dictionary corresponds to the following:
    umi: (number of reads, number of alignments)
    """

    pool = mp.Pool()

    chapter_names = cnt_table.keys()
    chapters = []
    for chapter in chapter_names:
        dct = cnt_table[chapter]
        chapters.append(dct)

    summary = dict(pool.map(summarize_cnt_chapter, zip(chapter_names, chapters)))

    pool.close()

    return summary


def prepare_dist_files(file_names):

    logger.info('This module prepares and saves umi distribution information.')
    logger.info('For each sample two separate distribution tables are generated.')
    logger.info('Tables record unique reads and multireads respectively.\n')

    for in_file, out_file1, out_file2, out_file3 in file_names:
        logger.info('Sample: %s.' % in_file)
        st_time = time.time()
        st = pysam.AlignmentFile(in_file, "rb")
        ud, md, cn_d = build_umi_dist_tables(st)
        pickle_it(ud, out_file1)
        pickle_it(md, out_file2)
        pickle_it(cn_d, out_file3)
        end_time = time.time()
        logger.info('Finished writing umi distribution tables in %d seconds:' % (end_time - st_time))
        logger.info('Saved umi distribution table for unique reads in %s.' % out_file1)
        logger.info('Saved umi distribution table for multi reads in %s.' % out_file2)
        logger.info('Saved multi read counter help table in %s.' % out_file3)


def summarize_dist_tables(file_names):

    logger.info('This module summarizes umi distribution information.')
    logger.info('Compressed summary tables for each file are stored on disk.\n')

    for in_file, out_file in file_names:
        logger.info('Sample: %s.' % in_file)
        st_time = time.time()
        dist_table = pickle.load(open(in_file, 'r'))
        summary = summarize_dist_table(dist_table)
        pickle_it(summary,out_file)
        end_time = time.time()
        logger.info('Finished writing summary table in %d seconds:' % (end_time - st_time))
        logger.info('Table saved in %s.' % out_file)


''' Functions for generating multiple views on hash tables'''


def region_based_generator(dist_table, reference='', includes=True):

    regions = sum(dist_table.values())

    for ref, start in regions:
        if reference in ref:
            if includes:
                yield ref, start
        else:
            if not includes:
                yield ref, start


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    fh = logging.FileHandler('report.log')
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)

    logger.info('Call to Distribution Analyzer module. \n')

    '''
    file_names = []
    for sample in samples:
        in_file = prepare_file_name(sample, pt, 'sample_')
        out_file1 = prepare_file_name(sample, path, 'ut_', '.pkl')
        out_file2 = prepare_file_name(sample, path, 'mt_', '.pkl')
        out_file3 = prepare_file_name(sample, path, 'cnt_', '.pkl')
        file_names.append((in_file, out_file1, out_file2, out_file3))

    prepare_dist_files(file_names)
    '''

    file_names = []
    for sample in samples:

        f1 = prepare_file_name(sample, path, 'ut_', '.pkl')
        f2 = prepare_file_name(sample, path, 'summ_ut_', '.pkl')
        f3 = prepare_file_name(sample, path, 'mt_', '.pkl')
        f4 = prepare_file_name(sample, path, 'summ_mt_', '.pkl')
        file_names.append((f1, f2))
        file_names.append((f3, f4))

    summarize_dist_tables(file_names)


    '''
    file_names = []
    for sample in samples:

        f1 = prepare_file_name(sample, path, 'cnt_', '.pkl')
        f2 = prepare_file_name(sample, path, 'summ_cnt_', '.pkl')
        file_names.append((f1, f2))

    counter = 1

    for in_file, out_file in file_names:

        start_time = time.time()

        cnt_table = pickle.load(open(in_file, 'r'))
        summ_cnt = summarize_cnt_table(cnt_table)
        pickle.dump(summ_cnt, open(out_file, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

        end_time = time.time()
        logger.info('Finished summarize operation %d in %d seconds. \n' %(counter, end_time - start_time))
        counter += 1
    '''




