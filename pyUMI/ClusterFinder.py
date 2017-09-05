#!/usr/bin/env python

#TODO : method documentation
#TODO : exception handling

import argparse
import pysam
import pandas as pd
import logging
import time
import cPickle as pickle
import os

logger = logging.getLogger(__name__)


def get_args():

    """Function for collecting command-line arguments.

    Returns:
        args: The return value. 'argparse.Namespace' dictionary.

    .. _PEP 484:
        https://www.python.org/dev/peps/pep-0484/

    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', '-i', type=str, help='bam/sam input file', required=True)
    parser.add_argument('--window', '-w', type=int, help='Cluster window size (default=1500nct)', default=1500)
    '''
    parser.add_argument('--xtag', '-x', type=str, help='UMI barcode tag', default='XM')
    parser.add_argument('--tmpdir', '-t', type=str, help='Output directory', default='.')
    '''
    args = parser.parse_args()

    return args


def generate_bam_view(pysam_iter):

    """This function accepts a pysam iterator over a bam file and returns a data frame.
    The returned data frame contains information of a subset of fields of each bam record.

    Parameters::

        pysam_ite(str):
            Iterator of type pysam.AlignmentFile.

    Returns:
        pandas data frame:

        The returned data frame contains the following fields of all bam records
        with 'NH' tags greater than 0 (mapped):

        'XM' tag content ('XM' column)
        query_name ('QName' column)
        reference_name ('Ref' column)
        reference_start ('Start' column)
        'NH' tag content ('NH' column)
    """

    pysam_iter.reset()
    reads = pysam_iter.fetch(until_eof=True)

    records = []
    for r in reads:
        if r.get_tag('NH') > 0:
            xm = r.get_tag('XM')
            qn = r.query_name
            ref = r.reference_name
            start = r.reference_start
            nh = r.get_tag('NH')
            records.append((xm, qn, ref, start, nh))
    df = pd.DataFrame(records, columns=['XM', 'QName', 'Ref', 'Start', 'NH'])

    return df


def gaps_of_size(l, length):
    x = sorted(list(set(l)))
    a = x[1:]
    b = x[:-1]

    gaps = []
    for i, j in zip(a, b):
        if i - j > length:
            gaps += [j, i]

    return gaps


def loci_clusters(l, length):
    interval = sorted(l)
    gaps = gaps_of_size(l, length)
    first = interval[0]
    last = interval[-1]
    gaps.append(first)
    gaps.append(last)
    gaps.sort()

    clusters = [(i, j + 50) for i, j in zip(*2 * [iter(gaps)]) if not i == j]

    if clusters:
        return clusters
    else:
        return None


def assign_cluster(clusters, loci, xm, ref):
    prefix = xm + ':' + ref + ':'
    cluster_map = [None] * len(loci)
    for idx, locus in enumerate(loci):
        for i, j in clusters:
            if i <= locus <= j:
                cluster_map[idx] = prefix + str(i) + ':' + str(j)
                break

    return cluster_map


def update_cluster_map(cluster_map, cluster_id, props):

    if cluster_id not in cluster_map:
        cluster_map.update({cluster_id: [props]})
    else:
        cluster_map[cluster_id].append(props)


def generate_cluster_map(df, length):

    cluster_map = {}
    for i, j in df.groupby(['XM', 'Ref']):

        if len(j) > 1:

            locs = list(j['Start'])
            nhs = list(j['NH'])
            clusters = loci_clusters(locs, length)

            if clusters:
                xm, ref = i
                c_map = assign_cluster(clusters, locs, xm, ref)
                names = list(j['QName'])

                for index, item in enumerate(c_map):
                    if item:
                        props = names[index], locs[index], nhs[index]
                        update_cluster_map(cluster_map, item, props)

    return cluster_map


def generate_umi_group_map(df):

    umi_group_map = {}
    for i, j in df.groupby(['XM']):
        size = len(set(j['QName']))
        umi_group_map.update({i: size})

    return umi_group_map


def generate_maximal_cluster_maps(cluster_map, umi_group_map):

    emc_map = {}
    emc_dup_map = {}
    mc_map = {}
    mc_dup_map = {}
    tmp = {}

    for key, value in cluster_map.items():

        a = key.split(':')[0]
        # special case to handle repetitive reads in close vicinity
        s = len(set([i for i, j, k in value]))
        if (a not in tmp) or (a in tmp and s > tmp[a]):
            tmp.update({a: s})

    for key, value in cluster_map.items():

        a = key.split(':')
        xm = a[0]
        gs = umi_group_map[xm]
        s = len(set([i for i, j, k in value]))
        s_max = len([i for i, j, k in value])

        if gs == s:
            if s_max == gs:
                emc_map.update({key: value})
            elif s_max > gs:
                emc_dup_map.update({key: value})
            else:
                logger.info('Strange case of exact maximal cluster: check your code!')
        else:
            m_size = tmp[xm]
            if m_size == s:
                if s_max == m_size:
                    mc_map.update({key: value})
                elif s_max > m_size:
                    mc_dup_map.update({key: value})
                else:
                    logger.info('Strange case of maximal cluster: check your code!')

    return emc_map, mc_map, emc_dup_map, mc_dup_map


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to Cluster Finder module.')

    params = vars(get_args())
    in_file = params['infile']
    window_size = params['window']

    prefix = os.path.splitext(in_file)[0]
    out_file1 = '%s.emc.pkl' % prefix
    out_file2 = '%s.mc.pkl' % prefix
    out_file3 = '%s.emc_dup.pkl' % prefix
    out_file4 = '%s.mc_dup.pkl' % prefix

    logger.info('Calculates exact maximal and maximal clusters for alignment file:')
    logger.info('%s' % in_file)
    logger.info('')

    st = pysam.AlignmentFile(in_file,"rb")

    stt = time.time()
    logger.info('Extracting alignment information ...')
    df = generate_bam_view(st)
    logger.info('-'*40)
    logger.info('%s alignments in total.' % format(len(df), ","))

    logger.info('Finding UMI groups ...')
    umi_group_map = generate_umi_group_map(df)
    logger.info('-'*40)
    logger.info('%s UMI groups in total.' % format(len(umi_group_map), ','))

    logger.info('Generating clusters ...')
    cluster_map = generate_cluster_map(df, window_size)
    logger.info('-'*40)
    logger.info('%s clusters in total.' % format(len(cluster_map), ','))

    logger.info('Generating maximal clusters ...')
    emc_map, mc_map, emc_dup_map, mc_dup_map = generate_maximal_cluster_maps(cluster_map, umi_group_map)
    logger.info('-'*40)
    logger.info('%s exact maximal clusters in total.' % format(len(emc_map), ','))
    logger.info('%s exact maximal clusters with duplicate reads in total.' % format(len(emc_dup_map), ','))
    logger.info('%s maximal clusters in total.' % format(len(mc_map), ','))
    logger.info('%s maximal clusters with duplicate in total.' % format(len(mc_dup_map), ','))
    logger.info('')
    edt = time.time()
    elapsed = edt - stt

    logger.info('Elapsed time: %f seconds.' % elapsed)  # ToDo stats

    pickle.dump(emc_map, open(out_file1,'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(mc_map, open(out_file2, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(emc_dup_map, open(out_file3, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(mc_dup_map, open(out_file4, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    logger.info('Exact maximal cluster map saved in %s' % out_file1)
    logger.info('Maximal cluster map saved in %s' % out_file2)
    logger.info('Exact maximal cluster map with duplicates saved in %s' % out_file3)
    logger.info('Maximal cluster map with duplicates saved in %s' % out_file4)
    logger.info('Done!')
    logger.info('')







