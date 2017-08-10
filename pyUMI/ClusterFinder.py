#!/usr/bin/env python

#TODO : method documentation
#TODO : exception handling

import argparse
import pysam
import pandas as pd
import logging
import time
import cPickle as pickle

logger = logging.getLogger(__name__)

'''
def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--in_file', '-i', type=str, help='Input BAM file name.', required=True)
    parser.add_argument(
        '--out_dir', '-o', type=str, default="~",
        help='Destination for output file(s).')

    args = parser.parse_args()

    return args
'''


def generate_bam_view(pysamIter):

    pysamIter.reset()
    reads = pysamIter.fetch(until_eof=True)

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
    gaps = gaps_of_size(l, 1500)
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


def generate_cluster_map(df, length=1500):

    cluster_map = {}
    for i, j in df.groupby(['XM', 'Ref']):

        if len(j) > 1:

            locs = list(j['Start'])
            clusters = loci_clusters(locs, length)

            if clusters:
                xm, ref = i
                c_map = assign_cluster(clusters, locs, xm, ref)
                names = list(j['QName'])

                for index, item in enumerate(c_map):
                    if item:
                        props = names[index], locs[index]
                        update_cluster_map(cluster_map, item, props)

    return cluster_map


def generate_umi_group_map(df):

    umi_group_map = {}
    for i, j in df.groupby(['XM']):
        size = len(set(j['QName']))
        umi_group_map.update({i: size})

    return umi_group_map


def generate_maximal_cluster_map(cluster_map, umi_gorup_map):

    maximal_cluster_map = {}
    for key, value in cluster_map.items():

        a = key.split(':')
        xm = a[0]
        gs = umi_gorup_map[xm]
        if gs - 1 <= len(value) <= gs:
            maximal_cluster_map.update({key: value})

    return maximal_cluster_map


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to Cluster Finder module.')

    #params = vars(get_args())
    #output_dir = params['out_dir']
    #in_file = params['in_file']

    in_file = '/data/parastou/Star-Lab/test/NStar25.Aligned.out.tagged.bam'
    st = pysam.AlignmentFile(in_file,"rb")
    logger.info('Alignment file loaded:')
    logger.info('%s' % in_file)

    stt = time.time()

    df = generate_bam_view(st)
    cluster_map = generate_cluster_map(df, 1500)
    logger.info('Cluster map generated')  # ToDo stats
    umi_group_map = generate_umi_group_map(df)
    max_cluster_map = generate_maximal_cluster_map(cluster_map, umi_group_map)

    edt = time.time()
    elapsed = edt - stt

    logger.info('Generated maximal clusters in %f seconds.' % elapsed)  # ToDo stats

    out_file = 'MaximalClusters.pkl'
    pickle.dump(max_cluster_map , open(out_file,'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    logger.info('Maximal clusters map saved in %s' % out_file)







