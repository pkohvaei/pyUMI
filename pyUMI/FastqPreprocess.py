#!/usr/bin/env python

import pysam
import numpy as np
import cPickle as pickle
from collections import Counter
from operator import itemgetter
from os import path
import pandas as pd
from Bio import SeqIO
import logging
import argparse

logger = logging.getLogger(__name__)


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--sequence_in_file', '-S', type=str, help='Reads fastq file', required=True)
    parser.add_argument(
        '--umi_in_file', '-U', type=str, help='UMIs fastq file', required=True)
    parser.add_argument(
        '--out_dir', '-o', type=str, default='', help='Directory for output file(s).')
    parser.add_argument(
        '--prefix', '-p', type=str, default='', help='Output files prefix.'
    )

    args = parser.parse_args()

    return args


def umi_fastq_to_df(in_file):

    records = []
    for record in SeqIO.parse(in_file, 'fastq'):
        r = str(record.seq)
        records.append(r)
    x = pd.DataFrame(records, columns=['SEQ'])

    return x


def seq_fastq_to_df(in_file):

    records = []
    for record in SeqIO.parse(in_file, 'fastq'):
        r = (str(record.seq), record.letter_annotations["phred_quality"])
        records.append(r)
    x = pd.DataFrame(records, columns=['SEQ', 'QUAL'])

    return x


def N_indices(record_df, column_name='SEQ'):

    idx = record_df[record_df[column_name].str.contains('N', na=False)].index.tolist()

    return idx


def write_fastq(records_iter, out_file):
    with open(out_file, 'w+') as f:
        SeqIO.write(records_iter, f, 'fastq')


'''
    records = []
    for read in reads:
        records.append((read.query_name, read.query_sequence, str(read.qual)))
    rdf = pd.DataFrame(records, columns=['N', 'S', 'Q'])
    rdf2 = rdf.drop_duplicates(subset=['N'])
    pickle.dump(rdf2, open('rdf.pkl', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    with open('/home/parastou/mm10/downsample.fastq', 'w+') as f:
        for index, row in rdf2.iterrows():
            n, s, q = row
            f.write('@%s\n%s\n+\n%s\n'  % (n, s, q))
'''



def remove_N_associated_records(umi_df, seq_df):

    idx = N_indices(umi_df)
    new_umi_df = umi_df.drop(umi_df.index[idx])
    new_seq_df = seq_df.drop(seq_df.index[idx])

    logger.info('Low quality UMIs marked successfully.')

    return new_umi_df, new_seq_df


def preprocess_UMIs(u_in_file, s_in_file, prefix='', out_dir=''):

    umi_df = umi_fastq_to_df(u_in_file)
    logger.info('UMI fastq file %s loaded into memory.' % u_in_file)
    seq_df = seq_fastq_to_df(s_in_file)
    logger.info('Reads fastq file %s loaded into memory.' % s_in_file)
    new_umi_df, new_seq_df = remove_N_associated_records(umi_df, seq_df)

    umi_records = new_umi_df['REC'].tolist()
    seq_records = new_seq_df['REC'].tolist()

    print len(umi_df)
    print len(seq_df)
    print len(umi_records)
    print len(seq_records)

    umi_out_file = path.join(out_dir, prefix, 'umi.fastq')
    seq_out_file = path.join(out_dir, prefix, 'seq.fastq')

    #write_fastq(umi_records, umi_out_file)
    #write_fastq(seq_records, seq_out_file)


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('This module pre-processes reads and their associated UMIs.')
    logger.info('Accepts as input separate fastq files for UMIs and reads.')
    logger.info('Outputs modified fastq files for UMIs and associated reads respectively')
    logger.info('-' * 40)

    params = vars(get_args())

    output_dir = params['out_dir']
    seq_in_file = params['sequence_in_file']
    umi_in_file = params['umi_in_file']
    prefix = params['prefix']

    preprocess_UMIs(umi_in_file, seq_in_file)



