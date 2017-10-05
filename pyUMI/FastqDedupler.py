#!/usr/bin/env python

import os
import logging
import argparse
import cPickle as pickle

import pysam
import numpy as np
import pandas as pd
import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


logger = logging.getLogger(__name__)


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--in_file', '-f', type=str, help='Reads input file (fastq).', required=True)
    parser.add_argument(
        '--umi_in_file', '-u', type=str, help='UMIs input file (txt).', required=True)
    parser.add_argument(
        '--out_dir', '-o', type=str, default='', help='Output directory (default . ).')
    parser.add_argument(
        '--prefix', '-p', type=str, default='', help='Output files prefix.'
    )

    args = parser.parse_args()

    return args


def fastq_umi_to_df(in_file, in_umi):

    fastq_records = []
    for record in SeqIO.parse(in_file, 'fastq'):
        r = (record.id, str(record.seq), record.letter_annotations["phred_quality"])
        fastq_records.append(r)

    umis = []
    with open(in_umi, 'r') as fu:
        for rec in fu:
            umis.append(rec.rstrip())

    records = []
    for index, (i, j, k) in enumerate(fastq_records):
        records.append((i, j, k, umis[index]))

    x = pd.DataFrame(records, columns=['NAME', 'SEQ', 'QUAL', 'XM'])

    return x


def modify_qual(quals):

    a = np.vstack(quals)
    new_qual = np.amax(a, axis=0)

    return new_qual


def to_IOSeq_rec(row):

    qname = row['NAME']
    seq = row['SEQ']
    qual = row['QUAL']
    xm = row['XM']

    record = SeqRecord(Seq(seq, generic_dna), id=qname, name=qname, description='', dbxrefs=[])
    record.letter_annotations["phred_quality"] = qual

    return record, xm


def dedup_df_records(df):

    new_records = []
    dup_count = 0
    rec_count = 0

    for (xm, seq), g in df.groupby(['XM', 'SEQ']):

        if len(g) > 1:
            name = list(g['NAME'])[0] + 'M'
            qual = modify_qual(list(g['QUAL']))
            new_records.append((name, seq, qual, xm))
            dup_count += (len(g) - 1)
        else:
            new_records.append((g['NAME'].item(), seq, g['QUAL'].item(), xm))

        rec_count += 1

    logger.info('')
    logger.info('Number of duplicate records:\t\t%s ' % format(dup_count, ','))
    logger.info('Total records in deduplicated file:\t%s' % format(rec_count, ','))
    logger.info('Checksum (duplicates + new):\t\t%s' % format(dup_count + rec_count, ','))
    logger.info('')

    return new_records


def prepare_out_names(in_file, in_umi, out_dir='.', prefix=''):

    out_file = os.path.join(out_dir, prefix + '%s.dedup.fastq' % os.path.basename(in_file).split('.')[0])
    out_umi = os.path.join(out_dir, prefix + '%s.dedup.umi' % os.path.basename(in_umi).split('.')[0])

    return out_file, out_umi


def df_to_file(df, out_file, out_umi):

    with open(out_file, 'w') as fq:
        with open(out_umi, 'w') as fu:
            for index, row in df.iterrows():
                record, xm = to_IOSeq_rec(row)

                SeqIO.write(record, fq, 'fastq')
                fu.write(xm + '\n')

    fu.close()
    fq.close()

    logger.info('De-duplicated fastq file saved in: %s' % out_file)
    logger.info('Modified umi file saved in: %s' % out_umi)
    logger.info('-' * 40)


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to fastq de-duplicator module.')
    logger.info('')
    logger.info('This module removes PCR duplicates from fastq file and modifies the associated UMI file.')
    logger.info('-' * 40)

    params = vars(get_args())

    output_dir = params['out_dir']
    in_file = params['in_file']
    in_umi = params['umi_in_file']
    prefix = params['prefix']

    x = fastq_umi_to_df(in_file, in_umi)
    y = pd.DataFrame(dedup_df_records(x), columns=['NAME', 'SEQ', 'QUAL', 'XM'])

    out_file, out_umi = prepare_out_names(in_file, in_umi, output_dir, prefix)
    df_to_file(y, out_file, out_umi)

    logger.info('De-dupler finished successfully!')

