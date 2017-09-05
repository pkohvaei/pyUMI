#!/usr/bin/env python

import os
import sys
import pysam
import logging
import Bio.SeqIO
import argparse


logger = logging.getLogger(__name__)


def get_args():
    # TODO : docstring

    parser = argparse.ArgumentParser()
    parser.add_argument('--read', '-r', type=str, help='.fastq file containing the reads', required=True)
    parser.add_argument('--umi', '-u', type=str, help='.fastq file containing UMIs', required=True)
    parser.add_argument('--bam', '-b', type=str, help='bam/sam file to be tagged', required=True)
    parser.add_argument('--xtag', '-x', type=str, help='UMI barcode tag', default='XM')
    parser.add_argument('--tmpdir', '-t', type=str, help='Output directory', default='.')
    args = parser.parse_args()

    return args


def rfq_to_names(in_file):
    # TODO : docstring

    names = []
    fq = Bio.SeqIO.parse(in_file, 'fastq')

    for rec in fq:
        names.append(rec.id)

    if not len(names) == len(set(names)):
        logger.info('Warning! Found duplicate query names in the .fastq file.')
        logger.info('This might cause inconsistencies between reads and their unique molecular barcode.')

    return names


def ufq_to_umis(in_file):
    # TODO : docstring

    umis = []
    fu = Bio.SeqIO.parse(in_file, 'fastq')

    for rec in fu:
        umis.append(rec.seq)

    if not len(umis) == len(set(umis)):
        logger.info('Warning! Found duplicate umi entries in the .fastq file.')
        logger.info('This might cause inconsistencies between reads and their unique molecular barcode.')

    return umis
  
    
def utxt_to_umis(in_file):
    # TODO : docstring

    umis = []
    with open(in_file, 'r') as fu:
        for rec in fu:
            umis.append(rec.rstrip())

    if not len(umis) == len(set(umis)):
        logger.info('Warning! Found duplicate umi entries in the .fastq file.')
        logger.info('This might cause inconsistencies between reads and their unique molecular barcode.')

    return umis


def umi_to_name_map(fq_file, umi_file):
    # TODO : docstring

    names = rfq_to_names(fq_file)
    #umis = ufq_to_umis(umi_file)
    umis = utxt_to_umis(umi_file)

    if not len(names)==len(umis):
        logger.info('Error: Size mismatch between input files.')
        logger.info('Program quits.')
        sys.exit()

    un_dict = dict(zip(names, umis))

    return un_dict


def xm_tag_reads(bam_file, fq_file, umi_file, tag='XM', tmp_dir='.'):
    # TODO : docstring

    umi_name_map = umi_to_name_map(fq_file, umi_file)

    st = pysam.AlignmentFile(bam_file, "rb")

    out_file = "%s.tagged.bam" % os.path.splitext(bam_file)[0]
    tmp_file = os.path.join(tmp_dir, 'tmp.bam')
    tagged_reads = pysam.AlignmentFile(tmp_file, 'wb', template=st)

    st.reset()
    reads = st.fetch(until_eof=True)
    for r in reads:
        # TODO : proper re-formatting of tags
        if r.has_tag('RG'):
            r.set_tag('RG', None)
        xm = umi_name_map[r.query_name]
        r.tags += [(tag, xm)]
        tagged_reads.write(r)

    tagged_reads.close()
    pysam.sort("-o", out_file, tmp_file)
    pysam.index(out_file)
    os.remove(tmp_file)
    
    return out_file


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to bam/sam umi-tagger')
    logger.info('This package tags alignments in bam/sam files with their corresponding UMI.')
    

    params = vars(get_args())

    bam_file = params['bam']
    fq_file = params['read']
    umi_file = params['umi']
    xm_tag = params['xtag']
    tmp_dir = params['tmpdir']
    
    logger.info('UMI barcode information is added with the tag : \'%s\'' %xm_tag)

    # TODO : test code

    out_file = xm_tag_reads(bam_file, fq_file, umi_file, tag=xm_tag, tmp_dir=tmp_dir)
    logger.info('Finished tagging the bam file. Generated sorted output : %s' %out_file)
