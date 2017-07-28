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
    # TODO : args definition and default assignment

    parser = argparse.ArgumentParser()
    parser.add_argument('--species', '-s', type=str, help='Name of UCSC database (e.g. hg38)', default=None)
    args = parser.parse_args()

    return args


def rfq_to_names(in_file):
    # TODO : docstring

    names = []
    fq = Bio.SeqIO.parse(in_file, 'fastq')

    for rec in fq:
        names.append(rec.id)

    return names


def ufq_to_umis(in_file):
    # TODO : docstring

    umis = []
    fu = Bio.SeqIO.parse(in_file, 'fastq')

    for rec in fu:
        umis.append(rec.seq)

    return umis


def umi_to_name_map(fq_file, umi_file):
    # TODO : docstring

    names = rfq_to_names(fq_file)
    umis = ufq_to_umis(umi_file)

    if not len(names)==len(umis):
        logger.info('Error: Size mismatch between input files.')
        logger.info('Program quits.')
        sys.exit()

    un_dict = dict(zip(names, umis))

    return un_dict


def xm_tag_reads(bam_file, fq_file, umi_file):
    # TODO : docstring

    umi_name_map = umi_to_name_map(fq_file, umi_file)

    st = pysam.AlignmentFile(bam_file, "rb")

    out_file = "%s.tagged.bam" % os.path.splitext(bam_file)[0]
    tagged_reads = pysam.AlignmentFile('tmp.bam', 'wb', template=st)

    st.reset()
    reads = st.fetch(until_eof=True)
    for r in reads:
        # TODO : proper re-formatting of tags
        if r.has_tag('RG'):
            r.set_tag('RG', None)
        xm = umi_name_map[r.query_name]
        r.tags += [('XM', xm)]
        tagged_reads.write(r)

    tagged_reads.close()
    pysam.sort("-o", out_file, 'tmp.bam')
    pysam.index(out_file)


if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)
    logger.info('Call BumiLabelMaker.')
    # TODO : logger help string

    params = vars(get_args())

    # TODO : vars assignments and checking
    # TODO : test code

    xm_tag_reads()
