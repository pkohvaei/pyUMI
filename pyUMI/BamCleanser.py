#!/usr/bin/env python

import argparse
import logging
import time

import pysam
import numpy as np
from util import pickle_it
import cPickle as pickle
import pandas as pd

from umiDistAnalyzer import get_umi_list


samples = ['AAAGTT']#,'ATATAG','ATCAAA','CCCTGG','CCGGAC','GACCGC','GGCGTC','GTTCGA','TAAAGT','TCTGCA','TTAATC']


logger = logging.getLogger(__name__)


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--in_file', '-i', type=str, help='Input BAM file name.', required=True)
    parser.add_argument(
        '--out_dir', '-o', type=str, default="~",
        help='Destination for output file(s).')

    args = parser.parse_args()

    return args


def low_phred_base(phred_array, threshold=10):

    if np.sum(phred_array < threshold) > 0:
        return True
    else:
        return False


def bam_to_deduplicated_fastq(in_file, out_file1, out_file2):

    st = pysam.AlignmentFile(in_file, "rb")
    reads = st.fetch(until_eof=True)

    fx_records = []
    fxk_records = []
    for r in reads:
        if low_phred_base(r.query_qualities):
            fxk_records.append((r.query_name, r.qual, r.get_tag('XM') + ':' + r.query_sequence, ''.join(map(str, r.query_qualities))))
        else:
            fx_records.append((r.query_name, r.qual, r.get_tag('XM') + ':' + r.query_sequence))

    df = pd.DataFrame(fx_records, columns=['N', 'Q', 'XS'])
    dfk = pd.DataFrame(fxk_records, columns=['N', 'Q', 'XS', 'F'])

    # df.drop_duplicates(subset=['XS']).drop_duplicates(subset=['Q'])
    # dfk.drop_duplicates(subset=['F'])

    with open(out_file1, 'a+') as fs:
        with open(out_file2, 'a+') as fu:
            for index, row in df.drop_duplicates(subset=['XS']).drop_duplicates(subset=['Q']).iterrows():

                n, q, xs = row
                x, s = str.split(xs, ':')
                fs.write('@%s\n%s\n+\n%s\n' % (n, s, q))
                fu.write('%s\n' % x)

            for index, row in dfk.drop_duplicates(subset=['XS', 'F']):

                n, q, xs, _ = row
                x, s = str.split(xs, ':')
                fs.write('@%s\n%s\n+\n%s\n' % (n, s, q))
                fu.write('%s\n' % x)


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to bam Cleanser module.')

    # params = vars(get_args())
    # output_dir = params['out_dir']
    # in_file = params['in_file']
    # TODO : Check the out_dir

    for sample in samples:

        in_file = '/data/UMI/data/MUS/sample_' + sample + '.bam'
        # out_file = '/data/UMI/data/MUS/sample_' + sample + 'umis.pkl'
        out_file1 = '/data/UMI/data/MUS/sample_' + sample + '.fastq'
        out_file2 = '/data/UMI/data/MUS/sample_' + sample + '.umi'

        # st = pysam.AlignmentFile(in_file, "rb")
        # umis = get_umi_list(st)
        # pickle_it(umis, out_file)

        bam_to_deduplicated_fastq(in_file, out_file1, out_file2)

    '''
    s_time1 = time.time()
    e_time1 = time.time()

    logger.info("Total number of samples: %d" %len(cell_barcodes))
    logger.info("Elapsed time: %d seconds" %(e_time1-s_time1))
    '''
    '''
    # TODO : work on this further for parallel computing version
    for cell_barcode in cell_barcodes:

        locals()["sample_%s" %cell_barcode] = []

    alignment_file.reset()
    aligned_segments = alignment_file.fetch(until_eof=True)
    s_time2 = time.time()
    for aligned_segment in aligned_segments:
        xc_tag = aligned_segment.get_tag("XC")
        locals()["sample_%s" % xc_tag].append(aligned_segment)

    for i, cell_barcode in enumerate(cell_barcodes):
        sample_file_name = 'sample_' + cell_barcode + ".bam"
        demux_reads = pysam.AlignmentFile(sample_file_name, 'wb', template=alignment_file)
        for s in locals()["sample_%s" %cell_barcode]:
            demux_reads.write(s)
        demux_reads.close()
    e_time2 = time.time()
    '''

    '''
    s_time2 = time.time()
    cb_counter = 1
    for cell_barcode in cell_barcodes:
        logger.info('Extracting sample no. %d with cell barcode %s :' %(cb_counter, cell_barcode))
        alignment_file.reset()
        aligned_segments = alignment_file.fetch(until_eof=True)
        sample_file_name = output_dir + 'sample_' + cell_barcode + ".bam"
        demux_reads = pysam.AlignmentFile(sample_file_name, 'wb', template=alignment_file)
        for aligned_segment in aligned_segments:
            xc_tag = aligned_segment.get_tag("XC")
            if xc_tag == cell_barcode:
                demux_reads.write(aligned_segment)
        demux_reads.close()
        pysam.index(sample_file_name)
        cb_counter += 1

    e_time2 = time.time()
    alignment_file.close()

    #TODO : extensive report on dir, name and size of the produced files

    logger.info("De-multiplexing finish")
    logger.info("Elapsed time: %d seconds" %(e_time2-s_time2))
    '''



