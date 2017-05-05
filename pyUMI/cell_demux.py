#!/usr/bin/env python

#TODO : method documentation
#TODO : exception handling
#TODO : consider the replacement of generators for datastructures

import argparse
import logging
import time

import pysam
import numpy as np

from bamStats import cellbarcode_info


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



if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to Cell DeMux module.')

    params = vars(get_args())

    output_dir = params['out_dir']
    in_file = params['in_file']

    #TODO : Check the out_dir

    alignment_file = pysam.AlignmentFile(in_file,"rb")
    aligned_segments = alignment_file.fetch(until_eof=True)

    s_time1 = time.time()
    cell_barcodes = cellbarcode_info(alignment_file)
    e_time1 = time.time()

    logger.info("Total number of samples: %d" %len(cell_barcodes))
    logger.info("Elapsed time: %d seconds" %(e_time1-s_time1))

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

    s_time2 = time.time()
    for cell_barcode in cell_barcodes:
        alignment_file.reset()
        aligned_segments = alignment_file.fetch(until_eof=True)
        sample_file_name = 'sample_' + cell_barcode + ".bam"
        demux_reads = pysam.AlignmentFile(sample_file_name, 'wb', template=alignment_file)
        for aligned_segment in aligned_segments:
            xc_tag = aligned_segment.get_tag("XC")
            if xc_tag == cell_barcode:
                demux_reads.write(aligned_segment)
        demux_reads.close()
        pysam.index(sample_file_name)

    e_time2 = time.time()
    alignment_file.close()

    #TODO : extensive report on dir, name and size of the produced files

    logger.info("De-multiplexing finish")
    logger.info("Elapsed time: %d seconds" %(e_time2-s_time2))




