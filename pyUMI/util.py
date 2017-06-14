#!/usr/bin/env python

import logging
import matplotlib.pyplot as plt
import numpy as np
import time
import pysam
from subprocess import check_output
import cPickle as pickle


logger = logging.getLogger(__name__)


def initialize_iterator(pysam_iterator):

    if type(pysam_iterator) == pysam.calignmentfile.AlignmentFile:
        pysam_iterator.reset()
        reads = pysam_iterator.fetch(until_eof=True)
    else:
        reads = pysam_iterator

    return reads


def prepare_file_name(sample, path='', prefix='', extn='.bam'):
    return path + prefix + sample + extn


def pickle_it(obj, out_file):

    # TODO : try and exception handling
    pickle.dump(obj, open(out_file, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)


def pcr_duplicates(in_file):

    flags = [1024, 1089, 1097, 1105, 1107, 1121, 1123, 1137, 1145, 1153, 1161, 1169, 1171, 1185, 1187, 1201, 1209, 4]
    duplicates = 0
    cmd = "samtools view -f "
    opt = " | wc -l"

    st_time = time.time()
    for dup_flag in flags:
        command = cmd + str(dup_flag) + ' ' + in_file + opt
        counts = int(check_output(command, shell=True))
        logger.info('Found %s PCR-duplicates with flag %d.' %(format(counts,','),dup_flag))
        duplicates += counts

    end_time = time.time()

    logger.info('Operation took %s seconds.' %format(end_time - st_time,","))
    logger.info("Total number of PCR-duplicates: %s" %format(duplicates,","))

    return duplicates


def rescale( x, oMin, oMax, nMin, nMax ):

    #range check
    if oMin == oMax:
        print "Warning: Zero input range"
        return None

    if nMin == nMax:
        print "Warning: Zero output range"
        return None

    #check reversed input range
    reverseInput = False
    oldMin = min( oMin, oMax )
    oldMax = max( oMin, oMax )
    if not oldMin == oMin:
        reverseInput = True

    #check reversed output range
    reverseOutput = False
    newMin = min( nMin, nMax )
    newMax = max( nMin, nMax )
    if not newMin == nMin :
        reverseOutput = True

    portion = (x-oldMin)*(newMax-newMin)/(oldMax-oldMin)
    if reverseInput:
        portion = (oldMax-x)*(newMax-newMin)/(oldMax-oldMin)

    result = portion + newMin
    if reverseOutput:
        result = newMax - portion

    return result


def plot_bar(data1, data2=None, bar_title=None, xlabel=None, ylabel=None, size=(20, 8), log_scale=True):

    fig = plt.figure(figsize=size)
    ax = fig.add_subplot(111)
    plt.grid()
    ax = plt.subplot(111)
    #plt.title(bar_title, y=1.08, fontsize=14)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.ylim(0, 50000)
    plt.xlim(data1[0], data1[-1])
    plt.xticks(np.arange(data1[0],data1[-1], 50))

    if data2 is not None:
        #data2_x = [3 * x - 2 for x in range(len(data2))]
        ax.bar(data1, data2, width=.5, color='b')
    if log_scale:
        plt.yscale('log')

    plt.show()


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to util package.')