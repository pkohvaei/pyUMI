#!/usr/bin/env python

import logging
from collections import Counter
from pylab import *
import pickle
import time
import multiprocessing as mp

import pysam


from util import prepare_file_name, initialize_iterator


samples = ['GACCGC', 'AAAACT', 'GGCGTC', 'AAAGTT', 'GTTCGA', 'ATATAG',
           'TAAAGT', 'ATCAAA', 'TCTGCA' , 'CCCTGG', 'TTAATC', 'CCGGAC']

path = '/data/UMI/data/MUS/DT/'
pt = '/data/UMI/data/MUS/'

umi_length = 10


logger = logging.getLogger(__name__)


'''Functions for visualization'''


def bar_prep(umi_dist):

    sorted_dist = sorted(umi_dist.items())
    x = [i for (i, j) in sorted_dist]
    y = [j for (i, j) in sorted_dist]
    return x, y


def umi_dist(pysam_iterator):

    reads = initialize_iterator(pysam_iterator)

    u1 = []
    u2 = {}

    for r in reads:
        if not r.is_unmapped:
            nh = r.get_tag('NH')
            xm = r.get_tag('XM')
            qname = r.query_name

            if nh == 1:
                u1.append(xm)
            else:
                if xm in u2:
                    u2[xm].append(qname)
                else:
                    u2.update({xm:[qname]})

    uniq_dist = Counter(Counter(u1).values())

    B = {}
    C = {}
    for item in u2:
        B.update({item: len(u2[item])})
        C.update({item: len(np.unique(u2[item]))})

    map_dist = Counter(Counter(B).values())
    multi_dist = Counter(Counter(C).values())

    return uniq_dist, multi_dist, map_dist


def prepare_plot_data(samples):

    bar_plot_data = {}

    for sample in samples:

        in_file = prepare_file_name(sample, pt, 'sample_','.bam')
        st = pysam.AlignmentFile(in_file, "rb")

        uniq_dist, multi_dist, _ = umi_dist(st)

        x1, y1 = bar_prep(uniq_dist)
        x2, y2 = bar_prep(multi_dist)

        bar_plot_data.update({sample:(x1, y1, x2, y2)})


def plot_bars(bar_plot_data):

    number_of_subplots = len(bar_plot_data)

    plt.close('all')

    l1 = plt.axhline(y=0.8, color='b')
    l2 = plt.axhline(y=0.8, color='r')
    #l3 = plt.axhline(y=0.5, color='g')
    l1.remove()
    l2.remove()
    #l3.remove()

    fig = plt.figure(figsize=(10,80))
    fig.suptitle('Unique Molecular Identifier distribution per single cell', fontsize=13, fontweight='bold')

    for i, sample in enumerate(bar_plot_data):
        data1, data2, data3, data4, _ = bar_plot_data[sample]
        i += 1
        plt.grid()
        plt.ylim(0, 30000)
        plt.xlim(0, 50)
        plt.legend([l1, l2], ["Unique Reads", "Multi Reads"], loc=4)
        plt.tight_layout(pad=4, w_pad=0.5, h_pad=10)
        ax1 = subplot(number_of_subplots, 1, i)
        ax1.bar(data1, data2, width=0.8, color='b', align='center')
        ax1.bar(data3, data4, width=0.8, color='r', align='edge')
        #ax1.bar(data5, data6, width=0.8, color='g')

        ax1.set_xlabel('Frequency of reads', fontsize=11)
        ax1.set_ylabel('Number of distinct UMIs (' + sample + ')', fontsize=11)

    plt.show()


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)

    fh = logging.FileHandler('report.log')
    fh.setLevel(logging.INFO)
    logger.addHandler(fh)

    logger.info('Call to Distribution Analyzer module. \n')

    '''
    pickle.dump(bar_plot_data, open('barPlotData.pkl', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    '''

    bar_plot_data = pickle.load(open('barPlotData.pkl', 'r'))

    plot_bars(bar_plot_data)


