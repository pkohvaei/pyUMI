#!/usr/bin/env python

import logging
import matplotlib.pyplot as plt
from collections import Counter

logger = logging.getLogger(__name__)


def count_reads(alignmentFile):
    alignmentFile.reset()
    return alignmentFile.count(until_eof=True)


def count_um(alignmentFile):
    uniques = 0
    multimaps = 0

    alignmentFile.reset()
    reads = alignmentFile.fetch(until_eof=True)
    for read in reads:
        if not(read.is_unmapped) and read.has_tag("NH"):
            if read.get_tag("NH") > 1:
                multimaps += 1
            if read.get_tag("NH") == 1:
                uniques += 1
    return multimaps,uniques


def multimap_generator(alignmentFile):

    alignmentFile.reset()
    reads = alignmentFile.fetch(until_eof=True)
    for read in reads:
        if not(read.is_unmapped) and read.has_tag("NH"):
            if read.get_tag("NH") > 1:
                yield read


def uniques_generator(alignmentFile):

    alignmentFile.reset()
    reads = alignmentFile.fetch(until_eof=True)
    for read in reads:
        if not(read.is_unmapped) and read.has_tag("NH"):
            if read.get_tag("NH") == 1:
                yield read


def umapped_generator(alignmentFile, mapped=True):

    alignmentFile.reset()
    reads = alignmentFile.fetch(until_eof=True)

    for read in reads:
        if not(read.is_unmapped) and mapped:
            yield read
        elif read.is_unmapped and not(mapped):
            yield read


def mapping_stats(alignmentFile, draw_pie=False):

    alignmentFile.reset()

    mapped = alignmentFile.mapped
    unmapped = alignmentFile.unmapped
    multimaps,uniques = count_um(alignmentFile)
    total =  mapped + unmapped

    logger.info("\t"+"Total number of reads:"+"\t"+"%s" %format(total,",").rjust(15) )
    logger.info("\t"+"Unmapped reads:"+"\t\t"+"%s" %format(unmapped,",").rjust(15))
    logger.info("\t"+"Uniquely mapped reads:"+"\t"+"%s" %format(uniques,",").rjust(15))
    logger.info("\t"+"Multimapped reads:"+"\t"+"%s" %format(multimaps,",").rjust(15))

    if draw_pie:
        labels = "unmapped", "uniquely mapped", "multimapped"
        plt.pie([unmapped,uniques,multimaps], labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)


def cellbarcode_info(alignmentFile, report_info=False, draw_pie=False):

    alignmentFile.reset()
    reads = alignmentFile.fetch(until_eof=True)

    cell_barcodes = []
    for r in reads:
        cell_barcodes.append(r.get_tag("XC"))

    cbs = Counter(cell_barcodes)

    if report_info:
        logger.info("  Total number of cellular barcodes: %d" % len(cbs))
        logger.info("  Barcode frequency:")
        for cb in cbs:
            logger.info("  %s: %s" %(cb,format(cbs[cb],",")))

    if draw_pie:
        labels = [label for label in cbs.keys()]
        plt.pie([x for x in cbs.values()], labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)

    return cbs


def flag_stats(alignmentFile, report_info=False, draw_bar=False):

    alignmentFile.reset()
    reads = alignmentFile.fetch(until_eof=True)

    flags = []
    for r in reads:
        flags.append(r.flag)

    cf = Counter(flags)

    if report_info:
        logger.info("  Present flags: %s" %cf.keys())

    if draw_bar:
        plt.bar(range(len(cf)), cf.values(), .2, color='g', align='center')
        plt.xticks(range(len(cf)), cf.keys());
        plt.yticks(cf.values(), [format(x, ",") for x in cf.values()])

    return cf


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to UMI statistics module.')
