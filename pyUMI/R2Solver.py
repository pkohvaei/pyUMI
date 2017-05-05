#!/usr/bin/env python

import logging
import numpy as np
from collections import Counter
from operator import itemgetter

import MSolverLib as mr
from util import plot_bar


def format_columns(text,width): return width - len(text)


logger = logging.getLogger(__name__)


def multimap_multiregion_generator(multimap_hashtable, variance, no_gene= True):

    """
    Generator for multi_mapped alignments with no gene annotation information in any of their occurrences.

    Parameters
    ----------
    multimap_hashtable : dictionary
    variance : positive integer

    return
    ----------
    tuple of form (query name, xc, xm, list of (ref,coordinate) tuples)
    """

    logger.info("\t" + "Generator for non-gene-tagged multi_mapped reads with at least %d distinct reference names." %(variance+1))
    for read in multimap_hashtable:
        genes = [x[0] for x in multimap_hashtable[read][2] if x[0]]
        if not genes and no_gene:
            refs = np.unique(np.array([x[1] for x in multimap_hashtable[read][2]]))
            if len(refs) > variance:
                xc = multimap_hashtable[read][0]
                xm = multimap_hashtable[read][1]
                regions = [(j,k) for (i,j,k) in multimap_hashtable[read][2]]
                yield read, xc, xm, regions
        elif not no_gene:
            refs = np.unique(np.array([x[1] for x in multimap_hashtable[read][2]]))
            if len(refs) > variance:
                xc = multimap_hashtable[read][0]
                xm = multimap_hashtable[read][1]
                regions = [(j,k) for (i,j,k) in multimap_hashtable[read][2]]
                yield read, xc, xm, regions


def find_region_hits(uniques, matches, regions, window):

    """
    Comment
    """

    hits = []
    resolved = []
    for match in matches:
        if not uniques[match][0]:
            _, ref, cord = uniques[match]
            region = (ref, cord)
            in_vicin = in_vicinity(region, regions, window)
            if in_vicin:
                hits.append((match,uniques[match]))
                resolved += in_vicin
    return hits, resolved


def in_vicinity(region, regions, window):

    resolved = []
    reference, coordinate = region
    for ref, cord in regions:
        if ref == reference and abs(cord - coordinate) <= window:
            resolved.append((ref, cord))
    return resolved


def coordinate_distance(region, regions, window, absolute=True):

    distances = []
    reference, coordinate = region
    for ref, cord in regions:
        if ref == reference and abs(cord - coordinate) <= window:
            if absolute:
                distances.append(abs(cord - coordinate))
            else:
                distances.append(cord - coordinate)
    return distances


def find_hit_distances(uniques, matches, regions, window, absolute=True):

    distances = []

    for match in matches:
        if not uniques[match][0]:
            _, ref, cord = uniques[match]
            region = (ref, cord)
            ds = coordinate_distance(region, regions, window, absolute)
            if ds:
                distances += ds
    return distances


def rank2_multimap_solver(mm_generator=None, uniques_hashtable=None, uniques=None, window=1500, report=False):

    """
    Read mapping method for multi_mapped alignments with no gene annotation information in any of their occurrences.
    Selects the occurrence of a multi_mapped alignment set which has the shortest distance to a uniquely mapped read
    with the same umi. The distance should be less than the window size.

    Parameters
    ----------
    mm_generator : generator of type multimap_multiregion_generator
    uniques_hashtable : dictionary
    uniques : dictionary
    window : positive integer
             Maximum acceptable coordinate distance.
    report : Boolean
             Enables information logging.

    return
    ----------
    pair of lists.
    list1 : tuple of form (query name, xc, xm, list of (ref,coordinate) tuples)
    list2 : list of lists in form of (query_name, (None,ref,coordinate))
    """

    rank2_solvers = []
    resolved_mm = []
    total = 0
    r_total = 0

    for mm in mm_generator:
        total += 1
        qname, xc, xm, regions = mm
        uniq_matches = mr.find_barcode_matches(uniques_hashtable, xc, xm)
        if len(uniq_matches) > 0:
            hits, resolved = find_region_hits(uniques, uniq_matches, regions, window)
            if len(hits) == 1:
                rank2_solvers.append(hits)
                resolved_mm.append((qname, xc, xm, resolved))
                r_total += len(regions)

    if report:
        b = "Number of acceptable multi read assignments to a unique genomic reference : %s"
        c = "Total number of alignments of assigned multi reads: %s"

        logger.info("\t" + "Finished multi read assignment using a vicinity threshold of " + str(window) + " bps.")
        logger.info("\t" + b %format(len(resolved_mm),",").rjust(format_columns(b, 85)))
        logger.info("\t" + c % format(r_total, ",").rjust(format_columns(c, 85)))

    return resolved_mm, rank2_solvers, r_total


def r2_solvers_distance_dist(mm_generator=None, uniques_hashtable=None, uniques=None, window=1500, absolute=True):

    total = 0
    distances = []

    for mm in mm_generator:
        total += 1
        qname, xc, xm, regions = mm
        uniq_matches = mr.find_barcode_matches(uniques_hashtable, xc, xm)
        if len(uniq_matches) > 0:
            ds = find_hit_distances(uniques, uniq_matches, regions, window, absolute)
            if ds:
                distances += ds

    return Counter(distances)


def draw_pie(dist,bar_title,absolute=True, xlabel=None, ylabel=None, size=(10, 4), log_scale=True):

    #TODO: Take care of absolute

    sorted_dist = sorted(dist.items())
    x = [i for (i, j) in sorted_dist]
    y = [j for (i, j) in sorted_dist]

    plot_bar(x, y, bar_title, xlabel, ylabel, size, log_scale)


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to rank2 multimap solver package.')