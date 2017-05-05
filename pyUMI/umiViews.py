#!/usr/bin/env python

import logging
import pysam

from util import initialize_iterator


logger = logging.getLogger(__name__)


def tag_based_generator(pysam_iterator, tag=None, values=None, include=True):

    reads = initialize_iterator(pysam_iterator)

    for r in reads:
        if r.has_tag(tag):
            if r.get_tag(tag) in values:
                if include:
                    yield r
            else:
                if not include:
                    yield r


def field_based_generator(pysam_iterator, field=None, values=None, include=True):

    reads = initialize_iterator(pysam_iterator)

    for r in reads:
        value = getattr(r, field)
        if value in values:
            if include:
                yield r
        else:
            if not include:
                yield r


def tag_annotated_generator(pysam_iterator, tag='', includes=True):

    """Generator of type pysam.alignedSegment.
       Takes an iterable of type pysam.alignedSegment as input.
       Yields tag annotated/non-annotated reads (annotate=True/False)."""

    reads = initialize_iterator(pysam_iterator)

    for r in reads:
        if r.has_tag(tag) and includes:
            yield r
        elif not r.has_tag(tag) and not includes:
                yield r


def total_alignments(pysam_iterator, include='mapped'):

    reads = initialize_iterator(pysam_iterator)
    total = None

    if include == 'unmapped':
        total = len([r for r in reads])
    elif include == 'mapped':
        total = len([r for r in reads if not r.is_unmapped])
    elif include == 'unique':
        total = len([r for r in reads if not r.is_unmapped if r.get_tag('NH') == 1])
    elif include == 'multi':
        total = len([r for r in reads if not r.is_unmapped if not r.get_tag('NH') == 1])

    return total


'''
def locus_pre_processor(pysam_iterator, nMin=0, nMax=3000):

    if type(pysam_iterator) == pysam.calignmentfile.AlignmentFile:
        pysam_iterator.reset()
        reads = pysam_iterator.fetch(until_eof=True)
    else:
        reads = pysam_iterator

    coordinate_map = [(item.reference_name, item.reference_start) for item in reads]

    references = np.unique([i for i,j in coordinate_map])
    if len(references) > 1:
        logger.info('Failed to pr_process records from different reference chromosomes.')
        exit()

    reference = references[0]
    oMin = 0
    oMax = 195471971
    #oMax = ref_lengths[reference]
    nMin = nMin
    nMax = nMax

    rescaled_coordinate_map = [(i,rescale(j, oMin, oMax, nMin, nMax)) for (i,j) in coordinate_map]

    return rescaled_coordinate_map
'''


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to umiViews module.')
    logger.info('This module provides a library of generators for viewing filtered presentations of a bam file.')