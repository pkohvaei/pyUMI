#!/usr/bin/env python

import logging

from util import initialize_iterator

logger = logging.getLogger(__name__)


def tag_based_generator(pysam_iterator, tag=None, values=None, has_values=True):

    reads = initialize_iterator(pysam_iterator)

    for r in reads:
        if r.has_tag(tag):
            if r.get_tag(tag) in values:
                if has_values:
                    yield r
            else:
                if not has_values:
                    yield r


def field_based_generator(pysam_iterator, field=None, values=None, has_values=True):

    reads = initialize_iterator(pysam_iterator)

    for r in reads:
        value = getattr(r, field)
        if value in values:
            if has_values:
                yield r
        else:
            if not has_values:
                yield r


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to umiViews module.')
    logger.info('This module provides a library of generators for viewing filtered presentations of a bam file.')