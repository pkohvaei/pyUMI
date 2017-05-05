#!/usr/bin/env python

import logging
from collections import Counter, Iterable
from itertools import chain

from util import initialize_iterator


logger = logging.getLogger(__name__)


'''Functions for generating statistical views on multi read and unique hash tables'''


def gene_tagged_multi_maps(multi_read_table):

    reads = len(multi_read_table)
    alignments = 0
    gene_tagged = 0
    all_gene_tagged = 0
    only_ref = 0

    for multi_read in multi_read_table:
        mappings = multi_read_table[multi_read][2]
        a = len(mappings)
        alignments += a
        g = len([gene for gene, ref, start in mappings if gene is not None])
        n = len([gene for gene, ref, start in mappings if gene is None])

        if g:
            gene_tagged += a

        if a == g:
            all_gene_tagged += g

        if g == 0:
            only_ref += n

    return reads, alignments, gene_tagged, all_gene_tagged, only_ref


def update_uniq_hashtable(uniq_hashtable, xc, xm, qname):
    if xm in uniq_hashtable[xc]:
        uniq_hashtable[xc][xm].append(qname)
    else:
        uniq_hashtable[xc].update({xm:[qname]})


def update_uniq_map(uniq_map, qname, locus):

    if qname in uniq_map:
        logger.info('Item already marked as uniquely mapped read!')
        quit()
    else:
        uniq_map.update({qname:locus})


def update_multimap(multi_map, qname, xc, xm, ge, qref, refpos):
    if not(qname in multi_map):
        multi_map.update({qname:(xc,xm,[(ge,qref,refpos)])})
    else:
        multi_map[qname][2].append((ge,qref,refpos))


def find_barcode_matches(uniq_map, xc, xm):
    if not (xc in uniq_map):
        print "key error!"
        return
    else:
        if xm in uniq_map[xc]:
            return uniq_map[xc][xm]
    return []


def build_multimapping_hashtable(pysam_iterator, cell_barcode=None, region='all'):

    reads = initialize_iterator(pysam_iterator)

    multi_maps = {}

    for read in reads:
        if not read.is_unmapped:

            qname = read.query_name
            xc = read.get_tag("XC")
            xm = read.get_tag("XM")
            refname = read.reference_name
            refpos = read.reference_start
            gene = None if not read.has_tag("GE") else read.get_tag("GE")

            if not read.get_tag("NH") == 1:
                if not cell_barcode or cell_barcode == xc:
                    if region == 'all' or (region == 'gene' and gene) or (region == 'ref' and not gene):
                        update_multimap(multi_maps, qname, xc, xm, gene, refname, refpos)

    return multi_maps


def resolve_multimap(multimap_hashtable, qname, gene):

    xc, xm, genes = multimap_hashtable[qname]
    rescued = [(i,j,k) for (i,j,k) in genes if i == gene]
    if len(rescued) == 1:
        multimap_hashtable.pop(qname, None)
        return (qname,xc,xm,rescued[0])
    else:
        return False


def update_resolved_multimaps(multimap_hashtable, uniq_hashtable, uniq_map, resolved_multimaps):

    uniques_updates = []
    for qname, gene in resolved_multimaps:
        result = resolve_multimap(multimap_hashtable, qname, gene)
        if result:
            qname, xc, xm, locus = result
            uniques_updates.append(result)
            update_uniq_hashtable(uniq_hashtable, xc, xm, qname)
            update_uniq_map(uniq_map,qname,locus)

    return uniques_updates


def build_uniques_associations(pysam_iterator, cell_barcodes=None):

    reads = initialize_iterator(pysam_iterator)

    uniques_hashtable = {}
    uniques = {}

    if not cell_barcodes or not isinstance(cell_barcodes,Iterable):
        logger.info('Cell barcode iterable not present! Quitting the program...')
        exit()

    for cb in cell_barcodes:
        uniques_hashtable.update({cb: {}})

    for read in reads:
        if not read.is_unmapped:
            if read.get_tag("NH") == 1:
                xc = read.get_tag("XC")
                if xc in cell_barcodes:
                    qname = read.query_name
                    xm = read.get_tag("XM")
                    refname = read.reference_name
                    refpos = read.reference_start
                    gene = None if not (read.has_tag("GE")) else read.get_tag("GE")
                    update_uniq_hashtable(uniques_hashtable, xc, xm, qname)
                    uniques.update({qname: (gene, refname, refpos)})

    return uniques_hashtable, uniques


def get_uniquesmap_genes(uniques, gene=True):

    switch = 0
    if not(gene):
        switch = 1

    result = [uniques[k][switch] for k in uniques if uniques[k][switch]]

    return result


def get_uniques_qname(uniques):

    q_names = [n for n in uniques]

    return q_names


def query_uniquesmap(uniques, keyword=None):

    result = []

    if keyword == 'gene':
        result = get_uniquesmap_genes(uniques)
    elif keyword == 'region':
        result = get_uniquesmap_genes(uniques, gene=False)
    elif keyword == 'qname':
        get_uniques_qname(uniques)
    else:
        logger.info("Keyword must be one of: 'qname', 'gene', 'region'")

    return result


def query_uniques_hashtable(uniques_hashtable, keyword='None'):

    result = []

    if keyword == 'umi':
        result = uniques_hashtable.items()
    else:
        logger.info("Keyword not known.")

    return result


def get_bamfile_genes(alignmentFile):

    alignmentFile.reset()
    reads = alignmentFile.fetch(until_eof=True)

    genes = []

    for read in reads:
        if read.has_tag("GE"):
            genes.append(read.get_tag("GE"))

    return genes


def get_multimap_genes(multimap_hashtable, gene=True):

    result = []

    switch = 0
    if not(gene):
        switch = 1

    r1 = [i[1][2] for i in multimap_hashtable.items()]
    r2 = list(chain.from_iterable(r1))
    result = [k[switch] for k in r2 if k[switch]]

    return result


def get_multimap_umis(multimap_hashtable):

    result = [i[1][1] for i in multimap_hashtable.items()]

    return result


def query_multimap_hashtable(multimap_hashtable, keyword=None):

    result = []

    if keyword == 'gene':
        result = get_multimap_genes(multimap_hashtable)
    elif keyword == 'region':
        result = get_multimap_genes(multimap_hashtable, gene=False)
    elif keyword == 'umi':
        result = get_multimap_umis(multimap_hashtable)
    elif keyword == 'qname':
        result = multimap_hashtable.keys()
    else:
        logger.info("Keyword must be one of: 'qname', 'umi', 'gene', 'region'")

    return result


if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logger.info('Call to multimap resolution toolkit.')