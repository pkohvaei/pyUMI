import os
import sys

import pysam
from Bio import SeqIO, Seq, SeqRecord
import StringIO
import glob


def bam_to_umi(in_file):

    """Generator to extract UMI barcodes from bam records.
    """
    bam_file = pysam.Samfile(in_file, "rb")
    for read in bam_file:
        xm = read.get_tag('XM')

        yield xm


def bam_to_rec(in_file):
    """Generator to convert BAM files into Biopython SeqRecords.
    """
    bam_file = pysam.Samfile(in_file, "rb")
    for read in bam_file:
        seq = Seq.Seq(read.seq)
        q = read.qual
        n = read.query_name

        fastq_string = "@%s\n%s\n+\n%s\n" % (n, seq, q)
        record = SeqIO.read(StringIO.StringIO(fastq_string), "fastq-sanger")

        yield record


def main(in_file):

    out_file = "%s.fastq" % os.path.splitext(in_file)[0]
    umi_file = "%s.umi" % os.path.splitext(in_file)[0]

    with open(out_file, "w") as out_handle:
        # Write records from the BAM file one at a time to the output file.
        # Works lazily as BAM sequences are read so will handle large files.
        SeqIO.write(bam_to_rec(in_file), out_handle, "fastq-sanger")

    with open(umi_file, "w") as out_handle:

        for umi in bam_to_umi(in_file):
            out_handle.write(umi + '\n')

    '''
        if read.is_reverse:
            seq = seq.reverse_complement()
        rec = SeqRecord.SeqRecord(seq, read.qname, "", "")
    '''


def test_output(in_file_bam, in_file_fq, in_file_umi):

    l1 = []
    l2 = []
    l = []

    for record in SeqIO.parse(in_file_fq, "fastq-sanger"):
        l1.append(record.id)

    with open(in_file_umi, 'r+') as fu:
        for line in fu:
            l2.append(line.rstrip())

    st = pysam.AlignmentFile(in_file_bam, "rb")
    st.reset()
    reads = st.fetch(until_eof=True)
    for r in reads:
        n = r.query_name
        xm = r.get_tag('XM')
        l.append((n, xm))

    if zip(l1,l2) == l:
        return 'Equal!'
    else:
        return 'Not equal!'


if __name__ == "__main__":

    '''
    os.chdir("/data/parastou/UMI/data/HG/AMLPDXdemux/")

    for in_file in glob.glob("*.bam"):

        # main(in_file)
        umi_file = "%s.umi" % os.path.splitext(in_file)[0]
        fq_file = "%s.fastq" % os.path.splitext(in_file)[0]

        print in_file, test_output(in_file, fq_file, umi_file)
    '''
    
    in_file = '/data/parastou/UMI/data/HG/AMLPDXdemux/HGsample_AATGTA.bam'
    main(in_file)	




