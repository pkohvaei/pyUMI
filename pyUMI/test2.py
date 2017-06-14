#!/usr/bin/env python

import pysam
import pandas as pd
import pickle

if __name__ == "__main__":

    sample = 'ATATAG'
    names = ['HISEQ:280:C9J9KANXX:1:2215:15146:97826', 'HISEQ:280:C9J9KANXX:1:2213:10613:32557', 'HISEQ:280:C9J9KANXX:1:2201:10903:83526']
    in_file = '/data/UMI/data/MUS/sample_' + sample + '.bam'
    st = pysam.AlignmentFile(in_file,"rb")

    st.reset()
    reads = st.fetch(until_eof=True)

    records = []
    for read in reads:
        if read.query_name in names:
            records.append((read.query_name, read.query_sequence, str(read.qual)))
    rdf = pd.DataFrame(records, columns=['N', 'S', 'Q'])

    with open('/home/parastou/mm10/sample.fastq', 'w+') as f:
        for index, row in rdf.iterrows():
            n, s, q = row
            f.write('@%s\n%s\n+\n%s\n'  % (n, s, q))

    print 'Success.'
