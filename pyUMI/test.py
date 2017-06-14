#!/usr/bin/env python

import pysam
import pandas as pd
import pickle

if __name__ == "__main__":

    sample = 'ATATAG'
    in_file = '/data/UMI/data/MUS/sample_' + sample + '.bam'
    st = pysam.AlignmentFile(in_file,"rb")

    st.reset()
    reads = st.fetch(until_eof=True)

    records = []
    for read in reads:
        records.append((read.query_name, read.query_sequence, str(read.qual)))
    rdf = pd.DataFrame(records, columns=['N', 'S', 'Q'])
    rdf2 = rdf.drop_duplicates(subset=['N'])
    pickle.dump(rdf2, open('rdf.pkl', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    with open('/home/parastou/mm10/downsample.fastq', 'w+') as f:
        for index, row in rdf2.iterrows():
            n, s, q = row
            f.write('@%s\n%s\n+\n%s\n'  % (n, s, q))

    print 'Success.'
