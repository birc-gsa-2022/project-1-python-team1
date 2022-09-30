import time
import matplotlib.pyplot as plt
from naive import *



def time_naive(fasta, fastq):
    n = []
    t = []
    ps, ps_names = parse_fastq(fastq)
    xs, xs_names = fasta_recs(fasta)
    for i in range(len(ps)):
        for j in range(len(xs)):
            t0 = time.time()
            idx = naive_idx(ps[i], xs[j])
            t1 = time.time()
            total = t1 - t0
            n.append(len(xs[i]))
            t.append(total)
    return n, t


result = time_naive("/home/mathilde/Documents/Kandidat/GSA/Failure_files/lin/__TEST__/tools/lin/genome-100-10.fa", "/home/mathilde/Documents/Kandidat/GSA/Failure_files/lin/__TEST__/tools/lin/genome-100-10.fa-reads-10-10-0.fq")
plt.plot(result[0],result[1], 'ro')
plt.show()


