from statistics import mean
import time
import matplotlib.pyplot as plt
from naive import *
from lin import *
import numpy as np




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
            n.append(len(xs[j]) * len(ps[i]))
            t.append(total)
    
    return [mean(n)], [mean(t)]

result1 = time_naive("test_data/fasta/genome-100-10.fa", "test_data/fastq/genome-100-10.fa-reads-10-10-0.fq")
result2 = time_naive("test_data/fasta/genome-100-10.fa", "test_data/fastq/genome-100-10.fa-reads-10-50-0.fq")
result3 = time_naive("test_data/fasta/genome-500-10.fa", "test_data/fastq/genome-500-10.fa-reads-10-10-0.fq")
result4 = time_naive("test_data/fasta/genome-500-10.fa", "test_data/fastq/genome-500-10.fa-reads-10-50-0.fq")
result5 = time_naive("test_data/fasta/genome-1000-10.fa", "test_data/fastq/genome-1000-10.fa-reads-10-10-0.fq")
result6 = time_naive("test_data/fasta/genome-1000-10.fa", "test_data/fastq/genome-1000-10.fa-reads-10-50-0.fq")
result7 = time_naive("test_data/fasta/worst-case.fa", "test_data/fastq/worst-case.fq")
result8 = time_naive("test_data/fasta/worst-case1.fa", "test_data/fastq/worst-case.fq")

x = result1[0] + result2[0] + result3[0] + result4[0] + result5[0] + result6[0] + result7[0] + result8[0]
y = result1[1] + result2[1] + result3[1] + result4[1] + result5[1] + result6[1] + result7[1] + result8[1]


def time_lin(fasta, fastq):
    n = []
    t = []
    ps, ps_names = parse_fastq(fastq)
    xs, xs_names = fasta_recs(fasta)
    for i in range(len(ps)):
        for j in range(len(xs)):
            t0 = time.time()
            idx = occurrence(ps[i], xs[j])
            t1 = time.time()
            total = t1 - t0
            n.append(len(xs[j]) * len(ps[i]))
            t.append(total)
    
    return [mean(n)], [mean(t)]

result12 = time_lin("test_data/fasta/genome-100-10.fa", "test_data/fastq/genome-100-10.fa-reads-10-10-0.fq")
result22 = time_lin("test_data/fasta/genome-100-10.fa", "test_data/fastq/genome-100-10.fa-reads-10-50-0.fq")
result32 = time_lin("test_data/fasta/genome-500-10.fa", "test_data/fastq/genome-500-10.fa-reads-10-10-0.fq")
result42 = time_lin("test_data/fasta/genome-500-10.fa", "test_data/fastq/genome-500-10.fa-reads-10-50-0.fq")
result52 = time_lin("test_data/fasta/genome-1000-10.fa", "test_data/fastq/genome-1000-10.fa-reads-10-10-0.fq")
result62 = time_lin("test_data/fasta/genome-1000-10.fa", "test_data/fastq/genome-1000-10.fa-reads-10-50-0.fq")
result72 = time_lin("test_data/fasta/worst-case.fa", "test_data/fastq/worst-case.fq")
result82 = time_lin("test_data/fasta/worst-case1.fa", "test_data/fastq/worst-case.fq")

y2 = result12[1] + result22[1] + result32[1] + result42[1] + result52[1] + result62[1] + result72[1] + result82[1]

fig, ax = plt.subplots()
ax.scatter(x, y, label = "naive")
ax.scatter(x, y2, label = "linear")
ax.set_title('Time complexity')
ax.legend()
ax.set_ylabel('Time')
ax.set_xlabel('n*m')
plt.savefig('/home/mathilde/Documents/Kandidat/GSA/Project/project1/project-1-python-team1/figs/time.png')