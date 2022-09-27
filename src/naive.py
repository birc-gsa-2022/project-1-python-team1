"""Implementation of the naive exact matching algorithm."""

import argparse
from utils import *

def main():
    argparser = argparse.ArgumentParser(
        description="Exact matching using the naive method")
    argparser.add_argument("genome", type=argparse.FileType('r'))
    argparser.add_argument("reads", type=argparse.FileType('r'))
    args = argparser.parse_args()
    print(f"Find every reads in {args.reads.name} " +
          f"in genome {args.genome.name}")


if __name__ == '__main__':
    main()


def naive_idx(p, x):
    # Naive algortihm to find the start index of matching.
    n = len(x)
    m = len(p)
    idx = []
    for i in range(n - m + 1):
        j = 0
        while j < m:
            if x[i + j] == p[j]:
                j += 1
            else:
                break
        else:
            idx.append(i)


    return idx

def naive(fasta, fastq):
    #Naive algortihm that outputs simple SAM format
    ps, ps_names = parse_fastq(fastq)
    xs, xs_names = fasta_recs(fasta)

    for i in range(len(ps)):
        for j in range(len(xs)):
            idx = naive_idx(ps[i], xs[j])
            for k in idx:
                e = get_edits(ps[i], xs[j][k :(k + len(ps[i])-1)])[2]
                l = edits_to_cigar(e)
            
                print(ps_names[i] + '\t' + xs_names[j] + '\t' + str(k + 1) + '\t' + l + '\t' + ps[i])
    return ''