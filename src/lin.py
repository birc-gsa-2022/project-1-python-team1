"""Implementation of a linear time exact matching algorithm."""

import argparse

from utils import *

def borderArray(p, x):
    # Step 1: Construct s = p$x and its borderarray b
    s = p + "$" + x
    ba = np.zeros(len(s), dtype=int)
    ba[0] = 0
    for i in range(1, len(s)):
        b = ba[i - 1]
        while b > 0 and s[i] != s[b]:
            b = ba[b - 1]
        if s[i] == s[b]:
            ba[i] = b + 1
        else:
            ba[i] = 0
    return ba, s


def occurrence(p, x):
    # Step 2: Report an occurence of p in x at position i-2m if b[i]=m
    idx = []
    ba = borderArray(p, x)[0]
    s = borderArray(p, x)[1]
    for i in range(len(ba)):
        if ba[i] == len(p):
            idx.append(i - 2 * len(p))
    return idx

def lin(fasta, fastq):
    #Linear border array that outputs simple SAM format
    ps, ps_names = parse_fastq(fastq)
    xs, xs_names = fasta_recs(fasta)

    for i in range(len(ps)):
        for j in range(len(xs)):
            idx = occurrence(ps[i], xs[j])
            for k in idx:
                e = get_edits(ps[i], xs[j][k :(k + len(ps[i])-1)])[2]
                l = edits_to_cigar(e)
            
                print(ps_names[i] + '\t' + xs_names[j] + '\t' + str(k + 1) + '\t' + l + '\t' + ps[i])
    return ''

def main():
    argparser = argparse.ArgumentParser(
        description="Exact matching in linear time")
    argparser.add_argument("genome", type=argparse.FileType('r'))
    argparser.add_argument("reads", type=argparse.FileType('r'))
    args = argparser.parse_args()
    print(f"Find every reads in {args.reads.name} " +
          f"in genome {args.genome.name}")
    lin(args.genome, args.reads)


if __name__ == '__main__':
    main()
