"""Implementation of a linear time exact matching algorithm."""

import argparse

from utils import *

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
            
                print(ps_names[i] + '   ' + xs_names[j] + '   ' + str(k + 1) + '   ' + l + '   ' + ps[i])
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
