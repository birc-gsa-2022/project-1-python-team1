import numpy as np

from Bio import SeqIO  # Bio.SeqIO is the sequence Input/Output interface for BioPython

# -----------------------------------------
# Project 1: Basic exact pattern matching
# -----------------------------------------


def naive_idx(p, x):
    # Naive algortihm to find the start index of matching.
    n = len(x)
    m = len(p)
    idx = []
    for i in range(n - m + 1):
        for j in range(m):
            if x[i + j] == p[j]:
                if j == m - 1:
                    idx.append(i)

    return idx


# Border arrays and exact pattern matching


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


def parse_fastq(file):
    # Function to return sequences and fastq names from fastq file:
    ps = []
    ps_names = []
    for line in file:
        if line.startswith("@"):
            newline = line.strip("@")
            newline1 = newline.strip()
            ps_names.append(newline1)
        if not line.startswith("@"):
            newline = str(line)
            newline1 = newline.strip()
            ps.append(newline1)
    return ps, ps_names



def fasta_recs(file):
    # Function to return sequences from fasta file:
    name = []
    sequence = []
    for record in SeqIO.parse(file, "fasta"):
        name.append(record.id)
        sequence.append(str(record.seq))

    return sequence, name





def get_edits(str1, str2):
    # Function that translates between two representations.
    # get_edits("aacgt-c", "a-attac") -> ("aacgtc", "aattac", "MDMMMIM")
    if len(str1) and len(str2) == 0:
        return ""

    elif not "-" in str1 and str2:
        return str1.replace("-", ""), str2.replace("-", ""), len(str1) * "M"

    else:
        str = ""
        for i in range(len(str1)):
            if str1[i] == "-" and str2[i] != "-":
                str = str + "I"

            elif str1[i] != "-" and str2[i] == "-":
                str = str + "D"

            else:
                str = str + "M"
        return str1.replace("-", ""), str2.replace("-", ""), str


def edits_to_cigar(e):
    # edits_to_cigar("MMIMMDDM") -> "2M1I2M2D1M"
    cigar = ""
    count = 0
    for i in range(len(e)):
        if e[i] == e[i - 1]:
            count += 1
            char = str(count) + e[i]
        else:
            count = 1
            cigar = cigar + char
            char = str(count) + e[i]
    cigar = cigar + char
    return cigar
