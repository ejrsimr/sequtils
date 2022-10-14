#!/usr/bin/env python3
# ejr: 2022-10-13
# INPUT = FASTA
# for each sequence compare ends and calculate identity
# we are going to assume gaps are not allowed
# OUTPUT = sequence name with inverted repeat score
import sys
import signal
import string

# MAIN
def main():
    fasta = {}
    filename = sys.argv[1]

    with open(filename) as my_file:
        fasta = read_fasta(my_file)

    ends = get_ends(fasta, 20)
    scores = score_ends(ends)

    print('\t'.join(("name", "left", "right", "score")))
    for header in scores:
        print('\t'.join((header, scores[header]['left'], scores[header]['right'],str(scores[header]['score']))))


# READ IN FASTA
def read_fasta(filename):
    header = ""
    fasta = {}

    for line in filename:
        line = line.rstrip()
        if (line[0] == ">"):
            header = line[1:]
            fasta[header] = []
        else:
            fasta[header].append(line)

    for header in fasta:
        fasta[header] = ''.join(fasta[header])

    return fasta

# reverse complement sequences
def rc(seq):

    complements = seq.maketrans('ATGCatgc', 'TACGtacg')
    rcseq = seq.translate(complements)[::-1]
    return rcseq


# count number of identical seqs
def count_identical_bp(seq1, seq2):
    score = 0
 
    # Traverse the string 1 char by char
    for i in range(len(seq1)) :
        if seq1[i] == seq2[i] :
            score += 1;

    return(score)

def score_ends(ends):

    for header in ends:
        left = ends[header]['left']
        right = ends[header]['right']
        rc_right = rc(right)

        score = count_identical_bp(left, rc_right)

        ends[header]['score'] = score

    return(ends)

# extract ends from fasta file
def get_ends(fasta, size):

    ends = {}
    for header in fasta:
        ends[header] = {}
        seq = fasta[header]
        ends[header]['left'] = seq[0:size]
        ends[header]['right'] = seq[(size) * -1:] 

    return(ends)

# MAIN
if __name__ == "__main__":

    # this catches sigpipe errors so you don't get an error message if you tail of head output
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    main()