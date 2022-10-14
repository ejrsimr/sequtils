#!/usr/bin/env python3
# ejr: 2022-10-13
# INPUT = FASTA
# OUTPUT = bed file of inverted termianl repeats.
# name = repeat
# score = num matches between pair
import sys
import signal
import string
window_size = 4000
min_itr = 16
min_score = 14

# MAIN
def main():
    fasta = {}
    filename = sys.argv[1]
    ## Open FASTA file
    with open(filename) as my_file:
        fasta = read_fasta(my_file)

    ## for each seq search for inverted repeat <= 4k upstream
    ## if there is a match return seq, score, position
    for header in fasta:
        seq = fasta[header]
        for i in range(len(seq)):
            # revcom the left end of window (cheaper than rc on window)
            left_rc = rc(seq[i:i+min_itr])
            # for the 4k window downstream of location
            for ii in range(i+min_itr,i+window_size-min_itr):
                score = get_score(left_rc, seq[ii:ii+min_itr])
                if score > min_score:
                    print("\t".join((header, str(i), str(ii+min_itr), seq[i:i+min_itr], str(score), ".")))


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
def get_score(seq1, seq2):
    score = 0
 
    print(seq1, seq2)
    # Traverse the string 1 char by char
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            score += 1;

    return(score)


# MAIN
if __name__ == "__main__":

    # this catches sigpipe errors so you don't get an error message if you tail of head output
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    main()