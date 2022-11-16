#!/usr/bin/env python3
# ejr: 2022-11-10
# lmd: 2022-11-10
# collection of functions for FASTA and sequence manipulation

import os
import sys
import re
from re import finditer

def main():
    print("It Ran")

###############################################################################
# Read FASTA filehandle into dictionary
###############################################################################
def read_fasta(fh):
    header = ""
    fasta = {}

    for line in fh:
        line = line.rstrip()
        if (line.startswith(">")):
            header = line[1:]
            fasta[header] = []
        else:
            fasta[header].append(line)
    for header in fasta:
        fasta[header] = ''.join(fasta[header])

    return fasta

###############################################################################
# Read seqtable filehandle into dictionary - seqtable = header\tseq
###############################################################################
def read_seqtable(fh):
    header = ""
    fasta = {}

    for line in fh:
        line = line.rstrip()
        fields = line.split("\t")
        fasta[fields[0]] = fields[1]

    return fasta

###############################################################################
# Add newlines every 80 characters for FASTA formatting
###############################################################################
def insert_newlines(string, every=80):
    lines = []

    for i in range(0, len(string), every):
        lines.append(string[i:i+every])

    return '\n'.join(lines)

###############################################################################
# Reverse Complement Sequence
###############################################################################
def reverse_complement(seq):
    # complement sequence
    bases = str.maketrans('AGCTagct','TCGAtcga')
    # reverse sequences; return
    return seq.translate(bases)[::-1]

###############################################################################
# GREP FASTA for REGEX - save into BED dictionary
###############################################################################
def grep_regex(fasta, regex):
    bed = {}
    for header in fasta:
        seqlen = len(fasta[header])
        # forward matches
        for match in finditer(regex, fasta[header]):
            name = header + ":" + start + "-" + end
            bed[name] = {}
            bed[name]['start'] = str(match.start())
            bed[name]['end'] = str(match.end())
            bed[name]['chr'] = header
            bed[name]['strand'] = "+"
            bed[name]['score'] = "0"

        # reverse complement matches
        rc = reverse_complement(fasta[header])
        for match in finditer(regex, rc):
            name = header + ":" + start + "-" + end
            bed[name] = {}
            bed[name]['start'] = str(seqlen - match.end())
            bed[name]['end'] = str(seqlen - match.start())
            bed[name]['chr'] = header
            bed[name]['strand'] = "-"
            bed[name]['score'] = "0"

    return(bed)

###############################################################################
# Print FASTA Dictionary as FASTA to STDOUT
###############################################################################
def print_fasta(fasta):
    for header in fasta:
        print(">", header, sep="")
        print(insert_newlines(fasta[header]))

###############################################################################
# Print FASTA Dictionary as seqtable to STDOUT
###############################################################################
def print_seqtable(fasta):
    for header in fasta:
        print(header, fasta[header], sep="\t")

###############################################################################
# Print BED Dictionary as BED to STDOUT
###############################################################################
def print_bed(bed): 
    for name in bed:
        print("\t".join([bed[name]['chr'], bed[name]['start'], bed[name]['end'], name, bed[name]['score'], bed[name]['strand']]))


###############################################################################
# Run MAIN
###############################################################################
if __name__ == "__main__":
    main()
