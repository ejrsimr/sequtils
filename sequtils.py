#!/usr/bin/env python3
# ejr: 2022-11-10
# lmd: 2022-11-16
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
        # starts with handles blank lines better than line[0]
        if (line.startswith(">")):
            header = line[1:]
            fasta[header] = []
        else:
            fasta[header].append(line)
    # append is more efficient that string concatenation
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
# Read file of sequence names or partial names and return list
###############################################################################
def file_to_list(fh):
    out_list = []

    for line in fh:
        line = line.rstrip()
        out_list.append(line)

    return out_list

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
# Filter FASTA File by length and gc content
###############################################################################
def filter_fasta(fasta, min_length, max_length, min_gc, max_gc):
    fasta_out = {}

    for header in fasta:
        s_gc = calc_gc(fasta[header])
        s_len = len(fasta[header])
        if s_gc >= min_gc & s_gc <= max_gc & s_len >= min_length & s_len <= max_length:
            fasta_out[header] = fasta[header]
    
    return fasta_out

###############################################################################
# calculate GC content of sequence - returns percentage in decimal (e.g. 0.21)
###############################################################################    
def calc_gc(seq):
    useq = seq.upper()
    g = useq.count('G')
    c = useq.count('C')
    a = useq.count('A')
    t = useq.count('T')
    perc_gc = (g + c) / (g + c + a + t)

    return(perc_gc)

###############################################################################
# Subset FASTA - require full header (fast) - returns FASTA dictionary
###############################################################################    
def subset_fasta_is(fasta, names_list, exclude):
    # this is faster than pattern matching, but on big FASTA files, this
    # will be a lot slower than indexed retrieval
    fasta_out = {}

    if exclude == True:
        for header in fasta:
            if header not in names_list: 
                fasta_out[header] = fasta[header]

    else:
        for header in fasta:
            if header in names_list:
                fasta_out[header] = fasta[header]

    return(fasta_out)

###############################################################################
# Subset FASTA - only startswith (slow) - returns FASTA dictionary
###############################################################################    
def subset_fasta_startswith(fasta, names_list, exclude):
    fasta_out = {}

    # for exclusions we need a full list and not just partials
    match_list = []
    for header in fasta:
        for pattern in names_list:
            if header.startswith(pattern):
                match_list.append(header)

    if exclude == True:
        for header in fasta:
            if header not in match_list:
                fasta_out[header] = fasta[header]
    else:
        for header in match_list:
            fasta_out[header] = fasta[header]

    return(fasta_out)

###############################################################################
# Subset FASTA - contains (very slow) - returns FASTA dictionary
###############################################################################    
def subset_fasta_contains(fasta, names_list, exclude):
    fasta_out = {}

    # for exclusions we need a full list and not just partials
    match_list = []
    for header in fasta:
        for pattern in names_list:
            if pattern in header:
                match_list.append(header)

    if exclude == True:
        for header in fasta:
            if header not in match_list:
                fasta_out[header] = fasta[header]
    else:
        for header in match_list:
            fasta_out[header] = fasta[header]

    return(fasta_out)


###############################################################################
# Run MAIN
###############################################################################
if __name__ == "__main__":
    main()
