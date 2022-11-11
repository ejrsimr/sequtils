#!/usr/bin/env python3
# ejr: 2022-11-10
# lmd: 2022-11-10
# generate histogram of sequences lengths from FASTA file

from sequtils import *
import argparse 
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sb


def main():
    args = get_args()
    fasta = read_fasta(args.file)
    length_histogram(fasta, 100)


def length_histogram(fasta, num_bins):
    
    df = pd.DataFrame(list(fasta.items()), columns = ['header','seq'])  
    df['length'] = df['seq'].apply(len)

    hist = sb.histplot(data = df, x='length', bins=num_bins)
    fig = hist.get_figure()
    fig.savefig('length_histogram.png')  



###############################################################################
# Get command-line arguments using argparse
###############################################################################
def get_args():
    parser = argparse.ArgumentParser(description="Create Histogram of Sequence Lengths from FASTA file")
    # file defaults to stdin 
    parser.add_argument('--file', type = argparse.FileType('r'), default = sys.stdin, help = 'Input FASTA - defaults to STDIN')
    args = parser.parse_args()

    return args

if __name__ == "__main__":
    main()