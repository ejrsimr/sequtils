#!/usr/bin/env python3
# ejr: 2022-11-10
# lmd: 2022-11-10
# generate histogram of sequences lengths from FASTA file

from sequtils import *
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sb

###############################################################################
### MAIN ######################################################################
###############################################################################
def main():
    fasta = read_fasta(sys.stdin)
    length_histogram(fasta, 100)

###############################################################################
### LENGTH HISTOGRAM ##########################################################
###############################################################################
def length_histogram(fasta, num_bins):
    
    df = pd.DataFrame(list(fasta.items()), columns = ['header','seq'])  
    df['length'] = df['seq'].apply(len)

    hist = sb.histplot(data = df, x='length', bins=num_bins)
    fig = hist.get_figure()
    fig.savefig('length_histogram.png')  
    
###############################################################################
### RUN MAIN ##################################################################
###############################################################################
if __name__ == "__main__":
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    main()