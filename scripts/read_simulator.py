# -*- coding: utf-8 -*-
"""Genomic read simulator.

This module simulates error-free genomic reads of arbitary length and coverage.

Author: 
    Yuchen (Peter) Ge

Email: 
    yge15@jhmi.edu
    
Usage:
    $ python read_simulator.py -l $readLength -c $coverage ref.fa ...

Attributes:
    parse_fasta(str):  parse one or more input FASTA reference(s)
    
    write_fasta(list): print a FASTA-like input to stdout
    
    read_simulation(str, int, int): simulate error-free reads by scanning 
        through the geome


"""
import sys
import argparse
#%%
def parse_fasta(fn, seqs):
    """Parse referece sequences from a FASTA file name. 
        For each entry, return a dictionary that maps its sequence identifier 
        to the nuleotide sequeneces.
    
    Args:
        fn (str):       file name of the FASTA file
        seqs (dict):    dictionary of sequences identifiers to nucleitode 
                        sequences

    Returns:
        None
    
    """ 
    with open(fn) as fh:
        while True:
            line = fh.readline()
            if len(line) == 0:
                break  # end of file
            if line.startswith('>'):
                name = line[1:].rstrip()
                seq = ''
            else:
                seq += line.rstrip()
            
            seqs[name] = seq
            
    return
#%%
def read_simulation(seqs, readLen, cov):
    """Simulate error-free reads from referece sequences, given a specific read
        length and coverage. 
        Return a list of tuples of read identifier with range in the original 
        reference and read sequence, in the alphabetically order of the 
        original sequence identifiers.
    
    Args:
        seqs (dict):        dictionary of sequences identifiers to nucleitode 
                            sequences
        readLen (int):      read length
        cov (float):        coverage, i.e. read depth

    Returns:
        reads (list):       list of tuples of (readID, readSeq) for each
                            error-free simulated read
    
    """ 
    reads = []
    for name in sorted(seqs):
        seq = seqs[name]
        for i in range(0, len(seq) - readLen + 1, int(readLen/cov)):
            readID = name + '|[%d,%d)' % (i, i + readLen)
            readSeq = seq[i: i + readLen]
            reads.append((readID, readSeq))
    
    return reads     
#%%
def write_fasta(reads):
    """Write to standard output the simulated reads.
    
    Args:
        reads (list):   list of tuples of (readID, readSeq) for each
                        error-free simulated read
    
    Returns:
        None
    
    """ 
    for readID, readSeq in reads:
        readID  =  '>' + readID
        sys.stdout.write(readID + '\n')
        sys.stdout.write(readSeq + '\n')
        
    return
        
#%%
def write_fastq(reads):
    pass
#%%
def main():
    """Main driver function of the program.
    
    """
    parser = argparse.ArgumentParser(prog='read_simulator', \
                                     description= 'This module simulates '
                                     'error-free genomic reads of arbitary '
                                     'length and coverage.')
    parser.add_argument('-l', nargs=1, help='read length')
    parser.add_argument('-c', nargs=1, help='coverage')
    parser.add_argument('refs', metavar='ref.fa', nargs='+', 
                        help='FASTA references')
    
    args = parser.parse_args()
    readLen = int(args.l[0])
    cov = float(args.c[0])
    refSeqs_fn = args.refs
    
    refSeqs = {}
    for name in refSeqs_fn:
        parse_fasta(name, refSeqs)

    reads = read_simulation(refSeqs, readLen, cov)
    write_fasta(reads)
    
    return    
#%%
if __name__ == "__main__": 
    main()
