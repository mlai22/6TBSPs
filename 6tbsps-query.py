# -*- coding: utf-8 -*-
'''
To build a compressed hashtable for protein databases.
Input: one or more multi-fasta protein database(s)
Output: compressed dictionary (hashtable) of kmers and reference proteins

Author:
	Yuchen Ge

E-mail:
	yge15@jhmi.edu

Usage:
	$ python 6tbsps.py-build [-k $kmer_len] -o $out_basename protein.faa

Attributes:
	Parse an input FASTA and translate the DNA sequences in 6-frame, 
	output a single FASTA file of compressed protein sequences.
'''
#%%
import os
import pickle
import argparse

# custom scripts
import scripts.file_io as fio
import scripts.six_frame_translation as sft
from scripts.local_alignment import LocalAlignment
#%%
def main():
	'''Main function of the 6tbsps-query'''
	parser = argparse.ArgumentParser(prog = '6tbsps-query', \
		description = '')
	# parser.add_argument('--db', help='database base name of k-mer indices')
	parser.add_argument('-o', help = 'output directory', required = True)
	parser.add_argument('reads', metavar = 'reads.fa', nargs = '+', \
		help = 'DNA reads in fasta format')
	
	args = parser.parse_args()
	out_dir = args.o
	in_files = args.reads

	# read in all reads
	reads = {}
	for name in in_files:
		fio.parse_fasta(name, reads)

	# query
	for read_id, seq in reads.items():
		# make directory
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)

		frames = sft.six_frames(seq)
		# for each of the 6 frames:
		for f in [-3, -2, -1, 1, 2, 3]:
			with open(os.path.join(out_dir, read_id+'_'+str(f)+'.out'), 'w') as out_file:
				frame = frames[f]
				peptide = sft.translation(sft.transcription(frame))
				print('Frame:', str(f), file=out_file)
				print('Query:', peptide, file=out_file)
				print('Length:', len(peptide), file=out_file)
	
	return 
#%%
if __name__ == "__main__":
	main()
