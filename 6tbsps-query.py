# -*- coding: utf-8 -*-
'''
To build a compressed hashtable for protein databases.
Input: one or more multi-fasta protein database(s)
Output: 2D local alignment graph 

Author:
	Yuchen Ge

E-mail:
	yge15@jhmi.edu

Usage:
	$ python 6tbsps-query --db $database -o $out_dir reads.fa

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
from scripts.seed_and_extend import naive_seed_and_extend
from scripts.local_alignment_affine import LocalAlignment
#%%
def main():
	'''Main function of the 6tbsps-query'''
	parser = argparse.ArgumentParser(prog = '6tbsps-query', \
		description = '')
	parser.add_argument('--db', required = True, \
		help='database base name of k-mer indices')
	parser.add_argument('-o', required = True, help = 'output directory')
	parser.add_argument('reads', metavar = 'reads.fa', nargs = '+', \
		help = 'DNA reads in fasta format')
	
	args = parser.parse_args()
	out_dir = args.o
	in_files = args.reads

	# get the pre-indexed protein database and sequences
	prot_db = fio.read_dict(args.db+'.kmer')
	prot_seq = fio.read_dict(args.db+'.prot')
	# get the constant k
	k = len(list(prot_db.keys())[0])

	# get all reads
	reads = {}
	for name in in_files:
		fio.parse_fasta(name, reads)

	# query
	for read_id, seq in reads.items():
		# make directory
		if not os.path.exists(out_dir):
			os.makedirs(out_dir)

		frames = sft.six_frames(seq)
		with open(os.path.join(out_dir, read_id.replace('/', '|')+'.out'), 'w') as out_file:
			# for each of the 6 frames:
			for f in [-3, -2, -1, 1, 2, 3]:
				frame = frames[f]
				# 6-frame translation
				query = sft.translation(sft.transcription(frame))
				# seed and extend
				regions = naive_seed_and_extend(query, prot_db, k)
				# local alignment
				for ref_id, (s, e) in regions:
					subject = prot_seq[ref_id][s:e]
					la = LocalAlignment(query, subject, None)
					la.fill_matrix()
					la.traceback()
					
					print('Frame:', str(f), file=out_file)
					print('Query:', query, file=out_file)
					print('Length:', len(query), file=out_file)
					print('Subject:', ref_id, file=out_file)
					print('Length:', len(prot_seq[ref_id]), file=out_file)
					print(file=out_file)

					la.display(out_file)
					print(file=out_file)

	return 
#%%
if __name__ == "__main__":
	main()
