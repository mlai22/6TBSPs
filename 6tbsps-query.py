# -*- coding: utf-8 -*-
'''
To query DNA sequences against a protein database.
Input: one or more multi-fasta DNA/RNA sequences
Output: 2D local alignment graph summary
Add multiprocess support

Author:
	Yuchen Ge, Zitong He

E-mail:
	yge15@jhmi.edu
	hezt@jhu.edu

Usage:
	$ python 6tbsps-query [-h] --db DB -o O [-t [T]] [--sm [SM]] reads.fa [reads.fa ...]

'''
#%%
import os
import pickle
import argparse
import multiprocessing
from itertools import product
import time

# custom scripts
import scripts.file_io as fio
import scripts.six_frame_translation as sft
from scripts.seed_and_extend import naive_seed_and_extend
from scripts.local_alignment_affine import LocalAlignment
from scripts.score_matrix import score_matrix, e_value_cal
#%%
def main():
	s_time = time.time()
	'''Main function of the 6tbsps-query'''
	parser = argparse.ArgumentParser(prog = '6tbsps-query', \
		description = 'Search DNA sequences against a pre-indexed protein database.')
	parser.add_argument('--db', required = True, \
		help='database base name of k-mer indices')
	parser.add_argument('-o', required = True, help = 'output directory')
	parser.add_argument('-t', default=1, nargs='?', type=int, help = 'number of threads')
	parser.add_argument('--sm', '--score-matrix', default='BLOSUM62', nargs='?', type=str, \
		help = 'scoring matrix: BLOSUM45, BLOSUM62 (default), BLOSUM80')
	parser.add_argument('reads', metavar = 'reads.fa', nargs = '+', \
		help = 'DNA reads in fasta format')

	args = parser.parse_args()
	sm = score_matrix(args.sm) # default BLOSUM62
	out_dir = args.o
	num_threads = args.t
	in_files = args.reads

	# get the pre-indexed protein database and sequences
	prot_db = fio.read_dict(args.db+'.kmer')
	prot_seq = fio.read_dict(args.db+'.prot')
	# get the constant k
	k = len(list(prot_db.keys())[0])

	# calculate total length of protein database
	n = sum([len(prot_seq[name]) for name in prot_seq])

	# get all reads
	reads = {}
	for name in in_files:
		fio.parse_fasta(name, reads)
	query_args = [[read_id, seq, out_dir, k, prot_db, sm, prot_seq, n] \
		for read_id, seq in reads.items()]
	with multiprocessing.Pool(processes=num_threads) as pool:
		pool.starmap(query, query_args)	
	print('running time: {}'.format(time.time() - s_time))

	return
#%%
def query(read_id, seq, out_dir, k, prot_db, sm, prot_seq, n):
	# query
	# make directory
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	frames = sft.six_frames(seq)
	with open(os.path.join(out_dir, read_id.replace('/', '|')+'.out'), 'w') as out_file, \
		open(os.path.join(out_dir, read_id.replace('/', '|')+'.summary'), 'w') as sum_file:
		output = [] # store output
		# for each of the 6 frames:
		for f in [-3, -2, -1, 1, 2, 3]:
			frame = frames[f]
			# 6-frame translation
			query = sft.translation(sft.transcription(frame))
			# calculate query length
			m = len(query)
			# seed and extend
			regions = naive_seed_and_extend(query, prot_db, k)
			# local alignment
			for ref_id, (s, e) in regions:
				subject = prot_seq[ref_id][s:e]
				la = LocalAlignment(query, subject, sm)
				S = la.fill_matrix()
				evalue = e_value_cal(m, n, S)
				
				output.append([f, read_id, query, ref_id, prot_seq, s, la, S, evalue])
		
		# sort by evalue, then by raw score
		output = sorted(output, key=lambda x: (x[-1], x[-2]))
		fio.align_out(output, out_file, sum_file)
	
	return
#%%
if __name__ == "__main__":
	main()
