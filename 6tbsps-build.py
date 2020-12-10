# -*- coding: utf-8 -*-
'''
To build a compressed hashtable for protein databases.
Input: one or more multi-fasta protein database(s)
Output: compressed dictionary (hashtable) of kmers and reference proteins

Author:
	Mei-Yu Lai
	Yuchen Ge

E-mail:
	mlai22@jhu.edu
	yge15@jhmi.edu

Usage:
	$ python 6tbsps-build [-h] [-k [KMER]] --db DB protein.faa [protein.faa ...]
	
'''
#%%
import os
import argparse

# custom src
import src.file_io as fio
#%%
def protein_kmer_table(seqs, k):
	'''
	Create a dictionary of k-mer protein sequence

	Args:
		seqs (dict): 		protein dictionary
		k (int): 			an interger k for k-mer

	Returns:
		table (dict): 		a dictionary that mapts each k-mer to the list of 
							referece names and positions
	'''
	table = {}
	for refID, seq in seqs.items():
		for loc in range(len(seq) - k + 1):
			kmer = seq[loc:loc+k]
			if kmer not in table:
				table[kmer] = []
			table[kmer].append((refID, loc))

	return table
#%%
def main():
	'''The main function for 6TBSPS-build'''

	parser = argparse.ArgumentParser(prog = '6tbsps-build', \
		description='Build a compressed hashtable for protein databases.')
	parser.add_argument('prot_faa', metavar='protein.faa', nargs = '+', \
		help='protein FASTA filename')
	parser.add_argument('-k', '--kmer', default=3, nargs='?', type=int, \
		help='k-mer length (default:3)')
	parser.add_argument('--db', '--database', \
		help='database base name of k-mer indices', required = True)

	args = parser.parse_args()
	in_files = args.prot_faa
	out_dir = os.path.dirname(args.db)
	out_base = os.path.basename(args.db)
	k = int(args.kmer)

	prot_seqs = {}
	for name in in_files:
		fio.parse_fasta(name, prot_seqs)
	fio.write_dict(prot_seqs, out_dir, out_base+'.prot')

	prot_kmer = protein_kmer_table(prot_seqs, k)
	fio.write_dict(prot_kmer, out_dir, out_base+'.kmer')

	return
#%%
if __name__ == '__main__':
	main()
