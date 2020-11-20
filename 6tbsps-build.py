# -*- coding: utf-8 -*-
'''
To build a compressed hashtable for protein databases.
Input: one or more multi-fasta protein database(s)
Output: compressed dictionary (hashtable) of kmers and reference proteins

Author:
	Mei-Yu Lai

E-mail:
	mlai22@jhu.edu

Usage:
	$ python 6tbsps.py-build protein.faa

Attributes:
	Parse an input FASTA and translate the DNA sequences in 6-frame, 
	output a single FASTA file of compressed protein sequences.
'''
#%%
import os
import argparse
import pickle
#%%
def parse_fasta(filename, ID_seq):
	'''
	Parse protein sequences from FASTA file

	Args:
        filename (str):   file name of the protein FASTA file
		ID_seq (dict):    a dictionary with keys reference IDs and values 
                            protein sequences

	Returns:
		None
	'''
	with open(filename, 'r') as fh:
		while True:
			line = fh.readline()	
			if len(line) == 0:
				break
			if line[0] == '>':
				refID = line[1:].rstrip()
				protein_seq = ''
			else:
				protein_seq += line.rstrip()
			ID_seq[refID] = protein_seq
			
	return
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
def write_dict(dictionary, dir_name, base_name):
	'''
	Write a dictionary to file, compressed using pickle

	Args:
		dictionary (dict): 		a python dictionary
		dir_name (str): 		output directory name
		base_name (str):		output base name

	Returns:
		None
	'''
	with open(os.path.join(dir_name, base_name)+'.pickle', 'wb') as handle:
		pickle.dump(dictionary, handle, protocol = pickle.HIGHEST_PROTOCOL)	
	
	return
#%%
def main():
	'''The main function for 6TBSPS-build'''

	parser = argparse.ArgumentParser(prog = '6tbsps-build', \
		description='Build a compressed hashtable for protein databases.')
	parser.add_argument('prot_faa', metavar='protein.faa', nargs = '+', \
		help='protein FASTA filename')
	parser.add_argument('-k', default=5, nargs='?', type=int, \
		help='k-mer length (default:5)')
	parser.add_argument('-o', help='base name of the k-mer indices')

	args = parser.parse_args()
	in_files = args.prot_faa
	out_path = os.path.join(os.getcwd(), args.o)
	out_dir = os.path.dirname(out_path)
	out_base = os.path.basename(out_path)
	k = int(args.k)

	prot_seqs = {}
	for name in in_files:
		parse_fasta(name, prot_seqs)
	write_dict(prot_seqs, out_dir, out_base+'.prot')

	prot_kmer = protein_kmer_table(prot_seqs, k)
	write_dict(prot_kmer, out_dir, out_base+'.kmer')

	return
#%%
if __name__ == '__main__':
	main()
