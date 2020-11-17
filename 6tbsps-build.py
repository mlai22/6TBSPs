# -*- coding: utf-8 -*-

'''
To build a compressed hashtable for protein databases.
Input: one or more multi-fasta protein database(s)
Output: compressed dictionary (hashtable)

Author:
	Mei-Yu Lai

E-mail:
	mlai22@jhu.edu

Usage:
	$ python 6tbsps.py protein.fa

Attributes:
	Parse an input FASTA and translate the DNA sequences in 6-frame, 
	output a single FASTA file of compressed protein sequences.
'''
#%%
import argparse
import pickle
#%%
def parse_fasta(filename):
	'''
	Parse protein sequences from FASTA file

	Args:
        filename (str):   file name of the FASTA file

	Returns:
		ID_seq (dict):    a dictionary with keys reference IDs and values 
                            protein sequences
	'''
	ID_seq = {}
	with open(filename, 'r') as fh:
		while True:
			line = fh.readline()	
			if len(line) == 0:
				break
			if line[0] == '>':
				readID = line.rstrip()
				protein_seq = ''
			else:
				protein_seq += line.rstrip()
			ID_seq[readID] = protein_seq
	return ID_seq
#%%
def protein_kmer_table(seqs, k):
	'''
	Create a dictionary of k-mer protein sequence

	Args:
		seqs (dict): read dictionary
		k (int) : an interger k for k-mer

	Returns:
		dict: a dictionary that mapts each k-mer to the set of names of reads and the position
	'''
	table = {}
	for readID, sequence in seqs.items():
		for loc in range(len(sequence) - k + 1):
			kmer = sequence[loc:loc+k]
			if kmer not in table:
				table[kmer] = []
			table[kmer].append((readID, loc))
	return table
#%%
def main():
	'''The main function for 6TBSPS-build'''

	parser = argparse.ArgumentParser(prog = '6tbsps-build',
									 description = 'To build a compressed hashtable for protein databases.')
	parser.add_argument('proteinDBS', metavar = 'protein.fa', help = 'protein FASTA filename')
	parser.add_argument('-i', nargs = 1, type = int, help = 'k-mer')
	args = parser.parse_args()
	filename = args.proteinDBS
	k = int(args.i[0])

	with open('protein_dictionary.pickle', 'wb') as handle:
		pickle.dump(parse_fasta(filename), handle, protocol = pickle.HIGHEST_PROTOCOL)	
	with open('protein_dictionary.pickle', 'rb') as handle:
		protein_seq = pickle.load(handle)	

	with open('kmer_protein_dictionary.pickle', 'wb') as handle:
		pickle.dump(protein_kmer_table(protein_seq, k), handle, protocol = pickle.HIGHEST_PROTOCOL)
	with open('kmer_protein_dictionary.pickle', 'rb') as handle:
		kmer_protein_seq = pickle.load(handle)

	return
#%%
if __name__ == '__main__':
	main()
