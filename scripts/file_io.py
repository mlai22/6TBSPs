# -*- coding: utf-8 -*-
'''
To build a compressed hashtable for protein databases.
Input: one or more multi-fasta protein database(s)
Output: compressed dictionary (hashtable) of kmers and reference proteins

Author:
	Mei-Yu Lai

E-mail:
	mlai22@jhu.edu

Attributes:
	Parse an input FASTA and translate the DNA sequences in 6-frame, 
	output a single FASTA file of compressed protein sequences.
'''
#%%
import os
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
def read_dict(db_name):
	'''
	Read in a dictionary from a compressed pickle

	Args:
		db_name (str):			pre-index kmer database name

	Returns:
		A decompressed dictionary
	'''
	with open(db_name+'.pickle', 'rb') as handle:
		return pickle.load(handle)
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
