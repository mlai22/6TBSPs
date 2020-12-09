# -*- coding: utf-8 -*-

'''
6-frame translation
Input: one or more multi-fasta DNA reads
Output: Compressed protein sequences
		a single fasta file of 6 peptides append [-3:+3:1] label to original read ID
		e.g. read1|+2, read3|-1

Author:
	Mei-Yu Lai

E-mail:
	mlai22@jhu.edu

Usage:
	$ python six-frame_translation.py ref.fa

Attributes:
	Parse an input FASTA and translate the DNA sequences in 6-frame, 
	output a single FASTA file of compressed protein sequences.
'''
#%%
import sys
#%%
def reverse_complement(dna):
	'''
	Reverse complement of DNA string

	Args:
		dna (str): a DNA string

	Returns:
		str: a string of reversed complement of the input DNA string
	'''
	complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	revComp = []
	for base in dna:
		if base in ['A', 'T', 'C', 'G']:
			revComp.append(complement[base])
	return ''.join(revComp[::-1])

def six_frames(dna):
	'''
	Build six frames for each read

	Args:
		dna (str):
	
	Returns:
		dict: a dictionary with keys: frame ID and values: DNA sequences
	'''
	frame = {}

	'''Frame +1'''
	frame_1 = dna[0:]

	'''Frame +2'''
	frame_2 = dna[1:]

	'''Frame +3'''
	frame_3 = dna[2:]

	'''Frame -1'''
	frame_4 = reverse_complement(dna)[0:]

	'''Frame -2'''
	frame_5 = reverse_complement(dna)[1:]

	'''Frame -3'''
	frame_6 = reverse_complement(dna)[2:]

	frame[1] = frame_1
	frame[2] = frame_2
	frame[3] = frame_3
	frame[-1] = frame_4
	frame[-2] = frame_5
	frame[-3] = frame_6

	return frame

def transcription(dna):
	'''
	Transcribe DNA to RNA that Thymine (T) is replaced by Uracil (U)
	
	Args:
		dna (str): input a DNA sequence

	Returns:
		str: trascibe a DNA sequence to a RNA sequence
	'''
	rna_seq = dna.replace('T', 'U')
	return rna_seq

def translation(rna):
	'''
	Translate RNA codons to amino acid sequence
	
	Args:
		rna (str): input a RNA sequence

	Returns:
		dict: a dictinoary of keys: codon and values: amino acids letters
	'''
	amino_acid_table = {'AUG':'M', 'UUG':'L', 'GUG':'V', 'CUG':'L', 
						'AUA':'I', 'UUA':'L', 'GUA':'V', 'CUA':'L',
						'AUC':'I', 'UUC':'F', 'GUC':'V', 'CUC':'L', 
						'AUU':'I', 'UUU':'F', 'GUU':'V', 'CUU':'L', 
						'AGG':'R', 'UGG':'W', 'GGG':'G', 'CGG':'R', 
				  		'AGA':'R', 'UGA':'.', 'GGA':'G', 'CGA':'R', 
				  		'AGC':'S', 'UGC':'C', 'GGC':'G', 'CGC':'R',
				  		'AGU':'S', 'UGU':'C', 'GGU':'G', 'CGU':'R', 
				  		'ACG':'T', 'UCG':'S', 'GCG':'A', 'CCG':'P', 
				  		'ACA':'T', 'UCA':'S', 'GCA':'A', 'CCA':'P', 
				  		'ACC':'T', 'UCC':'S', 'GCC':'A', 'CCC':'P', 
				  		'ACU':'T', 'UCU':'S', 'GCU':'A', 'CCU':'P',
				  		'AAG':'K', 'UAG':'.', 'GAG':'E', 'CAG':'Q', 
				  		'AAA':'K', 'UAA':'.', 'GAA':'E', 'CAA':'Q', 
				  		'AAC':'N', 'UAC':'Y', 'GAC':'D', 'CAC':'H', 
				  		'AAU':'N', 'UAU':'Y', 'GAU':'D', 'CAU':'H'}			  						  						  						  				  						  						  		 				  		
	protein = ''
	loc = 0
	while(loc+3 <= len(rna)):
		codon = rna[loc:loc+3]
		# neglect uncommon condon
		if codon not in amino_acid_table:
			continue
		# stop codon
		if amino_acid_table[codon] == '.':
			break
		protein += amino_acid_table[codon]
		loc += 3

	return protein