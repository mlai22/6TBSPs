# -*- coding: utf-8 -*-
"""Seed and extend module.

This module seeds the query and extends into local alignment regions within the 
    subject. CAUTION: Very Naive!!!

Author: 
    Yuchen (Peter) Ge

Email: 
    yge15@jhmi.edu
    
Attributes:
    seed(str, int, int):  break the input sequence into multiple small seeds
        for fast exact matching
    
    extend(list, dict): find the possible regions in subject for local
        alignment
        
    naive_seed_and_extend(list, dict, int, int): wrapper function of seed and
        extend, given a seed length that is compatible with the pre-built
        dictionary and a gap length

"""
#%%
def seed(seq, slen, sgap):
    """Seed a peptide sequence given a seed length and a gap length. 
        Return a list of strings of seeds.
    
    Args:
        seq (str):      file name of the FASTA file
        slen (int):     length of seeds
        sgap (int):     gap between seeds

    Returns:
        seeds (list):   list of strings of seeds in a reading frame
    
    """ 
    seeds = []
    for i in range(0, len(seq)-slen+1, sgap):
        seed = seq[i:i+slen]
        seeds.append(seed)
    
    return seeds
#%%
def extend(seeds, subject):
    """Extend the seeds of exact matching in the protein subject into a region
        for local alignment.
        Return a list of all possible reference ids followed by start and end
        positions.
    
    Args:
        seeds (list):       list of strings of seeds in a reading frame
        subject (dict):     pre-index kmer dictionary (protein database)

    Returns:
        targets (list):     list of tuples of (ref_id, (start, end)) in the 
                            protein reference
    
    """ 
    targets = {}
    for seed in seeds:
        for hits in subject[seed]:
            ref_id = hits[0]
            start = hits[1]
            end = hits[1] + len(seed)
            
            if ref_id not in targets:
                targets[ref_id] = (start, end)
            else:
                if targets[ref_id][0] > start:
                    targets[ref_id][0] = start
                if targets[ref_id][1] < end:
                    targets[ref_id][1] = end
    
    return [(rid, pos) for rid, pos in targets.items()]
#%%
def naive_seed_and_extend(query, subject, seed_len, seed_gap):
    """Wrap up seed and extend mechanism.
        Return a list of all possible reference ids followed by start and end
        positions for each frame.
    
    Args:
        query (list):       list of tuples of (frame, sequence) for each read
        subject (dict):     pre-index kmer dictionary (protein database)
        seed_len (int):     length of seeds
        seed_gap (int):     gap between seeds

    Returns:
        results (list):     list of tuples of (ref_id, (start, end)) in the 
                            protein reference
    
    """     
    results = []
    for frame, seq in query:
        seeds = seed(seq, seed_len, seed_gap)
        results.append((frame, extend(seeds, subject)))
        
    return results
