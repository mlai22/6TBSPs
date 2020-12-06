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
def hamming_dist(str1, str2):
    """Calculate hamming distance between two strings of the same length
    
    Args:
        str1 (str):     the first string
        str2 (str):     the second string
    
    Returns:
        dist (int):     the hamming distance between str1 and str2
    
    """
    assert len(str1) == len(str2)

    dist = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            dist += 1
    
    return dist
#%%
def seed(seq, k):
    """Seed a peptide sequence given a seed length and a gap length. 
        Return a list of strings of seeds.
    
    Args:
        seq (str):      protein sequence
        k (int):     length of seeds

    Returns:
        seeds (list):   list of strings of seeds in a reading frame
    
    """ 
    if len(seq) < k:
        return []
    
    head = seq[:k]
    tail = seq[len(seq)-k:]
    seeds = [head, tail]
    
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
        if seed not in subject:
            continue
        for hits in subject[seed]:
            ref_id = hits[0]
            start = hits[1]
            end = hits[1] + len(seed)
            
            if ref_id not in targets:
                targets[ref_id] = [start, end]
            else:
                if targets[ref_id][0] > start:
                    targets[ref_id][0] = start
                if targets[ref_id][1] < end:
                    targets[ref_id][1] = end
    
    return [(rid, tuple(pos)) for rid, pos in targets.items()]
#%%
def naive_seed_and_extend(query, subject, k):
    """Wrap up seed and extend mechanism.
        Return a list of all possible reference ids followed by start and end
        positions for each frame.
    
    Args:
        query (str):        protein sequence derived from each reading frame
        subject (dict):     pre-index kmer dictionary (protein database)
        k (int):            length of seeds

    Returns:
        results (list):     list of tuples of (ref_id, (start, end)) in the 
                            protein reference
    
    """     
    seeds = seed(query, k)
    return extend(seeds, subject)
