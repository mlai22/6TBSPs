""" Evaluator for 6TBSPs and BlastX output
​
Author: 
    Zitong He
​
Email: 
    hezt@jhu.edu
    
Usage:
    usage: python evalutaor.py [-h] -b file_path -6 directory_path
​
Attributes:
    None
​
"""

import argparse
import os
import pandas as pd


def ranking_loss(six_content, blastx_content):
    '''
    calculate ranking loss
    '''
    # get golden top result from blastx
    # print(blastx_content)
    print(six_content)
    blastx_dic = {} 
    for i in blastx_content:
        if i[0] not in blastx_dic.keys():
            blastx_dic[i[0]] = [[i[1], i[2]]]
        else:
            blastx_dic[i[0]].append([i[1], i[2]])
    for k in blastx_dic.keys():
        min_eval = min(blastx_dic[k], key=lambda x: x[1])[1]
        blastx_dic[k] = [i[0] for i in blastx_dic[k] if i[1] == min_eval]
    # put 6tbsps result to dict
    six_sid_dic = {} 
    six_evalue_dic = {}
    for i in six_content:
        if i[0] not in six_sid_dic.keys():
            six_sid_dic[i[0]] = [i[1]]
            six_evalue_dic[i[0]] = [i[2]]
        else:
            six_sid_dic[i[0]].append(i[1]) 
            six_evalue_dic[i[0]].append(i[2])
    # find blastx top results: blastx_dic[k][...] in 6tbsps results 
    ranking_loss_value = 0
    miss_counter = 0
    hit_counter = 0
    for qid, top_sid_list in blastx_dic.items():
        rank_candidates = []
        for top_sid in top_sid_list:
            if top_sid in six_sid_dic[qid]:
                # hit
                # there may exist tie: there are several same evalue good hits
                # break this tie by get the first hit's rank
                # among all same evalue good hits
                hit_loc = six_sid_dic[qid].index(top_sid)
                evalue = six_evalue_dic[qid][hit_loc]
                # .index() will return the first matched elem's index
                rank_candidates.append(six_evalue_dic[qid].index(evalue)) 
            else:
                # miss
                rank_candidates.append(len(six_sid_dic[qid]))
        min_rank = min(rank_candidates)
        print(rank_candidates)
        if min_rank == len(six_sid_dic[qid]):
            # miss
            miss_counter += 1
        else:
            hit_counter += 1 
        ranking_loss_value += min_rank
    print('Miss: ', miss_counter)
    print('Hit: ', hit_counter)
    return ranking_loss_value

def read_6tbsps_results(src_dir):
    '''
    read query results from 6tbsps output directory
    only read file with ext name: .summary
    clean the data
    '''
    six_content = []
    for file_name in os.listdir(src_dir):
        if os.path.splitext(file_name)[-1] == '.summary':
            with open(os.path.join(src_dir, file_name), 'r') as f:
                for line in f.readlines():
                    line_list = line.rstrip().split('\t')
                    qid = line_list[0].split(' ')[0]
                    sid = line_list[1].split(' ')[0]
                    evalue = float(line_list[-1])
                    line_list = [qid, sid, evalue]
                    six_content.append(line_list)
    return six_content

def read_blastx_results(src_file):
    '''
    read query results from 6tbsps output directory
    clean the data
    '''
    content = []
    with open(src_file, 'r') as f:
        for line in f.readlines():
            line_list = line.rstrip().split('\t')
            qid = line_list[0].split(' ')[0]
            sid = line_list[1].split(' ')[0]
            sid = sid.replace('@', ':')
            sid = sid.split('|')
            if len(sid) == 1:
                sid = sid[0]
            else:
                sid = sid[1]
            evalue = float(line_list[-2])
            line_list = [qid, sid, evalue]
            content.append(line_list)
    return content

def main():
    '''
    Main function of the evaluator
    Print ranking loss in console
    '''
    parser = argparse.ArgumentParser(prog = 'python evalutaor.py',
        description = '')
    parser.add_argument('-b', metavar='blastx_file_path',
        required = True,
        help='BLASTX query result file')
    parser.add_argument('-s', metavar='6tbsps_directory_path',
        required= True,
        help = '6TBSPs query result files directory')
    # parser.add_argument('-o', metavar='directory_path', \
        # required = True, help = 'output directory')
    args = parser.parse_args()
    blastx_file_path = args.b
    six_dir_path = args.s
    six_content = read_6tbsps_results(six_dir_path)
    blast_content = read_blastx_results(blastx_file_path)
    ranking_loss_value = ranking_loss(six_content, blast_content)
    print('Ranking Loss: ', ranking_loss_value)
    return ranking_loss_value


if __name__ == "__main__":
    main()