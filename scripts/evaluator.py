""" Evaluator for 6TBSPs and BlastX output
​
Author: 
    Zitong He
​
Email: 
    hezt@jhu.edu
    
Usage:
    TODO:
​
Attributes:
    TODO:
​
"""

import argparse
import os
import pandas as pd

def ranking_loss():
    
    return 0

def evaluate_6tbsps():
    print('6')
    
def evaluate_blastx(src_dir):
    for f in os.listdir(src_dir):
        src_f = os.path.join(src_dir, f)
        # dst_f = os.path.join(dst_dir, f + '_ms.fa')
        df = pd.read_csv(src_f, header=None)
        df = df[df.iloc[:, 2] == 1]
        f_content = []
        for index, row in df.iterrows():
            content = '>'
            for idx, c in enumerate(row):
                if idx == 6:
                    c = str(c).replace(' ', '')
                content += str(c) + ','
            f_content.append(content + '\n')
            f_content.append(row.iloc[1] + '\n')
        print(dst_f)
        with open(dst_f, 'w') as f:
            f.writelines(f_content)
    print('b') 

def main():
    '''
    Main function of the 6tbsps-query
    '''
    parser = argparse.ArgumentParser(prog = 'evalutaor', \
        description = '')
    parser.add_argument('-t', required = True, choices=['BLASTX', '6TBSPs'], \
        help='which tool you use for query, BLASTX or ')
    parser.add_argument('-o', metavar='evaluation_output_directory', \
        required = True, help = 'output directory')
    parser.add_argument('query_result_directory', nargs = '+', \
        help = 'query result files directory')
    args = parser.parse_args()
    out_dir = args.o
    in_dir = args.query_result_directory
    query_tool = args.t
    if query_tool == 'BLASTX':
        evaluate_blastx()
    else:
        evaluate_6tbsps()
    return

if __name__ == "__main__":
    main()