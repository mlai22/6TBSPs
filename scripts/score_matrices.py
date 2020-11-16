"""Value

This module store different score matrices for local alignment and a function for e-value calculate.

Author:
    Ruijing (Johnson) Zhang

Email:
    rzhang73@jh.edu


Attributes:
    Score matrix: BLOSUM45, BLOSUM62, BLOSUM80
    e_Value_cal(m,n,S)



"""

import math
import pandas

'''Cited from the original paper of BLOSUM, We used the same gap penalties for all matrices, 
-12 for the first residue in a gap, and -4 for subsequent residues in a gap. 
And in the matrix, '*' represents gap and the penalty is -5, -4 and -6 according to NCBI.

'''
BLOSUM45 = pandas.read_csv("BLOSUM45.csv")
assert (BLOSUM45 == BLOSUM45.T).all().all()
BLOSUM62 = pandas.read_csv("BLOSUM62.csv")
assert (BLOSUM62 == BLOSUM62.T).all().all()
BLOSUM80 = pandas.read_csv("BLOSUM80.csv")
assert (BLOSUM80 == BLOSUM80.T).all().all()


def e_Value_cal(m, n, S, name='BLOSUM62'):
    '''
    Actually, the parameters in the evalue calculating equation is quite hard to determine. We use data from this site
    https://bioinfo.lifl.fr/reblosum/all-matrices.html
    param m: length of the query
    param n: length of the reference
    param S: local alignment score
    return: e-value = kmne^-lambda*S

    '''
    if name == 'BLOSUM62':
        e = math.e
        E = 0.139042 * m * n * e ** (-0.320733 * S)
    if name == 'BLOSUM45':
        e = math.e
        E = 0.095168 * m * n * e ** (-0.231019 * S)
    if name == 'BLOSUM80':
        e = math.e
        E = 0.185160 * m * n * e ** (-0.350826 * S)
    return E


result = pandas.DataFrame(columns=('#', 'query', 'reference', 'score', 'range', 'e-value','alignment'))
'''Just a suggestion on the final dataframe.'''
