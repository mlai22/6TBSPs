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


def score_matrix(matrix_name='BLOSUM62'):
    """
    Given characters from seq_x and seq_y, compute the score based on the score_matrix.
    Args:
        xc (char): the character from seq_x
        yc (char): the character from seq_y
    Returns:
        score (int): score based on xc and yc
    """
    path = './score_matrices/' + matrix_name + '.csv'
    matrix = pandas.read_csv(path)
    assert (matrix == matrix.T).all().all()
    # matrix.rename(index={'*': '-'}, columns={'*': '-'}, inplace=True)
    return matrix
    # return Matrix.loc[xc, yc]


def e_value_cal(m, n, S, name='BLOSUM62'):
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


result = pandas.DataFrame(columns=('#', 'query', 'reference', 'score', 'range', 'e-value', 'alignment'))
'''Just a suggestion on the final dataframe.'''

# print(score_matrix('A','-','BLOSUM45'))
