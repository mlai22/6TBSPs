""" Local Alignment using affine gap, Alg: Gotch(Local)
​
Author: 
    Zitong He
​
Email: 
    hezt@jhu.edu
    
Usage:
    This is a class for internal pipline, no external usage.
    Internal usage see LocalAlignment.display()
​
Attributes:
    None
​
"""

import sys
import os
import numpy as np
from scripts.score_matrix import score_matrix

# MATCH = 5.
# MISMATCH = -1.
# GAPOPEN = -3.
# GAPEXT = -1.

def test_score_matrix(xc, yc):
    '''
    Given characters from seq_x and seq_y, compute the score based on the score_matrix.
    TODO: Need delete this one in the future, not used in this impelemtation
    Args:
        xc (char): the character from seq_x
        yc (char): the character from seq_y
​
    Returns:
        score (int): score based on xc and yc

    '''
    # if xc == yc: return MATCH # match
    # if xc == '-' or yc == '-': EOFError # gap
    # return MISMATCH
    pass

class LocalAlignment:
    def __init__(self, seq_x, seq_y, matrix_name='BLOSUM62', gap_open=-12., gap_ext=-4.):
        '''
        Init class parameters and score_matrix
        Parameters include the seq_x and seq_y, which need to be aligned.

        Args:
            seq_x (string): sequence x
            seq_y (string): sequence y
            matrix_name (string): name of the score matrix which we want to use
            gap_open (float): gap opening
            gap_ext (float): gap extension
        Returns:
            None

        self:
            X (numpy.array): upper matrix, _ will appear in seq_y
            Y (numpy.array): lower matrix, _ will appear in seq_x
            M (numpy.array): match matrix
            go (float): gap open panalty, for a score matrix, all go are same
            ge (float): gap extend panalty, for a score matrix, all ge are same
            score_matrix: score_matrix
        '''
        self.score_matrix = score_matrix(matrix_name)# only about match
        self.go = gap_open
        self.ge = gap_ext
        self.x = seq_x
        self.y = seq_y

        # initialize three matrix for affine sw
        dim_i = len(seq_x) + 1
        dim_j = len(seq_y) + 1
        self.X = np.zeros((dim_i, dim_j), dtype=float)
        self.Y = np.zeros((dim_i, dim_j), dtype=float)
        self.M = np.zeros((dim_i, dim_j), dtype=float)
        # initialize the first row/col of X, Y.
        # M's first row/col are zero, init already
        self.X[0,0], self.Y[0,0] = np.NaN, np.NaN
        for i in range(1, dim_i):
            self.X[i,0] = np.NaN
            self.Y[i,0] = -np.inf
        for j in range(1, dim_j):
            self.X[0,j] = -np.inf
            self.Y[0,j] = np.NaN

    def _match(self, i, j):
        '''
        i (int): index of current char in seq x
        j (int): index of current char in seq y
        return (float): score
        '''
        return float(self.score_matrix.loc[self.x[i-1], self.y[j-1]])
        # return float(score_matrix(self.x[i-1], self.y[j-1], self.matrix_name))
        # if self.x[i-1] == self.y[j-1]:
            # return MATCH
        # else:
            # return MISMATCH

    def fill_matrix(self):
        '''
        Calculate local alignment values of sequences x and y using
        dynamic programming. Return maximal local alignment value.

        self:
            X (numpy.array): upper matrix, _ will appear in seq_y
            Y (numpy.array): lower matrix, _ will appear in seq_x
            M (numpy.array): match matrix
            go (float): gap open panalty
            ge (float): gap extend panalt
            max_loc_list (list[tuple]): the list of maxs locs
            score: max value in self.M matrix, the best global alignment value
                   in local alignment
        return:
            self.score

        '''
        dim_i = len(self.x) + 1
        dim_j = len(self.y) + 1
        for j in range(1, dim_j):
            for i in range(1, dim_i):
                self.X[i][j] = max(self.M[i-1][j] + self.go + self.ge,
                                   self.X[i-1][j] + self.ge)
                self.Y[i][j] = max(self.M[i][j-1] + self.go + self.ge,
                                   self.Y[i][j-1] + self.ge)
                self.M[i][j] = max(self.M[i-1][j-1] + self._match(i, j),
                                   self.X[i][j],
                                   self.Y[i][j],
                                   0)
        argmax = np.where(self.M == self.M.max())
        self.max_loc_list = [(i, j) for i, j in zip(argmax[0], argmax[1])]
        self.score = int(self.M[self.max_loc_list[0]])
        return self.score
    
    def traceback(self):
        """
        Trace back from given cell in local-alignment matrix V
        
        self:
            X (numpy.array): upper matrix
            Y (numpy.array): lower matrix
            M (numpy.array): match matrix
            s (func): score_matrix function
            x (string): seq_x
            y (string): seq_y
            xscript_list (list[string]): transcript of alignment
            max_loc_list (list[tuple]): the list of maxs locs
            score (float): max value in self.M matrix, the best global alignment value
                    in local alignment
            align_seq_x_list (list[string]): aligned substring of sequence x, contains _ for gap
            align_seq_y_list (list[string]): aligned substring of sequence y, contains _ for gap
            max_loc_x_list (list[list[start, end]]): aligned location, start and end
            max_loc_y_list (list[list[start, end]]): aligned location, start and end
        return:
            self.score
            self.align_seq_x
            self.align_seq_y
            self.xscript

        """
        self.align_seq_x_list = []
        self.align_seq_y_list = []
        self.xscript_list = []
        self.max_loc_x_list = []
        self.max_loc_y_list = []

        for i, j in self.max_loc_list:
            align_seq_x = ''
            align_seq_y = ''
            xscript = ''
            max_loc_x = [0, 0]
            max_loc_y = [0, 0]
            max_loc_x[1] = i - 1
            max_loc_y[1] = j - 1
            v = self.M[i,j]
            while (v != 0):
                if (self.M[i][j] == self.Y[i][j]):
                    # horizontal first
                    # horizontal(j), ext
                    align_seq_x += '_'
                    align_seq_y += self.y[j-1]
                    xscript += ' '
                    if (self.Y[i,j] == self.Y[i, j-1] + self.ge):
                        # continue horizontal(j) ext
                        v = self.Y[i, j-1]
                    else:
                        # horizontal(j) open
                        v = self.M[i, j-1]
                    j -= 1
                elif (self.M[i][j] == self.X[i][j]):
                    # vertical then
                    # vertical(i) ext
                    align_seq_x += self.x[i-1]
                    align_seq_y += '_'
                    xscript += ' '
                    if (self.X[i,j] == self.X[i-1, j] + self.ge):
                        v = self.X[i-1, j]
                    else:
                        # horizontal(j) open
                        v = self.M[i-1, j]
                    i -= 1
                elif (self.M[i][j] == self.M[i-1][j-1] + self._match(i, j)):
                    # match last
                    # diagnal
                    align_seq_x += self.x[i-1]
                    align_seq_y += self.y[j-1]
                    if self.x[i-1] == self.y[j-1]:
                        xscript += '|'
                    else:
                        xscript += '*'
                    v = self.M[i-1][j-1]
                    i -= 1
                    j -= 1
            max_loc_x[0] = i
            max_loc_y[0] = j

            self.align_seq_x_list.append(align_seq_x[::-1])
            self.align_seq_y_list.append(align_seq_y[::-1])
            self.xscript_list.append(xscript[::-1])
            self.max_loc_x_list.append(max_loc_x)
            self.max_loc_y_list.append(max_loc_y)

        return self.score, \
               self.align_seq_x_list, self.align_seq_y_list, self.xscript_list, \
               self.max_loc_x_list, self.max_loc_y_list

    def display(self):
        '''
        Display some parameters of this function.
        Args:
            out_file: output file handle
        '''
        print('score: {}'.format(self.score))
        for i, loc in enumerate(self.max_loc_list):
            # print('x: {} y: {}'.format(loc[0], loc[1]))
            print('{:<4d} {} {:>4d}'.format(self.max_loc_x_list[i][0],
                                        self.align_seq_x_list[i],
                                        self.max_loc_x_list[i][1]))
            print('     {}     '.format(self.xscript_list[i]))
            print('{:<4d} {} {:>4d}'.format(self.max_loc_y_list[i][0],
                                        self.align_seq_y_list[i],
                                        self.max_loc_y_list[i][1]))
            print('')
    
    def display_file(self, out_file=None):
        '''
        Display some parameters of this function.
        '''
        print('score: {}'.format(self.score), file=out_file)
        for i, loc in enumerate(self.max_loc_list):
            # print('x: {} y: {}'.format(loc[0], loc[1]))
            print('{:<4d} {} {:>4d}'.format(self.max_loc_x_list[i][0],
                                            self.align_seq_x_list[i],
                                            self.max_loc_x_list[i][1]), file=out_file)
            print('    {}    '.format(self.xscript_list[i]), file=out_file)
            print('{:<4d} {} {:>4d}'.format(self.max_loc_y_list[i][0],
                                            self.align_seq_y_list[i],
                                            self.max_loc_y_list[i][1]), file=out_file)
            print('', file=out_file)


if __name__ == "__main__":
    '''
    Test codes
    '''
    # la = LocalAlignment('MISLIAALAVDRVIGMENAMPFNLPADLAWFKRNTLDKPVIMGRHTWESIG', 'SLNCIVAVSQNMGIGKNGDLPWPPLRNEFRYFQRMTTTSSVEGKQNLVIMGKKTWFSIPE')
    # la = LocalAlignment('SLNCIVAVSQNMGIGKNGDLPWPPLRNEFRYFQRMTTTSSVEGKQNLVIMGKKTWFSIPE', 'MISLIAALAVDRVIGMENAMPFNLPADLAWFKRNTLDKPVIMGRHTWESIG')
    # la = LocalAlignment('QRNTLDKPVIMGRHTWESI', 'QRMTTTSSVEGKQNLVIMGKKTWFSI')
    # la = LocalAlignment('NAMPFNL', 'NGDLPWPPL')
    la = LocalAlignment('SLIAALAVDRVIGMENAMPFNL', 'SLNCIVAVSQNMGIGKNGDLPWPPL')
    # la = LocalAlignment('GGTATGCTGGCGCTA', 'TATATGCGGCGTTT')
    # la = LocalAlignment('ACACACTA','AGCACACA')
    # la = LocalAlignment('ATTGAGC','ATGC')
    # la = LocalAlignment('ATGC','ATTGAGC')
    la.fill_matrix()
    la.traceback()
    la.display()