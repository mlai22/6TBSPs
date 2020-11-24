""" Local Alignment using affine gap, Alg: SW
​
Author: 
    Zitong He
​
Email: 
    hezt@jhu.edu
    
Usage:
    None
​
Attributes:
    None
​
"""

import sys
import os
import numpy as np

MATCH = 5.
MISMATCH = -1.
GAPOPEN = -3.
GAPEXT = -1.

def test_score_matrix(xc, yc):
    '''
    Given characters from seq_x and seq_y, compute the score based on the score_matrix.

    Args:
        xc (char): the character from seq_x
        yc (char): the character from seq_y
​
    Returns:
        score (int): score based on xc and yc

    '''
    if xc == yc: return MATCH # match
    if xc == '-' or yc == '-': EOFError # gap
    return MISMATCH

class LocalAlignment:
    def __init__(self, seq_x, seq_y, score_matrix):
        '''
        Init class parameters and score_matrix
        Parameters include the seq_x and seq_y, which need to be aligned.

        Args:
            seq_x (string): sequence x
            seq_y (string): sequence y
            score_matrix: for now, still use fixed match/mismatch/gap value​
        Returns:
            None

        self:
            X (numpy.array): upper matrix, _ will appear in seq_y
            Y (numpy.array): lower matrix, _ will appear in seq_x
            M (numpy.array): match matrix
            go (float): gap open panalty
            ge (float): gap extend panalty
        '''
        self.go = GAPOPEN
        self.ge = GAPEXT
        self.x = seq_x
        self.y = seq_y
        self.s = score_matrix # only about match

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
        TODO: future add socore matrix
        '''
        if self.x[i-1] == self.y[j-1]:
            return MATCH
        else:
            return MISMATCH

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
            xscript_list (list of string): transcript of alignment
            max_loc_list (list[tuple]): the list of maxs locs
            score (float): max value in self.M matrix, the best global alignment value
                    in local alignment
            align_seq_x_list (list of string): aligned substring of sequence x, contains _ for gap
            align_seq_y_list (list of string): aligned substring of sequence y, contains _ for gap
        return:
            self.score
            self.align_seq_x
            self.align_seq_y
            self.xscript

        """
        self.align_seq_x_list = []
        self.align_seq_y_list = []
        self.xscript_list = []

        for i, j in self.max_loc_list:
            align_seq_x = ''
            align_seq_y = ''
            xscript = ''
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
                        j -= 1
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
                        i -= 1
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

            self.align_seq_x_list.append(align_seq_x[::-1])
            self.align_seq_y_list.append(align_seq_y[::-1])
            self.xscript_list.append(xscript[::-1])

        return self.score, self.align_seq_x_list, self.align_seq_y_list, self.xscript_list

    def display(self):
        '''
        Display some parameters of this function.
        '''
        print('score: {}'.format(self.score))
        for i, loc in enumerate(self.max_loc_list):
            print('x: {} y: {}'.format(loc[0], loc[1]))
            print(self.align_seq_x_list[i])
            print(self.xscript_list[i])
            print(self.align_seq_y_list[i])
            print('')


if __name__ == "__main__":
    '''
    Test codes
    '''
    # la = LocalAlignment('MISLIAALAVDRVIGMENAMPFNLPADLAWFKRNTLDKPVIMGRHTWESIG', 'SLNCIVAVSQNMGIGKNGDLPWPPLRNEFRYFQRMTTTSSVEGKQNLVIMGKKTWFSIPE', test_score_matrix)
    # la = LocalAlignment('SLNCIVAVSQNMGIGKNGDLPWPPLRNEFRYFQRMTTTSSVEGKQNLVIMGKKTWFSIPE', 'MISLIAALAVDRVIGMENAMPFNLPADLAWFKRNTLDKPVIMGRHTWESIG', test_score_matrix)
    # la = LocalAlignment('QRNTLDKPVIMGRHTWESI', 'QRMTTTSSVEGKQNLVIMGKKTWFSI', test_score_matrix)
    # la = LocalAlignment('NAMPFNL', 'NGDLPWPPL', None)
    la = LocalAlignment('SLIAALAVDRVIGMENAMPFNL', 'SLNCIVAVSQNMGIGKNGDLPWPPL', None)
    # la = LocalAlignment('GGTATGCTGGCGCTA', 'TATATGCGGCGTTT', test_score_matrix)
    # la = LocalAlignment('ACACACTA','AGCACACA', test_score_matrix)
    # la = LocalAlignment('ATTGAGC','ATGC', None)
    # la = LocalAlignment('ATGC','ATTGAGC', None)
    la.fill_matrix()
    la.traceback()
    la.display()