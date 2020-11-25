""" Local Alignment using linear gap, Alg: SW
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
    if xc == yc: return 2 # match
    if xc == '-' or yc == '-': return -6 # gap
    return -4

class LocalAlignment:
    def __init__(self, seq_x, seq_y, score_matrix):
        '''
        Init class parameters and score_matrix
        Parameters include the seq_x and seq_y, which need to be aligned.

        Args:
            seq_x (string): sequence x
            seq_y (string): sequence y
    ​
        Returns:
            None

        '''
        self.V = np.zeros((len(seq_x)+1, len(seq_y)+1), dtype=int)
        self.x = seq_x
        self.y = seq_y
        self.s = score_matrix

    def fill_matrix(self):
        ''' Calculate local alignment values of sequences x and y using
            dynamic programming. Return maximal local alignment value.

        Args:
            self:
                x (string): seq_x
                y (string): seq_y
                s (func): score_matrix function
        
        Returns:
            self:
                V (numpy.array): matrix which needed to be filled
                V_max: the max value of every matrix elements
        '''
        for i in range(1, len(self.x)+1):
            for j in range(1, len(self.y)+1):
                self.V[i, j] = max(self.V[i-1, j-1] + self.s(self.x[i-1], self.y[j-1]), # diagonal
                                   self.V[i-1, j  ] + self.s(self.x[i-1], '-'),         # vertical
                                   self.V[i  , j-1] + self.s('-',         self.y[j-1]), # horizontal
                                   0)                                                   # empty
        argmax = np.where(self.V == self.V.max())
        self.V_max = int(self.V[argmax])
    
    def traceback(self):
        """
        Trace back from given cell in local-alignment matrix V
        
        Args:
            self:
                V (numpy.array): the filled matrix
                s (func): score_matrix function
                x (string): seq_x
                y (string): seq_y
        
        Returns:
            self:
                xscript (string): transcript of alignment
                alignment (strings): three lines, seq_x, alignment symbol, seq_y

        """
        # get i, j for maximal cell
        i, j = np.unravel_index(np.argmax(self.V), self.V.shape)
        # traceback
        xscript, alx, aly, alm = [], [], [], []
        while (i > 0 or j > 0) and self.V[i, j] != 0:
            diag, vert, horz = 0, 0, 0
            if i > 0 and j > 0:
                diag = self.V[i-1, j-1] + self.s(self.x[i-1], self.y[j-1])
            if i > 0:
                vert = self.V[i-1, j] +   self.s(self.x[i-1], '-')
            if j > 0:
                horz = self.V[i, j-1] +   self.s('-', self.y[j-1])
            if diag >= vert and diag >= horz:
                match = self.x[i-1] == self.y[j-1]
                xscript.append('M' if match else 'R')
                alm.append('|' if match else ' ')
                alx.append(self.x[i-1]); aly.append(self.y[j-1])
                i -= 1; j -= 1
            elif vert >= horz:
                xscript.append('D')
                alx.append(self.x[i-1]); aly.append('-'); alm.append(' ')
                i -= 1
            else:
                xscript.append('I')
                aly.append(self.y[j-1]); alx.append('-'); alm.append(' ')
                j -= 1
        self.xscript = (''.join(xscript))[::-1]
        self.alignment = '\n'.join(map(lambda l: ''.join(l), [alx[::-1], alm[::-1], aly[::-1]]))

    def display(self):
        '''
        Display some parameters of this function.
        '''
        print(self.V) 
        print(self.V_max)
        print(self.alignment)


if __name__ == "__main__":
    '''
    Test codes
    '''
    la = LocalAlignment('GGTATGCTGGCGCTA', 'TATATGCGGCGTTT', test_score_matrix)
    la.fill_matrix()
    la.traceback()
    la.display()