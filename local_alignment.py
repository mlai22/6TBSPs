import sys
import os
import numpy as np


class LocalAlignment:
    def __init__(self, seq_x, seq_y):
        '''init class paras and score_matrix'''
        self.V = np.zeros((len(seq_x)+1, len(seq_y)+1), dtype=int)
        self.x = seq_x
        self.y = seq_y
        self.s = self.test_score_matrix

    def test_score_matrix(self, xc, yc):
        ''' Cost function: 2 to match, -6 to gap, -4 to mismatch '''
        if xc == yc: return 2 # match
        if xc == '-' or yc == '-': return -6 # gap
        return -4

    def fill_matrix(self):
        ''' Calculate local alignment values of sequences x and y using
            dynamic programming. Return maximal local alignment value. '''
        for i in range(1, len(self.x)+1):
            for j in range(1, len(self.y)+1):
                self.V[i, j] = max(self.V[i-1, j-1] + self.s(self.x[i-1], self.y[j-1]), # diagonal
                                   self.V[i-1, j  ] + self.s(self.x[i-1], '-'),         # vertical
                                   self.V[i  , j-1] + self.s('-',         self.y[j-1]), # horizontal
                                   0)                                                   # empty
        argmax = np.where(self.V == self.V.max())
        self.V_max = int(self.V[argmax])
    
    def traceback(self):
        pass

    def display(self):
        print(self.V) 
        print(self.V_max)


if __name__ == "__main__":
    la = LocalAlignment('GGTATGCTGGCGCTA', 'TATATGCGGCGTTT')
    la.fill_matrix()
    la.display()