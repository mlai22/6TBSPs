import sys
import os
import numpy as np

def test_score_matrix(xc, yc):
        ''' Cost function: 2 to match, -6 to gap, -4 to mismatch '''
        if xc == yc: return 2 # match
        if xc == '-' or yc == '-': return -6 # gap
        return -4

class LocalAlignment:
    def __init__(self, seq_x, seq_y, score_matrix):
        '''init class paras and score_matrix'''
        self.V = np.zeros((len(seq_x)+1, len(seq_y)+1), dtype=int)
        self.x = seq_x
        self.y = seq_y
        self.s = score_matrix

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
        """ Trace back from given cell in local-alignment matrix V """
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
        print(self.V) 
        print(self.V_max)
        print(self.alignment)


if __name__ == "__main__":
    la = LocalAlignment('GGTATGCTGGCGCTA', 'TATATGCGGCGTTT', test_score_matrix)
    la.fill_matrix()
    la.traceback()
    la.display()