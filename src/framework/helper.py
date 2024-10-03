"""Some helper functions
"""

import random

import numpy as np


NUC = set(list('ACGU'))
Mutation = {k: list(NUC-set(k)) for k in NUC}


def mutateSeq(seq, hamming=1):
    """For given sequence, returns a uniform and random mutatation of given hamming distance
    """
    n = len(seq)
    seqLst = list(seq)
    # Get positions to mutate
    for pos in np.random.permutation(n)[:hamming]:
        seqLst[pos] = random.choice(Mutation[seqLst[pos]])
    return ''.join(seqLst)


def cc_from_dbn(dbn):
    cc = []
    tmp = []
    for ind, x in enumerate(dbn):
        if x == '(':
            tmp.append(ind)
        elif x == ')':
            cc.append([tmp.pop(), ind])
        else:
            cc.append([ind])
    return cc
