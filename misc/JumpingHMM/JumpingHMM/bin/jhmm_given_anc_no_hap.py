#!/usr/bin/env python
from __future__ import print_function
'''
Now developed with Eclipse
'''

import numpy as np
dna2i = {'A':0, 'C':1, 'G':2, 'T':3}
B = 2 # change to something automatic

def constellation(b, L, N, toprint=False):
    '''Returns all possible datasets of N vectors of
    length L, on the alphabet B=range(b)
    '''
    from itertools import product
    X = product(range(b), repeat=L)
    Y = product(X, repeat=N)
    Yl = np.array(list(Y))
    assert Yl.shape[0] == b**(L*N)
    if toprint:
        print('Generated %d datasets' % b**(L*N), file=sys.stderr)
    
    return Yl

def forward_reads_given_anc(r, a, rho, mu):
    '''This takes an array of reads r, an array of ancestors a, rho, mu
    and computes p(r|a) using the simplest forward algorithm
    '''
    # f[i, j, k]
    # a[k, j]
    
    K, L = a.shape
    N, Lr = r.shape
    assert L == Lr
    M = B**K # cardinality of ancestors
    f = np.zeros([N, L, K])
    ff = np.zeros([L])
    
    j = 0 # Init
    fj = 1.0
    for i in range(N):
        for k in range(K):
            if r[i, j] == a[k, j]: f[i, j, k] = 1.-(B-1)*mu
            else: f[i, j, k] = mu
            f[i, j, k] /= K # prior
        fij = sum(f[i, j, :])
        fj *= fij
    ff[j] = fj

    for j in range(1, L): # Recurs
        fj = 1.0
        for i in range(N):
            for k in range(K):
                for l in range(K):
                    if k == l: f[i, j, k] += (1.-(K-1)*rho) * f[i, j-1, l]
                    else: f[i, j, k] += rho * f[i, j-1, l]
                    
                if r[i, j] == a[k, j]: f[i, j, k] *= 1.-(B-1)*mu
                else: f[i, j, k] *= mu
                     
            fij = sum(f[i, j, :])
            fj *= fij
        ff[j] = fj
        
    return ff


def forward_single_read_given_anc(r, a, rho, mu):
    '''This takes a single read r, a set of ancestors a, rho, mu
    and computes p(r|a) using the simplest forward algorithm
    '''
    # f[j, k]
    # a[k, j]
    K, L = a.shape
    assert L == len(r)
    M = B**K # cardinality of ancestors
    f = np.zeros([L, K])
    ff = np.zeros([L])
    
    j = 0 # Init
    for k in range(K):
        if r[j] == a[k, j]: f[j, k] = 1.-(B-1)*mu
        else: f[j, k] = mu
        f[j, k] /= K # prior
                
    ff[j] = sum(f[j,:])

    for j in range(1, L): # Recurs
        for k in range(K):
            for l in range(K):
                if k == l: f[j, k] += (1.-(K-1)*rho) * f[j-1, l]
                else: f[j, k] += rho * f[j-1, l]
                
            if r[j] == a[k, j]: f[j, k] *= 1.-(B-1)*mu
            else: f[j, k] *= mu
                    
        ff[j] = sum(f[j,:])
    
    return ff

def forward_reads_1(r, rho, mu, K):
    ''' This calls forward_reads_given_anc
    '''
    n, L = r.shape
    all_a = constellation(B, L, K)
    p = np.zeros(L)
    for i, a in enumerate(all_a):
        p += forward_reads_given_anc(r, a, rho, mu)
    return p/B**(K*L)


def main():
    
    '''What does a main do?
    '''
    import sys
    
    rho = 0.272
    mu = 0.123
    
    r1 = np.array([0,1,1])
    r2 = np.array([1,1,0])
    r3 = np.array([1,0,0])
    r = np.array([r1, r2, r3])
    
    a1 = [0,0,0]
    a2 = [1,0,1]
    a3 = [1,0,0]
    a = np.array([a1, a2, a3])
    
    # check single read r given a set A of ancestors
    N, L = 1, r1.shape[0]
    pr = forward_single_read_given_anc(r1, a, rho, mu)
    print('single read is\t', r1, file=sys.stderr)
    print('', file=sys.stderr)
    for i, aa in enumerate(a):
        print('ancestor', i, 'is\t', aa, file=sys.stderr)    
    print('\nsingle prob is', pr, file=sys.stderr)
    
    # This computes \sum_r p(r|A)
    p = np.zeros(L)
    all_r = constellation(B, L, 1)
    for i, rc in enumerate(all_r):
        p_here = forward_single_read_given_anc(rc[0], a, rho, mu)
        p += p_here
        #print(rc[0], p_here, file=sys.stderr)
    print('sum over all reads gives', file=sys.stderr)
    print('----------> check \tp =', p, file=sys.stderr)
    cp = np.array([B**(L-j-1) for j in range(L)], dtype=np.float32)
    print('----------> should be \tp =', cp, file=sys.stderr)
    
    # Now for N reads
    print('\n##########\n', file=sys.stderr)
    N, L = r.shape
    # R is a set of reads
    for i, rr in enumerate(r):
        print('read', i, 'is\t', rr, file=sys.stderr)
    print('', file=sys.stderr)
    for i, aa in enumerate(a):
        print('ancestor', i, 'is\t', aa, file=sys.stderr)
    # P(R|A)
    pr = forward_reads_given_anc(r, a, rho, mu)
    print('\nprob is', pr, file=sys.stderr)
    # This computes \sum_R p(R|A)
    p = np.zeros(L)
    all_R = constellation(B, L, N)
    for i, rc in enumerate(all_R):
        p_here = forward_reads_given_anc(rc, a, rho, mu)
        #print(rc, p_here, file=sys.stderr)
        p += p_here
    print('----------> check \tp =', p, file=sys.stderr)
    cp = np.array([B**(N*(L-j-1)) for j in range(L)], dtype=np.float32)
    print('----------> should be \tp =', cp, file=sys.stderr)
 
    print('\n-----------------\n', file=sys.stderr)
    K = 3
    L = len(r1)
    all_a = constellation(B, L, K)
    p1 = 0.0
    p2 = 0.0
    pboth = 0.0
    count = 0
    for a in all_a:
        count += 1
        pa1 = forward_single_read_given_anc(r1, a, rho, mu)[-1]
        pa2 = forward_single_read_given_anc(r2, a, rho, mu)[-1]
        pboth += pa1*pa2
        p1 += pa1
        p2 += pa2
    print('It is\t', count, pboth, p1*p2, file=sys.stderr)
    print('Should\t', B**(K*L),'some.numbers', B**(2*L*(K-1)), file=sys.stderr)
    sys.exit()   
    K = 2
    N = 2
    # This computes
    all_r = constellation(B, L, N)
    p = np.zeros(L)
    for i, rs in enumerate(all_r):
        p += forward_reads_1(rs, rho, mu, K)
    print('-------> forw_reads_1 \tp =', p, file=sys.stderr)
    cp = np.array([B**(N*(L-j-1)) for j in range(L)], dtype=np.float32)
    print('-------> should be \tp =', cp, file=sys.stderr)
    #print('----------> B**(K*L) =', B**(K*L), file=sys.stderr)
    #print('count =', count, file=sys.stderr)
    
if __name__ == '__main__':
    main()
