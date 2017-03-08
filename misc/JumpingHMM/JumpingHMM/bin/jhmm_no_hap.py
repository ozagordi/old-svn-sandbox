#!/usr/bin/env python
from __future__ import print_function
'''
Now developed with Eclipse
'''

import numpy as np
dna2i = {'A':0, 'C':1, 'G':2, 'T':3}
B = 2 # change to something automatic

anc2int = {}
int2anc = []

def build_int2anc(K):
    '''Returns a map from integer to ancestor vector
    '''
    for x in range(B**K):
        x2 = x
        a = []
        for k in range(K):
            a.append(None)
        #a[0] = x2 % B
        for k in range(K-1, -1, -1):
            a[k] = x2 % B
            x2 /= B
        int2anc.append(a)

def build_anc2int(K):
    from itertools import product
    all_strings = product(''.join(map(str, range(B))), repeat=K)
    for als in all_strings:
        s = ''.join(als)
        ih = sum( [int(v)*(B**(i)) for i, v in enumerate(als)] )
        ih = 0
        for i in range(K):
            ih += int(s[i])*(B**(K-i-1))
        anc2int[s] = ih
        
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

def forward_single_read(r, K, rho, mu):
    '''This takes a single read r, K (number of haplotypes), rho,
    mu and computes p(r|a) using the simplest forward algorithm
    '''
    # f[j, k, x]
    L = len(r)
    M = B**K # cardinality of ancestors
    f = np.zeros([L, K, M])
    ff = np.zeros(L)
    
    j = 0 # Init
    for k in range(K):
        for x in range(M):
            if r[j] == int2anc[x][k]: f[j, k, x] = 1.-(B-1)*mu
            else: f[j, k, x] = mu
            f[j, k, x] /= K # prior on the initial state
            f[j, k, x] /= B**K # prior on the ancestors

    sh = np.apply_over_axes(np.sum, f[j, :, :], [0, 1])
    assert sh.size == 1            
    ff[j] = sh[0][0]

    for j in range(1, L): # Recurs
        for k in range(K):
            for x in range(M):# so far to cycle through all entries of f
                # now to compute the single entry
                for l in range(K):
                    psum = sum(f[j-1, l, :])
                    if k == l: psum *= (1.-(K-1)*rho)
                    else: psum *= rho
                    f[j, k, x] += psum                
                if r[j] == int2anc[x][k]: f[j, k, x] *= 1.-(B-1)*mu
                else: f[j, k, x] *= mu
                f[j, k, x] /= B**K # prior on the ancestors
        sh = np.apply_over_axes(np.sum, f[j, :, :], [0, 1])
        assert sh.size == 1
        ff[j] = sh[0][0]
    #ff /= B**K
    return ff

def forward_reads(r, K, rho, mu):
    '''This takes a set of reads r, K (number of haplotypes), rho,
    mu and computes p(r|a) using the forward algorithm
    '''
    # f[i, j, k, x]
    N, L = r.shape
    M = B**K # cardinality of ancestors
    f = np.zeros([N, L, K, M])
    ff = np.zeros(L)
    
    j = 0 # Init
    for x in range(M):
        p_jx = 1.0
        for i in range(N):
            for k in range(K):
                if r[i, j] == int2anc[x][k]: f[i, j, k, x] = 1.-(B-1)*mu
                else: f[i, j, k, x] = mu
                f[i, j, k, x] /= K # prior on the initial state
                #f[i, j, k, x] /= B**K # prior on the ancestors
                # entry computed    
            s_ijx = sum(f[i, j, :, x])
            p_jx *= s_ijx
        ff[j] += p_jx
    
    for j in range(1, L): # Recurs
        for x in range(M):
            p_jx = 1.0
            for i in range(N):
                for k in range(K):# so far to cycle through all entries of f
                    # now to compute the single entry
                    for l in range(K):
                        psum = sum(f[i, j-1, l, :])
                        if k == l: psum *= (1.-(K-1)*rho)
                        else: psum *= rho
                        f[i, j, k, x] += psum
                    if r[i, j] == int2anc[x][k]: f[i, j, k, x] *= 1.-(B-1)*mu
                    else: f[i, j, k, x] *= mu
                    f[i, j, k, x] /= B**K # prior on the ancestors
                    # entry computed
                s_ijx = sum(f[i, j, :, x])
                p_jx *= s_ijx
            ff[j] += p_jx
            
    ff /= B**K
    return ff

def forward_reads_invert(r, K, rho, mu):
    '''This takes a set of reads r, K (number of haplotypes), rho,
    mu and computes p(r|a) using the forward algorithm
    '''
    # f[i, j, k, x]
    N, L = r.shape
    M = B**K # cardinality of ancestors
    f = np.zeros([N, L, K, M])
    ff = np.zeros(L)
    
    j = 0 # Init
    for x in range(M):
        p_jx = 1.0
        for i in range(N):
            for k in range(K):
                if r[i, j] == int2anc[x][k]: f[i, j, k, x] = 1.-(B-1)*mu
                else: f[i, j, k, x] = mu
                f[i, j, k, x] /= K # prior on the initial state
                #f[i, j, k, x] /= B**K # prior on the ancestors
                # entry computed    
            s_ijx = sum(f[i, j, :, x])
            p_jx *= s_ijx
        ff[j] += p_jx
    
    for j in range(1, L): # Recurs
        for x in range(M):
            p_jx = 1.0
            for i in range(N):
                for k in range(K):# so far to cycle through all entries of f
                    # now to compute the single entry
                    for l in range(K):
                        psum = sum(f[i, j-1, l, :])
                        if k == l: psum *= (1.-(K-1)*rho)
                        else: psum *= rho
                        f[i, j, k, x] += psum
                    if r[i, j] == int2anc[x][k]: f[i, j, k, x] *= 1.-(B-1)*mu
                    else: f[i, j, k, x] *= mu
                    f[i, j, k, x] /= B**K # prior on the ancestors
                    # entry computed
                s_ijx = sum(f[i, j, :, x])
                p_jx *= s_ijx
            ff[j] += p_jx
            
    ff /= B**K
    return ff



def main():    
    '''What does a main do?
    '''
    import sys
    K = 2
    build_int2anc(K)
    build_anc2int(K)
    for i in range(B**K):
        assert i == anc2int[''.join(map(str, int2anc[i]))]
        
    N = 1
    rho = 0.1244
    mu = 0.2463
    
    r1 = np.array([0,1,1,0,1])        
    # check single read r
    N, L = 1, r1.shape[0]
    pr = forward_single_read(r1, K, rho, mu)
    print('single read is\t', r1, '\tsingle prob is', pr, file=sys.stderr)
    
    # This computes \sum_r p(r|A)
    p = np.zeros(L)
    all_r = constellation(B, L, N)
    for i, rc in enumerate(all_r):
        tp = forward_single_read(rc[0], K, rho, mu)
        #print(rc[0], tp, file=sys.stderr)
        p += tp
    print('sum over all single reads gives', file=sys.stderr)
    print('----------> check \tp =', p, file=sys.stderr)
    cp = np.array([B**(L-j-1) for j in range(L)], dtype=np.float32)
    print('----------> should be \tp =',cp, file=sys.stderr)
    
    # Now for N reads
    print('\n##########\n', file=sys.stderr)
    # check set of reads
    N, L = 3, 2# r1.shape[0]
    all_R = constellation(B, L, N)
    p = np.zeros(L)
    for r in all_R:
        pr = forward_reads(r, K, rho, mu)
        p += pr
    print('---> check \tp =', p, file=sys.stderr)
    cp = np.array([B**(N*(L-j-1)) for j in range(L)], dtype=np.float32)
    print('------> should be \tp =',cp, file=sys.stderr)
    
    
    
if __name__ == '__main__':
    main()
