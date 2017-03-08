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


def forward_reads_given_anc(r, a, rho, mu, eps):
    '''This takes an array of reads r, an array of ancestors a, rho, mu
    and computes p(r|a) using the simplest forward algorithm
    '''
    # f[i, j, k, v]
    # a[k, j]
    
    K, L = a.shape
    N, Lr = r.shape
    assert L == Lr
    M = B**K # cardinality of ancestors
    f = np.zeros([N, L, K, B])
    
    j = 0 # Init
    for i in range(N):
        for k in range(K):
            for v in range(B):
                if r[i, j] == v: f[i, j, k, v] = 1.-(B-1)*eps
                else: f[i, j, k, v] = eps
                if v == a[k, j]: f[i, j, k, v] *= 1.-(B-1)*mu
                else: f[i, j, k, v] *= mu
                f[i, j, k, v] /= K # prior

    for j in range(1, L): # Recurs
        for i in range(N):
            for k in range(K):
                for v in range(B):
                    for l in range(K):
                        if k == l: f[i, j, k, v] += (1.-(K-1)*rho) * sum(f[i, j-1, l, :])
                        else: f[i, j, k, v] += rho * sum(f[i, j-1, l, :])
                    
                    if r[i, j] == v: f[i, j, k, v] *= 1.-(B-1)*eps
                    else: f[i, j, k, v] *= eps
                    if v == a[k, j]: f[i, j, k, v] *= 1.-(B-1)*mu
                    else: f[i, j, k, v] *= mu

    f_ijk = np.sum(f, axis=3)    
    f_ij = np.sum(f_ijk, axis=2)        
    ff = np.prod(f_ij, axis=0)                
    
    return ff

def main():
    
    '''What does a main do?
    '''
    import sys
    from Bio import SeqIO
    from itertools import product
    
    rho = 0.272
    mu = 0.123
    eps = 0.1
    
    r1 = np.array([0,1,1])
    r2 = np.array([1,1,0])
    r3 = np.array([1,0,0])
    r = np.array([r1, r2, r3])
    
    a1 = [0,0,1]
    a2 = [1,1,0]
    a3 = [1,1,0]
    a = np.array([a1, a2, a3])
    
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
    pr = forward_reads_given_anc(r, a, rho, mu ,eps)
    print('\nprob is', pr, file=sys.stderr)
    # This computes \sum_R p(R|A)
    p = np.zeros(L)
    all_R = constellation(B, L, N)
    for i, rc in enumerate(all_R):
        p_here = forward_reads_given_anc(rc, a, rho, mu, eps)
        #print(rc, p_here, file=sys.stderr)
        p += p_here
    print('----------> check \tp =', p, file=sys.stderr)
    cp = np.array([B**(N*(L-j-1)) for j in range(L)], dtype=np.float32)
    print('----------> should be \tp =', cp, file=sys.stderr)
    
    
    
    
    
    sys.exit()
        
    
    
    
    
    filename = '../data/toy_reads.fasta'
    reads_in = list(SeqIO.parse(open(filename), 'fasta'))
    N = len(reads_in)
    L = len(reads_in[0].seq)
    K = 2
    rho = 0.0015
    mu = 0.17
    eps = 0.01
    r = np.empty([N, L], dtype=int)
    print >> sys.stderr, 'N=%d, L=%d' % (N, L)
    header = [s.id for s in reads_in]
    for i, s in enumerate(reads_in):
        r_str = s.seq.tostring()
        r[i] = [dna2i[b] for b in r_str]
    # r is the reads r[i, j] = read i at position j
    build_int2anc(K)
    #for i, a in enumerate(int2anc): print i, a
    #print ''
    build_anc2int(K)
    #for k, v in anc2int.items(): print k, v
    #print 'xxx'
    for i in range(B**K):
        assert i == anc2int[''.join(map(str, int2anc[i]))]
    N = 2
    L = 2
    method ='full_joint'
    method = 'forward'
    all_A = constellation(B, L, K)
    all_R = constellation(B, L, N)
    count = 0
    p = 0.0
    if method == 'full_joint':
        all_Z = constellation(K, L, N)
        all_H = constellation(B, L, N)
        for r in all_R:
            print (count)
            for z in all_Z:
                for a in all_A:
                    for h in all_H:
                        count += 1
                        p += full_joint(r, z, a, h, mu, eps, rho)

        
    elif method == 'full_joint_i':
        all_Zi = constellation(K, L, 1)
        all_Hi = constellation(B, L, 1)
        for r in all_R:
            print ('count=', count)
            for a in all_A:
                q = 1
                for ri in r:
                    s = 0.0
                    for z in all_Zi:
                        for h in all_Hi:
                            count += 1
                            s += full_joint_i(ri.reshape([1, L]), z, a,
                                              h.reshape([1, L]), K, L, 1,
                                              mu, eps, rho)
                    q *= s
                p += q

    elif method == 'forward':
        for r in all_R:
            count += 1
            p += forward(r, K, rho, mu, eps)

    else: pass
    
if __name__ == '__main__':
    main()
