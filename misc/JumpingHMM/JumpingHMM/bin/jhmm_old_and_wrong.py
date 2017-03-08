#!/usr/bin/env python
'''
Now developed with Aptana
'''

import numpy as np
dna2i = {'A':0, 'C':1, 'G':2, 'T':3}
B = 2 # change to something automatic

anc2int = {}
int2anc = []
def f_brute_force(r, i, j, k, x, v, N, K, mu, eps, rho):
    '''
    '''
    L = len(r)
    all_A = constellation(B, L, K)
    all_R = constellation(B, L, N)
    all_Z = constellation(K, L, N)
    all_H = constellation(B, L, N)
    p = 0.0
    count = 0
    for R in all_R:
        for z in all_Z:
            for a in all_A:
                for h in all_H:
                    a2i = ''.join(map(str, a[j]))
                    print anc2int
                    if h[i,j] == v and \
                        anc2int[a2i] == x and \
                        z[i, j] == k and \
                        r[i, :j+1] == R[i, :j+1]:
                        count += 1
                        p += full_joint(r, z, a, h, mu, eps, rho)
    return p

                     
def full_joint_i(R, Z, A, H, K, L, N, mu, eps, rho):
    '''
    '''
    j = 0
    p = 1.0
    for k in range(K): p /= 1
    for i in range(N):
        p /= K
        if H[i, j] == A[Z[i, j], j] : p *= 1.-(B-1)*mu
        else: p *= mu
        if R[i, j] == H[i, j]: p *= 1.-(B-1)*eps
        else: p *= eps

    for j in range(1, L):
        for k in range(K): p /= 1
        for i in range(N):
            if Z[i,j] == Z[i, j-1]: p *= 1.-(K-1)*rho
            else: p *= rho
            if H[i,j] == A[Z[i, j], j] : p *= 1.-(B-1)*mu
            else: p *= mu
            if R[i, j] == H[i, j]: p *= 1.-(B-1)*eps
            else: p *= eps
    return p


def full_joint(R, Z, A, H, mu, eps, rho):
    '''a[k,j] k-th component. j-th pos
    '''
    j = 0
    p = 1.0
    K = len(list(A))
    L = len(A[0])
    N = len(list(R))
    for k in range(K): p /= B
    for i in range(N):
        p /= K
        if H[i,j] == A[Z[i, j], j] : p *= 1.-(B-1)*mu
        else: p *= mu
        if R[i, j] == H[i, j]: p *= 1.-(B-1)*eps
        else: p *= eps

    for j in range(1, L):
        for k in range(K): p /= B
        for i in range(N):
            if Z[i,j] == Z[i, j-1]: p *= 1.-(K-1)*rho
            else: p *= rho
            if H[i,j] == A[Z[i, j], j] : p *= 1.-(B-1)*mu
            else: p *= mu
            if R[i, j] == H[i, j]: p *= 1.-(B-1)*eps
            else: p *= eps
    return p

def constellation(b, L, N, toprint=False):
    '''Returns all possible datasets of N vectors of
    length L, on the alphabet B=range(b)
    '''
    from itertools import product
    X = product(range(b), repeat=L)
    Y = product(X, repeat=N)
    ## for i, y in enumerate(Y):
    ##     if toprint:
    ##         for x in y:
    ##             print x
    ##         print ''
    ##     count = i+1
    print 'Generated %d datasets' % b**(L*N)
    ## assert count == b**(L*N)
    return np.array(list(Y))

def baseb(i, K, b):
    x2 = i
    a = []
    for k in range(K):
        a.append(None)
    for k in range(K-1, -1, -1):
        a[k] = x2 % b
            
        x2 /= b
    int2anc.append(a)
    return a
    
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

def forward(r, K, rho, mu, eps):
    '''Forward takes the reads r
    K number of haplotypes
    ...
    ...
    ...
    '''
    # f[i,j,k,x,v]
    N, L = r.shape
    M = B**K # cardinality of ancestors
    f = np.zeros([N, L, K, M, B])
    ff = np.zeros([L])
    
    j = 0 # Init
    for i in range(N):
        for k in range(K):
            for x in range(M):
                for v in range(B):
                    if r[i, j] == v: f[i, j, k, x, v] = 1.-(B-1)*eps
                    else: f[i, j, k, x, v] = eps
                    
                    if v == int2anc[x][k]: f[i, j, k, x, v] *= 1.-(B-1)*mu
                    else: f[i, j, k, x, v] *= mu
                    
                    f[i, j, k, x, v] /= K # prior
                    #f[i, j, k, x, v] /= B**K
                    # fbf = f_brute_force(r, i, j, k, x, v, N, K, mu, eps, rho)
                    # print f[i, j, k, x, v]/fbf

    for i in range(N):
        ff[j] = 1.0
        s_xi = 0.0
        for x in range(M):
            for k in range(K):
                for v in range(B):
                    s_xi += f[i, j, k, x, v]
        ff[j] *= s_xi

    for j in range(1, L): # Recurs
        for i in range(N):
            for k in range(K):
                for x in range(M):
                    for v in range(B):
                        for l in range(K):
                            if k == l: s_jikxv = 1.-(K-1)*rho
                            else: s_jikxv = rho
                            s_jikxvl = 0.0
                            for y in range(M):
                                for w in range(B):
                                    s_jikxvl += f[i, j-1, l, y, w]
                            s_jikxv *= s_jikxvl
                            f[i, j, k, x, v] += s_jikxv
                        
                        if r[i, j] == v: f[i, j, k, x, v] *= 1.-(B-1)*eps
                        else: f[i, j, k, x, v] *= eps

                        if v == int2anc[x][k]: f[i, j, k, x, v] *= 1.-(B-1)*mu
                        else: f[i, j, k, x, v] *= mu
                        #f[i, j, k, x, v] /= B**K

    
        ## for x in range(M):
        ##     p_x = 1.0
        ##     for i in range(N):
        ##         s_xi = 0.0
        ##         for k in range(K):
        ##             for v in range(B):
        ##                 s_xi += f[i, j, k, x, v]
        ##         p_x *= s_xi
        ##     ff[j] += p_x

        for i in range(N):
            ff[j] = 1.0
            s_xi = 0.0
            for x in range(M):
                for k in range(K):
                    for v in range(B):
                        s_xi += f[i, j, k, x, v]
            ff[j] *= s_xi
        #ff[j] /= B**(K*L)
    #print 'ff_%d = %f' % (j, ff[L-1])
    return ff[L-1]

def main():
    
    '''What does a main do?
    '''
    import sys
    from Bio import SeqIO
    from itertools import product
    
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
            print count
            for z in all_Z:
                for a in all_A:
                    for h in all_H:
                        count += 1
                        p += full_joint(r, z, a, h, mu, eps, rho)

        
    elif method == 'full_joint_i':
        all_Zi = constellation(K, L, 1)
        all_Hi = constellation(B, L, 1)
        for r in all_R:
            print 'count=', count
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

    print '\n----------> check p =', p
    print '----------> B**(K*L*N)=', B**(K*N)
    #print '\nB**(N*L)=', B**(N*L)
    print 'count =', count
    
if __name__ == '__main__':
    main()
