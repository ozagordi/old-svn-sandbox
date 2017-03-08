def Simpson(data):
    ''' Computes Simpson's index (Simpson. Measurement of diversity.
    Nature (1949) vol. 163 (4148) pp. 688) returns index and variance
    '''
    from math import sqrt
    N = sum(data)
    f = [float(i)/N for i in data]
    SI = sum([i*i for i in f])
    sc = sum([i*i*i for i in f])
    # var = 4*(sc - SI*SI)/N
    var_prec = (4*N*(N-1)*(N-2)*sc + 2*N*(N-1)*SI - 2*N*(N-1)*(2*N-3)*SI*SI ) / (N*N*(N-1)*(N-1))
    return SI, sqrt(var_prec)
