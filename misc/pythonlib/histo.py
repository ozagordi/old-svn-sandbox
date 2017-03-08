
import sys

def histogram(data, start, stop, mesh):
    """
    Returns the histogram vector for data
    """

    histo = [0]*mesh
    step = float(stop - start)/mesh
    
    for d in data:
        for i in range(mesh):
            if start + i*step < d and d <= start + (i+1)*step:
                histo[i] = histo[i] + 1
    
    histo_data = [ [start + (i+1)*step, histo[i]] for i in range(mesh) ]
    
    return histo_data

def int_histogram(data, start=None, stop=None, title_str='Histogram,'):
    """
    Returns the histogram vector for data, where data only contains integer
    """
    for i in data:
        if type(i) != type(1):
            sys.exit('data must be integer')
    if not start:
        start = min(data)
    if not stop:
        stop = max(data)
    
    histo = [0]*(stop + 1 - start)
    step  = 1
    total = len(data)
    for d in data:
        if start <= d and d <= stop:
            histo[d-start] += 1
    
    histo_data = [ [i+start, histo[i]] for i in range(stop + 1 - start) ]

    try:
        import Gnuplot
        g = Gnuplot.Gnuplot(debug=0)
        g.set(title=title_str + ' total count = %d' % total )
        g('set style data boxes')
        g('set style fill solid 0.3')
        g('set xrange [%f:%f]' % (float(start)-1., float(stop)+1.))
        g.plot(histo_data)
    except:
        print 'did not import Gnuplot'
    return histo_data

