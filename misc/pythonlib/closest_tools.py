def find_closest(hr):
    '''
    The diff_thresh has been set to 0.025 because even when aligning error-free reads
    to the original haplotypes, the distribution of differences of the best 2
    identities goes from 0.023 to 0.091 (~9%)
    '''
    from pythonlib import Alignment
    import tempfile
    import subprocess
    import heapq
    import operator
    
    diff_thresh = 0.0125
    abs_thresh = 0.1

    # ref_file = './ref.fas'
    out = tempfile.NamedTemporaryFile()
    outname = out.name
    hap, ref_file = hr
    cmline = 'needle -asequence asis:\'%s\' -bsequence %s \
              -gapopen 10.0 -gapextend 1.0 -auto -adesshow3 -out %s -aformat3 markx10' \
        % (hap, ref_file, outname)
    subprocess.call(cmline, shell=True)
    d = Alignment.alignfile2dict([outname], 'n', 6.0, 3.0, Verbose = False)['asis']
    out.close()
    this = {}
    mm = {}
    gaps = {}
    for k, v in d.items():
        v.summary()
        this[v.id_b] = float(v.mismatch)/(v.stop - v.start + 1) #float(v.ident)/(v.stop - v.start + 1)
        mm[v.id_b] = v.mismatch # v.stop - v.start + 1 - v.ident
        gaps[v.id_b] = v.int_gaps
        
    best2 = heapq.nsmallest(2, this.items(), operator.itemgetter(1))
    rel_delta = (best2[1][1] - best2[0][1])#/best2[0][1]
    
    if  rel_delta >= diff_thresh and best2[0][1] <= abs_thresh:
        return best2[0][0], best2[0][1], mm[best2[0][0]], gaps[best2[0][0]]
    else:
        return None, gaps[best2[0][0]]
    
