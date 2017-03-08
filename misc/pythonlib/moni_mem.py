#!/usr/bin/env python
__author__ = "Osvaldo Zagordi"
__version__ = "$Revision: 0.1 $"
__date__ = "$Date: 2008/06/13$"
__copyright__ = ""
__license__ = "Python"

import sys
import subprocess
import time
told = None
PIPE = subprocess.PIPE
sys.py3kwarning = False

def main(told, inc=1):
    """ Only called if run not interactively
    """
    
    if not sys.argv[1].isdigit():
        cmd = 'ps axc | grep %s' % sys.argv[1]

        process = subprocess.Popen(cmd, shell=True, stdout=PIPE).stdout
        try:
            PID = process.readline().split()[0]
        except IndexError:
            print >> sys.stderr, '# No process name contains:', sys.argv[1]
            print '#\t0'
            time.sleep(inc)
            return
        assert PID.isdigit(), 'PID must be a number'
    else:
        PID = sys.argv[1]
    cmd = 'ps -p %s -o rss -o time' % PID

    process = subprocess.Popen(cmd, shell=True, stdout=PIPE).stdout
    keyline = process.readline().split()
    valueline = process.readline().split()
    
    monitor = {}
    try:
        for i in range(len(keyline)):
            monitor[keyline[i]] = valueline[i]
    except:
        sys.exit('# process not found')
    tnow = monitor['TIME']
    mnow = monitor['RSS']
    if tnow != told:
        told = tnow
        print tnow, mnow, 'K'
    sys.stdout.flush()
    time.sleep(inc)
    
    return told


if __name__ == "__main__":
    while True:
      told =  main(told)
