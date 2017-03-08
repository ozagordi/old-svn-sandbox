#!/usr/bin/env python

import sys

file_1, file_2 = sys.argv[1:]

for l in open(file_1):
    count, hap = l.split()
    for i in range(int(count)):
        print '1', ' '.join(hap)

for l in open(file_2):
    count, hap = l.split()
    for i in range(int(count)):
        print '2', ' '.join(hap)
