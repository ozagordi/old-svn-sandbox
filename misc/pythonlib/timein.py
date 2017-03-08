#!/usr/bin/env python
import sys
import datetime
try:
    delta = int(sys.argv[1])
except:
    delta = 0
h = open('/Users/ozagordi/Desktop/timein.txt', 'a+')
d = datetime.timedelta(minutes=delta)
timein = datetime.datetime.now() + d
date = '%4.4d%2.2d%2.2d' % (timein.year, timein.month, timein.day)
time = '%2.2d:%2.2d' % (timein.hour, timein.minute)
h.write('%s\t%s\n' % (date, time))
