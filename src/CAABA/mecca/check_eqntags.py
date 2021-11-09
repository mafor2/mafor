#!/usr/bin/env python

# authors:
# original code: Martin Schultz, Juelich (2012)
# modified:      Rolf Sander, Mainz (2012-2017)

import sys
import os

if len(sys.argv) > 1:
    eqnfilename = sys.argv[1]
else:
    eqnfilename = 'gas.eqn'

print '%s checks for duplicate equation tags ("reaction numbers") in %s:' % \
    (os.path.basename(sys.argv[0]), eqnfilename)
EQNFILE = open(eqnfilename, 'r')
eqntag_list = []
errorstring = ''
for lineno,line in enumerate(EQNFILE):
    # remove leading spaces:
    line = line.lstrip(' ')
    # test for equation tags like <G1000>:
    if line[0] == '<':
        eqntag = line[0:line.index('>')+1]
        if eqntag in eqntag_list:
            errorstring = '%s  %s at lines %s and %s\n' % \
                (errorstring, eqntag, eqntag_list.index(eqntag)+1, lineno+1)
        else:
            eqntag_list.append(eqntag)
    else:
        # synchronize lineno and eqntag_list:
        eqntag_list.append('DUMMY') 
EQNFILE.close()

if (errorstring == ''):
    print 'OK: No duplicate equation tags in %s' % (eqnfilename)
else:
    print 'ERROR: Duplicate equation tags in %s' % (eqnfilename)
    print errorstring,
    sys.exit(1) # exit with error status 1
