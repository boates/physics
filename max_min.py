#!/usr/bin/env python
"""
max_min.py
Author: Brian Boates

Return the maximum and minimum values
from a data file for a given column.
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        c = int(sys.argv[2])
        p = int(sys.argv[3])
    except:
        print '\n usage: '+sys.argv[0]+' file.dat extrema_col(i.e. 1,2,...) print_col(i.e. 1,2,...)\n'
        sys.exit(0)

    # Read in the data column
    lines = f.readlines()
    f.close()
    cdat, pdat = [], []
    for line in lines:
        cdat.append( float(line.split()[c-1]) )
        pdat.append( float(line.split()[p-1]) )

    print 'MAX:', pdat[cdat.index(max(cdat))], '  MIN:', pdat[cdat.index(min(cdat))]


if __name__ == '__main__':
    main()
