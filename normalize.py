#!/usr/bin/env python
"""
normalize.py
Author: Brian Boates

Read in a data file and do a running avg on a given column
"""
import os, sys, commands, glob
import Numeric

def main():

    try:
        f = open(sys.argv[1],'r')            # Open data file
        lines = f.readlines()
        while '#' in lines[0].split():
            lines.pop(0)
        f.close()
        col = int(sys.argv[2]) - 1           # Index of column to be normalized
    except:
        print '\nusage: '+sys.argv[0]+' fname, column_index\n'
        sys.exit(0)

    # Read in the data
    C = [[] for i in range(len(lines[0].split()))]
    for i in range(len(lines)):
        row = lines[i].split()
        for j in range(len(row)):
            C[j].append(float(row[j]))

    RAW = Numeric.array(C[col])
    NORM = RAW / Numeric.sum(RAW)
    C[col] = list(NORM)

    # Write to the output file
    out = open(sys.argv[1]+'.norm','w')
    for i in range(len(C[0])):

        for j in range(len(C)):

            out.write(str(C[j][i])+' ')

        out.write('\n')

    out.close()
    

if __name__ == '__main__':
    main()
