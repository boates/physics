#!/usr/bin/env python
"""
Read in a data file and do a running avg on a given column
"""
import os, sys, commands, glob

def main():

    try:
        f = open(sys.argv[1],'r')            # Open data file
        lines = f.readlines()
        while '#' in lines[0].split():
            lines.pop(0)
        f.close()
        col = int(sys.argv[2]) - 1           # Index of column to be averaged
        delta = (int(sys.argv[3])/2) * 2     # Number of points to include in centered avg
    except:
        print '\nusage: '+sys.argv[0]+' fname, column_index, number_of_pts_per_avg\n'
        sys.exit(0)

    # Read in the data
    C = [[] for i in range(len(lines[0].split()))]
    for i in range(len(lines)):
        row = lines[i].split()
        for j in range(len(row)):
            C[j].append(float(row[j]))

    # Do the averaging
    AVG = []
    for i in range(delta/2,len(lines)-delta/2):

        s = 0.0

        for j in range(i-delta/2,i+delta/2+1):

            s += C[col][j]
        
        AVG.append(s/float(delta+1))
        
    # Write to the output file
    out = open(sys.argv[1]+'.avg','w')
    for i in range(len(AVG)):

        C[col][i+delta/2] = AVG[i]

        for j in range(len(C)):

            out.write(str(C[j][i+delta/2])+' ')

        out.write('\n')

    out.close()
    

if __name__ == '__main__':
    main()
