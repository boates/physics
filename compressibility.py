#!/usr/bin/env python
"""
compressibility.py
Author: Brian Boates

Calculate Kt = (-1/V)*dV/dP for a file with
column 1 = V, column 2 = P

Uses central difference method for all points less the
first and final points where the forward and backward
difference methods are used.
"""
import os, sys, commands, glob
import numpy

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        if '#' in lines[0]: lines.pop(0)
        f.close()
    except:
        print '\n usage: '+sys.argv[0]+' VP.dat\n'
        sys.exit(0)

    # Read in data
    V, P = [], []
    for line in lines:
        row = line.split()
        V.append(float(row[0]))
        P.append(float(row[1]))

    # Differentiate
    out = open('compressibility.dat','w')
    for i in range(len(V)):
        if i == 0:
            dV = V[i+1] - V[i]
            dP = P[i+1] - P[i]
        elif i == len(V)-1:
            dV = V[i] - V[i-1]
            dP = P[i] - P[i-1]
        else:
            dV = V[i+1] - V[i-1]
            dP = P[i+1] - P[i-1]
        out.write(str(P[i])+'  '+str( (-1/V[i]) * (dV/dP) )+'\n')
    out.close()


if __name__ == '__main__':
    main()
