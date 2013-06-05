#!/usr/bin/env python
"""
calculate the compressibility from PV data
"""
import os, sys, commands, glob
import numpy
from scipy import integrate

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        if '#' in lines[0]: lines.pop(0)
        f.close()
    except:
        print '\n usage: '+sys.argv[0]+'  file_with_rs_and_P_as_first_two_colums.dat.\n'
        sys.exit(0)

    # Read in data and calculate the volumes
    rs, P = [], []
    for line in lines:
        row = line.split()
        rs.append(float(row[0]))
        P.append(float(row[1]))
    V = (numpy.pi*4.0/3.0)*numpy.array(rs)**3

    # Use forward difference method
    out = open('compressibility.dat','w')
    for i in range(len(P)-1):
        dP = P[i+1] - P[i]
        dV = V[i+1] - V[i]
        kT =  -1.0/V[i] * dV/dP 
        out.write(str(rs[i])+'  '+str(kT)+'  '+str(P[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
