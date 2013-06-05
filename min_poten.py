#!/usr/bin/env python
"""
min_poten.py
Author: Brian Boates

Find the '1st minimum' in pair potential
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        f.close()
    except:
        print '\n usage: '+sys.argv[0]+'  Poten_XX.dat\n'
        sys.exit(0)

    # Get the PP data
    r, v = [], []
    for line in lines:
        row = line.split()
        r.append( float(row[0]) )
        v.append( float(row[1]) )

    # Find the 1st minimum
    vMin = v[0]+1
    i = 0
    while v[i] <= vMin:
        vMin = v[i]
        i += 1
    rMin = r[i-1]

    # Find the barrier height
    barrier = max(v[i-1:]) - vMin
    
    # Print the result: rMin, vMin(eV), vMin(K), barrier(eV), barrier(K)
    out = open('min_poten.dat','w')
    out.write(str(rMin)+' '+str(vMin)+' '+str(vMin/8.6173423E-5)+' '+str(barrier)+' '+str(barrier/8.6173423E-5)+'\n')
    out.close()


if __name__ == '__main__':
    main()
