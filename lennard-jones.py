#!/usr/bin/env python
"""
lennard-jones.py
Author: Brian Boates

Script to output a data file for a Lennard-Jones potential
hard-sphere potential of the form: V = 4e*[(s/r)^12 - (s/r)^6]
e: depth, s: finite distance to minimum
"""
import os, sys
import Numeric

global RMIN
global RMAX
global STEPSIZE
RMIN = 1.0
RMAX = 12.0
STEPSIZE = 0.01

def main():

    try:
        e = float(sys.argv[1])
        s = float(sys.argv[2])
    except IndexError:
        print '\nPlease provide the potential parameters for LJ-potential: epsilon and sigma\n'
        sys.exit(0)

    r = Numeric.arange(RMIN+STEPSIZE,RMAX+STEPSIZE,STEPSIZE)
    V = 4*e*( (s / r)**12 - (s / r)**6 )

    out = open('LJ-potential.dat','w')
    out.write('# LJ-potential with epsilon = '+str(e)+' and sigma = '+str(s)+'\n')
    out.write('# r\tV(r) = 4*e*( (s / r)**12 - (s / r)**6 )')

    for i in range(len(r)):
        out.write(str(r[i])+'\t'+str(V[i])+'\n')

    out.close()
    
if __name__ == '__main__':
    main()
