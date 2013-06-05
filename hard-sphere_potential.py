#!/usr/bin/env python
"""
hard-sphere_potential.py
Author: Brian Boates

Script to output a data file for an analytic repulsive
hard-sphere potential of the form: V = A/r^12
A is in energy*distance^12
"""
import os, sys
import Numeric

global RMIN
global RMAX
global STEPSIZE
RMIN = 3.0
RMAX = 12.0
STEPSIZE = 0.01

def main():

    try:
        A = float(sys.argv[1])
    except IndexError:
        print '\nPlease provide the potential parameter (A for A/r^12)\n'
        sys.exit(0)

    r = Numeric.arange(RMIN+STEPSIZE,RMAX+STEPSIZE,STEPSIZE)
    V = A / r**12

    out = open('hard-sphere_A'+str(A)+'.dat','w')
    out.write('# Hard-sphere repulsive potential with A = '+str(A)+'\n')
    out.write('# r\tV(r) = A/r^12\n')

    for i in range(len(r)):
        out.write(str(r[i])+'\t'+str(V[i])+'\n')

    out.close
    
if __name__ == '__main__':
    main()
