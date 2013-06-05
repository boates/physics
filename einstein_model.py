#!/usr/bin/env python
"""
Calculate the difference in free energy contributions given
phonon frequencies and temperature
"""
import os, sys, commands, glob, math
from math import log, sinh

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1])
        T = int(sys.argv[2])
    except:
        print '\n usage: '+sys.argv[0]+'  einstein.in  T(K)\n'
        sys.exit(0)

    # Define constants
    c  = 2.997925E+10   # cm/s
    kB = 8.617E-05      # eV/K
    hbar = 6.5822E-16   # eV*s

    # Read in the frequencies
    w1_list = f.readline().split()
    w2_list = f.readline().split()

    # Determine the number of modes
    n1 = float(len(w1_list))
    n2 = float(len(w2_list))

    # Convert from cm^-1 to s^-1
    w1_prod, w2_prod = 1.0, 1.0
    for i in range(len(w1_list)):
        w1_prod *= float(w1_list[i])*c
    for i in range(len(w2_list)):
        w2_prod *= float(w2_list[i])*c

    # Calculate free energy differences
    t1 = ( hbar/(kB*T) * w1_prod )**(1/n1)
    t2 = ( hbar/(kB*T) * w2_prod )**(1/n2)

    f1 = kB*T * log( t1 )
    f2 = kB*T * log( t2 )
    df = kB*T * log( t1 / t2 )

    # Print the results
    print 'df =', df*1000, 'meV'

    
if __name__ == '__main__':
    main()
