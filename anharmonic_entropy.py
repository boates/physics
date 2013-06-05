#!/usr/bin/env python
"""
anharmonic_entropy.py
Author: Brian Boates

Read in <U(T,V)> - Umin(T=0,V)
at constant volume in eV/atom

<U(T,V)> is the internal energy MD average
as a function of temperature done on a grid
of different temperatures at constant volume

Umin(T=0,V) is the internal energy of the
relaxed zero-temperature structure

[see Bonev's notes & Lacks & Shukla, PRB 54, 3266, 1996]
"""
import os, sys, commands, glob, numpy
from scipy import integrate

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        f.close()

        Emin = float(sys.argv[2])
        
        kB = 8.6173423E-05 # eV/K
        
    except:
        print '\n usage: '+sys.argv[0]+' E.dat(eV/atom) Emin(eV/atom)\n'
        sys.exit(0)

    # Read in data, perform calculation, and write to output
    out = open('S_A.dat','w')
    T_list, dE_list = [], []
    for i in range(len(lines)):

        # Read in T & E, Emin = -6.36578 for rs=1.22 128 atom Cmcm gamma-point w/ T=0 shift correction
        T_list.append( float(lines[i].split()[0]) )
        dE_list.append( float(lines[i].split()[1]) - Emin - 3.0/2.0*kB*T_list[i] )

        T = numpy.array(T_list)
        dE = numpy.array(dE_list)
        
        # Integrate dE/T^2
        integral = integrate.trapz( dE/T**2, T )

        # dF_cl / dT
        dF_cl = - (integral + dE[i] / T[i])

        # Anharmonic entropy is ~ -dF_cl / dT
        S_A = - dF_cl

        # Write to file
        out.write(str(T[i])+' '+str(S_A/kB)+'\n')

    out.close()


if __name__ == '__main__':
    main()
