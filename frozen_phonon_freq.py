#!/usr/bin/env python
"""
Calculate frequency (in cm^-1) of a mode based on
energy difference and atomic displacement
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        U     = float(sys.argv[1])
        U_dx  = float(sys.argv[2])
        dU    = U_dx - U
        dx    = float(sys.argv[3])
    except:
        print '\n usage: '+sys.argv[0]+' U(x) U(x+dx) dx   (in eV, eV, angstroms, respectively)\n'
        sys.exit(0)

    # Calculate 'spring constant'
    k = 2 / dx**2 * dU   # eV / angstrom^2
    k = k * 16.021765    # J / m^2
    
    # Calculate the mass (assuming nitrogen)
    amu = 14.0067         # g / mol
    avo = 6.0221415E+23   # atoms /mol

    m = (amu / avo) * (1.0/1000.0)   # kg

    m = m**2/(2.0*m)

    w = ( k / m )**0.5   # Hz
    w = w / 2.41796E+14 * 8.0655E+3   # cm^-1

    print 'w =', w
    

if __name__ == '__main__':
    main()
