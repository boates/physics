#!/usr/bin/env python
"""
van_der_waals.py
Author: Brian Boates

Compute the Van der Waals eq'n of state (P-V)
"""

import os, sys
import Numeric

def main():

    try:
        typat = sys.argv[1]
        natom = int(sys.argv[2])
        temp = float(sys.argv[3])
    except IndexError:
        print '\nusage: ./van_der_waals.py typat(i.e. "N") natom temp(K)\n'
        sys.exit(0)

    if typat == 'N':
        a = 4.0e-49   # J m^3  __ Found in Shroeder stat mech book (pg 181)
        b = 6.0e-29   # m^3    _/
    else:
        print '\nSorry, only info for nitrogen (N) is avaiable... exiting\n'
        sys.exit(0)

    # Define constants
    k = 1.381e-23   # in J/K

    # Create volume array (in m^3)
    vol = Numeric.arange(10,1000.1,0.1) / 6.7483346e+30

    pres = (natom*k*temp)/(vol-natom*b) - (a*natom**2)/vol**2

    # Convert P to GPa and vol to bohr^3
    pres *= -1e-09
    vol *= 6.7483346e+30
 
    out = open('vdw-'+typat+'-'+str(natom)+'.dat','w')
    out.write('# typat='+typat+' natom='+str(natom)+' temp='+str(temp)+'K\n')
    out.write('# a='+str(a)+' , b='+str(b)+'\n')
    out.write('# P (GPa)\tV (bohr^3)\n')
    for i in range(len(pres)):
        out.write(str(pres[i])+'\t'+str(vol[i])+'\n')

    out.close()

if __name__ == '__main__':
    main()
