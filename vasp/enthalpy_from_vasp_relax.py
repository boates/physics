#!/usr/bin/env python
"""
Calculate the relaxed enthalpy from vasp relaxation
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        f = open('OUTCAR','r')
    except:
        print '\n Make sure OUTCAR file is present. \n'
        sys.exit(0)

    # Grep for relaxed parameters
    U = float(commands.getoutput("tail -1000 OUTCAR | grep 'free  energy' | tail -n-1 | awk '{print $5}'"))        # in eV
    P = float(commands.getoutput("tail -1000 OUTCAR | grep pressure | tail -n-1 | awk '{print $4/10.0+$9/10.0}'")) # in GPa
    V = float(commands.getoutput("tail -1000 OUTCAR | grep volume | tail -n-1 | awk '{print $5}'"))                # in ang^3
    natom = int(commands.getoutput("head -10000 OUTCAR | grep NIONS").split()[-1])

    # Calculate the enthalpy
    H = U + P*V*0.0062415097   # in eV
    H = H / natom

    print 'H =', H, 'eV/atom'


if __name__ == '__main__':
    main()
