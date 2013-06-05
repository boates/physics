#!/usr/bin/env python
"""
band_structure_vasp.py
Author: Brian Boates

Extract the band structure from VASP EIGENVAL file
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        Ef = float(sys.argv[2])
    except:
        print '\n usage: '+sys.argv[0]+' EIGENVAL Efermi\n'
        sys.exit(0)

    # Read in header information
    natom = f.readline().split()[0]
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    not_nband, nkpt, nband = f.readline().split()

    # Read in the band energies
    out = open('band_structure.dat','w')
    for k in range(int(nkpt)):
        f.readline()
        kx, ky, kz, w = f.readline().split()
        out.write(kx+' '+ky+' '+kz)
        for b in range(int(nband)):
            band = float(f.readline().split()[-1]) - Ef
            out.write(' '+str(band))
        out.write('\n')
    out.close()


if __name__ == '__main__':
    main()
