#!/usr/bin/env python
"""
Script to convert a DL_POLY CONFIG file (with only coords) to an xyz file
"""
import os, sys

def main():

    try:
        f = open('CONFIG','r')
        # Remove header lines and lattice vectors
        for i in range(5): f.readline()
    except:
        print '\nMake sure the CONFIG file is present\n'
        sys.exit(0)

    # Open the output xyz file
    out = open('CONFIG.xyz','w')

    # Read in rest of CONFIG file
    lines = f.readlines()
    f.close()

    # Get the number of atoms
    NATOM = len(lines) / 2

    # Write the "header" of the xyz file
    out.write(str(NATOM)+'\n1\n')

    # Loop through the rest of the CONFIG file and write to the xyz file
    for i in range(0,len(lines),2):
        TYPAT = lines[i].split()[0]
        COORD = lines[i+1]
        out.write(TYPAT+COORD)

    out.close()
    
if __name__ == '__main__':
    main()
