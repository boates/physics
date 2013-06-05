#!/usr/bin/env python
"""
DOSCAR_plottable.py
Author: Brian Boates

Read in a DOSCAR file
Write a DOSCAR.dat file for east
plotting (shifted by E_fermi)
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
    except:
        print '\n usage: '+sys.argv[0]+' DOSCAR\n'
        sys.exit(0)

    # Read in header information
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    line = f.readline().split()
    nLines  = int( line[2] )
    E_fermi = float( line[-2] )

    # Read in the DOS and write to file
#    dos = open(sys.argv[1]+'.dat','w')
    dos = open('DOSCAR.dat','w')
    for i in range(nLines):
        row = f.readline().split()
        E    = float( row[0] ) - E_fermi
        DOS  = row[1]
        IDOS = row[2]
        dos.write(str(E)+' '+DOS+' '+IDOS+'\n')
    dos.close()


if __name__ == '__main__':
    main()
