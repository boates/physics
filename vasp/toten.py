#!/usr/bin/env python
"""
Script to extract energy (eV) data from a vasp md simulation
"""
import os, sys, commands
import glob

def main():

    # Retrieve user input
    try:
        outcar = sys.argv[1]
    except:
        print '\n usage: '+sys.argv[0]+' OUTCAR\n'
        sys.exit(0)

    # Read in energy (eV)
    toten = commands.getoutput("grep 'free  energy' OUTCAR | awk '{print $5}'").split()
    
    # Write data to file
    out = open('toten.dat','w')
    for i in range(len(toten)):
        out.write(str(float(toten[i]))+'\n')   # in eV
    out.close()


if __name__ == '__main__':
    main()
