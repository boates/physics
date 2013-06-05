#!/usr/bin/env python
"""
Script to extract energy (eV) data from a vasp md simulation
"""
import os, sys, commands
import glob

def main():

    # Read in energy (eV)
    E = commands.getoutput("grep 'free  energy' OUTCAR | awk '{print $5}'").split()
    
    # Write data to file
    out = open('energy.dat','w')
    out.write('# energy (eV)\n')
    for i in range(len(E)):
        out.write(str(float(E[i]))+'\n')
    out.close()

    os.system('blocker energy.dat > energy.blocker')

if __name__ == '__main__':
    main()
