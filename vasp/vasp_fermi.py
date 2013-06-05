#!/usr/bin/env python
"""
Script to extract fermi energies (eV) data from a vasp md simulation
"""
import os, sys, commands
import glob

def main():

    # Read in timesteps and energies (eV)
    FERMI = commands.getoutput("grep -A100 Iteration OUTCAR | grep -B50 '(   1)' | grep E-fermi | awk '{print $3}'").split()
    
    # Write data to file
    out = open('fermi.dat','w')
    out.write('# tstep, fermi energy (eV)\n')
    for i in range(len(FERMI)):
        out.write(str(i)+'    '+str(float(FERMI[i]))+'\n')
    out.close()

if __name__ == '__main__':
    main()
