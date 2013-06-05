#!/usr/bin/env python
"""
Script to extract energy (eV) data from a vasp md simulation for hugoniot calculations
"""
import os, sys, commands
import glob

def main():

    # Read in "energy without entropy" (eV) and the kinetic energy too
    E = commands.getoutput("grep 'energy  without entropy' OUTCAR | awk '{print $4}'").split()
    K = commands.getoutput("grep EKIN OUTCAR | awk '{print $5}'").split()

    if len(E) != len(K):
        min_len = min([len(E),len(K)])
        E = E[:min_len]
        K = K[:min_len]

    E_hug = []
    for i in range(len(E)):
        E_hug.append(float(E[i])+float(K[i]))
    
    # Write data to file
    out = open('energy_hugoniot.dat','w')
    out.write('# energy (eV)\n')
    for i in range(len(E_hug)):
        out.write(str(E_hug[i])+'\n')
    out.close()

if __name__ == '__main__':
    main()
