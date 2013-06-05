#!/usr/bin/env python
"""
Script to extract temperature (K) data from a vasp md simulation
"""
import os, sys, commands
import glob

def main():

    # Read in timesteps and temperatures (K)
    TEMP = commands.getoutput("grep 'kinetic Energy' OUTCAR | awk '{print $7}'").split()
    
    # Write data to file
    out = open('temperature.dat','w')
    out.write('# tstep, temp (K)\n')
    for i in range(len(TEMP)):
        out.write(str(i+1)+'    '+str(TEMP[i])+'\n')
    out.close()

if __name__ == '__main__':
    main()
