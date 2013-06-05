#!/usr/bin/env python
"""
Script to extract energy (Ry) data from a pwscf md simulation
"""
import os, sys, commands
import glob

def main():

    # Check to see how many output files there are (for awk column purposes)
    if len(glob.glob('output/*.out')) == 1:
        X = -1
    else:
        X = 0

    # Read in timesteps, times (ps), and energies (Ry)
    TIME = commands.getoutput("grep 'pico-seconds' output/*.out | awk '{print $"+str(4+X)+"}'").split()
    TSTEP = commands.getoutput("grep 'Entering Dynamics' output/*.out | awk '{print $"+str(6+X)+"}'").split()
    ENER = commands.getoutput("grep 'Ekin + Etot (const)' output/*.out | awk '{print $"+str(7+X)+"}'").split()
    
    # Write data to file
    out = open('energy.dat','w')
    out.write('#energy (Ry)\n')
    for i in range(len(ENER)):
        out.write(str(ENER[i])+'\n')
#        out.write(str(TSTEP[i])+'    '+str(TIME[i])+'    '+str(ENER[i])+'\n')
    out.close()

if __name__ == '__main__':
    main()
    
