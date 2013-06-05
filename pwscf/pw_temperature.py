#!/usr/bin/env python
"""
Script to extract temperature (K) data from a pwscf md simulation
"""
import os, sys, commands
import glob

def main():

    # Check to see how many output files there are
    if len(glob.glob('output/*.out')) == 1:
        X = -1
    else:
        X = 0

    # Read in timesteps, times (ps), and temperatures (K)
    TIME = commands.getoutput("grep 'pico-seconds' output/*.out | awk '{print $"+str(4+X)+"}'").split()
    TSTEP = commands.getoutput("grep 'Entering Dynamics' output/*.out | awk '{print $"+str(6+X)+"}'").split()
    TEMP = commands.getoutput("grep temperature output/*.out | grep -v Starting | awk '{print $"+str(4+X)+"}'").split()
    
    # Write data to file
    out = open('temperature.dat','w')
#    out.write('# tstep, time (ps), temp (K)\n')
    out.write('# temp (K)\n')
    for i in range(len(TEMP)):
        out.write(str(TEMP[i])+'\n')
#        out.write(str(TSTEP[i])+'    '+str(TIME[i])+'    '+str(TEMP[i])+'\n')
    out.close()

if __name__ == '__main__':
    main()
    
