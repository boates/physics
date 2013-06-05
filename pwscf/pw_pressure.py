#!/usr/bin/env python
"""
Script to extract pressure (GPa) data from a pwscf md simulation
"""
import os, sys, commands
import glob

def main():

    # Check to see how many output files there are
    if len(glob.glob('output/*.out')) == 1:
        X = -1
    else:
        X = 0

    # Read in timesteps, times (ps), and pressures (GPa)
    TIME = commands.getoutput("grep 'pico-seconds' output/*.out | awk '{print $"+str(4+X)+"}'").split()
    TSTEP = commands.getoutput("grep 'Entering Dynamics' output/*.out | awk '{print $"+str(6+X)+"}'").split()
    PRES = commands.getoutput("grep kbar output/*.out | awk '{print $"+str(7+X)+"/10.0}'").split()

    # Write data to file
    out = open('pressure.dat','w')
#    out.write('# tstep, time (ps), pres (GPa)\n')
    out.write('# pres (GPa)\n')
    for i in range(len(PRES)):
        out.write(str(PRES[i])+'\n')
#        out.write(str(TSTEP[i])+'    '+str(TIME[i])+'    '+str(PRES[i])+'\n')
    out.close()

if __name__ == '__main__':
    main()
    
