#!/usr/bin/env python
"""
Script to extract pressure (GPa) data from a vasp md simulation

The Kineteic part MUST BE and IS added in by hand in this script
"""
import os, sys, glob, commands

def main():

    # Read data from OUTCAR
    kB = 8.6173423E-05 # eV/K
    factor = 160.21765 # eV/Ang^3 to GPa
    N_over_V = 1.0 / float(commands.getoutput("grep 'volume\/ion' OUTCAR").split()[-2]) # atom/Ang^3
    pressures = commands.getoutput("grep pressure OUTCAR | awk '{print $4/10.0}'").split() # GPa
    temps = commands.getoutput("grep 'kinetic Energy' OUTCAR | sed s/'ture'/' '/g | awk '{print $7}'").split() # K

    # Add the kinetic part
    total_p = []
    for i in range(min(len(pressures),len(temps))):
        P = float(pressures[i])
        T = float(temps[i])
        kin_p = N_over_V*kB*T*factor
        total_p.append( P + kin_p )
#        print P, kin_p, total_p[i]
    
    # Write data to file
    out = open('pressure.dat','w')
    for i in range(len(total_p)):
        out.write(str(str(i+1))+' '+str(total_p[i])+'\n')
    out.close()

if __name__ == '__main__':
    main()
