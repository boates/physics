#!/usr/bin/env python
"""
Script to extract all thermodynamic variables from a VASP
MD simulation OSZICAR & OUTCAR files

toten.dat:        free energy of the electrons     [eV]
energy.dat:       toten + kinetic energy of ions   [eV]
temperature.dat:  temperature                      [K]
pressure.dat:     pressure + kinetic part          [GPa]
"""
import os, sys, glob, commands

def main():

    # Retrieve user input
    try:
        oszi = sys.argv[1]
        out  = sys.argv[2]
    except:
        print '\n usage: '+sys.argv[0]+' OSZICAR OUTCAR\n'
        sys.exit(0)

    # Extract toten in eV
    os.system("grep F "+oszi+" | awk '{print $7/1.0}' > toten.dat")
    
    # Extract total energy in eV
    os.system("grep F "+oszi+" | awk '{print $7+$11}' > energy.dat")

    # Extract temperature in K
    os.system("grep F "+oszi+" | awk '{print $3/1.0}' > temperature.dat")
    temperatures = commands.getoutput("grep F "+oszi+" | awk '{print $3/1.0}'").split()

    # Extract pressure in GPa
    kB = 8.6173423E-05  # eV/K
    factor = 160.21765  # eV/Ang^3 to GPa
    
    N_over_V = 1.0 / float(commands.getoutput("grep 'volume\/ion' OUTCAR").split()[-2])    # atom/Ang^3
    pressures = commands.getoutput("grep pressure OUTCAR | awk '{print $4/10.0}'").split() # GPa

    # Add the kinetic part
    total_P = []
    for i in range(min(len(pressures),len(temperatures))):
        P = float(pressures[i])
        T = float(temperatures[i])
        kin_P = N_over_V*kB*T*factor
        total_P.append( P + kin_P )
    
    # Write data to file
    out = open('pressure.dat','w')
    for i in range(len(total_P)):
        out.write(str(total_P[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
