#!/usr/bin/env python
"""
Get the phonon contribution to the free energy from anaddb
combine this with the enthalpy calculated from the abinit
OUTPUT file for an adjusted free energy.
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        N = int(sys.argv[1])
        anaddb = sys.argv[2]
        OUTPUT = sys.argv[3]
    except:
        print '\n usage: '+sys.argv[0]+'  natom_primitive_cell  anaddb.out  OUTPUT\n'
        sys.exit(0)

    factor = str( 6.2415097E18 / 6.0221415E23 / N )  # J/mol-c to eV/atom
        
    # Grab phonon Helmholtz free energies (assume 51 temperatures done)
    os.system("grep -A51 Avogadro "+anaddb+" | tail -51 | awk '{print $1,$2*"+factor+"}' > F_phonon.dat")
    os.system("grep -A51 Avogadro "+anaddb+" | tail -51 | awk '{print $1,$3*"+factor+"}' > E_phonon.dat")

    e = float(commands.getoutput("grep 'Total energy' "+OUTPUT+" | awk '{print $3}'"))       # eV
    p = float(commands.getoutput("grep Pressure "+OUTPUT+" | awk '{print $8}'"))             # GPa
    v = float(commands.getoutput("grep volume "+OUTPUT+" | awk '{print $5}'").split()[-1])   # bohr^3

    h = e + p*v*0.00092489631   # eV
    h = h / N

    format = "\"%4.1f %4.6f \\n\""
    os.system("awk '{printf "+format+",$1,$2+"+str(h)+"}' F_phonon.dat > G.dat")

if __name__ == '__main__':
    main()
