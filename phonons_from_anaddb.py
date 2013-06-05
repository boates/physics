#!/usr/bin/env python
"""
Get phonons from anaddb output file, print as a row of data
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        anaddb = sys.argv[1]
        natom = int(sys.argv[2])
    except:
        print '\n usage: '+sys.argv[0]+' anaddb.out natom_unit_cell\n'
        sys.exit(0)

    nmodes = 3*natom
    freqs  = commands.getoutput("grep -A10 'Phonon freq' "+anaddb+" | grep -B10 Eigendisplacements | grep -v e | awk '{print $2,$3,$4,$5,$6}'").split()
    try:
        nph1l = int(commands.getoutput("grep nph1l "+anaddb).split()[-1])
    except:
        nph1l = 1
 
    if len(freqs)/nph1l != nmodes:
        print '\n WARNING: Number of frequencies found not consistent with 3*natom_unit_cell --- exiting... \n'
        sys.exit(0)

    for i in range(nph1l):
        print
        for j in range(nmodes):
            print freqs[i*nmodes+j],


if __name__ == '__main__':
    main()
