#!/usr/bin/env python
"""
force_match_vasp.py
Author: Brian Boates

Perform force-matching analysis on VASP output
to obtain pair-potentials.

Assumptions:
 - Orthorhombic cell
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        f = sys.argv[1]
    except:
        print '\n usage: '+sys.argv[0]+' OUTCAR\n'
        sys.exit(0)

    # Remove any files from a previous analysis
    os.system('rm -rf FORCE_OUT/ position_force TRAJECTORYn1 Input.in fmv.log')

    # Get lattice constants from OUTCAR [Ang]
    a, b, c = commands.getoutput("grep -A1 'length of vectors' "+f+" | tail -2 | grep -v length | awk '{print $1,$2,$3}'").split()
#    a, b, c = float(a), float(b), float(c)
    a, b, c = float(a)/0.5291772, float(b)/0.5291772, float(c)/0.5291772
    
    # Get natom from OUTCAR
    natom = int(commands.getoutput("grep 'number of ions' "+f).split()[-1])

    # Determine half the shortest lattice vector [Ang]
    half = min([a, b, c]) / 2.0

    # Create the position_force file in [Ang], [eV/Ang]
#    os.system("grep -A"+str(natom+1)+" 'TOTAL-FORCE' "+f+" | grep -v TOTAL | grep -v '\-\-' > position_force")
    # in [Bohr], [Ha/Bohr]
    os.system("grep -A"+str(natom+1)+" 'TOTAL-FORCE' "+f+" | grep -v TOTAL | grep -v '\-\-' | awk '{print $1/0.5291772,$2/0.5291772,$3/0.5291772,$4*0.5291772/27.211383,$5*0.5291772/27.211383,$6*0.5291772/27.211383}' > position_force")
    

    # Determine the number of steps
    Nsteps = int(commands.getoutput("grep pressure "+f+" | wc -l"))

    # Create the TRAJECTORYn1 file & remove position_force
    os.system("awk '{if (n % "+str(natom)+" == 0) {print n/"+str(natom)+" + 1;print $0} else print $0;n+=1}' position_force > TRAJECTORYn1")
    os.system("echo "+str(Nsteps + 1)+" >> TRAJECTORYn1")
    os.system('rm -f position_force')

    # Create Input.in file
    os.system('sed s/"ACELL"/"'+str(a)+'"/g /home/boates/bin/TEMPLATE.in | sed s/"BCELL"/"'+str(b)+'"/g | sed s/"CCELL"/"'+str(c)+'"/g > tmp.in')
    os.system('sed s/"NSTEPS"/"'+str(Nsteps)+'"/g tmp.in | sed s/"HALFBOX"/"'+str(half)+'"/g | sed s/"NATOM"/"'+str(natom)+'"/g > Input.in')
    os.system('rm -f tmp.in')

    # Run force-matching
    os.system('m.e > fmv.log')
    

if __name__ == '__main__':
    main()

