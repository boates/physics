#!/usr/bin/env python
"""
Script to convert a VASP POSCAR file into a DL_POLY CONFIG file
"""
import os, sys, commands, glob

global TYPAT
TYPAT = 'NN'

def main():

    # Check for required POSCAR file
    try:
        f = open('POSCAR','r')
        lines = f.readlines()
        f.close()
    except:
        print '\n Make sure POSCAR file is present.\n'
        sys.exit(0)

    # Header info from POSCAR file
    lines.pop(0)
    alat = float(lines.pop(0).strip())
    ax, ay, az = lines.pop(0).split()
    bx, by, bz = lines.pop(0).split()
    cx, cy, cz = lines.pop(0).split()
    ax, ay, az, bx, by, bz, cx, cy, cz = float(ax), float(ay), float(az), float(bx), float(by), float(bz), float(cx), float(cy), float(cz)
    natom = lines.pop(0).strip()
    lines.pop(0)

    # Read in POSCAR coordinates
    X, Y, Z = [], [], []
    for i in range(int(natom)):
        row = lines[i].split()
        X.append(float(row[0]))
        Y.append(float(row[1]))
        Z.append(float(row[2]))

    # Write CONFIG file header
    out = open('CONFIG','w')
    out.write('DL_POLY CONFIG: typat='+TYPAT+' natom='+natom+'\n')
    out.write('      0        1\n')
    out.write('       '+str(alat*ax)+'        '+str(alat*ay)+'        '+str(alat*az)+'\n')
    out.write('       '+str(alat*bx)+'        '+str(alat*by)+'        '+str(alat*bz)+'\n')
    out.write('       '+str(alat*cx)+'        '+str(alat*cy)+'        '+str(alat*cz)+'\n')

    # Write coordinates to CONFIG file
    for i in range(int(natom)):
        out.write(TYPAT+'          '+str(i+1)+'\n')
        x = X[i]*ax*alat + Y[i]*bx*alat + Z[i]*cx*alat
        y = X[i]*ay*alat + Y[i]*by*alat + Z[i]*cy*alat
        z = X[i]*az*alat + Y[i]*bz*alat + Z[i]*cz*alat
        out.write('      '+str(x)+'       '+str(y)+'         '+str(z)+'\n')

    out.close()

    # Reformat CONFIG file properly
    os.system('config_spacer.py')
    os.system('mv -f CONFIG_fixed CONFIG')


if __name__ == '__main__':
    main()
