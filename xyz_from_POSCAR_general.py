#!/usr/bin/env python
"""
Convert a reduced coordinate POSCAR file to an xyz file in angstroms
"""
import os, sys, commands

def main():

    try:
        poscar = sys.argv[1]
    except:
        print '\n usage: '+sys.argv[0]+' POSCAR\n'
        sys.exit(0)

    f = open(poscar,'r')
    header = f.readline()
    alat = float(f.readline().strip())
    a = f.readline().split()
    ax, ay, az = float(a[0])*alat, float(a[1])*alat, float(a[2])*alat
    b = f.readline().split()
    bx, by, bz = float(b[0])*alat, float(b[1])*alat, float(b[2])*alat
    c = f.readline().split()
    cx, cy, cz = float(c[0])*alat, float(c[1])*alat, float(c[2])*alat
    natoms = f.readline().split()
    natom_total = 0
    for natom in natoms: natom_total += int(natom)
    f.readline()
    atoms = ['N','C','O','Si','H','B','He','Li','Be','F','Ne']

    lines = f.readlines()
    out = open('POSCAR.xyz','w')
    out.write(str(natom_total)+'\n1\n')
    k = -1
    for i in range(len(natoms)):
        for j in range(int(natoms[i])):
            k += 1
            line = lines[k].split()
            x = float(line[0])*ax + float(line[1])*bx + float(line[2])*cx
            y = float(line[0])*ay + float(line[1])*by + float(line[2])*cy
            z = float(line[0])*az + float(line[1])*bz + float(line[2])*cz
            out.write(atoms[i]+' '+str(x)+' '+str(y)+' '+str(z)+'\n')

    f.close()
    out.close()

if __name__ == '__main__':
    main()
