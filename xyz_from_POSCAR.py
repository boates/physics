#!/usr/bin/env python
"""
Convert a reduced coordinate POSCAR file to an xyz file in angstroms
"""
import os, sys, commands

def main():

    f = open('POSCAR','r')
    header = f.readline()
    alat = float(f.readline().strip())
    ax = float(f.readline().split()[0]) * alat
    ay = float(f.readline().split()[1]) * alat
    az = float(f.readline().split()[2]) * alat
    natom = int(f.readline().strip())
    f.readline()

    lines = f.readlines()
    out = open('POSCAR.xyz','w')
    out.write(str(natom)+'\n1\n')
    for i in range(natom):
        x = float(lines[i].split()[0]) * ax
        y = float(lines[i].split()[1]) * ay
        z = float(lines[i].split()[2]) * az
        out.write('N '+str(x)+' '+str(y)+' '+str(z)+'\n')

    f.close()
    out.close()

if __name__ == '__main__':
    main()
