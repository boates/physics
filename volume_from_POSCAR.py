#!/usr/bin/env python
"""
volume_from_POSCAR.py
Author: Brian Boates

Calculate cell volume from POSCAR file
"""
import os, sys, commands, glob
import numpy

def main():

    # Retrieve user input
    try:
        pos = open(sys.argv[1],'r')
    except:
        print '\n usage: '+sys.argv[0]+' POSCAR\n'
        sys.exit(0)

    # Read relevent information
    header = pos.readline()
    alat = float(pos.readline().strip())
    ax, ay, az = pos.readline().split()
    bx, by, bz = pos.readline().split()
    cx, cy, cz = pos.readline().split()
    natoms = pos.readline().split()
    try:
        natom = 0
        for i in range(len(natoms)):
            natom += int(natoms[i])
    except:
        natoms = pos.readline().split()
        natom = 0
        for i in range(len(natoms)):
            natom += int(natoms[i])

    # Create lattice vectors
    a = numpy.array([float(ax), float(ay), float(az)]) * alat
    b = numpy.array([float(bx), float(by), float(bz)]) * alat
    c = numpy.array([float(cx), float(cy), float(cz)]) * alat

    # Calculate and print volume in Ang^3 and Ang^3/atom
    V = numpy.dot(a, numpy.cross(b, c))
    v = V / natom

    print V, v


if __name__ == '__main__':
    main()
