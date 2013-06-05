#!/usr/bin/env python
"""
generate_R-3c.py
Author: Brian Boates

Return a POSCAR of molecules about the R-3c
molecular centers. Oreintations given in
Lesar, JCP, 1984.
"""
import os, sys, commands, glob, numpy

def main():

    # Create the lattice vectors
    alat = 4.0
    ax, ay, az =   1.0,  0.00000000,  0.80200732
    bx, by, bz =  -0.5,  0.86602540,  0.80200732
    cx, cy, cz =  -0.5, -0.86602540,  0.80200732

    # Determine the cell volume per atom (atom NOT center)
    a_vec = numpy.array([ax, ay, az])
    b_vec = numpy.array([bx, by, bz])
    c_vec = numpy.array([cx, cy, cz])
    V = numpy.dot(a_vec, numpy.cross(b_vec, c_vec))
    vol = V / 16.0

    # Create output POSCAR and write header
    out = open('POSCAR.R-3c','w')
    out.write('R-3c, vol='+str(vol)+' ang^3/atom\n')
    out.write('    '+str(alat)+'\n')
    out.write('     '+str(ax)+'  '+str(ay)+'  '+str(az)+'\n')
    out.write('     '+str(bx)+'  '+str(by)+'  '+str(bz)+'\n')
    out.write('     '+str(cx)+'  '+str(cy)+'  '+str(cz)+'\n')
    out.write('  16\n')
    out.write('Direct\n')

    # Scale lattice vector components
    ax, ay, az = ax*alat, ay*alat, az*alat
    bx, by, bz = bx*alat, by*alat, bz*alat
    cx, cy, cz = cx*alat, cy*alat, cz*alat    

    # Define the lattice space and get the inverse
    space = numpy.matrix([ [ax,ay,az], [bx,by,bz], [cx,cy,cz] ])
    inverse = space.I

    # Define the molecular centers in reduced coordinates
    a1 = [0.000, 0.000, 0.000]
    a2 = [0.500, 0.500, 0.500]
    a3 = [0.009, 0.491, 0.250]
    a4 = [0.991, 0.509, 0.750]
    a5 = [0.491, 0.250, 0.009]
    a6 = [0.509, 0.750, 0.991]
    a7 = [0.250, 0.009, 0.491]
    a8 = [0.750, 0.991, 0.509]
    centers = [a1, a2, a3, a4, a5, a6, a7, a8]

    # Create displacement vectors for molecular ions in cartesian coordinates
    bond_length = 1.10
    shift1 = (bond_length/2.0) * numpy.array([1.0/3.0**0.5, 1.0/3.0**0.5, 1.0/3.0**0.5])
    shift2 = (bond_length/2.0) * numpy.array([1.0/3.0**0.5, 1.0/3.0**0.5, 1.0/3.0**0.5])
    shift3 = (bond_length/2.0) * numpy.array([2.0/3.0, 2.0/3.0, -1.0/3.0])
    shift4 = (bond_length/2.0) * numpy.array([2.0/3.0, 2.0/3.0, -1.0/3.0])
    shift5 = (bond_length/2.0) * numpy.array([2.0/3.0, -1.0/3.0, 2.0/3.0])
    shift6 = (bond_length/2.0) * numpy.array([2.0/3.0, -1.0/3.0, 2.0/3.0])
    shift7 = (bond_length/2.0) * numpy.array([-1.0/3.0, 2.0/3.0, 2.0/3.0])
    shift8 = (bond_length/2.0) * numpy.array([-1.0/3.0, 2.0/3.0, 2.0/3.0])

    print shift1

    # Invert coordinates to lattice vector space
    shift1 = shift1   *inverse
    shift2 = shift2   *inverse
    shift3 = shift3   *inverse
    shift4 = shift4   *inverse
    shift5 = shift5   *inverse
    shift6 = shift6   *inverse
    shift7 = shift7   *inverse
    shift8 = shift8   *inverse

    # Reformat the inverted vectors
    s1 = [shift1.item(0),shift1.item(1),shift1.item(2)]
    s2 = [shift2.item(0),shift2.item(1),shift2.item(2)]
    s3 = [shift3.item(0),shift3.item(1),shift3.item(2)]
    s4 = [shift4.item(0),shift4.item(1),shift4.item(2)]
    s5 = [shift5.item(0),shift5.item(1),shift5.item(2)]
    s6 = [shift6.item(0),shift6.item(1),shift6.item(2)]
    s7 = [shift7.item(0),shift7.item(1),shift7.item(2)]
    s8 = [shift8.item(0),shift8.item(1),shift8.item(2)]
    shifts = [s1, s2, s3, s4, s5, s6, s7, s8]

    print s1
    s1 = numpy.array(s1)
    print (s1*space).item(0),(s1*space).item(1),(s1*space).item(2)

    # Create and write molecular coordinates from centers
    for i in range(len(centers)):
        out.write('  '+str(float(centers[i][0]) + shifts[i][0]))
        out.write('  '+str(float(centers[i][1]) + shifts[i][1]))
        out.write('  '+str(float(centers[i][2]) + shifts[i][2]))
        out.write('\n')
        out.write('  '+str(float(centers[i][0]) - shifts[i][0]))
        out.write('  '+str(float(centers[i][1]) - shifts[i][1]))
        out.write('  '+str(float(centers[i][2]) - shifts[i][2]))
        out.write('\n')

    out.close()
        

if __name__ == '__main__':
    main()
