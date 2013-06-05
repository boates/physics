#!/usr/bin/env python
"""
Author: Brian Boates
Date: May 5, 2009

Generate a random cell of random shape for given rs and natom
"""
import os, sys, commands, glob
import optparse, numpy
from math import pi, sin, cos
from random import random

def volume_from_rs(rs,Nel):
    """
    Return cell volume in angstroms^3.
    rs: density parameter for system
    Nel: number of valence electrons in system
    """
    a0 = 0.5291772   # Bohr radius (angstroms/bohr)
    volume = (4.0*pi/3.0)*Nel * (rs*a0)**3

    return volume


def init_parser():

    usage = 'usage: %prog [options]'
    version="%prog 1.0"
    parser = optparse.OptionParser(usage=usage,version=version)

    parser.add_option('-n', '--natom', type='int', dest='natom', default=2,
                      help='Number of atoms')
    parser.add_option('-r', '--rs', type='float', dest='rs', default=1.23,
                      help='Wigner-Seitz radius (rs) used to determine cell volume')
    parser.add_option('-z', '--zval', type='int', dest='zval', default=5,
                      help='Number of valence electrons per atom')

    parser.add_option('--amin', type='float', dest='a_min', default=0.5,
                      help='Minimum magnitude of lattice vector "a" (before scaling)')
    parser.add_option('--amax', type='float', dest='a_max', default=1.5,
                      help='Maximum magnitude of lattice vector "a" (before scaling)')
    parser.add_option('--bmin', type='float', dest='b_min', default=0.5,
                      help='Minimum magnitude of lattice vector "b" (before scaling)')
    parser.add_option('--bmax', type='float', dest='b_max', default=1.5,
                      help='Maximum magnitude of lattice vector "b" (before scaling)')
    parser.add_option('--cmin', type='float', dest='c_min', default=0.5,
                      help='Minimum magnitude of lattice vector "c" (before scaling)')
    parser.add_option('--cmax', type='float', dest='c_max', default=1.5,
                      help='Maximum magnitude of lattice vector "c" (before scaling)')

    parser.add_option('--abmin', type='float', dest='min_ab', default=40.0,
                      help='Minimum angle between lattice vectors "a" and "b"')
    parser.add_option('--abmax', type='float', dest='max_ab', default=140.0,
                      help='Maximum angle between lattice vectors "a" and "b"')
    parser.add_option('--acmin', type='float', dest='min_ac', default=40.0,
                      help='Minimum angle between lattice vectors "a" and "c"')
    parser.add_option('--acmax', type='float', dest='max_ac', default=140.0,
                      help='Maximum angle between lattice vectors "a" and "c"')
    parser.add_option('--bcmin', type='float', dest='min_bc', default=40.0,
                      help='Minimum angle between lattice vectors "b" and "c"')
    parser.add_option('--bcmax', type='float', dest='max_bc', default=140.0,
                      help='Maximum angle between lattice vectors "b" and "c"')

    return parser


def main():

    # Initiate option parser
    parser = init_parser()
    options, args = parser.parse_args()

    # Lame way of avoiding calling a dictionary the whole time
    natom  = options.natom
    zval   = options.zval
    rs     = options.rs
    a_min  = options.a_min
    a_max  = options.a_max
    b_min  = options.b_min
    b_max  = options.b_max
    c_min  = options.c_min
    c_max  = options.c_max
    min_ab = options.min_ab
    max_ab = options.max_ab
    min_ac = options.min_ac
    max_ac = options.max_ac
    min_bc = options.min_bc
    max_bc = options.max_bc
    
    # Generate random angles between lattice vectors in radians
    d_ab = max_ab - min_ab
    d_ac = max_ac - min_ac
    d_bc = max_bc - min_bc
    theta_ab = ( min_ab + d_ab*random() ) * pi/180.0
    theta_ac = ( min_ac + d_ac*random() ) * pi/180.0
    theta_bc = ( min_bc + d_bc*random() ) * pi/180.0

    # Define vector 'a' in x-direction, based off of V^(1/3)
    nel = natom*zval
    vol = volume_from_rs(rs,nel)

    # Generate lattice vectors magnitudes unscaled by density
    mag_a = a_min + (a_max-a_min)*random()
    mag_b = b_min + (b_max-b_min)*random()
    mag_c = c_min + (c_max-c_min)*random()

    ax, ay, az = mag_a, 0.0, 0.0

    bx = mag_b*cos(theta_ab)
    by = mag_b*sin(theta_ab)
    bz = 0.0

    cx = mag_c*cos(theta_ac)
    cy = mag_c*( cos(theta_bc) - cos(theta_ac)*cos(theta_ab) ) / sin(theta_ab)
    cz = ( mag_c**2 - cx**2 - cy**2 )**0.5

    # Rescale lattice vectors for given volume
    # Make sure vectors give positive volume
    a = [ax, ay, az]
    b = [bx, by, bz]
    c = [cx, cy, cz]
    vol_0 = numpy.dot(a,numpy.cross(b,c))
    scale = (vol / abs(vol_0))**(1.0/3.0)
    if vol_0 > 0.0:
        a = [ax*scale, ay*scale, az*scale]
        b = [bx*scale, by*scale, bz*scale]
        c = [cx*scale, cy*scale, cz*scale]
    else:
        a = [ax*scale, ay*scale, az*scale]
        b = [cx*scale, cy*scale, cz*scale]
        c = [bx*scale, by*scale, bz*scale]
    ax, ay, az = a
    bx, by, bz = b
    cx, cy, cz = c
    mag_a = ( ax**2 + ay**2 + az**2 )**0.5
    mag_b = ( bx**2 + by**2 + bz**2 )**0.5
    mag_c = ( cx**2 + cy**2 + cz**2 )**0.5

    atoms = []
    for i in range(natom):
        atoms.append( [random(), random(), random()] )

    # Write coordinates to xyz file and POSCAR file
    pos = open('random.POSCAR','w')
    pos.write('Random cell, rs='+str(rs)+'\n')
    pos.write('    1.0000000000000000\n')
    pos.write('     '+str(round(ax,16))+'    '+str(round(ay,16))+'    '+str(round(az,16))+'\n')
    pos.write('     '+str(round(bx,16))+'    '+str(round(by,16))+'    '+str(round(bz,16))+'\n')
    pos.write('     '+str(round(cx,16))+'    '+str(round(cy,16))+'    '+str(round(cz,16))+'\n')
    pos.write(' '+str(natom)+'\n')
    pos.write('Direct\n')
    for atom in atoms:
        pos.write('  '+str(atom[0])+'  '+str(atom[1])+'  '+str(atom[2])+'\n')
    pos.close()

    # Write relevant information to log file
    out = open('random.log','w')
    out.write(' # data generated based on following input:\n')
    out.write(' natom:                               '+str(natom)+'\n')
    out.write(' zval:                                '+str(zval)+'\n')
    out.write(' rs:                                  '+str(rs)+'\n')
    out.write(' volume:                              '+str(vol)+'\n')
    out.write(' a_min, a_max (starting components)   '+str(a_min)+', '+str(a_max)+'\n')
    out.write(' b_min, b_max (starting components)   '+str(a_min)+', '+str(a_max)+'\n')
    out.write(' c_min, c_max (starting components)   '+str(a_min)+', '+str(a_max)+'\n')
    out.write(' min_ab, max_ab                       '+str(min_ab)+', '+str(max_ab)+'\n')
    out.write(' min_ac, max_ac                       '+str(min_ac)+', '+str(max_ac)+'\n')
    out.write(' min_bc, max_bc                       '+str(min_bc)+', '+str(max_bc)+'\n\n')
    out.write(' # results from calculations:\n')
    out.write(' Nel =                '+str(nel)+'\n')
    out.write(' vol_0 =              '+str(vol_0)+'\n')
    out.write(' lattice vector a =   ('+str(ax)+', '+str(ay)+', '+str(az)+')\n')
    out.write(' lattice vector b =   ('+str(bx)+', '+str(by)+', '+str(bz)+')\n')
    out.write(' lattice vector c =   ('+str(cx)+', '+str(cy)+', '+str(cz)+')\n')
    out.write(' |a|, |b|, |c| =      '+str(mag_a)+', '+str(mag_b)+', '+str(mag_c)+'\n')
    out.write(' theta_ab =           '+str(theta_ab*180.0/pi)+'\n')
    out.write(' theta_ac =           '+str(theta_ac*180.0/pi)+'\n')
    out.write(' theta_bc =           '+str(theta_bc*180.0/pi)+'\n')
    out.close()


if __name__ == '__main__':
    main()
