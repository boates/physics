#!/usr/bin/env python
"""
grow_random_chains.py
Author: Brian Boates

Generate a random cell of random shape for given rs and natom
and grow polymeric chains inside
"""
import os, sys, commands, glob
import optparse, numpy
from math import pi, sin, cos
from random import random

def get_distance(d,dtol):
    """
    Return a distance randomly generated based on the
    given average and tolerance.
    d: average bond distance (angstroms)
    dtol: tolerance about average (angstroms)
    """
    return d + dtol* (2*random()-1.0)


def get_angle(a,atol):
    """
    Return angle randomly generated based on the
    given average and tolerance.
    a: average bonding angle (degrees)
    atol: tolerance about average (degrees)
    """
    return a + atol* (2*random()-1.0)


def get_dihedral(dihed):
    """
    Return a dihedral angle based on a
    given maximum value.
    dihed: maximum dihedral angle (degrees)
    """
    return dihed * random()


def make_3_coordinated(ratio):
    """
    Return bool on whether an atom will be
    3-coordinated.
    ratio: 
    """
    return random() < ratio


def random_sign():
    """
    Return either +/- 1 randomly
    """
    return int(2*random())*2 - 1

def sign(val):
    """
    Return either +/- 1 for pos/neg values of val
    """
    return val / abs(val)


def pbc_round(input_value):
    i = int(input_value)
    if (abs(input_value-i) >= 0.5):
        if (input_value > 0): i+= 1
        if (input_value < 0): i-= 1
    return i


def too_close(atom,atoms,exclude,a,b,c):
    """
    Check that the new 'atom' is not within
    a distance of 'exclude' of any other 'atoms'
    atom: [x,y,z] for new atom
    atoms: [[x1,y1,z1],[x2,y2,z2],...[xN,yN,zN]] for current atoms
    exclude: minimum distance allowed (angstroms)
    """
    # Create supercell images of atom
    a, b, c = numpy.array(a), numpy.array(b), numpy.array(c)
 
    supercell = []
    supercell.append(atom)
    supercell.append(atom+a)
    supercell.append(atom-a)
    supercell.append(atom+b)
    supercell.append(atom-b)
    supercell.append(atom+c)
    supercell.append(atom-c)

    supercell.append(atom+a+b)
    supercell.append(atom+a-b)
    supercell.append(atom-a+b)
    supercell.append(atom-a-b)

    supercell.append(atom+a+c)
    supercell.append(atom+a-c)
    supercell.append(atom-a+c)
    supercell.append(atom-a-c)

    supercell.append(atom+b+c)
    supercell.append(atom+b-c)
    supercell.append(atom-b+c)
    supercell.append(atom-b-c)
    
    supercell.append(atom+a+b+c)
    supercell.append(atom+a+b-c)
    supercell.append(atom+a-b+c)
    supercell.append(atom+a-b-c)
    supercell.append(atom-a+b+c)
    supercell.append(atom-a+b-c)
    supercell.append(atom-a-b+c)
    supercell.append(atom-a-b-c)

    supercell = numpy.array(supercell)

    tooclose = False
    for a in atoms:
        for s in supercell:
            dr = numpy.array(a) - numpy.array(s)
            r  = (dr[0]**2 + dr[1]**2 + dr[2]**2)**0.5

            if r < exclude: tooclose = True

    return tooclose


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

    # System parameters
    parser.add_option('-n', '--natom', type='int', dest='natom', default=2,
                      help='Number of atoms')
    parser.add_option('-z', '--zval', type='int', dest='zval', default=5,
                      help='Number of valence electrons per atom')
    parser.add_option('-r', '--rs', type='float', dest='rs', default=1.25,
                      help='Wigner-Seitz radius (rs) used to determine cell volume')

    # Distance parameters
    parser.add_option('-d', '--dist', type='float', dest='dist', default=1.25,
                      help='Average distance between neighboring atoms (angstroms)')
    parser.add_option('--dtol', type='float', dest='dtol', default=0.10,
                      help='Tolerance on the distance between neighboring atoms (angstroms)')

    parser.add_option('-x', '--exclude', type='float', dest='exclude', default=1.0,
                      help='Do not allow atoms within this distance of each other (angstroms)')

    # Angle Parameters
    parser.add_option('-a','--angle', type='float', dest='angle', default=111.0,
                      help='Average angle between neighboring atoms (degrees)')
    parser.add_option('--atol', type='float', dest='atol', default=10.0,
                      help='Tolerance on the angle between neighboring atoms (degrees)')
    parser.add_option('--dihed', type='float', dest='dihed', default=10.0,
                      help='Maximum allowed dihedral angle (degrees)')

    # Empirical parameters
    parser.add_option('--ratio', type='float', dest='ratio', default=0.5,
                      help='Average fraction of 3-coordinated atoms (2-coordinated fraction = 1 - ratio)')

    # Cell parameters
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
    natom   = options.natom
    zval    = options.zval
    rs      = options.rs
    dist    = options.dist
    dtol    = options.dtol
    exclude = options.exclude
    angle   = options.angle
    atol    = options.atol
    dihed   = options.dihed
    ratio   = options.ratio
    a_min   = options.a_min
    a_max   = options.a_max
    b_min   = options.b_min
    b_max   = options.b_max
    c_min   = options.c_min
    c_max   = options.c_max
    min_ab  = options.min_ab
    max_ab  = options.max_ab
    min_ac  = options.min_ac
    max_ac  = options.max_ac
    min_bc  = options.min_bc
    max_bc  = options.max_bc
    
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
    space = numpy.matrix( [a, b, c] )
    inverse = space.I

    alat = mag_a
    blat = mag_b
    clat = mag_c

    # Grow the chains and fill the box with natom atoms
    atoms = [ numpy.array([0.0, 0.0, 0.0]) ]
    while len(atoms) != natom:

        # Final routine to add third or more atoms
#        if len(atoms) >= 3:

            # Get the bond distance and angles
#            dist = get_distance(d,dtol)
#            ang  = get_angle(angle,atol)
#            dihed = get_dihedral(dihed)

        # Different routine to add the third atom
        if len(atoms) >= 2:

            # Get the bond distance and angle
            d   = get_distance(dist,dtol)
            ang = get_angle(angle,atol)

            # Determine the new bond vector
            oldx = previous_vec[0]
            oldy = previous_vec[1]
            oldz = previous_vec[2]
            newx = (random()+0.5)* random_sign()
            newy = (random()+0.5)* random_sign()
            newz = (1./oldz) * (previous_mag* d* cos((180.-ang)*pi/180.) - oldx*newx - oldy*newy)

            new  = numpy.array([newx,newy,newz])
            mag = (new[0]**2 + new[1]**2 + new[2]**2)**0.5

            new = atoms[-1] + new * (d/mag)

            # Wrap the new atom's coordinate
            reduced = new * inverse
            new = numpy.array( [reduced.item(0), reduced.item(1), reduced.item(2)] )
            for i in range(len(new)):
                while abs(new[i]) >= 1: new[i] -= sign(new[i])

            # Convert back to angstroms
            cart = new * space
            new = numpy.array( [cart.item(0), cart.item(1), cart.item(2)] )

            if not too_close(new,atoms,exclude,a,b,c):
                atoms.append(new)

                # Track the previous vector for next atom
                previous_vec = new
                previous_mag = d

                print len(atoms), d, ang, new


        # Different routine to add the second atom
        if len(atoms) == 1:
            
            # Get the bond distance
            d = get_distance(dist,dtol)

            # Determine the 'bond vector' (for orientation)
            new = numpy.array([random(), random(), random()])
            mag = (new[0]**2 + new[1]**2 + new[2]**2)**0.5

            # Rescale based on given bond distance
            new = atoms[-1] + ( new * (d/mag) )

            # Wrap the new atom's coordinate
            reduced = new * inverse
            new = numpy.array( [reduced.item(0), reduced.item(1), reduced.item(2)] )
            for i in range(len(new)):
                while abs(new[i]) >= 1: new[i] -= sign(new[i])

            # Convert back to angstroms
            cart = new * space
            new = numpy.array( [cart.item(0), cart.item(1), cart.item(2)] )

            # If atom isn't too close to close to another
            if not too_close(new,atoms,exclude,a,b,c):
                atoms.append(new)

                # Track the previous vector for next atom
                previous_vec = new
                previous_mag = d


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
    out = open('grow_random.log','w')
    out.write(' # data generated based on following input:\n')
    out.write(' natom:                               '+str(natom)+'\n')
    out.write(' zval:                                '+str(zval)+'\n')
    out.write(' rs:                                  '+str(rs)+'\n')
    out.write(' volume:                              '+str(vol)+'\n')
    out.write(' dist, dtol                           '+str(dist)+', '+str(dtol)+'\n')
    out.write(' exclude:                             '+str(exclude)+'\n')
    out.write(' angle, atol                          '+str(angle)+', '+str(atol)+'\n')
    out.write(' dihed:                               '+str(dihed)+'\n')
    out.write(' ratio:                               '+str(ratio)+'\n')
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
