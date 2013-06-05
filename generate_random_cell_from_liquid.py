#!/usr/bin/env python
"""
generate_random_cell_from_liquid.py

Sift through snapshots of liquid data from a wrapped xyz file
and create randomly shaped cells about atoms withing a snapshot
to generate a unit cell of desired number of atoms and density.
"""
import os, sys, commands, glob
import optparse, numpy
from math import pi, sin, cos
from random import random

def inside_polygon(x,y,poly):
    """
    Determine if a point is inside a given polygon or not
    Polygon is a list of (x,y) pairs. This fuction
    returns True or False.  The algorithm is called
    Ray Casting Method

    http://pseentertainmentcorp.com/smf/index.php?topic=545.0

    return True / False
    """
    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside


def volume_from_rs(rs,Nel):
    """
    Return cell volume in angstroms^3.
    rs: density parameter for system
    Nel: number of valence electrons in system
    """
    a0 = 0.5291772   # Bohr radius (angstroms/bohr)
    volume = (4.0*pi/3.0)*Nel * (rs*a0)**3

    return volume


def get_supercell(x,y,z,a):

    # Create the supercell array
    X, Y, Z = [], [], []
    for i in [-1,0,1]:
        for j in [-1,0,1]:
            for k in [-1,0,1]:
                X.append(x+i*a)
                Y.append(y+j*a)
                Z.append(z+k*a)

    return X, Y, Z

def init_parser():

    usage = 'usage: %prog [options]'
    version="%prog 1.0"
    parser = optparse.OptionParser(usage=usage,version=version)

    parser.add_option('-f', '--file', type='str', dest='fname', default='wrapped.xyz',
                      help='wrapped_liquid.xyz, default=wrapped.xyz')
    parser.add_option('-a', '--alat', type='float', dest='alat', default=None,
                      help='Liquid lattice parameter (Angtroms), default=None')
    parser.add_option('-c', '--nCells', type='int', dest='nCells', default=1,
                      help='Number of cells to generate, default=1')
    parser.add_option('-s', '--nSkip', type='int', dest='nSkip', default=100,
                      help='Number of configurations to skip in between generations, default=100')

    parser.add_option('-n', '--natom', type='int', dest='natom', default=2,
                      help='Number of atoms, default=2')
    parser.add_option('-r', '--rs', type='float', dest='rs', default=1.23,
                      help='Wigner-Seitz radius (rs) used to determine cell volume, default=1.23')
    parser.add_option('-z', '--zval', type='int', dest='zval', default=5,
                      help='Number of valence electrons per atom, default=5')

    parser.add_option('--amin', type='float', dest='a_min', default=0.5,
                      help='Minimum magnitude of lattice vector "a" (before scaling), default=0.5')
    parser.add_option('--amax', type='float', dest='a_max', default=1.5,
                      help='Maximum magnitude of lattice vector "a" (before scaling), default=1.5')
    parser.add_option('--bmin', type='float', dest='b_min', default=0.5,
                      help='Minimum magnitude of lattice vector "b" (before scaling), default=0.5')
    parser.add_option('--bmax', type='float', dest='b_max', default=1.5,
                      help='Maximum magnitude of lattice vector "b" (before scaling), default=1.5')
    parser.add_option('--cmin', type='float', dest='c_min', default=0.5,
                      help='Minimum magnitude of lattice vector "c" (before scaling), default=0.5')
    parser.add_option('--cmax', type='float', dest='c_max', default=1.5,
                      help='Maximum magnitude of lattice vector "c" (before scaling), default=1.5')

    parser.add_option('--abmin', type='float', dest='min_ab', default=40.0,
                      help='Minimum angle between lattice vectors "a" and "b", default=40.0')
    parser.add_option('--abmax', type='float', dest='max_ab', default=140.0,
                      help='Maximum angle between lattice vectors "a" and "b", default=140.0')
    parser.add_option('--acmin', type='float', dest='min_ac', default=40.0,
                      help='Minimum angle between lattice vectors "a" and "c", default=40.0')
    parser.add_option('--acmax', type='float', dest='max_ac', default=140.0,
                      help='Maximum angle between lattice vectors "a" and "c", default=140.0')
    parser.add_option('--bcmin', type='float', dest='min_bc', default=40.0,
                      help='Minimum angle between lattice vectors "b" and "c", default=40.0')
    parser.add_option('--bcmax', type='float', dest='max_bc', default=140.0,
                      help='Maximum angle between lattice vectors "b" and "c", default=140.0')

    return parser


def main():

    # Initiate option parser
    parser = init_parser()
    options, args = parser.parse_args()

    # Lame way of avoiding calling a dictionary the whole time
    fname  = options.fname
    alat   = options.alat
    nCells = options.nCells
    nSkip  = options.nSkip
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

    # Print useful information for user
    if alat == None:
        print '\n Parameter alat (the lattice parameter of the cubic liquid system) '+\
              'must be given using the -a option - exiting\n'
        sys.exit(0)
    if '-f' not in sys.argv:
        try:
            xyz = open(fname,'r')
            xyz.close()
        except IOError:
            print '\n Code assumes liquid xyz file is named wrapped.xyz - Could not detect file. '+\
                  'To change this use the -f option\n'
            sys.exit(0)
    print '\n To see all default values for tunable parameters use the -h option\n'

    # Open xyz file and determine nConfig and nLiq
    xyz = open(fname,'r')
    lines = xyz.readlines()
    xyz.close()
    nLiq = int(lines[0])
    nConfig = int( len(lines)/(nLiq+2.0) )
    xyz = open(fname,'r')

    # Check to make sure there are enough configurations
    if nCells > nConfig / nSkip * nLiq:
        print '\n You have asked for too many cells to be generated given'
        print ' the number of configs and atoms in this xyz file - exiting\n'
        sys.exit(0)

    # Create and change into a new directory where
    # the POSCAR and log files will be written
    dirname = 'RANDOM'
    os.system('rm -rf '+dirname)
    os.system('mkdir '+dirname)
    os.chdir(dirname)

    # Generate the cells about liquid atoms
    breaker = False
    cells = 0
    for i in range( int(nConfig/nSkip) ):

        if breaker == True: break

        # Read off natom and timestep
        xyz.readline()
        xyz.readline()
            
        x, y, z = [], [], []
        X, Y, Z = [], [], []
        for j in range(nLiq):
            row = xyz.readline().split()
            x.append(float(row[1]))
            y.append(float(row[2]))
            z.append(float(row[3]))

            f, g, h = get_supercell(float(row[1]),float(row[2]),float(row[3]),alat)
            X += f
            Y += g
            Z += h
        
        for j in range(nLiq):

            if cells >= nCells:
                print ' Current timestep: '+str(i)+' --- found '+str(cells)+' cells out of '+str(nCells)
                print '\n Found '+str(cells)+' cells as requested - exiting\n'
                breaker = True
                break

            # Rename current atom's coordinates
            x0 = x[j]
            y0 = y[j]
            z0 = z[j]
            atom = numpy.array([x0, y0, z0])

            # Generate random angles between lattice vectors in radians
            d_ab = max_ab - min_ab
            d_ac = max_ac - min_ac
            d_bc = max_bc - min_bc
            theta_ab = ( min_ab + d_ab*random() ) * pi/180.0
            theta_ac = ( min_ac + d_ac*random() ) * pi/180.0
            theta_bc = ( min_bc + d_bc*random() ) * pi/180.0

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
            try:
                proceed = True
                cz = ( mag_c**2 - cx**2 - cy**2 )**0.5
            except ValueError:
                proceed = False

            if proceed:
            
                # Calculate volume for desired density
                nel = natom*zval
                vol = volume_from_rs(rs,nel)
        
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
    
                # Now do the point in polygon tests:
                # The cell has 6 sides & 8 corners
    
                # Place the current atom at (0, 0, 0)
                c1 = atom
                c2 = atom + numpy.array(a)
                c3 = atom + numpy.array(b) 
                c4 = atom + numpy.array(a) + numpy.array(b)
                c5 = atom + numpy.array(c)
                c6 = atom + numpy.array(a) + numpy.array(c)
                c7 = atom + numpy.array(b) + numpy.array(c)
                c8 = atom + numpy.array(a) + numpy.array(b) + numpy.array(c)
    
                # Define the sides:
                left   = [ (c1[1],c1[2]), (c3[1],c3[2]), (c7[1],c7[2]), (c5[1],c5[2]) ]
                right  = [ (c2[1],c2[2]), (c4[1],c4[2]), (c8[1],c8[2]), (c6[1],c6[2]) ]
                front  = [ (c1[0],c1[2]), (c2[0],c2[2]), (c6[0],c6[2]), (c5[0],c5[2]) ]
                back   = [ (c3[0],c3[2]), (c4[0],c4[2]), (c8[0],c8[2]), (c7[0],c7[2]) ]
                bottom = [ (c1[0],c1[1]), (c2[0],c2[1]), (c4[0],c4[1]), (c3[0],c3[1]) ]
                top    = [ (c5[0],c5[1]), (c6[0],c6[1]), (c8[0],c8[1]), (c7[0],c7[1]) ]
    
                atoms = []
                for k in range(nLiq*3*3*3):
    
                    # Check if current atom (or its supercell replicated self) is inside the cell
                    in_left   = inside_polygon(Y[k],Z[k],left)
                    in_right  = inside_polygon(Y[k],Z[k],right)
                    in_front  = inside_polygon(X[k],Z[k],front)
                    in_back   = inside_polygon(X[k],Z[k],back)
                    in_bottom = inside_polygon(X[k],Y[k],bottom)
                    in_top    = inside_polygon(X[k],Y[k],top)
    
                    if in_left and in_right and in_front and in_back and in_bottom and in_top:
                        atoms.append( numpy.array([ X[k], Y[k], Z[k] ]) )
    
                # If the number of atoms found inside the cell is less than the
                # number of atoms wanted then skip this round
                still_proceed = True
                if len(atoms) != natom:
                    if len(atoms) < natom:
                        still_proceed = False
                    while len(atoms) > natom:
                        atoms.pop(-1)
    
                if still_proceed:
    
                    # Shift the atoms back to the origin and put in reduced coordinates
                    shifted_and_reduced = []
                    for ion in atoms:
                        shift = ion - atom
                        ion = shift * inverse
                        shifted_and_reduced.append( [ion.item(0), ion.item(1), ion.item(2)] )
                    atoms = shifted_and_reduced

#                    print atoms
            
                    # Write coordinates to xyz file and POSCAR file
                    pos = open('POSCAR.liq-'+str(natom)+'atoms-'+str(i)+'-'+str(j),'w')
                    pos.write('Random cell, rs='+str(rs)+'\n')
                    pos.write('    1.0000000000000000\n')
                    pos.write('     '+str(round(ax,16))+'    '+str(round(ay,16))+'    '+str(round(az,16))+'\n')
                    pos.write('     '+str(round(bx,16))+'    '+str(round(by,16))+'    '+str(round(bz,16))+'\n')
                    pos.write('     '+str(round(cx,16))+'    '+str(round(cy,16))+'    '+str(round(cz,16))+'\n')
                    pos.write(' '+str(natom)+'\n')
                    pos.write('Direct\n')
                    for ion in atoms:
                        pos.write('  '+str(ion[0])+'  '+str(ion[1])+'  '+str(ion[2])+'\n')
                    pos.close()

                    # Write relevant information to log file
                    out = open('random.log-'+str(natom)+'atoms-'+str(i)+'-'+str(j),'w')
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

                    # Count how many cells have been generated
                    cells += 1

        # Skip a given number of configurations to avoid reproducing
        # very similar configurations
        for j in range( (nLiq+2) *nSkip):
            xyz.readline()


if __name__ == '__main__':
    main()
