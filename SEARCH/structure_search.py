#!/usr/bin/env python
"""
structure_search.py
Version 1.0

Author: Brian Boates

Perform random ab initio structural searches
"""
import os, sys, commands, glob, copy
import numpy
from math import pi, sin, cos, asin, acos
from random import random

def masses():
    """
    Creates an atomic mass dictionary
    """
    mass = {'H' : 1.00794,'He': 4.002602,'Li': 6.941,'Be': 9.012182,'B' : 10.811,'C' : 12.0107,
            'N' : 14.0067,'O' : 15.9994,'F' : 18.9984032,'Ne': 20.1797,'Na': 22.98976928,
            'Mg': 24.3050,'Al': 26.9815386,'Si': 28.0855,'P' : 30.973762,'S' : 32.065,'Cl': 35.453,
            'Ar': 39.948,'K' : 39.0983,'Ca': 40.078,'Sc': 44.955912,'Ti': 47.867,'V' : 50.9415,
            'Cr': 51.9961,'Mn': 54.938045,'Fe': 55.845,'Co': 58.933195,'Ni': 58.6934,'Cu': 63.546,
            'Zn': 65.38,'Ga': 69.723,'Ge': 72.64,'As': 74.92160,'Se': 78.96,'Br': 79.904,'Kr': 83.798,
            'Rb': 85.4678,'Sr': 87.62,'Y' : 88.90585,'Zr': 91.224,'Nb': 92.90638,'Mo': 95.96,'Tc': 98.0,
            'Ru': 101.07,'Rh': 102.90550,'Pd': 106.42,'Ag': 107.8682,'Cd': 112.411,'In': 114.818,
            'Sn': 118.710,'Sb': 121.760,'Te': 127.60,'I' : 126.90447,'Xe': 131.293,'Cs': 132.9054519,
            'Ba': 137.327,'La': 6.145,'Ce': 6.77,'Pr': 6.773,'Nd': 7.007,'Pm': 7.26,'Sm': 7.52,
            'Eu': 5.243,'Gd': 7.895,'Tb': 8.229,'Dy': 8.55,'Ho': 8.795,'Er': 9.066,'Tm': 9.321,
            'Yb': 6.965,'Lu': 174.9668,'Hf': 178.49,'Ta': 180.94788,'W' : 183.84,'Re': 186.207,
            'Os': 190.23,'Ir': 192.217,'Pt': 195.084,'Au': 196.966569,'Hg': 200.59,'Tl': 204.3833,
            'Pb': 207.2,'Bi': 208.98040,'Po': 210.0,'At': 210.0,'Rn': 222.0,'Fr': 223.0,'Ra': 226.0}
    return mass


class lattice:
    """
    A 'lattice' object contains information about the lattice
    vectors, components, volume, etc.
    """
    def __init__(self,a,b,c):
        """
        a: lattice vector a (type=numpy.array (of 3 floats))
        b: lattice vector b (type=numpy.array (of 3 floats))
        c: lattice vector c (type=numpy.array (of 3 floats))

        Initialize lattice
        """
        # Define self lattice vectors
        self.a = a
        self.b = b
        self.c = c

        # Define self lattice vector components
        self.ax = a[0]
        self.ay = a[1]
        self.az = a[2]
        self.bx = b[0]
        self.by = b[1]
        self.bz = b[2]
        self.cx = c[0]
        self.cy = c[1]
        self.cz = c[2]

        # Define self lattice vector magnitudes
        self.amag = (a[0]**2 + a[1]**2 + a[2]**2)**0.5
        self.bmag = (b[0]**2 + b[1]**2 + b[2]**2)**0.5
        self.cmag = (c[0]**2 + c[1]**2 + c[2]**2)**0.5

        # Define volume
        self.vol = numpy.dot(a,numpy.cross(b,c))


    def rescale_vol(self,scale):
        """
        scale: factor to rescale volume by

        Randomly rescale the volume of lattice up to maxvol to
        attempt to accomodate space for new atoms/units

        return rescaled_lattice (type=lattice)
        """
        a_scale = self.a * scale**(1.0/3.0)
        b_scale = self.b * scale**(1.0/3.0)
        c_scale = self.c * scale**(1.0/3.0)

        scaled_lattice = lattice(a_scale, b_scale, c_scale)

        return scaled_lattice


class unit:
    """
    A unit is something that represents a molecule and 
    contains relevant info about itself such as number 
    of atoms, atomic positions, atom types, etc.
    """
    def __init__(self,input,begin,end):
        """
        input: String of input from SEARCH file
        begin: index to begin reading in unit data
        end: index to end reading in unit data

        Initialize 'unit'
        """
        # Read in the raw text as a string
        raw = input[begin:end].split('\n')

        # Unit name is the first line in unit definition
        try:
            name = raw[0].split()[2]
        # If no name is given, assign None
        except IndexError:
            name = None

        # Number of atoms in unit is second line in unit definition
        try:
            natom = int(raw[1].split()[2])
        # If issue reading natom
        except:
            print '\n Problem reading ATOMS from NEWUNIT '+name+'\n'
            print ' Check input file --- exiting...\n'
            sys.exit(1)

        # Loop over number of atoms in unit
        types, x, y, z = [], [], [], []
        for i in range(natom):

            # types gives each atom type in the unit
            types.append(raw[2+i].split()[0])

            # x, y, & z give the x, y, and z coordinates
            # for each atom in the unit (in Ang)
            try:
                x.append(float(raw[2+i].split()[1]))
                y.append(float(raw[2+i].split()[2]))
                z.append(float(raw[2+i].split()[3]))
            # If issue reading atom positions
            except:
                print '\n Problem reading coordinates for'
                print ' atom '+str(i+1)+' (type='+types[i]+') in NEWUNIT "'+name+'"\n'
                print ' Check input file --- exiting...\n'
                sys.exit(1)
            
        # Define the self variables
        self.name  = name
        self.natom = natom
        self.types = types
        self.x     = x
        self.y     = y
        self.z     = z

        # Compute center of mass
        mass = masses()
        xcom, ycom, zcom, mtot = 0.0, 0.0, 0.0, 0.0
        for i in range(natom):
            xcom += mass[types[i]]*x[i]
            ycom += mass[types[i]]*y[i]
            zcom += mass[types[i]]*z[i]
            mtot += mass[types[i]]

        # Create self variable for center of mass (type=numpy.array)
        self.com = numpy.array([xcom, ycom, zcom]) / mtot

        # Create self variable for total mass of unit (type=float)
        self.mtot = mtot
        

    def type(self,index):
        """
        index: atom index in unt (type=int)
               - begins at 1!

        return: atomic species (type=str)
        """
        i = index - 1
        return self.types[i]
        

    def atom(self,index):
        """
        index: atom index in unit (type=int)
               - begins at 1!

        return: numpy.array of [x, y, z] (in Ang)
        """
        i = index - 1
        return numpy.array([self.x[i], self.y[i], self.z[i]])


    def shift(self,x_shift,y_shift,z_shift):
        """
        Shifts unit coordinates

        x_shift: amount to shift unit coordinates by in x direction (float)
        y_shift: amount to shift unit coordinates by in x direction (float)
        z_shift: amount to shift unit coordinates by in x direction (float)
        """
        new_x = numpy.array( [ x+x_shift for x in self.x ] )
        new_y = numpy.array( [ y+y_shift for y in self.y ] )
        new_z = numpy.array( [ z+z_shift for z in self.z ] )

        self.x = new_x
        self.y = new_y
        self.z = new_z
        
        return 1


    def rotate(self,x_ang,y_ang,z_ang):
        """
        Rotates unit
        
        x_ang: angle to rotate about x-axis (degrees), (type=float)
        y_ang: angle to rotate about y-axis (degrees), (type=float)
        z_ang: angle to rotate about y-axis (degrees), (type=float)
        """
        # Convert each angle from degrees to radians
        x_ang = x_ang * pi/180.0
        y_ang = y_ang * pi/180.0
        z_ang = z_ang * pi/180.0

        # Define x rotation matrix
        Rx = numpy.matrix([ [1.0, 0.0, 0.0],
                            [0.0, cos(x_ang), -sin(x_ang)],
                            [0.0, sin(x_ang),  cos(x_ang)] ])
        
        # Define y rotation matrix
        Ry = numpy.matrix([ [cos(y_ang), 0.0, sin(y_ang)],
                            [0.0, 1.0, 0.0],
                            [-sin(y_ang), 0.0, cos(y_ang)] ])
        
        # Define z rotation matrix
        Rz = numpy.matrix([ [cos(z_ang), -sin(z_ang), 0.0],
                            [sin(z_ang),  cos(z_ang), 0.0],
                            [0.0, 0.0, 1.0] ])

        # Combine into one rotation matrix
        R = numpy.dot( Rz, numpy.dot(Ry,Rx) )

        # Rotate the coordinates of each atom in the unit
        unit_rot = []
        for i in range(self.natom):

            # Get vector of x,y,z for atom i (indexing starts at 1)
            vec = self.atom(i+1)

            # Rotate the coordinate
            vec_rot = numpy.array( numpy.dot(R, vec) )

            # Append rotated coordinate to new rotated unit array
            unit_rot.append( [v for v in vec_rot[0]] )  # round to avoid tiny/weird E-17 numbers

        # Reassign the self coordinates for the unit
        self.x = numpy.array( [c[0] for c in unit_rot] )
        self.y = numpy.array( [c[1] for c in unit_rot] )
        self.z = numpy.array( [c[2] for c in unit_rot] )
        
        return 1


    def reduce(self,a,b,c):
        """
        a: lattice vector a (type=numpy.array)
        b: lattice vector b (type=numpy.array)
        c: lattice vector c (type=numpy.array)

        return reduced_unit_coords (type=numpy.array of numpy.arrays of floats)
        """
        # Define lattice matrix in terms of cell vectors
        lat = numpy.matrix( [a, b, c] )
#        lat = numpy.matrix( [numpy.array([5.0,0.0,0.0]),
#                             numpy.array([5.0,5.0,0.0]),
#                             numpy.array([0.0,0.0,5.0])] )
        print lat
        # Take the inverse
        lat_inv = lat.I
        print lat_inv

        # Loop over each atom coord in unit
        unit_red = []
        print self.atom(1)
        for i in range(self.natom):

            # Get vector of x,y,z for atom i (indexing starts at 1)
            vec = self.atom(i+1)
            print '================================',i
            print vec

            # Convert to reduced coordinates
            vec_red = numpy.array( numpy.dot(numpy.array(lat_inv), vec) )
            print vec_red
            vec_orig = numpy.array( numpy.dot(numpy.array(lat), vec_red) )
            print vec_orig
            print [round(v,12) for v in vec_orig]

            # Append reduced coordinate to new reduced unit array
            unit_red.append( vec_red )  # round to avoid tiny/weird E-17 numbers

        # Update x, y, & z lists for unit
        self.x = numpy.array( [c[0] for c in unit_red] )
        self.y = numpy.array( [c[1] for c in unit_red] )
        self.z = numpy.array( [c[2] for c in unit_red] )
        print self.x
        print self.y
        print self.z

        # Compute reduced center of mass coordinates
        com_red = numpy.array( numpy.dot(lat_inv, self.com) )

        # Rename self.com variable to reduced value
        self.com = com_red[0]
                                                                                            
        return 1


def default(var):
    """
    var: Variable name (type=str)
    
    return: default value for var
    """
    if   var == 'VOLUME':      return None
    elif var == 'MAXVOL':      return 1.50
    elif var == 'PRESSURE':    return None
    elif var == 'ANGMIN':      return 40.0
    elif var == 'ANGMAX':      return 140.0
    elif var == 'ALATLIM':     return 4.0
    elif var == 'MINDIST':     return 0.7
    elif var == 'NCELLS':      return 100
    elif var == 'NTYPES':      return None
    elif var == 'NUNITS':      return None
    elif var == 'NCOARSE':     return 3
    elif var == 'NSTEPCOARSE': return 3
    elif var == 'NFINE':       return 4
    elif var == 'NSTEPFINE':   return 20
    elif var == 'POTCAR':      return './'
    else:                      return None
    

def get_var(input,var,unpack=1):
    """
    input: String of SEARCH input file (type=str)
    var: Variable name (type=str)
    unpack: Number of items to unpack for var (type=int)

    return: value given for var (type=str)
            values if unpack > 1 for var (type=list[str])
    """
    # Get beginning and end indices in input for var
    begin = input.find(var)
    end   = input.find('\n',begin)

    # If var is not found, retrieve default value from default() function
    if begin == -1:
        val = default(var)
        # If there is no default, alert user and prepare to exit
        if val == None:
            print '\n variable "'+var+'" has no default, please specify in input file --- exiting...\n'
        return val

    # Check for how many values need to be unpacked for var
    if unpack == 1:
        return input[begin:end].split()[2]
    else:
        return input[begin:end].split()[2:2+unpack]


def get_units(input,ntypes):
    """
    input: String of SEARCH input file (type=str)
    ntypes: Number of unit types (type=int)

    return: units list, each item in the units list
            contains the atom types and relative coordinates
            for each unit.

            i.e.: [ unit1, unit2, ..., unit_ntypes ]

            where unit1, etc. are defined by the unit class
    """
    # Initialize begin index
    b0 = 0

    # Loop over number of unit types
    units = []
    for i in range(ntypes):

        # Determine beginning and ending indices for current unit
        begin  = input.find('NEWUNIT',b0)
        end    = input.find('ENDUNIT',begin)

        # Make sure proper number of units are given
        if begin == -1:
            print '\n Problem with "NEWUNIT" and/or "ENDUNIT" directives in input file:'
            print ' Make sure there are as many NEWUNIT directives as value NTYPES \n'
            print ' exiting...\n'
            break
        
        b0 = end

        # Create and append a unit item to the units list
        units.append( unit(input,begin,end) )

    return units


def make_vectors(volume,angmin,angmax,alatlim):
    """
    volume: total cell volume in Ang^3
    angmin: min angle between lattice vectors
    angmax: max angle between lattice vectors
    alatlim: max ratio between lattice vector lengths
    
    return: lattice vectors a, b, & c (type=numpy.array)
    """
    # Generate lattice vector magnitudes before scaling by volume
    x = alatlim - 1.0
    a_mag = x*random() + 1.0
    b_mag = x*random() + 1.0
    c_mag = x*random() + 1.0

    # Set a vector in x-direction
    ax = a_mag
    ay = 0.0
    az = 0.0
    a = numpy.array( [ax, ay, az] )

    # Set a-b surface to x-y plane
    ang_ab = ( (angmax-angmin)*random() + angmin )  # degrees
    bx = b_mag*cos(ang_ab *pi/180.)
    by = b_mag*sin(ang_ab *pi/180.)
    bz = 0.0
    b = numpy.array( [bx, by, bz] )

    # Determine c vector from a & b vectors
    # Determine possible value for ang_bc within desired range
    # see http://journals.iucr.org/a/issues/2011/01/00/au5114/au5114bdy.html

    # Create fake ang_ac & ang_bc both = 0
    ang_ac, ang_bc = 0.0, 0.0

    # Loop over ang_bc values until a suitable ang_ac is also generated
    while (ang_ac < angmin or ang_ac > angmax or ang_bc < angmin or ang_bc > angmax):
        ang_bc = ( (angmax-angmin)*random() + angmin )  # degrees
#        print ang_bc
        cx = c_mag * cos(ang_bc) * cos(ang_ab)
        cy = c_mag * cos(ang_bc) * sin(ang_ab)
        cz = c_mag * sin(ang_bc)
        c = numpy.array( [cx, cy, abs(cz)] )  ## force z-component to be positive (doesn't change things)
        c_mag = ( cx**2 + cy**2 + cz**2 )**0.5
        ang_ac = acos(numpy.dot(a,c) / (a_mag*c_mag)) * (180./pi)
        ang_bc = acos(numpy.dot(b,c) / (b_mag*c_mag)) * (180./pi)

#    print 'theta_ab = ', acos(numpy.dot(a,b) / (a_mag*b_mag)) * (180./pi)
#    print 'theta_ac = ', acos(numpy.dot(a,c) / (a_mag*c_mag)) * (180./pi)
#    print 'theta_bc = ', acos(numpy.dot(b,c) / (b_mag*c_mag)) * (180./pi)

    # Compute the current volume of the cell
    v0 = numpy.dot(a,numpy.cross(b,c))

    # Compute the scaling factor for desired volume
    scale = (volume / v0)**(1.0/3.0)

    # Re-scale the lattice vectors and magnitudes to give desired volume
    a = a*scale
    b = b*scale
    c = c*scale
    a_mag = a_mag*scale
    a_mag = a_mag*scale
    a_mag = a_mag*scale

    return a, b, c
    

def main(fin='SEARCH'):
    """
    fin: Input filename for structure_search.py (default='SEARCH')

    Reads in and returns the input parameters from input file, SEARCH
    """
    # Read the input file
    f = open(fin,'r')
    input = f.read()
    f.close()

    # Read in each variable
    try:
        volume      = float(get_var(input,'VOLUME'))   ## total cell volume in Ang^3
        maxvol      = float(get_var(input,'MAXVOL'))   ## max rescaling of volume (i.e. 1.5)
        pressure    = float(get_var(input,'PRESSURE')) ## pressure for optimization in GPa
        angmin      = float(get_var(input,'ANGMIN'))   ## min angle between lattice vectors
        angmax      = float(get_var(input,'ANGMAX'))   ## max angle between lattice vectors
        alatlim     = float(get_var(input,'ALATLIM'))  ## max ratio of lattice vector lengths
        mindist     = float(get_var(input,'MINDIST'))  ## min distance between atoms in cell
        ncells      = int(get_var(input,'NCELLS'))     ## number of cells to generate
        ntypes      = int(get_var(input,'NTYPES'))     ## number of different "units"
        ncoarse     = int(get_var(input,'NCOARSE'))    ## number of coarse optimizations
        nstepcoarse = int(get_var(input,'NSTEPCOARSE')) ## number of steps for coarse optimiaztions
        nfine       = int(get_var(input,'NFINE'))       ## number of fine optimizations
        nstepfine   = int(get_var(input,'NSTEPFINE'))   ## number of steps for fine optimizations
        potcar      = get_var(input,'POTCAR')           ## filepath to POTCAR files
    except:
        sys.exit(1)
    try:
        # Length of ntypes, gives number of each unit for the cells
        nunits = [int(v) for v in get_var(input,'NUNITS',unpack=ntypes)]  ## list of how many of each unit to use
    except:
        print '\n Proble reading NUNITS - Make sure there are NTYPES specified'
        print ' exiting...\n'
        sys.exit(1)


    ######################
    ##  BEGIN ROUTINES  ##
    ######################

    # Retrieve the "units" from input (i.e. molecules, atoms, etc.)
    units = get_units(input,ntypes)

    # Make random lattice vectors and create cell lattice object
    a, b, c = make_vectors(volume, angmin, angmax, alatlim)
    print a
    print b
    print c
    cell = lattice(a,b,c)

    out = open('POSCAR','w')

    # Total number of atoms types
    N = numpy.dot( numpy.array([t.natom for t in units]) , nunits )

    # atomic species list
    raw = []
    for i in range(ntypes):
        raw += units[i].types
    species = []
    for s in raw:
        if s not in species:
            species.append(s)
#    print species

    raw = [t.types for t in units]
#    print raw
    the = []
    for i in range(ntypes):
        for j in range(nunits[i]):
            the += raw[i]
#    print the

    nspecies = []
    for s in species:
        nspecies.append(the.count(s))
#    print nspecies
        
    out.write('POSCAR\n  1.00000000\n')
    out.write('    '+str(cell.ax)+'  '+str(cell.ay)+'  '+str(cell.az)+'\n')
    out.write('    '+str(cell.bx)+'  '+str(cell.by)+'  '+str(cell.bz)+'\n')
    out.write('    '+str(cell.cx)+'  '+str(cell.cy)+'  '+str(cell.cz)+'\n')
    for s in species:
        out.write('  '+s)
    out.write('\n')
    for n in nspecies:
        out.write('  '+str(n))
    out.write('\n')
    out.write('C Direct\n')

    typat, x, y, z = [], [], [], [] 
    for i in range(ntypes):
        for j in range(nunits[i]):
            u = copy.copy(units[i])

            # Rotate the unit randomly
            u.rotate(random()*360., random()*360., random()*360.)

            # Shift the unit randomly
            ra = random()
            rb = random()
            rc = random()
            x_shift = ra*cell.ax + rb*cell.bx + rc*cell.cx
            y_shift = ra*cell.ay + rb*cell.by + rc*cell.cy
            z_shift = ra*cell.az + rb*cell.bz + rc*cell.cz
            u.shift(x_shift,y_shift,z_shift)

            # Wrap the unit coordinates
            #wrap

            # Convert to reduced coordinates
#            u.reduce(a,b,c)
            

            for k in range(u.natom):

                typat.append(u.types[k])
                x.append(u.x[k])
                y.append(u.y[k])
                z.append(u.z[k])

    typat_new, x_new, y_new, z_new = [], [], [], []
    for i in range(len(species)):
        s = species[i]
        for j in range(len(x)):
            if typat[j] == s:
                typat_new.append(typat[j])
                x_new.append(x[j])
                y_new.append(y[j])
                z_new.append(z[j])
#    print typat_new
    for i in range(N):
        out.write(' '+str(x_new[i])+'  '+str(y_new[i])+'  '+str(z_new[i])+'\n')
    out.close()






    
    # bottom











if __name__ == '__main__':
    main()
