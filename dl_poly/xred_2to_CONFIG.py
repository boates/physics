#!/usr/bin/env python

# Make a CONFIG file from an xred.dat file

import os, sys, commands

global atom_type
atom_type = 'NN'

def main():

    """Perform multiple shell commands to end up with a CONFIG file"""

    try:
        natom = sys.argv[1]
        a = sys.argv[2]
        xdim = int(sys.argv[3])
        ydim = int(sys.argv[4])
        zdim = int(sys.argv[5])
    except:
        print '\nPlease provide natom and the lattice constant (angstroms) and xdim, ydim, zdim\n'
        sys.exit(0)

    try:
        xred = open('xred.dat','r')
        xred.close()
    except:
        print '\nPlease make sure the xred.dat file is present\n'
        sys.exit(0)

    ax = str(float(a)*xdim)
    ay = str(float(a)*ydim)
    az = str(float(a)*zdim)

    os.system("awk '{print $1-0.5,$2-0.5,$3-0.5}' xred.dat > xred.half")
    os.system("awk '{print $1*"+ax+",$2*"+ay+",$3*"+az+"}' xred.half > xred.half2")

    # Turn it into an xyz file
    half2 = open('xred.half2','r')
    lines = half2.readlines()
    half2.close()
    xyz = open('xred.xyz','w')
    xyz.write(natom+'\n'+'0'+'\n')

    for line in lines:
        xyz.write(atom_type+' '+line)
    xyz.close()

    # Make an improper CONFIG file
    os.system('xyz_2to_CONFIG.py xred.xyz '+str(ax)+' '+str(ay)+' '+str(az))

    # Fix the formatting of the CONFIG file and replace old
    os.system('config_spacer.py')
    os.system('mv CONFIG_fixed CONFIG')

    # Remove unnecessary tmp files
    os.system('rm xred.half xred.half2')# xred.xyz')

if __name__ == '__main__':
    main()
