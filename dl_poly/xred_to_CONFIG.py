#!/usr/bin/env python

# Make a CONFIG file from an xred.dat file

import os, sys, commands

def main():

    """Perform multiple shell commands to end up with a CONFIG file"""

    try:
        natom = sys.argv[1]
        a = sys.argv[2]
    except:
        print '\nPlease provide natom and the lattice constant (angstroms)\n'
        sys.exit(0)

    try:
        xred = open('xred.dat','r')
        xred.close()
    except:
        print '\nPlease make sure the xred.dat file is present\n'
        sys.exit(0)

    os.system("awk '{print $1-0.5,$2-0.5,$3-0.5}' xred.dat > xred.half")
    os.system("awk '{print $1*"+a+",$2*"+a+",$3*"+a+"}' xred.half > xred.half2")

    # Turn it into an xyz file
    half2 = open('xred.half2','r')
    lines = half2.readlines()
    half2.close()
    xyz = open('xred.xyz','w')
    xyz.write(natom+'\n'+'0'+'\n')

    atom_type = 'LI'
    for line in lines:
        xyz.write(atom_type+' '+line)
    xyz.close()

    # Make an improper CONFIG file
    os.system('xyz_to_CONFIG.py xred.xyz '+a)

    # Fix the formatting of the CONFIG file and replace old
    os.system('config_spacer.py')
    os.system('mv CONFIG_fixed CONFIG')

    # Remove unnecessary tmp files
    os.system('rm xred.half xred.half2')# xred.xyz')

if __name__ == '__main__':
    main()
