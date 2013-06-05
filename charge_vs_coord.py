#!/usr/bin/env python
"""
Read in from a bader ACF.dat file, the respective charges on each atom
and compare with the coordination of each atom.
"""
import os, sys, commands, glob

def main():

    try:
        f = open('ACF.dat','r')
        g = open(glob.glob('coordination_rc_*.dat')[0],'r')
    except:
        print '\nMake sure a coordination file and an ACF.dat file are both present.\n'
        sys.exit(0)

    # Remove header info
    f.readline()
    f.readline()
    for i in range(10):
        g.readline()
    natom = int(commands.getoutput('wc -l coordination_rc_*.dat').split()[0]) - 10

    charge, coord = [], []
    for i in range(natom):
        charge.append(f.readline().split()[4])
        coord.append(g.readline().split()[0])

    out = open('cc.dat','w')
#    out.write('# coord, charge\n')
    for i in range(natom):
        out.write(coord[i]+' '+charge[i]+'\n')

    f.close()
    g.close()
    out.close()

if __name__ == '__main__':
    main()
