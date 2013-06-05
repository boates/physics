#!/usr/bin/env python

# Script to pull out atomic coords from pwscf output file
# and make an xyz file from them

import os, sys, commands

def main():
    
    try:
        fname = sys.argv[1]
        natom = int(sys.argv[2])
    except IndexError:
        print '\nPlease provide the name of the pwscf output file and natom\n'
        sys.exit(0)

    coords = commands.getoutput('grep -A'+str(natom)+' ATOMIC_POSITIONS '+fname)
    coords = coords.split('ATOMIC_POSITIONS (alat)')
    coords.pop(0)

    out = open(fname.split('.out')[0]+'.xyz','w')

    for i in range(len(coords)):
        out.write(str(natom)+'\n')
        out.write(str(i)+'\n')
        cur_coords = coords[i].split('\n')
        for cur_coord in cur_coords:
            if len(cur_coord.split()) == 4:
                out.write(cur_coord+'\n')

    out.close()


if __name__ == '__main__':
    main()