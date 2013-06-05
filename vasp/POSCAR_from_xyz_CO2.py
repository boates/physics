#!/usr/bin/env python
"""
Convert an xyz snapshot to a POSCAR file
"""
import os, sys, commands, glob

def main():

    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        f.close()
        alat = sys.argv[2]
        scale = sys.argv[3]
    except IndexError:
        print '\nusage: '+sys.argv[0]+' single_timestep.xyz, lattice_const (angstroms), reduce? (y/n)\n'
        sys.exit(0)

    if scale == 'y':
        pass
    else: alat = '1.00000000'

    natom = lines.pop(0).strip()
    tstep = lines.pop(0).strip()

    out = open('POSCAR','w')
    out.write('N\n')
    out.write('    '+alat+'\n')
    out.write('     1.0000000000000000    0.0000000000000000    0.0000000000000000\n')
    out.write('     0.0000000000000000    1.0000000000000000    0.0000000000000000\n')
    out.write('     0.0000000000000000    0.0000000000000000    1.0000000000000000\n')
    out.write(' 32  64\n')
    out.write('Direct\n')
    for line in lines:
        row = line.split()
        out.write(' '+str(float(row[1])/float(alat))+' '+str(float(row[2])/float(alat))+' '+\
                                                         str(float(row[3])/float(alat))+'\n')
    out.close()

if __name__ == '__main__':
    main()
