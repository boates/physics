#!/usr/bin/env python

""" Script to covnvert .xyz files to DL_POLY CONFIG files """

import os, sys, commands

def main():

    """ Create CONFIG file """

    try:
        fname_xyz = sys.argv[1]
        acell = float(sys.argv[2])
    except IndexError:
        print 'Please provide name of .xyz file and lattice constant (angstrom)'
        sys.exit(0)

    file_xyz = open(fname_xyz,'r')
    natom = int( file_xyz.readline().strip() )
    tstep = int( file_xyz.readline().strip() )

    out = open('CONFIG','w')
    out.write('DL_POLY: CONFIG file \n')
    out.write('\t'+'0'+'\t'+'1'+'\n')
    out.write('\t'+str(acell)+'\t'+'0.00000000'+'\t'+'0.00000000'+'\n')
    out.write('\t'+'0.00000000'+'\t'+str(acell)+'\t'+'0.00000000'+'\n')
    out.write('\t'+'0.00000000'+'\t'+'0.00000000'+'\t'+str(acell)+'\n')

    lines = file_xyz.readlines()
    for i in range(len(lines)):

        row = lines[i].split()
        typat = row[0]
        x = float(row[1])
        y = float(row[2])
        z = float(row[3])

        out.write(typat+'\t'+str(i+1)+'\n')
        out.write('\t'+str(x*acell)+'\t'+str(y*acell)+'\t'+str(z*acell)+'\n')

    out.close()

    os.system('config_spacer.py')
    os.system('mv CONFIG_fixed CONFIG')

if __name__ == '__main__':
    main()
