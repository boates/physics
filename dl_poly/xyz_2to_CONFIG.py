#!/usr/bin/env python

""" Script to convert .xyz files to DL_POLY CONFIG files """

import os, sys, commands

def main():

    """ Create CONFIG file """

    try:
        fname_xyz = sys.argv[1]
        ax = sys.argv[2]
        ay = sys.argv[3]
        az = sys.argv[4]
    except IndexError:
        print '\nPlease provide name of .xyz file & lattice consts: ax,ay,az\n'
        sys.exit(0)

    file_xyz = open(fname_xyz,'r')
    natom = int( file_xyz.readline().strip() )
    tstep = int( file_xyz.readline().strip() )

    out = open('CONFIG','w')
    out.write('DL_POLY: CONFIG file \n')
    out.write('\t'+'0'+'\t'+'1'+'\n')
    out.write('\t'+ax+'\t'+'0.00000000'+'\t'+'0.00000000'+'\n')
    out.write('\t'+'0.00000000'+'\t'+ay+'\t'+'0.00000000'+'\n')
    out.write('\t'+'0.00000000'+'\t'+'0.00000000'+'\t'+az+'\n')

    lines = file_xyz.readlines()
    for i in range(len(lines)):

        row = lines[i].split()
        typat = row[0]

        out.write(typat+'\t'+str(i+1)+'\n')
        out.write('\t'+row[1]+'\t'+row[2]+'\t'+row[3]+'\n')

if __name__ == '__main__':
    main()
