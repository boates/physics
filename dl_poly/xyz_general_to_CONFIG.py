#!/usr/bin/env python
"""
Script to covnvert .xyz files to DL_POLY CONFIG files
"""
import os, sys, commands

def main():
    """
    Create CONFIG file
    """
    try:
        fname_xyz = sys.argv[1]
        a = sys.argv[2]
        b = sys.argv[3]
        c = sys.argv[4]
    except IndexError:
        print 'Please provide name of .xyz file and lattice consts, a,b,c'
        sys.exit(0)

    file_xyz = open(fname_xyz,'r')
    natom = int( file_xyz.readline().strip() )
    tstep = int( file_xyz.readline().strip() )

    out = open('CONFIG','w')
    out.write('DL_POLY: CONFIG file \n')
    out.write('\t'+'0'+'\t'+'1'+'\n')
    out.write('\t'+a+'\t'+'0.00000000'+'\t'+'0.00000000'+'\n')
    out.write('\t'+'0.00000000'+'\t'+b+'\t'+'0.00000000'+'\n')
    out.write('\t'+'0.00000000'+'\t'+'0.00000000'+'\t'+c+'\n')

    lines = file_xyz.readlines()
    for i in range(len(lines)):

        row = lines[i].split()
        typat = row[0]

        out.write(typat+'\t'+str(i+1)+'\n')
        out.write('\t'+row[1]+'\t'+row[2]+'\t'+row[3]+'\n')

    file_xyz.close()

if __name__ == '__main__':
    main()
