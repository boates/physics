#!/usr/bin/env python
"""
RDF_extend.py
Author: Brian Boates

Take an xyz file and create a supercell using supercell.x
Then calculate the RDF using RDF
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        xyz  = sys.argv[1]
        alat = sys.argv[2]
    except:
        print '\n usage: '+sys.argv[0]+' TRAJEC.xyz alat(Ang)\n'
        sys.exit(0)

    f = open('supercell.in','w')
    f.write(xyz+'\n')
    f.write(alat+', '+alat+', '+alat+'\n')
    f.close()
    os.system('supercell.x < supercell.in > supercell.out')

    f = open('ERDF.in','w')
    f.write('SUPERCELL.xyz\n')
    f.write('ERDF.dat\n')
    f.write('N\nN\n')
    f.write(alat+'\n')
    f.write('0.02\n')
    f.write('0\n')
    f.write(str(float(alat)*2.0)+', '+str(float(alat)*2.0)+', '+str(float(alat)*2.0)+'\n')
    f.close()
    os.system('RDF < ERDF.in > ERDF.out')

    os.system('rm -f supercell.out ERDF.out')
    

if __name__ == '__main__':
    main()
