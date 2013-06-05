#!/usr/bin/env python
"""
Samples several configs from a wrapped.xyz file and creates directories and
runs SCF VASP calculations for many generated POSCAR files
"""
import os, sys, commands, glob
import random

def main():

    try:
        f = open('wrapped.xyz','r')
        f.close()
        alat = float(sys.argv[1])
        N = int(sys.argv[2])
    except:
        print '\nusage:   '+sys.argv[0]+'  alat(ang) number_configs_to_sample'
        print '\nMake sure the wrapped.xyz file is present.\n'
        sys.exit(0)

    natom = int(commands.getoutput('head -1 wrapped.xyz').split()[0])
    nConfig = int(commands.getoutput('wc -l wrapped.xyz').split()[0]) / (natom+2)

    for i in range(N):

        z = ''
        for j in range( 4 - len(str(i+1)) ):
            z += '0'
        dir = z+str(i+1)

        S = int(random.random() * nConfig)
        if nConfig > 2000:
            if S < 1000:
                S += 1000

        # File and data preparation
        os.system('mkdir '+dir)
        os.system('select_snapshot.py wrapped.xyz '+str(S))
        os.system('POSCAR_from_xyz_CO2.py snapshot.xyz '+str(alat)+' y')
        os.system('nn_distance_distributions_cnn.py 200 '+str(alat/0.529177)+' ' \
                                  +str(alat/0.529177)+' '+str(alat/0.529177)+' snapshot.xyz > nn.out')
        os.system('rm -f nn.out')
        os.system('mv -f nn_average.hist TRAJEC.cnn snapshot.xyz POSCAR '+dir)
        os.system('cp /work2/sbonev/boates/carbon_dioxide/vasp/bader/INPUTS/* '+dir)

        # Execute vasp
        cwd = os.getcwd()     
        os.chdir(dir)
        os.system('qsub bader.irm.s')
        os.chdir(cwd)
        

if __name__ == '__main__':
    main()
