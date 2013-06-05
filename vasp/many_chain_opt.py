#!/usr/bin/env python
"""
Generate random lattices of nitrogen chains and submit for
optimization in VASP
"""
import os, sys, commands, glob
import random

def main():

    try:
        input = sys.argv[1]
        N = int(sys.argv[2])
    except:
        print '\n usage:   '+sys.argv[0]+'  INPUT_file_for_generator  number_of_lattices_to_generate\n'
        sys.exit(0)

    os.system('rm -f chains.xyz chains.log chains.POSCAR')
    for i in range(N):

        z = ''
        for j in range( 4 - len(str(i+1)) ):
            z += '0'
        dir = z+str(i+1)

        # File and data preparation
        os.system('cp -r template/ '+dir)
        k = -1
        while k < 0:
            os.system('generate_random_chain_cell.py '+input+' > chains.out')
            if len(glob.glob('chains.xyz')) == 1:
                k = 1
        os.system('mv -f chains.POSCAR '+dir+'/POSCAR')
        os.system('mv -f chains.log '+dir)
        os.system('rm -f chains.xyz chains.out')

        # Execute vasp
        cwd = os.getcwd()     
        os.chdir(dir)
        os.system('qsub vasp.irm.s')
        os.chdir(cwd)
        

if __name__ == '__main__':
    main()
