#!/usr/bin/env python
"""
Loop over available pressure directories and create a
pressure vs. enthalpy per atom data file
"""
import os, sys, commands, glob, commands

def main():

    # Check for data
    dirs = glob.glob('*GPa')
    if len(dirs) == 0:
        print '\n No data detected. \n'
        sys.exit(0)

    cwd = os.getcwd()
    out = open('enthalpy.dat','w')
    for d in dirs:
        os.chdir(d)
        if 'OUTCAR' in glob.glob('*'):
            end = commands.getoutput('tail -n-1 OUTCAR').split()
            if 'Voluntary' in end:
                H = commands.getoutput('enthalpy_from_vasp_relax.py').split()[-2]
                out.write(d.rstrip('GPa')+' '+H+'\n')
                out.flush()
        os.chdir(cwd)
    out.close()
        

if __name__ == '__main__':
    main()
