#!/usr/bin/env python
"""
Loop enthalpy_from_vasp_relax.py over several numbered directories
"""
import os, sys, commands, glob

def main():

    # Detect present directories
    dirs = glob.glob('0*')
    if len(dirs) == 0:
        print '\n No directories of the form 0*/ detected - exiting... \n'
        sys.exit(0)

    cwd = os.getcwd()
    out = open('enthalpy.dat','w')
    for d in dirs:
        os.chdir(d)
        if 'OUTCAR' in glob.glob('*'):
            try:
                end = commands.getoutput('tail -n-1 OUTCAR').split()[0]
                if 'Voluntary' == end:
                    H = commands.getoutput('enthalpy_from_vasp_relax.py').split()[-2]
                    out.write(d+' '+H+'\n')
                    out.flush()
            except: pass
        os.chdir(cwd)
    out.close()


if __name__ == '__main__':
    main()
