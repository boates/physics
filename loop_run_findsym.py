#!/usr/bin/env python
"""
Loop run_findsym.py over several numbered directories
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        tol = sys.argv[1]
    except:
        print '\n usage: '+sys.argv[0]+' tolerance(0.01)\n'
        sys.exit(0)

    # Detect present directories
    dirs = glob.glob('0*')
    if len(dirs) == 0:
        print '\n No directories of the form 0*/ detected - exiting... \n'
        sys.exit(0)

    cwd = os.getcwd()
    for d in dirs:
        os.chdir(d)
        if ('OUTCAR' and 'CONTCAR') in glob.glob('*'):
            try:
                end = commands.getoutput('tail -n-1 OUTCAR').split()[0]
                if 'Voluntary' == end:
                    os.system('run_findsym.py CONTCAR '+tol+'')
            except: pass
        os.chdir(cwd)

    # Grep all of the resulting space groups to file
    os.system("grep 'Space Group' */findsym.out > space_groups_tol"+tol+".dat")

if __name__ == '__main__':
    main()
