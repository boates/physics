#!/usr/bin/env python
"""
Remove directories of random relaxations that
did not finish properly
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        natom = sys.argv[1]
    except:
        print '\n usage: '+sys.argv[0]+' natom\n'
        sys.exit(0)

    cwd = os.getcwd()
    os.chdir(natom+'_atoms')
    tails = commands.getoutput('tail -n-1 */OUTCAR').split("==>")
    print 'Directories removed:'
    for t in tails:
        d = t.split('/')[0]
        if 'Voluntary' not in t.split():
            os.system('rm -rf '+d)
            print d


if __name__ == '__main__':
    main()
