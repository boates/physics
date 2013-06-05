#!/usr/bin/env python
"""
Grab a given snapshot of an xyz file
"""
import os, sys, commands, glob

def main():

    try:
        f = sys.argv[1]
        natom = int(commands.getoutput('head -1 '+f))
        Nlines = int(commands.getoutput("wc -l "+sys.argv[1]+" | awk '{print $1}'"))
        print '\nThere are '+str(Nlines/(natom+2))+' snapshots in this file.\n'
        snapshot = int(sys.argv[2])
    except IndexError:
        print '\nusage: '+sys.argv[0]+' TRAJEC.xyz which_snapshot\n'
        sys.exit(0)

    os.system('tail -'+str( Nlines - (snapshot-1)*(natom+2) )+' '+f+' | head -'+str(natom+2)+' > snapshot.xyz')

if __name__ == '__main__':
    main()
