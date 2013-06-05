#!/usr/bin/env python
"""
HISTORY_select_snapshot
Author: Brian Boates

Select an arbitrary snapshot from a DL_POYL HISTORY file in CONFIG file form
"""
import os, sys, glob, commands

def main():

    # Retrieve user input
    try:
        fhistory = sys.argv[1]
        tstep = int(sys.argv[2])
    except:
        print '\n usage: '+sys.argv[0]+' HISTORY tstep \n'
        sys.exit(1)

    natom  = int(commands.getoutput("head -2 "+fhistory+" | tail -n-1 | awk '{print $3}'"))
    levcfg = int(commands.getoutput("head -2 "+fhistory+" | tail -n-1 | awk '{print $1}'"))

    fname = 'CONFIG.'+str(tstep)
    os.system("head -2 "+fhistory+" > "+fname)
    os.system("head -"+str((natom*(levcfg+2)+4)*(tstep+1)+2)+" "+fhistory+" | tail -"+str(natom*(levcfg+2)+3)+" >> "+fname)


if __name__ == '__main__':
    main()

