#!/usr/bin/env python
"""
Extract (Fx,Fy,Fz) from an OUTCAR file and write to a FORCES.xyz file
"""
import os, sys, commands

def main():

    try:
        f = open('OUTCAR','r')
        f.close()
    except:
        print '\nMake sure OUTCAR file is present\n'
        sys.exit(0)

    natom = int(commands.getoutput("grep NIONS OUTCAR | tail -1 | awk '{print $12}'"))
    os.system("grep -A"+str(natom+1)+" 'TOTAL-FORCE' OUTCAR | grep -v '\-\-' | grep -v 'TOTAL-FORCE' | awk '{print $4,$5,$6}' > FORCES.xyz")

if __name__ == '__main__':
    main()
