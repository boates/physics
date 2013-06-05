#!/usr/bin/env python
"""
A script to check for incomplete configurations in an XDATCAR file
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        for i in range(5): f.readline()
        natom = int(sys.argv[2])
    except:
        print '\n usage: '+sys.argv[0]+' XDATCAR natom\n'
        sys.exit(0)

    # Check to see if file is ok
    lines = int(commands.getoutput("wc -l "+sys.argv[1]).split()[0])
    if (lines - 5.0) % (natom+1) == 0.0:
        print '\n XDATCAR file is ok - exiting...\n'
        sys.exit(0)

    # If the file is not ok
    for i in range(int(lines/(natom+1.0))):
        j = 0
        k = f.readline()
        if "Konfig" not in k and k.strip() != "":
            print i, j, k, 'f.readline() did not return a blank line as expected'
            sys.exit(0)
        for j in range(natom):
            xyz = f.readline().split()
            if len(xyz) != 3:
                print i, j, 'f.readline() did not return x, y, z as expected'
                sys.exit(0)

    f.close()



if __name__ == '__main__':
    main()
