#!/usr/bin/env python

""" Script to convert the XDATCAR output file from VASP into xyz format. """

import os, sys, commands

global FNAME, ATOM, HEADER_LINES
FNAME = 'XDATCAR'
ATOM = 'N'
HEADER_LINES = 5

def main():
    
    """ Write the xyz file """

    try:
        acell = float(sys.argv[1])
        natom = int(sys.argv[2])
    except IndexError:
        print '\nPlease provide the lattice constant and natom.\n'
        sys.exit(0)
    
    f = open(FNAME,'r')
    for i in range(HEADER_LINES): f.readline()
    lines = f.readlines()
    f.close()
    
    out = open(FNAME+'.xyz','w')
    
    ntime = len(lines) / (natom + 1)  # +1 for the blank line in between each timestep
    
    for i in range(ntime):
    
        lines.pop(0) # Remove the blank line
        out.write(str(natom)+'\n')
        out.write(str(i+1)+'\n')
    
        for j in range(natom):
            
            row = lines.pop(0).split()
            for i in range(len(row)): row[i] = float(row[i])*acell
    
            out.write(ATOM+' '+str(row[0])+' '+str(row[1])+' '+str(row[2])+'\n')
    
    out.close()

if __name__ == "__main__":
    main()
