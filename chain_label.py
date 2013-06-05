#!/usr/bin/env python
"""
Read in a CHAIN.cbn file output from chain_angles.py
and create a chain_label.dat file for visualization in VMD.
"""
import os, sys, commands, glob

def main():

    try:
        f = open(sys.argv[1],'r')
        lines = int(commands.getoutput('wc -l '+sys.argv[1]).split()[0]) - 10
    except:
        print '\n usage '+sys.argv[0]+' CHAIN.cbn\n'
        sys.exit(0)

    # Read cbn file header
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    natom = int(f.readline().split()[-1])
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()

    out = open('chain_label.dat','w')
    steps = lines/natom
    k = 1
    for i in range(steps):

        labels = []
        
        # Monitor progress
        if float(k)/steps*100 % 10.0  == 0:
            print int(float(k)/steps*100),'%'
        k += 1

        for j in range(natom):
            
            nn = int(f.readline().split()[-1])

            # label = 0 for molecule, label = 1 for chain
            if nn >= 0:
                labels.append('0')
            if nn == -1:
                labels.append('1')
        
        for label in labels:
            out.write(label+' ')
        out.write('\n')

    out.close()

if __name__ == '__main__':
    main()
