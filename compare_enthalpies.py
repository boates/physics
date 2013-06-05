#!/usr/bin/env python
"""
Compare enthalpies of different structures over a pressure range
"""
import os, sys, commands, glob, numpy

def main():

    # Retrieve user input
    try:
        spgrp = sys.argv[1]
    except:
        print '\n usage: '+sys.argv[0]+' space_group_to_compare_with \n'
        sys.exit(0)

    # Get the data
    os.system('get_enthalpies.py')
    p1r = commands.getoutput("awk '{print $1}' enthalpy.dat").split()
    h1r = commands.getoutput("awk '{print $2}' enthalpy.dat").split()
    p2 = commands.getoutput("awk '{print $1}' ../"+spgrp+"/enthalpy.dat").split()
    h2 = commands.getoutput("awk '{print $2}' ../"+spgrp+"/enthalpy.dat").split()

    # Sort with pressure
    for i in range(len(p1r)):
        p1r[i] = int(p1r[i])
    args = numpy.argsort(p1r)
    p1, h1 = [], []
    for a in args:
        p1.append(str(p1r[a]))
        h1.append(h1r[a])

    # Write the data to file
    out = open('dH.dat','w')
    out.write('# dH (with respct to '+spgrp+') meV/atom\n')
    for i in range(len(p1)):
        if p1[i] in p2:
            dh = ( float(h1[i]) - float(h2[p2.index(p1[i])]) )*1000.0
            out.write(p1[i]+' '+str(dh)+'\n')
    out.close()

if __name__ == '__main__':
    main()
