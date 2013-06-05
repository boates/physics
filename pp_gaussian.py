#!/usr/bin/env python
"""
Create a gaussian to 'correct' a pair potential
"""
import os, sys, numpy, math

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        f.close()
        peak = float(sys.argv[2])
        a = float(sys.argv[3])
        b = float(sys.argv[4])
    except:
        print '\n usage: '+sys.argv[0]+'  Poten_XX.dat  peak  a  b  ( from a*e^(-b*(x-peak)^2) )\n'
        sys.exit(1)

    # Get r and V values from pair potential data file
    r, V = [], []
    for line in lines:
        r.append(float(line.split()[0]))
        V.append(float(line.split()[1]))        
    r = numpy.array(r)
    V = numpy.array(V)

    # Create gaussian correction function
    e = numpy.e
    gaussian = a*e**(-b*(r-peak)**2)

    # Create new pair potential as original - gaussian
    out = open('pp.new','w')
    for i in range(len(r)):
        out.write(str(r[i])+' '+str(V[i]-gaussian[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
