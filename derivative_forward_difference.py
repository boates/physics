#!/usr/bin/env python
"""
derivative_forward_difference.py
Author: Brian Boates

Calculate dy/dx for a file with x and y as columns 1 and
2, respectively using the forward difference method.
"""
import os, sys, commands, glob
import numpy

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        if '#' in lines[0]: lines.pop(0)
        f.close()
    except:
        print '\n usage: '+sys.argv[0]+' XY.dat\n'
        sys.exit(0)

    # Read in data
    x, y = [], []
    for line in lines:
        row = line.split()
        x.append(float(row[0]))
        y.append(float(row[1]))

    # Differentiate
    out = open('derivative.dat','w')
    for i in range(len(y)-1):
        dy = y[i+1] - y[i]
        dx = x[i+1] - x[i]
        out.write(str(x[i])+'  '+str(dy/dx)+'\n')
    out.close()


if __name__ == '__main__':
    main()
