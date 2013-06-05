#!/usr/bin/env python
"""
linear_interpolation_point.py
Author: Brian Boates

Return linearly interpolated value for y for
a given x value (input file xy.dat)
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        f.close()
        x0 = float(sys.argv[2])
    except:
        print '\n usage: '+sys.argv[0]+'  xy.dat  x0 \n'
        sys.exit(1)

    # read data in from file
    x, y = [], []
    for line in lines:
        row = line.split()
        x.append(float(row[0]))
        y.append(float(row[1]))

    # loop over data
    found = False
    for i in range(1, len(x)):

        # find the correct interpolation interval
        if x[i-1] < x0 < x[i] or x[i] < x0 < x[i-1]:

            found = True
            
            x1 = x[i-1]
            x2 = x[i]
            y1 = y[i-1]
            y2 = y[i]

            # linearly interpolate for y0
            dx = x2 - x1
            dy = y2 - y1
            w = (x0 - x1) / dx
            y0 = y1 + w*dy

            print x0, y0

    if found == False:
        print '\n x[i] < x0 < x[i-1] never satisifed for this x0 and data set \n'


if __name__ == '__main__':
    main()
