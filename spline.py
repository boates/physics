#!/usr/bin/env python
"""
Take a data set and spline it to a 
more dense set of points
"""
import os, sys, commands, glob, numpy
from scipy import interpolate

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        f.close()
        col1 = int(sys.argv[2])
        col2 = int(sys.argv[3])
        density = int(sys.argv[4])
    except:
        print '\n usage: '+sys.argv[0]+'  fname  col1  col2  density\n'
        print ' fname: file with data'
        print ' col1: data to be interpolated (x)'
        print ' col2: data to be interpolated (y)'
        print ' density: amount of new points per point (i.e. 10)\n'
        sys.exit(0)

    # Read in the data
    x, y = [], []
    for line in lines:
        x.append(float(line.split()[col1-1]))
        y.append(float(line.split()[col2-1]))

    # Assume constant spacing in x
    dx = x[1] - x[0]

    # Get max/min of x
    xmax = max(x)
    xmin = min(x)

    # Convert to numpy arrays
    x = numpy.array(x)
    y = numpy.array(y)

    # Create the new, more dense, x
    xnew = numpy.arange(xmin,xmax,dx/density)

    # Do the interpolation
    tck  = interpolate.splrep(x,y,s=0.0,k=1) # s=0: no smoothing, k=1: linear spline
    ynew = interpolate.splev(xnew,tck)

    # Write to file
    out = open('spline.dat','w')
    for i in range(len(xnew)):
        out.write(str(xnew[i])+' '+str(ynew[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
