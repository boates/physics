#!/usr/bin/env python
"""
3d_surface.py
Author: Brian Boates

Create a 3D surface color plot from data
"""
import os, sys, commands, glob, numpy
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        f.close()
        if len(sys.argv) >= 2:
            ybeg = 0
            yend = 1
        if len(sys.argv) >= 3:
            ybeg = int(sys.argv[4])
        if len(sys.argv) >= 4:
            yend = int(sys.argv[5])

    except:
        print '\n usage: '+sys.argv[0]+' xyz.dat ybeg(~70) yend(~250)\n'
        print 'x is the temperature'
        print 'y is the distance'
        print 'z is the potential'
        print 'ybeg / yend are the number of x data points to omit at the beginning / end\n'
        sys.exit(0)

    # Read in xyz data
    k = -1
    x0 = None #float(lines[0].split()[0])
    xa, ya, za = [], [], []
    for i in range(len(lines)):
        line = lines[i].split()
        x = float(line[0])
        y = float(line[1])
        z = float(line[2])
        if x != x0:
            k += 1
            xa.append( x )
            ya.append( [] )
            za.append( [] )
            x0 = x
        ya[k].append( y )
        za[k].append( z )

    # Convert lists to arrays
    x = numpy.array( xa )
    ytails = []
    for i in range(len(x)):
        ytails.append(ya[i][-1])
        ya[i] = numpy.array(ya[i])
        za[i] = numpy.array(za[i])
    y = numpy.array( ya )
    z = numpy.array( za )

    # Interpolate all y vs. z data to same y grid
    dy = y[0][1] - y[0][0]
    ymin = numpy.min(y)
    ymax = min(ytails)   # Length for all datasets is truncated to the shortest one
    y_new = numpy.arange(ymin,ymax,dy/2.)

    z_new = []
    for i in range(len(x)):
        tck = interpolate.splrep(y[i],z[i],s=0,k=1)   # s=0: no smoothing, k=1: linear spline
        z_new.append( interpolate.splev(y_new,tck) )
    y = y_new
    z = numpy.array( z_new )

#    print numpy.shape(x), numpy.shape(y), numpy.shape(z)
        
    # Create the x/y meshgrid
    x, y = numpy.meshgrid(x, y)

    # Create the figure
    fig = plt.figure()
    ax = Axes3D(fig)

#    ax.set_xlabel('r [A]')
#    ax.set_ylabel('T [K]')
#    ax.set_zlabel('V(T,r) [eV]')

    ax.set_xlabel('V [A^3/atom]')
    ax.set_ylabel('T [K]')
    ax.set_zlabel('P(T,V) [GPa]')

    reverse = False
    if reverse:
        x = x[::-1]
        z = numpy.transpose(z[::-1])
    else:
        z = numpy.transpose(z)

    ax.plot_surface(y, x, z, rstride=1, cstride=1, cmap=cm.hsv) # hsv

    plt.savefig('3d_surface.png')


if __name__ == '__main__':
    main()
