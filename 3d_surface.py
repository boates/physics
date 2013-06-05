#!/usr/bin/env python
"""
3d_surface.py
Author: Brian Boates

Create a 3D surface color plot from data
"""
import os, sys, commands, glob, numpy

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

def main():

    # Retrieve user input
    try:
        fx = open(sys.argv[1],'r')
        fy = open(sys.argv[2],'r')
        fz = open(sys.argv[3],'r')
        xlines = fx.readlines()
        ylines = fy.readlines()
        zlines = fz.readlines()
        fx.close()
        fy.close()
        fz.close()
        if len(sys.argv) >= 4:
            ybeg = 0
            yend = 1
        if len(sys.argv) >= 5:
            ybeg = int(sys.argv[4])
        if len(sys.argv) >= 6:
            yend = int(sys.argv[5])

    except:
        print '\n usage: '+sys.argv[0]+' x.dat y.dat z.dat ybeg(~70) yend(~250)\n'
        print 'x.dat is the temperature'
        print 'y.dat is the distance'
        print 'z.dat is the potential, each (row) column is for a different (x) y value\n'
        print 'ybeg / yend are the number of x data points to omit at the beginning / end\n'
        sys.exit(0)

    # Read in x data
    x, t = [], []
    for line in xlines:
        x.append( float(line) )
        t.append( float(line) )

    # Read in y data
    y = []
    for line in ylines:
        y.append( float(line) )
    y = y[ybeg:-yend]

    # Read in z data
    z = []
    for i in range(len(y)):
        z.append( [] )
        line = zlines[i+ybeg].split()
        for j in range(len(x)):
            z[i].append( float(line[j]) )
        z[i] = numpy.array( z[i] )
    z = numpy.array( z )

    # To plot the minima curve
    """
    f = open('1.42_min.dat','r')
    lines = f.readlines()
    f.close()
    r, m = [], []
    for line in lines:
        row = line.split()
        r.append(float(row[0]))
        m.append(float(row[1]))
    r = numpy.array(r)
    m = numpy.array(m)
    """

    # Create the x/y meshgrid
    x, y = numpy.meshgrid(x, y)

    # Create the figure
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.set_xlabel('r [A]')
    ax.set_ylabel('T [K]')
    ax.set_zlabel('V(T,r) [eV]')
    
    ax.plot_surface(y, x, z, rstride=1, cstride=1, cmap=cm.hsv) # hsv

#    ax.plot(r, t, m)

    plt.savefig('3d_surface.png')


if __name__ == '__main__':
    main()
