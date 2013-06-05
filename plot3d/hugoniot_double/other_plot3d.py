#!/usr/bin/env python
"""
plot3d.py
Author: Brian Boates

Make 3d contourf or surface plot
of data.XYZ files

interpolates Y and Z
slightly hackish, but it works.
"""
import os, sys, commands, glob, numpy
import matplotlib.pyplot as plt
from scipy import interpolate

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


def main():

    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        f.close()
    except IndexError:
        print '\nusage: '+sys.argv[0]+'  data.XYZ\n'
        sys.exit(0)

    X, Y, Z, count = [], [], [], -1
    for line in lines:
        row = line.split()
        if float(row[0]) not in X:
            X.append(float(row[0]))
            Y.append([])
            Z.append([])
            count += 1
        Y[count].append(float(row[1]))
        Z[count].append(float(row[2]))
    Y = Y[0]
 
    fig = plt.figure(figsize=(12,6),facecolor='w',edgecolor='k')
    plt.contour(Y,X,Z,100,linewidth=0.05,cmap=cm.jet)
    plt.contourf(Y,X,Z,100,cmap=cm.jet)
    
    plt.xlim( (min(Y),max(Y)) )
    plt.ylim( (2000,max(X)) )
    plt.colorbar()
#        plt.colorbar(ticks=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],drawedges=False)

    plt.xlabel('Pressure [GPa]',fontsize=18)
    plt.ylabel('Temperature [K]',fontsize=18)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    plt.savefig('contourf.pdf')

    
if __name__ == '__main__':
    main()
