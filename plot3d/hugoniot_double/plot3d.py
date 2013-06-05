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
        plot_type = sys.argv[2]
    except IndexError:
        print '\nusage: '+sys.argv[0]+'  data.XYZ  plot_type(surface or contourf)\n'
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

    # Interpolate
    dY   = 0.1    # Pressure grid of 1 GPa
    Ymin = 5.0   # From 0 GPa
    Ymax = 14.0  # To Pmax GPa
    Y_new = numpy.arange(Ymin,Ymax,dY)

#    for i in range(len(Y)):
#        Y[i] = Y[i][::-1]
#        Z[i] = Z[i][::-1]
#        j = 0
#        while j < len(Y[i]):
#            if j > 0 and j < len(Y[i]) - 1:
#                if not Y[i][j-1] < Y[i][j] < Y[i][j+1]:
#                    Y[i].pop(j)
#                    Z[i].pop(j)
#                else: j += 1
#            else: j += 1

    Z_new = []
    for i in range(len(Z)):
        tck = interpolate.splrep(Y[i],Z[i],s=0,k=1)   # s=0: no smoothing, k=1: linear spline
        Z_new.append( interpolate.splev(Y_new, tck) )

    if plot_type == 'contourf':

        fig = plt.figure(figsize=(12,6),facecolor='w',edgecolor='k')
        plt.contour(Y_new,X,Z_new,20,linewidth=0.05,cmap=cm.jet)
        plt.contourf(Y_new,X,Z_new,20,cmap=cm.jet)

#        plt.xlim( (min(Y_new),max(Y_new)) )
#        plt.ylim( (1500,max(X)) )
#        plt.colorbar(ticks=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],drawedges=False)
        plt.colorbar()

#        plt.title(r'made with plot3d.py',fontsize=20)
        plt.xlabel('Pressure [GPa]',fontsize=18)
        plt.ylabel('Temperature [K]',fontsize=18)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)

#        plt.show()
        plt.savefig('contourf.pdf')

    
    elif plot_type == 'surface':

        # Create the x/y meshgrid
        x, y = numpy.meshgrid(X, Y_new)

        # Create the figure
        fig = plt.figure()
        ax = Axes3D(fig)

#        ax.set_xlabel('V [A^3/atom]')
#        ax.set_ylabel('T [K]')
#        ax.set_zlabel('P(T,V) [GPa]')
        ax.set_xlabel('P [GPa]')
        ax.set_ylabel('T [K]')
        ax.set_zlabel('fraction')

        reverse = False
        if reverse:
            x = X[::-1]
            z = numpy.transpose(Z_new[::-1])
        else:
            z = numpy.transpose(Z_new)

        ax.plot_surface(y, x, z, rstride=1, cstride=1, cmap=cm.hsv) # hsv

        plt.savefig('surface.pdf')
        

if __name__ == '__main__':
    main()
