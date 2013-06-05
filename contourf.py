#!/usr/bin/env python
"""
Make contourf plot of TPZ.dat files

interpolates P and Z
"""
import os, sys, commands, glob
import pylab, numpy
from scipy import interpolate

def main():

    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        f.close()
    except IndexError:
        print '\nusage: '+sys.argv[0]+'  TPZ.dat\n'
        sys.exit(0)

    T, P, z, count = [], [], [], -1
    for line in lines:
        row = line.split()
        if float(row[0]) not in T:
            T.append(float(row[0]))
            P.append([])
            z.append([])
            count += 1
        P[count].append(float(row[1]))
        z[count].append(float(row[2]))

    # Interpolate
    dP = 10.0
    Pmin = 0.0
    Pmax = 200.0
    P_new = numpy.arange(Pmin,Pmax,dP/10.)

    z_new = []
    for i in range(len(z)):
        tck = interpolate.splrep(P[i],z[i],s=0,k=1)   # s=0: no smoothing, k=1: linear spline
        z_new.append( interpolate(P_new, tck) )
    z_new = numpy.array(z_new)
    

    fig = pylab.figure(figsize=(12,6),facecolor='w',edgecolor='k')
    pylab.contourf(P_new,T,z_new,10,linewidth=0.05)
    pylab.xlim( (max(y),1.22) )
    pylab.colorbar()
    pylab.savefig('test.png')

if __name__ == '__main__':
    main()
