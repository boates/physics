#!/usr/bin/env python
"""
dos_avg_vasp.py
Author: Brian Boates

Average multiple DOS calculations into one dataset

Interpolate each DOS file to ensure the averaging
takes place at same energies.
"""
import os, sys, commands, glob, numpy
from scipy import interpolate

def main():

    # Print usage warning to the user
    print '\n usage: '+sys.argv[0]+'\n'
    print 'Assumes files */nscf/DOSCAR are present \n'

    # Create the necessary DOSCAR.dat files
    dirs = glob.glob('*/nscf/')
    cwd = os.getcwd()
    for d in dirs:
        os.chdir(d)
        os.system('DOSCAR_plottable.py DOSCAR')
        os.chdir(cwd)

    # Grab list of DOS files to be averaged
    DOSs = glob.glob('*/nscf/DOSCAR.dat')
    nDOS = len(DOSs)

    # Retrieve data from detected files
    E, g = [], []
    for dos in DOSs:
        E.append(commands.getoutput("awk '{print $1}' "+dos).split())   # in eV
        g.append(commands.getoutput("awk '{print $2}' "+dos).split())

    # Convert all types to floats, and lists to arrays
    Emins, Emaxs = [], []
    for i in range(nDOS):
        for j in range(len(E[i])):
            E[i][j] = float(E[i][j])
            g[i][j] = float(g[i][j])
        Emins.append(min(E[i]))
        Emaxs.append(max(E[i]))
    for i in range(nDOS):
        E[i] = numpy.array(E[i])
        g[i] = numpy.array(g[i])

    # Assume constant energy spacing (i.e. E[i+1] - E[i] = const, for all i)
    dE = E[0][1] - E[0][0]

    # Select max and min energies
    Emin = max(Emins)
    Emax = min(Emaxs)

    # Create the new energy array using 10x the resolution of the original
    E_new = numpy.arange(Emin,Emax,dE/10.)

    # Do the interpolation
    g_new = []
    for i in range(nDOS):
        tck = interpolate.splrep(E[i],g[i],s=0,k=1)   # s=0: no smoothing, k=1: linear spline
        g_new.append( interpolate.splev(E_new,tck) )
    g_new = numpy.array(g_new)

    # Do the averaging (along appropriate axis)
    g_avg = numpy.mean(g_new,axis=0)

    # Write to file
    out = open('dos_avg.dat','w')
    out.write('# E (eV), DOS averaged over '+str(nDOS)+' configurations\n')
    for i in range(len(E_new)):
        out.write(str(E_new[i])+' '+str(g_new[0][i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
