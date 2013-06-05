#!/usr/bin/env python
"""
Integrate a PSD spectrum of a VACF within a chosen wavenumber range.
"""
import os, sys, commands, glob
import numpy
from scipy import integrate

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        w_min = float(sys.argv[2])
        w_max = float(sys.argv[3])
    except:
        print '\nusage: '+sys.argv[0]+'  psd.dat  wavenumber_min  wavenumber_max\n'
        sys.exit(0)

    # Read in PSD spectrum data
    lines = f.readlines()
    w, intensity = [], []
    norm = 0.0
    for line in lines:
        row = line.split()
        if w_min <= float(row[0]) <= w_max:
            w.append(float(row[0]))
            intensity.append(float(row[1]))
        norm += float(row[1])

    # Integrate using two methods
    simps = integrate.simps(numpy.array(intensity)/norm,numpy.array(w))
    trapz = integrate.trapz(numpy.array(intensity)/norm,numpy.array(w))

    # Report results
    print '\nIntegrating the spectrum in '+sys.argv[1]
    print 'from',w_min,'cm^-1 to',w_max,'cm^-1 gives:'
    print '    ',simps,' using scipy.integrate.simps'
    print '    ',trapz,' using scipy.integrate.trapz\n'


if __name__ == '__main__':
    main()
