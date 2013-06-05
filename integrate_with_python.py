#!/usr/bin/env python
"""
Integrate using scipy.integrate
"""
import os, sys, commands
import numpy
from scipy import integrate

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
    except:
        print '\n usage: '+sys.argv[0]+' file_to_integrate_first_two_colums_XY_no_header.dat \n'
        sys.exit(0)

    x, y = [], []
    for line in lines:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))

    x = numpy.array(x)
    y = numpy.array(y)

    trapz = integrate.trapz(y,x)
    simps = integrate.simps(y,x)

    print 'using trapz: ',trapz,'\nusing simps: ',simps


if __name__ == '__main__':
    main()
