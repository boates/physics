#!/usr/bin/env python

""" Fit a line using lregress.py module to enthalpy curves """

import os, sys
sys.path.append('/home/boates/software')
import lregress

def main():
    """
    Read in the data and fit the line
    Print out a,b,dy,da,db to screen
    """
    try:
        f = open(sys.argv[1])
    except:
        print '\nPlease give name of file with data to be fit (in columns: X Y)\n'
        sys.exit(0)

    header = f.readline()
    lines = f.readlines()
    f.close()
    X = []
    Y = []
    for line in lines:
        X.append(float(line.split()[0]))
        Y.append(float(line.split()[1]))

    a,b,dy,da,db = lregress.lregress(X,Y)

    print '\na =',a,'\tb =',b,'\nda =',da,'\tdb =',db,'\tdy =',dy
    print '\nY_fit = '+str(a)+' + '+str(b)+' * X\n'

if __name__ == '__main__':
    main()