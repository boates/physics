#!/usr/bin/env python
"""
Calculate the standard deviation of a given column from a data file.
"""
import os, sys, commands
import Numeric

def stdev(X):
    """
    Calculate standard deviation of array X.
    """
    mean = float(Numeric.average(X))
    sd = 0.0
    for x in X:
        sd += (float(x) - mean)**2
    sd = (sd/len(X))**0.5

    return sd

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        c = int(sys.argv[2])
    except:
        print '\n usage: '+sys.argv[0]+' f.dat column\n'
        sys.exit(0)

    # Read in data from file
    X = []
    for line in lines:
        X.append(float(line.split()[c-1]))

    # Calculate stdev
    sd = stdev(X)

    # Print results
    print '\nStandard deviation of column='+str(c)+' in '+sys.argv[1]+' is:',sd,'\n'

if __name__ == '__main__':
    main()
