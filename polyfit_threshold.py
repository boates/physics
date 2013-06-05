#!/usr/bin/env python
"""
Polynomial fitter
"""
import os, sys, commands, glob
import Numeric, math
from scipy import optimize
from scipy import linalg

global THRESHOLD
THRESHOLD = 10.0E-10000

def fitter(rs,P,n):

    # Build nth order basis
    x = rs
    y = P

    basis=[]
    for i in range(1,n):
        basis.append(x**i)

    # Basis array transposed
    M = Numeric.transpose(Numeric.array(basis))
    a, resid, rank, sigma = linalg.lstsq(M,y)

    # Build nth order fit
    fit = Numeric.zeros(len(rs),typecode=Numeric.Float)

    for i in range(1,n):
        if abs(a[i-1]) > THRESHOLD:
            fit += a[i-1]*rs**i

    return fit, a

def main():

    # Retrieve data from user and format
    try:
        f = open(sys.argv[1])
        header = f.readline()
        if '#' != header.strip()[0]:
            f.close()
            f = open(sys.argv[1])
        lines = f.readlines()
        X, Y = [], []
        for line in lines:
            X.append(float(line.split()[0]))
            Y.append(float(line.split()[1]))
        X = Numeric.array(X)
        Y = Numeric.array(Y)
        ORDER = int(sys.argv[2])+1
    except:
        print '\n usage: '+sys.argv[0]+' X_Y.dat order_of_polyfit xmin ymin\n'
        sys.exit(0)

    # Get the polynomial fit
    fit_Y, a = fitter(X,Y,ORDER)
    if len(sys.argv) == 5:
        min_X = float(sys.argv[3])
        max_X = float(sys.argv[4])
    else:
        min_X = min(X)
        max_X = max(X)
    stepsize = (max_X - min_X)/1000.
    new_X = Numeric.arange(min_X,max_X,stepsize,typecode=Numeric.Float)
    new_Y = Numeric.zeros(Numeric.shape(new_X),typecode=Numeric.Float)
    for i in range(len(a)):
        if abs(a[i]) > THRESHOLD:
            new_Y += a[i]*new_X**(i+1)

    X, Y = new_X, new_Y

    print "a[i] =", a

    # Write the fit to file
    out = open('polyfit.dat','w')
    for i in range(len(X)):
        out.write(str(X[i])+'  '+str(Y[i])+'\n')
    out.close()
    
if __name__ == '__main__':
    main()
