#!/usr/bin/env python
"""
Create a gaussian data set and its fourier transform
"""
import os, sys, Numeric, math

def main():

    # Create gaussian data set
    a = 50.0
    x = Numeric.arange(-3,3,0.001,typecode=Numeric.Float)
    y = Numeric.exp(-a*x**2.0)

    # Analytic solution for fourier transform of above gaussian
    N = float(len(x))
    dx = x[1] - x[0]
    T = N * dx
    df = 1.0 / T
#    f  = Numeric.arange(N,typecode=Numeric.Float)*df
#    Y = (math.pi/a)**0.5 * Numeric.exp(-4.0 * f**2.0 / a)    

    # Write to file
    out = open('gaussian.dat','w')
    for i in range(len(x)):
        out.write(str(x[i])+' '+str(y[i])+'\n')#' '+str(f[i])+' '+str(Y[i])+'\n')
    out.close()
#    f = open('ftgaussian.dat','w')
#    for i in range(len(x)):
#        f.write(str(f[i])+' '+str(Y[i])+'\n')
#    f.close()


if __name__ == '__main__':
    main()
