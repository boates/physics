#!/usr/bin/env python
"""
FFT testing script
"""
import os, sys, commands
import FFT, Numeric

def main():

    X = Numeric.arange(0,20*Numeric.pi,0.01,typecode=Numeric.Float)
    Y1 = Numeric.sin(X)
    Y2 = Numeric.sin(2*X)

    out = open('fncs.dat','w')
    for i in range(len(X)):
        out.write(str(X[i])+' '+str(Y1[i])+' '+str(Y2[i])+'\n')
    out.close()

    # Take FFT's
    FY1  = FFT.fft(Y1).real
    FY2  = FFT.fft(Y2).real
    N  = float(len(X)) # get # of data points
    dX = X[1]-X[0] # get distance spacing
    T  = N*dX  # define the period (total time)
    dQ = 1./T # define frequency step
    Q=Numeric.arange(N,typecode=Numeric.Float)*dQ

    print Q[list(FY1).index(max(FY1))], max(FY1)
    print Q[list(FY2).index(max(FY2))], max(FY2)
    print list(FY1).index(max(FY1))
    print list(FY2).index(max(FY2))

    out = open('FFTs.dat','w')
    for i in range(len(Q)):
        out.write(str(Q[i])+' '+str(FY1[i])+' '+str(FY2[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
