#!/usr/bin/env python
"""
Test to make sure fft works properly.
Done my fft'ing an fft to make sure it
is the same as the original data set.
"""
import os, sys, Numeric
from FFT import fft

def main():

    N = 3000
    x = Numeric.arange(N,typecode=Numeric.Float)*0.01
    y = Numeric.sin(3.0*x) + 4*Numeric.cos(x)
        
    # Calculate the FFT
    H = fft(y)
    h = fft(H).real / N
    psd = (H*Numeric.conjugate(y)).real / N
    dt = x[1]-x[0]
    

    # Write FFT_inverse to file
    out = open('fft.dat','w')
    for i in range(len(x)):
        out.write(str(x[i])+' '+str(y[i])+' '+str(psd[i])+' '+str(h[len(h)-i-1])+'\n')

    out.close()


if __name__ == '__main__':
    main()
