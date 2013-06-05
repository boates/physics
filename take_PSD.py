#!/usr/bin/env python
"""
Take FFT of data set
"""
import os, sys, Numeric
import FFT

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        f.close()
    except:
        print '\n usage: '+sys.argv[0]+' XY.dat\n'
        sys.exit(0)

    # Parse the data
    t, y = [], []
    for line in lines:
        t.append(float(line.split()[0]))
        y.append(float(line.split()[1]))

    # Calculate the FFT and get the frequencies
    N  = float(len(t))
    dt = t[1] - t[0]
    T  = N * dt
    df = 1.0 / T
    f  = Numeric.arange(N,typecode=Numeric.Float)*df
    H  = ( FFT.fft(y)*Numeric.conjugate(FFT.fft(y)) ).real / N

    # Write to file
    out = open('PSD.dat','w')
    for i in range(len(f)/2):
        out.write(str(f[i])+' '+str(H[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
