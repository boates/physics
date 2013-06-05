#!/usr/bin/env python
"""
Calculate the Power Spectral Density (PSD) of a data set
"""
import os, sys, Numeric
from FFT import fft

def psd(t,y):
    """
    Calculate and return frequencies and PSD of y(t)
    """
    N  = float(len(t))   # GET # OF DATA POINTS
    dt = t[1] - t[0]     # GET TIME INTERVAL
    T  = N * dt          # DEFINE THE PERIOD (TOTAL TIME)
    df = 1.0 / T         # DEFINE FREQUENCY STEP
    H  = fft(y)          # ,n=256) ADDITIONAL OPTION

    # Caculate frequencies and PSD
    f = Numeric.arange(N,typecode=Numeric.Float)*df
    PSD = (Numeric.conjugate(H)*H).real / N

    return f, PSD

def main():

    # Retrieve user input
    try:
        fin = sys.argv[1]
	f = open(fin,'r')
    except IndexError:
        print '\nusage '+sys.argv[0]+' file.dat\n'
	sys.exit(0)

    # Read in data from file
    lines = f.readlines()
    f.close()
    t, y = [], []
    for line in lines:
        row = line.split()
        t.append(float(row[0]))
	y.append(float(row[1]))

    # Calculate the PSD and write to file
    f, PSD = psd(t,y)
    out = open('PSD.dat','w')
    for i in range(len(f)/2):
        out.write(str(f[i])+'    '+str(PSD[i])+'\n')

    out.close()


if __name__ == '__main__':
    main()
