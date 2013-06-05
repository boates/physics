#!/usr/bin/env python
"""
Integrate a neighbor distribution curve from an nn_average.hist file
"""
import os, sys, commands, glob
import Numeric
from scipy import integrate

def main():

    # Retrieve user input
    try:
        neighbor = int(sys.argv[1])
        fin = 'nn_average.hist'
        f = open(fin,'r')
        f.close()
        os.system('running_avg.py '+fin+' '+str(neighbor+1)+' 4')
        f = open('avg.dat','r')
        lines = f.readlines()
        f.close()
    except:
        print '\n usage: '+sys.argv[0]+' neighbor_index(nn is 1)\n'
        sys.exit(0)

    ### Use smoothed data to locate extrema ###
        
    # Read in the data
    R, R_n = [], []
    for line in lines:
        row = line.split()
        R.append(float(row[0]))
        R_n.append(float(row[neighbor]))

    # Calculate how many peaks the neighbor distribution has
    maxima, minima = [], []
    for i in range(2,len(R_n)-2):
        if R_n[i-2] < R_n[i-1] < R_n[i] > R_n[i+1] > R_n[i+2]:
            maxima.append(i)
        if R_n[i-2] > R_n[i-1] > R_n[i] < R_n[i+1] < R_n[i+2]:
            minima.append(i)

    for i in range(len(maxima)):
        print 'maximum at:',R[maxima[i]]
    for i in range(len(minima)):
        print 'minimum at:',R[minima[i]]

    if len(minima) > 1:
        print '\nMore than one local mimumum found... try further smoothing, exiting...\n'
        sys.exit(0)

    # Open actual un-averaged data
    f = open(fin,'r')
    f.readline()
    lines = f.readlines()
    f.close()

    # Read in the data
    R, R_n = [], []
    for line in lines:
        row = line.split()
        R.append(float(row[0]))
        R_n.append(float(row[neighbor]))

    # Perform integration
    if len(minima) == 1:
        minimum = minima[0]
        X = R[:minimum]
        Y = R_n[:minimum]

        min_integral = integrate.simps(Y,X)
#        min_integral = integrate.trapz(Y,X)
    else:
        min_integral = False

    all_integral = integrate.simps(R_n,R)
#    all_integral = integrate.trapz(R_n,R)

    print '\nIntegral of entire neighbor distribution =',all_integral
    if min_integral:
        print 'Integral up to minimum =',min_integral
        print 'ratio =',min_integral/all_integral
    print


if __name__ == '__main__':
    main()
