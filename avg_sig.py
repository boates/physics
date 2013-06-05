#!/usr/bin/env python
"""
Average multiple conductivity calculations into one dataset
Assumes all use same energy axis
Also calculates standard deviations of averages
"""
import os, sys, commands, glob, numpy

def stdev(X):
    """
    Calculate standard deviation of list X.
    """
    mean = float(numpy.average(X))
    sd = 0.0
    for x in X:
        sd += (float(x) - mean)**2
    sd = (sd/len(X))**0.5

    return sd

def main():

    # Retrieve user input
    try:
        fname = sys.argv[1]
    except:
        print '\n usage: '+sys.argv[0]+' fname_beginning\n'
        print ' (i.e. all must be similarly named (presumeably numbered) and end with .dat'
        print ' such as, sig01.dat sig02.dat ... sig08.dat, you would hand this code "sig")\n'
        sys.exit(0)

    # Grab list of DOS files to be averaged
    SIGs = glob.glob(fname+'*.dat')
    nSIG = len(SIGs)
    C = []

    # Get the data from all detected files
    for sig in SIGs:
        E = commands.getoutput("grep -v '\#' "+sig+" | awk '{print $2}'").split()   # in eV
        C.append(commands.getoutput("grep -v '\#' "+sig+" | awk '{print $4}'").split())   # in (ohm*cm)^-1
        
    # Average the data and write to file
    out = open('avg_sig.dat','w')
    out.write('# Average of '+str(nSIG)+' sig data files (units: eV vs. (ohm*cm)^-1), stdev\n')
    for i in range(len(E)):
        sig_avg = 0.0
        stdev_list = []
        for j in range(nSIG):
            sig_avg += float(C[j][i]) / nSIG
            stdev_list.append(float(C[j][i]))
        sd = stdev(stdev_list)
        out.write(str(E[i])+' '+str(sig_avg)+' '+str(sd)+'\n')
    out.close()

if __name__ == '__main__':
    main()
