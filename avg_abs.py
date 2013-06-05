#!/usr/bin/env python
"""
Average multiple reflectivity calculations into one dataset
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
        print ' such as, abs01.dat abs02.dat ... abs08.dat, you would hand this code "abs")\n'
        sys.exit(0)

    # Grab list of DOS files to be averaged
#    ABSs = glob.glob(fname+'*.dat')
    ABSs = glob.glob('0*/CONOUT_abs')
    nABS = len(ABSs)
    R = []

    # Get the data from all detected files
    for a in ABSs:
        E = commands.getoutput("awk '{print $1}' "+a).split()[1:]   # in eV
        R.append(commands.getoutput("awk '{print $4}' "+a).split()[1:])   # 0 < R < 1
        
    # Average the data and write to file
    out = open('avg_abs.dat','w')
    out.write('# Average of '+str(nABS)+' abs data files (units: eV vs. 0 < R < 1), stdev\n')
    for i in range(len(E)):
        abs_avg = 0.0
        stdev_list = []
        for j in range(nABS):
            abs_avg += float(R[j][i]) / nABS
            stdev_list.append(float(R[j][i]))
        sd = stdev(stdev_list)
        out.write(str(E[i])+' '+str(abs_avg)+' '+str(sd)+'\n')
    out.close()

if __name__ == '__main__':
    main()
