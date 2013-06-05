#!/usr/bin/env python

# Script to average the results of several similar diffusion datasets

import os, sys, Numeric

def main():

    """Read in the files and average the results"""

    try:
        nruns = int(sys.argv[1])
    except:
        print '\nPlease given the number of runs to be averaged\n'
        sys.exit(0)

    runs = [i+1 for i in range(nruns)]

    time = []
    datasets = []
    for run in runs:
        datasets.append([])
        f = open('run'+str(run)+'/diffusion.dat','r')
        lines = f.readlines()
        f.close()
        for line in lines:
            time.append( line.split()[0] )
            datasets[run-1].append( float(line.split()[1]) )
        datasets[run-1] = Numeric.array(datasets[run-1])

    avgdiff = Numeric.zeros(len(datasets[0]))
    for dataset in datasets:
        avgdiff = avgdiff + dataset
    avgdiff = avgdiff/len(datasets)

    out = open('avg_diffusion.dat','w')
    for i in range(len(avgdiff)):
        out.write(time[i]+'      '+str(avgdiff[i])+'\n')
    out.close()

if __name__ == '__main__':
    main()
