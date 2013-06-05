#!/usr/bin/env python
"""
Find the '1st minimum' in VACF for correlator calculations
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        lines = f.readlines()
        f.close()
        tstep = float(sys.argv[2]) / 1000.0
    except:
        print '\n usage: '+sys.argv[0]+'  VACF.dat  tstep(fs)\n'
        sys.exit(0)

    # Get the VACF data in ps
    t, vacf = [], []
    for line in lines:
        row = line.split()
        t.append( float(row[0]) )
        vacf.append( float(row[2]) )

    # Find the 1st minimum
    vacfMin = 2.0
    i = 0
    while vacf[i] <= vacfMin:
        vacfMin = vacf[i]
        i += 1

    tMin = t[i-1]
    nSteps = tMin / tstep
    
    # Print the result
    out = open('min_vacf.dat','w')
    out.write(str(tMin)+'  '+str(int(round(nSteps)))+'\n')
    out.close()


if __name__ == '__main__':
    main()
