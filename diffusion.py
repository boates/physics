#!/usr/bin/env python

# Script to create a diffusion file from an MSD output file
# MSD file should be named msd.dat

import os, sys, commands

def main():

    """Extract the msd data and output diffusion"""

    try:
        f = open(sys.argv[1],'r')
        out = open('D.dat','w')
    except:
        print '\nusage: '+sys.argv[0]+' msd.dat\n'
        sys.exit(0)

    # Time in ps is given in the 1st column, and msd in the 2nd
    lines = f.readlines()
    time,msd=[],[]
    t = '      '
    for line in lines:
        time.append( float(line.split()[0]) )
        msd.append( float(line.split()[1]) )
    # Einstein's diffusion relation: 2*t*D = 1/3 * msd
    for i in range(len(msd)-1):
        diffusion = abs( (msd[i] / time[i])*(1/6.0) )*0.0001 # Convert to cm^2/s
        out.write( str(time[i])+t+str(diffusion)+'\n' )        

    out.close()

if __name__ == '__main__':
    main()
