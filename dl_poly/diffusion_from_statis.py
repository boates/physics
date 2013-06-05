#!/usr/bin/env python

# Script to create a diffusion file from a statis.dat file
# created by STATIS_extract.f / STATIS_extract_formatter.py

import os, sys, commands

def main():

    """Extract the msd data and output diffusion"""

    try:
        f = open('statis.dat','r')
        out = open('diffusion.dat','w')
    except:
        print '\nPlease make sure the statis.dat file is present\n'
        sys.exit(0)

    # Time in ps is given in the second column, and msd is in the eighth
    lines = f.readlines()
    t = '      '
    for line in lines:
        time = float( line.split()[1] )
        msd  = float( line.split()[7] )
        # Einstein's diffusion relation: 2*t*D = 1/3 * msd
        diffusion = ( ( msd / 6.0 ) / time )
        out.write( str(time)+t+str(diffusion)+'\n' )        

    out.close()

if __name__ == '__main__':
    main()
