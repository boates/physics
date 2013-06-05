#!/usr/bin/env python

# Script to analyze bond frequencies using PSD.py

import sys
import PSD

def main():

    # Read in the users data file as two columns
    # timestep, dr (ang)
    try:
        fin = sys.argv[1]
	step = float(sys.argv[2])
	f = open(fin,'r')
    except IndexError:
        print '\nusage '+sys.argv[0]+' f_bond.dat timestep_in_A.U.\n'
	sys.exit(0)

    lines = f.readlines()
    f.close()
    x, y = [], []
    for line in lines:
        row = line.split()
        x.append(float(row[0]))
	y.append(float(row[1]))

    # Calculate the PSD and save a figure
    psd_obj = PSD.psd(x,y)
    f, psd = psd_obj.getPSD()

    # Write PSD data to file
    out = open('PSD.dat','w')
    out.write('#f        psd\n')
    
    # Neglect the zeroth value
    for i in range(len(f)/2-1):
        out.write(str(f[i+1])+'    '+str(psd[i+1])+'\n')

    # Domiqnant frequency
    psd = list(psd)[1:len(psd)/2]
    f = list(f)[1:len(f)/2]
    max_f = f[ psd.index(max(psd)) ]

    # Corresponding period
    T = 1.0 / max_f

    # Convert period from number of timesteps to femtoseconds
    T = T * (step * 2.41880e-02)

    out.close()

    print 'The period is '+str(T)+' femtoseconds'

if __name__ == '__main__':
    main()
