#!/usr/bin/env python
"""
Calcaulate histograms of charge associated with C and O atoms
based on bader charge analysis.
"""
import os, sys, commands, glob, Numeric

def main():

    # Retrive user input
    try:
        Nbins = int(sys.argv[1])
        nC    = int(sys.argv[2])
        nO    = int(sys.argv[3])
    except:
        print '\n usage: '+sys.argv[0]+' Nbins(i.e. 100) number_of_carbons number_of_oxygens\n'
        sys.exit(0)

    # Glob for all available ACF.dat files
    ACF = glob.glob('*/ACF.dat')

    # Initiate C & O histograms
    Chist, Ohist = Numeric.zeros(Nbins), Numeric.zeros(Nbins)

    # Loop over all found ACF.dat files and analyze
    Cnorm, Onorm = 0.0, 0.0
    for acf in ACF:
        f = open(acf,'r')
        header = f.readline()
        header = f.readline()
        lines = f.readlines()
        lines.pop(-1)
        lines.pop(-1)
        for line in lines:
            row = line.split()
            i = int(row[0])
            c = float(row[4])
            if i <= nC:
                Cbin = int(c*Nbins/16.0 + 0.5)
                Chist[Cbin] += 1
                Cnorm += 1.0
            elif i > nC:
                Obin = int(c*Nbins/16.0 + 0.5)
                Ohist[Obin] += 1
                Onorm += 1.0

    # Write the histograms to file
    Cout = open('bader_C.hist','w')
    Oout = open('bader_O.hist','w')
    for i in range(Nbins):
        Cout.write( str(i*16.0/Nbins)+' '+str(Chist[i]/Cnorm)+'\n' )
        Oout.write( str(i*16.0/Nbins)+' '+str(Ohist[i]/Onorm)+'\n' )
    Cout.close()
    Oout.close()


if __name__ == '__main__':
    main()
