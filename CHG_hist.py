#!/usr/bin/env python
"""
CHG_hist.py
Author: Brian Boates

Read an CHGCAR file and return a histogram of
values from the entire grid, essentially removing
any spatial component from the analysis for a broad
overview of the charge density across an entire system.
"""
import os, sys, commands, glob, numpy

def main():

    # Retrieve user input
    try:
        chg = open(sys.argv[1],'r')
        nbins = int(sys.argv[2])
    except:
        print '\n usage: '+sys.argv[0]+' CHGCAR nbins\n'
        sys.exit(0)

    # Read header (POSCAR) information
    chg.readline()
    alat = float(chg.readline())
    ax, ay, az = chg.readline().split()
    bx, by, bz = chg.readline().split()
    cx, cy, cz = chg.readline().split()
    ax, ay, az = alat*float(ax), alat*float(ay), alat*float(az)
    bx, by, bz = alat*float(ax), alat*float(ay), alat*float(az)
    cx, cy, cz = alat*float(ax), alat*float(ay), alat*float(az)
    natom = int(chg.readline())
    chg.readline()
    x, y, z = [], [], []
    for i in range(natom):
        row = chg.readline().split()
        x.append(float(row[0]))
        y.append(float(row[0]))
        z.append(float(row[0]))
    chg.readline()

    # Now process the CHG data
    grid = chg.readline().split()

    # Number of elements to read from the CHGCAR file
    nGrid = int(grid[0])*int(grid[1])*int(grid[2])

    # Number of elements VASP prints per line to CHGCAR
    nRow = 10

    # Number of lines to read from CHGCAR
    nLines = nGrid / nRow
    if nGrid % nRow > 0: nLines += 1

    # Read in the CHG
    CHG = []
    for i in range(nLines):
        row = chg.readline().split()
        for r in row:
            CHG.append(float(r))

    # Create the histogram and write to file
    hist = numpy.histogram(numpy.array(CHG),bins=nbins,normed=False,range=(min(CHG),max(CHG)))

    out = open('hist.CHG','w')
    for i in range(len(hist[0])):
        out.write(str(hist[1][i])+' '+str(hist[0][i]/float(nGrid))+'\n')
    out.close()


if __name__ == '__main__':
    main()
