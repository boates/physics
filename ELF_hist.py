#!/usr/bin/env python
"""
ELF_hist.py
Author: Brian Boates

Read an ELFCAR file and return a histogram of
values from the entire grid, essentially removing
any spatial component from the analysis for a broad
overview of the ELF across an entire system.
"""
import os, sys, commands, glob, numpy

def main():

    # Retrieve user input
    try:
        elf = open(sys.argv[1],'r')
        nbins = int(sys.argv[2])
    except:
        print '\n usage: '+sys.argv[0]+' ELFCAR nbins\n'
        sys.exit(0)

    # Read header (POSCAR) information
    elf.readline()
    alat = float(elf.readline())
    ax, ay, az = elf.readline().split()
    bx, by, bz = elf.readline().split()
    cx, cy, cz = elf.readline().split()
    ax, ay, az = alat*float(ax), alat*float(ay), alat*float(az)
    bx, by, bz = alat*float(ax), alat*float(ay), alat*float(az)
    cx, cy, cz = alat*float(ax), alat*float(ay), alat*float(az)
    natom = int(elf.readline())
    elf.readline()
    x, y, z = [], [], []
    for i in range(natom):
        row = elf.readline().split()
        x.append(float(row[0]))
        y.append(float(row[0]))
        z.append(float(row[0]))
    elf.readline()

    # Now process the ELF data
    grid = elf.readline().split()

    # Number of elements to read from the ELFCAR file
    nGrid = int(grid[0])*int(grid[1])*int(grid[2])

    # Number of elements VASP prints per line to ELFCAR
    nRow = 10

    # Number of lines to read from ELFCAR
    nLines = nGrid / nRow
    if nGrid % nRow > 0: nLines += 1

    # Read in the ELF
    ELF = []
    for i in range(nLines):
        row = elf.readline().split()
        for r in row:
            ELF.append(float(r))

    # Create the histogram and write to file
    hist = numpy.histogram(numpy.array(ELF),bins=nbins,normed=False,range=(0.0,1.0))

    out = open('hist.ELF','w')
    for i in range(len(hist[0])):
        out.write(str(hist[1][i])+' '+str(hist[0][i]/float(nGrid))+'\n')
    out.close()


if __name__ == '__main__':
    main()
