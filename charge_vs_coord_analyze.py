#!/usr/bin/env python
"""
Reads in CC.dat file (concatenation of many cc.dat files) and produces
a histogram of coordination vs. charge
"""
import os, sys, commands, glob
import pylab, Numeric

global Nbins, Natom, Nvalence
Nbins, Natom, Nvalence = 200, 64, 5

def main():

    try:
        f = open('CC.dat','r')
        lines = f.readlines()
        f.close()
        Nconfig = int(commands.getoutput('wc -l CC.dat').split()[0])/float(Natom)
    except:
        print '\nMake sure CC.dat file is present.\n'
        sys.exit(0)

    coord, charge = [], []
    for line in lines:
        row = line.split()
        coord.append(int(row[0]))
        charge.append(float(row[1]))

    coord1_charge = [charge[i] for i in range(len(charge)) if coord[i] == 1]
    coord2_charge = [charge[i] for i in range(len(charge)) if coord[i] == 2]
    coord3_charge = [charge[i] for i in range(len(charge)) if coord[i] == 3]

    hist1 = Numeric.zeros(Nbins,typecode=Numeric.Float)
    hist2 = Numeric.zeros(Nbins,typecode=Numeric.Float)
    hist3 = Numeric.zeros(Nbins,typecode=Numeric.Float)
    HIST  = Numeric.zeros(Nbins,typecode=Numeric.Float)

    X = [((max(charge) - min(charge))/Nbins)*(i+1)+min(charge) for i in range(Nbins)]

    for i in range(len(coord1_charge)):
        bin = int( (coord1_charge[i] - min(charge))*Nbins/(max(charge)-min(charge)) )
        hist1[bin-1] += 1
        HIST[bin-1]  += 1
    for i in range(len(coord2_charge)):
        bin = int( (coord2_charge[i] - min(charge))*Nbins/(max(charge)-min(charge)) )
        hist2[bin-1] += 1
        HIST[bin-1]  += 1
    for i in range(len(coord3_charge)):
        bin = int( (coord3_charge[i] - min(charge))*Nbins/(max(charge)-min(charge)) )
        hist3[bin-1] += 1
        HIST[bin-1]  += 1

    norm = Natom*Nvalence*Nconfig

    out = open('coord_charge.dat','w')
    for i in range(len(X)):
        out.write(str(X[i])+'  '+str(HIST[i]/norm)+'  '+str(hist1[i]/norm)+'  '+str(hist2[i]/norm)+'  '+str(hist3[i]/norm)+'\n')

    out.close()

if __name__ == '__main__':
    main()
