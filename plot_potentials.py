#!/usr/bin/env python
"""
Make a contourf plot of pair potentials (Poten_NN.dat files)
NOTE: Assumes all potentials are created with same distance spacing in Poten files
"""
import os, sys, commands, glob
import pylab, Numeric

def main():

    potens = glob.glob('1.*/FMTW/FORCE_OUT/Poten_*.dat')
    Nlines, RS = [], []
    for poten in potens:
        Nlines.append(int(commands.getoutput('wc -l '+poten).split()[0]))
        RS.append(poten.split('/')[0])

    # Truncate all potentials to be the same length
    os.system('rm -rf Potens/')
    os.system('mkdir Potens/')
    for poten in potens:
        os.system('head -'+str(min(Nlines))+' '+poten+' > Potens/'+str(RS[potens.index(poten)])+'_'+poten.split('/')[-1])

    # Enter new potential directory and read in and format data
    cwd = os.getcwd()
    os.chdir('Potens/')
    potens = glob.glob('*Poten_*.dat')
    potens.sort()
    RS = [float(poten.split('_')[0]) for poten in potens]
    POTENS = []
    for poten in potens:
        R = []
        POTENS.append([])
        f = open(poten,'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            row = line.split()
            r = float(row[0])
            p = float(row[1])
            if p == 0.0:
                p = 1.0E-6
            R.append(r)
            POTENS[potens.index(poten)].append(p)

    # Convert RS values to floats
    for i in range(len(RS)):
        RS[i] = float(RS[i])

    # Create the plot
    pylab.figure(num=1,figsize=(8,4),facecolor='w',edgecolor='k')
    pylab.contourf(RS,R,Numeric.transpose(POTENS),200,antialiased=False)
    pylab.title('Pair Potential vs. Density',fontsize=14)
    pylab.xlabel('Rs',fontsize=12)
    pylab.ylabel('R',fontsize=12)
    
    pylab.savefig('POTENS.png')
    pylab.clf()

#    pylab.show()


if __name__ == '__main__':
    main()
