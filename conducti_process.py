#!/usr/bin/env python
"""
Script to process conducti calculations
"""
import os, sys, commands, glob, math, numpy

def main():

    # Find all conducti output files
    sig_files = glob.glob('*/CONOUT_sig')
    abs_files = glob.glob('*/abs.dat')

    # get the frequencies, energies, and corresponding wavelengths
    w = commands.getoutput("awk '{print $1}' "+sig_files[0]).split()
    E = commands.getoutput("awk '{print $2}' "+sig_files[0]).split()
    w.pop(0)
    E.pop(0)
    L = []
    for i in range(len(E)):
        E[i] = float(E[i])
        w[i] = float(w[i])
        L.append( (6.5822*29.979*2*math.pi) / E[i] )

    w = numpy.array(w)
    E = numpy.array(E)
    L = numpy.array(L)
    
    # Extract and average the conductivities
    cond = []
    for i in range(len(sig_files)):
        f = sig_files[i]
        c = commands.getoutput("awk '{print $4}' "+f).split()
        c.pop(0)
        cond.append(c)

    # Calculate the averaged conductivity
    cond_avg = numpy.zeros(len(cond[0]))
    for i in range(len(cond[0])):
        for j in range(len(cond)):
            cond_avg[i] += float(cond[j][i]) / len(cond)

    # Extract and average the reflectivities
    refl = []
    for i in range(len(abs_files)):
        f = abs_files[i]
        r = commands.getoutput("awk '{print $4}' "+f).split()
        refl.append(r)

    # Calculate the averaged reflectivity
    refl_avg = numpy.zeros(len(refl[0]))
    for i in range(len(refl[0])):
        for j in range(len(refl)):
            refl_avg[i] += float(refl[j][i]) / len(refl)

    # Write to output file
    out = open('CONDUCTI.dat','w')
    out.write('# w(ua), E(eV), lambda(nm), cond((ohm*cm)^-1), refl\n')
    for i in range(len(cond_avg)):
        out.write(str(w[i])+' '+str(E[i])+' '+str(L[i])+' '+str(cond_avg[i])+' '+str(refl_avg[i])+'\n')


if __name__ == '__main__':
    main()
