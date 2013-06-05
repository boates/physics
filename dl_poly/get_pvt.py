#!/usr/bin/env python

# Try retrieving P,V,T, etc. from the STATIS file, may be easier than from OUTPUT...

import os, sys, commands

def readlines(FILE,n):
    '''Read n lines from FILE'''
    for i in range(n):
        FILE.readline()
        
try:
    s = open('STATIS','r')
    header1 = s.readline()
    header2 = s.readline()
    
    c = open('CONTROL','r')
    lines = c.readlines()
    for line in lines:
        if len(line.split()) == 2:
            var, value = line.split()
            if var == 'steps':
                steps = int(value)
            elif var == 'stats':
                stats = int(value)
    c.close()
except:
    print 'Could not open STATIS and CONTROL files successfully--stopping'
    sys.exit(0)

# Total energy is row 1 value 1
# Temp is row 1, value 2
# Pres is row 6, value 2
# Vol is row 4, value 4

nblocks = int(steps)/int(stats)
out = open('pvt.dat','w')
out.write('# --Data extracted from STATIS file--\n')
out.write('#tstep\tpres (GPa)\tvol (ang^3)\ttemp (K)\tetot (eV)\t\tpot (eV)\n')
for i in range(nblocks):
    tstep, t, elements = s.readline().split()
    row1 = s.readline().split()
    Etot = str( float(row1[0]) * 1.036426865E-4 ) # convert unit to eV
    T = row1[1]
    s.readline()
    s.readline()
    V = s.readline().split()[3]
    s.readline()
    P = str( float(s.readline().split()[1]) * 0.016605402 ) # convert atm unit to GPa
    
    # Every line has 5 values, each line read is 5 elements gone
    leftover = int(elements) - 5*6
    if leftover % 5 == 0:
        extra_lines = leftover/5
    else:
        extra_lines = leftover/5 + 1
    readlines(s,extra_lines)

    # Calculate Etot - 3*k_b*T
    k_b = 8.617343E-5 # Boltzmann's const in eV/K
    pot = str( float(Etot) - 3*k_b*float(T) )

    out.write(tstep+'\t'+P+' \t'+V+'\t'+T+'\t'+Etot+'\t'+pot+'\n')

s.close()
out.close()
