#!/usr/bin/env python
"""
correlators_bond.py
Author: Brian Boates

To calculate the survival probability of C-C bonds
fin:
tstep 1
1 2
3 4
tstep 2
1 2
5 6
7 8
tstep 3
...

fin is a list of atomic indices for each
timestep, between which exist a C-C bond.

fin is a .bond file from molecules.py
"""
import sys, os, glob, commands, numpy

def main():

    # Retrieve user input
    try:
        fin    = sys.argv[1]
        window = int(sys.argv[2])
    except:
        print '\n usage: ' + sys.argv[0] + ' CC.bond window_size \n'
        print ' CC.bond is a .bond file from molecules.py'
        print ' window_size is how far in time you want to calculate the'
        print ' survival probability, output given in tstep vs. surv. prob.'
        sys.exit(0)

    # Read info from input file
    nsteps = int(commands.getoutput('grep tstep '+fin+' | tail -n-1').split()[-1])
    bonds = open(fin,'r')
    lines = bonds.readlines()
    bonds.close()

    # Read in the C-C pair indices
    all, pairs, t = [], [], -1
    for i in range(len(lines)):
        line = lines[i].split()
        if 'tstep' in line:
            pairs.append([])
            t += 1
        else:
            line = [int(line[0]),int(line[1])]
            line.sort()
            if line not in pairs[t]:
                pairs[t].append(line)
            if line not in all:
                all.append(line)

    # Create a list of lists indicating at each time steps if a given pair is 'present'
    present = []
    for i in range(len(all)):
        present.append([])
        for j in range(len(pairs)):
            if all[i] in pairs[j]:
                present[i].append(1)
            else:
                present[i].append(0)

    # Create windowed lists for an averaged result
    windowed = []
    for p in present:
        for start_time in range(0, len(p)-window):
            slice = p[start_time+1:start_time+window+1]
            if slice[0] == 1:  # Only include if pair is present at first time
                windowed.append(slice)

    # Create normlaized histogram
    hist = numpy.add.reduce( numpy.array(windowed) ) / float(len(windowed))

    # Write to output file
    out = open('correlator_bond.dat','w')
    for i in range(len(hist)):
        out.write(str(i)+' '+str(hist[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
