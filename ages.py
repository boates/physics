#!/usr/bin/env python
"""
ages.py
Author: Brian Boates

Reads a .cbn file
Creates a histogram of lifetimes
"""
import os, sys, commands, glob
import numpy

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        tstep = float(sys.argv[2])   # in fs
        tbin = int(sys.argv[3])      # in fs
    except:
        print '\n usage: '+sys.argv[0]+'  TRAJEC.cbn  tstep_in_fs  tbin_in_fs(i.e. 5)\n'
        sys.exit(0)

    # Read the cbn header, get natom
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    natom = int(f.readline().split()[-1])
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()

    # Get nsteps, check for consistency with natom
    nlines = int(commands.getoutput("wc -l "+sys.argv[1]+" | awk '{print $1}'"))
    if (nlines - 10) % natom != 0:
        print '\n (nlines-10)/natom is not an integer, check .cbn file --- exiting...\n'
        sys.exit(0)
    else:
        nsteps = (nlines - 10) / natom

    # Read in the bonding information
    NN = [[] for k in range(natom)]
    for i in range(nsteps):
        for j in range(natom):
            NN[j].append(f.readline().split()[-1])
    f.close()

    # Check to see if system remains purely molecular indefinitely
    pure = []
    for i in range(natom):
        # If an atoms bonding information is constant for the entire simulation
        if NN[i].count(NN[i][0]) == len(NN[i]):
            pure.append(True)
        else:
            pure.append(False)
    if pure == [True for i in range(natom)]:
        print '\n System remains purely molecular for the entire simulation --- exiting...'
        print ' Simulation length =',nsteps*tstep/1000.0,'ps\n'
        sys.exit(0)

    immortal = pure.count(True)

    # Track the "births" and "deaths" of molecules
    births, deaths, survivors = [], [], []
    for i in range(natom):
        for j in range(nsteps):
            # Special case: initial timestep
            if j == 0:
                # If atom is bonded at t=0: record a birth
                if NN[i][j] != -1:
                    births.append(j)
            # Special case: final timestep
            elif j == (nsteps-1):
                # If there remains a living molecule: exclude it
                # and count it as a "survivor" as it still lives
                # at the end of the simulation
                if NN[i][j] != -1:
                    survivors.append((nsteps - births.pop(-1))*tstep)
            else:
                # If atom is no longer with previous partner
                if NN[i][j] != NN[i][j-1]:
                    # and is now not bonded: record a death
                    if NN[i][j] == -1:
                        deaths.append(j)
                    # and is with a new partner: record a birth
                    elif NN[i][j] != -1:
                        births.append(j)
                        # If atom wasn't dead in between lives: record a death (for its 1st molecule)
                        if NN[i][j-1] != -1:
                            deaths.append(j)

    # Check that number of births = number of deaths
    if len(deaths) != len(births):
        print '\n Number of births = '+str(len(births))+' != Number of deaths = '+str(len(deaths))
        print ' Exiting...\n'
        sys.exit(0)

    # Determine the ages of the detected molecules in fs
    ages = []
    for i in range(len(births)):
        age = ( deaths[i] - births[i] )*tstep
        ages.append(age)

#    # For the immortal molecules, count as age=nsteps*tstep
#    for i in range(immortal):
#        ages.append(nsteps*tstep)

    # Create a histogram of molecular ages
    nbins = int(max(survivors) / tbin)  # Number of bins required for bin width ~ tbin
    hist = numpy.histogram(numpy.array(ages), bins=nbins, normed=False, range=(0.0,max(ages)))
    norm = numpy.histogram(numpy.array(ages), bins=nbins, normed=True, range=(0.0,max(ages)))

    # Write the histogram to file
    out = open('ages.dat','w')
    out.write('# t(fs), normalized_probability ### tracked '+str(len(deaths))+' molecules\n')
    for i in range(len(hist[0])):
        out.write(str(hist[1][i])+' '+str(norm[0][i])+' '+str(hist[0][i]/2.0)+'\n')
    out.close()

    # Create a histogram for the surviving molecules
    nbins = int(max(survivors) / tbin)  # Number of bins required for bin width ~ tbin
    hist = numpy.histogram(numpy.array(survivors), bins=nbins, normed=False, range=(0.0,max(survivors)))
    norm = numpy.histogram(numpy.array(survivors), bins=nbins, normed=True, range=(0.0,max(survivors)))

    # Write the histogram to file
    out = open('survivors.dat','w')
    out.write('# t(fs), normalized_probability ### '+str(len(survivors)/2.0)+' surviving molecules, '+str(immortal/2.0)+' immortal\n')
    for i in range(len(hist[0])):
        out.write(str(hist[1][i])+' '+str(norm[0][i])+' '+str(hist[0][i]/2.0)+'\n')
    out.close()


if __name__ == '__main__':
    main()
