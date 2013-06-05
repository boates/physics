#!/usr/bin/env python
"""
correlators.py
Author: Brian Boates

Adapted from Isaac's original code.
"""
import sys, os, glob, commands, numpy

def main():

    # Retrieve user input
    try:
        fin    = sys.argv[1]
        nneigh = int(sys.argv[2])
        window = int(sys.argv[3])
        snapshot_stride  = 1  # should be a user input
    except:
        print '\n usage: ' + sys.argv[0] + ' c*n_file nneighbours window_size \n'
        sys.exit(0)

    # Determine size of file
    nlines = int(commands.getoutput('wc -l '+fin).split()[0]) - 10

    # Read in the file header information
    trajec = open(fin,'r')
    trajec.readline()
    alat  = float( trajec.readline().split()[-1] )
    blat  = float( trajec.readline().split()[-1] )
    clat  = float( trajec.readline().split()[-1] )
    natom = int( trajec.readline().split()[-1] )
    nnmax = int( trajec.readline().split()[-1] )
    trajec.readline()
    trajec.readline()
    trajec.readline()
    trajec.readline()

    # Determine the number of timesteps
    nsteps = nlines / natom

    # Read in the neighbor information
    nn = numpy.zeros((natom,nsteps,nneigh), dtype=numpy.int)
    for i in range(nsteps):
        for j in range(natom):
            nn_tmp = trajec.readline().split()[4:4+nneigh]
            nn[j][i] = numpy.sort(nn_tmp)
    trajec.close()

    # Initialize the histogram
    hist = numpy.zeros(window, dtype=float)

    # Perform the analysis
    for atom_history in nn:
        for start_time in range(0, nsteps - window):
            bool_array = numpy.equal( atom_history[start_time+1:start_time+window+1], atom_history[start_time])
            int_array  = numpy.array( bool_array, numpy.int)
            correlator = numpy.add.reduce(int_array, 1) / nneigh
            hist += correlator

    # Write the histogram to output file
    out = open('correlator_' + str(nneigh) + 'nn.dat','w')
    time = 0
    for item in hist:
        out.write(str(time)+' '+str(item/hist[0])+'\n')
        time += 1
    out.close()


if __name__ == '__main__':
    main()
