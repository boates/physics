#!/usr/bin/env python
"""
This is a modification of Isaac's /usr/local/analysis/nn_angle_distribution.py
It is inteneded to include an extra criteria, namely to track angles for atoms
of a specific coordination.
"""
import os, sys, commands, glob
import numpy
import string

def main():

    # Retrieve user input
    try:
        fname = sys.argv[1]
        cnnFile  = open(fname,'r')
        for i in range(10): cnnFile.readline()	# skip 10 info lines
        f = open(glob.glob('coordination_rc_*.dat')[0],'r')
        for i in range(10): f.readline() # remove header line
        inner = int(sys.argv[2]) 
        outer = int(sys.argv[3])
        stride = int(sys.argv[4])
        Nbins  = int(sys.argv[5])
    except:
        print '\n usage: '+sys.argv[0]+' TRAJEC.cnn inner_neighbor, outer_neighbor, snapshot_stride(=1), Nbins\n'
        print ' Make sure coordination_rc_*.dat file is present.\n'
        sys.exit(0)

    # Check that inner and outer are set reasonably
    if inner > outer or inner == 0:
       print '\ninner_neighbor must be an integer greater than 0, but less than outer_neighbor, exiting...\n'
       sys.exit(0)

    # Determine length of cnn file
    try:
        Nlines = int(commands.getoutput('wc -l ' + fname).split()[0]) - 10 # 10 info lines
    except:
        print '\nError determining length of cnn file, exiting...\n'
        sys.exit(0)

    # Set lattice constants
    a_x = float(commands.getoutput('grep "# a = " ' + fname).split()[3])
    a_y = float(commands.getoutput('grep "# b = " ' + fname).split()[3])
    a_z = float(commands.getoutput('grep "# c = " ' + fname).split()[3])

    # Set system variables
    natom = int(commands.getoutput('grep number_of_particles ' + fname).split()[3])
    Nsteps = Nlines/(natom)         
    Nneighbors = int(natom - 1)

    # Read in coordination info to an array
    COORD_array = []
    for i in range(Nsteps):
        COORD_array.append([])
        for j in range(natom):
            COORD_array[i].append(int(f.readline().strip()))
    COORD_array = numpy.array(COORD_array)

    # Histgram details
    min_distance = 0.0
    max_distance = numpy.sqrt((a_x/2.0)**2 + (a_y/2.0)**2 + (a_z/2.0)**2)
    min_theta = 0.0
    max_theta = numpy.pi
    histogram = [numpy.zeros(Nbins, dtype=numpy.int) for i in range(3)]
    theta_step = (max_theta - min_theta)/float(Nbins)

    # Begin the analysis using Isaac's already written algorithm
    s = 0
    while s < Nsteps:

        SNAPSHOT_array = numpy.zeros((natom, 3), dtype=numpy.float)
        neighbour_list = numpy.zeros((natom, 1 + Nneighbors), dtype=numpy.int)

        for p in range(natom):

            line = cnnFile.readline()

            for i in range(3):
                SNAPSHOT_array[p][i] = float(line.split()[i + 1])

            for c in numpy.arange(inner, outer + 1):
                neighbour_list[p][c] = int(line.split()[c - 1 + 4])

        # At this point, SNAPSHOT_array an neighbourlist are ready to go
        for p in range(natom):

            px = SNAPSHOT_array[p][0]
            py = SNAPSHOT_array[p][1]
            pz = SNAPSHOT_array[p][2]

            pc = COORD_array[s][p]

            for o in range(natom):

                x = SNAPSHOT_array[o][0] - px
                y = SNAPSHOT_array[o][1] - py
                z = SNAPSHOT_array[o][2] - pz

                SNAPSHOT_array[o][0] = x - (int(x/a_x+Nsteps+0.5)-Nsteps)*a_x
                SNAPSHOT_array[o][1] = y - (int(y/a_y+Nsteps+0.5)-Nsteps)*a_y
                SNAPSHOT_array[o][2] = z - (int(z/a_z+Nsteps+0.5)-Nsteps)*a_z

            for o in numpy.arange(inner, outer + 1):
            
                o_vector = numpy.array([SNAPSHOT_array[neighbour_list[p][o]][0], SNAPSHOT_array[neighbour_list[p][o]][1], SNAPSHOT_array[neighbour_list[p][o]][2]])
            
                for oo in numpy.arange(o + 1, outer + 1):

                    oo_vector = numpy.array([SNAPSHOT_array[neighbour_list[p][oo]][0], SNAPSHOT_array[neighbour_list[p][oo]][1], SNAPSHOT_array[neighbour_list[p][oo]][2]])

                    numerator = numpy.dot(o_vector,oo_vector)
                    denomenator = numpy.sqrt(numpy.dot(o_vector, o_vector))*numpy.sqrt(numpy.dot(oo_vector, oo_vector))

                    theta = numpy.arccos(numerator/denomenator)
                    bin = int((theta - min_theta)/theta_step)
#                    print o, oo
                    histogram[pc-1][bin] += 1

        outputFile_nn_angle = open ('nn_angle.' + string.zfill(inner,3) + '_' + string.zfill(outer,3) + '.hist', 'w')
        outputFile_nn_angle.write('# bin (deg), p(theta): coord=1, coord=2, coord=3, nsteps_used = ' + str(Nsteps) + '\n')

        for theta_bin in range(Nbins):

#            volume_element = numpy.cos(theta_bin*theta_step) - numpy.cos((theta_bin+1)*theta_step)
            volume_element = 1.0
            outputFile_nn_angle.write(repr( (theta_bin*theta_step + theta_step/2.0 + min_theta)*180.0/numpy.pi ) + '  ')
            outputFile_nn_angle.write(str(float(outer + 1 - inner)*(histogram[0][theta_bin]/float(histogram[0].sum()))/volume_element) + ' ')
            outputFile_nn_angle.write(str(float(outer + 1 - inner)*(histogram[1][theta_bin]/float(histogram[1].sum()))/volume_element) + ' ')
            outputFile_nn_angle.write(str(float(outer + 1 - inner)*(histogram[2][theta_bin]/float(histogram[2].sum()))/volume_element) + ' ')
            outputFile_nn_angle.write('\n')

        outputFile_nn_angle.close()

#        print s    

        s += stride


if __name__ == '__main__':
    main()
