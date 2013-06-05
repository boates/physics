#!/usr/bin/env python
"""
Revised Isaac code, made to work more generally (i.e. for CO2)
"""
# This is the second generation of this code. In the original version, all 
# angles were calculed, up to some max. This means that angles between the 
# 1st and 2nd shell were always calculated, which isn't necessarily what we i
# want to recover the behaviour of the old code, see the following:
#
# old:		nn_angle_distributions.py nbins 8 TRAJEC.cnn
# current:	nn_angle_distributions.py nbins 1 8 TRAJEC.cnn

import os, sys, commands, string, glob
import numpy

def main():

    # Retrieve user input
    try:
        number_of_bins  = int(sys.argv[1])
        inner_neighbour = int(sys.argv[2]) 
        outer_neighbour = int(sys.argv[3])
        input_filename  = sys.argv[4]
    except:
        print '\nusage: '+sys.argv[0]+' number_of_bins inner_neighbour outer_neighbour FILE.cnn\n'
        sys.exit(0)

    if inner_neighbour > outer_neighbour or inner_neighbour == 0:
        print '\n inner_neighbour must be an integer greater than 0, but less than outer_neighbour... exiting \n'
        sys.exit(0)

    # Should be a user input
    snapshot_stride = 1

    # Check size of file
    command_line_counter = commands.getoutput('wc -l ' + input_filename).split()
    if len(command_line_counter) != 2:
        print '\nError determining file size... exiting\n'
        sys.exit(0)
    else:
        number_of_lines = int(command_line_counter[0]) - 10 # there are 10 header lines
    
    size_of_box_x = float(commands.getoutput('grep "# a = " ' + input_filename).split()[3])
    size_of_box_y = float(commands.getoutput('grep "# b = " ' + input_filename).split()[3])
    size_of_box_z = float(commands.getoutput('grep "# c = " ' + input_filename).split()[3])
    
    cnnFile = open(input_filename, 'r')

    # Skip 10 info lines
    for i in range(10):
        cnnFile.readline()
    
    number_of_particles = int(commands.getoutput('grep number_of_particles ' + input_filename).split()[3])
    number_of_snapshots = number_of_lines/(number_of_particles)         
    number_of_neighbours = int(number_of_particles - 1)
    
    ###### histgram stuff
    
    min_distance = 0.0
    max_distance = numpy.sqrt((size_of_box_x/2.0)**2 + (size_of_box_y/2.0)**2 + (size_of_box_z/2.0)**2)
    
    min_theta = 0.0
    max_theta = numpy.pi
    
    histogram = numpy.zeros(number_of_bins, dtype=numpy.int)
    theta_step = (max_theta - min_theta)/float(number_of_bins)
    
    s = 0
    
    while s < number_of_snapshots:
    
        SNAPSHOT_array = numpy.zeros((number_of_particles, 3), dtype=numpy.float)
        neighbour_list = numpy.zeros((number_of_particles, 1 + number_of_neighbours), dtype=numpy.int)
    
        for p in range(number_of_particles):
    
            line = cnnFile.readline()
    
            for i in range(3):
                SNAPSHOT_array[p][i] = float(line.split()[i + 1])
    
            for c in numpy.arange(inner_neighbour, outer_neighbour + 1):
                neighbour_list[p][c] = int(line.split()[c - 1 + 4])
    
        # At this point, SNAPSHOT_array an neighbourlist are ready to go
        for p in range(number_of_particles):
    
            px = SNAPSHOT_array[p][0]
            py = SNAPSHOT_array[p][1]
            pz = SNAPSHOT_array[p][2]
    
            for o in range(number_of_particles):
    
                x = SNAPSHOT_array[o][0] - px
                y = SNAPSHOT_array[o][1] - py
                z = SNAPSHOT_array[o][2] - pz
    
                SNAPSHOT_array[o][0] = x - (int(x/size_of_box_x+number_of_snapshots+0.5)-number_of_snapshots)*size_of_box_x
                SNAPSHOT_array[o][1] = y - (int(y/size_of_box_y+number_of_snapshots+0.5)-number_of_snapshots)*size_of_box_y
                SNAPSHOT_array[o][2] = z - (int(z/size_of_box_z+number_of_snapshots+0.5)-number_of_snapshots)*size_of_box_z
    
            for o in numpy.arange(inner_neighbour, outer_neighbour + 1):
                
                o_vector = numpy.array([SNAPSHOT_array[neighbour_list[p][o]][0], SNAPSHOT_array[neighbour_list[p][o]][1], SNAPSHOT_array[neighbour_list[p][o]][2]])
                
                for oo in numpy.arange(o + 1, outer_neighbour + 1):
    
                    oo_vector = numpy.array([SNAPSHOT_array[neighbour_list[p][oo]][0], SNAPSHOT_array[neighbour_list[p][oo]][1], SNAPSHOT_array[neighbour_list[p][oo]][2]])
    
                    numerator = numpy.dot(o_vector,oo_vector)
                    denomenator = numpy.sqrt(numpy.dot(o_vector, o_vector))*numpy.sqrt(numpy.dot(oo_vector, oo_vector))
    
                    theta = numpy.arccos(numerator/denomenator)
                    bin = int((theta - min_theta)/theta_step)
    #                print o, oo
                    histogram[bin] += 1
    
        outputFile_nn_angle = open ('nn_angle.' + string.zfill(inner_neighbour,3) + '_' + string.zfill(outer_neighbour,3) + '.hist', 'w')
        outputFile_nn_angle.write('# bin (deg), p(theta), nsteps_used = ' + str(number_of_snapshots) + '\n') 
    
        for theta_bin in range(number_of_bins):
    
    #       volume_element = numpy.cos(theta_bin*theta_step) - numpy.cos((theta_bin+1)*theta_step)
            volume_element = 1.0
            outputFile_nn_angle.write(repr( (theta_bin*theta_step + theta_step/2.0 + min_theta)*180.0/numpy.pi ) + '  ')
            outputFile_nn_angle.write(str(float(outer_neighbour + 1 - inner_neighbour)*(histogram[theta_bin]/float(histogram.sum()))/volume_element) + ' ')
            outputFile_nn_angle.write('\n')
    
        outputFile_nn_angle.close()
    
    #    print s    
    
        s += snapshot_stride
    

if __name__ == '__main__':
    main()
