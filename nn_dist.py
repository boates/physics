#!/usr/bin/env python
"""
Revised Isaac code: calculate nn-histograms and create .cnn file
"""
import os, sys, commands, glob
import numpy

def pbc_round(input_value):
    i = int(input_value)
    if (abs(input_value-i) >= 0.5):
        if (input_value > 0): i+= 1
        if (input_value < 0): i-= 1
    return i

def main():

    # Retrieve user input
    try:
        number_of_bins = int(sys.argv[1])    
        size_of_box_x = float(sys.argv[2])
        size_of_box_y = float(sys.argv[3])
        size_of_box_z = float(sys.argv[4])
        fname = sys.argv[5]
        try:
            label = sys.argv[6]
        except: label = ''
    except:
        print '\n usage: '+sys.argv[0]+' number_of_bins alat_x alat_y alat_z(bohr) FILE.xyz label(optional)\n'
        sys.exit(0)

    # Should be a user input
    snapshot_stride = 1 
    if snapshot_stride < 1:
        snapshot_stride = 1
    
    # Check for size of file
    command_line_counter = commands.getoutput('wc -l ' + fname).split()
    if len(command_line_counter) != 2:
        print 'Error determining file size... exiting'
        sys.exit(0)
    else:
        number_of_lines = int(command_line_counter[0])

    if label != '':
        fout = 'TRAJEC_'+label+'.cnn'
    else:
        fout = 'TRAJEC.cnn'

    # Open needed files
    xyzFile = open(fname, 'r')
    cnnFile = open(fout, 'w')

    # Read info from xyz file
    number_of_particles = int(xyzFile.readline())
    number_of_snapshots = number_of_lines/(number_of_particles + 2)
    number_of_neighbours = int(number_of_particles - 1)

    # Write cnn file header
    cnnFile.write('# comment =  \n')
    cnnFile.write('# a = ' + str(size_of_box_x) + '\n')
    cnnFile.write('# b = ' + str(size_of_box_y) + '\n')
    cnnFile.write('# c = ' + str(size_of_box_z) + '\n')
    cnnFile.write('# number_of_particles = '  + str(number_of_particles) + '\n')
    cnnFile.write('# number_of_neighbours = ' + str(number_of_neighbours) + '\n')
    cnnFile.write('#\n')
    cnnFile.write('#\n')
    cnnFile.write('#\n')
    cnnFile.write('# units = bohr\n')

    #==============================#
    # BEGIN HISTOGRAM CALCULATIONS #
    #==============================#
    
    min_distance = 0.0
    max_distance = max(size_of_box_x, size_of_box_y, size_of_box_z)/2.0

    histogram = numpy.zeros((int(number_of_neighbours), int(number_of_bins)), dtype=numpy.int)
    bin_size = (max_distance - min_distance)/float(number_of_bins)    
    
    SNAPSHOT_array = numpy.zeros((number_of_particles, 3), dtype=numpy.float)

    # Return to beginning of xyzFile
    xyzFile.seek(0)
    
    for s in range(number_of_snapshots):
    
        distance_array, typat = [], []
    
        xyzFile.readline()    # Skip number of particles
        xyzFile.readline()    # Skip MD step or cell parameters

        for p in range(number_of_particles):
            line = xyzFile.readline()
            typat.append(line.split()[0])
            SNAPSHOT_array[p][0] = float(line.split()[1])/0.529177
            SNAPSHOT_array[p][1] = float(line.split()[2])/0.529177
            SNAPSHOT_array[p][2] = float(line.split()[3])/0.529177
    
        # At this point, SNAPSHOT_array is ready to go
    
        for p in range(number_of_particles):
             
             cnnFile.write(typat[p] + ' % .8e' % SNAPSHOT_array[p][0] + ' % .8e' % SNAPSHOT_array[p][1] + ' % .8e' % SNAPSHOT_array[p][2])
    
             neighbour_distances = [] 
    
             for o in range(number_of_particles):
    
                  dx = SNAPSHOT_array[p][0] - SNAPSHOT_array[o][0]
                  dy = SNAPSHOT_array[p][1] - SNAPSHOT_array[o][1]
                  dz = SNAPSHOT_array[p][2] - SNAPSHOT_array[o][2]
    
                  dx -= size_of_box_x*pbc_round(dx/size_of_box_x)
                  dy -= size_of_box_y*pbc_round(dy/size_of_box_y)
                  dz -= size_of_box_z*pbc_round(dz/size_of_box_z)
    
                  distance = (dx**2 + dy**2 + dz**2)**(0.5)
    
                  neighbour_distances.append([distance, o])
    
             list.sort(neighbour_distances)
    
             distance_array.append([])
    
             for i in range(number_of_neighbours):
                  distance_array[-1].append(neighbour_distances[i + 1][0])
                  cnnFile.write(' ' + str.rjust(str(neighbour_distances[i + 1][1]), 3) ) 
    
             cnnFile.write('\n')
    
        # At this point, distance_array is ready to go
    
        for distance_list in distance_array:
            for neighbour_index in range(number_of_neighbours):
                if distance_list[neighbour_index] < max_distance:
                    bin = int((distance_list[neighbour_index] - min_distance)/bin_size)
                    histogram[neighbour_index][bin] += 1

    if label != '':
        f_hist = 'nn_average_'+label+'.hist'
    else:
        f_hist = 'nn_average.hist'
    cnnFile_nn_average = open (f_hist, 'w')
    cnnFile_nn_average.write('# bin (Angst), nn1, nn2, ... \n')
    
    for bin_index in range(number_of_bins):
    
         cnnFile_nn_average.write(repr( (bin_index*bin_size + bin_size/2.0 + min_distance)*0.529177 ) + '  ')
         for neighbour_index in range(number_of_neighbours):
              cnnFile_nn_average.write(str(histogram[neighbour_index][bin_index]/float(number_of_particles*number_of_snapshots)) + ' ')
         cnnFile_nn_average.write('\n')
    
    cnnFile_nn_average.close()


if __name__ == '__main__':
     main()
