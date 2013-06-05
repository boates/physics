#!/usr/bin/env python
"""
correlators.py
Author: Brian Boates

Adapted from Isaac's original code.
"""
import sys, os, glob, commands, numpy

def main():

    # Retrieve user input
    

    if len(sys.argv) != 2:
        print 'usage: ' + sys.argv[0] + ' c*n_file nneighbours window_size'

    input_filename  = sys.argv[1]
    number_of_neighbours  = int(sys.argv[2])
    window_size  = int(sys.argv[3])

    snapshot_stride  = 1 # should be a user input

    command_line_counter = commands.getoutput('wc -l ' + input_filename).split()

    if len(command_line_counter) != 2:
        print 'Error determining file size'
    else:
        number_of_lines = int(command_line_counter[0]) - 10 # there are 10 info lines

    size_of_box_x = float(commands.getoutput('grep "# a = " ' + input_filename).split()[3])
    size_of_box_y = float(commands.getoutput('grep "# b = " ' + input_filename).split()[3])
    size_of_box_z = float(commands.getoutput('grep "# c = " ' + input_filename).split()[3])

    cbnFile = open(input_filename, 'r')

    for i in range(10):
        cbnFile.readline()  # skip 10 info lines

    number_of_particles = int(commands.getoutput('grep natom ' + input_filename).split()[3])
    number_of_snapshots = number_of_lines/(number_of_particles)         # the 2 is due to xyz format

    BONDING_array  = numpy.zeros((number_of_particles, number_of_snapshots, number_of_neighbours), dtype=numpy.int)

    for s in range(number_of_snapshots):
        for p in range(number_of_particles):
             string = cbnFile.readline()
             neighbour_list = []
             for n in range(number_of_neighbours):
                 neighbour_list.append(int(string.split()[4 + n]))
             neighbour_list = numpy.sort(neighbour_list)
             BONDING_array[p][s] = neighbour_list
    HISTOGRAM = numpy.zeros(window_size ,dtype=float)

    for particle_history in BONDING_array:                                           # loops over particles
        for start_time in numpy.arange(0, number_of_snapshots - window_size):        # loops over start time
                correlator = numpy.add.reduce(numpy.array(numpy.equal(particle_history[start_time + 1:start_time + window_size + 1], particle_history[start_time]), numpy.int), 1)/number_of_neighbours
                HISTOGRAM += correlator

    histFile = open('correlator_' + str(number_of_neighbours) + 'nn.dat', 'w')

    time = 1

    for element in HISTOGRAM:
        histFile.write(str(time) + ' ' + str(element/HISTOGRAM[0]) + '\n')
        time += 1

    histFile.close()


if __name__ == '__main__':
    main()
