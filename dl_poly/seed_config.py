#!/usr/bin/env python

# Create a cubic seed inside a lattice.

import os, sys, commands

def main():

    f = open('CONFIG','r')
    out = open('CONFIG_SEED','w')

    out.write(f.readline())
    out.write(f.readline())
    line = f.readline()
    out.write(line)
    ax = float(line.split()[0])
    line = f.readline()
    out.write(line)
    ay = float(line.split()[1])
    line = f.readline()
    out.write(line)
    az = float(line.split()[2])

    lines = f.readlines()
    coords = [line for line in lines if lines.index(line) % 2 == 1]

    coord_seed = []
    coord_remaining = []
    for coord in coords:
        x = float(coord.split()[0])
        y = float(coord.split()[1])
        z = float(coord.split()[2])
        if x < -(ax/2)/2 and y < -(ay/2)/2 and z < -(az/2)/2:
            coord_seed.append(coord)
        else:
            coord_remaining.append(coord)

    for i in range(len(coord_seed)):
        out.write('SS                '+str(i+1)+'\n')
        out.write(coord_seed[i])
    for i in range(len(coord_remaining)):
        out.write('LI                '+str(i+len(coord_seed)+1)+'\n')
        out.write(coord_remaining[i])

    f.close()
    out.close()

if __name__ == '__main__':
    main()
    
