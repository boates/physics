#!/usr/bin/env python

# Format the HISTORY data file to obtain coords of atoms in
# desired format (natom x y z)
# NOTE: this only works for history files generated with keytrj=0
# i.e. only the coords in the file
# Augmentation of this code would be quite simple to include keytrj>0 cases

import os, sys, commands

# Open & prepare the HISTORY file for processing
try:
    f = open('HISTORY','r')
except IOError:
    print '\nHISTORY file does not exist--stopping\n'
    sys.exit(0)
header1 = f.readline()
header2 = f.readline()

# Get necessary info from the CONTROL file
g = open('CONTROL','r')
lines = g.readlines()
for line in lines:
    if len(line.split()) == 2:
        var, value = line.split()
        if var == 'steps':
            steps = int(value)
        elif var == 'equilibration':
            equilibration = int(value)
        elif var == 'print' and value != 'rdf':
            interval = int(value)

try:
    nsteps = int((steps-equilibration)/interval)+1
except NameError:
    print '\n--Could not retrieve all required information from CONTROL file--'
    print 'Need: steps, equilibration, & print\n'
    sys.exit(0)

h = open('traj.xyz','w')
k = open('traj_force.dat','w')

for i in range(nsteps):
    info = f.readline().split()
    timestep, natom, keytrj, pbkey, tstep = info[1:]
    a1 = f.readline().split()
    a2 = f.readline().split()
    a3 = f.readline().split()
    h.write(natom+'\n'+timestep+'\n')

    for j in range(int(natom)):
        atom,iatom,m,q = f.readline().split()
        if int(keytrj) >= 0:
            coords = f.readline().split()
            h.write(atom+'\t'+str(float(coords[0]))+'    '+str(float(coords[1]))+'    '+str(float(coords[2]))+'\n')
            k.write(coords[0]+'\t'+coords[1]+'\t'+coords[2]+'\t')
        if int(keytrj) >= 1:
            vels = f.readline().split()
#            h.write(atom+'\t'+vels[0]+'\t'+vels[1]+'\t'+vels[2]+'\n')
        if int(keytrj) >= 2:
            forces = f.readline().split()
#            h.write(atom+'\t'+forces[0]+'\t'+forces[1]+'\t'+forces[2]+'\n')
            k.write(forces[0]+'\t'+forces[1]+'\t'+forces[2]+'\n')

h.close()
k.close()
