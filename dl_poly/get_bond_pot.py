#!/usr/bin/env python

# Program to do get_traj.py, get_pvt.py and create a bond_pot file
# CAREFUL: currently only works for DL_POLY simulations with only 2 atoms!

import os, sys, commands

os.system('get_traj.py')
os.system('get_pvt.py')

traj = open('traj.xyz','r')
traj_lines = traj.readlines()
traj.close()

g = open('CONTROL','r')
CONTROL_lines = g.readlines()
for line in CONTROL_lines:
    if len(line.split()) == 4:
        var, start, step, trjkey = line.split()
    if var == 'traj':
        trjkey = int(trjkey)
        break

bonds = []
bond_initial = int(traj_lines[1])
if trjkey == 0:
    for i in range(len(traj_lines)):
        if i%4 == 3:
            bonds.append(str(float(traj_lines[i].split()[1])*2.))
        if i == len(traj_lines) - 3:
            bond_final = int(traj_lines[i])

elif trjkey == 2:
    for j in range(5): traj_lines.pop(0)
    for i in range(len(traj_lines)):
        if i%8 == 0:
            bonds.append(str(float(traj_lines[i].split()[1])*2.))
        if i == len(traj_lines) - 7:
            bond_final = int(traj_lines[i])
                                        
pvt = open('pvt.dat','r')
pvt.readline() # Remove
pvt.readline() # Headers
pvt_lines = pvt.readlines()
pvt.close()

tstep = []
P = []
vol = []
T = []
etot = []
pot = []
for line in pvt_lines:
    row = line.split()
    tstep.append(row[0])
    P.append(row[1])
    vol.append(row[2])
    T.append(row[3])
    etot.append(row[4])
    pot.append(row[5])

pot_initial = int(tstep[0])
pot_final = int(tstep[-1])
t = int(tstep[1]) - int(tstep[0])

if bond_initial == pot_initial:
    pass
elif bond_initial < pot_initial:
    while bond_initial != pot_initial:
        bonds.pop(0)
        bond_initial += t
elif bond_initial > pot_initial:
    while bond_initial != pot_initial:
        pot.pop(0)
        pot_initial += t

if bond_final == pot_final:
    pass
elif bond_final < pot_final:
    while bond_final != pot_final:
        pot.pop(-1)
        pot_final -= t
elif bond_final > pot_final:
    while bond_final != pot_final:
        bonds.pop(-1)
        bond_final -= t

# Uncomment for a good check of time synchronized data sets
#print bond_initial, pot_initial, bond_final, pot_final, t

# Now write the bond_pot file
out = open('bond_pot.dat','w')
out.write('# Potential (check CONFIG for units) as a function of bond length (angstroms)\n')

#if bond_initial == pot_initial and bond_final == pot_final and len(bonds) == len(pot):
for i in range(len(bonds)):
    out.write(bonds[i]+'\t'+pot[i]+'\n')

out.close()

