#!/usr/bin/env python
import os

# Bonding lifetime analysis
os.system('nice -15 correlators.py TRAJEC.cnn 1  1000 >  cluster1.out')
os.system('nice -15 correlators.py TRAJEC.cnn 2  1000 >  cluster2.out')
os.system('nice -15 correlators.py TRAJEC.cnn 3  1000 >  cluster3.out')
os.system('nice -15 correlators.py TRAJEC.cnn 4  1000 >  cluster4.out')
os.system('nice -15 correlators.py TRAJEC.cnn 5  1000 >  cluster5.out')
os.system('nice -15 correlators.py TRAJEC.cnn 6  1000 >  cluster6.out')
os.system('nice -15 correlators.py TRAJEC.cnn 7  1000 >  cluster7.out')
os.system('nice -15 correlators.py TRAJEC.cnn 8  1000 >  cluster8.out')
os.system('rm -f cluster*.out')

# Convert from tstep to picoseconds
os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_1nn.dat  > c01.dat")
os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_2nn.dat  > c02.dat")
os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_3nn.dat  > c03.dat")
os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_4nn.dat  > c04.dat")
os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_5nn.dat  > c05.dat")
os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_6nn.dat  > c06.dat")
os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_7nn.dat  > c07.dat")
os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_8nn.dat  > c08.dat")
os.system('rm -f correlator_*nn.dat')

# Grab the correlators at given timesteps
os.system('grab_correlators_at_timestep_CO2.sh 100')
