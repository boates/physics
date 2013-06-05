#!/usr/bin/env python
"""
loop_kpt_snapshots.py
Author: Brian Boates

Select snapshots randomly from a wrapped xyz file and convert
to POSCAR's, create number directories and submit jobs on brasdor

run from within a TEMPK/1.rs/ directory (i.e. 1000K/1.40/)
"""
import os, sys, commands, glob, random

def main():

    # Retrieve user input
    try:
        fname = sys.argv[1]
        f = open(fname,'r')
        lines = f.readlines()
        natom = int(lines[0].strip())
        nConfigs = int( len(lines) / (natom+2.0) )

        alat = sys.argv[2]
        ay_over_ax = sys.argv[3]
        az_over_ax = sys.argv[4]
        nruns = int(sys.argv[5])
        launch = int(sys.argv[6])
    except:
        print '\n usage: '+sys.argv[0]+' wrapped.xyz alat(Ang) ay/ax az/ax nruns starting_index_for_dirname\n'
        sys.exit(0)

    # Generate list of new directories for calculations
    dir_nums = range(launch,launch+nruns)
    dirs = ['0'+str(d) for d in dir_nums]

    # Create existing list of directories for checks
    old = glob.glob('0*')

    # Proceed with random selections and submissions
    cwd = os.getcwd()
    for d in dirs:
        if d not in old:
            os.mkdir(d)
            os.system('cp -rf ../template/* '+d)
            os.chdir(d)
            k = 1
            while k != 0:
                snapshot = int( (0.10 + 0.90*random.random()) * nConfigs )
                k = os.system('select_snapshot.py ../'+fname+' '+str(snapshot)+' >& select_snapshot.out')
            os.system('rm -f select_snapshot.out')
            os.system('POSCAR_from_xyz_general.py snapshot.xyz '+alat+' '+ay_over_ax+' '+az_over_ax+' y')
            os.chdir(cwd)
        else:
            print '\n Directory', d, 'already exists - skipping...\n'
    os.system('sleep 61s')
    for d in dirs:
        if d not in old:
            os.chdir(d)
            os.system('qsub vasp.brasdor.s')
            os.chdir(cwd)


if __name__ == '__main__':
    main()
