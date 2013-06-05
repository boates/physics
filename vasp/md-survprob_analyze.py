#!/usr/bin/env python
"""
md-survprob_analyze.py
Author: Brian Boates

Run analysis on a monatomic system to get
survival probabilities
"""
import os, sys, commands, glob, time

global tFinal
tFinal = 5  # in picoseconds

def main():
    """
    Execute several vasp md analysis scripts
    """
    # Take final "tFinal" picoseconds of snapshots
    fnameFinal = 'FINAL-'+str(tFinal)+'ps.xyz'

    os.system('cp analysis/FINAL-5ps.xyz ./')
    
    # Check to make sure all the needed vasp files are present
    if ('POSCAR' and 'INCAR') not in glob.glob('*'):
        print '\nCould not find one or more of: XDATCAR, POSCAR, & INCAR\n'
	sys.exit(0)

    # Retrieve timestep (fs), natom, & lattice constant
    step_in_fs = float(commands.getoutput("grep POTIM INCAR | awk '{print $3}'"))
    step_in_au = step_in_fs / 2.41880e-02
    natom = int(commands.getoutput("head -6 POSCAR | tail -1"))
    alat = float(commands.getoutput("head -2 POSCAR | tail -1"))  ### In angstroms
    a = float(commands.getoutput("head -3 POSCAR | tail -1").split()[0])*alat
    b = float(commands.getoutput("head -4 POSCAR | tail -1").split()[1])*alat
    c = float(commands.getoutput("head -5 POSCAR | tail -1").split()[2])*alat

    # Create cnn and nn_average.hist files
    print 'Extracting cnn and nn_dist file...'
    out = open('cnn.in','w')
    out.write(fnameFinal+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write('200\n')
    out.write('16\n')
    out.close()
    os.system('nice -15 nn_dist.x < cnn.in > nn_dist.out')
    os.system('rm -f cnn.in nn_dist.out')
    print 'done'

    # Determine nn angle distribution
    print 'Calculating angle distributions...'
    out = open('angles.in','w')
    out.write('TRAJEC.cnn\n')
    out.write('1,2\n')
    out.write('200\n')
    out.write('y\n')
    out.close()
    os.system('nn_angles.x < angles.in > angles.out')
    os.system('rm -f angles.in angles.out')
    print 'done'

    # Bonding lifetime analysis
    print 'Calculating cluster lifetimes...'
    os.system('nice -15 correlators.py TRAJEC.cnn 1  500 >  cluster1.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 2  500 >  cluster2.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 3  500 >  cluster3.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 4  500 >  cluster4.out')
    os.system('rm -f cluster*.out')

    # Convert from tstep to picoseconds
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_1nn.dat  > c01.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_2nn.dat  > c02.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_3nn.dat  > c03.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_4nn.dat  > c04.dat")
    os.system('rm -f correlator_*nn.dat')

    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')
#    os.system('mkdir analysis')
    os.system('mv -f *.dat *.xyz *.cnn *.hist date analysis/')

    os.system('rm -f XDATCAR')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
