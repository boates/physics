#!/usr/bin/env python
"""
md-monatomic_analyze.py VERSION 1.0
Author: Brian Boates

Run analysis on a monatomic system
"""
import os, sys, commands, glob, time

global tFinal
tFinal = 5  # in picoseconds

def main():
    """
    Execute several vasp md analysis scripts
    """
    print 'Concatenating VASP files...'

    # Move most recent MD files (if present)
    if 'OUTCAR' in glob.glob('*'):
        os.system('vasp_mover.sh')
    os.system('vasp_cleanup')

    # Concatenate XDATCAR files
    cwd = os.getcwd()
    os.chdir('output/')
    os.system('nice -15 XDATCAR_cat.py')
    os.system('mv -f XDATCAR.cat ../XDATCAR')
        
    # Concatenate the OUTCAR files
    os.system('nice -15 cat OUTCAR.* > OUTCAR')
    os.system('mv -v OUTCAR ../OUTCAR')

    # Concatenate the OSZICAR files
    os.system('nice -15 cat OSZICAR.* > OSZICAR')
    os.system('mv -v OSZICAR ../OSZICAR')

    os.chdir(cwd)
    print 'done...'

    # Check to make sure all the needed vasp files are present
    if ('OUTCAR' and 'XDATCAR' and 'OSZICAR' and 'POSCAR' and 'INCAR') not in glob.glob('*'):
        print '\nCould not find one or more of: OUTCAR, XDATCAR, OSZICAR, POSCAR, INCAR\n'
	sys.exit(0)

    # Remove any previous analysis
    if 'analysis/' in glob.glob('*/'):
        os.system('rm -rf analysis/')

    # Retrieve timestep (fs), natom, & lattice constant
    step_in_fs = float(commands.getoutput("grep POTIM INCAR | awk '{print $3}'"))
    step_in_au = step_in_fs / 2.41880e-02
    natom = int(commands.getoutput("head -6 POSCAR | tail -1"))
    alat = float(commands.getoutput("head -2 POSCAR | tail -1"))  ### In angstroms
    a = float(commands.getoutput("head -3 POSCAR | tail -1").split()[0])*alat
    b = float(commands.getoutput("head -4 POSCAR | tail -1").split()[1])*alat
    c = float(commands.getoutput("head -5 POSCAR | tail -1").split()[2])*alat

    # Determine nFinalConfigs
    nFinalConfigs = int( round( tFinal / (step_in_fs/1000.) ) )
            
    # Extract thermodynamic quantities
    print 'Extracting thermodynamic variables and averaging...'
    os.system('nice -15 thermodynamic_vasp.py OSZICAR OUTCAR')
    
    os.system('tail -'+str(nFinalConfigs)+' toten.dat    > e.dat')
    os.system('tail -'+str(nFinalConfigs)+' energy.dat   > E.dat')
    os.system('tail -'+str(nFinalConfigs)+' pressure.dat > p.dat')

    os.system('blocker e.dat > toten.blocker')
    os.system('blocker E.dat > energy.blocker')
    os.system('blocker p.dat > pressure.blocker')

    os.system('rm -f e.dat E.dat p.dat')
    print 'done'

    # Extract the xyz file
    print 'Extracting xyz file... '
    os.system('nice -15 xyz_from_vasp.x')
    print 'done'

    # Write an input file for unwrap_PBC.x
    print 'Unwrapping xyz file...'
    out = open('unwrap.in','w')
    out.write('TRAJEC.xyz\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.close()
    os.system('nice -15 unwrap_PBC.x < unwrap.in > unwrap.out')
    os.system('mv TRAJEC.xyz wrapped.xyz')
    os.system('mv unwrapped.xyz TRAJEC.xyz')
    os.system('rm -f unwrap.in unwrap.out')
    print 'done'

    # Determine the number of steps
    nConfigs = int(commands.getoutput('wc -l TRAJEC.xyz').split()[0]) / (natom+2)

    # Take final "tFinal" picoseconds of snapshots
    fnameFinal = 'FINAL-'+str(tFinal)+'ps.xyz'
    nFinalLines = (natom + 2)*nFinalConfigs
    os.system('tail -'+str(nFinalLines)+' TRAJEC.xyz > '+fnameFinal)

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
    os.system('nice -15 correlators.py TRAJEC.cnn 1  1000 >  cluster1.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 2  1000 >  cluster2.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 3  1000 >  cluster3.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 4  1000 >  cluster4.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 5  1000 >  cluster5.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 6  1000 >  cluster6.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 7  1000 >  cluster7.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 8  1000 >  cluster8.out')
#    os.system('nice -15 correlators.py TRAJEC.cnn 9  1000 >  cluster9.out')
#    os.system('nice -15 correlators.py TRAJEC.cnn 10 1000 > cluster10.out')
#    os.system('nice -15 correlators.py TRAJEC.cnn 11 1000 > cluster11.out')
#    os.system('nice -15 correlators.py TRAJEC.cnn 12 1000 > cluster12.out')
#    os.system('nice -15 correlators.py TRAJEC.cnn 13 1000 > cluster13.out')
#    os.system('nice -15 correlators.py TRAJEC.cnn 14 1000 > cluster14.out')
#    os.system('nice -15 correlators.py TRAJEC.cnn 15 1000 > cluster15.out')
#    os.system('nice -15 correlators.py TRAJEC.cnn 16 1000 > cluster16.out')
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
#    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_9nn.dat  > c09.dat")
#    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_10nn.dat > c10.dat")
#    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_11nn.dat > c11.dat")
#    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_12nn.dat > c12.dat")
#    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_13nn.dat > c13.dat")
#    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_14nn.dat > c14.dat")
#    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_15nn.dat > c15.dat")
#    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_16nn.dat > c16.dat")
    os.system('rm -f correlator_*nn.dat')

    # Grab the correlators at given timesteps
    os.system('grab_correlators_at_timestep_CO2.sh 100')
    print 'done'

    # Calculate RDF
    print 'Calculating RDF...'
    out = open('RDF.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF.dat\n')
    out.write('N\nN\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF.in > RDF.out')
    os.system('rm -f RDF.in RDF.out')
    print 'done'

    # Calculate atomic coordination
    print 'Calculating coordinations...'
    min_gr = commands.getoutput('min_gr.py RDF.dat')
    out = open('coord.in','w')
    out.write(fnameFinal+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write(min_gr+'\n')
    out.close()
    os.system('nice -15 coordination.x < coord.in > coord.out')
    os.system('rm -f coord.out')
    print 'done'

    # Calculate VACF
    print 'Calculating VACF...'
    out = open('VACF.in','w')
    out.write(fnameFinal+'\n')
    out.write('VACF.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('ALL\n')
    out.close()
    os.system('nice -15 VACF < VACF.in > VACF.out')
    os.system('mv -f DiffCoef.data DiffCoef.dat')
    os.system('rm -f VACF.in VACF.out')
    print 'done'

    # Construct VDOS input file
    print 'Calculating VDOS...'
    out = open('VDOS.in','w')
    out.write('VACF.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS.in > VDOS.out')
    os.system('vdos_to_wavenumber.sh VDOS.dat')
    os.system('running_avg.py VDOS.dat 2 10')
    os.system('rm -f VDOS.in VDOS.out')

    # Calculate msd
    print 'Calculating msd, MSD, and diffusion...'
    out = open('msd.in','w')
    out.write('TRAJEC.xyz\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')

    # Calculate MSD & diffusion
    out = open('MSD.in','w')
    out.write('TRAJEC.xyz\n')
    out.write('MSD.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('100\n')
    out.write('5\n')
    out.write(' \n')
    out.close()
    os.system('nice -15 MSD < MSD.in > MSD.out')
    os.system('rm -f MSD.in MSD.out')
    os.system('nice -15 diffusion.py MSD.dat')
    os.system('mv -f D.dat diffusion.dat')
    os.system('rm -f fort.* diss.dat')
    print 'done'

    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')
    os.system('mkdir analysis')
    os.system('mv -f nohup.out analysis/log')
    os.system('mv -f *.dat *.xyz *.cnn *.mol *.in *.blocker *.hist *.avg date analysis/')

    os.system('rm -f OUTCAR OSZICAR XDATCAR')

    # Copy any relevant plotting scripts to analysis/
    os.system('cp /home/boates/software/gnuplot_scripts/c1-16.plt analysis/')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
