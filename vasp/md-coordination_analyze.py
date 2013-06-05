#!/usr/bin/env python
"""
md-hugoniot-analyze.py
Author: Brian Boates

Run analysis on a monatomic system to get
TD variables to compute Hugoniot
"""
import os, sys, commands, glob, time

global tFinal
tFinal = 5  # in picoseconds

def main():
    """
    Execute several vasp md analysis scripts
    """
    print 'Concatenating VASP files...'

    # Concatenate XDATCAR files
    cwd = os.getcwd()
    os.chdir('output/')
    os.system('nice -15 XDATCAR_cat.py')
    os.system('mv -f XDATCAR.cat ../XDATCAR')
    os.chdir(cwd)
    print 'done...'

    # Check to make sure all the needed vasp files are present
    if ('XDATCAR' and 'POSCAR' and 'INCAR') not in glob.glob('*'):
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

    # Determine nFinalConfigs
    nFinalConfigs = int( round( tFinal / (step_in_fs/1000.) ) )
            
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

    # Calculate msd
    print 'Calculating msd, MSD, and diffusion...'
    out = open('msd.in','w')
    out.write('TRAJEC.xyz\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')

    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')
    os.system('mkdir analysis')
    os.system('mv -f *.dat *.xyz *.in  *.avg date analysis/')

    os.system('rm -f XDATCAR')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
