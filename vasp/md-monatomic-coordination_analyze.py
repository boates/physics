#!/usr/bin/env python
"""
md-monatomic-coordination_analyze.py VERSION 1.0
Author: Brian Boates

Run coordination analysis on a monatomic system
"""
import os, sys, commands, glob, time

global tFinal
tFinal = 3  # in picoseconds

def main():
    """
    Execute several vasp md analysis scripts
    """
    print 'Concatenating VASP files...'

    # Concatenate XDATCAR files
    cwd = os.getcwd()
    os.chdir('output/')
    os.system('gunzip XDATCAR.*') # unzip
    os.system('nice -15 XDATCAR_cat.py')
    os.system('mv -f XDATCAR.cat ../XDATCAR')
    os.system('gzip XDATCAR.*') # rezip
    os.chdir(cwd)
    print 'done...'

    # Check to make sure all the needed vasp files are present
    if ('XDATCAR' and 'POSCAR' and 'INCAR') not in glob.glob('*'):
        print '\nCould not find one or more of:, XDATCAR, POSCAR, INCAR\n'
	sys.exit(1)

    # Remove any previous analysis
    if 'coord/' in glob.glob('*/'):
        os.system('rm -rf coord/')

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
    print 'Calculating coordination...'
    min_gr = commands.getoutput('min_gr.py RDF.dat')
    out = open('coord.in','w')
    out.write(fnameFinal+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write(min_gr+'\n')
    out.close()
    os.system('nice -15 coordination.x < coord.in > coord.out')
    os.system('rm -f coord.out')
    print 'done'

    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')
    os.system('mkdir coord')
    os.system('rm -f TRAJEC.xyz')
    os.system('gzip *xyz')
    os.system('mv -f *.dat *xyz.gz *.in *.out *.avg date coord/')

    os.system('rm -f XDATCAR')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
