#!/usr/bin/env python
"""
md-TP-vasp_analyze.py VERSION 1.0
Author: Brian Boates
ANALYZE TWO PHASE MD DATA FROM VASP
  ---> Concatenate all OUTCAR and XDATCAR files together
  ---> Extract T, P, & E using vasp_temperature.py, vasp_pressure.py, & vasp_energy.py
  ---> Extract .xyz file using xyz_from_vasp.x and unwrap using unwrap_PBC.x
  ---> Extract .cbn file using cbn_from_xyz.x
  ---> Separate species into respective .xyz files
  ---> Create .cnn files and calculate nn_average.hist files using nn_dist.py
  ---> Calculate carbon centered angle distributions using CO2_angles.x
  ---> Calculate cluster lifetimes for clusters of 3,4,5,6,7 using bonding_lifetime_n_body.py
  ---> Calculate RDF's using RDF
  ---> Calculate MSD using msd.x
  ---> Calculate VACF's using VACF
  ---> Calculate VDOS's using VDOS.x and convert to cm^-1
  ---> Move all created analysis files to analysis/ directory
"""
import os, sys, commands, glob

global nFinalConfigs, nSkipConfigs
nFinalConfigs = 2500
nSkipConfigs = 1000

def main():
    """
    Execute several vasp md analysis scripts
    """
    if 'output' in glob.glob('*'):

        if 'output/OUTCAR.100' not in glob.glob('output/*'):

            if 'OUTCAR' in glob.glob('*'):
                os.system('mv OUTCAR output/OUTCAR.100')
                os.system('mv XDATCAR output/XDATCAR.100')
                os.system('mv OSZICAR output/OSZICAR.100')

            else:
                print '\nNo output files found ---- exiting...\n'
                sys.exit(0)

        else:
            if 'OUTCAR' in glob.glob('*'):
                os.system('mv OUTCAR output/OUTCAR.199')
                os.system('mv XDATCAR output/XDATCAR.199')
                os.system('mv OSZICAR output/OSZICAR.199')

        print 'Concatenating VASP files...'

        os.system('cp output/XDATCAR.100 ./XDATCAR.100')
        xdats = glob.glob('output/XDATCAR*')
        xdats.pop(xdats.index('output/XDATCAR.100'))
        xdats.sort()
        for xdat in xdats:
            os.system('tail -n+6 '+xdat+' > XDATCAR.'+xdat[-3:])

        # Concatenate the XDATCAR files
        xdats = glob.glob('XDATCAR.*')
        xdats.sort()
        cmd = 'cat '
        for xdat in xdats:
            cmd += xdat+' '
        cmd += '> XDATCAR'
        os.system('nice -15 '+cmd)
        
        # Concatenate the OUTCAR files
        outs = glob.glob('output/OUTCAR*')
        outs.sort()
        cmd = 'cat '
        for out in outs:
            cmd += out+' '
        cmd += '> OUTCAR'
        os.system('nice -15 '+cmd)

        # Remove intermediate data files
        os.system('rm -f XDATCAR.*')

        print 'done...'

        if 'job_files' not in glob.glob('*'):
            os.system('mkdir job_files')    
        os.system('mv -f CO2*K.* job_files/')

    # Check to make sure all the needed vasp files are present
    if ('OUTCAR' and 'INCAR' and 'XDATCAR' and 'POSCAR') not in glob.glob('*'):
        print '\nCould not find file(s): OUTCAR, INCAR, XDATCAR, OSZICAR, POSCAR\n'
	sys.exit(0)

    # Move any previous analysis
    if 'analysis/' in glob.glob('*/'):
        print 'Moving previous analysis... '
        if 'previous_analysis/' in glob.glob('*/'):
            os.system('rm -rf previous_analysis/')
        os.system('mv -f analysis/ previous_analysis')
        print 'done'

    # Retrieve timestep (au), natom, & lattice constant
    step_in_au = int(round(float(commands.getoutput("grep POTIM INCAR | awk '{print $3}'"))/2.41880e-05))/1000
    nC = int(commands.getoutput("head -6 POSCAR | tail -1").split()[0])
    nO = int(commands.getoutput("head -6 POSCAR | tail -1").split()[1])
    natom = nC + nO
    alat = float(commands.getoutput("head -2 POSCAR | tail -1"))  ### In angstroms
    a = float(commands.getoutput("head -3 POSCAR | tail -1").split()[0])*alat
    b = float(commands.getoutput("head -4 POSCAR | tail -1").split()[1])*alat
    c = float(commands.getoutput("head -5 POSCAR | tail -1").split()[2])*alat
            
    # Extract thermodynamic quantities
    print 'Extracting T, P, & E and averaging...'
    os.system('nice -15 vasp_pressure.py')
    os.system('tail -5000 pressure.dat > p.dat')
    os.system('blocker p.dat 2 > pressure.blocker')
    os.system('rm -f p.dat')
    os.system('nice -15 vasp_temperature.py')
    os.system('nice -15 vasp_energy.py')
    print 'done'

    # Create the xyz file from the pwscf output files
    print 'Extracting xyz file... '
    os.system('nice -15 xyz_from_vasp_CO2.x')
    print 'done'

    # Write an input file for unwrap_PBC.x
    print 'Unwrapping xyz file...'
    out = open('unwrap.in','w')
    out.write('TRAJEC.xyz\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.close()
    os.system('nice -15 unwrap_PBC.x < unwrap.in > unwrap.out')
    os.system('rm -f unwrap.out')
    os.system('mv TRAJEC.xyz wrapped.xyz')
    os.system('mv unwrapped.xyz TRAJEC.xyz')

    # Determine the number of steps
    nConfigs = int(commands.getoutput('wc -l TRAJEC.xyz').split()[0])/(natom+2)

    # Take the last nFinalConfigs configurations for further analysis
    fnameFinal = 'FINAL-'+str(nFinalConfigs)+'.xyz'
    nFinalLines = (natom + 2)*nFinalConfigs
    os.system('tail -'+str(nFinalLines)+' TRAJEC.xyz > '+fnameFinal)

    fnameSkip = 'SKIP-'+str(nSkipConfigs)+'.xyz'
    if nConfigs < 10000 + nSkipConfigs:
        if nConfigs > nSkipConfigs*2:
            nSkipLines = (nConfigs - nSkipConfigs)*(natom+2)
        else:
            nSkipLines = nConfigs*(natom+2)
        os.system('tail -'+str(nSkipLines)+' TRAJEC.xyz > '+fnameSkip)
    elif nConfigs >= 10000 + nSkipConfigs:
        fnameSkip = 'FINAL-10000.xyz'
        os.system('tail -'+str(int(10000*(natom+2)))+' TRAJEC.xyz > '+fnameSkip)
    print 'done'

    # Write an input file for cbn_from_xyz.x
    print 'Extracting cbn file...'
    out = open('cbn.in','w')
    out.write(fnameSkip+'\n')
    out.write(str(a/0.529177)+','+str(b/0.529177)+','+str(c/0.529177)+'\n')
    out.close()

    # Create the cbn file
    os.system('nice -15 cbn_from_xyz.x < cbn.in > cbn.out')
    os.system('rm -f cbn.out')
    print 'done'

    # Create cnn file
    print 'Extracting cnn file...'
    os.system('nice -15 nn_dist.py 200 '+str(a/0.529177)+' '+str(a/0.529177)+' '+str(a/0.529177)+' FINAL-'+str(nFinalConfigs)+'.xyz')
    print 'done'

    # Extract separate species and perform separate cnn analysis
    print 'Performing separate species cnn analysis...'
    os.system('separate_CO2_xyz.py FINAL-'+str(nFinalConfigs)+'.xyz '+str(nFinalConfigs))
    os.system('nice -15 nn_dist.py 200 '+str(a/0.529177)+' '+str(a/0.529177)+' '+str(a/0.529177)+' FINAL-'+str(nFinalConfigs)+'_C.xyz C')
    os.system('nice -15 nn_dist.py 200 '+str(a/0.529177)+' '+str(a/0.529177)+' '+str(a/0.529177)+' FINAL-'+str(nFinalConfigs)+'_O.xyz O')
    os.system('nice -15 nn_dist_two_species.py 200 '+str(a/0.529177)+' '+str(a/0.529177)+' '+str(a/0.529177)+' FINAL-'+str(nFinalConfigs)+'.xyz '+str(nO)+' CO')
    print 'done'

    # Calculate all orientational information
    print 'Calculating only carbon cenetered bond angles...'
    out = open('angles.in','w')
    out.write('TRAJEC.cnn\n200\n')
    out.close()
    os.system('CO2_angles.x < angles.in')
    os.system('rm -f angles.in')
    print 'done'

    print 'Performing CO2_orientation.x calculations...'
    os.system('CO2_orientation_run.sh')
    os.system('CO2_angles_norm.py CO2_angles.dat')
    print 'done'
                
    # Calculate lifetimes of clusters of 3, 4, 5, 6, & 7 atoms
    print 'Calculating cluster lifetimes...'
    if nFinalConfigs/2 >= nConfigs:
        var = str(nConfigs - 1)
    else:
        var = str(nFinalConfigs/2)
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC_CO.cnn 2 '+var+' > cluster2.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC_CO.cnn 3 '+var+' > cluster3.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC_CO.cnn 4 '+var+' > cluster4.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC_CO.cnn 5 '+var+' > cluster5.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC_CO.cnn 6 '+var+' > cluster6.out')
    os.system('rm -f cluster*.out')
    print 'done'

    # Construct RDF input file
    print 'Calculating RDF, MSD, VACF, & VDOS... '
    out = open('RDF_CC.in','w')
    out.write(fnameSkip+'\n')
    out.write('RDF_CC.dat\n')
    out.write('C\nC\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()

    out = open('RDF_CO.in','w')
    out.write(fnameSkip+'\n')
    out.write('RDF_CO.dat\n')
    out.write('C\nO\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()

    out = open('RDF_OO.in','w')
    out.write(fnameSkip+'\n')
    out.write('RDF_OO.dat\n')
    out.write('O\nO\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()

    # Construct MSD input file
    out = open('msd.in','w')
    out.write('TRAJEC.xyz\n')
    out.close()

    # Construct VACF input file
    out = open('VACF.in','w')
    out.write(fnameSkip+'\n')
    out.write('VACF.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('ALL\n')
    out.close()

    # Construct VDOS input file
    out = open('VDOS.in','w')
    out.write('VACF.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()

    # Calculate RDF, MSD, VACF, & VDOS
    os.system('nice -15 RDF < RDF_CC.in > RDF_CC.out')
    os.system('nice -15 RDF < RDF_CO.in > RDF_CO.out')
    os.system('nice -15 RDF < RDF_OO.in > RDF_OO.out')
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('nice -15 VACF < VACF.in > VACF.out')
    os.system('VDOS.x < VDOS.in > VDOS.out')
    os.system('vdos_to_wavenumber.sh VDOS.dat')
    os.system('running_avg.py VDOS.dat 2 10')
    os.system('rm -f VDOS.in msd.in fort.* DiffCoef.data diss.dat')
    print 'done'

    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')
    os.system('mkdir analysis')
    os.system('mv -f *.dat *.xyz *.cbn *.cnn *.in *.blocker *.hist *.avg *.norm date analysis/')

    os.system('rm -f OUTCAR XDATCAR')
    os.system('rm -f mini_analysis/')
    os.system('rm -f *.out')

    # Copy any relevant plotting scripts to analysis/
    os.system('cp /home/boates/software/gnuplot_scripts_CO2/* analysis/')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
