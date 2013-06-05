#!/usr/bin/env python
"""
md-MgSiO3-vasp_analyze.py VERSION 1.0
Author: Brian Boates

Run all applicable scripts to analyze a MgSiO3 system
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
    if 'OUTCAR' in glob.glob('*') and commands.getoutput('head CONTCAR') != '':
        os.system('vasp_mover.sh')
    os.system('vasp_cleanup; rm -f CONTCAR')

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
    nMg = int(commands.getoutput("grep -B1 Direct POSCAR | grep -v Direct").split()[0])
    nSi = int(commands.getoutput("grep -B1 Direct POSCAR | grep -v Direct").split()[1])
    nO  = int(commands.getoutput("grep -B1 Direct POSCAR | grep -v Direct").split()[2])
    natom = nMg + nSi + nO
    typat1, typat2, typat3 = 'Mg', 'Si', 'O'
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
#    out = open('xyz.in','w')
#    out.write(typat1+'\n'+typat2+'\n'+typat3+'\n')
#    out.close()
    os.system('nice -15 xyz_from_vasp5_3types.x') # < xyz.in > xyz.out')
#    os.system('rm -f xyz.in xyz.out')
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
    print 'done'

    # Determine the number of steps
    nConfigs = int(commands.getoutput('wc -l TRAJEC.xyz').split()[0])/(natom+2)

    # Slice the xyz lengths for analysis
    fnameFinal = 'FINAL-'+str(tFinal)+'ps.xyz'
    nFinalLines = (natom + 2)*nFinalConfigs
    os.system('tail -'+str(nFinalLines)+' TRAJEC.xyz > '+fnameFinal)

    # Create cnn and nn_average.hist files
    print 'Extracting cnn files...'
    """
    out = open('cnn_MgSi.in','w')
    out.write(fnameFinal+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write('200\n')
    out.write('16\n')
    out.write('Mg,Si\n')
    out.write(str(nMg)+','+str(nSi)+'\n')
    out.close()
    os.system('nice -15 nn_dist_two_species.x < cnn_MgSi.in > nn_dist.out')
    os.system('mv TRAJEC.cnn TRAJEC_MgSi.cnn')
    os.system('mv nn_dist.hist nn_dist_MgSi.hist')

    out = open('cnn_MgO.in','w')
    out.write(fnameFinal+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write('200\n')
    out.write('16\n')
    out.write('Mg,O\n')
    out.write(str(nMg)+','+str(nO)+'\n')
    out.close()
    os.system('nice -15 nn_dist_two_species.x < cnn_MgO.in > nn_dist.out')
    os.system('mv TRAJEC.cnn TRAJEC_MgO.cnn')
    os.system('mv nn_dist.hist nn_dist_MgO.hist')

    out = open('cnn_SiO.in','w')
    out.write(fnameFinal+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write('200\n')
    out.write('16\n')
    out.write('Si,O\n')
    out.write(str(nSi)+','+str(nO)+'\n')
    out.close()
    os.system('nice -15 nn_dist_two_species.x < cnn_SiO.in > nn_dist.out')
    os.system('mv TRAJEC.cnn TRAJEC_SiO.cnn')
    os.system('mv nn_dist.hist nn_dist_SiO.hist')
    """

    out = open('cnn.in','w')
    out.write(fnameFinal+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write('200\n')
    out.write('16\n')
    out.close()
    os.system('nice -15 nn_dist.x < cnn.in > nn_dist.out')
                                
    os.system('rm -f cnn.in nn_dist.out')
#    os.system('rm -f cnn.in cnn_MgSi.in cnn_MgO.in cnn_SiO.in nn_dist.out')
    print 'done'

    # Bonding lifetime analysis
    print 'Calculating cluster lifetimes...'
    if 500 >= nConfigs:
        var = str(nConfigs - 1)
    else:
        var = '500'
    os.system('nice -15 correlators.py TRAJEC.cnn 1 '+var+' > cluster1.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 2 '+var+' > cluster2.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 3 '+var+' > cluster3.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 4 '+var+' > cluster4.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 5 '+var+' > cluster5.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 6 '+var+' > cluster6.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 7 '+var+' > cluster7.out')
    os.system('nice -15 correlators.py TRAJEC.cnn 8 '+var+' > cluster8.out')
    os.system('rm -f cluster*.out')

    # Convert from tstep to picoseconds
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_1nn.dat > c01.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_2nn.dat > c02.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_3nn.dat > c03.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_4nn.dat > c04.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_5nn.dat > c05.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_6nn.dat > c06.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_7nn.dat > c07.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_8nn.dat > c08.dat")
    os.system('rm -f correlator_*nn.dat')
    print 'done'

    # Calculate RDF's
    print 'Calculating RDFs...'
    out = open('RDF_MgSi.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF_MgSi.dat\n')
    out.write('Mg\nSi\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_MgSi.in > RDF_MgSi.out')

    out = open('RDF_MgO.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF_MgO.dat\n')
    out.write('Mg\nO\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_MgO.in > RDF_MgO.out')

    out = open('RDF_SiO.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF_SiO.dat\n')
    out.write('Si\nO\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_SiO.in > RDF_SiO.out')

    out = open('RDF_MgMg.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF_MgMg.dat\n')
    out.write('Mg\nMg\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_MgMg.in > RDF_MgMg.out')

    out = open('RDF_SiSi.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF_SiSi.dat\n')
    out.write('Si\nSi\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_SiSi.in > RDF_SiSi.out')

    out = open('RDF_OO.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF_OO.dat\n')
    out.write('O\nO\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_OO.in > RDF_OO.out')

    os.system('rm -f RDF_MgSi.in RDF_MgO.in RDF_SiO.in RDF_MgMg.in RDF_SiSi.in RDF_OO.in')
    os.system('rm -f RDF_MgSi.out RDF_MgO.out RDF_SiO.out RDF_MgMg.out RDF_SiSi.out RDF_OO.out')
    print 'done'

    # Calculating "two_types" coordination: C wrt O only and O wrt C only
    """
    print 'Calculating coordinations...'
    min_gr = commands.getoutput('min_gr.py RDF_NX.dat')
    out = open('coord_NX.in','w')
    out.write(fnameFinal+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write(min_gr+'\n')
    out.write('N,X\n')
    out.write(str(nN)+','+str(nX)+'\n')
    out.close()
    os.system('nice -15 coordination_two_types.x < coord_NX.in > coord.out')
    os.system('mv -f coordination_two_types.xyz coordination_NwrtX.xyz')
    os.system('mv -f coordination_two_types.dat coordination_NwrtX.dat')

    out = open('coord_XN.in','w')
    out.write(fnameFinal+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write(min_gr+'\n')
    out.write('X,N\n')
    out.write(str(nX)+','+str(nN)+'\n')
    out.close()
    os.system('nice -15 coordination_two_types.x < coord_XN.in > coord.out')
    os.system('mv -f coordination_two_types.xyz coordination_XwrtN.xyz')
    os.system('mv -f coordination_two_types.dat coordination_XwrtN.dat')

    os.system('rm -f coord_NX.in coord_XN.in coord.out')
    print 'done'
    """

    # Calculate VACF's
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
    os.system('rm -f VACF.in VACF.out')
    os.system('mv -f DiffCoef.data DiffCoef.dat')
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
    print 'done'

    # Calculate msd's
    print 'Calculating msd, MSD, and diffusion...'
    out = open('msd.in','w')
    out.write('TRAJEC.xyz\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')

    # Calculate MSD's & diffusion
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
    os.system('mv -f *.dat *.xyz *.cnn *.in *.blocker *.hist *.avg date analysis/')

    os.system('rm -f OUTCAR OSZICAR XDATCAR')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
