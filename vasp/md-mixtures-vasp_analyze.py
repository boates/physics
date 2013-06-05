#!/usr/bin/env python
"""
md-mixtures-vasp_analyze.py VERSION 3.0
Author: Brian Boates

Run all applicable scripts to analyze a binary mixture
"""
import os, sys, commands, glob, time

global massN, massX
massN = 14.00674 # in amu
#massX = 30.973762 # P
#massX = 12.0107 # C
massX = 28.0855 # Si
#massX = 72.64 # Ge


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
    os.system('cp output/CONTCAR.100 ./POSCAR')  # Ensure a POSCAR is present

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
    nX = int(commands.getoutput("head -6 POSCAR | tail -1").split()[0])
    nN = int(commands.getoutput("head -6 POSCAR | tail -1").split()[1])
    natom = nN + nX
    typat1, typat2 = 'X', 'N'
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

    os.system("blocker e.dat | grep Final | awk '{print $4}' > toten")
    os.system("blocker E.dat | grep Final | awk '{print $4/"+str(natom)+"}' > E")  # eV/atom
    os.system("blocker p.dat | grep Final | awk '{print $4}' > P")  # GPa
    os.system("volume_from_POSCAR.py POSCAR | awk '{print $2}' > V")  # A^3 per atom
    os.system("grep TEBEG INCAR | awk '{print $3}' > T")  # K
                    
    os.system('rm -f e.dat E.dat p.dat')
    print 'done'

    # Extract the xyz file
    print 'Extracting xyz file... '
    out = open('xyz.in','w')
    out.write(typat1+'\n'+typat2+'\n')
    out.close()
    os.system('nice -15 xyz_from_vasp_2types.x < xyz.in > xyz.out')
    os.system('rm -f xyz.in xyz.out')
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
    out = open('cnn.in','w')
    out.write(fnameFinal+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write('200\n')
    out.write('16\n')
    out.close()
    os.system('nice -15 nn_dist.x < cnn.in > nn_dist.out')

    out = open('cnn_XN.in','w')
    out.write(fnameFinal+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write('200\n')
    out.write('16\n')
    out.write('X,N\n')
    out.write(str(nX)+','+str(nN)+'\n')
    out.close()
    os.system('nice -15 nn_dist_two_species.x < cnn_XN.in > nn_dist.out')
    os.system('mv TRAJEC.cnn TRAJEC_XN.cnn')
    os.system('mv nn_dist.hist nn_dist_XN.hist')

    os.system('rm -f cnn.in cnn_XN.in nn_dist.out')
    print 'done'

    # Bonding lifetime analysis
    print 'Calculating cluster lifetimes for XN_x...'
    os.system('nice -15 correlators.py TRAJEC_XN.cnn 1 500 > cluster1.out')
    os.system('nice -15 correlators.py TRAJEC_XN.cnn 2 500 > cluster2.out')
    os.system('nice -15 correlators.py TRAJEC_XN.cnn 3 500 > cluster3.out')
    os.system('nice -15 correlators.py TRAJEC_XN.cnn 4 500 > cluster4.out')
    os.system('nice -15 correlators.py TRAJEC_XN.cnn 5 500 > cluster5.out')
    os.system('nice -15 correlators.py TRAJEC_XN.cnn 6 500 > cluster6.out')
    os.system('nice -15 correlators.py TRAJEC_XN.cnn 7 500 > cluster7.out')
    os.system('nice -15 correlators.py TRAJEC_XN.cnn 8 500 > cluster8.out')
    os.system('rm -f cluster*.out')

    # Convert from tstep to picoseconds
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_1nn.dat > XN1.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_2nn.dat > XN2.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_3nn.dat > XN3.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_4nn.dat > XN4.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_5nn.dat > XN5.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_6nn.dat > XN6.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_7nn.dat > XN7.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_8nn.dat > XN8.dat")
    os.system('rm -f correlator_*nn.dat')
    print 'done'

    # Calculate RDF's
    print 'Calculating RDFs...'
    out = open('RDF_NN.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF_NN.dat\n')
    out.write('N\nN\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_NN.in > RDF_NN.out')

    out = open('RDF_NX.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF_NX.dat\n')
    out.write('N\nX\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_NX.in > RDF_NX.out')

    out = open('RDF_XX.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF_XX.dat\n')
    out.write('X\nX\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_XX.in > RDF_XX.out')

    os.system('rm -f RDF_NN.in RDF_NX.in RDF_XX.in')
    os.system('rm -f RDF_NN.out RDF_NX.out RDF_XX.out')
    print 'done'

    # Calculating "two_types" coordination: A wrt B only and B wrt O only
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

    # Calculate VACF's
    print 'Calculating VACF...'
    out = open('VACF_N.in','w')
    out.write(fnameFinal+'\n')
    out.write('VACF_N.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('4\n')
    out.write('N\n')
    out.close()
    os.system('nice -15 VACF < VACF_N.in > VACF_N.out')
    os.system('rm -f VACF_N.in VACF_N.out')
    os.system('mv -f DiffCoef.data DiffCoef_N.dat')

    out = open('VACF_X.in','w')
    out.write(fnameFinal+'\n')
    out.write('VACF_X.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('4\n')
    out.write('X\n')
    out.close()
    os.system('nice -15 VACF < VACF_X.in > VACF_X.out')
    os.system('rm -f VACF_X.in VACF_X.out')
    os.system('mv -f DiffCoef.data DiffCoef_X.dat')
    print 'done'

    # Construct VDOS input file
    print 'Calculating VDOS...'
    out = open('VDOS_N.in','w')
    out.write('VACF_N.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS_N.in > VDOS_N.out')
    os.system('mv VDOS.dat VDOS_N.dat')

    # create dos.dat in THz normalized to 3N for entropy
    os.system("awk '{print $1,$2*"+str(nN)+"}' VDOS_N.dat > dos_N.dat")

    os.system('vdos_to_wavenumber.sh VDOS_N.dat')
    os.system('rm -f VDOS_N.in VDOS_N.out')

    out = open('VDOS_X.in','w')
    out.write('VACF_X.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS_X.in > VDOS_X.out')
    os.system('mv VDOS.dat VDOS_X.dat')

    # create dos.dat in THz normalized to 3N for entropy
    os.system("awk '{print $1,$2*"+str(nX)+"}' VDOS_X.dat > dos_X.dat")

    os.system('vdos_to_wavenumber.sh VDOS_X.dat')
    os.system('rm -f VDOS_X.in VDOS_X.out')

    fN = nN / float(natom)
    fX = nX / float(natom)
    os.system("paste VDOS_N.dat VDOS_X.dat | awk '{print $1,$2+$4}' > VDOS_NX.dat")

    os.system('running_avg.py VDOS_N.dat 2 10')
    os.system('running_avg.py VDOS_X.dat 2 10')
    os.system('running_avg.py VDOS_NX.dat 2 10')
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

    # Now get the entropy and compute free energies
    out = open('fluidicity_N.in','w')
    out.write('dos_N.dat\n')
    vol = str(float(commands.getoutput('cat V').strip())*int(natom))
    out.write(vol+'\n')
    out.write(str(nN)+'\n')
    temp = commands.getoutput('cat T').strip()
    out.write(temp+'\n')
    out.write(str(massN)+'\n')
    out.close()
    os.system("fluidicity.x < fluidicity_N.in > fluidicity_N.out")

    out = open('entropy_N.in','w')
    out.write('dos_N.dat\n')
    out.write(str(natom)+'\n')
    out.write(str(nN)+'\n')
    out.write(vol+'\n')
    out.write(temp+'\n')
    out.write(str(massN)+'\n')
    fluidicity = commands.getoutput('tail -n-1 fluidicity_N.out').split()[-1]
    out.write(fluidicity+'\n')
    out.close()
    os.system('entropy_liquid.x < entropy_N.in > entropy_N.out')

    os.system("grep 'Entropy per atom in eV/K' entropy_N.out | awk '{print $6}' > S_N") # in eV/K
    os.system("grep 'Entropy per atom in k_B'  entropy_N.out | grep -v id | awk '{print $6}' > S_N_in_kB") # in kB

    out = open('fluidicity_X.in','w')
    out.write('dos_X.dat\n')
    vol = str(float(commands.getoutput('cat V').strip())*int(natom))
    out.write(vol+'\n')
    out.write(str(nX)+'\n')
    temp = commands.getoutput('cat T').strip()
    out.write(temp+'\n')
    out.write(str(massX)+'\n')
    out.close()
    os.system("fluidicity.x < fluidicity_X.in > fluidicity_X.out")

    out = open('entropy_X.in','w')
    out.write('dos_X.dat\n')
    out.write(str(natom)+'\n')
    out.write(str(nX)+'\n')
    out.write(vol+'\n')
    out.write(temp+'\n')
    out.write(str(massX)+'\n')
    fluidicity = commands.getoutput('tail -n-1 fluidicity_X.out').split()[-1]
    out.write(fluidicity+'\n')
    out.close()
    os.system('entropy_liquid.x < entropy_X.in > entropy_X.out')

    os.system("grep 'Entropy per atom in eV/K' entropy_X.out | awk '{print $6}' > S_X") # in eV/K
    os.system("grep 'Entropy per atom in k_B'  entropy_X.out | grep -v id | awk '{print $6}' > S_X_in_kB") # in kB
    ######
    os.system("paste S_N S_X | awk '{print $1+$2}' > S")
    os.system("paste S_N_in_kB S_X_in_kB | awk '{print $1+$2}' > S_in_kB")
    ######

    os.system("paste T S | awk '{print $1*$2}' > TS")  # eV/atom
    os.system("paste P V | awk '{print $1*$2*0.0062415097}' > PV")  # eV/atom
    os.system("paste E PV | awk '{print $1+$2}' > H")  # eV/atom
    os.system("paste E TS | awk '{print $1-$2}' > F")  # eV/atom
    os.system("paste E PV TS | awk '{print $1+$2-$3}' > G")  # eV/atom
    
    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')
    os.system('mkdir analysis')
    os.system('mv -f nohup.out analysis/log')
    os.system('mv -f *.dat *.xyz *.cnn *.in *.hist *.avg *.out date analysis/')
    os.system('mv -f toten E P V T S S_* PV TS H F G analysis/')

    os.system('rm -f OUTCAR OSZICAR XDATCAR')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
