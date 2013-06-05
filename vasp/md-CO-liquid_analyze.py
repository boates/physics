#!/usr/bin/env python
"""
md-CO2-vasp_analyze.py VERSION 4.0
Author: Brian Boates

Run all applicable scripts to analyze a CO2 system
"""
import os, sys, commands, glob, time

global massC, massO
massC = 12.0107  # in amu
massO = 15.9994

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
    nC = int(commands.getoutput("head -6 POSCAR | tail -1").split()[0])
    nO = int(commands.getoutput("head -6 POSCAR | tail -1").split()[1])
    natom = nC + nO
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
    os.system('mv TRAJEC.xyz wrapped_CO2.xyz')
    os.system('mv unwrapped.xyz TRAJEC_CO2.xyz')
    print 'done'

    # Extract seperate xyz files for C and O
    print 'Separating xyz files for C and O...'
    out = open('separate_wrapped.in','w')
    out.write('wrapped_CO2.xyz\n')
    out.write(str(nC)+','+str(nO)+'\n')
    out.write('C,O\n')
    out.close()

    os.system('separate_xyz_two_types.x < separate_wrapped.in > separate.out')
    os.system('mv TRAJEC_type1.xyz wrapped_C.xyz')
    os.system('mv TRAJEC_type2.xyz wrapped_O.xyz')

    out = open('separate.in','w')
    out.write('TRAJEC_CO2.xyz\n')
    out.write(str(nC)+','+str(nO)+'\n')
    out.write('C,O\n')
    out.close()

    os.system('separate_xyz_two_types.x < separate.in > separate.out')
    os.system('mv TRAJEC_type1.xyz TRAJEC_C.xyz')
    os.system('mv TRAJEC_type2.xyz TRAJEC_O.xyz')

    os.system('rm -f separate_wrapped.in separate.in separate.out')
    print 'done'

    # Determine the number of steps
    nConfigs = int(commands.getoutput('wc -l TRAJEC_CO2.xyz').split()[0])/(natom+2)

    # Slice the xyz lengths for analysis for CO2
    fnameFinal_CO2 = 'FINAL-'+str(tFinal)+'ps_CO2.xyz'
    nFinalLines_CO2 = (natom + 2)*nFinalConfigs
    os.system('tail -'+str(nFinalLines_CO2)+' TRAJEC_CO2.xyz > '+fnameFinal_CO2)

    # Slice the xyz lengths for analysis for C
    fnameFinal_C = 'FINAL-'+str(tFinal)+'ps_C.xyz'
    nFinalLines_C = (nC + 2)*nFinalConfigs
    os.system('tail -'+str(nFinalLines_C)+' TRAJEC_C.xyz > '+fnameFinal_C)

    # Slice the xyz lengths for analysis for O
    fnameFinal_O = 'FINAL-'+str(tFinal)+'ps_O.xyz'
    nFinalLines_O = (nO + 2)*nFinalConfigs
    os.system('tail -'+str(nFinalLines_O)+' TRAJEC_O.xyz > '+fnameFinal_O)

    # Create cnn and nn_average.hist files
    print 'Extracting cnn files...'
    out = open('cnn_CO2.in','w')
    out.write(fnameFinal_CO2+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write('200\n')
    out.write('16\n')
    out.close()
    os.system('nice -15 nn_dist.x < cnn_CO2.in > nn_dist.out')
    os.system('mv TRAJEC.cnn TRAJEC_CO2.cnn')
    os.system('mv nn_dist.hist nn_dist_CO2.hist')

    out = open('cnn_C.in','w')
    out.write(fnameFinal_C+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write('200\n')
    out.write('16\n')
    out.close()
    os.system('nice -15 nn_dist.x < cnn_C.in > nn_dist.out')
    os.system('mv TRAJEC.cnn TRAJEC_C.cnn')
    os.system('mv nn_dist.hist nn_dist_C.hist')

    out = open('cnn_O.in','w')
    out.write(fnameFinal_O+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write('200\n')
    out.write('16\n')
    out.close()
    os.system('nice -15 nn_dist.x < cnn_O.in > nn_dist.out')
    os.system('mv TRAJEC.cnn TRAJEC_O.cnn')
    os.system('mv nn_dist.hist nn_dist_O.hist')

    out = open('cnn_CO.in','w')
    out.write(fnameFinal_CO2+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write('200\n')
    out.write('16\n')
    out.write('C,O\n')
    out.write(str(nC)+','+str(nO)+'\n')
    out.close()
    os.system('nice -15 nn_dist_two_species.x < cnn_CO.in > nn_dist.out')
    os.system('mv TRAJEC.cnn TRAJEC_CO.cnn')
    os.system('mv nn_dist.hist nn_dist_CO.hist')

    os.system('rm -f cnn_CO2.in cnn_C.in cnn_O.in cnn_CO.in nn_dist.out')
    print 'done'

    # Bonding lifetime analysis
#    print 'Calculating cluster lifetimes for COx...'
#    if nFinalConfigs/2 >= nConfigs:
#        var = str(nConfigs - 1)
#    else:
#        var = str(nFinalConfigs/2)
#    os.system('nice -15 correlators.py TRAJEC_CO.cnn 1 '+var+' > cluster1.out')
#    os.system('nice -15 correlators.py TRAJEC_CO.cnn 2 '+var+' > cluster2.out')
#    os.system('nice -15 correlators.py TRAJEC_CO.cnn 3 '+var+' > cluster3.out')
#    os.system('nice -15 correlators.py TRAJEC_CO.cnn 4 '+var+' > cluster4.out')
#    os.system('rm -f cluster*.out')

    # Convert from tstep to picoseconds
#    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_1nn.dat > CO1.dat")
#    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_2nn.dat > CO2.dat")
#    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_3nn.dat > CO3.dat")
#    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_4nn.dat > CO4.dat")
#    os.system('rm -f correlator_*nn.dat')

    # Grab the correlators at given timesteps
#    os.system('grab_correlators_at_timestep_CO2.sh 100')
#    os.system('grab_correlators_at_timestep_CO2.sh 1000')
#    print 'done'

    # Calculate RDF's
    print 'Calculating RDFs...'
    out = open('RDF_CC.in','w')
    out.write(fnameFinal_CO2+'\n')
    out.write('RDF_CC.dat\n')
    out.write('C\nC\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_CC.in > RDF_CC.out')

    out = open('RDF_CO.in','w')
    out.write(fnameFinal_CO2+'\n')
    out.write('RDF_CO.dat\n')
    out.write('C\nO\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_CO.in > RDF_CO.out')

    out = open('RDF_OO.in','w')
    out.write(fnameFinal_CO2+'\n')
    out.write('RDF_OO.dat\n')
    out.write('O\nO\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_OO.in > RDF_OO.out')

    os.system('rm -f RDF_CC.in RDF_CO.in RDF_OO.in')
    os.system('rm -f RDF_CC.out RDF_CO.out RDF_OO.out')
    print 'done'

    # Calculate coordination's
    print 'Calculating coordinations...'
    min_gr = commands.getoutput('min_gr.py RDF_CO.dat')
    out = open('coord_CO2.in','w')
    out.write(fnameFinal_CO2+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write(min_gr+'\n')
    out.close()
    os.system('nice -15 coordination.x < coord_CO2.in > coord.out')
    os.system('mv -f coordination.xyz coordination_CO2.xyz')
    os.system('mv -f coordination.dat coordination_CO2.dat')

    min_gr = commands.getoutput('min_gr.py RDF_CC.dat')
    out = open('coord_C.in','w')
    out.write(fnameFinal_C+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write(min_gr+'\n')
    out.close()
    os.system('nice -15 coordination.x < coord_C.in > coord.out')
    os.system('mv -f coordination.xyz coordination_C.xyz')
    os.system('mv -f coordination.dat coordination_C.dat')

    min_gr = commands.getoutput('min_gr.py RDF_OO.dat')
    out = open('coord_O.in','w')
    out.write(fnameFinal_O+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write(min_gr+'\n')
    out.close()
    os.system('nice -15 coordination.x < coord_O.in > coord.out')
    os.system('mv -f coordination.xyz coordination_O.xyz')
    os.system('mv -f coordination.dat coordination_O.dat')

    os.system('rm -f coord_CO2.in coord_C.in coord_O.in coord.out')

    # Calculating "two_types" coordination: C wrt O only and O wrt C only
    min_gr = commands.getoutput('min_gr.py RDF_CO.dat')
    out = open('coord_CO.in','w')
    out.write(fnameFinal_CO2+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write(min_gr+'\n')
    out.write('C,O\n')
    out.write(str(nC)+','+str(nO)+'\n')
    out.close()
    os.system('nice -15 coordination_two_types.x < coord_CO.in > coord.out')
    os.system('mv -f coordination_two_types.xyz coordination_CwrtO.xyz')
    os.system('mv -f coordination_two_types.dat coordination_CwrtO.dat')

    out = open('coord_OC.in','w')
    out.write(fnameFinal_CO2+'\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.write(min_gr+'\n')
    out.write('O,C\n')
    out.write(str(nO)+','+str(nC)+'\n')
    out.close()
    os.system('nice -15 coordination_two_types.x < coord_OC.in > coord.out')
    os.system('mv -f coordination_two_types.xyz coordination_OwrtC.xyz')
    os.system('mv -f coordination_two_types.dat coordination_OwrtC.dat')

    os.system('rm -f coord_CO.in coord_OC.in coord.out')
    print 'done'

    # Calculate VACF's
    print 'Calculating VACFs...'
    out = open('VACF.in','w')
    out.write(fnameFinal_CO2+'\n')
    out.write('VACF_CO2.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('ALL\n')
    out.close()
    os.system('nice -15 VACF < VACF.in > VACF.out')
    os.system('rm -f VACF.in VACF.out')
    os.system('mv -f DiffCoef.data DiffCoef_CO2.dat')

    out = open('VACF.in','w')
    out.write(fnameFinal_C+'\n')
    out.write('VACF_C.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('ALL\n')
    out.close()
    os.system('nice -15 VACF < VACF.in > VACF.out')
    os.system('rm -f VACF.in VACF.out')
    os.system('mv -f DiffCoef.data DiffCoef_C.dat')

    out = open('VACF.in','w')
    out.write(fnameFinal_O+'\n')
    out.write('VACF_O.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('ALL\n')
    out.close()
    os.system('nice -15 VACF < VACF.in > VACF.out')
    os.system('rm -f VACF.in VACF.out')
    os.system('mv -f DiffCoef.data DiffCoef_O.dat')
    print 'done'

    # Construct VDOS input file
    print 'Calculating VDOSs...'
    out = open('VDOS.in','w')
    out.write('VACF_C.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS.in > VDOS.out')
    os.system('mv -f VDOS.dat VDOS_C.dat')
    os.system("awk '{print $1,$2*"+str(nC)+"}' VDOS_C.dat > dos_C.dat")
    os.system('vdos_to_wavenumber.sh VDOS_C.dat')
    os.system('running_avg.py VDOS_C.dat 2 10')
    os.system('rm -f VDOS.in VDOS.out')

    out = open('VDOS.in','w')
    out.write('VACF_O.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS.in > VDOS.out')
    os.system('mv -f VDOS.dat VDOS_O.dat')
    os.system("awk '{print $1,$2*"+str(nO)+"}' VDOS_O.dat > dos_O.dat")
    os.system('vdos_to_wavenumber.sh VDOS_O.dat')
    os.system('running_avg.py VDOS_O.dat 2 10')
    os.system('rm -f VDOS.in VDOS.out')
    print 'done'

    # Now get the entropy and compute free energies
    out = open('fluidicity_C.in','w')
    out.write('dos_C.dat\n')
    vol = str(float(commands.getoutput('cat V').strip())*int(natom))
    out.write(vol+'\n')
    out.write(str(nC)+'\n')
    temp = commands.getoutput('cat T').strip()
    out.write(temp+'\n')
    out.write(str(massC)+'\n')
    out.close()
    os.system("fluidicity.x < fluidicity_C.in > fluidicity_C.out")

    out = open('entropy_C.in','w')
    out.write('dos_C.dat\n')
    out.write(str(natom)+'\n')
    out.write(str(nC)+'\n')
    out.write(vol+'\n')
    out.write(temp+'\n')
    out.write(str(massC)+'\n')
    fluidicity = commands.getoutput('tail -n-1 fluidicity_C.out').split()[-1]
    out.write(fluidicity+'\n')
    out.close()
    os.system('entropy_liquid.x < entropy_C.in > entropy_C.out')

    os.system("grep 'Entropy per atom in eV/K' entropy_C.out | awk '{print $6}' > S_C") # in eV/K
    os.system("grep 'Entropy per atom in k_B'  entropy_C.out | grep -v id | awk '{print $6}' > S_C_in_kB") # in kB

    out = open('fluidicity_O.in','w')
    out.write('dos_O.dat\n')
    vol = str(float(commands.getoutput('cat V').strip())*int(natom))
    out.write(vol+'\n')
    out.write(str(nO)+'\n')
    temp = commands.getoutput('cat T').strip()
    out.write(temp+'\n')
    out.write(str(massO)+'\n')
    out.close()
    os.system("fluidicity.x < fluidicity_O.in > fluidicity_O.out")

    out = open('entropy_O.in','w')
    out.write('dos_O.dat\n')
    out.write(str(natom)+'\n')
    out.write(str(nO)+'\n')
    out.write(vol+'\n')
    out.write(temp+'\n')
    out.write(str(massO)+'\n')
    fluidicity = commands.getoutput('tail -n-1 fluidicity_O.out').split()[-1]
    out.write(fluidicity+'\n')
    out.close()
    os.system('entropy_liquid.x < entropy_O.in > entropy_O.out')

    os.system("grep 'Entropy per atom in eV/K' entropy_O.out | awk '{print $6}' > S_O") # in eV/K
    os.system("grep 'Entropy per atom in k_B'  entropy_O.out | grep -v id | awk '{print $6}' > S_O_in_kB") # in kB
    ######
    os.system("paste S_C S_O | awk '{print $1+$2}' > S")
    os.system("paste S_C_in_kB S_O_in_kB | awk '{print $1+$2}' > S_in_kB")
    ######

    os.system("paste T S | awk '{print $1*$2}' > TS")  # eV/atom
    os.system("paste P V | awk '{print $1*$2*0.0062415097}' > PV")  # eV/atom
    os.system("paste E PV | awk '{print $1+$2}' > H")  # eV/atom
    os.system("paste E TS | awk '{print $1-$2}' > F")  # eV/atom
    os.system("paste E PV TS | awk '{print $1+$2-$3}' > G")  # eV/atom
    
    # Calculate msd's
    print 'Calculating msds, MSDs, and diffusion...'
    out = open('msd.in','w')
    out.write('TRAJEC_CO2.xyz\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')
    os.system('mv -f msd.dat msd_CO2.dat')
    os.system('rm -f msd.in msd.out')

    out = open('msd.in','w')
    out.write('TRAJEC_C.xyz\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')
    os.system('mv -f msd.dat msd_C.dat')
    os.system('rm -f msd.in msd.out')

    out = open('msd.in','w')
    out.write('TRAJEC_O.xyz\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')
    os.system('mv -f msd.dat msd_O.dat')
    os.system('rm -f msd.in msd.out')

    # Calculate MSD's & diffusion
    out = open('MSD.in','w')
    out.write('TRAJEC_CO2.xyz\n')
    out.write('MSD_CO2.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('100\n')
    out.write('5\n')
    out.write(' \n')
    out.close()
    os.system('nice -15 MSD < MSD.in > MSD.out')
    os.system('rm -f MSD.in MSD.out')
    os.system('nice -15 diffusion.py MSD_CO2.dat')
    os.system('mv -f D.dat diffusion_CO2.dat')

    out = open('MSD.in','w')
    out.write('TRAJEC_C.xyz\n')
    out.write('MSD_C.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('100\n')
    out.write('5\n')
    out.write(' \n')
    out.close()
    os.system('nice -15 MSD < MSD.in > MSD.out')
    os.system('rm -f MSD.in MSD.out')
    os.system('nice -15 diffusion.py MSD_C.dat')
    os.system('mv -f D.dat diffusion_C.dat')

    out = open('MSD.in','w')
    out.write('TRAJEC_O.xyz\n')
    out.write('MSD_O.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('100\n')
    out.write('5\n')
    out.write(' \n')
    out.close()
    os.system('nice -15 MSD < MSD.in > MSD.out')
    os.system('rm -f MSD.in MSD.out')
    os.system('nice -15 diffusion.py MSD_O.dat')
    os.system('mv -f D.dat diffusion_O.dat')

    os.system('rm -f fort.* diss.dat')
    print 'done'

    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')
    os.system('mkdir analysis')
    os.system('mv -f nohup.out analysis/log')
    os.system('mv -f *.dat *.xyz *.cnn *.in *.hist *.avg *.out date analysis/')
    os.system('mv -f toten E P V T S S_* PV TS H F G analysis/')

    os.system('rm -f OUTCAR OSZICAR XDATCAR')

    # Copy any relevant plotting scripts to analysis/
    os.system('cp /home/boates/software/CO2/gnuplot_scripts_CO2/* analysis/')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
