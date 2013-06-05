#!/usr/bin/env python
"""
md-basic-CO2-vasp_analyze.py VERSION 1.0
Author: Brian Boates

Run some applicable scripts to analyze a CO2 system
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
    os.system('nice -15 vasp_hugoniot_energy.py')
    
    os.system('tail -'+str(nFinalConfigs)+' toten.dat    > e.dat')
    os.system('tail -'+str(nFinalConfigs)+' energy.dat   > E.dat')
    os.system('tail -'+str(nFinalConfigs)+' energy_hugoniot.dat   > hugE.dat')
    os.system('tail -'+str(nFinalConfigs)+' pressure.dat > p.dat')

    os.system('blocker e.dat > toten.blocker')
    os.system('blocker E.dat > energy.blocker')
    os.system('blocker hugE.dat > energy_hugoniot.blocker')
    os.system('blocker p.dat > pressure.blocker')

    os.system('rm -f e.dat E.dat hugE.dat p.dat')
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

    # Determine the number of steps
    nConfigs = int(commands.getoutput('wc -l TRAJEC_CO2.xyz').split()[0])/(natom+2)

    # Slice the xyz lengths for analysis for CO2
    fnameFinal_CO2 = 'FINAL-'+str(tFinal)+'ps_CO2.xyz'
    nFinalLines_CO2 = (natom + 2)*nFinalConfigs
    os.system('tail -'+str(nFinalLines_CO2)+' TRAJEC_CO2.xyz > '+fnameFinal_CO2)

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

    # Calculate msd's
    print 'Calculating msds, MSDs, and diffusion...'
    out = open('msd.in','w')
    out.write('TRAJEC_CO2.xyz\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')
    os.system('mv -f msd.dat msd_CO2.dat')

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
    os.system('rm -f MSD.in MSD.out fort.*')
    os.system('nice -15 diffusion.py MSD_CO2.dat')
    os.system('mv -f D.dat diffusion_CO2.dat')
    print 'done'

    # Calculate VACF's
    print 'Calculating VACFs...'
    out = open('VACF_C.in','w')
    out.write(fnameFinal_CO2+'\n')
    out.write('VACF_C.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('C\n')
    out.close()
    os.system('nice -15 VACF < VACF_C.in > VACF_C.out')
    os.system('rm -f VACF_C.in VACF_C.out')
    os.system('mv -f DiffCoef.data DiffCoef_C.dat')

    out = open('VACF_O.in','w')
    out.write(fnameFinal_CO2+'\n')
    out.write('VACF_O.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('O\n')
    out.close()
    os.system('nice -15 VACF < VACF_O.in > VACF_O.out')
    os.system('rm -f VACF_O.in VACF_O.out')
    os.system('mv -f DiffCoef.data DiffCoef_O.dat')
    print 'done'

    # Calculate VDOS's
    print 'Calculating VDOSs and entropies/free energies...'
    out = open('VDOS_C.in','w')
    out.write('VACF_C.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS_C.in > VDOS_C.out')
    os.system('mv VDOS.dat VDOS_C.dat')

    out = open('VDOS_O.in','w')
    out.write('VACF_O.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS_O.in > VDOS_O.out')
    os.system('mv VDOS.dat VDOS_O.dat')

    os.system('rm -f VDOS_*.in VDOS_*.out')

    # Combine C and O into CO2 VDOS
    os.system("paste VDOS_C.dat VDOS_O.dat | awk '{print $1,$2/3.+2*$4/3.}' > VDOS_CO2.dat")
    os.system("awk '{print $1,$2*"+str(natom)+"}' VDOS_CO2.dat > dos.dat")

    # Get entropy
    out = open('gibbs.in','w')
    out.write('dos.dat\n')
    out.write('POSCAR\n')
    temp = commands.getoutput("grep TEBEG INCAR | awk '{print $3}'")
    p_avg = commands.getoutput("grep Final pressure.blocker | awk '{print $4}'")
    e_avg = commands.getoutput("grep Final energy_hugoniot.blocker | awk '{print $4}'")
    out.write(temp+', '+p_avg+', '+e_avg+'\n')
    out.close()

    os.system('gibbs-solid.x < gibbs.in > gibbs.dat')

    os.system('vdos_to_wavenumber.sh VDOS_C.dat')
    os.system('running_avg.py VDOS_C.dat 2 10')

    os.system('vdos_to_wavenumber.sh VDOS_O.dat')
    os.system('running_avg.py VDOS_O.dat 2 10')

    os.system('vdos_to_wavenumber.sh VDOS_CO2.dat')
    os.system('running_avg.py VDOS_CO2.dat 2 10')

    os.system('mv entropy.dat S')
    os.system("grep Final energy_hugoniot.blocker | awk '{print $4/"+str(natom)+"}' > E")
    os.system("grep Final pressure.blocker | awk '{print $4}' > P")
    os.system("volume_from_POSCAR.py POSCAR | awk '{print $2}' > V")
    os.system("grep TEBEG INCAR | awk '{print $3}' > T")

    os.system("paste P V | awk '{print $1*$2*0.0062415097}' > PV")
    os.system("paste E PV | awk '{print $1+$2}' > H")
    os.system("paste T S | awk '{print $1*$2}' > TS")
    os.system("paste H TS | awk '{print $1-$2}' > G")
                    
    print 'done'

    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')
    os.system('mkdir analysis')
    os.system('mv -f nohup.out analysis/log')
    os.system('mv -f *.dat *.xyz *.in *.blocker *.avg date analysis/')
    os.system('mv -f E P V T S H G PV TS analysis/')

    os.system('rm -f OUTCAR OSZICAR XDATCAR')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
