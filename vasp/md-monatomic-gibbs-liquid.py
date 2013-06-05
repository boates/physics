#!/usr/bin/env python
"""
md-monatomic-gibbs-liquid.py VERSION 1.0
Author: Brian Boates

compute Gibbs for a monatomic liquid
"""
import os, sys, commands, glob, time

global mass
mass = 14.00674  # in amu

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
    if 'gibbs/' in glob.glob('*/'):
        os.system('rm -rf gibbs/')

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

    os.system("blocker e.dat | grep Final | awk '{print $4}' > toten")
    os.system("blocker E.dat | grep Final | awk '{print $4/"+str(natom)+"}' > E")  # eV/atom
    os.system("blocker p.dat | grep Final | awk '{print $4}' > P")  # GPa
    os.system("volume_from_POSCAR.py POSCAR | awk '{print $2}' > V")  # A^3 per atom
    os.system("grep TEBEG INCAR | awk '{print $3}' > T")  # K

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

    # create dos.dat in THz normalized to 3N for entropy
    os.system("awk '{print $1,$2*"+str(natom)+"}' VDOS.dat > dos.dat")
    
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
    print 'done'

    # Now get the entropy and compute free energies
    out = open('fluidicity.in','w')
    out.write('dos.dat\n')
    vol = str(float(commands.getoutput('cat V').strip())*int(natom))
    out.write(vol+'\n')
    out.write(str(natom)+'\n')
    temp = commands.getoutput('cat T').strip()
    out.write(temp+'\n')
    out.write(str(mass)+'\n')
    out.close()
    os.system("fluidicity.x < fluidicity.in > fluidicity.out")

    out = open('entropy.in','w')
    out.write('dos.dat\n')
    out.write(str(natom)+'\n')
    out.write(str(natom)+'\n')
    out.write(vol+'\n')
    out.write(temp+'\n')
    out.write(str(mass)+'\n')
    fluidicity = commands.getoutput('tail -n-1 fluidicity.out').split()[-1]
    out.write(fluidicity+'\n')
    out.close()
    os.system('entropy_liquid.x < entropy.in > entropy.out')

    os.system("grep 'Entropy per atom in eV/K' entropy.out | awk '{print $6}' > S") # in eV/K
    os.system("grep 'Entropy per atom in k_B'  entropy.out | awk '{print $6}' > S_in_kB") # in kB

    os.system("paste T S | awk '{print $1*$2}' > TS")  # eV/atom
    os.system("paste P V | awk '{print $1*$2*0.0062415097}' > PV")  # eV/atom
    os.system("paste E PV | awk '{print $1+$2}' > H")  # eV/atom
    os.system("paste E TS | awk '{print $1-$2}' > F")  # eV/atom
    os.system("paste E PV TS | awk '{print $1+$2-$3}' > G")  # eV/atom
    
    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')
    os.system('mkdir gibbs')
    os.system('mv -f nohup.out gibbs/log')
    os.system('mv -f *.dat *.xyz *.in *.avg *.out date toten E P V T S S_in_kB PV TS H F G gibbs/')

    os.system('rm -f OUTCAR OSZICAR XDATCAR')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
