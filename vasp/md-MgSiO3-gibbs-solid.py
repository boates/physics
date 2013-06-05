#!/usr/bin/env python
"""
md-MgO-gibbs-solid.py VERSION 1.0
Author: Brian Boates

compute Gibbs for MgO solid
"""
import os, sys, commands, glob, time

global massMg, massSi, massO
massMg = 24.3050  # in amu
massSi = 28.0855
massO  = 15.9994

global tFinal
tFinal = 3  # in picoseconds

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
    nMg = int(commands.getoutput("head -6 POSCAR | tail -1").split()[0])
    nSi = int(commands.getoutput("head -6 POSCAR | tail -1").split()[1])    
    nO = int(commands.getoutput("head -6 POSCAR | tail -1").split()[2])    
    natom = nMg + nSi + nO
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
    out.write('Mg\nSi\nO\n')
    out.close()
    os.system('nice -15 xyz_from_vasp4_3types.x < xyz.in >& xyz.out')
    os.system('rm -f xyz.in xyz.out')
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
    out = open('VACF_Mg.in','w')
    out.write(fnameFinal+'\n')
    out.write('VACF_Mg.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('Mg\n')
    out.close()
    os.system('nice -15 VACF < VACF_Mg.in > VACF_Mg.out')
    os.system('mv -f DiffCoef.data DiffCoef_Mg.dat')
    os.system('rm -f VACF_Mg.in VACF_Mg.out')

    out = open('VACF_Si.in','w')
    out.write(fnameFinal+'\n')
    out.write('VACF_Si.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('Si\n')
    out.close()
    os.system('nice -15 VACF < VACF_Si.in > VACF_Si.out')
    os.system('mv -f DiffCoef.data DiffCoef_Si.dat')
    os.system('rm -f VACF_Si.in VACF_Si.out')

    out = open('VACF_O.in','w')
    out.write(fnameFinal+'\n')
    out.write('VACF_O.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('O\n')
    out.close()
    os.system('nice -15 VACF < VACF_O.in > VACF_O.out')
    os.system('mv -f DiffCoef.data DiffCoef_O.dat')
    os.system('rm -f VACF_O.in VACF_O.out')
    print 'done'

    # Construct VDOS input file
    print 'Calculating VDOS...'
    out = open('VDOS_Mg.in','w')
    out.write('VACF_Mg.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS_Mg.in > VDOS_Mg.out')
    os.system('mv VDOS.dat VDOS_Mg.dat')
    # create dos.dat in THz normalized to 3N for entropy
    os.system("awk '{print $1,$2*"+str(nMg)+"}' VDOS_Mg.dat > dos_Mg.dat")
    os.system('vdos_to_wavenumber.sh VDOS_Mg.dat')
    os.system('running_avg.py VDOS_Mg.dat 2 10')
    os.system('rm -f VDOS_Mg.in VDOS_Mg.out')

    out = open('VDOS_Si.in','w')
    out.write('VACF_Si.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS_Si.in > VDOS_Si.out')
    os.system('mv VDOS.dat VDOS_Si.dat')
    # create dos.dat in THz normalized to 3N for entropy
    os.system("awk '{print $1,$2*"+str(nSi)+"}' VDOS_Si.dat > dos_Si.dat")
    os.system('vdos_to_wavenumber.sh VDOS_Si.dat')
    os.system('running_avg.py VDOS_Si.dat 2 10')
    os.system('rm -f VDOS_Si.in VDOS_Si.out')

    out = open('VDOS_O.in','w')
    out.write('VACF_O.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS_O.in > VDOS_O.out')
    os.system('mv VDOS.dat VDOS_O.dat')
    # create dos.dat in THz normalized to 3N for entropy
    os.system("awk '{print $1,$2*"+str(nO)+"}' VDOS_O.dat > dos_O.dat")
    os.system('vdos_to_wavenumber.sh VDOS_O.dat')
    os.system('running_avg.py VDOS_O.dat 2 10')
    os.system('rm -f VDOS_O.in VDOS_O.out')
    print 'done'
    
    # Calculate msd
    print 'Calculating msd...'
    out = open('msd.in','w')
    out.write('TRAJEC.xyz\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')
    print 'done'

    # Now get the entropy and compute free energies
    out = open('entropy_Mg.in','w')
    out.write('dos_Mg.dat\n')
    out.write('POSCAR\n')
    temp = commands.getoutput('cat T').strip()
    out.write(temp+'\n')
    out.close()
    os.system('entropy_solid_three_types.x < entropy_Mg.in > entropy_Mg.out')
            
    os.system("grep 'Entropy per atom in eV/K' entropy_Mg.out | awk '{print $6}' > S_Mg") # in eV/K
    os.system("grep 'Entropy per atom in k_B'  entropy_Mg.out | awk '{print $6}' > S_Mg_in_kB") # in kB

    out = open('entropy_Si.in','w')
    out.write('dos_Si.dat\n')
    out.write('POSCAR\n')
    temp = commands.getoutput('cat T').strip()
    out.write(temp+'\n')
    out.close()
    os.system('entropy_solid_three_types.x < entropy_Si.in > entropy_Si.out')
            
    os.system("grep 'Entropy per atom in eV/K' entropy_Si.out | awk '{print $6}' > S_Si") # in eV/K
    os.system("grep 'Entropy per atom in k_B'  entropy_Si.out | awk '{print $6}' > S_Si_in_kB") # in kB

    out = open('entropy_O.in','w')
    out.write('dos_O.dat\n')
    out.write('POSCAR\n')
    temp = commands.getoutput('cat T').strip()
    out.write(temp+'\n')
    out.close()
    os.system('entropy_solid_three_types.x < entropy_O.in > entropy_O.out')
            
    os.system("grep 'Entropy per atom in eV/K' entropy_O.out | awk '{print $6}' > S_O") # in eV/K
    os.system("grep 'Entropy per atom in k_B'  entropy_O.out | awk '{print $6}' > S_O_in_kB") # in kB

    os.system("paste S_Mg S_Si S_O | awk '{print $1+$2+$3}' > S")
    os.system("paste S_Mg_in_kB S_Si_in_kB S_O_in_kB | awk '{print $1+$2+$3}' > S_in_kB")

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
    os.system('mv -f *.dat *.xyz *.in *.avg *.out date toten E P V T S S_* PV TS H F G gibbs/')

    os.system('rm -f OUTCAR OSZICAR XDATCAR')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
