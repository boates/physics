#!/usr/bin/env python
"""
md-hugoniot-analyze.py
Author: Brian Boates

Run analysis on a monatomic system to get
TD variables to compute Hugoniot
"""
import os, sys, commands, glob, time

global tFinal
tFinal = 3  # in picoseconds

def main():
    """
    Execute several vasp md analysis scripts
    """
    print 'Concatenating VASP files...'

    os.system('vasp_cleanup')

    cwd = os.getcwd()
    os.chdir('output/')        
    # Concatenate the OUTCAR files
    os.system('nice -15 cat OUTCAR.* > OUTCAR')
    os.system('mv -v OUTCAR ../OUTCAR')
    # Concatenate the OSZICAR files
    os.system('nice -15 cat OSZICAR.* > OSZICAR')
    os.system('mv -v OSZICAR ../OSZICAR')
    os.chdir(cwd)
    print 'done...'

    # Check to make sure all the needed vasp files are present
    if ('OUTCAR' and 'OSZICAR' and 'POSCAR' and 'INCAR') not in glob.glob('*'):
        print '\nCould not find one or more of: OUTCAR, OSZICAR, POSCAR, INCAR\n'
	sys.exit(0)

    # Retrieve natom and tstep
    natom = int(commands.getoutput("head -6 POSCAR | tail -1"))
    step_in_fs = float(commands.getoutput("grep POTIM INCAR | awk '{print $3}'"))
    
    # Determine nFinalConfigs
    nFinalConfigs = int( round( tFinal / (step_in_fs/1000.) ) )
            
    # Extract thermodynamic quantities
    os.system('rm -f V T P E VTPE') # remove old
    
    print 'Extracting thermodynamic variables and averaging...'

    # VOLUME #
    os.system("volume_from_POSCAR.py POSCAR | awk '{print $2}' > V")

    # TEMPERATURE #
    os.system("grep TEBEG INCAR | awk '{print $3}' > T")

    # PRESSURE #
    os.system('nice -15 thermodynamic_vasp.py OSZICAR OUTCAR')
    os.system('tail -'+str(nFinalConfigs)+' pressure.dat > P.dat')
    os.system("blocker P.dat | grep Final | awk '{print $4}' > P")

    # ENERGY #
    os.system('nice -15 vasp_hugoniot_energy.py')
    os.system('mv energy_hugoniot.dat Ehug.dat')
    os.system('tail -'+str(nFinalConfigs)+' Ehug.dat > E.dat')
    os.system("blocker E.dat | grep Final | awk '{print $4/"+str(natom)+"}' > E")

    # ALL #
    os.system('paste V T P E > VTPE')

    os.system('rm -f E.dat P.dat pressure.dat temperature.dat Ehug.dat toten.dat energy.dat')
    print 'done'

    os.system('rm -f OUTCAR OSZICAR')

    print 'Analysis complete'

if __name__ == '__main__':
    main()
