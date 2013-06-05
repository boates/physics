#!/usr/bin/python
"""
md-CO2-dlpoly-analyze.py VERSION 1.0
Author: Brian Boates

Run analysis on a CO2 system run in DL_POLY
"""
import os, sys, commands, glob, time

global tFinal
tFinal = 5  # in picoseconds

def main():
    """
    Execute several vasp md analysis scripts
    """
    # Check to make sure all the needed vasp files are present
    if ('OUTPUT' and 'STATIS' and 'HISTORY' and 'CONFIG' and 'CONTROL') not in glob.glob('*'):
        print '\nCould not find one or more of: OUTPUT, STATIS, HISTORY, CONFIG, CONTROL'
	sys.exit(0)

    # Remove any previous analysis
    if 'analysis/' in glob.glob('*/'):
        os.system('rm -rf analysis/')

    # Retrieve timestep (fs), natom, & lattice constant
    step_in_fs = float(commands.getoutput("grep timestep CONTROL | awk '{print $2}'"))*1000.
    step_in_au = step_in_fs / 2.41880e-02
    nC = int(commands.getoutput("grep 'number of molecules' OUTPUT | head -1").split()[-1])
    nO = int(commands.getoutput("grep 'number of molecules' OUTPUT | tail -n-1").split()[-1])
    natom = nC + nO
    a = float(commands.getoutput("head -3 CONFIG | tail -n-1").split()[0])  #
    b = float(commands.getoutput("head -4 CONFIG | tail -n-1").split()[1])  # In Angstroms
    c = float(commands.getoutput("head -5 CONFIG | tail -n-1").split()[2])  #

    # Determine nFinalConfigs
    nFinalConfigs = int( round( tFinal / (step_in_fs/1000.) ) )
            
    # Extract thermodynamic quantities
    print 'Extracting thermodynamic variables and averaging...'
    os.system('run_STATIS_twotypes.sh')
    os.system("awk '{print $2,$3}' statis.dat > pressure.dat")
    os.system("awk '{print $2,$5}' statis.dat > temperature.dat")
    os.system("awk '{print $2,$6}' statis.dat > energy.dat")
    
    os.system('tail -'+str(nFinalConfigs)+' pressure.dat    > P.dat')
    os.system('tail -'+str(nFinalConfigs)+' temperature.dat > T.dat')
    os.system('tail -'+str(nFinalConfigs)+' energy.dat      > E.dat')

    os.system('blocker P.dat > pressure.blocker')
    os.system('blocker T.dat > temperature.blocker')
    os.system('blocker E.dat > energy.blocker')

    os.system('rm -f P.dat T.dat E.dat')
    print 'done'

    # Extract the xyz file
    print 'Extracting xyz file... '
    os.system('nice -15 convert_HISTORY_xyz.x')
    print 'done'

    # Write an input file for unwrap_PBC.x
    print 'Unwrapping xyz file...'
    out = open('unwrap.in','w')
    out.write('HISTORY.xyz\n')
    out.write(str(a)+','+str(b)+','+str(c)+'\n')
    out.close()
    os.system('nice -15 unwrap_PBC.x < unwrap.in > unwrap.out')
    os.system('rm -f unwrap.out')
    os.system('mv HISTORY.xyz wrapped.xyz')
    os.system('mv unwrapped.xyz TRAJEC.xyz')
    print 'done'

    # Slice the xyz lengths for analysis for CO2
    fnameFinal = 'FINAL-'+str(tFinal)+'ps.xyz'
    nFinalLines = (natom + 2)*nFinalConfigs
    os.system('tail -'+str(nFinalLines)+' TRAJEC.xyz > '+fnameFinal)

    # Calculate RDF's
    print 'Calculating RDFs...'
    out = open('RDF_CC.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF_CC.dat\n')
    out.write('C\nC\n')
    out.write(str(min([a,b,c])/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_CC.in > RDF_CC.out')

    out = open('RDF_CO.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF_CO.dat\n')
    out.write('C\nO\n')
    out.write(str(min([a,b,c])/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_CO.in > RDF_CO.out')

    out = open('RDF_OO.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF_OO.dat\n')
    out.write('O\nO\n')
    out.write(str(min([a,b,c])/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_OO.in > RDF_OO.out')

    os.system('rm -f RDF_*.in RDF_*.out')
    print 'done'

    # Calculate VACF's
    print 'Calculating VACFs...'
    out = open('VACF_C.in','w')
    out.write(fnameFinal+'\n')
    out.write('VACF_C.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('C\n')
    out.close()
    os.system('nice -15 VACF < VACF_C.in > VACF_C.out')
    os.system('mv -f DiffCoef.data DiffCoef_C.dat')

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

    os.system('rm -f VACF_*.in VACF_*.out')
    print 'done'

    # Construct VDOS input file
    print 'Calculating VDOSs...'
    out = open('VDOS_C.in','w')
    out.write('VACF_C.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS_C.in > VDOS_C.out')
    os.system('vdos_to_wavenumber.sh VDOS.dat')
    os.system('mv -f VDOS.dat VDOS_C.dat')
    os.system('running_avg.py VDOS_C.dat 2 10')

    out = open('VDOS_O.in','w')
    out.write('VACF_O.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS_O.in > VDOS_O.out')
    os.system('vdos_to_wavenumber.sh VDOS.dat')
    os.system('mv -f VDOS.dat VDOS_O.dat')
    os.system('running_avg.py VDOS_O.dat 2 10')

    xC = float(nC) / float(natom)
    xO = float( 1.0 - xC )
    os.system("paste VDOS_C.dat VDOS_O.dat | awk '{print $1,"+str(xC)+"*$2+"+str(xO)+"*$4}' > VDOS_CO2.dat")
    os.system('running_avg.py VDOS_CO2.dat 2 10')

    os.system('rm -f VDOS_*.in VDOS_*.out')
    print 'done'

    # Calculate msd's
    print 'Calculating msd...'
    out = open('msd.in','w')
    out.write('TRAJEC.xyz\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')
    os.system('rm -f msd.in msd.out')
    print 'done'

    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')
    os.system('mkdir analysis')
    os.system('mv -f nohup.out analysis/log')
    os.system('mv -f *.dat *.xyz *.in *.blocker *.avg date analysis/')

    # Copy any relevant plotting scripts to analysis/
#    os.system('cp /home/boates/software/CO2/gnuplot_scripts_CO2/* analysis/')
    print 'Analysis complete'


if __name__ == '__main__':
    main()
