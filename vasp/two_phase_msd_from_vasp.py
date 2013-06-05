#!/usr/bin/env python
"""
Separate two phase trajectories respectively and calculate separate
msd curves for each 'half'
"""
import os, sys, commands

def main():

    # Retrieve user input
    try:
        f = open('XDATCAR','r')
        f.close()
        f = open('POSCAR','r')
        f.readline()
        alat = f.readline().strip()
        f.close()
    except:
        print '\n Make sure both XDATCAR and POSCAR are present. \n'
        sys.exit(0)

    # Extract xyz file from XDATCAR
    print 'Extracting xyz file...'
    os.system('xyz_from_vasp.x')
    
    # Unwrap the xyz file
    print 'Unwrapping xyz file...'
    out = open('unwrap.in','w')
    out.write('TRAJEC.xyz\n')
    out.write(alat+','+alat+','+alat+'\n')
    out.close()
    os.system('unwrap_PBC.x < unwrap.in > unwrap.out')
    os.system('rm -f unwrap.*')

    # Relabel each of the two phases
    print 'Renaming atoms in xyz file...'
    out = open('rename.in','w')
    out.write('unwrapped.xyz\n')
    out.close()
    os.system('two_phase_rename_half_xyz.x < rename.in > rename.out')
    os.system('rm -f rename.*')

    # Change the labelling of the two phases to resemble CO2 naming
    os.system("sed s/N/O/g TWO_PHASE.xyz > RENAMED.xyz")
    os.system('mv RENAMED.xyz TWO_PHASE.xyz')

    # Separate into separate xyz files
    print 'Separating into separate species xyz files...'
    os.system('separate_CO2_xyz.py TWO_PHASE.xyz')

    # Calculate individual species MSD curves
    print 'Calculating MSD...'
    natom = float(commands.getoutput('head -1 TRAJEC.xyz').strip())
    nConfig = str( int( int(commands.getoutput('wc -l TRAJEC.xyz').split()[0]) / (natom+2) ) )
    out = open('msd1.in','w')
    out.write('FINAL-'+nConfig+'_C.xyz\n')
    out.close()
    os.system('msd.x < msd1.in > msd1.out')
    os.system('mv msd.dat msd_cgN.dat')
    out = open('msd2.in','w')
    out.write('FINAL-'+nConfig+'_O.xyz\n')
    out.close()
    os.system('msd.x < msd2.in > msd2.out')
    os.system('mv msd.dat msd_polyliquid.dat')
    os.system('rm -f msd1.* msd2.*')


if __name__ == '__main__':
    main()
