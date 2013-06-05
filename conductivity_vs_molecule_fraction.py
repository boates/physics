#!/usr/bin/env python
"""
Correlate the percentage of molecules with the w=0 conductivity
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        fcond =  sys.argv[1]
        fxyz  =  sys.argv[2]
        alat  =  str( float(sys.argv[3]) / 0.5291772 )   # Convert to Bohr
    except:
        print '\n usage: '+sys.argv[0]+'  sigXX.dat  snapshot.xyz  alat(angstroms)\n'
        sys.exit(0)

    # Get the w=0 conductivity from the file
    column = '4'
    if 'avg' in fcond:
        column = '2'
    cond = commands.getoutput("head -2 "+fcond+" | grep -v '\#' | awk '{print $"+column+"}'")
    cond = str( float(cond) )

    # Create the TRAJEC.cbn file using cbn_from_xyz.x
    inp = open('cbn.in','w')
    inp.write(fxyz+'\n')
    inp.write(alat+','+alat+','+alat+'\n')
    inp.close()
    os.system('cbn_from_xyz.x < cbn.in > cbn.out')
    os.system('rm -f cbn.in cbn.out')

    # Calculate the fraction of molecules using dissociation_from_cbn.x
    inp = open('diss.in','w')
    inp.write('TRAJEC.cbn \n')
    inp.close()
    os.system('dissociation_from_cbn.x < diss.in > diss.out')
    fmol = commands.getoutput("tail -n-1 dissociation.dat | awk '{print $3}'")
    fmol = str( float(fmol) )
    os.system('rm -f TRAJEC.cbn diss.in diss.out dissociation.dat')

    # Print the results to the user
    print ' fraction_of_molecules =\t', fmol
    print ' conductivity_in_(ohm*cm)^-1 =\t', cond


if __name__ == '__main__':
    main()
