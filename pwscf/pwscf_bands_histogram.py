#!/usr/bin/env python
"""
Create a histogram of band energies, and shift them back with
respect to each respective fermi energy value.
"""
import os, sys, commands, glob
import numpy, matplotlib, pylab

def main():
    """
    Execute from simulation directory (the one with an output/ directory in it)
    Creates a .png file of a histogram of the band energies shifted with respect
    to each respective calculated Fermi energy.
    """
    if len( glob.glob('input.*.pw') ) == 1:
        f_in = glob.glob('input.*.pw')[0]
        print '\nNBANDS & NATOMS found in the input file: '+f_in+'\n'
    else:
        # Ask the user for the name of their pwscf input file if not given as sys.argv[1]
        f_in = input('\nWhat is the name of your pwscf input file? '+\
                                   '(please enter in quotations: \'input.dat\')\n')

    # Get the number of bands & atoms from the input file
    NBANDS = int(commands.getoutput("grep nbnd "+f_in+" | awk '{print $3}'"))
    NATOMS = int(commands.getoutput("grep nat "+f_in+" | awk '{print $3}'"))

    # PWSCF output files give 8 band energies per line of output
    if NBANDS % 8 == 0:
        NLINES = int( NBANDS / 8 )
    else:
        NLINES = int( NBANDS / 8 ) + 1

    # Retrieve a list of all pwscf output files for the simulation
    # NOTE: This assumes you run bands_histogram.py from the simulation directory
    #       with an output/ directory
    output_files = commands.getoutput('ls -1 output/*.out').split() + glob.glob('*.out')

    # Get the number of valence electrons from one of the output files
    try: 
        NVALENCE = int(float(commands.getoutput("head -1000 "+output_files[0]+" | grep -A1 valence "+
                                                "| tail -1 | awk '{print $2}'")))
    except ValueError:
        NVALENCE = input('\nCould not find number of valence electrons from output file.\n'+
                         'Please give the number of valence electrons:   ')

    # Indicate to the user how many bands, atoms and valence electrons were found
    print 'Found NBANDS = '+str(NBANDS)+', NATOMS = '+str(NATOMS)+', & NVALENCE = '+str(NVALENCE)

    # Read in the band energies (in eV) as a 1D array (the +1 is for additional non-band lines)
    BAND_array = []
    HOMO_array = []
    LUMO_array = []
    for out in output_files:
        BAND_array += commands.getoutput("grep -A"+str(NLINES+1)+" 'bands (ev)' "+out+\
                                                                 " | grep -v 'bands (ev)'").split()
        HOMO_array += commands.getoutput("grep 'highest occupied' "+out+" | awk '{print $7}'").split() # in eV
        LUMO_array += commands.getoutput("grep 'highest occupied' "+out+" | awk '{print $8}'").split() # in eV

    # Remove grep artifacts
    for i in range(BAND_array.count('--')):
        BAND_array.remove('--')
    for i in range(HOMO_array.count('--')):
        HOMO_array.remove('--')
    for i in range(LUMO_array.count('--')):
        LUMO_array.remove('--')

#    print ( len(BAND_array) / NBANDS ), len(HOMO_array), len(LUMO_array)

    # Perform a consistency check (between BAND_array and FERMI_array)
    CONSISTENT = ( len(BAND_array) / NBANDS ) == len(HOMO_array) == len(LUMO_array)
    if CONSISTENT:
        print '\nThe number of band energies is consistent with the number of '+\
              'bands and HOMO/LUMO energies.\n'
    if not CONSISTENT:
        print '\n******** CRASH ********'
        print 'The number of band energies is inconsistent with the number of '+\
              'bands and HOMO/LUMO energies... exiting\n'
        sys.exit(0)

    # Subtract off the respective fermi energy for each set of band energies
    for i in range( len(BAND_array) ):
        BAND_array[i] = float(BAND_array[i]) - float(HOMO_array[i/NBANDS])
    for i in range( len(LUMO_array) ):
        LUMO_array[i] = float(LUMO_array[i]) - float(HOMO_array[i/NBANDS])
    AVG_LUMO = numpy.average(LUMO_array)

    # Check to see how many states are filled below the Fermi level
    NELEC = NATOMS*NVALENCE
    NFILLED = NELEC / 2.0
    NBELOW = []
#    for i in ( range(0,len(BAND_array),NBANDS) ):
#        count_below = 0
#        for j in range(NBANDS):
#            if BAND_array[j+i] < 0:
#                count_below += 1
#        NBELOW.append(count_below)

#    # Print relevant filled bands information
#    MIN_BELOW = min(NBELOW)
#    MAX_BELOW = max(NBELOW)
#    AVG_BELOW = numpy.average(NBELOW)
#    print 'The number of filled bands in the system is:',NFILLED
#    print 'The minimum number of bands below the HOMO level was:',MIN_BELOW
#    print 'The maximum number of bands below the HOMO level was:',MAX_BELOW
#    print 'The average number of bands below the HOMO level is:',AVG_BELOW,'\n'

    # Create the figure window to a reasonable size and format the plot
    pylab.figure(num=1,figsize=(6,4),facecolor='w',edgecolor='k')
    pylab.title('Histogram of Band Energies',fontsize=10)
    pylab.xlabel('Band energy shifted wrt HOMO (eV)',fontsize=8)
    pylab.ylabel('Normalized to one',fontsize=8)
    pylab.xticks(fontsize=6)
    pylab.yticks(fontsize=6)

    # Create the histogram with 250 bins, normalized to 1
    pylab.hist(BAND_array, 250, 1, label='_nolegend_')
    
    # Get the original axis limits
    axis = pylab.axis()

    # Plot a vertical line at the HOMO energies
    pylab.plot([0.0,0.0],[0.0,1.0],'r--',linewidth=1,label='HOMO')

    # Plot a vertical line at the AVERAGE LUMO energy
    pylab.plot([AVG_LUMO,AVG_LUMO],[0.0,1.0],'g--',linewidth=1,label='AVG_LUMO')

    # Create the legend
    pylab.legend()

    # Reset the original axis limits
    pylab.axis(axis)

    # Save the histogram to file
    pylab.savefig('bands_histogram.png')

if __name__ == '__main__':
    main()