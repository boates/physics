#!/usr/bin/env python
"""
md-pwscf_analyze.py VERSION 6.0
Author: Brian Boates
ANALYZE MD DATA FROM PWSCF
  ---> Extract thermodynamic variables using pw_temperature.py, pw_pressure.py, & pw_energy.py
  ---> Calculate the averages of the above thermodynamic variables using blocker
  ---> Extract xyz file using xyz_from_pwscf.py to TRAJEC.xyz
  ---> Create TRAJEC.cbn file form xyz file using cbn_from_xyz.x
  ---> Create TRAJEC.cnn file and nn_average.hist usinng Isaac's nn_distance_distributions_cnn.py
  ---> Calculate angle distributions for neighbors (1,3), (1,2), & (2,3) with nn_angle_distributions.py
  ---> Calculate cluster lifetimes for clusters of 2,3,4,5,6 using bonding_lifetime_n_body.py
  ---> (removed) Use bond_freq.py to create histograms of bond lengths, periods, & frequencies
  ---> Calculate coordination of atoms using Bonev's get_coord.pl
  ---> Calculate VACF using VACF
  ---> Calculate RDF using RDF
  ---> Calculate MSD using MSD and my msd.x code
  ---> Calculate the power spectral density of the VACF
  ---> Calculate Diffusion from MSD.dat and msd.dat with diffusion.py
  ---> Move all created analysis files to analysis/ directory
"""
import os, sys, commands, glob

global nFinalConfigs
nFinalConfigs = 1000

def TRAJEC_trim(fname,nskip=0):
        """
	Trim an xyz file (remove incomplete final configuration if present)
        nskip: Number of initial configurations to skip
	"""
        # Get the number of lines from the xyz file
        nlines = int(commands.getoutput('wc -l '+fname).split()[0])

        # Get the number of atoms in the system
        f = open(fname,'r')
        natom = int( f.readline().strip() )

        # Number of initial lines to skip
        nskiplines = nskip * (natom+2)

        # Get the number of lines in the incomplete configuration
        nleftover = int( round( ( float(nlines)/float(natom+2) - nlines/(natom+2) )*(natom+2) ) )

        # Calculate the number of lines for the trimmed file
        N = nlines - nleftover

        # Replace the original with the trimmed xyz file
        os.system('head -n+'+str(N)+' '+fname+' | tail -'+str(N-nskiplines)+' > '+fname+'.trimmed')
        os.system('mv '+fname+'.trimmed'+' '+fname)

def gr_min(rs):
    """
    Return the min of g(r) (in bohr) for a given rs
    """
    gr_mins = {'1.10':'2.87','1.12':'3.02','1.13':'3.00','1.14':'3.06','1.15':'3.08','1.16':'3.12','1.17':'3.12',\
	       '1.18':'3.14','1.19':'3.16','1.20':'3.16','1.21':'3.16','1.22':'3.17','1.23':'3.17','1.24':'3.17',\
	       '1.25':'3.12','1.26':'3.02','1.27':'3.10','1.28':'3.06','1.29':'2.97','1.30':'2.93','1.31':'2.83',\
	       '1.32':'2.93','1.35':'2.82','1.36':'2.83','1.37':'2.91','1.38':'2.91','1.39':'2.83','1.40':'2.72',\
	       '1.41':'2.83','1.43':'2.89','1.45':'2.97','1.46':'2.97','1.47':'2.93','1.48':'2.97','1.49':'2.97',\
	       '1.50':'2.83','1.55':'2.83','1.60':'3.21','1.70':'3.14','1.80':'3.21'}

    return gr_mins[str(rs)]


def main():
    """
    Execute several pwscf md analysis scripts
    """
    # Efforts to obtain name of pwscf input file for given md run
    try:
        pw_input = sys.argv[1]
    except IndexError:
        try:
            pw_input = glob.glob('input.*.pw')[0]
        except:
            print '\nMake sure you are in the proper directory, exiting now...\n'
	    sys.exit(0)

    # Move any previous analysis
    if 'analysis/' in glob.glob('*/'):
        print 'Moving previous analysis... '
        if 'previous_analysis/' in glob.glob('*/'):
            os.system('rm -rf previous_analysis/')
        os.system('mv -f analysis/ previous_analysis')
        print 'done'
    os.system('rm -f output/fake.out output/output.out')

    # Put output.out into the output directory
    try:
        f = open('output.out','r')
        f.close()
        os.system('mv -f output.out output/')
    except:
        pass

    # Retrieve timestep (au) and lattice constants from input file
    rs = os.getcwd().split('/')[-1]
    step_in_au = int(commands.getoutput("grep dt "+pw_input+" | awk '{print $3}'")) * 2
    celldm = commands.getoutput("grep celldm "+pw_input+" | awk '{print $3}'").split()
    natom = int(commands.getoutput("grep nat "+pw_input+" | awk '{print $3}'"))
    a = str( float(celldm[0]) * 0.5291772 )
    b, c = a, a

    # Extract band energy information (w/ or w/o smearing)
    print 'Analyzing bands...'
    if len(commands.getoutput("grep smearing "+pw_input).split()) > 0:
        os.system('fermi_bands_histogram.py > bands.out')
    else:
        os.system('bands_histogram.py > bands.out')    
    print 'done'
    
#    # Retrieve run information (w/ or w/o smearing)
#    print 'Running pw_monitor.py... '
#    if len(commands.getoutput("grep smearing "+pw_input).split()) > 0:
#        os.system('fermi_pw_monitor.py')
#    else:
#        os.system('pw_monitor.py')
#    print 'done'

    # Extract thermodynamic quantities
    print 'Extracting T, P, & E and averaging...'
    os.system('pw_temperature.py')
    os.system('tail -'+str(nFinalConfigs)+' temperature.dat > T.dat')
    os.system('blocker T.dat 3 > temperature.blocker')
    os.system('pw_pressure.py')
    os.system('tail -'+str(nFinalConfigs)+' pressure.dat > P.dat')
    os.system('blocker P.dat 3 > pressure.blocker')
    os.system('pw_energy.py')
    os.system('tail -'+str(nFinalConfigs)+' energy.dat > E.dat')
    os.system('blocker E.dat 3 > energy.blocker')
    os.system('rm -f T.dat P.dat E.dat')
    print 'done'

    # Create the xyz file from the pwscf output files
    print 'Extracting xyz file... '
    os.system('xyz_from_pwscf.py '+pw_input)

    # Remove trailing incomplete configuration
    TRAJEC_trim('TRAJEC.xyz')
    print 'done'

    # Take the last nFinalConfigs configurations for further analysis
    fnameFinal = 'FINAL-'+str(nFinalConfigs)+'.xyz'
    nFinalLines = (natom + 2)*nFinalConfigs
    os.system('tail -'+str(nFinalLines)+' TRAJEC.xyz > '+fnameFinal)

    # Write an input file for cbn_from_xyz.x
    out = open('cbn.in','w')
    out.write(fnameFinal+'\n')
    out.write(celldm[0]+','+celldm[0]+','+celldm[0]+'\n')
    out.close()

    # Create the cbn file
    print 'Extracting cbn file...'
    os.system('cbn_from_xyz.x < cbn.in > cbn.out')
    os.system('rm -f cbn.out')
    print 'done'

    # Create TRAJEC.cnn and nn_average.hist files
    print 'Extracting cnn file and calculating nn_average.hist...'
    Nbins = 200
    os.system('nn_distance_distributions_cnn.py '+str(Nbins)+' '+celldm[0]+' '+celldm[0]+' ' \
	                                                        +celldm[0]+' '+fnameFinal+' > nn.out')
    os.system('rm -f nn.out')
    print 'done'

    # Calculate angle distributions out to 3rd neighbors
    print 'Calculating angle distributions...'
    os.system('nn_angle_distributions.py 100 1 2 TRAJEC.cnn > angle1.out &')
    os.system('nn_angle_distributions.py 100 1 3 TRAJEC.cnn > angle2.out &')
    os.system('nn_angle_distributions.py 100 2 3 TRAJEC.cnn > angle3.out')
    os.system('rm -f angle*.out')
    print 'done'

    # Calculate lifetimes of clusters of 2, 3, & 4 atoms
    print 'Calculating cluster lifetimes...'
    os.system('bonding_lifetime_n_body.py TRAJEC.cnn 1 '+str(nFinalConfigs/2)+' > cluster1.out &')
    os.system('bonding_lifetime_n_body.py TRAJEC.cnn 2 '+str(nFinalConfigs/2)+' > cluster2.out &')
    os.system('bonding_lifetime_n_body.py TRAJEC.cnn 3 '+str(nFinalConfigs/2)+' > cluster3.out &')
    os.system('bonding_lifetime_n_body.py TRAJEC.cnn 4 '+str(nFinalConfigs/2)+' > cluster4.out &')
    os.system('bonding_lifetime_n_body.py TRAJEC.cnn 5 '+str(nFinalConfigs/2)+' > cluster5.out')
    os.system('rm -f cluster*.out')
    print 'done'
    
    # Calculate the coordination of the atoms
    print 'Calculating coordination of atoms...'
    os.system('get_coord.pl TRAJEC.cnn '+gr_min(rs))
    os.system('av_coor_to_column.py av_coordination_rc_'+gr_min(rs)+'.dat')
    print 'done'

#    # Calculate bond frequencies/periods, and lengths from cbn file
#    print 'Performing bond analysis...'
#    os.system('bond_freq.py TRAJEC.cbn '+str(step_in_au))
#    print 'done'

    # Construct RDF input file
    print 'Calculating RDF, MSD, & VACF... '
    out = open('RDF.in','w')
    out.write(fnameFinal+'\n')
    out.write('RDF.dat\n')
    out.write('N\nN\n')
    out.write(str(float(a)/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()

    # Construct MSD input files
    out = open('MSD.in','w')
    out.write('TRAJEC.xyz\n')
    out.write('MSD.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('1\n')
    out.write('N\n')
    out.close()
    out = open('msd.in','w')
    out.write('TRAJEC.xyz\n')
    out.close()

    # Construct VACF input file
    out = open('VACF.in','w')
    out.write(fnameFinal+'\n')
    out.write('VACF.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('1\n')
    out.write('ALL\n')
    out.close()

    # Calculate RDF, MSD & VACF
    os.system('RDF < RDF.in > RDF.out')
    os.system('MSD < MSD.in > MSD.out')
    os.system('msd.x < msd.in > msd.out')
    os.system('VACF < VACF.in > VACF.out')
    os.system('rm -f RDF.out MSD.out msd.out msd.in VACF.out fort.*')
    print 'done'

    # Calculate PSD of VACF
    print 'Taking PSD of VACF...'
    os.system('PSD.py VACF.dat')
    os.system('mv PSD.dat VACF.psd')
    print 'done'

    # Calculate diffusion
    print 'Calculating diffusion...'
    os.system('diffusion.py MSD.dat DIFFUSION.dat')
    os.system('diffusion.py msd.dat diffusion.dat')
    print 'done'

    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')
    os.system('rm -f temp.xyz')
    os.system('mkdir analysis')
    if 'tmp/pwscf.msd.dat' in glob.glob('tmp/pwscf.msd.dat'):
        os.system('mv -f tmp/pwscf.msd.dat analysis/')
    os.system('mv -f RDF.in RDF.dat MSD.in MSD.dat msd.dat VACF.in VACF.dat VACF.psd DiffCoef.data \
                     DIFFUSION.dat diffusion.dat *.blocker temperature.dat pressure.dat energy.dat \
		     bands* *.hist cbn.in diss.dat correlator_*nn.dat TRAJEC.xyz TRAJEC.cbn TRAJEC.cnn \
	             av_coordination_rc* coordination_rc_*.dat date '+fnameFinal+' analysis/')
    # Copy any relevant plotting scripts to analysis/
    os.system('cp /home/boates/data/nitrogen/pwscf/remote/gnuplot_scripts/* analysis/')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
