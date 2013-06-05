#!/usr/bin/env python
"""
md-CO2-vasp_analyze.py VERSION 4.0
Author: Brian Boates

Run all applicable scripts to analyze a CO2 system
"""
import os, sys, commands, glob, time

global nFinalConfigs, nSkipConfigs
nFinalConfigs = 4000  # Final 3 ps (w/ 0.75 fs tstep)
nSkipConfigs = 7000   # Skip initial 5 ps (w/ 0.75 fs tstep)
#nFinalConfigs = 1000
#nSkipConfigs = 500

def main():
    """
    Execute several vasp md analysis scripts
    """
    if 'output' in glob.glob('*'):

        if 'output/OUTCAR.100' not in glob.glob('output/*'):

            if 'OUTCAR' in glob.glob('*'):
                os.system('mv OUTCAR output/OUTCAR.100')
                os.system('mv XDATCAR output/XDATCAR.100')
                os.system('mv OSZICAR output/OSZICAR.100')

            else:
                print '\nNo output files found ---- exiting...\n'
                sys.exit(0)

        else:
            if 'OUTCAR' in glob.glob('*'):
                os.system('mv OUTCAR output/OUTCAR.199')
                os.system('mv XDATCAR output/XDATCAR.199')
                os.system('mv OSZICAR output/OSZICAR.199')

        print 'Concatenating VASP files...'

        os.system('cp output/XDATCAR.100 ./XDATCAR.100')
        xdats = glob.glob('output/XDATCAR*')
        xdats.pop(xdats.index('output/XDATCAR.100'))
        xdats.sort()
        for xdat in xdats:
            os.system('tail -n+6 '+xdat+' > XDATCAR.'+xdat[-3:])

        # Concatenate the XDATCAR files
        xdats = glob.glob('XDATCAR.*')
        xdats.sort()
        cmd = 'cat '
        for xdat in xdats:
            cmd += xdat+' '
        cmd += '> XDATCAR'
        os.system('nice -15 '+cmd)
        
        # Concatenate the OUTCAR files
        outs = glob.glob('output/OUTCAR*')
        outs.sort()
        cmd = 'cat '
        for out in outs:
            cmd += out+' '
        cmd += '> OUTCAR'
        os.system('nice -15 '+cmd)

        # Remove intermediate data files
        os.system('rm -f XDATCAR.*')

        print 'done...'

    # Check to make sure all the needed vasp files are present
    if ('OUTCAR' and 'INCAR' and 'XDATCAR' and 'POSCAR') not in glob.glob('*'):
        print '\nCould not find file(s): OUTCAR, INCAR, XDATCAR, OSZICAR, POSCAR\n'
	sys.exit(0)

    # Move any previous analysis
    if 'analysis/' in glob.glob('*/'):
        print 'Moving previous analysis... '
        if 'previous_analysis/' in glob.glob('*/'):
            os.system('rm -rf previous_analysis/')
        os.system('mv -f analysis/ previous_analysis')
        print 'done'

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
            
    # Extract thermodynamic quantities
    print 'Extracting T, P, & E and averaging...'
    os.system('nice -15 vasp_pressure.py')
    os.system('tail -5000 pressure.dat > p.dat')
    os.system('blocker p.dat 2 > pressure.blocker')
    os.system('rm -f p.dat')
    os.system('nice -15 vasp_temperature.py')
    os.system('nice -15 vasp_energy.py')
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
    os.system('separate_CO2_xyz.py wrapped_CO2.xyz')
    os.system('mv TRAJEC_C.xyz wrapped_C.xyz')
    os.system('mv TRAJEC_O.xyz wrapped_O.xyz')

    os.system('separate_CO2_xyz.py TRAJEC_CO2.xyz')
    print 'done'

    # Determine the number of steps
    nConfigs = int(commands.getoutput('wc -l TRAJEC_CO2.xyz').split()[0])/(natom+2)

    # Slice the xyz lengths for analysis for CO2
    fnameFinal_CO2 = 'FINAL-'+str(nFinalConfigs)+'_CO2.xyz'
    nFinalLines_CO2 = (natom + 2)*nFinalConfigs
    os.system('tail -'+str(nFinalLines_CO2)+' TRAJEC_CO2.xyz > '+fnameFinal_CO2)
    fnameSkip_CO2 = 'SKIP-'+str(nSkipConfigs)+'_CO2.xyz'
    nSkipLines_CO2 = (nConfigs - nSkipConfigs)*(natom+2)
    os.system('tail -'+str(nSkipLines_CO2)+' TRAJEC_CO2.xyz > '+fnameSkip_CO2)

    # Slice the xyz lengths for analysis for C
    fnameFinal_C = 'FINAL-'+str(nFinalConfigs)+'_C.xyz'
    nFinalLines_C = (nC + 2)*nFinalConfigs
    os.system('tail -'+str(nFinalLines_C)+' TRAJEC_C.xyz > '+fnameFinal_C)
    fnameSkip_C = 'SKIP-'+str(nSkipConfigs)+'_C.xyz'
    nSkipLines_C = (nConfigs - nSkipConfigs)*(nC+2)
    os.system('tail -'+str(nSkipLines_C)+' TRAJEC_C.xyz > '+fnameSkip_C)

    # Slice the xyz lengths for analysis for O
    fnameFinal_O = 'FINAL-'+str(nFinalConfigs)+'_O.xyz'
    nFinalLines_O = (nO + 2)*nFinalConfigs
    os.system('tail -'+str(nFinalLines_O)+' TRAJEC_O.xyz > '+fnameFinal_O)
    fnameSkip_O = 'SKIP-'+str(nSkipConfigs)+'_O.xyz'
    nSkipLines_O = (nConfigs - nSkipConfigs)*(nO+2)
    os.system('tail -'+str(nSkipLines_O)+' TRAJEC_O.xyz > '+fnameSkip_O)

    # Create cnn and nn_average.hist files
    print 'Extracting cnn files...'
    os.system('nice -15 nn_dist.py 200 '+str(a/0.529177)+' '+str(b/0.529177)+' '+str(c/0.529177)+' '+fnameFinal_CO2+' ALL')
    os.system('nice -15 nn_dist_two_species.py 200 '+str(a/0.529177)+' '+str(b/0.529177)+' '+str(c/0.529177)+' '+fnameFinal_CO2+' '+str(nO)+' CO')
    os.system('nice -15 nn_dist.py 200 '+str(a/0.529177)+' '+str(b/0.529177)+' '+str(c/0.529177)+' '+fnameFinal_C+' C')
    os.system('nice -15 nn_dist.py 200 '+str(a/0.529177)+' '+str(b/0.529177)+' '+str(c/0.529177)+' '+fnameFinal_O+' O')
    print 'done'

    # Create TRAJEC.mol file
    print 'Extracting mol file...'
    out = open('mol.in','w')
    out.write('TRAJEC_ALL.cnn\n')
    out.write(str(nC)+'\n')
    out.write(str(nO)+'\n')
    out.close()
    os.system('nice -15 CO2_molecule_detector.x < mol.in > mol.out')
    os.system('rm -f mol.in mol.out')
    print 'done'

    # Create COM xyz file
    print 'Extracting COM xyz and cnn files...'
    out = open('com.in','w')
    out.write('TRAJEC.mol\n')
    out.close()
    os.system('CO2_COM.x < com.in > com.out')
    os.system('rm -f com.in com.out')
    print 'done'

    # Calculate all orientational information
    print 'Calculating only carbon cenetered bond angles...'
    out = open('angles.in','w')
    out.write('TRAJEC.mol\n')
    out.close()
    os.system('CO2_angles.x < angles.in > angles.out')
    os.system('rm -f angles.in angles.out')
    print 'done'

    # Run CO2_track_angles.x to track individual molecule angles over time
    print 'Running CO2_track_angles.x...'
    out = open('track.in','w')
    out.write('TRAJEC.mol\n')
    out.write(str(step_in_fs)+'\n')
    out.close()
    os.system('CO2_track_angles.x < track.in > track.out')
    os.system('rm -f track.in track.out')
    print 'done'

    # Run angluar velocity analysis to find equilibrium CO2 angles
    # CO2_tracked_angles.dat file must be present.
    print 'Running CO2 angular velocity analysis'
    out = open('ang_vel.in','w')
    out.write(str(nC)+'\n')
    out.close()
    os.system('CO2_angular_velocity_analysis.x < ang_vel.in > ang_vel.out')
    os.system('rm -f ang_vel.in ang_vel.out')
    print 'done'

    # Calculate "Angular distribution function" (ADF) w/ CO2_ADF.x
    # CO2_tracked_angles.dat file must be present.
    print 'Running CO2_ADF.x...'
    out = open('adf.in','w')
    out.write(str(nC)+'\n')
    out.close()
    os.system('CO2_ADF.x < adf.in > adf.out')
    os.system('rm -f adf.in adf.out')
    print 'done'

    # Run CO2 orientation analysis 
    print 'Running CO2 orientation analysis...'
    out = open('all.in','w')
    out.write('TRAJEC_ALL.cnn\n')
    out.write('TRAJEC_C.cnn\n')
    out.write('0\n')
    out.write('200\n')
    out.close()
    os.system('CO2_orientation.x < all.in > all.out')
    os.system('rm -f all.in all.out')

    out = open('one.in','w')
    out.write('TRAJEC_ALL.cnn\n')
    out.write('TRAJEC_C.cnn\n')
    out.write('1\n')
    out.write('200\n')
    out.close()
    os.system('CO2_orientation.x < one.in > one.out')
    os.system('rm -f one.in one.out')
    print 'done'

    # Bonding lifetime analysis
    print 'Calculating cluster lifetimes for COx...'
    if nFinalConfigs/2 >= nConfigs:
        var = str(nConfigs - 1)
    else:
        var = str(nFinalConfigs/2)
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC_CO.cnn 1 '+var+' > cluster1.out &')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC_CO.cnn 2 '+var+' > cluster2.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC_CO.cnn 3 '+var+' > cluster3.out &')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC_CO.cnn 4 '+var+' > cluster4.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC_CO.cnn 5 '+var+' > cluster5.out &')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC_CO.cnn 6 '+var+' > cluster6.out')

    # Stall time in seconds (2 minutes)
    time.sleep(120)
    os.system('rm -f cluster*.out')

    # Convert from tstep to picoseconds
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_1nn.dat > CO1.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_2nn.dat > CO2.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_3nn.dat > CO3.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_4nn.dat > CO4.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_5nn.dat > CO5.dat")
    os.system("awk '{print $1*"+str(step_in_fs/1000.0)+",$2}' correlator_6nn.dat > CO6.dat")
    os.system('rm -f correlator_*nn.dat')

    # Grab the correlators at given timesteps
    os.system('grab_correlators_at_timestep_CO2.sh 100')
    os.system('grab_correlators_at_timestep_CO2.sh 1000')
    print 'done'

    # Calculate RDF's
    print 'Calculating RDFs...'
    out = open('RDF_CC.in','w')
    out.write(fnameSkip_CO2+'\n')
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
    out.write(fnameSkip_CO2+'\n')
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
    out.write(fnameSkip_CO2+'\n')
    out.write('RDF_OO.dat\n')
    out.write('O\nO\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_OO.in > RDF_OO.out')

    out = open('RDF_com.in','w')
    out.write('CO2_COM.xyz\n')
    out.write('RDF_com.dat\n')
    out.write('H\nH\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()
    os.system('nice -15 RDF < RDF_com.in > RDF_com.out')

    os.system('rm -f RDF_CC.in RDF_CO.in RDF_OO.in RDF_com.in')
    os.system('rm -f RDF_CC.out RDF_CO.out RDF_OO.out RDF_com.out')
    print 'done'

    # Calculate coordination's
    print 'Calculating coordinations...'
#    os.system('min_gr.py RDF_CC.dat')
#    min_gr = commands.getoutput("tail min_gr.dat | awk '{print $1}'")
#    os.system('rm -f min_gr.dat RDF_CC.dat.avg')
#    os.system('nice -15 get_coord.pl TRAJEC_C.cnn '+min_gr)
#    os.system('av_coor_to_column.py av_coordination_rc_'+min_gr+'.dat')
#    os.system('mv -f av_coordination_rc_'+min_gr+'.dat av_coordination_rc_'+min_gr+'_C.dat')
#    os.system('mv -f coordination_rc_'+min_gr+'.dat coordination_rc_'+min_gr+'_C.dat')
#    os.system('mv -f coordination.dat coordination_C.dat')

    os.system('min_gr.py RDF_CO.dat')
    min_gr = commands.getoutput("tail min_gr.dat | awk '{print $1}'")
    os.system('rm -f min_gr.dat RDF_CO.dat.avg')
    os.system('nice -15 get_coord.pl TRAJEC_CO.cnn '+min_gr)
    os.system('av_coor_to_column.py av_coordination_rc_'+min_gr+'.dat')
    os.system('mv -f av_coordination_rc_'+min_gr+'.dat av_coordination_rc_'+min_gr+'_CO.dat')
    os.system('mv -f coordination_rc_'+min_gr+'.dat coordination_rc_'+min_gr+'_CO.dat')
    os.system('mv -f coordination.dat coordination_CO.dat')

    os.system('min_gr.py RDF_OO.dat')
    min_gr = commands.getoutput("tail min_gr.dat | awk '{print $1}'")
    os.system('rm -f min_gr.dat RDF_OO.dat.avg')
    os.system('nice -15 get_coord.pl TRAJEC_O.cnn '+min_gr)
    os.system('av_coor_to_column.py av_coordination_rc_'+min_gr+'.dat')
    os.system('mv -f av_coordination_rc_'+min_gr+'.dat av_coordination_rc_'+min_gr+'_O.dat')
    os.system('mv -f coordination_rc_'+min_gr+'.dat coordination_rc_'+min_gr+'_O.dat')
    os.system('mv -f coordination.dat coordination_O.dat')
    print 'done'

    # Calculate VACF's
    print 'Calculating VACFs...'
    out = open('VACF.in','w')
    out.write(fnameSkip_CO2+'\n')
    out.write('VACF_ALL.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('ALL\n')
    out.close()
    os.system('nice -15 VACF < VACF.in > VACF.out')
    os.system('rm -f VACF.in VACF.out')
    os.system('mv -f DiffCoef.data DiffCoef_ALL.dat')

    out = open('VACF.in','w')
    out.write(fnameSkip_C+'\n')
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
    out.write(fnameSkip_O+'\n')
    out.write('VACF_O.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('ALL\n')
    out.close()
    os.system('nice -15 VACF < VACF.in > VACF.out')
    os.system('rm -f VACF.in VACF.out')
    os.system('mv -f DiffCoef.data DiffCoef_O.dat')

    out = open('VACF.in','w')
    out.write('CO2_COM.xyz\n')
    out.write('VACF_com.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('ALL\n')
    out.close()
    os.system('nice -15 VACF < VACF.in > VACF.out')
    os.system('rm -f VACF.in VACF.out')
    os.system('mv -f DiffCoef.data DiffCoef_com.dat')
    print 'done'

    # Construct VDOS input file
    print 'Calculating VDOSs...'
    out = open('VDOS.in','w')
    out.write('VACF_ALL.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS.in > VDOS.out')
    os.system('mv -f VDOS.dat VDOS_ALL.dat')
    os.system('vdos_to_wavenumber.sh VDOS_ALL.dat')
    os.system('running_avg.py VDOS_ALL.dat 2 10')
    os.system('rm -f VDOS.in VDOS.out')

    out = open('VDOS.in','w')
    out.write('VACF_C.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS.in > VDOS.out')
    os.system('mv -f VDOS.dat VDOS_C.dat')
    os.system('vdos_to_wavenumber.sh VDOS_C.dat')
    os.system('running_avg.py VDOS_C.dat 2 10')
    os.system('rm -f VDOS.in VDOS.out')

    out = open('VDOS.in','w')
    out.write('VACF_O.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS.in > VDOS.out')
    os.system('mv -f VDOS.dat VDOS_O.dat')
    os.system('vdos_to_wavenumber.sh VDOS_O.dat')
    os.system('running_avg.py VDOS_O.dat 2 10')
    os.system('rm -f VDOS.in VDOS.out')

    out = open('VDOS.in','w')
    out.write('VACF_com.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()
    os.system('VDOS.x < VDOS.in > VDOS.out')
    os.system('mv -f VDOS.dat VDOS_com.dat')
    os.system('vdos_to_wavenumber.sh VDOS_com.dat')
    os.system('running_avg.py VDOS_com.dat 2 10')
    os.system('rm -f VDOS.in VDOS.out')
    print 'done'
    
    # Calculate msd's
    print 'Calculating msds, MSDs, and diffusion...'
    out = open('msd.in','w')
    out.write(fnameSkip_CO2+'\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')
    os.system('mv -f msd.dat msd_ALL.dat')
    os.system('rm -f msd.in msd.out')

    out = open('msd.in','w')
    out.write(fnameSkip_C+'\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')
    os.system('mv -f msd.dat msd_C.dat')
    os.system('rm -f msd.in msd.out')

    out = open('msd.in','w')
    out.write(fnameSkip_O+'\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')
    os.system('mv -f msd.dat msd_O.dat')
    os.system('rm -f msd.in msd.out')

    out = open('msd.in','w')
    out.write('CO2_COM.xyz\n')
    out.close()
    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('rm -f msd.in msd.out')
    os.system('mv -f msd.dat msd_com.dat')
    os.system('rm -f msd.in msd.out')

    # Calculate MSD's & diffusion
    out = open('MSD.in','w')
    out.write(fnameSkip_CO2+'\n')
    out.write('MSD_ALL.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('100\n')
    out.write('5\n')
    out.write(' \n')
    out.close()
    os.system('nice -15 MSD < MSD.in > MSD.out')
    os.system('rm -f MSD.in MSD.out')
    os.system('nice -15 diffusion.py MSD_ALL.dat')
    os.system('mv -f D.dat diffusion_ALL.dat')

    out = open('MSD.in','w')
    out.write(fnameSkip_C+'\n')
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
    out.write(fnameSkip_O+'\n')
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

    out = open('MSD.in','w')
    out.write('CO2_COM.xyz\n')
    out.write('MSD_com.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('100\n')
    out.write('5\n')
    out.write(' \n')
    out.close()
    os.system('nice -15 MSD < MSD.in > MSD.out')
    os.system('rm -f MSD.in MSD.out')
    os.system('nice -15 diffusion.py MSD_com.dat')
    os.system('mv -f D.dat diffusion_com.dat')

    os.system('rm -f fort.* diss.dat')
    print 'done'

    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')
    os.system('mkdir analysis')
    os.system('mv -f nohup.out analysis/log')
    os.system('mv -f *.dat *.xyz *.cnn *.mol *.in *.blocker *.hist *.avg date analysis/')

    os.system('rm -f OUTCAR XDATCAR')

    # Copy any relevant plotting scripts to analysis/
    os.system('cp /home/boates/software/CO2/gnuplot_scripts_CO2/* analysis/')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
