#!/usr/bin/env python
import os, sys, glob, commands, time

def main():

    # Retrieve timestep (au), natom, & lattice constant
    step_in_fs = float(commands.getoutput("grep POTIM INCAR | awk '{print $3}'"))
    step_in_au = step_in_fs / 2.41880e-02
    natom = int(commands.getoutput("head -6 POSCAR | tail -1"))
    alat = float(commands.getoutput("head -2 POSCAR | tail -1"))  ### In angstroms
    a = float(commands.getoutput("head -3 POSCAR | tail -1").split()[0])*alat
    b = float(commands.getoutput("head -4 POSCAR | tail -1").split()[1])*alat
    c = float(commands.getoutput("head -5 POSCAR | tail -1").split()[2])*alat
    nConfigs = int(commands.getoutput('wc -l unwrapped.xyz').split()[0])/(natom+2)

    # Average thermodynamic variables
    os.system("awk '{print $2}' pressure.dat > p.dat")
    os.system("awk '{print $2}' temperature.dat > t.dat")
    os.system("mv -f p.dat pressure.dat")
    os.system("mv -f t.dat temperature.dat")
    os.system("blocker pressure.dat > pressure.blocker")
    os.system("blocker temperature.dat > temperature.blocker")
    os.system("blocker energy.dat > energy.blocker")
    
    # Name the xyz files
    fnameSkip = 'unwrapped.xyz'
    fnameFinal = 'unwrapped.xyz'

    # Write an input file for cbn_from_xyz.x
    print 'Extracting cbn file...'
    out = open('cbn.in','w')
    out.write(fnameSkip+'\n')
    out.write(str(a/0.529177)+','+str(b/0.529177)+','+str(c/0.529177)+'\n')
    out.close()

    # Create the cbn file
    os.system('nice -15 cbn_from_xyz.x < cbn.in > cbn.out')
    os.system('rm -f cbn.out')
    print 'done'

    # Create TRAJEC.cnn and nn_average.hist files
    print 'Extracting cnn file and calculating nn_average.hist...'
    os.system('nice -15 nn_dist.py 200 '+str(a/0.529177)+' '+str(a/0.529177)+' '+str(a/0.529177)+' '+fnameFinal)
    print 'done'

    # Calculate angle distributions for nearest neighbors
    print 'Calculating angle distributions...'
    os.system('nice -15 nn_angles.py 100 1 2 TRAJEC.cnn > angles.out')
    os.system('rm -f angles.out')
    print 'done'

    # Calculate lifetimes of clusters of 2-16 atoms
    print 'Calculating cluster lifetimes...'
    var = str(nConfigs/2)
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 1 '+var+' > cluster01.out &')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 2 '+var+' > cluster02.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 3 '+var+' > cluster03.out &')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 4 '+var+' > cluster04.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 5 '+var+' > cluster05.out &')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 6 '+var+' > cluster06.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 7 '+var+' > cluster07.out &')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 8 '+var+' > cluster08.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 9 '+var+' > cluster09.out &')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 10 '+var+' > cluster10.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 11 '+var+' > cluster11.out &')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 12 '+var+' > cluster12.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 13 '+var+' > cluster13.out &')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 14 '+var+' > cluster14.out')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 15 '+var+' > cluster15.out &')
    os.system('nice -15 bonding_lifetime_n_body.py TRAJEC.cnn 16 '+var+' > cluster16.out')

    # Stall time in seconds (2 minutes)
    time.sleep(120)

    os.system('rm -f cluster*.out')
    os.system('mv -f correlator_1nn.dat c01.dat')
    os.system('mv -f correlator_2nn.dat c02.dat')
    os.system('mv -f correlator_3nn.dat c03.dat')
    os.system('mv -f correlator_4nn.dat c04.dat')
    os.system('mv -f correlator_5nn.dat c05.dat')
    os.system('mv -f correlator_6nn.dat c06.dat')
    os.system('mv -f correlator_7nn.dat c07.dat')
    os.system('mv -f correlator_8nn.dat c08.dat')
    os.system('mv -f correlator_9nn.dat c09.dat')
    os.system('mv -f correlator_10nn.dat c10.dat')
    os.system('mv -f correlator_11nn.dat c11.dat')
    os.system('mv -f correlator_12nn.dat c12.dat')
    os.system('mv -f correlator_13nn.dat c13.dat')
    os.system('mv -f correlator_14nn.dat c14.dat')
    os.system('mv -f correlator_15nn.dat c15.dat')
    os.system('mv -f correlator_16nn.dat c16.dat')
    print 'done'

    # Grab the correlators at given timesteps
    os.system('grab_correlators_at_timestep.sh 10')
    os.system('grab_correlators_at_timestep.sh 20')
    os.system('grab_correlators_at_timestep.sh 50')
    os.system('grab_correlators_at_timestep.sh 100')
    os.system('grab_correlators_at_timestep.sh 500')
    os.system('grab_correlators_at_timestep.sh 1000')

    # Construct RDF input file
    print 'Calculating RDF, MSD, VACF, & VDOS... '
    out = open('RDF.in','w')
    out.write(fnameSkip+'\n')
    out.write('RDF.dat\n')
    out.write('N\nN\n')
    alat_min = min([a,b,c])
    out.write(str(alat_min/2.0)+'\n')
    out.write('0.02\n')
    out.write('0\n')
    out.write(str(a)+', '+str(b)+', '+str(c)+'\n')
    out.close()

    # Construct msd.x input file
    out = open('msd.in','w')
    out.write(fnameSkip+'\n')
    out.close()

    # Construct MSD input file
    out = open('MSD.in','w')
    out.write(fnameSkip+'\n')
    out.write('MSD.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('100\n')
    out.write('5\n')
    out.write(' \n')
    out.close()

    # Construct VACF input file
    out = open('VACF.in','w')
    out.write(fnameSkip+'\n')
    out.write('VACF.dat\n')
    out.write(str(step_in_au)+'\n')
    out.write('1\n')
    out.write('5\n')
    out.write('ALL\n')
    out.close()

    # Construct VDOS input file
    out = open('VDOS.in','w')
    out.write('VACF.dat\n')
    out.write(str(step_in_au)+'\n')
    out.close()

    out = open('gibbs.in','w')
    out.write('dos.dat\n')
    out.write('POSCAR\n')
    temp = commands.getoutput("grep TEBEG INCAR | awk '{print $3}'")
    p_avg = commands.getoutput("grep Final pressure.blocker | awk '{print $4}'")
    e_avg = commands.getoutput("grep Final energy.blocker | awk '{print $4}'")
    out.write(temp+', '+p_avg+', '+e_avg+'\n')
    out.close()
                    
    # Calculate RDF, MSD, diffusion, VACF, VDOS, etc
    os.system('nice -15 RDF < RDF.in > RDF.out')
    os.system('min_gr.py RDF.dat')
    min_gr = commands.getoutput("tail min_gr.dat | awk '{print $1}'")
    os.system('nice -15 get_coord.pl TRAJEC.cnn '+min_gr)
    os.system('av_coor_to_column.py av_coordination_rc_'+min_gr+'.dat')

    os.system('nice -15 VACF < VACF.in > VACF.out')
    os.system('min_vacf.py VACF.dat '+str(step_in_fs))
    min_vacf = commands.getoutput("tail min_vacf.dat | awk '{print $2}'")
    os.system('grab_correlators_at_timestep.sh '+min_vacf)

    os.system('VDOS.x < VDOS.in > VDOS.out')
    os.system("awk '{print $1,$2*"+str(natom)+"}' VDOS.dat > dos.dat")
    os.system('gibbs-solid.x < gibbs.in > gibbs.dat')
    os.system('vdos_to_wavenumber.sh VDOS.dat')
    os.system('running_avg.py VDOS.dat 2 10')

    os.system('nice -15 msd.x < msd.in > msd.out')
    os.system('nice -15 MSD < MSD.in > MSD.out')
    os.system('nice -15 diffusion.py MSD.dat')

    os.system('rm -f VDOS.in msd.in fort.* diss.dat MSD.out')
    os.system('mv -f DiffCoef.data DiffCoef.dat')
    print 'done'

    # Put everything in a separate analysis sub-directory
    print 'Finalizing... '
    os.system('date > date')

    # Copy any relevant plotting scripts to analysis/
    os.system('cp /home/boates/software/gnuplot_scripts/c1-16.plt ./')
    print 'Analysis complete'

if __name__ == '__main__':
    main()
