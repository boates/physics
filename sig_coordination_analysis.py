#!/usr/bin/env python
"""
Loop over the different configurations for which conductivities
have been calculated and determine their coordination fractions
"""
import os, sys, commands, glob, math

def alat_from_rs(rs,natom,nvalence=5):
    """
    Give rs and natom and returns lattice parameter in angstroms
    """
    a0 = 5.291772108E-09
    volume = (1.0/3.0)*(4.0*natom*nvalence*math.pi)*(rs*a0)**3
    length = volume**(1.0/3.0)

    return length/1E-8


def main():

    # Retrieve user input
    try:
        rs = float(sys.argv[1])
    except:
        print '\n usage: '+sys.argv[0]+'  rs\n'
        sys.exit(0)

    # Locate individual xyz-snapshot and conductivity files
    xyzs = glob.glob('snapshot*.xyz')
    sigs = glob.glob('sig*.dat')

    # Match the files (in case more of one than the other are present)
    xyzs_stripped, xyzs_stripped_int = [], []
    for xyz in xyzs:
        xyzs_stripped.append(xyz.lstrip('snapshot').rstrip('.xyz'))
        xyzs_stripped_int.append(int(xyz.lstrip('snapshot').rstrip('.xyz')))
    sigs_stripped, sigs_stripped_int = [], []
    for sig in sigs:
        sigs_stripped.append(sig.lstrip('sig').rstrip('.dat'))
        sigs_stripped_int.append(int(sig.lstrip('sig').rstrip('.dat')))

    # It is always assumed that if a sig file exists, its snapshot is there also
    xyzs_new, sigs_new = [], []
    for i in range(len(sigs_stripped)):
        if sigs_stripped[i] in xyzs_stripped:
            xyzs_new.append(xyzs[xyzs_stripped.index(sigs_stripped[i])])
            sigs_new.append(sigs[i])

    # Determine system variables
    natom = int(commands.getoutput("head -1 "+xyzs[0]))
    alat_ang = str(alat_from_rs(rs,natom))
    alat_bohr = str( float(alat_ang) / 0.5291772 )

    ### xyzs_new and sigs_new are now ordered with eachother ###

    # Perform analysis and write to file
    N = len(sigs_new)
    c1_sum, c2_sum, c3_sum = 0.0, 0.0, 0.0
    out = open('cond_vs_coordination.dat','w')
    out.write('# 1coord, 2coord, 3coord, conductivity(ohm*cm)^-1\n')
    for i in range(len(sigs_new)):

        # Grab w=0 conductivity
        sig_zero = commands.getoutput("head -2 "+sigs_new[i]+" | tail -n-1 | awk '{print $4}'")

        # Calculate the '1st minimum' in the cooresponding snapshot's g(r)
        min_g = str( float(commands.getoutput("min_gr.py "+xyzs_new[i]+" "+alat_ang).split()[0]) / 0.5291772 )

        # Generate a .cnn file for get_coord.pl
        os.system('nn_dist.py 200 '+alat_bohr+' '+alat_bohr+' '+alat_bohr+' '+xyzs_new[i])

        # Caclculate the coordination fractions, keep only for 1, 2, & 3 coordinated atoms
        os.system('get_coord.pl TRAJEC.cnn '+min_g)
        coords = commands.getoutput("tail -n-1 av_coordination_rc_*.dat").split()[1:4]

        # Remove ALL intermediate files
        os.system('rm -f av_coordination_rc_*.dat coordination_rc_*.dat')
        os.system('rm -f nn_average.hist TRAJEC.cnn')

        # Normalize coordination to natom
        c1, c2, c3 = float(coords[0])/natom, float(coords[1])/natom, float(coords[2])/natom
        c1_sum += c1
        c2_sum += c2
        c3_sum += c3

        # Write to file
        out.write(str(c1)+' '+str(c2)+' '+str(c3)+' '+sig_zero+'\n')

    out.close()

    # Run sig averaging and write overall average for the density to file
    os.system('avg_sig.py sig')
    sig_zero_avg = commands.getoutput("head -2 avg_sig.dat | tail -n-1 | awk '{print $2}'")

    out = open('avg_cond_vs_coordination.dat','w')
    out.write('# avg_1coord, avg_2coord, avg_3coord, avg_conductivity (ohm*cm)^-1\n')
    out.write(str(c1_sum/N)+' '+str(c2_sum/N)+' '+str(c3_sum/N)+' '+sig_zero_avg+'\n')
    out.close()

if __name__ == '__main__':
    main()
