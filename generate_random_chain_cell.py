#!/usr/bin/env python
"""
Randomly (within criteria) generate a supercell of N chains
in a POSCAR file for relaxation.
"""
import os, sys, commands, glob
import math, Numeric
from random import random

global pi
pi = math.pi

def make_chain(L,rNN,eta):
    """
    Generate the coords for a chain of given length (L) in angstroms
    
    L: number of atoms in chain
    rNN: fixed N-N bond length (1.25 ang)
    min_eta: minimum N-N-N angle (in degrees)
    max_eta: maximum N-N-N angle 
    R_eta: a random variable between 0 and 1 for determining eta angle
    """
    ystep = 2*rNN*math.sin(eta) / 2.0
    zstep = ( rNN**2 - ystep**2 )**0.5
    units = L/2
    chain = []
    for i in range(units):
        chain.append([0.0,2*i*ystep,0.0])
        chain.append([0.0,(2*i+1)*ystep,zstep])

    return chain
    
def volume_from_rs(rs,Nel):
    """
    Return cell volume in angstroms^3.
    rs: density parameter for system
    Nel: number of valence electrons in system
    """
    a0 = 0.5291772  # Bohr radius (angstroms/bohr)
    volume = (4.0*pi/3.0)*Nel * (rs*a0)**3

    return volume

def pbc_round(input_value):
    i = int(input_value)
    if (abs(input_value-i) >= 0.5):
        if (input_value > 0): i+= 1
        if (input_value < 0): i-= 1
    return i
            

def main():

    # Retrieve user input
    try:
        INPUT = sys.argv[1]
        f = open(INPUT,'r')
        L = int(f.readline())
        Nchain = int(f.readline())
        zval = int(f.readline())
        rNN = float(f.readline())
        rCC = float(f.readline())
        etas = f.readline().split()
        eta_min, eta_max = float(etas[0])*pi/180.0/2.0, float(etas[1])*pi/180.0/2.0
        rot_max = float(f.readline())*pi/180.0
        theta_acs = f.readline().split()
        theta_ac_min, theta_ac_max = float(theta_acs[0])*pi/180.0, float(theta_acs[1])*pi/180.0
        rss = f.readline().split()
        rs_min, rs_max = float(rss[0]), float(rss[1])
        ca_ratios = f.readline().split()
        ca_ratio_min, ca_ratio_max = float(ca_ratios[0]), float(ca_ratios[1])
    except:
        print '\n usage: '+sys.argv[0]+' INPUT\n'
        print ' INPUT should be formatted as follows:'
        print ' number_of_atoms_in_each_chain'                   # 8
        print ' number_of_chains'                                # 9
        print ' zval'                                            # 5
        print ' N-N_distance (angstroms)'                        # 1.25
        print ' min_chain-chain_distance (angstroms)'            # 2.50
        print ' min_N-N-N_angle (2*eta from notes)   max_2*eta'  # 90 120
        print ' max_rotation_angle_of_chains_in_ac_plane (deg)'  # 20
        print ' min_a-c_angle   max_a-c_angle'                   # 60 90
        print ' min_rs   max_rs'                                 # 1.20 1.25
        print ' min_c/a_ratio   max_c/a_ratio\n'                 # 0.5 2.0
        sys.exit(0)

    # Create chain via random number
    R_eta = random()
    eta = eta_min + (eta_max-eta_min)*R_eta
    chain = make_chain(L,rNN,eta)

    # Determine lattice vector, b, from chain length
    # Chains run parallel to b
    bx = 0.0
    by = L*chain[1][1]
    bz = 0.0
    b = Numeric.array([bx, by, bz])
    mag_b = (bx**2 + by**2 + bz**2)**0.5

    # Determine rs via random number method
    R_rs = random()
    rs = rs_min + (rs_max-rs_min)*R_rs

    # Calculate volume from generated rs
    natom = Nchain*L
    Nel = natom*zval
    volume = volume_from_rs(rs,Nel)

    # Determine area of ac-plane
    A_ac = volume / by

    # Generate c/a lattice vector ratio
    R_ratio = random()
    ca_ratio = 0.5 + 1.5*R_ratio

    # Generate angle between a and c lattice vectors
    R_theta_ac = random()
    theta_ac = theta_ac_min + (theta_ac_max-theta_ac_min)*R_theta_ac

    # Determine ax and thus a
    ax = ( A_ac / (ca_ratio*math.sin(theta_ac)) )**0.5
    ay = 0.0
    az = 0.0
    a = Numeric.array([ax, ay, az])
    mag_a = (ax**2 + ay**2 + az**2)**0.5    

    # Determine c from a
    mag_c = ca_ratio * mag_a
    cx = mag_c*math.cos(theta_ac)
    cy = 0.0
    cz = mag_c*math.sin(theta_ac)
    c = Numeric.array([cx, cy, cz])

    #==========================================#
    # NOW HAVE a, b, c, rs, chain, A_ac, and V #
    #                                          #
    # Now to start randomly positioning chains #
    #                                          #
    # - Chains must be at least rCC apart      #
    # - Must lie within supercell boundaries   #
    # - Define ac-plane atom positions and     #
    #   build chains from those                #
    #==========================================#

    ac_plane = [[0.0,0.0,0.0]]
    loops = 0
    while len(ac_plane) < Nchain:
        if loops == 1E5:
            print '\nLooped',loops,'times - not going to happen at this point - exiting...\n'
            sys.exit(0)
        loops += 1
        k = 0
        R_a, R_c = random(), random()
        xn = R_a * ax + R_c * cx
        zn = R_c * cz
        xn_prime = R_a
        zn_prime = R_c
        for atom in ac_plane:

            # Change coordinate system to cell vectors
            xo, yo, zo = atom
            zo_prime = zo / cz
            xo_prime = (xo - zo_prime*cx) / ax

            # Calculate distance w/ minimum image convention
            dx_prime = xn_prime - xo_prime
            dz_prime = zn_prime - zo_prime
            
            dx_prime -= pbc_round(dx_prime)
            dz_prime -= pbc_round(dz_prime)

            dx = dx_prime*ax + dz_prime*cx
            dz = dz_prime*cz
            
            d = (dx**2 + dz**2)**0.5

            # Check distance criteria for chains
            if d < rCC:
                k = 1

        # If passes criteria, append.
        if k == 0:
            ac_plane.append([xn, 0.0, zn])
            print round(len(ac_plane) / float(Nchain) * 100.0, 1),'% complete'

    # Write coordinates to xyz file and POSCAR file
    out = open('chains.xyz','w')
    out.write(str(natom)+'\n1\n')
    pos = open('chains.POSCAR','w')
    pos.write('N chains, rs='+str(rs)+'\n')
    pos.write('    '+str(ax)+'\n')
    pos.write('     1.0000000000000000    0.0000000000000000    0.0000000000000000\n')
    pos.write('     0.0000000000000000    '+str(by/mag_a)+'    0.0000000000000000\n')
    pos.write('     '+str(cx/mag_a)+'    0.0000000000000000    '+str(cz/mag_a)+'\n')
    pos.write(' '+str(natom)+'\n')
    pos.write('Direct\n')
    for shift in ac_plane:
        R_rotate = 2*(random()-0.5)                         # Number between -1 and 1
        rot = R_rotate * rot_max # + int(random()+0.5)*pi   # 'up' or 'down'
        for ion in chain:
            x = ion[0] + shift[0]
            y = ion[1] + shift[1]
            z = ion[2] + shift[2]
            # Randomly rotate chain in the ac-plane
            if ion[2] != chain[0][2]:
                x = (chain[0][0] + shift[0]) - rNN*math.sin(rot)
                z = (chain[0][2] + shift[2]) + rNN*math.cos(eta)*math.cos(rot)
            out.write('N '+str(x)+' '+str(y)+' '+str(z)+'\n')
            pos.write('  '+str((x - z*cx/cz)/ax)+'  '+str(y/by)+'  '+str(z/cz)+'\n')
    out.close()
    pos.close()

    # Write relevant information to log file
    out = open('chains.log','w')
    out.write(' # data generated based on following input:\n')
    out.write(' number_of_atoms_in_each_chain:                   '+str(L)+'\n')
    out.write(' number_of_chains:                                '+str(Nchain)+'\n')
    out.write(' zval:                                            '+str(zval)+'\n')
    out.write(' N-N_distance (angstroms):                        '+str(rNN)+'\n')
    out.write(' min_chain-chain_distance (angstroms):            '+str(rCC)+'\n')
    out.write(' min_N-N-N_angle (2*eta from notes)   max_2*eta:  '+str(eta_min*180.0/pi*2.0)+'   '+str(eta_max*180.0/pi*2.0)+'\n')
    out.write(' max_rotation_angle_of_chains_in_ac_plane (deg)   '+str(rot_max*180.0/pi)+'\n')
    out.write(' min_a-c_angle   max_a-c_angle:                   '+str(theta_ac_min*180.0/pi)+'   '+str(theta_ac_max*180.0/pi)+'\n')
    out.write(' min_rs   max_rs:                                 '+str(rs_min)+'   '+str(rs_max)+'\n')
    out.write(' min_c/a_ratio   max_c/a_ratio:                   '+str(ca_ratio_min)+'   '+str(ca_ratio_max)+'\n\n')
    out.write(' # results from calculations (finished after '+str(loops)+' loops):\n')
    out.write(' lattice vector a =  ('+str(ax)+', '+str(ay)+', '+str(az)+')\n')
    out.write(' lattice vector b =  ('+str(bx)+', '+str(by)+', '+str(bz)+')\n')
    out.write(' lattice vector c =  ('+str(cx)+', '+str(cy)+', '+str(cz)+')\n')
    out.write(' |a|, |b|, |c| =     '+str(mag_a)+', '+str(mag_b)+', '+str(mag_c)+'\n')
    out.write(' c/a ratio =         '+str(ca_ratio)+'\n')
    out.write(' theta_ac =          '+str(theta_ac*180.0/pi)+'\n')
    out.write(' eta =               '+str(eta*180.0/pi)+'\n')
    out.write(' 2*eta =             '+str(2*eta*180.0/pi)+'\n')
    out.write(' rs =                '+str(rs)+'\n')
    out.write(' natom =             '+str(natom)+'\n')
    out.write(' Nel =               '+str(Nel)+'\n')
    out.write(' volume =            '+str(volume)+'\n')
    out.write(' A_ac =              '+str(A_ac)+'\n')
    out.close()


if __name__ == '__main__':
    main()
