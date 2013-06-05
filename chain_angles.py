#!/usr/bin/env python
"""
chain angles code test in python
"""
import os, sys, commands, glob
import numpy

def pbc_round(x1,y1,z1,x2,y2,z2,ax,ay,az,verbose=False):

    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2

    s = [dx/ax,dy/ay,dz/az]
    p = [int(s[0]),int(s[1]),int(s[2])]

    for i in range(len(s)):
        if abs(s[i]-p[i]) >= 0.5:
            if s[i] > 0: p[i] += 1
            if s[i] < 0: p[i] -= 1

    dx -= ax*p[0]
    dy -= ay*p[1]
    dz -= az*p[2]

    r = (dx**2+dy**2+dz**2)**0.5

    if verbose:
        return dx, dy, dz
    else:
        return r

def main():

    # Retrieve input
    try:
        f = open(sys.argv[1],'r')
        Nlines = len(f.readlines()) - 10
        f.close()
        f = open(sys.argv[1],'r')
        rcut = float(sys.argv[2])
        NbinsT = int(sys.argv[3])
        NbinsP = int(sys.argv[4])
    except:
        print '\n usage: '+sys.argv[0]+' TRAJEC.cbn, rcut (angstroms), NbinsT, NbinsP\n'
        sys.exit(0)

    # Read file header
    f.readline() 
    ax = float(f.readline().split()[3])*0.529177
    ay = float(f.readline().split()[3])*0.529177
    az = float(f.readline().split()[3])*0.529177
    natom = int(f.readline().split()[3])
    f.readline() 
    f.readline() 
    f.readline() 
    f.readline() 
    f.readline()

    # Open CHAIN.cbn file and write header
    cbn = open('CHAIN.cbn','w')
    cbn.write('# cbn file created with chain_angles.py\n')
    cbn.write('# a =   '+str(ax/0.529177)+'\n')
    cbn.write('# b =   '+str(ay/0.529177)+'\n')
    cbn.write('# c =   '+str(az/0.529177)+'\n')
    cbn.write('# number_of_particles =        '+str(natom)+'\n')
    cbn.write('# number_of_neighbors =        '+str(1)+'\n')
    cbn.write('#\n')
    cbn.write('#\n')
    cbn.write('#\n')
    cbn.write('# units = bohr\n')

    # Initialize the histograms
    theta_hist = numpy.zeros(NbinsT,dtype=numpy.float)
    phi_hist = numpy.zeros(NbinsP,dtype=numpy.float)

    for i in range(Nlines/natom):

        X, Y, Z, nn, typat = [], [], [], [], []
        chains, chain_list = [], []

        for j in range(natom):

            line = f.readline().split()
            typat.append(line[0])
            X.append(float(line[1])*0.529177)
            Y.append(float(line[2])*0.529177)
            Z.append(float(line[3])*0.529177)
            nn.append(int(line[4]))

        # Loop over atoms to find neighbors within cutoff
        for j in range(natom):

            chain_list = []

            for k in range(natom):
                
                r = pbc_round(X[j],Y[j],Z[j],X[k],Y[k],Z[k],ax,ay,az)
                
                if r <= rcut:

                    chain_list.append(k)

            # Involved chain sorting and appending algorithm
            attached = False
            for chain in chains:

                for atom in chain_list:

                    if atom in chain:

                        repeats = [c for c in chain_list if c in chain]

                        for r in repeats:

                            chain_list.pop(chain_list.index(r))

                        chains[chains.index(chain)] += chain_list

                        attached = True
                            
            if not attached:

                chains.append(chain_list)
            
        for chain in chains:

            for ch in chains:

                for atom in chain:

                    if atom in ch and chains.index(chain) != chains.index(ch):

                        repeats = [c for c in chain if c in ch]

                        for r in repeats:

                            chain.pop(chain.index(r))

                        chains[chains.index(chain)] += chains[chains.index(ch)]

        for chain in chains:

            if chains.count(chain) > 1:

                chains.pop(chains.index(chain))

        CHAINS = [chain for chain in chains if len(chain) > 2]

        # A strange final double check seems to be necessary...
        for CHAIN in CHAINS:

            for CH in CHAINS:

                for atom in CHAIN:

                    if atom in CH and CHAINS.index(CHAIN) != CHAINS.index(CH):

                        repeats = [C for C in CHAIN if C in CH]

                        for r in repeats:

                            CHAIN.pop(CHAIN.index(r))

                        CHAINS[CHAINS.index(CHAIN)] += CHAINS[CHAINS.index(CH)]

        k = 0
        while k < len(CHAINS):

            if CHAINS.count(CHAINS[k]) == 1:

                k += 1

            else:

                CHAINS.pop(CHAINS.index(CHAINS[k]))

                k = 0
                
        for CHAIN in CHAINS:

            if CHAINS.count(CHAIN) > 1:

                CHAINS.pop(CHAINS.index(CHAIN))

        # Warn of any repeating chain units
        ALL = []
        for CHAIN in CHAINS:
            ALL += CHAIN
        F = []
        for a in ALL:
            if ALL.count(a) > 1 and a not in F:
                F.append(a)

        if len(F) > 0:
            print '\n Warning the following atoms have been found in multiple chains:'
            print F,'\n exiting...\n'
            print 'CHAINS ='
            for CHAIN in CHAINS: print CHAIN
            print
            sys.exit(0)

        # Write the new chain cbn file
        for j in range(natom):
            if j not in ALL:
                cbn.write(typat[j]+'  '+str(X[j]/0.529177)+'  '+str(Y[j]/0.529177)+'  '+str(Z[j]/0.529177)+'  '+str(nn[j])+'\n')
            else:
                cbn.write(typat[j]+'  '+str(X[j]/0.529177)+'  '+str(Y[j]/0.529177)+'  '+str(Z[j]/0.529177)+'  '+str(-1)+'\n')

        # Do the angle calculations
        triples = []
        quads = []
        for chain in chains:

            for p in chain:

                for o in chain:

                    r_p_o = pbc_round(X[p],Y[p],Z[p],X[o],Y[o],Z[o],ax,ay,az)
                    if r_p_o <= rcut:

                        for oo in chain:

                            r_o_oo = pbc_round(X[o],Y[o],Z[o],X[oo],Y[oo],Z[oo],ax,ay,az)
                            r_p_oo = pbc_round(X[p],Y[p],Z[p],X[oo],Y[oo],Z[oo],ax,ay,az)

                            if r_o_oo <= rcut and r_p_oo > rcut:

                                arg = ( (r_p_o**2+r_o_oo**2-r_p_oo**2) / (2.0*r_p_o*r_o_oo) )
                                theta = numpy.arccos(arg)*(180.0/numpy.pi)

                                triple = [p,o,oo]
                                triple.sort()

                                if triple not in triples:

                                    triples.append(triple)
                                    theta_hist[int(theta*(NbinsT/180.0))] += 1

                                for ooo in chain:

                                    r_oo_ooo = pbc_round(X[oo],Y[oo],Z[oo],X[ooo],Y[ooo],Z[ooo],ax,ay,az)
                                    r_o_ooo = pbc_round(X[o],Y[o],Z[o],X[ooo],Y[ooo],Z[ooo],ax,ay,az)
                                    r_p_ooo = pbc_round(X[p],Y[p],Z[p],X[ooo],Y[ooo],Z[ooo],ax,ay,az)

                                    if r_oo_ooo <= rcut and r_o_ooo > rcut and r_p_ooo > rcut:

                                        Ax, Ay, Az = pbc_round(X[oo],Y[oo],Z[oo],X[p],Y[p],Z[p],ax,ay,az,verbose=True)
                                        Bx, By, Bz = pbc_round(X[oo],Y[oo],Z[oo],X[o],Y[o],Z[o],ax,ay,az,verbose=True)
                                        Cx, Cy, Cz = pbc_round(X[oo],Y[oo],Z[oo],X[ooo],Y[ooo],Z[ooo],ax,ay,az,verbose=True)

                                        AXB = numpy.array([Ay*Bz-Az*By,Az*Bx-Ax*Bz,Ax*By-Ay*Bx])
                                        AXBdC = AXB * numpy.array([Cx,Cy,Cz])
                                        numerator = AXBdC[0] + AXBdC[1] + AXBdC[2]
                                        denominator = (AXB[0]**2+AXB[1]**2+AXB[2]**2)**0.5 * (Cx**2+Cy**2+Cz**2)**0.5
                                        phi = 90.0 - numpy.arccos(numerator/denominator)*(180.0/numpy.pi)

                                        if phi > 90.0: phi -= 180.0
                                        quad = [p,o,oo,ooo]
                                        quad.sort()

                                        if quad not in quads:

                                            quads.append(quad)
                                            phi_hist[int(round(phi*(NbinsP/180.0)))] += 1

    # Write theta info to file
    out = open('chain_theta.dat','w')
    out.write('# theta, hist\n')
    for i in range(len(theta_hist)):
        out.write(str(i*180.0/NbinsT)+'   '+str(theta_hist[i]/numpy.sum(theta_hist))+'\n')
    out.close()

    # Write phi info to file
    out = open('chain_phi.dat','w')
    out.write('# phi, hist\n')
    phi_list = numpy.array(range(len(phi_hist)))*180.0/NbinsP
    for i in range(len(phi_list)):
        if phi_list[i] > 90.0: phi_list[i] -= 180.0
    while phi_list[-1] < phi_list[0]:
        phi_list = numpy.array([phi_list[-1]] + list(phi_list[:-1]))
        phi_hist = numpy.array([phi_hist[-1]] + list(phi_hist[:-1]))
    for i in range(len(phi_hist)):
        out.write(str(phi_list[i])+'   '+str(phi_hist[i]/numpy.sum(phi_hist))+'\n')
    out.close()

    cbn.close()

    # Generate the CHAIN_LABEL.dat file for visualization of chains
#    os.system('echo "CHAIN.cbn" > chain.in')
#    os.system('chain_label.x < chain.in > chain_label.out')
#    os.system('rm -f chain_label.out')
                            

if __name__ == '__main__':
    main()
