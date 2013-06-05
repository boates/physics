#!/usr/bin/env python
"""
pdos.py
Author Brian Boates

Average a VASP PROCAR file by band, ion, and s/p/d character
"""
import sys, commands
import numpy
from scipy.special import erf

def main():

    # retrieve user input
    try:
        f = open(sys.argv[1],'r')
        atom_i = int( sys.argv[2] )
        atom_f = int( sys.argv[3] )
        band_i = int( sys.argv[4] )
        band_f = int( sys.argv[5] )
        nbins  = int( sys.argv[6] )
        sigma  = float( sys.argv[7] )
        if atom_i == atom_f == 0:
            atom_list = commands.getoutput('paste '+sys.argv[8]).split()
            atoms = [int(a) for a in atom_list]
        else:
            atoms = range(atom_i,atom_f+1)
    except:
        print '\n usage: '+sys.argv[0]+' PROCAR atom_i atom_f band(s) nbins sigma\n'
        print ' atom_i:  first atom index in range to average'
        print ' atom_f:  final atom index in range to average'
        print ' band_i: first band in range to consider'
        print ' band_f: final band in range to consider'
        print ' nbins:  number of bins for histograms'
        print ' sigma:  smearing parameter \n'
        print ' auxillary use: if atom_i and atom_f are both set to zero'
        print ' give 8th kwarg as filepath to a single column file of indices for'
        print ' which atoms to use (NOTE: indexing starts at atom 1, NOT zero \n'
        sys.exit(0)

    # determine range of band energies (max/min)
    band_energies = commands.getoutput("grep band "+sys.argv[1]+" | grep -v points | awk '{print $5}'").split()
    ebands = [float(e) for e in band_energies]
    emin = round( min(ebands) -1 )
    emax = round( max(ebands) +1 )
    de = emax - emin

    # header info
    f.readline() # header
    line   = f.readline().split()
    nkpts  = int( line[3] )
    nbands = int( line[7] )
    natom  = int( line[-1] )

    # create necessary arrays
    s = numpy.zeros(nbins+1)
    p = numpy.zeros(nbins+1)
    d = numpy.zeros(nbins+1)
    tot = numpy.zeros(nbins+1)
    kpts, weights, ibands, ebands, obands = [], [], [], [], []
    
    # loop over kpts
    for i in range(nkpts):

        # append new lists to band related arrays
        ibands.append( [] )
        ebands.append( [] )
        obands.append( [] )

        # blank line
        f.readline()

        # kpt info (dumb vasp formatting fix too)
        line = f.readline().replace('00-0.','00 -0.').split()
        kraw, weight = line[3:6], float(line[-1])
        kpt = [float(k) for k in kraw]
        kpts.append(kpt)
        weights.append(weight)

        # loop over bands
        for j in range(nbands):

            # blank line
            f.readline()

            # band info
            line = f.readline().split()
            iband = int(line[1])
            eband = float(line[4])
            oband = float(line[-1])

            # whether or not to consider this band
            if band_i <= iband <= band_f:
                ibands[i].append( iband )
                ebands[i].append( eband )
                obands[i].append( oband )

            # blank line
            f.readline()

            # header line (ion, s, p, d, tot)
            f.readline()

            print 'kpt: '+str(i+1)+'/'+str(nkpts),'band: '+str(j)+'/'+str(nbands)  ##################

            # loop over atoms
            for a in range(natom):

                # read in info for current atom
                line = f.readline().split()

                # if considering this band
                if band_i <= iband <= band_f:

                    # if considering this atom (shift indices by one for python)
                    if a+1 in atoms:

                        # loop over energy bins
                        for b in range(nbins):

                            # simplify some variable names
                            E = emin + b*(de/nbins)
                            w = weight * 2.0  # two e per state
                            arg1 = (E + de/nbins - eband) / (2**0.5 * sigma)
                            arg2 = (E - eband) / (2**0.5 * sigma)

                            # build each DOS with smearing
                            s[b]   += w/2. *( erf(arg1) - erf(arg2) ) *float(line[1])
                            p[b]   += w/2. *( erf(arg1) - erf(arg2) ) *float(line[2])
                            d[b]   += w/2. *( erf(arg1) - erf(arg2) ) *float(line[3])
                            tot[b] += w/2. *( erf(arg1) - erf(arg2) ) *float(line[4])

###########################################################################
#                        # determine current bin based on band energy
#      OLD WAY           bin = int( (eband-emin)*(nbins/de) +0.5 )
#      W/O SMEARING      s[bin] += float( line[1] )
#                        p[bin] += float( line[2] )
#                        d[bin] += float( line[3] )
#                        tot[bin] += float( line[4] )
###########################################################################

            # tot line
            f.readline()

        # blank line
        f.readline()

    # create energy array
    e = []
    for i in range(nbins+1):
        e.append(emin+i*de/nbins)

    # write s, p, d, and tot DOS for desired bands/atoms to file
    out = open('pdos.dat','w')
    out.write('# from file '+sys.argv[1]+', using atoms '+str(atom_i)+'-'+str(atom_f))
    out.write(' and bands '+str(band_i)+'-'+str(band_f)+'\n')
    out.write('# E(eV), tot, s, p, d\n')
    for i in range(nbins+1):
        out.write(str(e[i])+' '+str(tot[i])+' '+str(s[i])+' '+str(p[i])+' '+str(d[i])+'\n')
    out.close()
        


if __name__ == '__main__':
    main()
