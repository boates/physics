#!/usr/bin/env python
"""
Calculate bond lengths and frequencies/periods
"""
import sys, Numeric
import pylab, PSD

def pbc_round(input_value):
     i = int(input_value)
     if (abs(input_value-i) >= 0.5):
         if (input_value > 0): i+= 1
         if (input_value < 0): i-= 1
     return i

def nearest_neighbor(fcbn):
    """
    Open a cnn file and return a list of nearest neighbors
    from the initial configuration
    """
    # Open file and read header info
    f = open(fcbn,'r')
    for i in range(9):
        if i == 4:
            natom = int(f.readline().split()[-1])
        f.readline()
    
    # Read in the nn info for the first configuration
    nn = []
    for i in range(natom):
        nn.append( int(f.readline().split()[4]) )

    f.close()

    return nn

def main():

    # Check for required user input
    try:
        fname = sys.argv[1]
        step = int(sys.argv[2])
    except IndexError:
        print '\nusage: '+sys.argv[0]+' TRAJEC.cbn step_in_au\n'
        sys.exit(0)

    # Get the neighbor list from cbn file
    nn = nearest_neighbor(fname)

    # Open and read header from cbn file
    f = open(fname,'r')
    for i in range(8):
        if i == 3:
            alat = float(f.readline().split()[-1])*0.5291772 # Convert to angstroms
            natom = int(f.readline().split()[-1])
        f.readline()
    lines = f.readlines()
    f.close()

    # Calculate all the bond lengths over time
    dr = [[] for i in range(len(nn))]
    Nconfigs = len(lines)/len(nn)
    for j in range(len(lines)/(len(nn))):

        # Report progress to user
        if j % 200 == 0 and j != 0: print int(round((float(j)/Nconfigs)*100)), '% ...'
        
        for i in range(natom):

            # Current atom's coords and its partner's
            typat1,x1,y1,z1,dummy1 = lines[i].split()
            typat2,x2,y2,z2,dummy2 = lines[nn[i]].split()

            # Caclulate bond length using the minimum image convention
            dx = (float(x1) - float(x2)) * 0.5291772
            dy = (float(y1) - float(y2)) * 0.5291772
            dz = (float(z1) - float(z2)) * 0.5291772

            dx -= alat*pbc_round(dx/alat)
            dy -= alat*pbc_round(dy/alat)
            dz -= alat*pbc_round(dz/alat)

            # Calculate bond length in angstroms
            r = ( dx**2 + dy**2 + dz**2 )**0.5
            dr[i].append(abs(r))

        for k in range(natom): lines.pop(0)

    frequencies, periods, bond_lengths = [], [], []

    for dr_array in dr:
        
        # Calculate the PSD
        psd_obj = PSD.psd([i for i in range(len(dr_array))],dr_array)
        f, psd = psd_obj.getPSD()
        f, psd = list(f), list(psd)

        # Dominant frequency (in region of interest)
        max_f = 0
        while max_f == 0:
            max_f = f.pop( psd[1:len(psd)/2].index(max(psd[1:len(psd)/2])) )

        # Corresponding period in femtoseconds
        T = (1.0 / max_f) * (step * 2.4188000000000003e-02)

        frequencies.append(max_f)
        periods.append(T)

        bond_lengths += dr_array

    """
    # Write PSD data to file
    out = open('PSD.dat','w')
    out.write('#f        psd\n')

    # Neglect the zeroth value and the second half
    for i in range(len(f)/2-1):
    out.write(str(f[i+1])+'    '+str(psd[i+1])+'\n')
    out.close()
    """

    # Histograms of values
    pylab.figure(num=1,figsize=(5,3),facecolor='w',edgecolor='k')
    pylab.hist(frequencies, normed=1, label='_nolegend_')
    pylab.savefig('bond_frequencies_hist.png')
    pylab.clf()
    pylab.hist(periods, normed=1, label='_nolegend_')
    pylab.savefig('bond_periods_hist.png')
    pylab.clf()
    pylab.hist(bond_lengths, bins=250, normed=1, label='_nolegend_')
    pylab.savefig('bond_length_hist.png')
    pylab.clf()


if __name__ == '__main__':
    main()
