#!/usr/bin/env python
"""
BLA.py

Create a histogram of the difference in bond length of
each atoms 1st and 2nd nearest neighbors.
(Originally for investigations of Peierls distortions)
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

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        f.readline()
        a = float(f.readline().split()[-1])
        b = float(f.readline().split()[-1])
        c = float(f.readline().split()[-1])
        natom = int(f.readline().split()[-1])
        nneighbours = int(f.readline().split()[-1])
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        lines = f.readlines()
        f.close()

        nbins = int(sys.argv[2])
    except:
        print '\n usage: '+sys.argv[0]+' TRAJEC.cnn nbins\n'
        sys.exit(0)

    # Define variables
    nConfigs = len(lines) / natom
    BLA = []
    typat, x, y, z, nn1, nn2 = [], [], [], [], [], []

    # Loop over all snapshots
    for i in range(nConfigs):

        # Read in data for each snapshot
        for j in range(natom):
            row = lines[(i*natom+j)].split()
            typat.append( row[0] )
            x.append( float(row[1]) )
            y.append( float(row[2]) )
            z.append( float(row[3]) )
            nn1.append( int(row[4]) )
            nn2.append( int(row[5]) )

        # For each atom calculate 1st and 2nd neighbour distances
        for j in range(natom):
            r1 = pbc_round(x[j],y[j],z[j],x[nn1[j]],y[nn1[j]],z[nn1[j]],a,b,c)
            r2 = pbc_round(x[j],y[j],z[j],x[nn2[j]],y[nn2[j]],z[nn2[j]],a,b,c)

           # Append these difference to the bond length alternation list
            BLA.append( r2 - r1 )

    # Create histogram from list of values
    hist = numpy.histogram(BLA,bins=nbins)

    # Write histogram to file
    out = open('BLA.hist','w')
    for i in range(len(hist[0])):
        out.write(str(hist[1][i])+' '+str(hist[0][i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
