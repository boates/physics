#!/usr/bin/env python
"""
Find the '1st minimum' in g(r) for coordination calculations
- Calculates g(r) from RDF
- Does a running average using 10 points and my running_avg.py code
- Returns '1st minimum' in g(r) in angstroms
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        xyz = sys.argv[1]
        alat = float(sys.argv[2])
    except:
        print '\n usage: '+sys.argv[0]+'  TRAJEC.xyz  alat(angstroms)\n'
        sys.exit(0)

    # Create the RDF input file
    out = open('rdf.in','w')
    out.write(xyz+'\n')
    out.write('rdf.dat\nN\nN\n')
    out.write(str(alat/2.0)+'\n')
    out.write('0.02\n0\n')
    out.write(str(alat)+','+str(alat)+','+str(alat)+'\n')
    out.close()

    # Calculate g(r) (remove tmp files)
    os.system('RDF < rdf.in > rdf.out')
    os.system('rm -f rdf.in rdf.out')

    # Calculate the running average of g(r) (for smoothing purposes)
    os.system('running_avg.py rdf.dat 2 10')

    # Open the averaged g(r) file
    f = open('rdf.dat.avg','r')
    lines = f.readlines()
    f.close()
    r, g = [], []
    for line in lines:
        row = line.split()
        r.append( float(row[0]) )
        g.append( float(row[1]) )

    # Remove the beginning of g(r) so as not to confuse with '1st minimum'
    gMax = g.index(max(g))
    gNew = g[gMax:-1]
    rNew = r[gMax:-1]

    # Find the '1st minimum' in g(r)
    gMin = gNew.index(min(gNew))
    rMin = rNew[gMin]

    # Print the result
    print '\n', rMin, 'angstroms is the location of the 1st minimum in g(r)\n'

    # Remove rdf.dat and rdf.dat.avg
    os.system('rm -f rdf.dat rdf.dat.avg')


if __name__ == '__main__':
    main()
