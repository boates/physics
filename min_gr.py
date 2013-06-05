#!/usr/bin/env python
"""
Find the '1st minimum' in g(r) for coordination calculations
- Does a running average using 10 points and my running_avg.py code
- Returns '1st minimum' in g(r) in angstroms
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        rdf = sys.argv[1]
    except:
        print '\n usage: '+sys.argv[0]+'  RDF.dat\n'
        sys.exit(0)

    # Calculate the running average of g(r) (for smoothing purposes)
    os.system('running_avg.py '+rdf+' 2 6')

    # Open the averaged g(r) file
    f = open(sys.argv[1]+'.avg','r')
    lines = f.readlines()
    f.close()
    r, g = [], []
    for line in lines:
        row = line.split()
        r.append( float(row[0]) )
        g.append( float(row[1]) )

    # Remove the beginning of g(r) so as not to confuse with '1st minimum'
    g0 = g[0]
    while g[0] == g0:
        g.pop(0)
        r.pop(0)
    gNew = g
    rNew = r

    # Find the '1st minimum' in g(r)
    for i in range(2,len(g)+-2):
        if g[i-3] > g[i-2] > g[i-1] > g[i]:
            if g[i] < g[i+1] < g[i+2] < g[i+3]:
                gMin = i
                break
        
#    gMin = gNew.index(min(gNew))
    rMin = r[gMin]

    # Print the result
#    out = open('min_gr.dat','w')
#    out.write(str(rMin)+'\n')
#    out.close()
    print rMin

    # Remove rdf.dat and rdf.dat.avg
#    os.system('rm -f '+sys.argv[1]+'.avg')


if __name__ == '__main__':
    main()
