#!/usr/bin/env python
"""
BZ_paths.py

Create paths in the Brillouin Zone for anaddb
dispersion calculations
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
        stepsize = float(sys.argv[2])
        lines = f.readlines()
        f.close()
    except:
        print '\n usage: '+sys.argv[0]+' file_w_path_endpoints  stepsize(e.g. 0.01)\n'
        sys.exit(0)

    """
    A file_w_path_endpoints should look like:
    0.00 0.00 0.00
    -0.50 0.50 0.50
    -0.75 0.25 0.25
    -1.00 0.00 0.00
    -1.00 0.00 0.50
    """

    EP = []
    for line in lines:
        row = line.split()
        EP.append([ float(row[0]), float(row[1]), float(row[2]) ])
#    EP.append(EP[-1])
    EP.append([0.1,0.2,0.3])

    paths = open('BZ_paths.dat','w')
    for i in range(len(EP)-1):
        x, y, z = EP[i][0], EP[i][1], EP[i][2]
        x_n, y_n, z_n = EP[i+1][0], EP[i+1][1], EP[i+1][2]
        dx = x_n - x
        dy = y_n - y
        dz = z_n - z

        Nx = int(abs(dx / stepsize))
        Ny = int(abs(dy / stepsize))
        Nz = int(abs(dz / stepsize))

#        if Nx != Ny or Ny != Nz or Nz != Nx:
#            print '\n WARNING distance to endpoints not equal in x, y, z - exiting'
#            print 'From', EP[i], 'to', EP[i+1]
#            sys.exit(0)

        Ns = [Nx, Ny, Nz]
#        print Ns
        Nmax = max([n for n in Ns if n != 0.0])
        Sx = float(Nx) / Nmax * stepsize
        Sy = float(Ny) / Nmax * stepsize
        Sz = float(Nz) / Nmax * stepsize
#        print [Sx,Sy,Sz]
        for i in range(Nmax):
            paths.write(str(x)+' '+str(y)+' '+str(z)+' 1\n')
            if dx != 0:
                x += Sx * round(dx/abs(dx))
            if dy != 0:
                y += Sy * round(dy/abs(dy))
            if dz != 0:
                z += Sz * round(dz/abs(dz))

    paths.write(str(x)+' '+str(y)+' '+str(z)+' 1\n')

    paths.close()


if __name__ == '__main__':
    main()
