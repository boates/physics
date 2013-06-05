#!/usr/bin/env python
"""
finite-temperature_POSCAR.py
Author: Brian Boates

Loosely based on Isaac's old finite-temperature.py
"""
import sys, random, numpy

def main():

    # Retrieve user input
    try:
        pos = sys.argv[1]
        nSteps  = int(sys.argv[2])
        smear   = float(sys.argv[3]) / 1000000.0
    except:
        print '\n usage: '+sys.argv[0]+'  supercell.POSCAR  nSteps  smear'
        print '\n Using nSteps >= 1000 and smear ~ TEMP in K works well \n'
        sys.exit(0)

    # Read POSCAR file header
    p = open(pos,'r')
    p.readline()
    alat = float(p.readline())
    ax, ay, az = p.readline().split()
    bx, by, bz = p.readline().split()
    cx, cy, cz = p.readline().split()
    natom = int(p.readline())
    p.readline()

    # Build the lattice vectors
    ax, ay, az = float(ax)*alat, float(ay)*alat, float(az)*alat
    bx, by, bz = float(bx)*alat, float(by)*alat, float(bz)*alat
    cx, cy, cz = float(cx)*alat, float(cy)*alat, float(cz)*alat
    a = numpy.array( [ax, ay, az] )
    b = numpy.array( [bx, by, bz] )
    c = numpy.array( [cx, cy, cz] )

    # Retrieve the reduced atomic coordinates
    atoms = []
    for i in range(natom):
        x, y, z = p.readline().split()
        atoms.append( [float(x), float(y), float(z)] )
    atoms = numpy.array(atoms)
        
    # Set the thermal smearing parameters
    x_smear = smear  
    y_smear = smear 
    z_smear = smear 

    # Open the output xyz file
    xyz = 'SMEAR_'+sys.argv[3]+'K.xyz'
    out = open(xyz,'w')

    # Write the smeared coordinates to file
    for tstep in range(nSteps):

         out.write(str(natom) + '\n')
         out.write(str(tstep + 1) + '\n')
     
         for atom in atoms:
         
              x_r = random.gauss(0.0, x_smear)
              y_r = random.gauss(0.0, y_smear)
              z_r = random.gauss(0.0, z_smear)

              x = ax*(atom[0]+x_r) + bx*(atom[1]+y_r) + cx*(atom[2]+z_r)
              y = ay*(atom[0]+x_r) + by*(atom[1]+y_r) + cy*(atom[2]+z_r)
              z = az*(atom[0]+x_r) + bz*(atom[1]+y_r) + cz*(atom[2]+z_r)

              out.write('H '+str(x)+' '+str(y)+' '+str(z)+'\n')
               
    out.close()


if __name__ == '__main__':
     main()
