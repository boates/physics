#!/usr/bin/env python
"""
Create input files for and run findsym based on
POSCAR / CONTCAR format
"""
import os, sys, commands, glob, math

def main():

    # Retrieve user input
    try:
        car = open(sys.argv[1],'r')
        tol = sys.argv[2]
    except:
        print '\n usage: '+sys.argv[0]+' CONTCAR(or POSCAR) tolerance(0.01)\n'
        sys.exit(0)

    # Read CONTCAR file
    header = car.readline()
    alat = float(car.readline())
    ax, ay, az = car.readline().split()
    ax, ay, az = float(ax)*alat, float(ay)*alat, float(az)*alat
    bx, by, bz = car.readline().split()
    bx, by, bz = float(bx)*alat, float(by)*alat, float(bz)*alat
    cx, cy, cz = car.readline().split()
    cx, cy, cz = float(cx)*alat, float(cy)*alat, float(cz)*alat
    natoms = car.readline().split()
    car.readline()

    # Determine magnitude of lattice vectors
    mag_a = ( ax**2 + ay**2 + az**2 )**0.5
    mag_b = ( bx**2 + by**2 + bz**2 )**0.5
    mag_c = ( cx**2 + cy**2 + cz**2 )**0.5

    # Determine angles between lattice vectors
    theta_ab = math.acos( (ax*bx + ay*by + az*bz) / (mag_a*mag_b) )
    theta_ac = math.acos( (ax*cx + ay*cy + az*cz) / (mag_a*mag_c) )
    theta_bc = math.acos( (bx*cx + by*cy + bz*cz) / (mag_b*mag_c) )

    # Create aptly named string variables for output
    a, b, c = str(mag_a), str(mag_b), str(mag_c)
    alpha = str(theta_bc*180.0/math.pi)
    beta  = str(theta_ac*180.0/math.pi)
    gamma = str(theta_ab*180.0/math.pi)

    # Find total number of atoms
    N = 0
    for i in range(len(natoms)):
        N += int(natoms[i])

    # Write findsym input file
    out = open('findsym_tol'+tol+'.in','w')
    out.write('Comment\n')
    out.write(tol+'\n')
    out.write('2\n')
    out.write(a+' '+b+' '+c+' '+alpha+' '+beta+' '+gamma+'\n')
    out.write('2\n')
    out.write('P\n')
    out.write(str(N)+'\n')
    for i in range(len(natoms)):
        out.write(str(natoms[i])+'*'+str(i+1)+'\n')
    for i in range(N):
        out.write(car.readline())
    car.close()
    out.close()

    # Run findsym
    os.system('findsym < findsym_tol'+tol+'.in > findsym_tol'+tol+'.out')


if __name__ == '__main__':
    main()
