#!/usr/bin/env python
"""
POSCAR_angles.py
"""
import sys, math, numpy

def main():

    # Retrieve user input
    try:
        f = open(sys.argv[1],'r')
    except:
        print '\n usage: '+sys.argv[0]+' POSCAR\n'
        sys.exit(1)

    f.readline()
    alat = float(f.readline().strip())

    a = f.readline().split()
    ax = float(a[0])*alat
    ay = float(a[1])*alat
    az = float(a[2])*alat
    a = numpy.array([ax,ay,az])
    a_mag = ( ax**2 + ay**2 + az**2 )**0.5

    b = f.readline().split()
    bx = float(b[0])*alat
    by = float(b[1])*alat
    bz = float(b[2])*alat
    b = numpy.array([bx,by,bz])
    b_mag = ( bx**2 + by**2 + bz**2 )**0.5

    c = f.readline().split()
    cx = float(c[0])*alat
    cy = float(c[1])*alat
    cz = float(c[2])*alat
    c = numpy.array([cx,cy,cz])
    c_mag = ( cx**2 + cy**2 + cz**2 )**0.5

    f.close()

    theta_ab = math.acos( numpy.dot(a,b) / (a_mag*b_mag) )*180./math.pi
    theta_ac = math.acos( numpy.dot(a,c) / (a_mag*c_mag) )*180./math.pi
    theta_bc = math.acos( numpy.dot(b,c) / (b_mag*c_mag) )*180./math.pi

    print
    print 'theta_ab, theta_ac, theta_bc = ', theta_ab, theta_ac, theta_bc
    print 'magnitudes (a, b, c): ', a_mag, b_mag, c_mag
    print 'relative magnitudes (a/b, a/c, b/c): ', a_mag/b_mag, a_mag/c_mag, b_mag/c_mag
    print 

if __name__ == '__main__':
    main()
