#!/usr/bin/env python
"""
Give the primitive vectors for the rhombohedral R3c structure for a given
lattice constant and angle.
"""
import sys
from math import sin, cos, pi

def main():

    """ Print out the rprim vectors """

    try:
        a = float(sys.argv[1])
        alpha = float(sys.argv[2])
        P = None
    except IndexError:
        print '\nPlease give the lattice constant AND the angle\n'
        sys.exit(0)
    except ValueError:
        if sys.argv[1] == 'N':
            try: P = int(sys.argv[2])
            except: P = -1
            while P not in [0,2,4,6,8,10]:
                print '\nPlease choose from the available choices for P\n'
                try: P = input('Please provide the pressure (either 0, 2, 4, 6, 8, or 10 GPa (LeSar, 1984)): ')
                except: P = -1
        else:
            print '\nMake sure the lattice constant and angle are numbers... obviously.\n'
            sys.exit(0)
            
    lesar_data = {0:(13.32,85.02),2:(12.15,84.87),4:(11.74,84.81),6:(11.48,84.85),\
                  8:(11.28,84.88),10:(11.11,84.97)}
    if P != None:
        a, alpha = lesar_data[P]

    # Convert alpha to radians
    alpha *= pi/180.

    rprim1 = (a*sin(pi/2-alpha)*cos(pi-alpha/2),a*sin(pi/2-alpha)*sin(pi-alpha/2),a*cos(pi/2-alpha))
    rprim2 = (a*sin(alpha)*cos(3*pi/2-alpha),a*sin(alpha)*sin(3*pi/2-alpha),a*cos(alpha))
    rprim3 = (0.0,a*sin(alpha),a*cos(alpha))

    print '\nUsing a lattice constant of %f bohr and angle of %f degrees,' % lesar_data[P]
    print 'the primitive vectors (in bohr) for the R3c rhombohedral lattice are:\n'
    print '%f\t%f\t%f' % rprim1
    print '%f\t%f\t%f' % rprim2
    print '%f\t%f\t%f\n' % rprim3

if __name__ == '__main__':
    main()
