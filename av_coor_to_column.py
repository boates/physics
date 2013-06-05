#!/usr/bin/env python
"""
Convert an av_coor... file from bonev's output of get_coord.pl to a column
"""
import sys

def main():

    try:
        f = open(sys.argv[1],'r')
    except IndexError:
        print '\nusage: '+sys.argv[0]+' av_coordination_rc_###.dat\n'
        sys.exit(0)

    # Read off the header
    for i in range(11): f.readline()

    # Read the average coordination information
    crd = f.readline().split()

    out = open('coordination.dat','w')

    for i in range(len(crd)):
        out.write(str(i)+'    '+crd[i]+'\n')

    out.close()

if __name__ == '__main__':
    main()
