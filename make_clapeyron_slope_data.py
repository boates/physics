#!/usr/bin/env python
"""
Give slope, P, and T, and it will make the 2-point data file of a
line on that point of length 250 K/Gpa.
"""
import sys

def main():

    # Retrieve user input
    try:
        slope = float(sys.argv[1])
        P = float(sys.argv[2])
        T = float(sys.argv[3])
    except:
        print '\n usage: '+sys.argv[0]+'  dT/dP  P(GPa)  T(K)\n'
        sys.exit(0)

    # Create output file
    dP = 200.0/slope
    out = open('cc_'+str(int(T))+'K_new.dat','w')
    out.write(str(P-dP)+' '+str(T-dP*slope)+'\n')
    out.write(str(P+dP)+' '+str(T+dP*slope)+'\n')


if __name__ == '__main__':
    main()
