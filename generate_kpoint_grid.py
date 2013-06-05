#!/usr/bin/env python
"""
Create even grid of kpoint coordinates.
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        x_trans = int(sys.argv[1])
        y_trans = int(sys.argv[2])
        z_trans = int(sys.argv[3])
    except:
        print '\n usage: '+sys.argv[0]+' xtrans, y_trans, z_trans\n'
        sys.exit(0)

    # Get kpoint stepsize
#    x_step = 0.5/x_trans
#    y_step = 0.5/y_trans
#    z_step = 0.5/z_trans
    x_step = 1.0/x_trans
    y_step = 1.0/y_trans
    z_step = 1.0/z_trans
    
    weight = 1.0/(x_trans*y_trans*z_trans)

    # Open and write output file
    out = open('kpoints.dat','w')
    for i in range(x_trans):
        for j in range(y_trans):
            for k in range(z_trans):
                out.write(str(round(i*x_step,3))+' '+str(round(j*y_step,3))+' '+str(round(k*z_step,3))+' '+str(weight)+'\n')
    out.close()
    

if __name__ == '__main__':
    main()
