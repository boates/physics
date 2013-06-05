#!/usr/bin/env python
"""
Extract average correlator info, for lifetimes of 1/5 a picosecond (assuming 32au tsteps)
"""
import os, sys, commands, glob

# RS is a list of all the names of the rs directories
global RS
RS = commands.getoutput('ls -1 | grep 1 | grep -v c').split()

def main():

    try:
        c = int(sys.argv[1])
    except IndexError:
        print '\nusage: '+sys.argv[0]+'  N_correlator_neighbors (i.e. 1,2,3,etc...)\n'
        sys.exit(0)

    if c < 10:
        c_str = '0'+str(c)
    else:
        c_str = str(c)
        
    # Open the output file
    out = open('c_'+c_str+'.dat','w')
    out.write('# rs, P (GPa), dP (GPa), correlator_at_0.2 picoseconds\n')

    for rs in RS:
        
        # Get pressure and correlator at 280th step (10x 14fs N2 oscillations (for tstep=0.50fs))
        try:
            P = float(commands.getoutput("grep Final "+rs+"/analysis/pressure.blocker | awk '{print $4}'"))
        except:
            P = '--------'
        try:
            f = float(commands.getoutput("tail -280 "+rs+"/analysis/c"+c_str+".dat | tail -n-1 | awk '{print $2}'"))
        except:
            f = '--------'
            
        # Write to the output file
        if (P or f) == '--------':
            out.write('#')
        out.write(rs+'   '+str(P)+'   '+str(f)+'\n')

    out.close()


if __name__ == '__main__':
    main()
