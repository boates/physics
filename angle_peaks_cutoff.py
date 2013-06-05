#!/usr/bin/env python
"""
Create a data set of nn_dist peak distances vs rs & P, for a given neighbor
"""
import os, sys, commands, glob

# RS is a list of all the names of the rs directories
global RS
RS = commands.getoutput('ls -1 | grep "1\." | grep -v c').split()

def main():

    # Retrieve user input
    try:
        cutoff = float(sys.argv[1])
    except:
        print '\n usage: '+sys.argv[0]+' cutoff\n'
        sys.exit(0)

    # Open the output file
    out = open('angle_EV_cutoff'+str(cutoff)+'.dat','w')
    out.write('# rs, <angle>, P(GPa)\n')

    for rs in RS:
        
        # Get pressure
        try:
            P = commands.getoutput("tail -2 "+rs+"/analysis/pressure.blocker | head -1 | awk '{print $4}'").strip()
        except:
            P = '--------'

        # Get location of peak
        try:
            EV = float(commands.getoutput("expectation_value_cutoff.py "+rs+"/analysis/CO2_angles.dat 1 2 "+str(cutoff)).split()[-1])
        except:
            EV = '--------'
            
        # Write to the output file
        if '--' in P or '--' in str(EV):
            out.write('#')
        out.write(rs+'  '+str(EV)+'  '+P+'\n')

    out.close()


if __name__ == '__main__':
    main()
