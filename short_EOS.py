#!/usr/bin/env python
"""
Extract EOS curve for present densities.
"""
import os, sys, commands, glob

# RS is a list of all the names of the rs directories
global RS

def main():

    RS = commands.getoutput('ls -1 | grep "1\." | grep -v c').split()
    RS += commands.getoutput('ls -1 | grep "2\." | grep -v c').split()

    # Open the output file
    out = open('EOS.dat','w')
    out.write('# rs, P (GPa), dP (GPa)\n')

    for rs in RS:
        
        # Get pressure and cell information and assume cubic cell
        try:
            P, dP = commands.getoutput("tail -2 "+rs+"/short_analysis/pressure.blocker | head -1 | awk '{print $4,$6}'").split()
        except:
            P, dP = '--------', '--------'

        # Write to the output file
        if (P or dP) == '--------':
            out.write('#')
        out.write(rs+'   '+P+'   '+dP+'\n')

    out.close()


if __name__ == '__main__':
    main()
