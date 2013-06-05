#!/usr/bin/env python
"""
Create a file of rs, P, avg_dissociation
"""
import os, sys, commands, glob

# RS is a list of all the names of the rs directories
global RS
RS = commands.getoutput('ls -1 | grep 1 | grep -v c').split()

def main():

    # Open the output file
    out1 = open('rs_diss.dat','w')
    out2 = open('P_diss.dat','w')
    out1.write('# rs, avg_dissociation, N_diss\n')
    out2.write('# P (GPa), avg_dissociation, N_diss\n')

    natom = commands.getoutput("head */output/XDATCAR.100 | head -2 | tail -1 | awk '{print $1}'")

    for rs in RS:
        
        # Get pressure and cell information and assume cubic cell
        try:
            P = commands.getoutput("tail -2 "+rs+"/analysis/pressure.blocker | head -1 | awk '{print $4,$6}'").split()[0]
        except:
            P = '--------'

        try:
            d = commands.getoutput("blocker "+rs+"/analysis/dissociation.dat 2 | tail -2 | head -1 | awk '{print $4}'")
            N_diss = str( int( (int(natom)*float(d)) + 0.5 ) )
        except:
            d, N_diss = '--------', '--------'

        # Write to the output file
        if '----' in d:
            out1.write('#')
            out2.write('#')
        if 'tail' in P: out2.write('#')
        out1.write(rs+'   '+d+'   '+N_diss+'\n')
        out2.write(P+'   '+d+'   '+N_diss+'\n')

    out1.close()
    out2.close()


if __name__ == '__main__':
    main()
