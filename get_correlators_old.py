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

    # Open the output file
    out = open('correlator_'+str(c)+'.dat','w')
    out.write('# rs, P (GPa), dP (GPa), correlator_at_0.2 picoseconds\n')

    for rs in RS:
        
        # Get pressure and cell information and assume cubic cell
        try:
            P, dP = commands.getoutput("tail -2 "+rs+"/analysis/pressure.blocker | head -1 | awk '{print $4,$6}'").split()
        except:
            P, dP = '--------', '--------'

        try:
            os.system("tail -400 "+rs+"/analysis/correlator_"+str(c)+"nn.dat | awk '{print $2}' > tmp.dat")
            f = commands.getoutput("blocker tmp.dat | tail -2 | head -1 | awk '{print $4}'")
            os.system('rm -f tmp.dat')
#            f = commands.getoutput("tail -259 "+rs+"/analysis/correlator_"+str(c)+"nn.dat | head -1 | awk '{print $2}'")
        except:
            f = '--------'
            
        # Write to the output file
        if (P or dP or f) == '--------':
            out.write('#')
        out.write(rs+'   '+P+'   '+dP+'   '+f+'\n')

    out.close()


if __name__ == '__main__':
    main()
