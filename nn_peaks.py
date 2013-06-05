#!/usr/bin/env python
"""
Create a data set of nn_dist peak distances vs rs & P, for a given neighbor
"""
import os, sys, commands, glob

# RS is a list of all the names of the rs directories
global RS
RS = commands.getoutput('ls -1 | grep "1\." | grep -v c').split()

def main():

    try:
        f_nn = sys.argv[1]
        nn = int(sys.argv[2])
        scale = sys.argv[3]
    except IndexError:
        print '\nusage: '+sys.argv[0]+'  nn_average.hist  neighbor(i.e. 1,2,3,etc...) scale_with_rs?(y/n)\n'
        sys.exit(0)

    # Open the output file
    out = open('peaks_'+str(nn)+'nn.dat','w')
    out.write('# rs, P (GPa), nn='+str(nn)+' peak_distance\n')

    for rs in RS:
        
        # Get pressure
        try:
            P = commands.getoutput("tail -2 "+rs+"/analysis/pressure.blocker | head -1 | awk '{print $4}'").strip()
        except:
            P = '--------'

        # Get location of peak
        try:
            peak = float(commands.getoutput("expectation_value.py "+rs+"/analysis/"+f_nn+" 1 "+str(int(nn)+1)).split()[-1])
            if scale == 'y':
                peak /= float(rs)
        except:
            peak = '--------'
            
        # Write to the output file
        if '--' in P or '--' in str(peak):
            out.write('#')
        out.write(rs+'  '+str(peak)+'  '+P+'\n')

    out.close()


if __name__ == '__main__':
    main()
