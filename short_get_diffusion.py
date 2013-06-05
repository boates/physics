#!/usr/bin/env python
"""
Grab Diffusion constants from the two different methods
"""
import os, sys, commands, glob

# RS is a list of all the names of the rs directories
global RS
RS = commands.getoutput('ls -1 | grep 1 | grep -v c').split()

def main():

    # Retrieve user input
    try:
        NLINES1 = sys.argv[1]
        NLINES2 = sys.argv[2]
    except:
        print '\n usage: '+sys.argv[0]+' NLINES_D NLINES_DIFFCOEF\n'
        sys.exit(0)

    # Open the output file
    out = open('diffusion.dat','w')
    out.write('# rs, P(GPa), diffusion(cm^2/s) from MSD-avg, from VACF-avg, min, max, NLINES1='+NLINES1+', NLINES2='+NLINES2+'\n')

    for rs in RS:
        
        # Get pressure and cell information and assume cubic cell
        try:
            P = commands.getoutput("tail -2 "+rs+"/short_analysis/pressure.blocker | head -1 | awk '{print $4}'")
        except:
            P = '----'

        try:
            os.system("tail -"+NLINES1+" "+rs+"/short_analysis/D.dat | awk '{print $2}' > D.tmp")
            d1 = commands.getoutput("blocker D.tmp | tail -2 | head -1 | awk '{print $4}'")
            os.system("tail -"+NLINES2+" "+rs+"/short_analysis/DiffCoef.dat | awk '{print $2}' > DiffCoef.tmp")
            d2 = commands.getoutput("blocker DiffCoef.tmp | tail -2 | head -1 | awk '{print $4}'")
            f = open('DiffCoef.tmp','r')
            lines = f.readlines()
            f.close()
            data = []
            for line in lines:
                data.append(float(line))
            dmax = max(data)
            dmin = min(data)
            diff = max(data) - min(data)
            os.system('rm -f D.tmp DiffCoef.tmp')
        except:
            d1, d2, dmax, dmin, diff = '----', '----', '----', '----', '----'
            
        # Write to the output file
        if (P or d1 or d2) == '----':
            out.write('#')
        out.write(rs+'  '+P+'  '+d1+'  '+d2+'  '+str(dmin)+'  '+str(dmax)+'\n')

    out.close()


if __name__ == '__main__':
    main()
