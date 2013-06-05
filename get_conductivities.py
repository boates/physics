#!/usr/bin/env python
"""
Get the w=0 conductivities for all densities @ a given T
and their errorbars via a standard deviation calculation
"""
import os, sys, commands, glob

def main():

    # Locate density directories
    rss_rough = glob.glob('1.*')
    rss = [rs for rs in rss_rough if len(glob.glob(rs+'/sig*.dat')) != 0]
    rss.sort()

#    exclude = ['1.31','1.30','1.29','1.28','1.26']
    exclude = []
    
    # Parse for values and write to output
    out = open('conductivity.dat','w')
    out.write('# P(GPa), cond(ohm*cm)^-1, stdev\n')
    cwd = os.getcwd()
    for rs in rss:
        if rs not in exclude:
            P = commands.getoutput("grep '"+rs.replace('.','\.')+" ' rs_pressure.dat | awk '{print $2}'")
            os.chdir(rs)
            os.system('avg_sig.py sig')
            sig = commands.getoutput("head -2 avg_sig.dat | tail -n-1 | awk '{print $2}'")
            sd  = commands.getoutput("head -2 avg_sig.dat | tail -n-1 | awk '{print $3}'")
#            print rs, P, sig, sd
            os.chdir(cwd)
            out.write(str(P)+' '+str(sig)+' '+str(sd)+'\n')
    out.close()


if __name__ == '__main__':
    main()
