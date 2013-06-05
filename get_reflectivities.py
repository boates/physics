#!/usr/bin/env python
"""
Get the w=0 reflectivities for all densities @ a given T
and their errorbars via a standard deviation calculation
"""
import os, sys, commands, glob

def main():

    # Retrieve user input
    try:
        wavelength = float(sys.argv[1])
    except:
        print '\n usage: '+sys.argv[0]+' wavelength(nm) \n'
        sys.exit(0)

    # Locate density directories
    rss_rough = glob.glob('1.*')
    rss = [rs for rs in rss_rough if len(glob.glob(rs+'/abs*.dat')) != 0]
    rss.sort()

#    exclude = ['1.31','1.30','1.29','1.28','1.26']
    exclude = []

    # Determine what line to use based on given wavelength
    hbar =  6.5822E-16     # eV * s
    c    =  2.997925E+08   # m / s
    E    =  hbar*c/(wavelength*1E-09)
    energies = commands.getoutput("awk '{print $1}' "+glob.glob(rss[0]+'/abs*.dat')[0]).split()
    dE = []
    for i in range(len(energies)):
        dE.append( abs(float(energies[i]) - E) )
    index = dE.index(min(dE))
    
    # Parse for values and write to output
    out = open('reflectivity.dat','w')
    out.write('# P(GPa), 0 < R < 1, stdev (using lambda = '+str(hbar*c/float(energies[index])*1E9)+' nm) \n')
    cwd = os.getcwd()
    
    for rs in rss:
        if rs not in exclude:
            P = commands.getoutput("grep '"+rs.replace('.','\.')+"' rs_pressure.dat | awk '{print $2}'")
            os.chdir(rs)
            os.system('avg_abs.py abs')
            # index + 2 accounts for the header in the avg file and shifting python zero indexing
            ref = commands.getoutput("head -"+str(int(index+2))+" avg_abs.dat | tail -n-1 | awk '{print $2}'")
            sd  = commands.getoutput("head -"+str(int(index+2))+" avg_abs.dat | tail -n-1 | awk '{print $3}'")
            os.chdir(cwd)
            out.write(str(P)+' '+str(ref)+' '+str(sd)+'\n')
    out.close()


if __name__ == '__main__':
    main()
