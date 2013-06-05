#!/usr/bin/env python

# Analyze all of the data from the enthalpy per atom simulations
# run_STATIS.sh to get H as column 7 of the statis file
# Use blocker to average enthalpy per atom

import os, sys, commands, glob

global NLINES
NLINES = False # Number of lines to include in blocker averaging

def run_blocker(fname,col,nlines=False):
    """
    Give the specified data fname and column number
    returns average and statistical error
    """
    if nlines:
        os.system("tail -"+str(nlines)+" "+fname+" > tail_statis.dat")
        fname = 'tail_statis.dat'
        
    block = commands.getoutput("blocker "+fname+" "+str(col)+" | grep 'Final result' | awk '{print $4,$6}'")
    avg = float(block.split()[0])
    error = float(block.split()[1]) 
    
    return avg, error

def main():
    """
    Analyze all of the data and obtain a unified result
    i.e. a file of natom vs H/atom
    """
    dirs_list = open('dirs.txt','r')
    dirs = dirs_list.readlines()
    dirs_list.close()
    out = open('enthalpy_vs_natom.dat','w')

    for dir in dirs:

        print round((float(int(dir) - int(dirs[0]))/float(int(dirs[-1]) - int(dirs[0])))*100.,2),'%'

        dir = dir.strip()
        natom = int(dir)
        cwd = os.getcwd()
        os.chdir(dir)

        if 'statis.dat' not in glob.glob('*'):
            os.system('run_STATIS.sh')
        if 'diffusion.dat' not in glob.glob('*'):
            os.system('diffusion_from_statis.py')
        if 'RDF_LI_LI.dat' not in glob.glob('*'):
            os.system('get_RDF.py all all n')
        
        H_avg, dH = run_blocker('statis.dat',7,nlines=NLINES)

        os.chdir(cwd)

        out.write(str(natom)+'    '+str(H_avg/natom)+'    '+str(dH/natom)+'\n')

    out.close()

if __name__ == '__main__':
    main()
