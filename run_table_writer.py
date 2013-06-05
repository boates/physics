#!/usr/bin/env python
"""
Script to run force-matching and table writing code automatically
"""
import os, sys, commands, glob

def main():

        # Get needed parameters
        rcut = float(commands.getoutput("tail -1 FORCE_OUT/Poten_NN.dat | awk '{print $1}'"))*0.529177
        rcut = str( round(int(rcut*100)/100.,2) )
        poten = commands.getoutput("awk '{print $2}' FORCE_OUT/Poten_NN.dat").split()
        for i in range(poten.count(poten[0])):
            poten.pop(0)
        core_poten = poten[4]
        core = str(round(float(commands.getoutput("grep "+core_poten+" FORCE_OUT/Poten_NN.dat | awk '{print $1}'")),2))

        # Prepare input files and run table writing code
        os.system('sed -e s/CORERADIUSINBOHR/'+core[0:4]+'/g ~/bin/INPUT_GROM.template > INPUT_GROM.tmp')
        os.system('sed -e s/CUTOFFINANGSTROMS/'+rcut[0:4]+'/g INPUT_GROM.tmp > INPUT_GROM')
        os.system('f_table_writer.e > tw.log')
	

if __name__ == '__main__':
    main()
