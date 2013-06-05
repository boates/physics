#!/usr/bin/env python
"""
Script to run force-matching  \\ and table writing code automatically
                              \\  - (not currently)

Must execute in output/ directory, i.e. the one with the OUTCAR.* files
"""
import os, sys, commands, glob

def main():

        # Grab OUTCAR's & skip OUTCAR.100 always
	outs = glob.glob('OUTCAR.*')
	outs.sort()
	outs.pop(0)

        # Use last three OUTCAR's
	if len(outs) <= 3: os.system('cat OUTCAR.* > OUTCAR')
	else: os.system('cat '+outs[-3]+' '+outs[-2]+' '+outs[-1]+' > OUTCAR')

        os.system('force_match_vasp.py OUTCAR > fm.log')
	
        os.system('mv -f Input.in fm.log fmv.log FORCE_OUT/')
        os.system('rm -f TRAJECTORYn1')

	os.system('mv FORCE_OUT ../analysis/')

	os.system('rm -rf OUTCAR')


####################################
#        os.system('cp /home/boates/bin/INPUT_GROM.template ./')

        # Get needed parameters
#        rcut = float(commands.getoutput("tail -1 FORCE_OUT/Poten_NN.dat | awk '{print $1}'"))*0.529177
#        rcut = str( round(int(rcut*100)/100.,2) )
#        poten = commands.getoutput("awk '{print $2}' FORCE_OUT/Poten_NN.dat").split()
#        for i in range(poten.count(poten[0])):
#            poten.pop(0)
#        core_poten = poten[4]
#        core = str(round(float(commands.getoutput("grep "+core_poten+" FORCE_OUT/Poten_NN.dat | awk '{print $1}'")),2))

#        # Prepare input files and run table writing code
#        os.system('sed -e s/CORERADIUSINBOHR/'+core[0:4]+'/g INPUT_GROM.template > INPUT_GROM.tmp')
#        os.system('sed -e s/CUTOFFINANGSTROMS/'+rcut[0:4]+'/g INPUT_GROM.tmp > INPUT_GROM')
#        os.system('f_table_writer.e > tw.log')
#        os.chdir(pwd)


if __name__ == '__main__':
    main()
