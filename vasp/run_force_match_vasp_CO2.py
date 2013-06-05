#!/usr/bin/env python

import os, sys, commands, glob

dirs = glob.glob('OUTCAR.*')
dirs.sort()
dirs = dirs[1:]

natom = commands.getoutput("grep 'number of ions' "+dirs[0]).split()[-1]

# Only keep the OUTCAR's that don't end in the middle of printing forces
good = []
for d in dirs:
    end = commands.getoutput("tail -n-"+natom+" "+d+" | grep 'TOTAL-FORCE'")
    if 'POSITION' not in end:
        good.append(d)

OUTCARS = ''
for g in good:
    OUTCARS += g+' '

os.system('mkdir force_match')
os.system('cat '+OUTCARS+' > force_match/OUTCAR')

os.chdir('force_match')

os.system('force_match_vasp_CO2.py OUTCAR 32 64 > fm.log &')


