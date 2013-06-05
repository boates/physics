#!/usr/bin/env python

# Script to automate the submission of several jobs
# with varying natom to check the effects on enthalpy/atom

import os, sys, commands
import math

global RS
global NATOM_i
global NATOM_f
global dN
RS = 3.0575351449264874  # "3.06"
NATOM_i = 432
NATOM_f = 108
dN = 1

def calculate_box_from_rs(rs,natom):
        """
        Calculate the lattice parameter of a cell based on
        given values of rs & natom. (cubic structure is assumed)
	reutrns acell in angstroms
        """
        A = 0.5291772108  # Bohr to ang conversion
        pi = math.pi

        volume = (1.0/3)*(4*int(natom)*pi)*(float(rs))**3 # in bohr^3
        length = volume**(1.0/3) # in bohr
        acell = length * A

        return acell

def main():
    """
    Submit many dl_poly jobs with different natoms.
    """
    natom = NATOM_i

    for i in range((NATOM_i - NATOM_f) / dN + 1):

        dirname = str(natom)

        os.system('mkdir '+dirname)
	os.system('cp CONTROL TABLE dl.par '+dirname+'/')
        os.system('sed s/NATOM/'+str(natom)+'/g FIELD.template > '+dirname+'/FIELD')
	
	acell = calculate_box_from_rs(RS,natom)
	acell = str( round(acell,8) )[:11]
	
	n_CONFIG_cur_lines = 5 + 2*natom
	
        os.system('sed s/ACELL/'+str(acell)+'/g CONFIG.template > '+dirname+'/CONFIG.tmp')
        os.system('head -'+str(n_CONFIG_cur_lines)+' '+dirname+'/CONFIG.tmp > '+dirname+'/CONFIG')
	os.system('rm '+dirname+'/CONFIG.tmp')

        cwd = os.getcwd()
	os.chdir(dirname+'/')
	os.system('qsub dl.par')
	os.chdir(cwd)

        natom -= dN

if __name__ == '__main__':
    main()
