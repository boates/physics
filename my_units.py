#!/usr/bin/env python

# VERSION 1.6
# Author: Brian "Completely Awesome" Boates

# Useful unit conversion tool
# $./units.py value unit_a unit_b
# The above command will attempt to convert the 
# given value from unit_a to unit_b

import os, sys, commands

def convert(a,b,value,a_to_b,factors,old_a=None,old_value=None,results=False):
    """
    Tries the quick and easy one step conversion from unit a to b.
    """
    if old_a:
        unit_1 = old_a
    if not old_a:
        unit_1 = a
    if (a,b) in a_to_b:
        results = str(value)+' '+unit_1+' = '+str(value*factors[a_to_b.index((a,b))])+' '+b
    elif (b,a) in a_to_b:
        results = str(value)+' '+unit_1+' = '+str(value/factors[a_to_b.index((b,a))])+' '+b
    if results and old_a:
        results = results.replace(str(value),str(old_value))
        
    return results

def converter(a,b,value,a_to_b,factors,results=False,old_a=None,old_value=None):
    """
    Sort through the UNITS file and perform conversion though more than
    one conversion if necessary to achieve desired conversion.
    """
    results = convert(a,b,value,a_to_b,factors,old_a=old_a,old_value=old_value)    

    if not results:
        possibilities = [pair[1] for pair in a_to_b if a == pair[0]] + \
                        [pair[0] for pair in a_to_b if a == pair[1]]
        a_to_c = []
        for p in possibilities:
            a_to_c.append((a,p))

        cfactors = []
        for p in a_to_c:
            try:
                cfactors.append(factors[a_to_b.index(p)])
            except ValueError:
                cfactors.append(1./factors[a_to_b.index((p[1],p[0]))])

        for p in a_to_c:
            tmp_results = convert(p[0],p[1],value,a_to_c,cfactors)
            if tmp_results:
                input = tmp_results.split()
                new_value = float(input[3])/float(input[0])
                new_a = input[-1]
                if a == first_a and value == first_value:
                    results = convert(new_a,b,value*new_value,a_to_b,factors,old_a=a,old_value=value)
                else:
                    results = convert(new_a,b,value*new_value,a_to_b,factors,old_a=first_a,old_value=first_value)
                if results:
                    break
        
    if not results:
        try:
            results = converter(new_a,b,value*new_value,a_to_b,factors,old_a=first_a,old_value=first_value)
        except RuntimeError:
            print 'Desired conversion not found in UNITS file'

    return results


  #==============#
  # MAIN PROGRAM #
  #==============#

##### Read in unit information from the UNITS file and organize into proper variables.
##### Run the above function "converter" which attempts to execute the function
##### "convert" immediately but will perform further measures if not successful.

UNITS_PATH = commands.getoutput('echo $HOME')+'/util/'
f = open(UNITS_PATH+'UNITS','r')
header = f.readline()

factors = []
a_to_b = []
check = []
for conversion in f.readlines():
    row = conversion.split()
    if len(row) == 5:
        factors.append(float(row[-2]))
        a_to_b.append((row[1],row[-1]))
        check.append(row[1])
        check.append(row[-1])
    
try:
    value = float(sys.argv[1])
    a = str(sys.argv[2])
    b = str(sys.argv[3])
except IndexError:
    print 'Please provide necessary arguments'
    print 'i.e. "$ ./units.py value unit_a unit_b" would convert the value from unit_a to unit_b'    
    sys.exit(0)

if a not in check and b not in check:
    print 'The units "'+a+'" and "'+b+'" were both not found in the UNITS file'
    sys.exit(0)
if a not in check:
    print 'The unit "'+a+'" was not be found in the UNITS file'
    sys.exit(0)
if b not in check:
    print 'The unit "'+b+'" was not be found in the UNITS file'
    sys.exit(0)

global first_a
global first_value

first_a = a
first_value = value

results = converter(a,b,value,a_to_b,factors)
if results:
    print results
