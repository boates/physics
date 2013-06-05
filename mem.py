#!/usr/bin/env python

# Check how much memory is being used by a directory
# Give directory name as first argument

import os, sys

home = os.getcwd()

try:
    dir = sys.argv[1]
except IndexError:
    print '\nPlease give the name of the directory you wish to check.'
    print '--Current working directory was used as default--\n'
    dir = home

os.chdir(dir)
os.system('du -h | sort -nu > f.tmp')
#os.chdir(home)
f = open('f.tmp','r')
lines = f.readlines()
for line in lines:
    if line.split()[-1] == '.':
        print 'Memory usage in', dir, '=', line.split()[0]+'\n'
os.remove('f.tmp')
