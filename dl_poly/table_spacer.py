#!/usr/bin/env python

# Fix crappily formatted TABLE files (hopefully)

import os, sys, commands

t = open('TABLE.orig','r')
lines = t.readlines()
t.close()

out = open('TABLE_fixed','w')
out.write(lines.pop(0).strip()+'\n')
out.write('  '+lines.pop(0).strip()+'\n')
out.write(lines.pop(0).strip()+'\n')

for line in lines:
    row = line.split()
    for i in range(len(row)):
        if float(row[i]) >= 0:
            out.write('  '+row[i][:8]+row[i][-4:])
        elif float(row[i]) < 0:
            out.write(' '+row[i][:9]+row[i][-4:])
    out.write('\n')

out.close()
