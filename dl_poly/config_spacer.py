#!/usr/bin/env python

# Fix crappily formatted CONFIG files (hopefully)

import os, sys, commands

t = open('CONFIG','r')
lines = t.readlines()
t.close()
t = '        '

out = open('CONFIG_fixed','w')
out.write(lines.pop(0)) # header
imcon, levcfg = lines.pop(0).split()
out.write(t+imcon+t+'  '+levcfg+'\n')
a1 = lines.pop(0).split()
a2 = lines.pop(0).split()
a3 = lines.pop(0).split()
#' + ' % .8e' % SN'
out.write(t+' % .8f' % float(a1[0])+t+' % .8f' % float(a1[1])+t+' % .8f' % float(a1[2])+'\n')
out.write(t+' % .8f' % float(a2[0])+t+' % .8f' % float(a2[1])+t+' % .8f' % float(a2[2])+'\n')
out.write(t+' % .8f' % float(a3[0])+t+' % .8f' % float(a3[1])+t+' % .8f' % float(a3[2])+'\n')

for line in lines:
    if lines.index(line) % 2 != 0:
        row = line.split()
        for i in range(len(row)):
            if float(row[i]) >= 0:
                out.write('       '+'% .8f' % float(row[i]))
            elif float(row[i]) < 0:
                out.write('      '+'% .8f' % float(row[i]))
        out.write('\n')
    else:
        row = line.split()
        for entry in row:
            out.write(entry)
            if entry == row[0]:
                out.write(t)
            if entry != row[-1]:
                out.write(t)
            else:
                out.write('\n')
out.close()
