#!/usr/bin/env python

"""Brief example of how to use optparse"""

import optparse

usage = 'usage: %prog [options] filename'
version="%prog 1.0"
parser = optparse.OptionParser(usage=usage,version=version)

parser.add_option('-f','--file',dest='filename',# Type defaults to str
                  help='Name of file to read in')
parser.add_option('-d', '--date', type='string', dest='date', default='', 
                  help='Date stamp YYYYMMDD of the file to retrieve')
parser.add_option('--noplot', action='store_true', dest='noplot',
                  default=False, help='Disables plotting [%default]')
parser.add_option('--secondord', action='store_true',dest='secondord',
                  default=None, help='Considers second order corrections')


options, args = parser.parse_args()
if len(args)>1:
    parser.error(usage)

if options.filename or args:
    filename = options.filename or args[0]
else:
    parser.error('A filename must be specified')

print '#================================================#'
print '# This is an example code and actually does      #'
print '# nothing!!                                      #'
print '#================================================#'
print
print 'The filename you specified is',filename
print
if options.date:
    print 'The date you passed is',options.date
    print
print 'The noplot flag is',options.noplot
if not options.secondord==None:
    print 'Considering second order corrections'
print
print
print '#================================================#'
print '# The program is complete, thanks for playing :) #'
print '#================================================#'
