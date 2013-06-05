#!/bin/bash
awk '{print $1/1E-12*(8.0655E+03/2.41796E+14),$2*1E-12/(8.0655E+03/2.41796E+14)}' $1 > vdos.dat
mv -f vdos.dat $1
