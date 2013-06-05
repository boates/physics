#!/bin/bash

t=$1
rm -f C_$t.dat
list=" \
01 \
02 \
03 \
04 \
05 \
06 \
07 \
08 \
09 \
10 \
11 \
12 \
13 \
14 \
15 \
16"

for i in $list; do
    c=`head -$t c$i.dat | tail -n-1 | awk '{print $2}'`
    echo $i $c >> C_$t.dat
done

exit 0
