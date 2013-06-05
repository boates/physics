#!/bin/bash

t=$1
rm -f CO$t.dat
list=" \
1 \
2 \
3 \
4 \
5 \
6 \
"

for i in $list; do
    c=`head -$t CO$i.dat | tail -n-1 | awk '{print $2}'`
    echo $i $c >> C_$t.dat
done

exit 0
