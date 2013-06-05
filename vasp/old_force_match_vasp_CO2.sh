#!/bin/tcsh

rm -f position force position_force TRAJECTORYn1
rm -rf FORCE_OUT

set acell=`grep -A1 'length of vectors' OUTCAR | head -2 | tail -1 | awk '{print $1/0.529177}'`
echo "acell: " $acell

set natoms=`grep 'NIONS =' OUTCAR | head -1 | awk '{print $12}'`
echo "natoms: " $natoms

set natomp1=`echo "scale=10; $natoms+1" | bc`

set half_box=`echo "scale=10; $acell/2.0" | bc`
echo "half_box: " $half_box

grep -A$natomp1 "POSITION                                       TOTAL-FORCE (eV/Angst)" OUTCAR | grep -v "\-\-" | grep -v POSITION | awk '{print $1/0.529177,$2/0.529177,$3/0.529177,$4*0.019446906,$5*0.019446906,$6*0.019446906}' > position_force

set number_of_lines=`wc -l position_force | awk '{print $1/'$natoms'}'`
echo "number_of_Lines: " $number_of_lines

awk '{if (n % '$natoms' == 0) {print n/'$natoms' + 1;print $0} else print $0;n+=1}' position_force > TRAJECTORYn1
set linecount=`wc -l TRAJECTORYn1 | awk '{print $1}'`
echo "linecount: " $linecount
set denominator=`expr $natoms + 1`
set number_of_snapshots=`expr $linecount / $denominator`
echo "number_of_snapshots: " $number_of_snapshots
set number_to_append=`expr $number_of_snapshots + 1`
echo $number_to_append >> TRAJECTORYn1

rm position_force

sed s/acell/"$acell"/g /home/boates/bin/TEMPLATE_CO2.in | sed s/number_of_snapshots/"$number_of_snapshots"/g | sed s/half_box/"$half_box"/g | sed s/natoms/"$natoms"/g > Input.in

/home/boates/bin/m.e
