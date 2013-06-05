#!/bin/bash

### Script to perform gradual melting simulations with DL_POLY
### Make sure that the correct TABLE, FIELD and CONFIG files
### as well as CONTROL-orig, and CONTROL-restart files and
### CONTROL template files are present.

temp=400
temp_STEP=10
temp_LIMIT=600

sed s/DL_TEMP/$temp.0/g CONTROL-orig.template > CONTROL-orig
cp CONTROL-orig CONTROL
qsub -sync y dl.par
mkdir 0$temp
cp CONTROL CONFIG TABLE FIELD REV* 0$temp/
mv HISTORY STATIS OUTPUT RDFDAT dl_poly.* 0$temp
$HOME/software/dl_poly/copy
temp=`expr $temp + $temp_STEP`

while [ "$temp" -le "$temp_LIMIT" ]
  do
    
  sed s/DL_TEMP/$temp.0/g CONTROL-restart.template > CONTROL-restart
  cp CONTROL-restart CONTROL
  qsub -sync y dl.par
  mkdir 0$temp
  cp CONTROL CONFIG TABLE FIELD REV* 0$temp/
  mv HISTORY STATIS OUTPUT RDFDAT dl_poly.* 0$temp
  $HOME/software/dl_poly/copy
  temp=`expr $temp + $temp_STEP`

  done
