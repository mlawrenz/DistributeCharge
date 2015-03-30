#!/bin/bash

file=$1
base=`basename $file`
sed "s/LYS    11/LYN    11/g" < $file | sed "/ HZ1 LYN/d" > charge-${base}
