#!/bin/bash

file=$1
sed "s/GLU/GLH/g" < $file | sed "s/ASP/ASH/g" |  sed "s/LYS/LYN/g" | sed "s/ARG/ARN/g" > netural-$file

