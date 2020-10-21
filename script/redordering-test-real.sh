#!/bin/bash
​  
arr=($( ls *))

for i in $files ; do
  echo Next: $i
  let counter=$counter+1
  echo $counter