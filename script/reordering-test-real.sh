#!/bin/bash
​  
for f in "../"*
do
  echo $(basename $f)
done
