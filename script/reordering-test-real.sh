#!/bin/bash

echo "The script you are running has basename `basename "$0"`, dirname `dirname "$0"`"

for f in "../"*
do
  echo $(basename $f)
done
