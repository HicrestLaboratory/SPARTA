#!/bin/bash

for f in "data/"*
do
  echo $(basename $f)
done
