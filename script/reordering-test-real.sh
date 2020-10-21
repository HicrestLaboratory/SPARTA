#!/bin/bash

for file in "data/"*
do
  printf "%s\n" "$file" | cut -d"/" -f8
done
