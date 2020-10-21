#!/bin/bash
​  
for file in "../"*
do
  printf "%s\n" "$file" | cut -d"/" -f8
done
