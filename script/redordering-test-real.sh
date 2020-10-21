#!/bin/bash
​  
for file in "../"*
do
  echo "%s\n" "$file" | cut -d"/" -f8
done
