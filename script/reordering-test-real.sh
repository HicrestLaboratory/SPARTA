#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
DATADIR="$(DIR)/../data"

for f in "$(DATADIR)"/*
do
  echo $(basename $f)
done
