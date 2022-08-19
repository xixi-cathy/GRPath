#!/bin/sh

## User-defined parameter
# r: number of repeats (default: 3)
while getopts ":r:" opt
do
   case "$opt" in
      r ) parameterR=$OPTARG ;;
   esac
done

if [ ! -n "$parameterR" ]; then
  parameterR=3
  printf "repeats: (default) $parameterR\n"
else
  printf "repeats: $parameterR\n"
fi

## Run
python link_cell_type.py -repeats $parameterR
Rscript regulation_path.R