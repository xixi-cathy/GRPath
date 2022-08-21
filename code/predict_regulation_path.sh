#!/bin/sh

## User-defined parameter
# r: number of repeats (default: 3)
while getopts ":f:r:" opt
do
   case "$opt" in
      f ) parameterF="$OPTARG" ;;
      r ) parameterR=$OPTARG ;;
   esac
done

if [ ! -n "$parameterF" ]; then
  parameterF="demo_data"
  printf "data folder: (default) $parameterF\n"
else
  printf "data folder: $parameterF\n"
fi

if [ ! -n "$parameterR" ]; then
  parameterR=3
  printf "repeats: (default) $parameterR\n"
else
  printf "repeats: $parameterR\n"
fi

## Run
python link_cell_type.py -datafolder $parameterF -repeats $parameterR
Rscript regulation_path.R