#!/bin/sh
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1

## User-defined parameters
# k: number of working kernels (default: 8)
# i: number of iterations (default: 5)
# r: number of repeats (default: 2)
# p: adjusted p-value threshold (default: 0.02)
while getopts ":k:i:r:p:" opt
do
   case "$opt" in
      k ) parameterK=$OPTARG ;;
      i ) parameterA=$OPTARG ;;
      r ) parameterB=$OPTARG ;;
      p ) parameterC=$OPTARG ;;
   esac
done

if [ ! -n "$parameterK" ]; then
  parameterK=8
  printf "kernels: (default) $parameterK\n"
else
  printf "kernels: $parameterK\n"
fi

if [ ! -n "$parameterA" ]; then
  parameterA=10
  printf "iterations: (default) $parameterA\n"
else
  printf "iterations: $parameterA\n"
fi

if [ ! -n "$parameterB" ]; then
  parameterB=2
  printf "repeats: (default) $parameterB\n"
else
  printf "repeats: $parameterB\n"
fi

if [ ! -n "$parameterC" ]; then
  parameterC=1e-3
  printf "adjusted p-value threshold: (default) $parameterC\n"
else
  printf "adjusted p-value threshold: $parameterC\n"
fi

# Run
python causal_regions_variants.py -kernels $parameterK -iters $parameterA -repeats $parameterB
Rscript causal_egenes.R $parameterC