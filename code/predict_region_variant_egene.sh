#!/bin/sh
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
python causal_regions_variants.py
Rscript causal_egenes.R