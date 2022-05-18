#!/bin/bash

. ~/libs/anaconda3/etc/profile.d/conda.sh
conda activate pyharm

# Set paths
KHARMADIR=../..

python3 $KHARMADIR/scripts/compare.py torus.out0.final.first.phdf torus.out0.final.second.phdf init_vs_restart

h5diff --exclude-path Info torus.out0.final.first.phdf torus.out0.final.second.phdf
