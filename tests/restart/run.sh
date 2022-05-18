#!/bin/bash

# Bash script testing initialization vs restart of a torus problem

# Set paths
KHARMADIR=../..

$KHARMADIR/run.sh -i $KHARMADIR/pars/sane.par parthenon/time/nlim=5 debug/archive_parameters=false perturbation/u_jitter=0

mv torus.out0.final.phdf torus.out0.final.init.phdf

python3 $KHARMADIR/scripts/fix_restart.py torus.out1.00000.rhdf torus.out1.fixed.rhdf

$KHARMADIR/run.sh -r torus.out1.fixed.rhdf parthenon/time/nlim=5 debug/archive_parameters=false perturbation/u_jitter=0

mv torus.out0.final.phdf torus.out0.final.restart.phdf
