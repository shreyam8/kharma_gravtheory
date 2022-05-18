#!/bin/bash

# Bash script testing initialization vs restart of a torus problem

# Set paths
KHARMADIR=../..

$KHARMADIR/run.sh -i $KHARMADIR/pars/sane.par parthenon/time/nlim=5 debug/archive_parameters=false perturbation/u_jitter=0

mv torus.out0.final.phdf torus.out0.final.first.phdf

$KHARMADIR/run.sh -i $KHARMADIR/pars/sane.par parthenon/time/nlim=5 debug/archive_parameters=false perturbation/u_jitter=0

mv torus.out0.final.phdf torus.out0.final.second.phdf
