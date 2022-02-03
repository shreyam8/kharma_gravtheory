#!/bin/bash

NAME=$1
VAR=$2
mkdir -p frames_${NAME}_${VAR}
cd frames_${NAME}_${VAR}

SCRIPT_DIR="$(dirname $0)"
parallel -P 8 python $SCRIPT_DIR/quick_plot.py {} $SCRIPT_DIR/../pars/${NAME}.par $VAR frame_{#} ::: ../${NAME}.*.phdf
