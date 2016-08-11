#!/bin/bash

export OMP_NUM_THREADS=4

# $1 is the input file

./build/main.x $1
