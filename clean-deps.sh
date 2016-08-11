#!/bin/bash

devnull=/dev/null

flibdir=./deps/fortran-lib
tridiagdir=./deps/tridiag

mkdir -p ./deps

pushd $flibdir > $devnull
scons --clean
popd > $devnull

pushd $tridiagdir > $devnull
scons --clean
popd > $devnull
