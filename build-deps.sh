#!/bin/bash

devnull=/dev/null

flibdir=./deps/fortran-lib
tridiagdir=./deps/tridiag

mkdir -p ./deps

# Pull and build generic Fortran helper library
git clone https://github.com/kramer314/fortran-lib.git $flibdir
pushd $flibdir > $devnull
scons
popd > $devnull

# Pull and build tridiagonal matrix library
git clone https://github.com/kramer314/tridiag.git $tridiagdir
pushd $tridiagdir > $devnull
scons
popd > $devnull
