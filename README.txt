Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
See the LICENSE.txt file at the top-level of this distribution.

THIS IS BETA SOFTWARE. FEATURES and FUNCTIONALITY ARE SUBJECT TO CHANGE WITHOUT
PRIOR NOTICE.

Description
===========
Implementation of a multi-dimensional flux-based method for extracting momentum
distributions from traveling wavepackets in time-dependent quantum mechanical
calculations. Both a semi-classical method and a method incorporating
first-order quantum corrections are implemented.

Example usage
=============
The following will run a simulation with the parameters given in the default
input input file `run.inp`:

    ./run.sh input/run.inp

After execution, the output directory in the input file (by default `output/`)
will have five updated files present:

* run.inp: Copy of input file
* run.log: Runtime logging
* vd_p.dat: Extracted momentum spectrum

[TODO]
* Time-dependent output of wavefunction
* Normalization checks

Required dependencies
=====================
Software to install:

* Fortran compiler (most recently tested on GNU Fortran 6.1)
* SCons (most recently tested on SCons 2.5)

External libraries to build from source (build scripts included):

* Fortran generics library, available at:
  https://github.com/kramer314/fortran-lib
  Note that this library is still in beta status, so the API is technically
  still in flux. Always use the most recent version.

* Fortran tridiagonal matrix library, available at:
  https://github.com/kramer314/tridiag
  Note that this library is still in beta status; always use the latest
  version.

* Wavefunction math library, available at:
  https://github.com/kramer314/wfmath
  Note that this library is still in beta status; always use the latest
  version.

Recommended Dependencies
========================
* OpenMP library (most recently tested on GCC OpenMP 5.3)

Basic setup
===========
The following will download and build this program and its library
dependencies:

    git clone https://github.com/kramer314/vd-2d-tdse.git
    cd ./vd-2d-tdse/
    ./make-deps.sh -j [N]
    scons -j [N]

where ``[N]`` is the number of build threads (omit the ``-j`` flag to build
everything in serial.

Manual setup
============
Distribution of Fortran libraries is pretty nasty until submodules from F2015
are supported. What the `./build-deps.sh` script does is clone the Fortran
generics, tridiagonal matrix, and wavefunction math library projects that this
program uses into the directory `deps/` and then builds it into the folder
`deps/[lib-name]/build/`. You can have this generics library project at some
other location; you just need to build the library as normal and create
symlinks as follows:

    ln -s /path/to/fortran-lib-project ./deps/fortran-lib
    ln -s /path/to/tridiag-matrix-project ./deps/tridiag
    ln -s /path/to/wfmath ./deps/wfmath

Note that this is a link to the project directory itself, *not* the directory
where the build folder where the compiled library is actually located. These
symlinks are necessary because the build script for this program automatically
adds the paths x`./deps/[lib-name]/buiild` to the compiler path for linking.
