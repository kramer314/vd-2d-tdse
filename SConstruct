# Copyright (c) 2016 Alex Kramer <kramer.alex.kramer@gmail.com>
# See the LICENSE.txt file at the top-level directory of this distribution.

import os

# Uncomment to have SCons handle updating and building dependencies as well.
# if GetOption("clean"):
#     os.system("./make_deps.sh --clean")
# else:
#     n_jobs = GetOption("num_jobs")
#     os.system("./make_deps.sh -j " + str(n_jobs))

def absdir(path):
    return os.path.abspath(path) + "/"

root_dir = absdir("./")
source_dir = root_dir + "src/"
build_dir = root_dir + "build/"
deps_dir = root_dir + "deps/"

libs = ["flib", "tridiag", "wfmath"]
lib_path = [deps_dir + "fortran-lib/build/", deps_dir + "tridiag/build", deps_dir + "wfmath/build"]

env = DefaultEnvironment(ENV = os.environ, TOOLS = ['default', "gfortran"])

IEEE_flags = "-fno-unsafe-math-optimizations -frounding-math -fsignaling-nans "
general_flags = "-frecursive "
openmp_flags = "-fopenmp "
debug_flags = "-Og -g3 -Wall -Wextra -Wconversion -Wunused-parameter " + \
    "-pedantic -std=f2008 -fcheck=all -fbacktrace "
prod_flags = "-O3 -march=native -fexternal-blas -lblas "

flags = general_flags + IEEE_flags + openmp_flags + debug_flags
flags = general_flags + IEEE_flags + openmp_flags + prod_flags

env.Replace(F90FLAGS = flags)
env.Replace(LINKFLAGS = flags)
env.Replace(FORTRANMODDIRPREFIX = "-J ")
env.Replace(FORTRANMODDIR = build_dir)
env.Replace(F90PATH = [lib_path, build_dir])

Export("env")
Export("libs")
Export("lib_path")

SConscript(source_dir+"SConscript", variant_dir=build_dir, duplicate=1)

# For whatever reason, we can't use duplicate=0 and have *.mod files in the
# build directory. But, if we duplicate the source tree into the build
# directory SCons doesn't automatically clean the source files, so we have to
# manually define the entire build directory as a cleaning target.
Clean(".", build_dir)
