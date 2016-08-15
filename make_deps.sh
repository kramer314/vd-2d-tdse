#!/bin/bash

devnull="/dev/null"

depsdir="./deps/"
flibdir=${depsdir}"fortran-lib/"
tridiagdir=${depsdir}"tridiag/"

threads=1

function usage {
    echo "Usage: make_deps [[-c | --clean] | [-h | --help]] |"\
	 "[-j | --jobs [N] ]]"
    echo ""
    echo "--clean: removes dependency source and build subdirectories."
    echo "--help: print this usage statement"
    echo ""
    echo "--jobs: number of build threads"
    echo ""
    echo "If no arguments are specified, dependencies will be updated"\
	 "and built in serial."
}

function clean {
    echo "Removing dependency subdirectories."
    rm -rf $depsdir
}

function build {
    echo "Building dependencies."
    mkdir -p $depsdir

    # Pull and build generic Fortran helper library
    echo "Updating generic Fortran helper library source."
    echo "-----------------------------------------------"
    git clone https://github.com/kramer314/fortran-lib.git $flibdir
    pushd $flibdir > $devnull
    git pull origin master
    echo "Building generic Fortran helper library."
    echo "----------------------------------------"
    echo $threads
    scons -j${threads}
    popd > $devnull

    # Pull and build tridiagonal matrix library
    echo "Updating tridiagonal matrix library source."
    echo "-------------------------------------------"
    git clone https://github.com/kramer314/tridiag.git $tridiagdir
    pushd $tridiagdir > $devnull
    git pull origin master
    echo "Building tridiagonal matrix library."
    echo "------------------------------------"
    scons -j$threads
    popd > $devnull
}

# Parse arguments
if [ -z "$1" ]; then
    build
    exit
fi

while [ "$1" != "" ]; do
    case $1 in
	-c | --clean)
	    clean
	    exit
	    ;;
	-h | --help)
	    usage
	    exit
	    ;;
	-j )
	    shift
	    threads=$1
	    build
	    exit
	    ;;
	* )
	    usage
	    exit
    esac
    shift
done
