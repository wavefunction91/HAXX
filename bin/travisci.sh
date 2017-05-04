#!/bin/sh

set -e

if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
else
    # no OpenMP support in clang, will use C++11 threads
    export CC=/usr/bin/clang-3.8
    export CXX=/usr/bin/clang++-3.8
fi
export CXXFLAGS="-std=c++11 -O3"
export FC=gfortran


CMAKE_URL="https://cmake.org/files/v3.7/cmake-3.7.2-Linux-x86_64.tar.gz"

wget --no-check-certificate --quiet -O - $CMAKE_URL | tar --strip-components=1 -xz -C cmake
PATH=$PWD/cmake/bin":"$PATH


mkdir build && cd build
cmake --version
cmake -DBOOST_LIBRARYDIR='/usr/lib' -DCMAKE_Fortran_FLAGS='-O3' ..
make -j2
make test
