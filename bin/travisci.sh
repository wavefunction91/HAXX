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

mkdir build && cd build
cmake -DBOOST_LIBRARYDIR='/usr/lib' ..
make -j2
make test
