#!/bin/bash

set -e # Abort on error.

git submodule init
git submodule update

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_GENERIC=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX -DUSE_JEMALLOC=OFF

make -j ${CPU_COUNT} ${VERBOSE_AT}
make install
