#!/bin/bash

set -e # Abort on error.

git submodule init
git submodule update

CMAKE_EXTRA_ARGS=""
if [[ "$(uname)" == "Darwin" ]]; then
  CMAKE_EXTRA_ARGS="-DUSE_SYSTEM_ZLIB=ON -DUSE_LIBCPP=ON"
fi

mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DBUILD_GENERIC=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX -DUSE_JEMALLOC=OFF $CMAKE_EXTRA_ARGS

make -j ${CPU_COUNT} ${VERBOSE_AT}
make install
