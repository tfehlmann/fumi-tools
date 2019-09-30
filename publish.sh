#!/usr/bin/env bash

set -e

SOURCE_DIR=$(realpath "$1")
VERSION="$2"

# build static binaries
mkdir /tmp/build && cd /tmp/build && cmake "${SOURCE_DIR}" -DCMAKE_INSTALL_PREFIX="${SOURCE_DIR}"/dist -DBUILD_SHARED_LIBS=OFF && make -j4 install
cd "${SOURCE_DIR}"/dist && tar cfv "fumi_tools_${VERSION}_static_linux_x86_64.tar.gz" bin

# build source tgz that includes version.hpp generated previously
mkdir /tmp/"${SOURCE_DIR}$$"
cp -r "${SOURCE_DIR}" /tmp/"${SOURCE_DIR}$$/fumi_tools-${VERSION}"
find /tmp/"${SOURCE_DIR}$$" -type d -name '.git' | xargs -n1 rm -rf
cd /tmp/"${SOURCE_DIR}$$"
tar cfv fumi_tools-${VERSION}.tar.gz fumi_tools-${VERSION}
mv fumi_tools-${VERSION}.tar.gz "${SOURCE_DIR}"/dist