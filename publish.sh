#!/usr/bin/env bash

set -e

SOURCE_DIR=$(realpath "$1")
VERSION="$2"

# generate version.hpp
mkdir /tmp/build && cd /tmp/build && cmake "${SOURCE_DIR}" -DBUILD_SHARED_LIBS=OFF -DBUILD_GENERIC=ON

# build source tgz that includes version.hpp generated previously
mkdir /tmp/$$
cp -r "${SOURCE_DIR}" /tmp/"$$/fumi_tools-${VERSION}"
find /tmp/$$ -type d -name '.git' | xargs -n1 rm -rf
cd /tmp/$$
tar czfv fumi_tools-${VERSION}.tar.gz fumi_tools-${VERSION}
mv fumi_tools-${VERSION}.tar.gz "${SOURCE_DIR}"/dist
