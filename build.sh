#!/usr/bin/env bash

set -e

SOURCE_DIR=$(realpath "$1")
VERSION="$2"

# build static binaries
mkdir /tmp/build && cd /tmp/build
${CMAKE} -S "${SOURCE_DIR}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="${SOURCE_DIR}/dist/${PLATFORM}" -DBUILD_SHARED_LIBS=OFF -DBUILD_GENERIC=ON "${CMAKE_ARGS}"
make -j4 install

cd "${SOURCE_DIR}/dist/${PLATFORM}" && tar czfv "../fumi_tools_${VERSION}_static_${PLATFORM}.tar.gz" bin && cd .. && rm -rf "${PLATFORM}"

# remove the tag here so that the tag will only be created if all builds pass
curl --silent --request DELETE --header "PRIVATE-TOKEN: ${GITLAB_TOKEN}" "${CI_API_V4_URL}/projects/${CI_PROJECT_ID}/repository/tags/${VERSION}"

git tag -d ${VERSION} || true
