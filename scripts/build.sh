#!/usr/bin/env bash
#
# Build static fumi_tools binaries for a given platform.
#
# Usage: scripts/build.sh <platform> [extra cmake args...]
#
# Example:
#   scripts/build.sh linux_x86_64 -DUSE_LIBCPP=OFF
#   scripts/build.sh apple_darwin_arm64 -DUSE_LIBCPP=ON -DBUILD_GENERIC=OFF

set -euo pipefail

PLATFORM="${1:?Usage: build.sh <platform> [cmake args...]}"
shift

SOURCE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="${SOURCE_DIR}/build"
DIST_DIR="${SOURCE_DIR}/dist"

# Use Ninja on macOS to avoid Make "multiple target patterns" bug with
# CMake-generated ExternalProject rules.
GENERATOR_ARGS=()
if [[ "$(uname)" == "Darwin" ]] && command -v ninja &>/dev/null; then
    GENERATOR_ARGS=(-G Ninja)
fi

mkdir -p "$BUILD_DIR" && cd "$BUILD_DIR"

cmake "$SOURCE_DIR" \
  "${GENERATOR_ARGS[@]}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${DIST_DIR}/${PLATFORM}" \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_GENERIC=ON \
  "$@"

cmake --build . --target install -j "$(getconf _NPROCESSORS_ONLN)"

mkdir -p "$DIST_DIR"
tar czf "${DIST_DIR}/fumi_tools_static_${PLATFORM}.tar.gz" -C "${DIST_DIR}/${PLATFORM}" bin
