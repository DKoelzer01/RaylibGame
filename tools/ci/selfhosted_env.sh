#!/bin/bash
# The self-hosted jobs' shared paths (#1016, #1028): appends OUT, RUN, ccache and the version to $GITHUB_ENV.
# Outside the checkout (actions/checkout cleans it): GNS, FetchContent, ccache and the run folder (the erosion bake
# cache in saves/) stay warm between runs. Needs REF (the branch) in the environment.
set -eu
root="$RUNNER_WORKSPACE/ci-persist"
mkdir -p "$root/out" "$root/run" "$root/ccache"
branch="${REF#refs/heads/}"
{
    echo "OUT=$root/out"
    echo "RUN=$root/run"
    echo "CCACHE_DIR=$root/ccache"
    echo "CCACHE_MAXSIZE=2G"
    echo "RAYLIBGAME_VERSION=$branch $(git rev-parse --short HEAD)"
} >> "$GITHUB_ENV"
if command -v ccache > /dev/null; then
    # CMake (3.17+) takes the launcher from the environment: tools/build_linux.sh needs no change.
    { echo "CMAKE_C_COMPILER_LAUNCHER=ccache"; echo "CMAKE_CXX_COMPILER_LAUNCHER=ccache"; } >> "$GITHUB_ENV"
fi
