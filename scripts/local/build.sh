#!/usr/bin/env bash

set -e

CWD=$(pwd)
# BUILD_TYPE=Release
BUILD_TYPE=RelWithDebInfo
fresh_build=""
cmake_gen=0

function usage {
    echo "usage: build.sh [-f] [-g] [-h]"
    echo "   ";
    echo "  -f | --fresh             : Fresh build";
    echo "  -g | --generate          : Build with cmake generation";
    echo "  -h | --help              : This message";
}

function parse_args {
    # positional args
    args=()

    # named args
    while [ "$1" != "" ]; do
        case "$1" in
            -f | --an_arg )               fresh_build="--fresh";;
            -g | --some_more_args )       cmake_gen=1;;
            -h | --help )                 usage;                   exit;; # quit and show usage
            * )                           args+=("$1")             # if no match, add it to the positional args
        esac
        shift # move to next kv pair
    done

    set -- "${args[@]}"

    if [ ${#args[@]} -ne 0 ]; then
        echo "not support args: ${args[@]}"
        exit;
    fi
}


function run {
    parse_args "$@"

    # Set OpenMP root for macOS (Homebrew libomp)
    if [[ "$(uname)" == "Darwin" ]]; then
        export OpenMP_ROOT=/opt/homebrew/opt/libomp
    fi

    if [[ $cmake_gen == 1 ]]; then
        cmake ${fresh_build} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
            -DCMAKE_INSTALL_PREFIX=${CWD}/build/${BUILD_TYPE}/installed/usr \
            -DBUILD_WITH_MKL=ON \
            -DMKL_DIR="/opt/intel/oneapi/mkl/latest/lib/cmake/mkl" \
            -DSMR_VERBOSE_CONFIGURE=ON \
            -B build/${BUILD_TYPE} -S . -G Ninja
    fi

    cmake --build build/${BUILD_TYPE}
}


run "$@";
