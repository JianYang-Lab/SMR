#!/usr/bin/env bash

set -e

CWD=$(pwd)
APP_NAME=smr

# BUILD_TYPE=Release
BUILD_TYPE=RelWithDebInfo

fresh_build=""
install_app=0
cmake_gen=0

function usage {
    echo "usage: build.sh [-f] [-g] [-h]"
    echo "   ";
    echo "  -f | --fresh             : Fresh build";
    echo "  -g | --generate          : Build with cmake generation";
    echo "  -i | --install           : Install";
    echo "  -h | --help              : This message";
}

function parse_args {
    # positional args
    args=()

    # named args
    while [ "$1" != "" ]; do
        case "$1" in
            -f | --fresh )     fresh_build="--fresh";;
            -g | --generate )  cmake_gen=1;;
            -i | --install )   install_app=1;;
            -h | --help )      usage;                   exit;; # quit and show usage
            * )                args+=("$1")             # if no match, add it to the positional args
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

    echo "Current working directory: ${CWD}"

    module load gcc cmake intelmkl

    if [[ $cmake_gen == 1 ]]; then
        # Generate cmake
        cmake ${fresh_build} -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
            -DCMAKE_TOOLCHAIN_FILE=cmake/hpc-toolchain.cmake \
            -DCMAKE_INSTALL_PREFIX=${CWD}/build/${BUILD_TYPE}/installed/usr \
            -B build/${BUILD_TYPE} -S .
    fi

    cmake --build build/${BUILD_TYPE}

    if [[ $install_app == 1 ]]; then
        # Install
        cmake --install build/${BUILD_TYPE}

        installed_lib=build/${BUILD_TYPE}/installed/usr/lib
        if [[ ! -d ${installed_lib} ]]; then
            mkdir -p installed_lib
        fi
        cp /soft/compiler/intel/oneapi-2022.2/mkl/2022.1.0/lib/intel64/libmkl_avx512.so.2 ${installed_lib}/
        cp /soft/compiler/intel/oneapi-2022.2/mkl/2022.1.0/lib/intel64/libmkl_avx2.so.2 ${installed_lib}/
        cp /soft/compiler/intel/oneapi-2022.2/mkl/2022.1.0/lib/intel64/libmkl_def.so.2 ${installed_lib}/
        strip ${installed_lib}/libmkl_avx512.so.2
        strip ${installed_lib}/libmkl_avx2.so.2
        strip ${installed_lib}/libmkl_def.so.2

        # You should set the linuxdeploy and appimagetool path.
        echo "Path of linuxdeploy: ${LINUX_DEPLOY_BIN}"
        echo "Path of appimagetool: ${APP_IMAGE_TOOL_BIN}"
        echo "Path of appimagetool: ${APP_IMAGE_RUNTIME_FILE}"
        # Packaging
        ${LINUX_DEPLOY_BIN} --appdir build/${BUILD_TYPE}/installed \
            --executable build/${BUILD_TYPE}/installed/usr/bin/${APP_NAME} \
            --desktop-file build/${BUILD_TYPE}/installed/usr/share/applications/${APP_NAME}.desktop \
            --icon-file build/${BUILD_TYPE}/installed/usr/share/icons/hicolor/256x256/apps/${APP_NAME}.png

        ${APP_IMAGE_TOOL_BIN} build/${BUILD_TYPE}/installed --runtime-file ${APP_IMAGE_RUNTIME_FILE}
    fi
}

run "$@";
