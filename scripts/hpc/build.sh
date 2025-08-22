#!/usr/bin/env sh

set -e

# You should cd to the project root
CWD=$(pwd)

# BUILD_TYPE=Release
# BUILD_TYPE=Debug
BUILD_TYPE=RelWithDebInfo
APP_NAME=smr

# You should set the linuxdeploy and appimagetool path.
echo "Current working directory: ${CWD}"
echo "Path of linuxdeploy ${LINUX_DEPLOY_BIN}"
echo "Path of appimagetool ${APP_IMAGE_TOOL_BIN}"
echo "Path of appimagetool ${APP_IMAGE_RUNTIME_FILE}"

module load gcc cmake intelmkl

# Generate cmake
cmake --fresh -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
    -DCMAKE_TOOLCHAIN_FILE=cmake/hpc-toolchain.cmake \
    -DCMAKE_INSTALL_PREFIX=${CWD}/build/${BUILD_TYPE}/installed/usr \
    -B build/${BUILD_TYPE} -S .

# Build
cmake --build build/${BUILD_TYPE}

# Install
# cmake --install build/${BUILD_TYPE}

# cp -v /soft/compiler/intel/oneapi-2022.2/mkl/2022.1.0/lib/intel64/libmkl_*.so.* build/${BUILD_TYPE}/installed/usr/lib/

# # Packaging
# ${LINUX_DEPLOY_BIN} --appdir build/${BUILD_TYPE}/installed \
#     --executable build/${BUILD_TYPE}/installed/usr/bin/${APP_NAME} \
#     --desktop-file build/${BUILD_TYPE}/installed/usr/share/applications/${APP_NAME}.desktop \
#     --icon-file build/${BUILD_TYPE}/installed/usr/share/icons/hicolor/256x256/apps/${APP_NAME}.png

# ${APP_IMAGE_TOOL_BIN} build/${BUILD_TYPE}/installed --runtime-file ${APP_IMAGE_RUNTIME_FILE}
