# SMR CMakeLists.txt

cmake_minimum_required(VERSION 3.16...3.28)

# Version: src/config.h is the single source of truth
file(STRINGS "${CMAKE_SOURCE_DIR}/src/config.h" SMR_VERSION_LINE REGEX "#define SMR_VERSION")
string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" SMR_VERSION "${SMR_VERSION_LINE}")
project(smr VERSION ${SMR_VERSION} LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

option(BUILD_WITH_MIMALLOC "Use mimalloc for improved memory allocation performance" ON)

# ---- Dependencies ----
include(cmake/CPM.cmake)
set(FETCHCONTENT_QUIET ON)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)

# Silence deprecation warnings from third-party projects with old
# cmake_minimum_required (e.g. Eigen 3.4.0). Only effective on CMake >= 4.0.
set(CMAKE_POLICY_VERSION_MINIMUM 3.10)

include(cmake/SMRBlas.cmake)
smr_configure_blas(smr_blas_ok)

find_package(OpenMP REQUIRED)

CPMAddPackage(NAME eigen GITHUB_REPOSITORY eigen-mirror/eigen GIT_TAG 3.4.0
              EXCLUDE_FROM_ALL YES SYSTEM YES OPTIONS "BUILD_TESTING OFF" "EIGEN_BUILD_DOC OFF")
# Fetch zlib as a release tarball (not git): zlib's CMakeLists renames
# zconf.h inside its own source tree, which permanently dirties CPM's
# git cache and triggers a 'cache is dirty' warning on every configure.
CPMAddPackage(NAME zlib
              URL https://github.com/madler/zlib/archive/refs/tags/v1.3.1.tar.gz
              URL_HASH SHA256=17e88863f3600672ab49182f217281b6fc4d3c762bde361935e436a95214d05c
              EXCLUDE_FROM_ALL YES SYSTEM YES OPTIONS "ZLIB_BUILD_EXAMPLES OFF")
CPMAddPackage(NAME spdlog GITHUB_REPOSITORY gabime/spdlog VERSION 1.15.3
              EXCLUDE_FROM_ALL YES SYSTEM YES OPTIONS "SPDLOG_BUILD_EXAMPLE OFF" "SPDLOG_FMT_EXTERNAL ON")
CPMAddPackage(NAME fmt GITHUB_REPOSITORY fmtlib/fmt GIT_TAG 12.1.0
              EXCLUDE_FROM_ALL YES SYSTEM YES OPTIONS "FMT_DOC OFF" "FMT_TEST OFF")
if(BUILD_WITH_MIMALLOC)
    CPMAddPackage(NAME mimalloc GITHUB_REPOSITORY microsoft/mimalloc VERSION 2.2.7
                  EXCLUDE_FROM_ALL YES SYSTEM YES
                  OPTIONS "MI_OVERRIDE ON" "MI_BUILD_SHARED OFF" "MI_BUILD_STATIC ON" "MI_BUILD_TESTS OFF" "MI_BUILD_OBJECT OFF")
endif()

# ---- Target ----
set(SRCS src/bfile.cpp src/CommFunc.cpp src/mmap_read.cpp src/perf_timer.cpp src/SMR.cpp
         src/SMR_data.cpp src/SMR_data_p1.cpp src/SMR_data_p2.cpp src/SMR_data_p3.cpp
         src/SMR_plot.cpp src/StatFunc.cpp src/StrFunc.cpp)
set(DCDFLIB_SRCS third_party/dcdflib/dcdflib.c third_party/dcdflib/ipmpar.c)

add_executable(${PROJECT_NAME} ${SRCS} ${DCDFLIB_SRCS})
target_include_directories(${PROJECT_NAME} PRIVATE third_party/dcdflib)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX Eigen3::Eigen zlib spdlog::spdlog fmt::fmt)

include(CheckCXXCompilerFlag)
set(smr_arch_flag "none")
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    check_cxx_compiler_flag("-march=x86-64-v3" SMR_HAS_MARCH_X86_64_V3)
    if(SMR_HAS_MARCH_X86_64_V3)
        target_compile_options(${PROJECT_NAME} PRIVATE -march=x86-64-v3)
        set(smr_arch_flag "-march=x86-64-v3")
    endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64|ARM64")
    check_cxx_compiler_flag("-march=armv8.5-a" SMR_HAS_MARCH_ARMV85A)
    if(SMR_HAS_MARCH_ARMV85A)
        target_compile_options(${PROJECT_NAME} PRIVATE -march=armv8.5-a)
        set(smr_arch_flag "-march=armv8.5-a")
    endif()
endif()

if(smr_blas_ok)
    smr_target_enable_blas(${PROJECT_NAME})
endif()
if(BUILD_WITH_MIMALLOC)
    target_link_libraries(${PROJECT_NAME} mimalloc-static)  # global malloc override at link time
endif()

# Keep compile_commands.json in the source root for clangd
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" "${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json")
endif()

# ---- Install ----
include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Required by AppImage packaging (scripts/hpc/build.sh, linuxdeploy)
install(FILES resources/smr.desktop DESTINATION share/applications)
install(FILES resources/smr.png DESTINATION share/icons/hicolor/256x256/apps)
install(FILES resources/smr.svg DESTINATION share/icons/hicolor/scalable/apps)

# ---- Summary ----
message(STATUS "smr ${PROJECT_VERSION} | ${CMAKE_BUILD_TYPE} | BLAS: ${SMR_BLAS_BACKEND} | mimalloc: ${BUILD_WITH_MIMALLOC} | arch: ${smr_arch_flag}")
