#SMR CMakeLists.txt

cmake_minimum_required(VERSION 3.16)
project(smr C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(BUILD_WITH_MKL ON)
option(SMR_VERBOSE_CONFIGURE "Print extended CMake configuration diagnostics" OFF)

if(UNIX AND NOT APPLE)
    set(LINUX TRUE)
endif()

#information of CPU architecture
message(STATUS "CPU_ARCH: " ${CMAKE_HOST_SYSTEM_PROCESSOR})

#information of operation system
message(STATUS ${CMAKE_HOST_SYSTEM})
if (UNIX)
    set(system "UNIX-like")
elseif(APPLE)
    set(system "macOS")
elseif(WIN32)
    set(system "windows")
else()
    message(FATAL_ERROR "Operation system not recognized, only support unix apple and windows.")
endif()
message(STATUS "OS: " ${system})

#information of compilation toolchain
message(STATUS "CMAKE_C_COMPILER:  " ${CMAKE_C_COMPILER} " " ${CMAKE_C_COMPILER_ID} ":" ${CMAKE_C_COMPILER_VERSION} " " ${CMAKE_C_STANDARD})
message(STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} " " ${CMAKE_CXX_COMPILER_ID} ":" ${CMAKE_CXX_COMPILER_VERSION} " " ${CMAKE_CXX_STANDARD})
message(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

if(SMR_VERBOSE_CONFIGURE)
    message(STATUS "CUSTORM_INCLUDE_SERACH_PATH: ${CUSTORM_INCLUDE_SERACH_PATH}")
    message(STATUS "CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES: ${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
    message(STATUS "CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES: ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES}")
    message(STATUS "CMAKE_C_IMPLICIT_LINK_DIRECTORIES: ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
    message(STATUS "CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES: ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES}")
    message(STATUS "CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
    message(STATUS "CMAKE_AR: ${CMAKE_AR}")
    message(STATUS "CMAKE_RANLIB: ${CMAKE_RANLIB}")
    message(STATUS "CMAKE_COMMAND: ${CMAKE_COMMAND}")
    message(STATUS "CMAKE_VERSION: ${CMAKE_VERSION}")
    message(STATUS "CMAKE_GENERATOR: ${CMAKE_GENERATOR}")
    message(STATUS "CMAKE_BUILD_TOOL: ${CMAKE_BUILD_TOOL}")
    message(STATUS "CMAKE_SHARED_LIBRARY_PREFIX: ${CMAKE_SHARED_LIBRARY_PREFIX}")
    message(STATUS "CMAKE_SHARED_LIBRARY_SUFFIX: ${CMAKE_SHARED_LIBRARY_SUFFIX}")
    message(STATUS "CMAKE_STATIC_LIBRARY_PREFIX: ${CMAKE_STATIC_LIBRARY_PREFIX}")
    message(STATUS "CMAKE_STATIC_LIBRARY_SUFFIX: ${CMAKE_STATIC_LIBRARY_SUFFIX}")
    message(STATUS "CMAKE_FIND_LIBRARY_PREFIXES: ${CMAKE_FIND_LIBRARY_PREFIXES}")
    message(STATUS "CMAKE_FIND_LIBRARY_SUFFIXES: ${CMAKE_FIND_LIBRARY_SUFFIXES}")
    message(STATUS "CMAKE_PROJECT_NAME: ${CMAKE_PROJECT_NAME}")
    message(STATUS "CMAKE_PROJECT_VERSION: ${CMAKE_PROJECT_VERSION}")
    message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
    message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
endif()

if(LINUX AND NOT OMP_LIBRARY)
    set(OMP_LIBRARY "/usr/lib/libgomp.so")
endif()

if(LINUX AND BUILD_WITH_MKL)
    if(NOT MKL_DIR)
        set(MKL_DIR "/opt/intel/oneapi/mkl/latest/lib/cmake/mkl")
    endif()
    set(MKL_THREADING "gnu_thread")
    set(MKL_INTERFACE "lp64")
endif()


# APPLE use Apple Accelerate
if(NOT APPLE AND BUILD_WITH_MKL)
    find_package(MKL REQUIRED)
endif()

find_package(OpenMP REQUIRED)

# Fetch 3rd parties
include(cmake/CPM.cmake)
set(FETCHCONTENT_QUIET ON)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)

CPMAddPackage(
    NAME eigen
    GITHUB_REPOSITORY eigen-mirror/eigen
    GIT_TAG 3.4.0
    EXCLUDE_FROM_ALL YES
    SYSTEM YES
    OPTIONS
        "BUILD_TESTING OFF"
        "EIGEN_BUILD_DOC OFF"
)
CPMAddPackage(
    NAME zlib
    GITHUB_REPOSITORY madler/zlib
    VERSION 1.3.1
    EXCLUDE_FROM_ALL YES
    SYSTEM YES
    OPTIONS
        "ZLIB_BUILD_EXAMPLES OFF"
)
CPMAddPackage(
    NAME spdlog
    GITHUB_REPOSITORY gabime/spdlog
    VERSION 1.15.3
    EXCLUDE_FROM_ALL YES
    SYSTEM YES
    OPTIONS
        "SPDLOG_BUILD_EXAMPLE OFF"
        "SPDLOG_BUILD_EXAMPLE_HO OFF"
        "SPDLOG_BUILD_TESTS OFF"
        "SPDLOG_BUILD_TESTS_HO OFF"
        "SPDLOG_BUILD_BENCH OFF"
)
CPMAddPackage(
    NAME fmt
    GITHUB_REPOSITORY fmtlib/fmt
    GIT_TAG 12.1.0
    EXCLUDE_FROM_ALL YES
    SYSTEM YES
    OPTIONS
        "FMT_DOC OFF"
        "FMT_TEST OFF"
        "FMT_FUZZ OFF"
        "FMT_CUDA_TEST OFF"
)

file(GLOB SRCS "${CMAKE_SOURCE_DIR}/src/*.cpp")
file(GLOB DCDFLIB_SRCS "${CMAKE_SOURCE_DIR}/third_party/dcdflib/*.c")
add_executable(${PROJECT_NAME} ${SRCS} ${DCDFLIB_SRCS})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/third_party/dcdflib)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX Eigen3::Eigen zlib spdlog::spdlog fmt::fmt)

# -march=native require rebuild it on AMD/Intel
target_compile_options(${PROJECT_NAME} PRIVATE -march=x86-64)

if(NOT APPLE AND BUILD_WITH_MKL)
    target_link_libraries(${PROJECT_NAME} MKL::MKL)
    target_compile_definitions(${PROJECT_NAME} PRIVATE EIGEN_USE_MKL_ALL)
    target_compile_options(${PROJECT_NAME} PRIVATE -fpermissive)
endif()


IF(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
        ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
    )
ENDIF()

if(SMR_VERBOSE_CONFIGURE)
    get_property(imported_targets DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" PROPERTY IMPORTED_TARGETS)
    foreach(target IN LISTS imported_targets)
        message(STATUS "imported target: ${target}")
    endforeach()

    get_property(local_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
    foreach(target IN LISTS local_targets)
        message(STATUS "built target: ${target}")
    endforeach()
endif()

# Installation
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)

# Targets
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)

# Resource
install(FILES ${CMAKE_SOURCE_DIR}/resources/smr.desktop
    DESTINATION share/applications)
install(FILES ${CMAKE_SOURCE_DIR}/resources/smr.png
    DESTINATION share/icons/hicolor/256x256/apps)
install(FILES ${CMAKE_SOURCE_DIR}/resources/smr.svg
    DESTINATION share/icons/hicolor/scalable/apps)
