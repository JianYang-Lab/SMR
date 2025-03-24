# smr CMakeLists.txt
cmake_minimum_required(VERSION 3.13)
project(smr C CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")


#set project source directory
set(PRJ_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(CUSTORM_INCLUDE_SERACH_PATH "" CACHE PATH "")
set(CUSTORM_LIB_SEARCH_PATH "" CACHE PATH "")
set(BUILD_STATIC FALSE CACHE BOOL "")


#information of CPU architecture
message(STATUS "CPU_ARCH: " ${CMAKE_HOST_SYSTEM_PROCESSOR})
message("")

#information of operation system
message(STATUS ${CMAKE_HOST_SYSTEM})
if (UNIX)
    set(system "UNIX-like")
elseif(APPLE)
    set(system "MacOS")
elseif(WIN32)
    set(system "Windows")
else()
    message(FATAL_ERROR "Operation system not recognized, only support unix apple and windows.")
endif()
message(STATUS "OS: " ${system})
message("")

#information of compilation toolchain
message(STATUS "CMAKE_C_COMPILER:  " ${CMAKE_C_COMPILER} " " ${CMAKE_C_COMPILER_ID} ":" ${CMAKE_C_COMPILER_VERSION} " " ${CMAKE_C_STANDARD})
message(STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} " " ${CMAKE_CXX_COMPILER_ID} ":" ${CMAKE_CXX_COMPILER_VERSION} " " ${CMAKE_CXX_STANDARD})
message(STATUS "CMAKE_AR: " ${CMAKE_AR})
message(STATUS "CMAKE_RANLIB: " ${CMAKE_RANLIB})
message(STATUS "CMAKE_COMMAND: " ${CMAKE_COMMAND})
message(STATUS "CMAKE_VERSION: " ${CMAKE_VERSION})
message(STATUS "CMAKE_GENERATOR: " ${CMAKE_GENERATOR})
message(STATUS "CMAKE_BUILD_TOOL: " ${CMAKE_BUILD_TOOL})

message(STATUS "CMAKE_SHARED_LIBRARY_PREFIX: " ${CMAKE_SHARED_LIBRARY_PREFIX})
message(STATUS "CMAKE_SHARED_LIBRARY_SUFFIX: " ${CMAKE_SHARED_LIBRARY_SUFFIX})
message(STATUS "CMAKE_STATIC_LIBRARY_PREFIX: " ${CMAKE_STATIC_LIBRARY_PREFIX})
message(STATUS "CMAKE_STATIC_LIBRARY_SUFFIX: " ${CMAKE_STATIC_LIBRARY_SUFFIX})

message(STATUS "CMAKE_FIND_LIBRARY_PREFIXES: " ${CMAKE_FIND_LIBRARY_PREFIXES})
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
message(STATUS "CMAKE_FIND_LIBRARY_SUFFIXES: " ${CMAKE_FIND_LIBRARY_SUFFIXES})

message(STATUS "CMAKE_PROJECT_NAME: " ${CMAKE_PROJECT_NAME})
message(STATUS "CMAKE_PROJECT_VERSION: " ${CMAKE_PROJECT_VERSION})
message(STATUS "CMAKE_BINARY_DIR: " ${CMAKE_BINARY_DIR})
message(STATUS "CMAKE_SOURCE_DIR: " ${CMAKE_SOURCE_DIR})
message("")


# detect zlib
set(zlib_path FALSE CACHE PATH "zlib package location")
set(zlib_include_path FALSE CACHE PATH "zlib head file location")
set(zlib_lib_path FALSE CACHE PATH "zlib library location")
if (zlib_path)
    set(zlib_include_dir ${zlib_path}/include)
    set(zlib_lib_dir ${zlib_path}/lib)
elseif (zlib_include_path AND zlib_lib_path)
    set(zlib_include_dir ${zlib_include_path})
    set(zlib_lib_dir ${zlib_lib_path})
else()
    find_path(zlib_include_dir zlib.h ${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES} ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
    find_library(zlib_lib_dir z ${CMAKE_C_IMPLICIT_LINK_DIRECTORIES} ${CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES})
    if (zlib_include_dir AND zlib_lib_dir)
        get_filename_component(zlib_lib_dir ${zlib_lib_dir} DIRECTORY)
        message(STATUS "zlib founded.")
    else()
        message(FATAL_ERROR "Can not find zlib lib, please specific them by -Dzlib_path or -Dzlib_include_path and -Dzlib_lib_path")
    endif()
endif()


# detect Eigen
set(eigen_path FALSE CACHE PATH "path where Eigen located.")
# set(eigen_path FALSE)
if (eigen_path)
    message(STATUS "You set eigen path to: " ${eigen_path})
    set(eigen_dir ${eigen_path}/eigen3)
else()
    find_path(eigen_dir "eigen3/Eigen" ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES} ${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES})
    if (eigen_dir)
        set(eigen_dir ${eigen_dir}/eigen3)
        message(STATUS "Find Eigen: " ${eigen_dir})
    else()
        message(FATAL_ERROR "can not find eigen, please using -Deigen_path=/where/your/eigen/located to set eigen path")
    endif()
endif()

# build
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
add_executable(smr ${PRJ_SOURCE_DIR}/bfile.cpp ${PRJ_SOURCE_DIR}/CommFunc.cpp
                   ${PRJ_SOURCE_DIR}/dcdflib.cpp ${PRJ_SOURCE_DIR}/SMR_data_p1.cpp
                   ${PRJ_SOURCE_DIR}/SMR_data_p2.cpp ${PRJ_SOURCE_DIR}/SMR_data_p3.cpp
                   ${PRJ_SOURCE_DIR}/SMR_data.cpp ${PRJ_SOURCE_DIR}/SMR_plot.cpp
                   ${PRJ_SOURCE_DIR}/SMR.cpp ${PRJ_SOURCE_DIR}/StatFunc.cpp
                   ${PRJ_SOURCE_DIR}/StrFunc.cpp)
    if (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64")
    if (system STREQUAL "UNIX-like")
        if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
            target_compile_options(smr PRIVATE -fopenmp)
            target_include_directories(smr PRIVATE ${zlib_include_dir} PRIVATE ${eigen_dir} PRIVATE ${CMAKE_SOURCE_DIR}/include)
            target_link_directories(smr PRIVATE ${zlib_lib_dir})
            if (NOT BUILD_STATIC)
                target_link_libraries(smr z gomp m pthread)
            else()
                message(WARNING "You now try to link libs staticlly, You need make sure all static lib was installed.")
                target_link_libraries(smr -static z gomp m dl pthread)
            endif()
        elseif (CMAKE_C_COMPILER_ID STREQUAL "Clang")
            target_include_directories(smr PRIVATE ${zlib_include_dir} PRIVATE ${eigen_dir} PRIVATE ${CMAKE_SOURCE_DIR}/include)
            target_link_directories(smr PRIVATE ${zlib_lib_dir})
            if (NOT BUILD_STATIC)
                target_link_libraries(smr z gomp m pthread)
            else()
                message(WARNING "You now try to link libs staticlly, You need make sure all static lib was installed.")
                target_link_libraries(smr -static z gomp m dl pthread)
            endif()
        endif()
    elseif (system STREQUAL "MacOS")
    elseif (system STREQUAL "Windows")
    else()
        message(FATAL_ERROR "OS not recgnized.")
    endif()
elseif (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64")
    if (system STREQUAL "UNIX-like")
        target_include_directories(smr PRIVATE ${zlib_include_dir} PRIVATE ${eigen_dir} PRIVATE ${CMAKE_SOURCE_DIR}/include)
        target_link_directories(smr PRIVATE ${zlib_lib_dir})
        target_link_libraries(smr m
            -force_load /opt/homebrew/Cellar/zlib/1.3/lib/libz.a
            -force_load /opt/homebrew/Cellar/libomp/17.0.6/lib/libomp.a
        )
    elseif (system STREQUAL "Windows")
    else()
        message(FATAL_ERROR "OS not recgnized.")
    endif()
else()
    message(STATUS "Until now, only x86_64 and arm64 architecture was tested.")
endif()
