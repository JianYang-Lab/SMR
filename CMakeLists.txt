# smr CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(smr C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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

find_package(OpenMP REQUIRED)

# Fetch 3rd parties
include(FetchContent)

set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
set(FETCHCONTENT_QUIET OFF)

FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE)

set(ZLIB_USE_STATIC_LIBS OFF)
set(ZLIB_ROOT "")
FetchContent_Declare(
    ZLIB
    GIT_REPOSITORY https://github.com/madler/zlib.git
    GIT_TAG v1.3.1
    GIT_SHALLOW TRUE)

FetchContent_MakeAvailable(Eigen ZLIB)

file(GLOB SRCS "${CMAKE_SOURCE_DIR}/src/*.cpp")
add_executable(${PROJECT_NAME} ${SRCS})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX Eigen3::Eigen zlib)

# List targets
get_property(imported_targets DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" PROPERTY IMPORTED_TARGETS)
foreach(target IN LISTS imported_targets)
    message(STATUS "imported target: ${target}")
endforeach()

get_property(local_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
foreach(target IN LISTS local_targets)
    message(STATUS "built target: ${target}")
endforeach()

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
