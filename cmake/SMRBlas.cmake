# BLAS backend dispatcher for SMR.
#
# SMR_BLAS_BACKEND selects the BLAS backend Eigen uses:
#   mkl      - Intel oneMKL (default; uses cmake/SMRMKL.cmake)
#   openblas - OpenBLAS (system install preferred, CPM fetch as fallback)
#   none     - Eigen's native backend (no external BLAS)
#
# Backward compatibility: BUILD_WITH_MKL=OFF maps to backend "none".
# OPENBLAS_ROOT (cache or environment) points at a system OpenBLAS install.

include_guard(GLOBAL)

set(SMR_BLAS_BACKEND "mkl" CACHE STRING "BLAS backend for Eigen: mkl, openblas, none")
set_property(CACHE SMR_BLAS_BACKEND PROPERTY STRINGS mkl openblas none)
set(OPENBLAS_ROOT "" CACHE PATH "OpenBLAS installation root (optional)")

# Legacy flag: BUILD_WITH_MKL=OFF maps to backend "none", but only when the
# backend was not explicitly chosen (an explicit SMR_BLAS_BACKEND wins).
if(DEFINED BUILD_WITH_MKL AND NOT BUILD_WITH_MKL AND SMR_BLAS_BACKEND STREQUAL "mkl")
    set(SMR_BLAS_BACKEND "none" CACHE STRING "BLAS backend for Eigen: mkl, openblas, none" FORCE)
endif()

function(smr_configure_openblas out_var)
    set(smr_use_openblas FALSE)

    # 1. Prefer a system OpenBLAS (e.g. distro package or HPC module)
    find_library(OPENBLAS_LIBRARY NAMES openblas
                 HINTS ${OPENBLAS_ROOT} ${OPENBLAS_ROOT}/lib $ENV{OPENBLAS_ROOT} $ENV{OPENBLAS_ROOT}/lib)
    if(OPENBLAS_LIBRARY)
        set(OPENBLAS_TARGET "${OPENBLAS_LIBRARY}")
        set(smr_use_openblas TRUE)
        message(STATUS "OpenBLAS: using system library ${OPENBLAS_LIBRARY}")
    else()
        # 2. Fallback: fetch and build OpenBLAS (BLAS only via NO_LAPACK; no Fortran needed,
        #    and avoids a 'lapack' target name collision with Eigen's lapack subdir)
        message(STATUS "OpenBLAS: no system library found; fetching and building from source (one-time)")
        CPMAddPackage(
            NAME openblas
            GITHUB_REPOSITORY OpenMathLib/OpenBLAS
            GIT_TAG v0.3.30
            EXCLUDE_FROM_ALL YES
            OPTIONS
                "BUILD_SHARED_LIBS OFF"
                "NO_LAPACK ON"
                "NO_LAPACKE ON"
                "BUILD_TESTING OFF"
                "DYNAMIC_ARCH ON"
                "USE_THREAD ON"
                "USE_OPENMP OFF"
                "CMAKE_POSITION_INDEPENDENT_CODE ON"
        )
        if(TARGET openblas_static)
            set(OPENBLAS_TARGET openblas_static)
        elseif(TARGET OpenBLAS::OpenBLAS)
            set(OPENBLAS_TARGET OpenBLAS::OpenBLAS)
        else()
            set(OPENBLAS_TARGET openblas)
        endif()
        set(smr_use_openblas TRUE)
        message(STATUS "OpenBLAS: building target ${OPENBLAS_TARGET} (DYNAMIC_ARCH=ON)")
    endif()

    # CACHE INTERNAL so the value is visible at the top-level scope and in
    # smr_target_enable_openblas (functions' local scopes are not enough here)
    set(OPENBLAS_TARGET "${OPENBLAS_TARGET}" CACHE INTERNAL "")
    set(${out_var} "${smr_use_openblas}" CACHE INTERNAL "")
endfunction()

function(smr_target_enable_openblas target_name)
    if(NOT TARGET "${target_name}")
        message(FATAL_ERROR "smr_target_enable_openblas: target '${target_name}' does not exist.")
    endif()
    target_compile_definitions(${target_name} PRIVATE EIGEN_USE_BLAS)
    target_link_libraries(${target_name} ${OPENBLAS_TARGET})
endfunction()

function(smr_configure_blas out_var)
    if(SMR_BLAS_BACKEND STREQUAL "mkl")
        include(cmake/SMRMKL.cmake)
        smr_configure_mkl(${out_var})
        # sub-configure set ${out_var} in our scope; re-export it to the top level
        set(${out_var} "${${out_var}}" PARENT_SCOPE)
    elseif(SMR_BLAS_BACKEND STREQUAL "openblas")
        smr_configure_openblas(${out_var})
        set(${out_var} "${${out_var}}" PARENT_SCOPE)
    else()
        message(STATUS "BLAS: using Eigen's native backend (SMR_BLAS_BACKEND=none)")
        set(${out_var} FALSE PARENT_SCOPE)
    endif()
endfunction()

function(smr_target_enable_blas target_name)
    if(SMR_BLAS_BACKEND STREQUAL "mkl")
        include(cmake/SMRMKL.cmake)
        smr_target_enable_mkl(${target_name})
    elseif(SMR_BLAS_BACKEND STREQUAL "openblas")
        smr_target_enable_openblas(${target_name})
    endif()
endfunction()
