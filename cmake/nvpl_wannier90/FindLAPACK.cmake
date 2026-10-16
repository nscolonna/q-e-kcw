# Shim FindLAPACK.cmake, active only while configuring the internal Wannier90
# submodule with -DBLA_VENDOR=NVPL (see external/wannier90.cmake, which
# prepends this directory to CMAKE_MODULE_PATH just for that add_subdirectory
# call).
#
# Wannier90 v4's own CMakeLists.txt calls find_package(LAPACK REQUIRED)
# itself. CMake's stock FindLAPACK module (as shipped with the CMake versions
# QE supports) has no notion of NVPL as a BLA_VENDOR, so letting that call
# run unshimmed would either fail outright or silently resolve a second,
# different LAPACK than the rest of QE just linked -- exactly the
# inconsistent-LAPACK problem QE_LAPACK_INTERNAL is refused for above.
#
# QE's own top-level CMakeLists.txt Lapack section already builds a
# LAPACK::LAPACK INTERFACE IMPORTED GLOBAL target wrapping the right
# nvpl::lapack_lp64_omp/seq target, and it runs before add_subdirectory(external)
# is reached, so that target is already visible here (GLOBAL imported targets
# are visible from any directory scope). Reuse it instead of re-deriving it.
if(NOT TARGET LAPACK::LAPACK)
    message(FATAL_ERROR
        "NVPL FindLAPACK shim: expected LAPACK::LAPACK to already exist "
        "(built by the top-level Lapack section for BLA_VENDOR=NVPL) before "
        "Wannier90's own find_package(LAPACK REQUIRED) runs.")
endif()
set(LAPACK_FOUND TRUE)
