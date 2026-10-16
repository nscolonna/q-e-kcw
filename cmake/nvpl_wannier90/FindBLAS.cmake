# Shim FindBLAS.cmake, active only while configuring the internal Wannier90
# submodule with -DBLA_VENDOR=NVPL (see external/wannier90.cmake, which
# prepends this directory to CMAKE_MODULE_PATH just for that add_subdirectory
# call). Companion to FindLAPACK.cmake in this same directory -- see there
# for why a shim is needed at all.
#
# Unlike LAPACK::LAPACK, QE's own top-level CMakeLists.txt never builds a
# standalone BLAS::BLAS target for the NVPL case (QE only ever links
# LAPACK::LAPACK, which pulls in nvpl_blas transitively). Wannier90 links
# BLAS::BLAS and LAPACK::LAPACK as two separate targets though, so this shim
# has to define BLAS::BLAS itself, wrapping the matching nvpl::blas_lp64_omp/
# seq target the same way the top-level CMakeLists.txt wraps LAPACK's.
find_package(nvpl REQUIRED COMPONENTS blas)
if(QE_ENABLE_OPENMP)
    set(_qe_nvpl_blas_target nvpl::blas_lp64_omp)
else()
    set(_qe_nvpl_blas_target nvpl::blas_lp64_seq)
endif()
if(NOT TARGET BLAS::BLAS)
    add_library(BLAS::BLAS INTERFACE IMPORTED GLOBAL)
    target_link_libraries(BLAS::BLAS INTERFACE ${_qe_nvpl_blas_target})
endif()
unset(_qe_nvpl_blas_target)
set(BLAS_FOUND TRUE)
