############################################################
# This is based on experience from thea Odyssey cluster
# at the University of Tokyo.
#
# Note that this compiler option will overwrite some of 
# the choices from CMakeLists.txt.
#
# It seems that FFTW is detected automatically so this is
# handled by the functionality in CMakeLists.txt
#
# But the BLAS, LAPACK and SCALAPACK options should be set
# using the convinient compiler flags provided by the Fujitsu
# compiler.
#
# Given that this compiler is mainly used in HPC environments
# we will here enforce SCALAPACK even without setting the
# appropriate cmake flags.
#
############################################################

# Stop if OFFLOAD enabled, as this is not tested
if(QE_ENABLE_OFFLOAD)
	message(FATAL_ERROR "The current setup for the Fujitsu fortran compiler does NOT "
		            "support the use of offload devices.")
endif()

############################################################
# Set flag to allow allocating assignment, 
# this is not default for the Fujitsu compiler.
############################################################

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Nalloc_assign -Kfast")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Kfast")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Kfast")


############################################################
# Setup OpenMP, BLAS, LAPACK and SCALAPACK libraries
############################################################

# Always compile with SCALAPACK ENABLED
qe_add_global_compile_definitions(__SCALAPACK)

if(QE_ENABLE_OPENMP)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Kopenmp -SCALAPACK -SSL2BLAMP")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Kopenmp -SCALAPACK -SSL2BLAMP")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Kopenmp -SCALAPACK -SSL2BLAMP")
else()
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -SCALAPACK -SSL2")
	set(CMAKE_C_FLAGS "${CMAKE_Fortran_FLAGS} -SCALAPACK -SSL2")
	set(CMAKE_CXX_FLAGS "${CMAKE_Fortran_FLAGS} -SCALAPACK -SSL2")
endif()


