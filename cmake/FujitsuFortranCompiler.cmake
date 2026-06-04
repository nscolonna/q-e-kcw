############################################################
# This is based on experience from thea Odyssey cluster
# at the University of Tokyo.
############################################################

# Stop if OFFLOAD enabled, as this is not tested
if(QE_ENABLE_OFFLOAD)
	message(FATAL_ERROR "The current setup for the Fujitsu fortran compiler does NOT "
		            "support the use of offload devices.")
endif()

############################################################
# Set flag to allow allocating assignment,
# and -Kfast for code optiomizations.
# this is not default for the Fujitsu compiler.
############################################################

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Nalloc_assign -Kfast")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Kfast")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Kfast")


############################################################
# Setup  BLAS and LAPACK libraries
############################################################

if(QE_ENABLE_OPENMP AND NOT BLA_VENDOR)
   set(BLA_VENDOR Fujitsu_SSL2BLAMPSVE)
elseif(NOT BLA_VENDOR)
   set(BLA_VENDOR Fujitsu_SSL2SVE)
endif()


