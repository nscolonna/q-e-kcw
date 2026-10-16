# Makefile for head (head.x)
# 
# Author G. Stenuit
#
include ../../make.inc

VPATH = @srcdir@

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
	  $(MOD_FLAG)../../PW/src \
	  $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../../dft-d3 \

HEADOBJS = \
bcast_gwl_input.o \
close_gwl.o \
openfilq_gwl.o \
gwl_readin.o \
lanczos_k.o \
solve_head.o 

QEMODS = $(BASEMODS)
LIBPWPH =  ../../PHonon/PH/libph.a ../../LR_Modules/liblrmod.a \
	   ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a


LIBMIN= ../minpack/minpacklib.a

# TLDEPS should contain the same list of dependencies as in the main Makefile
TLDEPS= phlibs
# Ensure that "make" run in this subdirectory builds the needed dependencies
include ../../install/tldeps

all : head.x  

head.x : head.o $(HEADOBJS) $(LIBPWPH) $(LIBMIN)
	$(LD) $(LDFLAGS) -o head.x head.o \
		 $(HEADOBJS) $(LIBPWPH)  $(QEMODS) $(LIBMIN) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../GWW/head/head.x . )

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L


ifneq ($(TOLERATE_MISSING_DEPEND),)
-include make.depend
else
ifeq ($(wildcard make.depend),)
$(error make.depend not found in $(CURDIR) -- run "make depend" first, or "make veryclean"/"make distclean" to clean without it)
endif
include make.depend
endif


