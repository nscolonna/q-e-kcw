# Makefile for school

include ../../make.inc

VPATH = @srcdir@

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3 

FDOBJS = \
        stop_pp.o 


PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a 
QEMODS = $(BASEMODS)

# see comments in file tldeps for explanations
TLDEPS=pwlibs
include ../../install/tldeps

all : fd.x fd_ifc.x fd_ef.x

fd_ifc.x : fd_ifc.o libfd.a  $(FDOBJS) $(PWOBJS) $(QEMODS)  
	$(LD) $(LDFLAGS) -o $@ \
	fd_ifc.o libfd.a $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PHonon/FD/$@ . )

fd.x : fd.o libfd.a  $(FDOBJS) $(PWOBJS) $(QEMODS)  
	$(LD) $(LDFLAGS) -o $@ \
	fd.o libfd.a $(PWOBJS) $(QEMODS) $(QELIBS) 
	- ( cd ../../bin ; ln -fs ../PHonon/FD/$@ . )

fd_ef.x : fd_ef.o libfd.a  $(FDOBJS) $(PWOBJS) $(QEMODS)  
	$(LD) $(LDFLAGS) -o $@ \
	fd_ef.o libfd.a $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PHonon/FD/$@ . )

libfd.a : $(FDOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *~ *_tmp.f90 *.d *.mod *.i *.L libfd.a

ifneq ($(TOLERATE_MISSING_DEPEND),)
-include make.depend
else
ifeq ($(wildcard make.depend),)
$(error make.depend not found in $(CURDIR) -- run "make depend" first, or "make veryclean"/"make distclean" to clean without it)
endif
include make.depend
endif
