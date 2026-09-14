# Makefile for ParO

include ../../make.inc

VPATH = @srcdir@

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG) ../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG).

PARO = \
paro_gamma.o paro_gamma_new.o \
paro_k.o paro_k_new.o \
pcg_gamma.o \
pcg_k.o \
bpcg_gamma.o \
bpcg_k.o 

all : libparo.a

libparo.a:  $(PARO)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x


ifneq ($(TOLERATE_MISSING_DEPEND),)
-include make.depend
else
ifeq ($(wildcard make.depend),)
$(error make.depend not found in $(CURDIR) -- run "make depend" first, or "make veryclean"/"make distclean" to clean without it)
endif
include make.depend
endif


