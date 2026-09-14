# Makefile for DAVID_RCI

include ../../make.inc

VPATH = @srcdir@

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG) ../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG).

DAVID_RCI = \
david_rci.o 

all : libdavid_rci.a


libdavid_rci.a:  $(DAVID_RCI)
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


