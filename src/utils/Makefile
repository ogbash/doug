# src/uitls/Makefile

include ../Make.def

#---------------------------------------------------------------------

include Make.deps

SRC = \
	RealKind.f90   \
	parameters.f90 \
	DOUG_utils.f90 \
	globals.f90

OBJS       = ${SRC:.f90=.o}
LOCAL_LIB  = utils.a
GLOBAL_LIB = ../../lib/$(DOUG_LIB).a

.PHONY: deps

deps:
	@../updatedeps.py -d

lib: $(LOCAL_LIB)($(OBJS))
	$(RANLIB) $(LOCAL_LIB)
	@echo Lib $(LOCAL_LIB) is now up-to-date
#	rm -rf $(OBJS)

globalLib: $(GLOBAL_LIB)($(OBJS))
# not to waste time on this - main Makefile will do it on the whole library 
#	$(RANLIB) $(GLOBAL_LIB) 
#	rm -rf *.o

%.o: %.f90
	$(PURIFY) $(F95) -cpp -DDVMAJOR=$(DVMAJOR) -DDVMINOR=$(DVMINOR) \
	$(FCOMPFLAGS) -c $<

%.a : %.f90
	$(PURIFY) $(F95) $(FCOMPFLAGS) -c $<
	$(AR) $(ARFLAGS) $@ $*.o
	rm -f $*.o

.PHONY: clean
clean:
	rm -f *.a *.o *.mod core* *~
