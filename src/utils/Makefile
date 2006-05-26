# src/utils/Makefile

include ../Make.def

#---------------------------------------------------------------------

include Make.deps

F_SRC = \
	RealKind.f90   \
	parameters.f90 \
	DOUG_utils.f90 \
	globals.f90

C_SRC = \
	mpi_ex.c

OBJS       = ${C_SRC:.c=.o} ${F_SRC:.f90=.o} 
LOCAL_LIB  = utils.a
GLOBAL_LIB = ../../lib/$(DOUG_LIB).a

.PHONY: deps

deps:
	@../updatedeps.py -d

all: lib

lib: $(LOCAL_LIB)($(OBJS))
	$(RANLIB) $(LOCAL_LIB)
	@echo Lib $(LOCAL_LIB) is now up-to-date
#	rm -rf $(OBJS)

globalLib: $(GLOBAL_LIB)($(OBJS))
# not to waste time on this - main Makefile will do it on the whole library 
#	$(RANLIB) $(GLOBAL_LIB) 
#	rm -rf *.o

DOUG_utils.o: mpi_ex.o

%.o: %.f90
	$(PURIFY) $(F95) -cpp -DDVMAJOR=$(DVMAJOR) -DDVMINOR=$(DVMINOR) \
	$(FCOMPFLAGS) -c $<

%.o: %.c
	$(PURIFY) $(CC) $(COPTFLAGS) $(INCL) -c $<

%.a : %.f90
	$(PURIFY) $(F95) $(FCOMPFLAGS) -c $<
	$(AR) $(ARFLAGS) $@ $*.o
	rm -f $*.o

%.a : %.c
	$(PURIFY) $(CC) $(COPTFLAGS) $(INCL) -c $<
	$(AR) $(ARFLAGS) $@ $*.o
	rm -f $*.o

.PHONY: clean
clean:
	rm -f *.a *.o *.mod core* *~
