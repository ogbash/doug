# src/solvers/Makefile

include ../Make.def

# --------------------------------------------

include Make.deps

INCL = -I../utils -I../datatypes $(INCL_UMFPACK4) $(INCL_MUMPS)

F_SRC = \
	solvers.f90 \
	ConvInf.f90 \
	Fact.f90    \
	cg.f90      \
	pcg.f90     \
        subsolvers.f90

C_SRC = \
        umfpack4wrap.c

#SOLVERS      = ${SRC:.f90=}
F_OBJS = ${F_SRC:.f90=.o} 
C_OBJS = ${C_SRC:.c=.o}
LOCAL_LIB   = solvers.a
GLOBAL_LIB  = ../../lib/$(DOUG_LIB).a

.PHONY: all utils datatypes deps
all: lib
utils:
	cd ../utils  && $(MAKE) lib
datatypes:
	cd ../datatypes  && $(MAKE) lib
deps:
	@../updatedeps.py -d

Fact.o: umfpack4wrap.o

lib: utils datatypes $(LOCAL_LIB)($(F_OBJS) $(C_OBJS))
	$(RANLIB) $(LOCAL_LIB)
# TEMORARY SOLUTION>
# move subsolvers.o with its undefined "umf4*" symbols to the 
# end of the archive; umfpack4wrap.o has the symbols in "T" area.
	ar m $(LOCAL_LIB) subsolvers.o
#	@echo "f: " $(F_OBJS) "c: " $(C_OBJS)
#	rm -f *.o

globalLib: $(GLOBAL_LIB)($(F_OBJS) $(C_OBJS))
# TEMORARY SOLUTION>
# move subsolvers.o with its undefined "umf4*" symbols to the 
# end of the archive; umfpack4wrap.o has them in "T" area.
	ar m $(GLOBAL_LIB) subsolvers.o
#	@echo "f: " $(F_OBJS) "c: " $(C_OBJS)
# not to waste time on this - main Makefile will do it on the whole library 
#	$(RANLIB) $(GLOBAL_LIB)
#	rm -f *.o

all_drivers:
	cd drivers; make drivers

##%.o: %.f90 utils

%.o: %.f90
	$(PURIFY) $(F95) $(FCOMPFLAGS) $(INCL) -c $<

%.o: %.c
	$(PURIFY) $(CC) $(COPTFLAGS) $(INCL) -c $<

.PHONY: clean clean_all
clean:
	rm -f *.a *.o *.mod core* *~

clean_all: clean
	cd drivers/ && $(MAKE) clean
