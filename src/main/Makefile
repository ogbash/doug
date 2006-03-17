# src/main/Makefile

include ../Make.def

#---------------------------------------------------------------

BINDIR = ../../bin
INCL = -I../utils -I../datatypes -I../solvers

LIBS = 	../../lib/$(DOUG_LIB).a \
		$(LIBS_UMFPACK4) \
		$(LIBS_MUMPS) \
		$(LIBS_METIS) \
		$(LIBS_PL) \
		$(LIBS_MPI)

include Make.deps

SRC = \
	main.f90          \
	main_drivers.f90  \
	DOUG.f90

AGGRSRC = \
	aggr.f90          \
	main_drivers.f90  \
	DOUG.f90

OBJS = ${SRC:.f90=.o}

AGGROBJS = ${AGGRSRC:.f90=.o}

.PHONY: all deps
all: aggr

deps:
	@../updatedeps.py -d

lib:
	cd ../../  &&  $(MAKE) lib

main: lib $(OBJS)
	$(F95) $(FCOMPFLAGS) $(INCL) -o $@ $(OBJS) $(LIBS)
	mv -f $@  $(BINDIR)/$(DOUG_EXE) && \
	chmod 755 $(BINDIR)/$(DOUG_EXE)

aggr: lib $(AGGROBJS)
	$(F95) $(FCOMPFLAGS) $(INCL) -o $@ $(AGGROBJS) $(LIBS)
	mv -f $@  $(BINDIR)/$(AGGR_EXE) && \
	chmod 755 $(BINDIR)/$(AGGR_EXE)

%.o: %.f90
	$(F95) $(FCOMPFLAGS) $(INCL) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.a *.o *.mod core* *~ $(BINDIR)/$(DOUG_EXE) \
        $(BINDIR)/$(AGGR_EXE)
