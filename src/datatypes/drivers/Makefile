# src/datatypes/drivers/Makefile

include ../../Make.def

#---------------------------------------------------------------------

INCL = -I../ -I../../utils
LIBS = $(LIBS_METIS) $(LIBS_PL) $(LIBS_UMFPACK2) $(LIBS_UMFPACK4)

COMMON_OBJS = ../datatypes.a ../../utils/utils.a

DRIVERS = \
	test_Graph    \
	test_Mesh     \
	test_Polygon  \
	test_PipeMesh \
	test_ElemMtxs \
	test_assembleSpMtx

.PHONY: drivers all_deps
drivers : $(DRIVERS) test_SpMtx_symmetry_at_pmvm
all_deps:
	cd ../../utils && ${MAKE} lib
	cd ../         && ${MAKE} lib

%: %.f90 all_deps
	$(F95) $(FCOMPFLAGS) $(INCL) -o $@ $< $(COMMON_OBJS) $(LIBS)

.PHONY: clean
clean:
	rm -f *.a *.o *.mod core* *~ $(DRIVERS) test_SpMtx_symmetry_at_pmvm

#
test_SpMtx_symmetry_at_pmvm: all_deps
	cd ../../main && \
	$(F95) $(FCOMPFLAGS) -I../datatypes -I../utils -c main_drivers.f90
	$(F95) $(FCOMPFLAGS) $(INCL) -I../../main -o $@ $@.f90 \
		../../main/main_drivers.o $(COMMON_OBJS) $(LIBS)
