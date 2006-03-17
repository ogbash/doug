# src/datatypes/Makefile

include ../Make.def

#---------------------------------------------------------------

include Make.deps

INCL = -I../utils/ -IElemMtxs/

SRC = \
	Points2D.f90 \
	Polygon.f90  \
	Graph.f90    \
	Aggregate_mod.f90 \
	Mesh.f90     \
	DenseMtx.f90 \
	Vect.f90     \
	Hash.f90     \
	BHeap.f90    \
	IdxMap.f90

#DATATYPES      = ${SRC:.f90=}
OBJS = ${SRC:.f90=.o}
LOCAL_LIB   = datatypes.a
GLOBAL_LIB  = ../../lib/$(DOUG_LIB).a

.PHONY: all utils deps
all: lib
utils:
	cd ../utils  && $(MAKE) lib
deps:
	@../updatedeps.py -d
	cd SpMtx/    && $(MAKE) deps
	cd ElemMtxs/ && $(MAKE) deps

lib: utils $(LOCAL_LIB)($(OBJS))
	$(RANLIB) $(LOCAL_LIB)
	cd SpMtx/    && $(MAKE) lib
	cd ElemMtxs/ && $(MAKE) lib
#	rm -f *.o

globalLib: $(GLOBAL_LIB)($(OBJS))
# not to waste time on this - main Makefile will do it on the whole library 
#	$(RANLIB) $(GLOBAL_LIB)
	cd SpMtx/    && $(MAKE) globalLib
	cd ElemMtxs/ && $(MAKE) globalLib
#	rm -f *.o

all_drivers:
	cd drivers; make drivers

##%.o: %.f90 utils

%.o: %.f90 
	$(PURIFY) $(F95) $(FCOMPFLAGS) $(INCL) -c $<


.PHONY: clean clean_all
clean:
	rm -f *.a *.o *.mod core* *~
	cd SpMtx/    && $(MAKE) clean
	cd ElemMtxs/ && $(MAKE) clean

clean_all:
	rm -f *.a *.o *.mod core* *~
	cd SpMtx/    && $(MAKE) clean	
	cd ElemMtxs/ && $(MAKE) clean	
	cd drivers/  && $(MAKE) clean	
