# src/datatypes/SpMtx/Makefile

include ../../Make.def

#---------------------------------------------------------------

include Make.deps

INCL = -I../../utils -I.. -I../ElemMtxs

SRC = \
	SpMtx_util.f90        \
	SpMtx_operation.f90   \
	SpMtx_generator.f90   \
	SpMtx_arrangement.f90 \
	SpMtx_permutation.f90 \
	SpMtx_op_AB.f90       \
	SpMtx_op_Ax.f90       \
	SpMtx_class.f90       \
	SpMtx.f90             \
        CoarseMtx.f90         \
        CoarseGrid.f90        \
        CreateRestrict.f90    \
        CreateCoarse.f90      \
        TransmitCoarse.f90    \
        CoarseAllgathers.f90  \
        GeomInterp.f90

OBJS = ${SRC:.f90=.o}
LOCAL_LIB   = ../datatypes.a
GLOBAL_LIB  = ../../../lib/$(DOUG_LIB).a

.PHONY: all utils deps
all: lib
utils:
	cd ../../utils && $(MAKE) lib
deps:
	@../../updatedeps.py -d

lib : utils $(LOCAL_LIB)($(OBJS))
	$(RANLIB) $(LOCAL_LIB)
	cd .. ; ln -sf SpMtx/*.mod . ; cd SpMtx

##	mv *.mod ../
#	rm -f *.o ; mv *.mod ../

globalLib : $(GLOBAL_LIB)($(OBJS))
# not to waste time on this - main Makefile will do it on the whole library 
#	$(RANLIB) $(GLOBAL_LIB)
	cd .. ; ln -sf SpMtx/*.mod . ; cd SpMtx

##	mv *.mod ../

#	rm -f *.o ; mv -f *.mod ../


SpMtx: utils SpMtx.o 
	$(F95) $(FCOMPFLAGS) $(INCL) -c $@.f90 $(LIBS)
SpMtx_arrangement: utils SpMtx_arrangement.o
	$(F95) $(FCOMPFLAGS) $(INCL) -c $@.f90 $(LIBS)
# ...

##%.o : utils %.f90

%.o : %.f90
	$(PURIFY) $(F95) $(FCOMPFLAGS) $(INCL) -c $*.f90


.PHONY: clean 
clean:
	rm -f *.a *.o *.mod core* *~ 

