# src/datatypes/ElemMtxs/Makefile

include ../../Make.def

#---------------------------------------------------------------

include Make.deps

INCL = -I../../utils -I.. -I../SpMtx

SRC = \
	ElemMtxs_assemble.f90    \
	ElemMtxs_distribute.f90  \
	ElemMtxs_base.f90        \
	ElemMtxs.f90

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
	cd .. ; ln -sf ElemMtxs/*.mod . ; cd ElemMtxs

##	mv *.mod ../
#	rm -f *.o ; mv *.mod ../

globalLib : $(GLOBAL_LIB)($(OBJS))
# not to waste time on this - main Makefile will do it on the whole library 
#	$(RANLIB) $(GLOBAL_LIB)
	cd .. ; ln -sf ElemMtxs/*.mod . ; cd ElemMtxs

##	mv *.mod ../

#	rm -f *.o ; mv -f *.mod ../



##%.o : utils %.f90

%.o : %.f90
	$(PURIFY) $(F95) $(FCOMPFLAGS) $(INCL) -c $*.f90


.PHONY: clean 
clean:
	rm -f *.a *.o *.mod core* *~ 

