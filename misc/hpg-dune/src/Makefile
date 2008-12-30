
#ALUGRID_PREFIX=~/DUNE/ALUGrid-1.12
#ALUGRID_PREFIX=/usr/local/mpi/gfortran64/openmpi-1.2.5
ALUGRID_PREFIX=/usr/local/src/ALUGrid-1.12
# ALUGRID_PREFIX=/usr/local/

DEBUG=-g
CPPFLAGS+=-DENABLE_ALUGRID
LIBS+=-ldunegrid -ldunecommon -lalugrid
CXXFLAGS+=$(DEBUG)
LDFLAGS+=$(DEBUG)

CPPFLAGS+=-I$(ALUGRID_PREFIX)/include/serial -I$(ALUGRID_PREFIX)/include 
LDFLAGS+=-L$(ALUGRID_PREFIX)/lib

all: gen

gen.o: accmatrix.hh  functors.hh  integrateentity.hh genio.hh info.hh

gen: gen.o
	g++ -o gen gen.o ${CXXFLAGS} ${CPPFLAGS} ${LDFLAGS} ${LIBS}

clean:
	rm -f gen
	rm -f *.o

