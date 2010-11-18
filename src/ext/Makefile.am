lib_LTLIBRARIES = libdoug-ext.la
libdoug_ext_la_SOURCES = doug_ext.f90 ops.c
libdoug_ext_la_FCFLAGS=-I@top_builddir@/src -I@top_builddir@/src/main
noinst_PROGRAMS=c_doug_aggr c_doug_cg c_doug_pcg1

c_doug_aggr_SOURCES=main_aggr.c
# this file is compiled with 'mpicc' and in static linking does not find MPI library for Fortran,
#  thus we hack it by providing required libraries, however it only works for OpenMPI
c_doug_aggr_LDADD=libdoug-ext.la ../libdoug.la ../main/doug_aggr-main_drivers.o -lmpi_f90 -lmpi_f77

c_doug_cg_SOURCES=main_cg.c
c_doug_cg_LDADD=libdoug-ext.la ../libdoug.la ../main/doug_aggr-main_drivers.o -lmpi_f90 -lmpi_f77

c_doug_pcg1_SOURCES=main_pcg1.c
c_doug_pcg1_LDADD=libdoug-ext.la ../libdoug.la ../main/doug_aggr-main_drivers.o -lmpi_f90 -lmpi_f77
