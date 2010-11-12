lib_LTLIBRARIES = libdoug-ext.la
libdoug_ext_la_SOURCES = DOUG_utils.f90
libdoug_ext_la_FCFLAGS=-I@top_builddir@/src -I@top_builddir@/src/main
noinst_PROGRAMS=c_doug_aggr
c_doug_aggr_SOURCES=main_aggr.c
# this file is compiled with 'mpicc' and in static linking does not find MPI library for Fortran,
#  thus we hack it by providing required libraries, however it only works for OpenMPI
c_doug_aggr_LDADD=libdoug-ext.la ../libdoug.la ../main/doug_aggr-main_drivers.o -lmpi_f90 -lmpi_f77
