
@MAKE_INCLUDE@ Make.deps

lib_LTLIBRARIES = libdoug.la
libdoug_la_SOURCES = 
libdoug_la_FCFLAGS=-I@srcdir@/utils -I@srcdir@/datatypes/SpMtx @FCFLAGS_F95@
EXTRA_DIST = 

include utils/Makefile.am
include datatypes/Makefile.am
include solvers/Makefile.am

SUBDIRS = main

mostlyclean-local:
	rm -f *.mod

maintainer-clean-local:
	rm -f Make.deps

all-recursive: all-am

Make.deps: Makefile.am
	python @top_srcdir@/scripts/amdeps.py --srcdir=@srcdir@
