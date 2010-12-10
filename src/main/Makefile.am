
@MAKE_INCLUDE@ Make.deps

# doug main executables
bin_PROGRAMS = doug_geom doug_aggr

doug_geom_SOURCES = \
	main_drivers.F90 \
	DOUG.f90 \
	geom.F90

doug_geom_FCFLAGS=-I.. @FCFLAGS_F95@
doug_geom_LDADD=../libdoug.la

doug_aggr_SOURCES = \
	main_drivers.F90 \
	DOUG.f90 \
	aggr.F90

doug_aggr_FCFLAGS=-I.. @FCFLAGS_F95@
doug_aggr_LDADD=../libdoug.la

Make.deps:
	python @top_srcdir@/scripts/amdeps.py -i../Makefile --srcdir=@srcdir@

mostlyclean-local:
	rm -f *.mod
	rm -f *.obj *.o
	rm -f Make.deps

maintainer-clean-local:
	rm -f Make.deps
