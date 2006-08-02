
@MAKE_INCLUDE@ Make.deps

# doug_main
bin_PROGRAMS = doug_main doug_aggr

doug_main_SOURCES = \
	main_drivers.F90 \
	DOUG.f90 \
	main.F90

doug_main_FCFLAGS=-I.. @FCFLAGS_F95@
doug_main_LDADD=../libdoug.la

doug_aggr_SOURCES = \
	main_drivers.F90 \
	DOUG.f90 \
	aggr.F90

doug_aggr_FCFLAGS=-I.. @FCFLAGS_F95@
doug_aggr_LDADD=../libdoug.la

Make.deps:
	python @top_srcdir@/scripts/amdeps.py -i../Makefile

mostlyclean-local:
	rm -f *.mod

maintainer-clean-local:
	rm -f Make.deps
