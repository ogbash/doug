
libdoug_la_SOURCES += \
	utils/mpi_ex.c \
	utils/RealKind.F90 \
	utils/parameters.F90 \
	utils/globals.F90 \
	utils/DOUG_utils.F90

EXTRA_DIST += \
	utils/controls.F90

noinst_PROGRAMS = matrix_gen
matrix_gen_SOURCES = utils/MatrixGen.F90
matrix_gen_LDADD=libdoug.la
