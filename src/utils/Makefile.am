
libdoug_la_SOURCES += \
	utils/mpi_ex.c \
	utils/RealKind.F90 \
	utils/parameters.F90 \
	utils/globals.F90 \
	utils/DOUG_utils.F90

EXTRA_DIST += \
	utils/controls.F90

noinst_PROGRAMS = matrix_gen vector_gen txt2xdr
matrix_gen_SOURCES = utils/MatrixGen.F90
matrix_gen_LDADD=libdoug.la
vector_gen_SOURCES = utils/VectorGen.F90
vector_gen_LDADD=libdoug.la
txt2xdr_SOURCES = utils/txt2xdr.F90
txt2xdr_LDADD=libdoug.la
