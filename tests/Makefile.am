INCFLAGS=-I@top_srcdir@/include -I@top_builddir@/include

AM_CXXFLAGS=$(INCFLAGS) @SSEXX_AUX_CFLAGS@ @QMT_CFLAGS@
AM_LDFLAGS=-L@top_builddir@/lib @LDFLAGS@ @SSEXX_AUX_LDFLAGS@ @QMT_LDFLAGS@
LDADD = -llevel3 @SSEXX_AUX_LIBS@ @LIBS@ @QMT_LIBS@

check_PROGRAMS = 

if BUILD_QDP_TESTS
check_PROGRAMS += t_dslash
endif

TESTS=$(check_PROGRAMS)

t_dslash_SOURCES = unittest.h \
	dslashm_w.h \
	dslashm_w.cc \
	reunit.h \
	reunit.cc \
	testDslashFull.h \
	testDslashFull.cc

if SCALAR_SINGLE
t_dslash_SOURCES += t_dslash.cc
endif

if SCALAR_DOUBLE
t_dslash_SOURCES += t_dslash.cc
endif


if PARSCALAR_SINGLE
t_dslash_SOURCES += t_dslash_parscalar.cc \
	testDecomp.h \
	testDecompHvv.h \
	testMvvRecons.h \
	testRecons.h \
	testDecomp32bit.cc \
	testDecompHvv32bit.cc \
	testMvvRecons32bit.cc \
	testRecons32bit.cc
endif

if PARSCALAR_DOUBLE
t_dslash_SOURCES += \
	t_dslash_parscalar.cc \
	testDecomp.h \
	testDecompHvv.h \
	testMvvRecons.h \
	testRecons.h \
	testDecomp64bit.cc \
	testDecompHvv64bit.cc \
	testMvvRecons64bit.cc \
	testRecons64bit.cc
endif


t_dslash_DEPENDENCIES = @top_builddir@/lib/liblevel3.a
