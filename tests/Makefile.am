INCFLAGS=-I@top_srcdir@/include -I@top_builddir@/include

AM_CXXFLAGS=$(INCFLAGS) @SSEXX_AUX_CFLAGS@ @QMT_CFLAGS@
AM_LDFLAGS=-L@top_builddir@/lib @LDFLAGS@ @SSEXX_AUX_LDFLAGS@ @QMT_LDFLAGS@
LDADD = -llevel3 @SSEXX_AUX_LIBS@ @LIBS@ @QMT_LIBS@

check_PROGRAMS = t_dslash 

if BUILD_SITE_TESTS
check_PROGRAMS += t_site_ops
endif

TESTS=
if RUN_TESTS
TESTS += $(check_PROGRAMS)
endif


FRAMEWORK_SRCS =  unittest.h \
	dslashm_w.h \
	dslashm_w.cc \
	reunit.h \
	reunit.cc
	testDslashFull.h \
	testDslashFull.cc

t_dslash_SOURCES = $(FRAMEWORK_SRCS) \
	testDslashFull.h \
	testDslashFull.cc \
	t_dslash.cc


if SCALAR_SINGLE
t_site_ops_SOURCES = $(FRAMEWORK_SRCS) \
	t_site_ops_scalar.cc \
	testSiteDslash.h \
	sse32.h \
	testSiteDslash32bit.cc
endif

if SCALAR_DOUBLE
t_site_ops_SOURCES = $(FRAMEWORK_SRCS) \
	t_site_ops_scalar.cc \
	testSiteDslash.h \
	sse64.h \
	testSiteDslash64bit.cc
endif


if PARSCALAR_SINGLE
t_site_ops_SOURCES = $(FRAMEWORK_SRCS) \
	t_site_ops_parscalar.cc \
	sse32.h \
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
t_site_ops_SOURCES = $(FRAMEWORK_SRCS) \
	t_site_ops_parscalar.cc \
	sse64.h \
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
t_site_ops_DEPENDENCIES = @top_builddir@/lib/liblevel3.a
