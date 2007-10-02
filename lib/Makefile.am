# $Id: Makefile.am,v 1.6 2007-10-02 20:40:21 bjoo Exp $
#
# Path for all files used in building the library
# This file was automatically generated

TOPSRCDIR=@top_srcdir@
#
# Include Path
#
INCFLAGS = -I$(TOPSRCDIR)/include
#
# Any other compiler flags
# 
AM_CFLAGS = $(INCFLAGS) @CFLAGS@ @DEFS@ @SSEXX_AUX_CFLAGS@ @QMT_CFLAGS@
AM_CXXFLAGS = $(INCFLAGS) @SSEXX_AUX_CFLAGS@ @QMT_CFLAGS@
AM_CC  = $(CC)
AM_CXX = $(CXX)

lib_LIBRARIES = liblevel3.a


liblevel3_a_SOURCES = 

#
# An AM_CONDITIONAL if statement to decide which fermion type
#

if SCALAR_SINGLE
liblevel3_a_SOURCES += shift_tables_scalar.c packer_nopad.c sse_su3dslash_32bit_scalar.c dispatch_scalar.c site_dslash_32bit_scalar.c
if BUILD_QDP_PACKERS
liblevel3_a_SOURCES += qdp_packer_nopad.cc
endif
endif

if SCALAR_DOUBLE
liblevel3_a_SOURCES += shift_tables_scalar.c packer_nopad.c sse_su3dslash_64bit_scalar.c dispatch_scalar.c site_dslash_64bit_scalar.c
if BUILD_QDP_PACKERS
liblevel3_a_SOURCES += qdp_packer_nopad.cc
endif
endif

if PARSCALAR_SINGLE
liblevel3_a_SOURCES += shift_tables_parscalar.c \
	packer_p4_pad.c \
	sse_su3dslash_32bit_parscalar.c \
	dispatch_parscalar.c \
	decomp_32bit_parscalar.c \
	decomp_hvv_32bit_parscalar.c \
	mvv_recons_32bit_parscalar.c \
	recons_32bit_parscalar.c

if BUILD_QDP_PACKERS
liblevel3_a_SOURCES += qdp_packer_p4_pad.cc
endif
endif

if PARSCALAR_DOUBLE
liblevel3_a_SOURCES += shift_tables_parscalar.c \
	packer_nopad.c \
	sse_su3dslash_64bit_parscalar.c \
	dispatch_parscalar.c \
	decomp_64bit_parscalar.c \
	decomp_hvv_64bit_parscalar.c \
	mvv_recons_64bit_parscalar.c \
	recons_64bit_parscalar.c

if BUILD_QDP_PACKERS
liblevel3_a_SOURCES += qdp_packer_nopad.cc
endif
endif
