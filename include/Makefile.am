nobase_include_HEADERS=sse_align.h \
	sse_config.h \
	sse_dslash.h \
	sse_dslash_qdp_packer.h \
	decomp.h \
	decomp_hvv.h \
	mvv_recons.h \
	recons.h

if SCALAR_SINGLE
nobase_include_HEADERS += sse32.h types32.h shift_tables_scalar.h dispatch_scalar.h
endif

if SCALAR_DOUBLE
nobase_include_HEADERS += sse64.h types64.h shift_tables_scalar.h dispatch_scalar.h
endif

if PARSCALAR_SINGLE
nobase_include_HEADERS += sse32.h types32.h shift_tables_parscalar.h dispatch_parscalar.h
endif

if PARSCALAR_DOUBLE
nobase_include_HEADERS += sse64.h types64.h shift_tables_parscalar.h  dispatch_parscalar.h
endif


nodist_include_HEADERS=sse_config_internal.h
