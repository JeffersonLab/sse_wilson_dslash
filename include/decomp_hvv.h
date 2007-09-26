#ifndef DECOMP_HVV
#define DECOMP_HVV

#include "sse_config.h"

#if SSE_PRECISION == 32
#include "types32.h"
#else 
#include "types64.h"
#endif

#ifdef __cplusplus
extern "C" { 
#endif

void decomp_hvv_gamma0_plus(const spinor_array src, 
			    const u_mat_array u,
			    halfspinor_array dst);

#ifdef __cplusplus
};
#endif

#endif
