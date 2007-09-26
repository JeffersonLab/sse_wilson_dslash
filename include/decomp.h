#ifndef DECOMP_H
#define DECOMP_H

#include "sse_config.h"

#if SSE_PRECISION == 32
#include "types32.h"
#else 
#include "types64.h"
#endif

#ifdef __cplusplus
extern "C" { 
#endif

void decomp_gamma0_minus(const spinor_array src, halfspinor_array dst);
void decomp_gamma1_minus(const spinor_array src, halfspinor_array dst);
void decomp_gamma2_minus(const spinor_array src, halfspinor_array dst);  
void decomp_gamma3_minus(const spinor_array src, halfspinor_array dst);

void decomp_gamma0_plus(const spinor_array src, halfspinor_array dst);
void decomp_gamma1_plus(const spinor_array src, halfspinor_array dst);
void decomp_gamma2_plus(const spinor_array src, halfspinor_array dst);  
void decomp_gamma3_plus(const spinor_array src, halfspinor_array dst);


#ifdef __cplusplus
};
#endif

#endif
