#include "decomp.h"
#include "sse_align.h"

#include <xmmintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef  union { 
  unsigned int a[4];
  __m128d vector;
} SSEMask;

void decomp_gamma0_minus(const spinor_array src, halfspinor_array dst) 
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;

  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  
  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0, 0x0}};

  /* Projection: Munge components 0 & 3 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  /* Shuffle the spinor components */
  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Components 1 & 2 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);
   
}


void decomp_gamma1_minus(const spinor_array src, halfspinor_array dst)
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;

  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  
  /* Projection: Munge components 0 & 3 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  /* Shuffle the spinor components */

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Components 1 & 2 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);

}

void decomp_gamma2_minus(const spinor_array src, halfspinor_array dst) 
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;

  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  
  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0, 0x0 }};

  /* Projection: Munge components 0 & 2 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  /* Shuffle the spinor components */
  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Components 1 & 3 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);

}

void decomp_gamma3_minus(const spinor_array src, halfspinor_array dst) 
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;

  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  
   /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);

}



void decomp_gamma0_plus(const spinor_array src, halfspinor_array dst) 
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;

  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  
  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0, 0x0 }};

  /* Projection: Munge components 0 & 3 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  /* Shuffle the spinor components */
  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Components 1 & 2 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);


}

void decomp_gamma1_plus(const spinor_array src, halfspinor_array dst) 
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;

  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  
  /* Projection: Munge components 0 & 3 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  /* Shuffle the spinor components */

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Components 1 & 2 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);


}

void decomp_gamma2_plus(const spinor_array src, halfspinor_array dst) 
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;

  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  
  SSEMask sse_sgn = {{0x0, 0x80000000, 0x0, 0x0 }};

  /* Projection: Munge components 0 & 2 */

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  /* Shuffle the spinor components */
  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Components 1 & 3 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  xmm3 = _mm_shuffle_pd( xmm3, xmm3, 0x1);
  xmm4 = _mm_shuffle_pd( xmm4, xmm4, 0x1);
  xmm5 = _mm_shuffle_pd( xmm5, xmm5, 0x1);
  
  xmm3 = _mm_xor_pd(sse_sgn.vector, xmm3);
  xmm4 = _mm_xor_pd(sse_sgn.vector, xmm4);
  xmm5 = _mm_xor_pd(sse_sgn.vector, xmm5);

  xmm0 = _mm_sub_pd(xmm0, xmm3);
  xmm1 = _mm_sub_pd(xmm1, xmm4);
  xmm2 = _mm_sub_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);

}

void decomp_gamma3_plus(const spinor_array src, halfspinor_array dst) 
{
  __m128d xmm0 ALIGN;
  __m128d xmm1 ALIGN;
  __m128d xmm2 ALIGN;

  __m128d xmm3 ALIGN;
  __m128d xmm4 ALIGN;
  __m128d xmm5 ALIGN;
  
   /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[0][0][0]);
  xmm1 = _mm_load_pd(&src[0][1][0]);
  xmm2 = _mm_load_pd(&src[0][2][0]);

  xmm3 = _mm_load_pd(&src[2][0][0]);
  xmm4 = _mm_load_pd(&src[2][1][0]);
  xmm5 = _mm_load_pd(&src[2][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[0][0][0], xmm0);
  _mm_store_pd(&dst[0][1][0], xmm1);
  _mm_store_pd(&dst[0][2][0], xmm2);

  /* Load spinor component 0 */
  xmm0 = _mm_load_pd(&src[1][0][0]);
  xmm1 = _mm_load_pd(&src[1][1][0]);
  xmm2 = _mm_load_pd(&src[1][2][0]);

  xmm3 = _mm_load_pd(&src[3][0][0]);
  xmm4 = _mm_load_pd(&src[3][1][0]);
  xmm5 = _mm_load_pd(&src[3][2][0]);

  xmm0 = _mm_add_pd(xmm0, xmm3);
  xmm1 = _mm_add_pd(xmm1, xmm4);
  xmm2 = _mm_add_pd(xmm2, xmm5);

  /* Store component */
  _mm_store_pd(&dst[1][0][0], xmm0);
  _mm_store_pd(&dst[1][1][0], xmm1);
  _mm_store_pd(&dst[1][2][0], xmm2);


}


#ifdef __cplusplus
};
#endif