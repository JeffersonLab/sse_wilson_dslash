#include "decomp.h"
#include "sse_align.h"

#include <xmmintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef union { 
    unsigned int a[4];
    __m128 vector;
  } SSESign;

  static SSESign signs24 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  static SSESign signs34 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  static SSESign signs23 __attribute__((unused)) ALIGN = {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};


  static SSESign signs13 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  static SSESign signs12 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  static SSESign signs14 __attribute__((unused)) ALIGN = {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};


void decomp_gamma0_minus( spinor_array src, halfspinor_array dst) 
{

  /* c <-> color, s <-> spin */

  /* Space for upper components */
  __m128 c0_s01;
  __m128 c1_s01;
  __m128 c2_s01;

  /* Space for lower components */
  __m128 c0_s23;
  __m128 c1_s23;
  __m128 c2_s23;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 c0_s32;
  __m128 c1_s32;
  __m128 c2_s32;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 ic0_s32;
  __m128 ic1_s32;
  __m128 ic2_s32;


  /* Load up the spinors */
  c0_s01 = _mm_loadl_pi(c0_s01, (__m64 *)&src[0][0][0]);
  c1_s01 = _mm_loadl_pi(c1_s01, (__m64 *)&src[0][1][0]);
  c2_s01 = _mm_loadl_pi(c2_s01, (__m64 *)&src[0][2][0]);

  c0_s01 = _mm_loadh_pi(c0_s01, (__m64 *)&src[1][0][0]);
  c1_s01 = _mm_loadh_pi(c1_s01, (__m64 *)&src[1][1][0]);
  c2_s01 = _mm_loadh_pi(c2_s01, (__m64 *)&src[1][2][0]);

  c0_s23 = _mm_loadl_pi(c0_s23, (__m64 *)&src[2][0][0]);
  c1_s23 = _mm_loadl_pi(c1_s23, (__m64 *)&src[2][1][0]);
  c2_s23 = _mm_loadl_pi(c2_s23, (__m64 *)&src[2][2][0]);

  c0_s23 = _mm_loadh_pi(c0_s23, (__m64 *)&src[3][0][0]);
  c1_s23 = _mm_loadh_pi(c1_s23, (__m64 *)&src[3][1][0]);
  c2_s23 = _mm_loadh_pi(c2_s23, (__m64 *)&src[3][2][0]);

 
  /* Swap the lower components  and multiply by -i*/
  c0_s32 = _mm_shuffle_ps(c0_s23, c0_s23, 0x1b);
  c1_s32 = _mm_shuffle_ps(c1_s23, c1_s23, 0x1b);
  c2_s32 = _mm_shuffle_ps(c2_s23, c2_s23, 0x1b);

  ic0_s32 = _mm_xor_ps(c0_s32, signs24.vector);
  ic1_s32 = _mm_xor_ps(c1_s32, signs24.vector);
  ic2_s32 = _mm_xor_ps(c2_s32, signs24.vector);

  /* Add */
  c0_s01 = _mm_add_ps(c0_s01, ic0_s32);
  c1_s01 = _mm_add_ps(c1_s01, ic1_s32);
  c2_s01 = _mm_add_ps(c2_s01, ic2_s32);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);
   
}


void decomp_gamma1_minus( spinor_array src, halfspinor_array dst)
{
  /* Space for upper components */
  __m128 c0_s01;
  __m128 c1_s01;
  __m128 c2_s01;

  /* Space for lower components */
  __m128 c0_s23;
  __m128 c1_s23;
  __m128 c2_s23;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 c0_s32;
  __m128 c1_s32;
  __m128 c2_s32;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 sc0_s32;
  __m128 sc1_s32;
  __m128 sc2_s32;


  /* Load up the spinors */
  c0_s01 = _mm_loadl_pi(c0_s01, (__m64 *)&src[0][0][0]);
  c1_s01 = _mm_loadl_pi(c1_s01, (__m64 *)&src[0][1][0]);
  c2_s01 = _mm_loadl_pi(c2_s01, (__m64 *)&src[0][2][0]);

  c0_s01 = _mm_loadh_pi(c0_s01, (__m64 *)&src[1][0][0]);
  c1_s01 = _mm_loadh_pi(c1_s01, (__m64 *)&src[1][1][0]);
  c2_s01 = _mm_loadh_pi(c2_s01, (__m64 *)&src[1][2][0]);

  c0_s23 = _mm_loadl_pi(c0_s23, (__m64 *)&src[2][0][0]);
  c1_s23 = _mm_loadl_pi(c1_s23, (__m64 *)&src[2][1][0]);
  c2_s23 = _mm_loadl_pi(c2_s23, (__m64 *)&src[2][2][0]);

  c0_s23 = _mm_loadh_pi(c0_s23, (__m64 *)&src[3][0][0]);
  c1_s23 = _mm_loadh_pi(c1_s23, (__m64 *)&src[3][1][0]);
  c2_s23 = _mm_loadh_pi(c2_s23, (__m64 *)&src[3][2][0]);

 
  /* Swap the lower components */
  c0_s32 = _mm_shuffle_ps(c0_s23, c0_s23, 0x4e);
  c1_s32 = _mm_shuffle_ps(c1_s23, c1_s23, 0x4e);
  c2_s32 = _mm_shuffle_ps(c2_s23, c2_s23, 0x4e);

  sc0_s32 = _mm_xor_ps(c0_s32, signs34.vector);
  sc1_s32 = _mm_xor_ps(c1_s32, signs34.vector);
  sc2_s32 = _mm_xor_ps(c2_s32, signs34.vector);

  /* Add */
  c0_s01 = _mm_add_ps(c0_s01, sc0_s32);
  c1_s01 = _mm_add_ps(c1_s01, sc1_s32);
  c2_s01 = _mm_add_ps(c2_s01, sc2_s32);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);

  

}

void decomp_gamma2_minus( spinor_array src, halfspinor_array dst) 
{

  /* Space for upper components */
  __m128 c0_s01;
  __m128 c1_s01;
  __m128 c2_s01;

  /* Space for lower components */
  __m128 c0_s23;
  __m128 c1_s23;
  __m128 c2_s23;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 c0_s32;
  __m128 c1_s32;
  __m128 c2_s32;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 sc0_s32;
  __m128 sc1_s32;
  __m128 sc2_s32;


  /* Load up the spinors */
  c0_s01 = _mm_loadl_pi(c0_s01, (__m64 *)&src[0][0][0]);
  c1_s01 = _mm_loadl_pi(c1_s01, (__m64 *)&src[0][1][0]);
  c2_s01 = _mm_loadl_pi(c2_s01, (__m64 *)&src[0][2][0]);

  c0_s01 = _mm_loadh_pi(c0_s01, (__m64 *)&src[1][0][0]);
  c1_s01 = _mm_loadh_pi(c1_s01, (__m64 *)&src[1][1][0]);
  c2_s01 = _mm_loadh_pi(c2_s01, (__m64 *)&src[1][2][0]);

  c0_s23 = _mm_loadl_pi(c0_s23, (__m64 *)&src[2][0][0]);
  c1_s23 = _mm_loadl_pi(c1_s23, (__m64 *)&src[2][1][0]);
  c2_s23 = _mm_loadl_pi(c2_s23, (__m64 *)&src[2][2][0]);

  c0_s23 = _mm_loadh_pi(c0_s23, (__m64 *)&src[3][0][0]);
  c1_s23 = _mm_loadh_pi(c1_s23, (__m64 *)&src[3][1][0]);
  c2_s23 = _mm_loadh_pi(c2_s23, (__m64 *)&src[3][2][0]);

 
  /* Swap the lower components */
  c0_s32 = _mm_shuffle_ps(c0_s23, c0_s23, 0xb1);
  c1_s32 = _mm_shuffle_ps(c1_s23, c1_s23, 0xb1);
  c2_s32 = _mm_shuffle_ps(c2_s23, c2_s23, 0xb1);

  sc0_s32 = _mm_xor_ps(c0_s32, signs23.vector);
  sc1_s32 = _mm_xor_ps(c1_s32, signs23.vector);
  sc2_s32 = _mm_xor_ps(c2_s32, signs23.vector);

  /* Add */
  c0_s01 = _mm_add_ps(c0_s01, sc0_s32);
  c1_s01 = _mm_add_ps(c1_s01, sc1_s32);
  c2_s01 = _mm_add_ps(c2_s01, sc2_s32);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);

  

}

void decomp_gamma3_minus( spinor_array src, halfspinor_array dst) 
{

  /* Space for upper components */
  __m128 c0_s01;
  __m128 c1_s01;
  __m128 c2_s01;

  /* Space for lower components */
  __m128 c0_s23;
  __m128 c1_s23;
  __m128 c2_s23;

  /* Load up the spinors */
  c0_s01 = _mm_loadl_pi(c0_s01, (__m64 *)&src[0][0][0]);
  c1_s01 = _mm_loadl_pi(c1_s01, (__m64 *)&src[0][1][0]);
  c2_s01 = _mm_loadl_pi(c2_s01, (__m64 *)&src[0][2][0]);

  c0_s01 = _mm_loadh_pi(c0_s01, (__m64 *)&src[1][0][0]);
  c1_s01 = _mm_loadh_pi(c1_s01, (__m64 *)&src[1][1][0]);
  c2_s01 = _mm_loadh_pi(c2_s01, (__m64 *)&src[1][2][0]);

  c0_s23 = _mm_loadl_pi(c0_s23, (__m64 *)&src[2][0][0]);
  c1_s23 = _mm_loadl_pi(c1_s23, (__m64 *)&src[2][1][0]);
  c2_s23 = _mm_loadl_pi(c2_s23, (__m64 *)&src[2][2][0]);

  c0_s23 = _mm_loadh_pi(c0_s23, (__m64 *)&src[3][0][0]);
  c1_s23 = _mm_loadh_pi(c1_s23, (__m64 *)&src[3][1][0]);
  c2_s23 = _mm_loadh_pi(c2_s23, (__m64 *)&src[3][2][0]);

 
  /* sub */
  c0_s01 = _mm_sub_ps(c0_s01, c0_s23);
  c1_s01 = _mm_sub_ps(c1_s01, c1_s23);
  c2_s01 = _mm_sub_ps(c2_s01, c2_s23);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);

}



void decomp_gamma0_plus( spinor_array src, halfspinor_array dst) 
{

  /* c <-> color, s <-> spin */

  /* Space for upper components */
  __m128 c0_s01;
  __m128 c1_s01;
  __m128 c2_s01;

  /* Space for lower components */
  __m128 c0_s23;
  __m128 c1_s23;
  __m128 c2_s23;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 c0_s32;
  __m128 c1_s32;
  __m128 c2_s32;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 ic0_s32;
  __m128 ic1_s32;
  __m128 ic2_s32;


  /* Load up the spinors */
  /* Color 0 */
  c0_s01 = _mm_loadl_pi(c0_s01, (__m64 *)&src[0][0][0]);
  c1_s01 = _mm_loadl_pi(c1_s01, (__m64 *)&src[0][1][0]);
  c2_s01 = _mm_loadl_pi(c2_s01, (__m64 *)&src[0][2][0]);

  c0_s01 = _mm_loadh_pi(c0_s01, (__m64 *)&src[1][0][0]);
  c1_s01 = _mm_loadh_pi(c1_s01, (__m64 *)&src[1][1][0]);
  c2_s01 = _mm_loadh_pi(c2_s01, (__m64 *)&src[1][2][0]);

  c0_s23 = _mm_loadl_pi(c0_s23, (__m64 *)&src[2][0][0]);
  c1_s23 = _mm_loadl_pi(c1_s23, (__m64 *)&src[2][1][0]);
  c2_s23 = _mm_loadl_pi(c2_s23, (__m64 *)&src[2][2][0]);

  c0_s23 = _mm_loadh_pi(c0_s23, (__m64 *)&src[3][0][0]);
  c1_s23 = _mm_loadh_pi(c1_s23, (__m64 *)&src[3][1][0]);
  c2_s23 = _mm_loadh_pi(c2_s23, (__m64 *)&src[3][2][0]);

 
  /* Swap the lower components  and multiply by +i*/
  c0_s32 = _mm_shuffle_ps(c0_s23, c0_s23, 0x1b);
  c1_s32 = _mm_shuffle_ps(c1_s23, c1_s23, 0x1b);
  c2_s32 = _mm_shuffle_ps(c2_s23, c2_s23, 0x1b);

  ic0_s32 = _mm_xor_ps(c0_s32, signs13.vector);
  ic1_s32 = _mm_xor_ps(c1_s32, signs13.vector);
  ic2_s32 = _mm_xor_ps(c2_s32, signs13.vector);

  /* Add */
  c0_s01 = _mm_add_ps(c0_s01, ic0_s32);
  c1_s01 = _mm_add_ps(c1_s01, ic1_s32);
  c2_s01 = _mm_add_ps(c2_s01, ic2_s32);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);
   

}

void decomp_gamma1_plus( spinor_array src, halfspinor_array dst) 
{
  /* Space for upper components */
  __m128 c0_s01;
  __m128 c1_s01;
  __m128 c2_s01;

  /* Space for lower components */
  __m128 c0_s23;
  __m128 c1_s23;
  __m128 c2_s23;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 c0_s32;
  __m128 c1_s32;
  __m128 c2_s32;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 sc0_s32;
  __m128 sc1_s32;
  __m128 sc2_s32;


  /* Load up the spinors */
  c0_s01 = _mm_loadl_pi(c0_s01, (__m64 *)&src[0][0][0]);
  c1_s01 = _mm_loadl_pi(c1_s01, (__m64 *)&src[0][1][0]);
  c2_s01 = _mm_loadl_pi(c2_s01, (__m64 *)&src[0][2][0]);

  c0_s01 = _mm_loadh_pi(c0_s01, (__m64 *)&src[1][0][0]);
  c1_s01 = _mm_loadh_pi(c1_s01, (__m64 *)&src[1][1][0]);
  c2_s01 = _mm_loadh_pi(c2_s01, (__m64 *)&src[1][2][0]);

  c0_s23 = _mm_loadl_pi(c0_s23, (__m64 *)&src[2][0][0]);
  c1_s23 = _mm_loadl_pi(c1_s23, (__m64 *)&src[2][1][0]);
  c2_s23 = _mm_loadl_pi(c2_s23, (__m64 *)&src[2][2][0]);

  c0_s23 = _mm_loadh_pi(c0_s23, (__m64 *)&src[3][0][0]);
  c1_s23 = _mm_loadh_pi(c1_s23, (__m64 *)&src[3][1][0]);
  c2_s23 = _mm_loadh_pi(c2_s23, (__m64 *)&src[3][2][0]);

 
  /* Swap the lower components */
  c0_s32 = _mm_shuffle_ps(c0_s23, c0_s23, 0x4e);
  c1_s32 = _mm_shuffle_ps(c1_s23, c1_s23, 0x4e);
  c2_s32 = _mm_shuffle_ps(c2_s23, c2_s23, 0x4e);

  sc0_s32 = _mm_xor_ps(c0_s32, signs12.vector);
  sc1_s32 = _mm_xor_ps(c1_s32, signs12.vector);
  sc2_s32 = _mm_xor_ps(c2_s32, signs12.vector);

  /* Add */
  c0_s01 = _mm_add_ps(c0_s01, sc0_s32);
  c1_s01 = _mm_add_ps(c1_s01, sc1_s32);
  c2_s01 = _mm_add_ps(c2_s01, sc2_s32);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);


}

void decomp_gamma2_plus( spinor_array src, halfspinor_array dst) 
{
  /* Space for upper components */
  __m128 c0_s01;
  __m128 c1_s01;
  __m128 c2_s01;

  /* Space for lower components */
  __m128 c0_s23;
  __m128 c1_s23;
  __m128 c2_s23;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 c0_s32;
  __m128 c1_s32;
  __m128 c2_s32;

  /* Swap upper and lower components */
  /* Compiler should spill, or use 64 bit extras */
  __m128 sc0_s32;
  __m128 sc1_s32;
  __m128 sc2_s32;


  /* Load up the spinors */
  c0_s01 = _mm_loadl_pi(c0_s01, (__m64 *)&src[0][0][0]);
  c1_s01 = _mm_loadl_pi(c1_s01, (__m64 *)&src[0][1][0]);
  c2_s01 = _mm_loadl_pi(c2_s01, (__m64 *)&src[0][2][0]);

  c0_s01 = _mm_loadh_pi(c0_s01, (__m64 *)&src[1][0][0]);
  c1_s01 = _mm_loadh_pi(c1_s01, (__m64 *)&src[1][1][0]);
  c2_s01 = _mm_loadh_pi(c2_s01, (__m64 *)&src[1][2][0]);

  c0_s23 = _mm_loadl_pi(c0_s23, (__m64 *)&src[2][0][0]);
  c1_s23 = _mm_loadl_pi(c1_s23, (__m64 *)&src[2][1][0]);
  c2_s23 = _mm_loadl_pi(c2_s23, (__m64 *)&src[2][2][0]);

  c0_s23 = _mm_loadh_pi(c0_s23, (__m64 *)&src[3][0][0]);
  c1_s23 = _mm_loadh_pi(c1_s23, (__m64 *)&src[3][1][0]);
  c2_s23 = _mm_loadh_pi(c2_s23, (__m64 *)&src[3][2][0]);

 
  /* Swap the lower components */
  c0_s32 = _mm_shuffle_ps(c0_s23, c0_s23, 0xb1);
  c1_s32 = _mm_shuffle_ps(c1_s23, c1_s23, 0xb1);
  c2_s32 = _mm_shuffle_ps(c2_s23, c2_s23, 0xb1);

  sc0_s32 = _mm_xor_ps(c0_s32, signs14.vector);
  sc1_s32 = _mm_xor_ps(c1_s32, signs14.vector);
  sc2_s32 = _mm_xor_ps(c2_s32, signs14.vector);

  /* Add */
  c0_s01 = _mm_add_ps(c0_s01, sc0_s32);
  c1_s01 = _mm_add_ps(c1_s01, sc1_s32);
  c2_s01 = _mm_add_ps(c2_s01, sc2_s32);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);


}

void decomp_gamma3_plus( spinor_array src, halfspinor_array dst) 
{
  /* Space for upper components */
  __m128 c0_s01;
  __m128 c1_s01;
  __m128 c2_s01;

  /* Space for lower components */
  __m128 c0_s23;
  __m128 c1_s23;
  __m128 c2_s23;

  /* Load up the spinors */
  c0_s01 = _mm_loadl_pi(c0_s01, (__m64 *)&src[0][0][0]);
  c1_s01 = _mm_loadl_pi(c1_s01, (__m64 *)&src[0][1][0]);
  c2_s01 = _mm_loadl_pi(c2_s01, (__m64 *)&src[0][2][0]);

  c0_s01 = _mm_loadh_pi(c0_s01, (__m64 *)&src[1][0][0]);
  c1_s01 = _mm_loadh_pi(c1_s01, (__m64 *)&src[1][1][0]);
  c2_s01 = _mm_loadh_pi(c2_s01, (__m64 *)&src[1][2][0]);

  c0_s23 = _mm_loadl_pi(c0_s23, (__m64 *)&src[2][0][0]);
  c1_s23 = _mm_loadl_pi(c1_s23, (__m64 *)&src[2][1][0]);
  c2_s23 = _mm_loadl_pi(c2_s23, (__m64 *)&src[2][2][0]);

  c0_s23 = _mm_loadh_pi(c0_s23, (__m64 *)&src[3][0][0]);
  c1_s23 = _mm_loadh_pi(c1_s23, (__m64 *)&src[3][1][0]);
  c2_s23 = _mm_loadh_pi(c2_s23, (__m64 *)&src[3][2][0]);

 
  /* sub */
  c0_s01 = _mm_add_ps(c0_s01, c0_s23);
  c1_s01 = _mm_add_ps(c1_s01, c1_s23);
  c2_s01 = _mm_add_ps(c2_s01, c2_s23);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);


}


#ifdef __cplusplus
};
#endif
