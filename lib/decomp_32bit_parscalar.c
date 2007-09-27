#include "decomp.h"
#include "sse_align.h"

#include <xmmintrin.h>

#ifdef __cplusplus
extern "C" {
#endif

void decomp_gamma0_minus(const spinor_array src, halfspinor_array dst) 
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

  union { 
    float a[4];
    __m128 vector;
  } signs ALIGN = {{1,-1,1,-1}};

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

  ic0_s32 = _mm_mul_ps(c0_s32, signs.vector);
  ic1_s32 = _mm_mul_ps(c1_s32, signs.vector);
  ic2_s32 = _mm_mul_ps(c2_s32, signs.vector);

  /* Add */
  c0_s01 = _mm_add_ps(c0_s01, ic0_s32);
  c1_s01 = _mm_add_ps(c1_s01, ic1_s32);
  c2_s01 = _mm_add_ps(c2_s01, ic2_s32);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);
   
}


void decomp_gamma1_minus(const spinor_array src, halfspinor_array dst)
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

  union { 
    float a[4];
    __m128 vector;
  } signs ALIGN = {{1,1,-1,-1}};

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

  sc0_s32 = _mm_mul_ps(c0_s32, signs.vector);
  sc1_s32 = _mm_mul_ps(c1_s32, signs.vector);
  sc2_s32 = _mm_mul_ps(c2_s32, signs.vector);

  /* Add */
  c0_s01 = _mm_add_ps(c0_s01, sc0_s32);
  c1_s01 = _mm_add_ps(c1_s01, sc1_s32);
  c2_s01 = _mm_add_ps(c2_s01, sc2_s32);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);

  

}

void decomp_gamma2_minus(const spinor_array src, halfspinor_array dst) 
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

  union { 
    float a[4];
    __m128 vector;
  } signs ALIGN = {{1,-1,-1,1}};

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

  sc0_s32 = _mm_mul_ps(c0_s32, signs.vector);
  sc1_s32 = _mm_mul_ps(c1_s32, signs.vector);
  sc2_s32 = _mm_mul_ps(c2_s32, signs.vector);

  /* Add */
  c0_s01 = _mm_add_ps(c0_s01, sc0_s32);
  c1_s01 = _mm_add_ps(c1_s01, sc1_s32);
  c2_s01 = _mm_add_ps(c2_s01, sc2_s32);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);

  

}

void decomp_gamma3_minus(const spinor_array src, halfspinor_array dst) 
{

  /* Space for upper components */
  __m128 c0_s01;
  __m128 c1_s01;
  __m128 c2_s01;

  /* Space for lower components */
  __m128 c0_s23;
  __m128 c1_s23;
  __m128 c2_s23;


  union { 
    float a[4];
    __m128 vector;
  } signs ALIGN = {{1,-1,-1,1}};

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



void decomp_gamma0_plus(const spinor_array src, halfspinor_array dst) 
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

  union { 
    float a[4];
    __m128 vector;
  } signs ALIGN = {{-1,1,-1,1}};

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

  ic0_s32 = _mm_mul_ps(c0_s32, signs.vector);
  ic1_s32 = _mm_mul_ps(c1_s32, signs.vector);
  ic2_s32 = _mm_mul_ps(c2_s32, signs.vector);

  /* Add */
  c0_s01 = _mm_add_ps(c0_s01, ic0_s32);
  c1_s01 = _mm_add_ps(c1_s01, ic1_s32);
  c2_s01 = _mm_add_ps(c2_s01, ic2_s32);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);
   

}

void decomp_gamma1_plus(const spinor_array src, halfspinor_array dst) 
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

  union { 
    float a[4];
    __m128 vector;
  } signs ALIGN = {{-1,-1,1,1}};

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

  sc0_s32 = _mm_mul_ps(c0_s32, signs.vector);
  sc1_s32 = _mm_mul_ps(c1_s32, signs.vector);
  sc2_s32 = _mm_mul_ps(c2_s32, signs.vector);

  /* Add */
  c0_s01 = _mm_add_ps(c0_s01, sc0_s32);
  c1_s01 = _mm_add_ps(c1_s01, sc1_s32);
  c2_s01 = _mm_add_ps(c2_s01, sc2_s32);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);


}

void decomp_gamma2_plus(const spinor_array src, halfspinor_array dst) 
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

  union { 
    float a[4];
    __m128 vector;
  } signs ALIGN = {{-1,1,1,-1}};

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

  sc0_s32 = _mm_mul_ps(c0_s32, signs.vector);
  sc1_s32 = _mm_mul_ps(c1_s32, signs.vector);
  sc2_s32 = _mm_mul_ps(c2_s32, signs.vector);

  /* Add */
  c0_s01 = _mm_add_ps(c0_s01, sc0_s32);
  c1_s01 = _mm_add_ps(c1_s01, sc1_s32);
  c2_s01 = _mm_add_ps(c2_s01, sc2_s32);

  /* Store */
  _mm_store_ps(&dst[0][0][0],c0_s01);
  _mm_store_ps(&dst[1][0][0],c1_s01);
  _mm_store_ps(&dst[2][0][0],c2_s01);


}

void decomp_gamma3_plus(const spinor_array src, halfspinor_array dst) 
{
  /* Space for upper components */
  __m128 c0_s01;
  __m128 c1_s01;
  __m128 c2_s01;

  /* Space for lower components */
  __m128 c0_s23;
  __m128 c1_s23;
  __m128 c2_s23;


  union { 
    float a[4];
    __m128 vector;
  } signs ALIGN = {{1,-1,-1,1}};

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