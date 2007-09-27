#include "mvv_recons.h"

#include "xmmintrin.h"

#ifdef __cplusplus
extern "C" {
#endif

void mvv_recons_gamma0_plus(const halfspinor_array src, 
			    const u_mat_array u,
			    halfspinor_array upper_sum, halfspinor_array lower_sum)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  union { 
  float a[4];
    __m128 vector;
  } signs13 ALIGN = {{-1,1,-1,1}};

  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][2][0]);
  xmm7 = _mm_load_ss(&u[2][0][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_mul_ps(signs13.vector, xmm0);
  xmm1 = _mm_mul_ps(signs13.vector, xmm1);
  xmm2 = _mm_mul_ps(signs13.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[0][1][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss(&u[0][2][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[2][0][1] );
  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */
  _mm_store_ps(&upper_sum[0][0][0],xmm3);
  _mm_store_ps(&upper_sum[1][0][0],xmm4);
  _mm_store_ps(&upper_sum[2][0][0],xmm5);

  /* Recons */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);
  
  xmm3 = _mm_mul_ps(signs13.vector, xmm3);
  xmm4 = _mm_mul_ps(signs13.vector, xmm4);
  xmm5 = _mm_mul_ps(signs13.vector, xmm5);
  
  /* Store up */
  _mm_store_ps(&lower_sum[0][0][0],xmm3);
  _mm_store_ps(&lower_sum[1][0][0],xmm4);
  _mm_store_ps(&lower_sum[2][0][0],xmm5);
  
}

void mvv_recons_gamma1_plus_add(const halfspinor_array src, 
				const u_mat_array u,
				halfspinor_array upper_sum, 
				halfspinor_array lower_sum)
{

  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  union { 
  float a[4];
    __m128 vector;
  } signs13 ALIGN = {{-1,1,-1,1}};

  union { 
  float a[4];
    __m128 vector;
  } signs12 ALIGN = {{-1,-1,1,1}};

  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][2][0]);
  xmm7 = _mm_load_ss(&u[2][0][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_mul_ps(signs13.vector, xmm0);
  xmm1 = _mm_mul_ps(signs13.vector, xmm1);
  xmm2 = _mm_mul_ps(signs13.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[0][1][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss(&u[0][2][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[2][0][1] );
  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);

  _mm_store_ps( &upper_sum[0][0][0],xmm0 );
  _mm_store_ps( &upper_sum[1][0][0],xmm1 );
  _mm_store_ps( &upper_sum[2][0][0],xmm2 );

  /* Load lower sum project and accumulate */
  xmm0 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm1 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm2 = _mm_load_ps( &lower_sum[2][0][0] );

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);
  
  xmm3 = _mm_mul_ps(signs12.vector, xmm3);
  xmm4 = _mm_mul_ps(signs12.vector, xmm4);
  xmm5 = _mm_mul_ps(signs12.vector, xmm5);

  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  _mm_store_ps( &lower_sum[0][0][0],xmm0 );
  _mm_store_ps( &lower_sum[1][0][0],xmm1 );
  _mm_store_ps( &lower_sum[2][0][0],xmm2 );
  

}

void mvv_recons_gamma2_plus_add(const halfspinor_array src, 
				const u_mat_array u,
				halfspinor_array upper_sum, 
				halfspinor_array lower_sum)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  union { 
  float a[4];
    __m128 vector;
  } signs13 ALIGN = {{-1,1,-1,1}};


  union { 
  float a[4];
    __m128 vector;
  } signs14 ALIGN = {{-1,1,1,-1}};


  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][2][0]);
  xmm7 = _mm_load_ss(&u[2][0][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_mul_ps(signs13.vector, xmm0);
  xmm1 = _mm_mul_ps(signs13.vector, xmm1);
  xmm2 = _mm_mul_ps(signs13.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[0][1][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss(&u[0][2][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[2][0][1] );
  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);

  _mm_store_ps( &upper_sum[0][0][0],xmm0 );
  _mm_store_ps( &upper_sum[1][0][0],xmm1 );
  _mm_store_ps( &upper_sum[2][0][0],xmm2 );

  /* Load lower sum project and accumulate */
  xmm0 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm1 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm2 = _mm_load_ps( &lower_sum[2][0][0] );

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);
  
  xmm3 = _mm_mul_ps(signs14.vector, xmm3);
  xmm4 = _mm_mul_ps(signs14.vector, xmm4);
  xmm5 = _mm_mul_ps(signs14.vector, xmm5);

  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  _mm_store_ps( &lower_sum[0][0][0],xmm0 );
  _mm_store_ps( &lower_sum[1][0][0],xmm1 );
  _mm_store_ps( &lower_sum[2][0][0],xmm2 );
  


}

void mvv_recons_gamma3_plus_add_store(const halfspinor_array src, 
			    const u_mat_array u,
			    const halfspinor_array upper_sum, 
			    const halfspinor_array lower_sum,
			    spinor_array dst)
{

  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  union { 
  float a[4];
    __m128 vector;
  } signs13 ALIGN = {{-1,1,-1,1}};


  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][2][0]);
  xmm7 = _mm_load_ss(&u[2][0][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_mul_ps(signs13.vector, xmm0);
  xmm1 = _mm_mul_ps(signs13.vector, xmm1);
  xmm2 = _mm_mul_ps(signs13.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[0][1][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss(&u[0][2][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[2][0][1] );
  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);

  /* Scatter out into the spinor */
  _mm_storel_pi((__m64 *)&dst[0][0][0],xmm0);
  _mm_storel_pi((__m64 *)&dst[0][1][0],xmm1);
  _mm_storel_pi((__m64 *)&dst[0][2][0],xmm2);

  _mm_storeh_pi((__m64 *)&dst[1][0][0],xmm0);
  _mm_storeh_pi((__m64 *)&dst[1][1][0],xmm1);
  _mm_storeh_pi((__m64 *)&dst[1][2][0],xmm2);


  /* Load lower sum and accumulate */
  xmm0 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm1 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm2 = _mm_load_ps( &lower_sum[2][0][0] );

  /* Recons -- sse_vector sub */
  xmm0 = _mm_sub_ps( xmm0, xmm3 );
  xmm1 = _mm_sub_ps( xmm1, xmm4 );
  xmm2 = _mm_sub_ps( xmm2, xmm5 );

  _mm_storel_pi((__m64 *)&dst[2][0][0],xmm0);
  _mm_storel_pi((__m64 *)&dst[2][1][0],xmm1);
  _mm_storel_pi((__m64 *)&dst[2][2][0],xmm2);

  _mm_storeh_pi((__m64 *)&dst[3][0][0],xmm0);
  _mm_storeh_pi((__m64 *)&dst[3][1][0],xmm1);
  _mm_storeh_pi((__m64 *)&dst[3][2][0],xmm2);

}



void mvv_recons_gamma0_minus(const halfspinor_array src, 
			    const u_mat_array u,
			    halfspinor_array upper_sum, halfspinor_array lower_sum)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  union { 
  float a[4];
    __m128 vector;
  } signs13 ALIGN = {{-1,1,-1,1}};

  union { 
  float a[4];
    __m128 vector;
  } signs24 ALIGN = {{1,-1,1,-1}};

  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][2][0]);
  xmm7 = _mm_load_ss(&u[2][0][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_mul_ps(signs13.vector, xmm0);
  xmm1 = _mm_mul_ps(signs13.vector, xmm1);
  xmm2 = _mm_mul_ps(signs13.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[0][1][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss(&u[0][2][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[2][0][1] );
  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */
  _mm_store_ps(&upper_sum[0][0][0],xmm3);
  _mm_store_ps(&upper_sum[1][0][0],xmm4);
  _mm_store_ps(&upper_sum[2][0][0],xmm5);

  /* Recons */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);
  
  xmm3 = _mm_mul_ps(signs24.vector, xmm3);
  xmm4 = _mm_mul_ps(signs24.vector, xmm4);
  xmm5 = _mm_mul_ps(signs24.vector, xmm5);
  
  /* Store up */
  _mm_store_ps(&lower_sum[0][0][0],xmm3);
  _mm_store_ps(&lower_sum[1][0][0],xmm4);
  _mm_store_ps(&lower_sum[2][0][0],xmm5);
  
}

void mvv_recons_gamma1_minus_add(const halfspinor_array src, 
				const u_mat_array u,
				halfspinor_array upper_sum, 
				halfspinor_array lower_sum)
{

  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  union { 
  float a[4];
    __m128 vector;
  } signs13 ALIGN = {{-1,1,-1,1}};

  union { 
  float a[4];
    __m128 vector;
  } signs34 ALIGN = {{1,1,-1,-1}};

  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][2][0]);
  xmm7 = _mm_load_ss(&u[2][0][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_mul_ps(signs13.vector, xmm0);
  xmm1 = _mm_mul_ps(signs13.vector, xmm1);
  xmm2 = _mm_mul_ps(signs13.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[0][1][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss(&u[0][2][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[2][0][1] );
  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);

  _mm_store_ps( &upper_sum[0][0][0],xmm0 );
  _mm_store_ps( &upper_sum[1][0][0],xmm1 );
  _mm_store_ps( &upper_sum[2][0][0],xmm2 );

  /* Load lower sum project and accumulate */
  xmm0 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm1 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm2 = _mm_load_ps( &lower_sum[2][0][0] );

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);
  
  xmm3 = _mm_mul_ps(signs34.vector, xmm3);
  xmm4 = _mm_mul_ps(signs34.vector, xmm4);
  xmm5 = _mm_mul_ps(signs34.vector, xmm5);

  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  _mm_store_ps( &lower_sum[0][0][0],xmm0 );
  _mm_store_ps( &lower_sum[1][0][0],xmm1 );
  _mm_store_ps( &lower_sum[2][0][0],xmm2 );
  

}

void mvv_recons_gamma2_minus_add(const halfspinor_array src, 
				const u_mat_array u,
				halfspinor_array upper_sum, 
				halfspinor_array lower_sum)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  union { 
  float a[4];
    __m128 vector;
  } signs13 ALIGN = {{-1,1,-1,1}};


  union { 
  float a[4];
    __m128 vector;
  } signs23 ALIGN = {{1,-1,-1,1}};


  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][2][0]);
  xmm7 = _mm_load_ss(&u[2][0][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_mul_ps(signs13.vector, xmm0);
  xmm1 = _mm_mul_ps(signs13.vector, xmm1);
  xmm2 = _mm_mul_ps(signs13.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[0][1][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss(&u[0][2][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[2][0][1] );
  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);

  _mm_store_ps( &upper_sum[0][0][0],xmm0 );
  _mm_store_ps( &upper_sum[1][0][0],xmm1 );
  _mm_store_ps( &upper_sum[2][0][0],xmm2 );

  /* Load lower sum project and accumulate */
  xmm0 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm1 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm2 = _mm_load_ps( &lower_sum[2][0][0] );

  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);
  
  xmm3 = _mm_mul_ps(signs23.vector, xmm3);
  xmm4 = _mm_mul_ps(signs23.vector, xmm4);
  xmm5 = _mm_mul_ps(signs23.vector, xmm5);

  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  _mm_store_ps( &lower_sum[0][0][0],xmm0 );
  _mm_store_ps( &lower_sum[1][0][0],xmm1 );
  _mm_store_ps( &lower_sum[2][0][0],xmm2 );
  


}

void mvv_recons_gamma3_minus_add_store(const halfspinor_array src, 
			    const u_mat_array u,
			    const halfspinor_array upper_sum, 
			    const halfspinor_array lower_sum,
			    spinor_array dst)
{

  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2; 
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;

  union { 
  float a[4];
    __m128 vector;
  } signs13 ALIGN = {{-1,1,-1,1}};


  /* Load Halfvector xmm0-xmm2 */
  xmm0 = _mm_load_ps( &src[0][0][0] );
  xmm1 = _mm_load_ps( &src[1][0][0] );
  xmm2 = _mm_load_ps( &src[2][0][0] );

  /* SU3 * 3 vector */

  xmm3 = _mm_load_ss(&u[0][0][0]);
  xmm6 = _mm_load_ss(&u[1][0][0]);
  xmm4 = _mm_load_ss(&u[0][1][0]);
  xmm7 = _mm_load_ss(&u[2][1][0]);
  xmm5 = _mm_load_ss(&u[0][2][0]);
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x0);
  xmm3 = _mm_mul_ps(xmm0,xmm3);
  xmm7 = _mm_shuffle_ps(xmm7,xmm7,0x0);
  xmm6 = _mm_mul_ps(xmm1,xmm6);
  xmm5 = _mm_shuffle_ps(xmm5,xmm5,0x0);
  xmm4 = _mm_mul_ps(xmm0, xmm4);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_mul_ps(xmm0, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][2][0]);
  xmm7 = _mm_load_ss(&u[2][0][0]);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm3 = _mm_add_ps(xmm7, xmm3);
  xmm6 = _mm_load_ss(&u[1][1][0]);
  xmm7 = _mm_load_ss(&u[2][2][0]);
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm4 = _mm_add_ps(xmm6, xmm4);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm6 = _mm_load_ss( &u[0][0][1] );
  xmm7 = _mm_load_ss( &u[1][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0xb1);
  xmm1 = _mm_shuffle_ps(xmm1, xmm1, 0xb1);
  xmm2 = _mm_shuffle_ps(xmm2, xmm2, 0xb1);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0 );
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0 );
  xmm0 = _mm_mul_ps(signs13.vector, xmm0);
  xmm1 = _mm_mul_ps(signs13.vector, xmm1);
  xmm2 = _mm_mul_ps(signs13.vector, xmm2);
  xmm6 = _mm_mul_ps(xmm0,xmm6);
  xmm7 = _mm_mul_ps(xmm1,xmm7);
  xmm3 = _mm_add_ps(xmm6,xmm3);
  xmm4 = _mm_add_ps(xmm7,xmm4);
  xmm6 = _mm_load_ss( &u[2][2][1] );
  xmm7 = _mm_load_ss( &u[0][1][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm2, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);
  xmm6 = _mm_load_ss(&u[1][0][1] );
  xmm7 = _mm_load_ss(&u[0][2][1] );
  xmm6 = _mm_shuffle_ps( xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps( xmm7, xmm7, 0x0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm0, xmm7);
  xmm3 = _mm_add_ps(xmm6, xmm3);
  xmm5 = _mm_add_ps(xmm7, xmm5);
  xmm0 = _mm_load_ss( &u[2][0][1] );
  xmm6 = _mm_load_ss( &u[1][2][1] );
  xmm7 = _mm_load_ss( &u[2][1][1] );
  xmm0 = _mm_shuffle_ps(xmm0, xmm0, 0x0);
  xmm6 = _mm_shuffle_ps(xmm6, xmm6, 0x0);
  xmm7 = _mm_shuffle_ps(xmm7, xmm7, 0x0);
  xmm0 = _mm_mul_ps(xmm2, xmm0);
  xmm6 = _mm_mul_ps(xmm1, xmm6);
  xmm7 = _mm_mul_ps(xmm2, xmm7);
  xmm3 = _mm_add_ps(xmm0, xmm3);
  xmm5 = _mm_add_ps(xmm6, xmm5);
  xmm4 = _mm_add_ps(xmm7, xmm4);

  /* Result in      xmm3,4,5 */
  /* END MVV */

  /* Load upper sum and accumulate */
  xmm0 = _mm_load_ps( &upper_sum[0][0][0] );
  xmm1 = _mm_load_ps( &upper_sum[1][0][0] );
  xmm2 = _mm_load_ps( &upper_sum[2][0][0] );

  xmm0 = _mm_add_ps(xmm3,xmm0);
  xmm1 = _mm_add_ps(xmm4,xmm1);
  xmm2 = _mm_add_ps(xmm5,xmm2);

  /* Scatter out into the spinor */
  _mm_storel_pi((__m64 *)&dst[0][0][0],xmm0);
  _mm_storel_pi((__m64 *)&dst[0][1][0],xmm1);
  _mm_storel_pi((__m64 *)&dst[0][2][0],xmm2);

  _mm_storeh_pi((__m64 *)&dst[1][0][0],xmm0);
  _mm_storeh_pi((__m64 *)&dst[1][1][0],xmm1);
  _mm_storeh_pi((__m64 *)&dst[1][2][0],xmm2);


  /* Load lower sum and accumulate */
  xmm0 = _mm_load_ps( &lower_sum[0][0][0] );
  xmm1 = _mm_load_ps( &lower_sum[1][0][0] );
  xmm2 = _mm_load_ps( &lower_sum[2][0][0] );

  /* Recons -- sse_vector sub */
  xmm0 = _mm_add_ps( xmm0, xmm3 );
  xmm1 = _mm_add_ps( xmm1, xmm4 );
  xmm2 = _mm_add_ps( xmm2, xmm5 );

  _mm_storel_pi((__m64 *)&dst[2][0][0],xmm0);
  _mm_storel_pi((__m64 *)&dst[2][1][0],xmm1);
  _mm_storel_pi((__m64 *)&dst[2][2][0],xmm2);

  _mm_storeh_pi((__m64 *)&dst[3][0][0],xmm0);
  _mm_storeh_pi((__m64 *)&dst[3][1][0],xmm1);
  _mm_storeh_pi((__m64 *)&dst[3][2][0],xmm2);

}



#ifdef __cplusplus
};
#endif