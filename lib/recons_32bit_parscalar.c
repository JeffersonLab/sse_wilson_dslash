#include "recons.h"
#include "xmmintrin.h"

#ifdef __cplusplus
extern "C" {
#endif

 typedef union { 
    unsigned int a[4];
    __m128 vector;
  } SSESign;




void recons_4dir_plus( halfspinor_array hs0,
		       halfspinor_array hs1,
		       halfspinor_array hs2,
		       halfspinor_array hs3,
		      spinor_array spinor)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;
  __m128 xmm8;

  __m128 xmm9;
  __m128 xmm10;


  SSESign signs24 __attribute__((unused)) ALIGN= {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN= {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN= {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};


  /* The partial sum is stored in right order */
  xmm0 = _mm_load_ps(&spinor[0][0][0]);
  xmm1 = _mm_load_ps(&spinor[0][2][0]);
  xmm2 = _mm_load_ps(&spinor[1][1][0]);

  xmm6 = _mm_load_ps(&spinor[2][0][0]);
  xmm7 = _mm_load_ps(&spinor[2][2][0]);
  xmm8 = _mm_load_ps(&spinor[3][1][0]);

  /* Top components. Don't need to reruct just
     accumulate */

  /* half spinor from dir 0 */
  /* index order is [color][spin][reim] */
  xmm3 = _mm_load_ps(&hs0[0][0][0]);
  xmm4 = _mm_load_ps(&hs0[1][0][0]); 
  xmm5 = _mm_load_ps(&hs0[2][0][0]);
 
  /* Add it */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

  xmm3 = _mm_xor_ps(signs24.vector, xmm3);
  xmm4 = _mm_xor_ps(signs24.vector, xmm4);
  xmm5 = _mm_xor_ps(signs24.vector, xmm5);
  
  xmm6 = _mm_add_ps( xmm3, xmm6 );
  xmm7 = _mm_add_ps( xmm4, xmm7 );
  xmm8 = _mm_add_ps( xmm5, xmm8 );

  /* Half  spinor from dir 1 */
  xmm3 = _mm_load_ps(&hs1[0][0][0]);
  xmm4 = _mm_load_ps(&hs1[1][0][0]);
  xmm5 = _mm_load_ps(&hs1[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

  xmm3 = _mm_xor_ps(signs34.vector, xmm3);
  xmm4 = _mm_xor_ps(signs34.vector, xmm4);
  xmm5 = _mm_xor_ps(signs34.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Half spinor from dir 2 */
  xmm3 = _mm_load_ps(&hs2[0][0][0]);
  xmm4 = _mm_load_ps(&hs2[1][0][0]);
  xmm5 = _mm_load_ps(&hs2[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(signs23.vector, xmm3);
  xmm4 = _mm_xor_ps(signs23.vector, xmm4);
  xmm5 = _mm_xor_ps(signs23.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);


  /* Half Spinor from dir 3 */
  xmm3 = _mm_load_ps(&hs3[0][0][0]);
  xmm4 = _mm_load_ps(&hs3[1][0][0]);
  xmm5 = _mm_load_ps(&hs3[2][0][0]);

  /* Add it */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Accumulate */
  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Deswizzle and store */
  /* Store top half of result spinor */
  _mm_storel_pi((__m64*)&spinor[0][0][0],xmm0);
  _mm_storel_pi((__m64*)&spinor[0][1][0],xmm1);
  _mm_storel_pi((__m64*)&spinor[0][2][0],xmm2);
  _mm_storeh_pi((__m64*)&spinor[1][0][0],xmm0);
  _mm_storeh_pi((__m64*)&spinor[1][1][0],xmm1);
  _mm_storeh_pi((__m64*)&spinor[1][2][0],xmm2);

  _mm_storel_pi((__m64*)&spinor[2][0][0],xmm6);
  _mm_storel_pi((__m64*)&spinor[2][1][0],xmm7);
  _mm_storel_pi((__m64*)&spinor[2][2][0],xmm8);
  _mm_storeh_pi((__m64*)&spinor[3][0][0],xmm6);
  _mm_storeh_pi((__m64*)&spinor[3][1][0],xmm7);
  _mm_storeh_pi((__m64*)&spinor[3][2][0],xmm8);

  
}


void recons_3dir_plus( halfspinor_array hs0,
		       halfspinor_array hs1,
		       halfspinor_array hs2,
		       spinor_array spinor)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;
  __m128 xmm8;

  __m128 xmm9;
  __m128 xmm10;

  SSESign signs24 __attribute__((unused)) ALIGN= {{ 0x00000000, 0x80000000, 0x00000000, 0x80000000 }};
  SSESign signs34 __attribute__((unused)) ALIGN= {{ 0x00000000, 0x00000000, 0x80000000, 0x80000000 }};
  SSESign signs23 __attribute__((unused)) ALIGN= {{ 0x00000000, 0x80000000, 0x80000000, 0x00000000 }};

  xmm0 = _mm_load_ps(&spinor[0][0][0]);
  xmm1 = _mm_load_ps(&spinor[0][2][0]);
  xmm2 = _mm_load_ps(&spinor[1][1][0]);
  xmm6 = _mm_load_ps(&spinor[2][0][0]);
  xmm7 = _mm_load_ps(&spinor[2][2][0]);
  xmm8 = _mm_load_ps(&spinor[3][1][0]);

  /* Top components. Don't need to reruct just
     accumulate */

  /* half spinor from dir 0 */
  /* index order is [color][spin][reim] */
  xmm3 = _mm_load_ps(&hs0[0][0][0]);
  xmm4 = _mm_load_ps(&hs0[1][0][0]); 
  xmm5 = _mm_load_ps(&hs0[2][0][0]);
 
  /* Add it */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

  xmm3 = _mm_xor_ps(signs24.vector, xmm3);
  xmm4 = _mm_xor_ps(signs24.vector, xmm4);
  xmm5 = _mm_xor_ps(signs24.vector, xmm5);
  
  xmm6 = _mm_add_ps( xmm3, xmm6 );
  xmm7 = _mm_add_ps( xmm4, xmm7 );
  xmm8 = _mm_add_ps( xmm5, xmm8 );


  /* Half  spinor from dir 1 */
  xmm3 = _mm_load_ps(&hs1[0][0][0]);
  xmm4 = _mm_load_ps(&hs1[1][0][0]);
  xmm5 = _mm_load_ps(&hs1[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

  xmm3 = _mm_xor_ps(signs34.vector, xmm3);
  xmm4 = _mm_xor_ps(signs34.vector, xmm4);
  xmm5 = _mm_xor_ps(signs34.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Half spinor from dir 2 */
  xmm3 = _mm_load_ps(&hs2[0][0][0]);
  xmm4 = _mm_load_ps(&hs2[1][0][0]);
  xmm5 = _mm_load_ps(&hs2[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(signs23.vector, xmm3);
  xmm4 = _mm_xor_ps(signs23.vector, xmm4);
  xmm5 = _mm_xor_ps(signs23.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Deswizzle and store */
  /* Store top half of result spinor */
  _mm_storel_pi((__m64*)&spinor[0][0][0],xmm0);
  _mm_storel_pi((__m64*)&spinor[0][1][0],xmm1);
  _mm_storel_pi((__m64*)&spinor[0][2][0],xmm2);
  _mm_storeh_pi((__m64*)&spinor[1][0][0],xmm0);
  _mm_storeh_pi((__m64*)&spinor[1][1][0],xmm1);
  _mm_storeh_pi((__m64*)&spinor[1][2][0],xmm2);


  _mm_storel_pi((__m64*)&spinor[2][0][0],xmm6);
  _mm_storel_pi((__m64*)&spinor[2][1][0],xmm7);
  _mm_storel_pi((__m64*)&spinor[2][2][0],xmm8);
  _mm_storeh_pi((__m64*)&spinor[3][0][0],xmm6);
  _mm_storeh_pi((__m64*)&spinor[3][1][0],xmm7);
  _mm_storeh_pi((__m64*)&spinor[3][2][0],xmm8);

  
}


void recons_4dir_minus( halfspinor_array hs0,
		        halfspinor_array hs1,
		        halfspinor_array hs2,
		        halfspinor_array hs3,
		       spinor_array spinor)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;
  __m128 xmm8;

  __m128 xmm9;
  __m128 xmm10;

  SSESign signs13 __attribute__((unused)) ALIGN= {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN= {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN= {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};  


  xmm0 = _mm_load_ps(&spinor[0][0][0]);
  xmm1 = _mm_load_ps(&spinor[0][2][0]);
  xmm2 = _mm_load_ps(&spinor[1][1][0]);
  xmm6 = _mm_load_ps(&spinor[2][0][0]);
  xmm7 = _mm_load_ps(&spinor[2][2][0]);
  xmm8 = _mm_load_ps(&spinor[3][1][0]);


  /* Top components. Don't need to reconstruct just
     accumulate */

  /* half spinor from dir 0 */
  /* index order is [color][spin][reim] */
  xmm3 = _mm_load_ps(&hs0[0][0][0]);
  xmm4 = _mm_load_ps(&hs0[1][0][0]); 
  xmm5 = _mm_load_ps(&hs0[2][0][0]);
 
  /* Add it */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

  xmm3 = _mm_xor_ps(signs13.vector, xmm3);
  xmm4 = _mm_xor_ps(signs13.vector, xmm4);
  xmm5 = _mm_xor_ps(signs13.vector, xmm5);
  
  xmm6 = _mm_add_ps( xmm3, xmm6 );
  xmm7 = _mm_add_ps( xmm4, xmm7 );
  xmm8 = _mm_add_ps( xmm5, xmm8 );


  /* Half  spinor from dir 1 */
  xmm3 = _mm_load_ps(&hs1[0][0][0]);
  xmm4 = _mm_load_ps(&hs1[1][0][0]);
  xmm5 = _mm_load_ps(&hs1[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

  xmm3 = _mm_xor_ps(signs12.vector, xmm3);
  xmm4 = _mm_xor_ps(signs12.vector, xmm4);
  xmm5 = _mm_xor_ps(signs12.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);


  /* Half spinor from dir 2 */
  xmm3 = _mm_load_ps(&hs2[0][0][0]);
  xmm4 = _mm_load_ps(&hs2[1][0][0]);
  xmm5 = _mm_load_ps(&hs2[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(signs14.vector, xmm3);
  xmm4 = _mm_xor_ps(signs14.vector, xmm4);
  xmm5 = _mm_xor_ps(signs14.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Half Spinor from dir 3 */
  xmm3 = _mm_load_ps(&hs3[0][0][0]);
  xmm4 = _mm_load_ps(&hs3[1][0][0]);
  xmm5 = _mm_load_ps(&hs3[2][0][0]);

  /* Add it */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Accumulate */
  xmm6 = _mm_sub_ps(xmm6, xmm3);
  xmm7 = _mm_sub_ps(xmm7, xmm4);
  xmm8 = _mm_sub_ps(xmm8, xmm5);

  /* Deswizzle and store */
  /* Store top half of result spinor */
  _mm_storel_pi((__m64*)&spinor[0][0][0],xmm0);
  _mm_storel_pi((__m64*)&spinor[0][1][0],xmm1);
  _mm_storel_pi((__m64*)&spinor[0][2][0],xmm2);
  _mm_storeh_pi((__m64*)&spinor[1][0][0],xmm0);
  _mm_storeh_pi((__m64*)&spinor[1][1][0],xmm1);
  _mm_storeh_pi((__m64*)&spinor[1][2][0],xmm2);

  _mm_storel_pi((__m64*)&spinor[2][0][0],xmm6);
  _mm_storel_pi((__m64*)&spinor[2][1][0],xmm7);
  _mm_storel_pi((__m64*)&spinor[2][2][0],xmm8);
  _mm_storeh_pi((__m64*)&spinor[3][0][0],xmm6);
  _mm_storeh_pi((__m64*)&spinor[3][1][0],xmm7);
  _mm_storeh_pi((__m64*)&spinor[3][2][0],xmm8);


}

void recons_3dir_minus( halfspinor_array hs0,
		        halfspinor_array hs1,
		        halfspinor_array hs2,
			spinor_array spinor)
{
  __m128 xmm0;
  __m128 xmm1;
  __m128 xmm2;
  __m128 xmm3;
  __m128 xmm4;
  __m128 xmm5;
  __m128 xmm6;
  __m128 xmm7;
  __m128 xmm8;

  __m128 xmm9;
  __m128 xmm10;
  
  SSESign signs13 __attribute__((unused)) ALIGN= {{ 0x80000000, 0x00000000, 0x80000000, 0x00000000 }};
  SSESign signs12 __attribute__((unused)) ALIGN= {{ 0x80000000, 0x80000000, 0x00000000, 0x00000000 }};
  SSESign signs14 __attribute__((unused)) ALIGN= {{ 0x80000000, 0x00000000, 0x00000000, 0x80000000 }};  


  xmm0 = _mm_load_ps(&spinor[0][0][0]);
  xmm1 = _mm_load_ps(&spinor[0][2][0]);
  xmm2 = _mm_load_ps(&spinor[1][1][0]);
  xmm6 = _mm_load_ps(&spinor[2][0][0]);
  xmm7 = _mm_load_ps(&spinor[2][2][0]);
  xmm8 = _mm_load_ps(&spinor[3][1][0]);

  /* Top components. Don't need to reconstruct just
     accumulate */

  /* half spinor from dir 0 */
  /* index order is [color][spin][reim] */
  xmm3 = _mm_load_ps(&hs0[0][0][0]);
  xmm4 = _mm_load_ps(&hs0[1][0][0]); 
  xmm5 = _mm_load_ps(&hs0[2][0][0]);
 
  /* Add it */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);


  /* recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x1b);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x1b);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x1b);

  xmm3 = _mm_xor_ps(signs13.vector, xmm3);
  xmm4 = _mm_xor_ps(signs13.vector, xmm4);
  xmm5 = _mm_xor_ps(signs13.vector, xmm5);
  
  xmm6 = _mm_add_ps( xmm3, xmm6 );
  xmm7 = _mm_add_ps( xmm4, xmm7 );
  xmm8 = _mm_add_ps( xmm5, xmm8 );

  /* Half  spinor from dir 1 */
  xmm3 = _mm_load_ps(&hs1[0][0][0]);
  xmm4 = _mm_load_ps(&hs1[1][0][0]);
  xmm5 = _mm_load_ps(&hs1[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0x4e);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0x4e);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0x4e);

  xmm3 = _mm_xor_ps(signs12.vector, xmm3);
  xmm4 = _mm_xor_ps(signs12.vector, xmm4);
  xmm5 = _mm_xor_ps(signs12.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Half spinor from dir 2 */
  xmm3 = _mm_load_ps(&hs2[0][0][0]);
  xmm4 = _mm_load_ps(&hs2[1][0][0]);
  xmm5 = _mm_load_ps(&hs2[2][0][0]);

  /* Accumulate */
  xmm0 = _mm_add_ps(xmm3, xmm0);
  xmm1 = _mm_add_ps(xmm4, xmm1);
  xmm2 = _mm_add_ps(xmm5, xmm2);

  /* Recons + add */
  xmm3 = _mm_shuffle_ps(xmm3, xmm3, 0xb1);
  xmm4 = _mm_shuffle_ps(xmm4, xmm4, 0xb1);
  xmm5 = _mm_shuffle_ps(xmm5, xmm5, 0xb1);

  xmm3 = _mm_xor_ps(signs14.vector, xmm3);
  xmm4 = _mm_xor_ps(signs14.vector, xmm4);
  xmm5 = _mm_xor_ps(signs14.vector, xmm5);

  xmm6 = _mm_add_ps(xmm3, xmm6);
  xmm7 = _mm_add_ps(xmm4, xmm7);
  xmm8 = _mm_add_ps(xmm5, xmm8);

  /* Deswizzle and store */
  _mm_storel_pi((__m64*)&spinor[0][0][0],xmm0);
  _mm_storel_pi((__m64*)&spinor[0][1][0],xmm1);
  _mm_storel_pi((__m64*)&spinor[0][2][0],xmm2);
  _mm_storeh_pi((__m64*)&spinor[1][0][0],xmm0);
  _mm_storeh_pi((__m64*)&spinor[1][1][0],xmm1);
  _mm_storeh_pi((__m64*)&spinor[1][2][0],xmm2);

  _mm_storel_pi((__m64*)&spinor[2][0][0],xmm6);
  _mm_storel_pi((__m64*)&spinor[2][1][0],xmm7);
  _mm_storel_pi((__m64*)&spinor[2][2][0],xmm8);
  _mm_storeh_pi((__m64*)&spinor[3][0][0],xmm6);
  _mm_storeh_pi((__m64*)&spinor[3][1][0],xmm7);
  _mm_storeh_pi((__m64*)&spinor[3][2][0],xmm8);


}


#ifdef __cplusplus
};
#endif
