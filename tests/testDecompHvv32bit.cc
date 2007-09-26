#include "testDecompHvv.h"
#include "sse32.h"
#include "types32.h"
#include "decomp_hvv.h"
#include "unittest.h"

using namespace QDP;
using namespace Assertions;
using namespace std;

#define _sse_42_gamma0_minus()     _sse_vector_xch_i_sub()
#define _sse_42_gamma0_plus()      _sse_vector_xch_i_add()
#define _sse_24_gamma0_minus_set() _sse_vector_xch_i_mul_up()
#define _sse_24_gamma0_plus_set()  _sse_vector_xch_i_mul_neg_up()
#define _sse_24_gamma0_minus_add() _sse_vector_xch_i_add()
#define _sse_24_gamma0_plus_add()  _sse_vector_xch_i_sub()


/* gamma 1 */
#define _sse_42_gamma1_minus()     _sse_vector_xch();\
				   _sse_vector_addsub()
#define _sse_42_gamma1_plus()      _sse_vector_xch();\
				   _sse_vector_subadd()
#define _sse_24_gamma1_minus()     _sse_vector_xch();\
				   _sse_vector_subadd()
#define _sse_24_gamma1_plus()      _sse_vector_xch();\
		 		   _sse_vector_addsub()

/* gamma 2 */
#define _sse_42_gamma2_minus()     _sse_vector_i_subadd()
#define _sse_42_gamma2_plus()      _sse_vector_i_addsub()
#define _sse_24_gamma2_minus()     _sse_vector_i_addsub()
#define _sse_24_gamma2_plus()      _sse_vector_i_subadd()

/* gamma 3 */
#define _sse_42_gamma3_minus()     _sse_vector_sub()
#define _sse_42_gamma3_plus()      _sse_vector_add()
#define _sse_24_gamma3_minus()     _sse_vector_sub()
#define _sse_24_gamma3_plus()      _sse_vector_add()
#define _sse_24_gamma3_minus_rows12() _sse_vector_add()
#define _sse_24_gamma3_plus_rows12() _sse_vector_add()


void testDecompHvv0Minus::run(void) 
{
  fail("Deliberate failure");
}


void testDecompHvv0Plus::run(void) 
{
  spinor_array spinor;
  halfspinor_array hspinor1, hspinor2;
  u_mat_array matrix;

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor1[col][spin2][reim] = 0;
	hspinor2[col][spin2][reim] = 0;
      }
    }
  }

  for(int col=0; col < 3; col++) { 
    for(int spin4=0; spin4 < 4; spin4++) {
      for(int reim=0; reim < 2; reim++) { 

	spinor[col][spin4][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	
      }
    }
  } // Color

  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  spinor_array* sm1 = &spinor;
  u_mat_array*  u = &matrix;
  halfspinor_array *s3 = &hspinor1;

  // Copy from code 
  _sse_pair_load((*sm1)[0],(*sm1)[1]);
  _sse_pair_load_up((*sm1)[2],(*sm1)[3]);
  _sse_42_gamma0_plus();
  _sse_su3_inverse_multiply((*u));
  _sse_vector_store_up(*s3);

  // My version
  decomp_hvv_gamma0_plus(spinor, matrix, hspinor2);

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = hspinor1[col][spin2][reim] - hspinor2[col][spin2][reim];
	assertion( diff < 1.0e-6 );
      }
    }
  }
}

void testDecompHvv1Minus::run(void) 
{
  fail("Deliberate failure");
}

void testDecompHvv1Plus::run(void) 
{
  fail("Deliberate failure");
}

void testDecompHvv2Minus::run(void) 
{
  fail("Deliberate failure");
}

void testDecompHvv2Plus::run(void) 
{
  fail("Deliberate failure");
}


void testDecompHvv3Minus::run(void) 
{
  fail("Deliberate failure");
}

void testDecompHvv3Plus::run(void) 
{
  fail("Deliberate failure");
}
