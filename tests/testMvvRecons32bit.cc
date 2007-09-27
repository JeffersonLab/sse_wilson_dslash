#include "testMvvRecons.h"
#include "sse32.h"
#include "types32.h"
#include "mvv_recons.h"
#include "unittest.h"

using namespace QDP;
using namespace Assertions;
using namespace std;

#define _sse_24_gamma0_minus_set() _sse_vector_xch_i_mul_up()
#define _sse_24_gamma0_plus_set()  _sse_vector_xch_i_mul_neg_up()
#define _sse_24_gamma0_minus_add() _sse_vector_xch_i_add()
#define _sse_24_gamma0_plus_add()  _sse_vector_xch_i_sub()

#define _sse_24_gamma1_minus()     _sse_vector_xch();\
				   _sse_vector_subadd()
#define _sse_24_gamma1_plus()      _sse_vector_xch();\
		 		   _sse_vector_addsub()

#define _sse_24_gamma2_minus()     _sse_vector_i_addsub()
#define _sse_24_gamma2_plus()      _sse_vector_i_subadd()

#define _sse_24_gamma3_minus()     _sse_vector_sub()
#define _sse_24_gamma3_plus()      _sse_vector_add()
#define _sse_24_gamma3_minus_rows12() _sse_vector_add()
#define _sse_24_gamma3_plus_rows12() _sse_vector_add()

void testMvvRecons0Plus::run(void) 
{
  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_store_up(r12_1);
  _sse_24_gamma0_minus_set(); 
  _sse_vector_store_up(r34_1);

  mvv_recons_gamma0_plus(hspinor, matrix, r12_2, r34_2);
  
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = r12_1[col][spin2][reim] - r12_2[col][spin2][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << endl;
#endif
	assertion( diff < 1.0e-6 );

	float diff2 = r34_1[col][spin2][reim] - r34_2[col][spin2][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff lower = " << diff2 
		    << endl;
#endif 
	assertion( diff2 < 1.0e-6 );

      }
    }
  }

}

void testMvvRecons1PlusAdd::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_vector_store(r12_1);
  _sse_vector_load(r34_1);
  _sse_24_gamma1_minus();
  _sse_vector_store(r34_1);




  mvv_recons_gamma1_plus_add(hspinor, matrix, r12_2, r34_2);
  

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = r12_1[col][spin2][reim] - r12_2[col][spin2][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << endl;
#endif

	assertion( diff < 1.0e-6 );

	float diff2 = r34_1[col][spin2][reim] - r34_2[col][spin2][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff lower = " << diff2 
		    << endl;
#endif
	assertion( diff2 < 1.0e-6 );

      }
    }
  }

}

void testMvvRecons2PlusAdd::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_vector_store(r12_1);
  _sse_vector_load(r34_1);
  _sse_24_gamma2_minus();
  _sse_vector_store(r34_1);




  mvv_recons_gamma2_plus_add(hspinor, matrix, r12_2, r34_2);
  

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = r12_1[col][spin2][reim] - r12_2[col][spin2][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << endl;
#endif

	assertion( diff < 1.0e-6 );

	float diff2 = r34_1[col][spin2][reim] - r34_2[col][spin2][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff lower = " << diff2 
		    << endl;
#endif
	assertion( diff2 < 1.0e-6 );

      }
    }
  }

}

void testMvvRecons3PlusAddStore::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  spinor_array spinor1;
  spinor_array spinor2;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_pair_store(spinor1[0], spinor1[1]);
  _sse_vector_load(r34_1);
  _sse_24_gamma3_minus();
  _sse_pair_store(spinor1[2], spinor1[3]);



  mvv_recons_gamma3_plus_add_store(hspinor, matrix, r12_2, r34_2, spinor2);
  

  for(int col=0; col < 3; col++) { 
    for(int spin=0; spin < 4; spin++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = spinor1[spin][col][reim] - spinor2[spin][col][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << endl;
#endif
	assertion( diff < 1.0e-6 );

      }
    }
  }

}

void testMvvRecons0Minus::run(void) 
{
  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_store_up(r12_1);
  _sse_24_gamma0_plus_set(); 
  _sse_vector_store_up(r34_1);

  mvv_recons_gamma0_minus(hspinor, matrix, r12_2, r34_2);
  
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = r12_1[col][spin2][reim] - r12_2[col][spin2][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << endl;
#endif
	assertion( diff < 1.0e-6 );

	float diff2 = r34_1[col][spin2][reim] - r34_2[col][spin2][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff lower = " << diff2 
		    << endl;
#endif 
	assertion( diff2 < 1.0e-6 );

      }
    }
  }

}

void testMvvRecons1MinusAdd::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_vector_store(r12_1);
  _sse_vector_load(r34_1);
  _sse_24_gamma1_plus();
  _sse_vector_store(r34_1);




  mvv_recons_gamma1_minus_add(hspinor, matrix, r12_2, r34_2);
  

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = r12_1[col][spin2][reim] - r12_2[col][spin2][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << endl;
#endif

	assertion( diff < 1.0e-6 );

	float diff2 = r34_1[col][spin2][reim] - r34_2[col][spin2][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff lower = " << diff2 
		    << endl;
#endif
	assertion( diff2 < 1.0e-6 );

      }
    }
  }

}

void testMvvRecons2MinusAdd::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_vector_store(r12_1);
  _sse_vector_load(r34_1);
  _sse_24_gamma2_plus();
  _sse_vector_store(r34_1);




  mvv_recons_gamma2_minus_add(hspinor, matrix, r12_2, r34_2);
  

  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = r12_1[col][spin2][reim] - r12_2[col][spin2][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << endl;
#endif

	assertion( diff < 1.0e-6 );

	float diff2 = r34_1[col][spin2][reim] - r34_2[col][spin2][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin2
		    << " re=" << reim 
		    << " diff lower = " << diff2 
		    << endl;
#endif
	assertion( diff2 < 1.0e-6 );

      }
    }
  }

}

void testMvvRecons3MinusAddStore::run(void) 
{

  halfspinor_array hspinor;
  u_mat_array matrix;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN;
  halfspinor_array r12_2 ALIGN,r34_2 ALIGN;

  spinor_array spinor1;
  spinor_array spinor2;


  /* Random numbers in halfspinors */
  for(int col=0; col < 3; col++) { 
    for(int spin2=0; spin2 < 2; spin2++) { 
      for(int reim=0; reim < 2; reim++) { 
	hspinor[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r12_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
	r34_1[col][spin2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;

	/* Set up 'existing partial sum' */
	r12_2[col][spin2][reim] = r12_1[col][spin2][reim];
	r34_2[col][spin2][reim] = r34_1[col][spin2][reim];
	

      }
    }
  }
  

  // Random matrix 
  for(int col1=0; col1 < 3; col1++) { 
    for(int col2=0; col2 < 3; col2++) { 
      for(int reim=0; reim < 2; reim++) { 
	matrix[col1][col2][reim] = (float)(rand() - RAND_MAX/2)/(float)(RAND_MAX/2)*2.0;
      }
    }
  }

  _sse_vector_load(hspinor);
  _sse_su3_multiply((matrix));
  _sse_vector_load(r12_1);
  _sse_vector_add();
  _sse_pair_store(spinor1[0], spinor1[1]);
  _sse_vector_load(r34_1);
  _sse_24_gamma3_plus();
  _sse_pair_store(spinor1[2], spinor1[3]);



  mvv_recons_gamma3_minus_add_store(hspinor, matrix, r12_2, r34_2, spinor2);
  

  for(int col=0; col < 3; col++) { 
    for(int spin=0; spin < 4; spin++) { 
      for(int reim=0; reim < 2; reim++) { 
	float diff = spinor1[spin][col][reim] - spinor2[spin][col][reim];

#if 0
	QDPIO::cout << "  col=" << col 
		    << " sp=" << spin
		    << " re=" << reim 
		    << " diff upper = " << diff 
		    << endl;
#endif
	assertion( diff < 1.0e-6 );

      }
    }
  }

}
