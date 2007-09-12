/*******************************************************************************
 * $Id: sse_su3dslash_32bit_parscalar.c,v 1.2 2007-09-12 21:00:50 bjoo Exp $
 * 
 * Action of the 32bit parallel Wilson-Dirac operator D_w on a given spinor field
 *
 * The externally accessible function is   sse_su3dslash_wilson
 * 
 * void sse_su3dslash_wilson(float *u, float *psi, float *res, int isign, int cb)
 *
 * NOTE:
 * u: base pointer to gauge field
 * psi: base pointer to input spinor field on FULL lattice
 * res: base pointer to output spinor field on FULL lattice
 * isign: 1-->normal, -1--> swaps 1 - gamma(mu^) for 1 + gamma(mu^)
 * cb: checkerboard (0/1) of input fields
 *
 * Oringal author: Chris McClendon <cmcclend@jlab.org>
 * Acknowledgements to: Martin Luescher <luscher@mail.desy.de>
 * Date: 9/15/2001
 *
 *******************************************************************************/

/* requires gauge fields packed by pack_gauge_field of intpar_table.c */

/* externally callable function: wnxtsu3dslash */
/* This routine applies the operator D' to Psi, putting the result in Res. */

/*	       Nd-1 */
/*	       --- */
/*	       \ */
/*   res(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu) */
/*	       /     mu                    mu */
/*	       --- */
/*	       mu=0 */

/*	             Nd-1 */
/*	             --- */
/*	             \    + */
/*                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu) */
/*	             /     mu                       mu */
/*	             --- */
/*	             mu=0 */

/* Arguments: */

/* U	     Gauge field					(Read) */
/* Psi	     Pseudofermion field				(Read) */
/* Res	     Pseudofermion field				(Write) */
/*		      + */
/* ISign     D' or D'  ( +1 | -1 ) respectively	                (Read) */
/* CB	     Checkerboard of input vector			(Read) */

/* xchi1     Pseudo-halffermion field with 2 tail buffers       (Read) */
/* xchi2     Pseudo-halffermion field with 2 tail buffers       (Read)  */
/* shift     Shift table to use..see intpar_table.c                  */
/* bound     an array with the number of boundaries per direction...see intpar_table.c (Read) */

#include <sse_config.h>

#ifdef __cplusplus
extern "C" {
#endif


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <qmp.h>

#ifdef DSLASH_USE_QMT_THREADS
#include <qmt.h>
#endif

#include <sse32.h>
extern void make_shift_tables(int *shift, int icolor_start[2], int bound[2][4][4]);

#define BASE 0x3f


#include <sse_align.h>

static int subgrid_vol = 0;
static int subgrid_vol_cb = 0;
static int initP=0;

/* iup idn shift table overlays with packed/unpacked gauge lookerupper*/
static int *newshift;
static int icolor_start[2];    /* starting site for each coloring (cb) */
static int icolor_end[2];      /* end site for each coloring (cb) */

#define BACKWARD 0
#define FORWARD 1


/* here are the macros for the accessing the shift tables, etc...
 * iup means x + mu^, idn means x - mu^, and most of the
 * mvv_shift and rec_shift values end up being identity values, except 
 * those that pull in from the "recieving" tail buffer (TAIL2)...it may 
 * be more efficient to compute
 * whether or not sites are on the boundary or not on the fly using MMX registers 
 */
#define iup(mysite,mymu) newshift[mymu+4*(mysite+subgrid_vol*(1))]
#define idn(mysite,mymu) newshift[mymu+4*(mysite+subgrid_vol*(0))]
#define mvv_shift(mysite,mymu) newshift[mymu+4*(mysite+subgrid_vol*(2))]
#define rec_shift(mysite,mymu) newshift[mymu+4*(mysite+subgrid_vol*(3))]

#define a_chia(mymu,mysite) (chia+mysite+3*subgrid_vol_cb*mymu)
#define a_chib(mymu,mysite) (chib+mysite+3*subgrid_vol_cb*mymu)



/* here's some redefines if non_temporal storing is on...non temporal stores 
 * bypass the cache on their way to memory, 
 * resulting in less cache pollution */
/* the non temporal prefetches just target one way of lvl2 cache..might be 
 * good to use, might not */

#define _sse_vector_store_NTA _sse_vector_store
#define _sse_vector_store_up_NTA _sse_vector_store

#define _prefetch_single_NTA _prefetch_single_nta
#define _prefetch_spinor_NTA _prefetch_spinor_nta

/* the default is packed...that's why the statements look non-intuitive */
/* here's how packing works
	u <-- gauge_field
   in lexical order of u(site, mu) the packed gauge fields are laid out as follows:
	u(site,0)
	u(site+1,0)
	u(site,1)
	u(site+1,1)
	u(site,2)
	u(site+1,2)
	u(site,3)
	u(site+1,3)
*/

/* for SZIN we use packed because we have unrolled the loops because we can take
 * advantage of the long cahce line on the P4...so pack_gauge_field in intpar_table.c
 * needs to be called */
  
#define _gauge_field0_0(mysite) gauge_field[mysite][0]
#define _gauge_field0_1(mysite) gauge_field[mysite][1] 
#define _gauge_field0_2(mysite) gauge_field[mysite][2]
#define _gauge_field0_3(mysite) gauge_field[mysite][3]
#define _gauge_field1_0(mysite) gauge_field[mysite+1][0]
#define _gauge_field1_1(mysite) gauge_field[mysite+1][1]
#define _gauge_field1_2(mysite) gauge_field[mysite+1][2]
#define _gauge_field1_3(mysite) gauge_field[mysite+1][3]
#define _gauge_field0_0_opp(mysite) gauge_field_opp[mysite][0]
#define _gauge_field0_1_opp(mysite) gauge_field_opp[mysite][1]
#define _gauge_field0_2_opp(mysite) gauge_field_opp[mysite][2]
#define _gauge_field0_3_opp(mysite) gauge_field_opp[mysite][3]
#define _gauge_field1_0_opp(mysite) gauge_field_opp[mysite+1][0]
#define _gauge_field1_1_opp(mysite) gauge_field_opp[mysite+1][1]
#define _gauge_field1_2_opp(mysite) gauge_field_opp[mysite+1][2]
#define _gauge_field1_3_opp(mysite) gauge_field_opp[mysite+1][3]
#define a_gauge_field0_0(mysite) &gauge_field[mysite][0]
#define a_gauge_field0_1(mysite) &gauge_field[mysite][1]
#define a_gauge_field0_2(mysite) &gauge_field[mysite][2]
#define a_gauge_field0_3(mysite) &gauge_field[mysite][3]
#define a_gauge_field1_0(mysite) &gauge_field[mysite+1][0]
#define a_gauge_field1_1(mysite) &gauge_field[mysite+1][1]
#define a_gauge_field1_2(mysite) &gauge_field[mysite+1][2]
#define a_gauge_field1_3(mysite) &gauge_field[mysite+1][3]
#define a_gauge_field0_0_opp(mysite) &gauge_field_opp[mysite][0]
#define a_gauge_field0_1_opp(mysite) &gauge_field_opp[mysite][1]
#define a_gauge_field0_2_opp(mysite) &gauge_field_opp[mysite][2]
#define a_gauge_field0_3_opp(mysite) &gauge_field_opp[mysite][3]
#define a_gauge_field1_0_opp(mysite) &gauge_field_opp[mysite+1][0]
#define a_gauge_field1_1_opp(mysite) &gauge_field_opp[mysite+1][1]
#define a_gauge_field1_2_opp(mysite) &gauge_field_opp[mysite+1][2]
#define a_gauge_field1_3_opp(mysite) &gauge_field_opp[mysite+1][3]


/* now overlays for spinors as arrays or structs */
typedef float chi_float[2];
typedef chi_float chi_two[2];
typedef float u_mat_array[3][3][2]  ALIGN;  /* color color re/im */ 
typedef float spinor_array[4][3][2] ALIGN; /* Nspin4 color re/im */
typedef chi_two chi_array[3]    ALIGN; /*..color Nspin2 re/im ::note:: color slowest varying */
typedef u_mat_array (*my_mat_array)[4] ALIGN;  


#define MY_SPINOR spinor_array
#define MY_SSE_VECTOR chi_array
#define MY_SSE_FLOAT sse_float  
#define _c1__ [0] 
#define _c2__ [1] 
#define _c3__ [2]
#define _c4__ [3]


#define MY_GAUGE u_mat_array
#define MY_GAUGE_ARRAY my_mat_array


/* end overlays */

/* macros for spin basis note: assume first two rows are linearly independent 
 * except for gamma3 */


/* macros for spin basis note: assume first two rows are linearly independent 
 * except for gamma3 */
/* it should be possible to change the spin basis by just correctly modifying 
 * these macros. However, some non-chiral spin basis may need additional
 * modifications inside the dslash routine*/


/* gamma 0 */
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


/* end spin basis section */


#define min(a,b)  ((a < b) ? a : b)

/* declare stuff...again, the vector constants must be 
 * declared in each routine to assure 16-byte alignment */


#define DECL_COMMON_ALIASES_TEMPS \
int ix1,iy1,iy2,iz1; \
MY_GAUGE* up1 ALIGN;       \
MY_GAUGE* up2 ALIGN;       \
MY_GAUGE* um1 ALIGN;       \
MY_GAUGE* um2 ALIGN;       \
MY_GAUGE* um3 ALIGN;       \
MY_SPINOR* s1 ALIGN;       \
MY_SPINOR* sp1 ALIGN;      \
MY_SPINOR* sp2 ALIGN;      \
MY_SPINOR* sm2 ALIGN;      \
MY_SPINOR* sm1 ALIGN;      \
MY_SPINOR* sn1 ALIGN;      \
sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};\
sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};\
sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};\
sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};\
sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};\
sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};\
sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f};\
MY_SSE_FLOAT fact1,fact2; \
MY_SSE_VECTOR r12_1 ALIGN,r34_1 ALIGN,r12_2 ALIGN,r34_2 ALIGN   /*so user can put semicolon in */


/* if we are allocating arrays on the fly or at compile time */ 
/* note down below things will need to be changed as well in each proc */
/* I assume multidimensional arrays are just a big single dimensional 
 * array at the site level, but inside working with colors and spins I 
 * use array indices to handle the offsetting for me :P */


typedef struct {
  MY_SPINOR* spinfun;
  MY_SSE_VECTOR *chifun;  /*note this must be allocated ..direction varying more slowly than vol */
  MY_GAUGE       (*u)[4];
  MY_GAUGE       (*u2)[4];
  int cb;
  /* the output array                         */
} Arg_s;



/****************************isign corresponding to +1  **************************/


/* the basic operation here is load up a spinor, do the spin decomposition, and store the halfspinor
to a lattice temporary */

void decomp_plus(size_t lo,size_t hi, int id, const void *ptr) /*need to fix decomp_minus */
{
  DECL_COMMON_ALIASES_TEMPS;
  const Arg_s *a = (Arg_s *)ptr;

  MY_SSE_VECTOR* chia = a->chifun; /* needs to be changed to MY_SSE_VECTOR and be an array*/
  MY_SSE_VECTOR* s3;
  MY_SSE_VECTOR* s4;
  int cb = a->cb;
  int  low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;

  MY_SPINOR* spinor_field= a->spinfun;
   
  iy1=iup(low,0);
  sp1=&spinor_field[low];


/************************ loop over all lattice sites *************************/
  for (ix1=low;ix1<high;ix1+=2) 
  {
    s1=&spinor_field[ix1+2];
    _prefetch_spinor(s1);

    /* prefetched input spinor for next two sites */
      
/******************************* direction +0 *********************************/
    /* ...(1-isign*gamma(0))... */
	  
    /* first of two sites */
    /* _load loads into xmm0-2, while _load_up loads into xmm3-5 */
    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    s3 = a_chia(0,idn(ix1,0));
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma0_minus();
    
    /* the halfspinor is now in xmm0-2 , so _store*/
      
    /* note: if the communications hardware likes its buffers out to memory instead of in cache, then non-temporal stores may be in order...check one of the #define switches above for more details */
    _sse_vector_store_NTA(*s3);
      
     /* second of two sites */
    _sse_pair_load((*(sp1+1))_c1__,(*(sp1+1))_c2__);
    s3 = a_chia(0,idn(ix1+1,0));
    _sse_pair_load_up((*(sp1+1))_c3__,(*(sp1+1))_c4__);
    _sse_42_gamma0_minus();

    _sse_vector_store_NTA(*s3);


/******************************* direction +1 *********************************/
    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    s3 = a_chia(1,idn(ix1,1));
    _sse_42_gamma1_minus();

    _sse_vector_store_NTA(*s3);
      
    _sse_pair_load((*(sp1+1))_c1__,(*(sp1+1))_c2__);
    s3 = a_chia(1,idn(ix1+1,1));
    _sse_pair_load_up((*(sp1+1))_c3__,(*(sp1+1))_c4__);
    _sse_42_gamma1_minus();

    _sse_vector_store_NTA(*s3);

/******************************* direction +2 *********************************/
    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    s3 = a_chia(2,idn(ix1,2));
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma2_minus();
      
    _sse_vector_store_NTA(*s3);

    _sse_pair_load((*(sp1+1))_c1__,(*(sp1+1))_c2__);
    s3 = a_chia(2,idn(ix1+1,2));
    _sse_pair_load_up((*(sp1+1))_c3__,(*(sp1+1))_c4__);
    _sse_42_gamma2_minus();

    _sse_vector_store_NTA(*s3);   	
	
    sp2=sp1+1;

/******************************* direction +3 *********************************/
    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    s3 = a_chia(3,idn(ix1,3));
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma3_minus();

    _sse_vector_store_NTA(*s3);

    _sse_pair_load((*sp2)_c1__,(*sp2)_c2__);
    s3 = a_chia(3,idn(ix1+1,3));
    _sse_pair_load_up((*sp2)_c3__,(*sp2)_c4__);
    _sse_42_gamma3_minus();
      
    _sse_vector_store_NTA(*s3);
	 
    iz1=ix1+2;

    sp1=&spinor_field[iz1];

/******************************** end of loop *********************************/
  }
}


/* the basic operations in this routine include loading a spinor, doing 
 * the spin projection, and multiplying the halfspinor by the appropriate 
 * gauge field, and saving the resulting halfspinor to a lattice temporary */

/* need gauge fields on opposite cb */
void decomp_hvv_plus(size_t lo,size_t hi, int id, const void *ptr)
{
  DECL_COMMON_ALIASES_TEMPS;
  const Arg_s *a = (const Arg_s *)ptr;
  MY_SPINOR* spinor_field = a->spinfun;

  MY_SSE_VECTOR* chib = a->chifun; /* a 1-d map of a 2-d array */
  MY_GAUGE_ARRAY gauge_field = a->u;
  MY_SSE_VECTOR* s3;
  MY_SSE_VECTOR* s4;

  int cb = a->cb;
  int  low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;

/************************ loop over all lattice sites *************************/
  for (ix1=low;ix1<high;ix1+=2) 
  {
/******************************* direction +0 *********************************/
    /* ...(1+gamma(0))... */
    sm1=&spinor_field[ix1];
    um1=a_gauge_field0_0(ix1); 
    um2=a_gauge_field0_1(ix1);
    
/******************************* direction -0 *********************************/
    /* ...(1-gamma(0))... */
    /* load the input spinor */  
    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    um3=a_gauge_field0_2(ix1);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);

    /* prefetch the next direction's worth of gauge fields */
    _prefetch_su3(um3);

    /* do the spin projection */
    _sse_42_gamma0_plus();

    s3 = a_chib(0,iup(ix1,0));
    /*_prefetch_single(s3);*/
    /* do the Hermitian conjugate multiplication */
    _sse_su3_inverse_multiply((*um1));
    /* store_up since the result is in xmm3-5 */
    _sse_vector_store_up(*s3);
      
    _sse_pair_load((*(sm1+1))_c1__,(*(sm1+1))_c2__);
    _sse_pair_load_up((*(sm1+1))_c3__,(*(sm1+1))_c4__);
     
    _sse_42_gamma0_plus();
    s4 = a_chib(0,iup(ix1+1,0));
    /*_prefetch_single(s4);*/

    _sse_su3_inverse_multiply((*(um2)));
    _sse_vector_store_up(*s4);

/******************************* direction +1 *********************************/
    um1=um3;
    um2=a_gauge_field0_3(ix1);
   
/******************************* direction -1 *********************************/
    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
	  
    um3=a_gauge_field1_0(ix1);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _prefetch_su3(um3);
    _sse_42_gamma1_plus();
	  
    s3 = a_chib(1,iup(ix1,1));
    /* _prefetch_single(s3);*/
    _sse_su3_inverse_multiply((*um1));
    _sse_vector_store_up(*s3);

    _sse_pair_load((*(sm1+1))_c1__,(*(sm1+1))_c2__);
    _sse_pair_load_up((*(sm1+1))_c3__,(*(sm1+1))_c4__);
    _sse_42_gamma1_plus();
    s4 = a_chib(1,iup(ix1+1,1));
    /*_prefetch_single(s4);*/
    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_store_up(*s4);


/******************************* direction +2 *********************************/
    um1=um3;
    um2=a_gauge_field1_1(ix1);
    
/******************************* direction -2 *********************************/
    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    um3=a_gauge_field1_2(ix1);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _prefetch_su3(um3);

    _sse_42_gamma2_plus();      

    s3 = a_chib(2,iup(ix1,2));
    /*_prefetch_single(s3);*/
    _sse_su3_inverse_multiply((*um1));

    _sse_vector_store_up(*s3);
      
    _sse_pair_load((*(sm1+1))_c1__,(*(sm1+1))_c2__);
    _sse_pair_load_up((*(sm1+1))_c3__,(*(sm1+1))_c4__);
    _sse_42_gamma2_plus();      

    s4 = a_chib(2,iup(ix1+1,2));
    /*_prefetch_single(s4);*/
    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_store_up(*s4);
     

/******************************* direction +3 *********************************/
    um1=um3;
    um2=a_gauge_field1_3(ix1);
    sm2=sm1+1;
	  
/******************************* direction -3 *********************************/
    iz1=ix1+2;
      
    s3 = a_chib(3,iup(ix1,3));
    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    sm1 = &spinor_field[iz1];

    _sse_42_gamma3_plus();
	  
    /*_prefetch_single(s3);*/
    _prefetch_spinor(sm1);
    _sse_su3_inverse_multiply((*um1));
    _sse_vector_store_up(*s3);
      
    um1=a_gauge_field0_0(iz1);  /* gauge packed or not this is the same */
    _prefetch_su3(um1);
      
    _sse_pair_load((*sm2)_c1__,(*sm2)_c2__);
    _sse_pair_load_up((*sm2)_c3__,(*sm2)_c4__);
    s4 = a_chib(3,iup(ix1+1,3));
    /*_prefetch_single(s4);*/

    _sse_42_gamma3_plus();
    _sse_su3_inverse_multiply((*um2));
    _sse_vector_store_up(*s4);

/******************************** end of loop *********************************/
  }
}
/***************end of decomp_hvv****************/


/* the basic operations in this routine include loading the halfspinor 
 * from memory, multiplying it by the appropriate gauge field, doing the 
 * spin reconstruction, and summing over directions, and saving the partial 
 * sum over directions */

void mvv_recons_plus(size_t lo,size_t hi, int id, const void *ptr)
{
  DECL_COMMON_ALIASES_TEMPS;

  const Arg_s *a =(Arg_s *)ptr;

  MY_SPINOR* spinor_field = a->spinfun;

  MY_SSE_VECTOR* chia = a->chifun; /* a 1-d map of a 2-d array */
  MY_GAUGE_ARRAY gauge_field = a->u;
  MY_SSE_VECTOR* s3;
  MY_SSE_VECTOR* s4;
  int cb = a->cb;
  int  low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;

  s3 = a_chia(0,mvv_shift(low,0));
  _prefetch_single_NTA(s3);
  iy1=iup(low,0);
  sp1=&spinor_field[iy1];
  up1=a_gauge_field0_0(low);
  up2=a_gauge_field0_1(low);

/************************ loop over all lattice sites *************************/
  for (ix1=low;ix1<high;ix1+=2) 
  {
    /* s1=&spinor_field[ix1];
       _prefetch_spinor(s1);*/
      
/******************************* direction +0 *********************************/
    /* ...(1-gamma(0))... */
    /* load from the temporary */
    _sse_vector_load(*s3);
    s4 = a_chia(0,mvv_shift(ix1+1,0));
    
    /*prefetch the next temp into one way of lvl2 cache, and prefetch the next gague field */
    _prefetch_single_NTA(s4);
    _prefetch_su3(up1+2); 
      /* do the SU(3) normal multiplication */
    _sse_su3_multiply((*up1));
       
    /* the first two rows of that 4x2 reconstruction matrix are the identity 2x2 matrix, 
     * so just _store_up */
    /* notation: r12_1 -- first two spin components, first of the two  adjacent sites 
     * r34_1 -- last two spin components, first site
     * r12_2 -- first two spin components, adjacent (ix1+1) site
     * r34_2 -- last two spin components, adjacent site
     */
    _sse_vector_store_up(r12_1);
      
    /* do some kind of multiplication to form the bottom two components from the top two */
    _sse_24_gamma0_minus_set(); 
    /*answer is still in xmm3-5, so store_up */
    _sse_vector_store_up(r34_1);

    /* now do the other sites's worth */
    _sse_vector_load(*s4);
    s3 = a_chia(1,mvv_shift(ix1,1));
    _prefetch_single_NTA(s3);
    _sse_su3_multiply((*(up2)));
    _sse_vector_store_up(r12_2);

    _sse_24_gamma0_minus_set(); 
    _sse_vector_store_up(r34_2);

/***************************** direction +1 ***********************************/
    up1=a_gauge_field0_2(ix1);
    up2=a_gauge_field0_3(ix1);   
      
    /* same kind of thing again */
    _sse_vector_load(*s3);
    s4 = a_chia(1,mvv_shift(ix1+1,1));

    _prefetch_single_NTA(s4);
    _prefetch_su3(up1+2); 

    _sse_su3_multiply((*up1));

    /* ok here's where things are different..now we load the partial sum over 
     * directions of the output spinor, then add on to them the spin 
     * reconstructed terms (using the bottom two rows of the 4x2 reconstruction
     * matrix for the spin reconstruction), and then save the temp
     * and it should stay in cache */

    /* remember, we load and store partial sums with xmm0-2, with the multiplied 
     * halfspinor in xmm3-5 */
    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);
      
    _sse_vector_load(r34_1);
    _sse_24_gamma1_minus();
    _sse_vector_store(r34_1);

    /* same kind of stuff from here on out until direction 3*/  
    _sse_vector_load(*s4);
    s3 = a_chia(2,mvv_shift(ix1,2));
    _prefetch_single_NTA(s3);

    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);
      
    _sse_vector_load(r34_2);
    _sse_24_gamma1_minus();
    _sse_vector_store(r34_2);

    up1=a_gauge_field1_0(ix1); /* default is packed gauge fields, thats why the 
				* statements are non-intuitive */
    up2=a_gauge_field1_1(ix1);

/******************************* direction +2 *********************************/
    _sse_vector_load(*s3);
    s4 = a_chia(2,mvv_shift(ix1+1,2));
    _prefetch_single_NTA(s4);
    _prefetch_su3(up1+2); 

    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);       

    _sse_vector_load(r34_1);
    _sse_24_gamma2_minus();
    _sse_vector_store(r34_1);       

    _sse_vector_load(*s4);
    s3 = a_chia(3,mvv_shift(ix1,3));
    _prefetch_single_NTA(s3);
     
    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);       

    _sse_vector_load(r34_2);
    _sse_24_gamma2_minus();
    _sse_vector_store(r34_2); 

    up1=a_gauge_field1_2(ix1);
    up2=a_gauge_field1_3(ix1);

/******************************* direction +3 *********************************/
    sn1=&spinor_field[ix1];

    _sse_vector_load(*s3);
    s4 = a_chia(3,mvv_shift(ix1+1,3));
    
    _prefetch_single_NTA(s4);
    /* prefetch the gague field for direction 0, site = ix1+2 */
    _prefetch_su3(up1+2); 
      
    _sse_su3_multiply((*up1));
    
    /* ok things get tricky here again...the rows12 thing is for compatibility 
     * with other spin basis in which the first two rows of that 4x2 
     * reconstruction matrix may be something other than the 2x2 identity matrix */
    _sse_vector_load(r12_1);
    _sse_24_gamma3_minus_rows12();
    _sse_pair_store((*sn1)_c1__,(*sn1)_c2__);

    _sse_vector_load(r34_1);
    _sse_24_gamma3_minus();
    _sse_pair_store((*sn1)_c3__,(*sn1)_c4__);   

    iz1=ix1+2;
    if (iz1==high)
      iz1=0;

    _sse_vector_load(*s4);
    s3 = a_chia(0,mvv_shift(iz1,0));
    _prefetch_single_NTA(s3);
      
    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_24_gamma3_minus_rows12();
    _sse_pair_store((*(sn1+1))_c1__,(*(sn1+1))_c2__);

    _sse_vector_load(r34_2);
    _sse_24_gamma3_minus();
    _sse_pair_store((*(sn1+1))_c3__,(*(sn1+1))_c4__); 

    up1=a_gauge_field0_0(iz1);
	 up2=a_gauge_field0_1(iz1);

/******************************** end of loop *********************************/
  }
}
/******************end of mvv_recons*************************/

   
/* this routine takes the partial sum from mvv_recons() and loops 
 * over the output spin components, 2 at a time doing a sum over directions 
 * for each set, accumulating in xmm0-2 and loading the halfspinor 
 * temporaries into xmm3-5 */

void recons_plus(size_t lo,size_t hi, int id, const void *ptr )	
{
  DECL_COMMON_ALIASES_TEMPS;
  const Arg_s *a = (Arg_s *)ptr;
  MY_SPINOR* spinor_field = a->spinfun;

  MY_SSE_VECTOR* hs0, *hs1, *hs2,*hs3,*hs4,*hs5,*hs6,*hs7,*hs8;
  MY_SSE_VECTOR* chib = a->chifun; /* a 1-d map of a 2-d array */
  int cb = a->cb;
  int low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;

  /* printf("\nlo:%i, hi:%i, id:%i, chib:%x", lo, hi, id, chib);*/

/************************ loop over all lattice sites *************************/
  for (ix1=low;ix1<high;ix1+=2) 
  {
#define PREFDIST 4
    hs0 = a_chib(0,rec_shift(ix1,0));
    sn1=&spinor_field[ix1];   
   
    _prefetch_spinor_NTA(sn1+PREFDIST);
    /* _prefetch_single(&iup(ix1,0)+16);*/
     
    /***** psi 1&2 site 1 ******/
    _sse_vector_load_up(*(hs0));   /* vector in xmm3-5 */
	                             
    /* accumulate in xmm0-2 */
    _sse_pair_load((*sn1)_c1__, (*sn1)_c2__); /*load in partial sum */
	  
    hs1 = a_chib(1,rec_shift(ix1,1));
    _sse_vector_add();
    _sse_vector_load_up(*(hs1)); /*direction +1 */ 
	 
    hs2 = a_chib(2,rec_shift(ix1,2));
    _sse_vector_add();          /* accumulating in xmm0-2 */
    _sse_vector_load_up(*(hs2));  /* direction +2 */
	
    hs3 = a_chib(3,rec_shift(ix1,3));
    _sse_vector_add();
    _sse_vector_load_up(*(hs3)); /* direction +3 */
   
    _sse_24_gamma3_plus_rows12();
    _sse_pair_store((*sn1)_c1__,(*sn1)_c2__);

    /***** psi's 3 and 4 sites 1*****/
    /*hs2= a_chib(0,rec_shift(ix1+PREFDIST,0));
    _prefetch_single(hs2);*/
    _sse_vector_load_up(*(hs0));       /* vector in xmm3-5 */
   
    /* accumulate in xmm0-2 */
    /*hs2 = a_chib(1,rec_shift(ix1+PREFDIST,1));
    _prefetch_single(hs2);*/
    _sse_pair_load((*sn1)_c3__,(*sn1)_c4__); /*load in partial sum */

    _sse_24_gamma0_plus_add();

    /*hs2 = a_chib(2,rec_shift(ix1+PREFDIST,2));
    _prefetch_single(hs2);*/
    _sse_vector_load_up(*(hs1));
    
    _sse_24_gamma1_plus();

    /*hs2 = a_chib(3,rec_shift(ix1+PREFDIST,3));
    _prefetch_single(hs2);*/
    _sse_vector_load_up(*(hs2));
   
    _sse_24_gamma2_plus();
      
    _sse_vector_load_up(*(hs3));
    _sse_24_gamma3_plus();
       
    _sse_pair_store((*sn1)_c3__,(*sn1)_c4__);

    _prefetch_single(hs2);


    /****** psi 1&2 site 2 ******/
    hs0 = a_chib(0,rec_shift(ix1+1,0)); 
    _sse_vector_load_up(*(hs0));   /* vector in xmm3-5 */
    hs1 = a_chib(1,rec_shift(ix1+1,1));
	  
    /* accumulate in xmm0-2 */
    _sse_pair_load((*(sn1+1))_c1__,(*(sn1+1))_c2__); /*load in partial sum */
   
    _sse_vector_add();
    _sse_vector_load_up(*(hs1)); /*direction +1 */
    hs2 = a_chib(2,rec_shift(ix1+1,2));
      
    _sse_vector_add();          /* accumulating in xmm0-2 */
	  
    _sse_vector_load_up(*(hs2)); /* direction +2 */
    hs3 = a_chib(3,rec_shift(ix1+1,3));

    _sse_vector_add();

    _sse_vector_load_up(*(hs3)); /* direction +3 */
     
    _sse_24_gamma3_plus_rows12();

    _sse_pair_store((*(sn1+1))_c1__,(*(sn1+1))_c2__);
	 
      
    /***** psi's 3 and 4 site 2 *****/
    /*hs2 = a_chib(0,rec_shift(ix1+1+PREFDIST,0));
    _prefetch_single(hs2);*/

    _sse_vector_load_up(*(hs0));       /* vector in xmm3-5 */
     
    _sse_pair_load((*(sn1+1))_c3__,(*(sn1+1))_c4__); /*load in partial sum */
    _sse_24_gamma0_plus_add();

    /*hs2 = a_chib(1,rec_shift(ix1+1+PREFDIST,1));
    _prefetch_single(hs2);*/
    _sse_vector_load_up(*(hs1));
    
    _sse_24_gamma1_plus();
    /*hs2 = a_chib(2,rec_shift(ix1+1+PREFDIST,2));
    _prefetch_single(hs2);*/
    _sse_vector_load_up(*(hs2));
    
    _sse_24_gamma2_plus();
    /*hs2 = a_chib(3,rec_shift(ix1+1+PREFDIST,3));
    _prefetch_single(hs2);*/
    _sse_vector_load_up(*(hs3));
    _sse_24_gamma3_plus();
   
    _sse_pair_store((*(sn1+1))_c3__,(*(sn1+1))_c4__); 
  
    /*************************end of loop ****************************/
  }
}
/*****************end of recons**************/





/*************** now for isign corresponding to -1  ****************************************/

void decomp_minus(size_t lo,size_t hi, int id, const void *ptr ) /*need to fix decomp_minus */
{
  DECL_COMMON_ALIASES_TEMPS;
  const Arg_s *a =(Arg_s *)ptr;

  MY_SSE_VECTOR* chia = a->chifun; /* needs to be changed to MY_SSE_VECTOR and be an array*/
  MY_SSE_VECTOR* s3;
  MY_SSE_VECTOR* s4;
  int cb = a->cb;
  int  low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;

  MY_SPINOR* spinor_field= a->spinfun;
   
  /*printf("\nlo:%i, hi:%i, id:%i, chia[0][0]:%x", lo, hi, id, chia[0][0]);*/
  iy1=iup(low,0);
  sp1=&spinor_field[low];
 
/************************ loop over all lattice sites *************************/

  for (ix1=low;ix1<high;ix1+=2) 
  {
    s1=&spinor_field[ix1+2];
    _prefetch_spinor(s1);
      
/******************************* direction +0 *********************************/
    /* ...(1-gamma(0))... */
    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    s3 = a_chia(0,idn(ix1,0));
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma0_plus();
    _sse_vector_store(*s3);
      
    _sse_pair_load((*(sp1+1))_c1__,(*(sp1+1))_c2__);
    s4 = a_chia(0,idn(ix1+1,0));
    _sse_pair_load_up((*(sp1+1))_c3__,(*(sp1+1))_c4__);
    _sse_42_gamma0_plus();
    _sse_vector_store(*s4);

/******************************* direction -0 *********************************/
   

/******************************* direction +1 *********************************/
    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    s3 = a_chia(1,idn(ix1,1));
    _sse_42_gamma1_plus();
    _sse_vector_store(*s3);
      
    _sse_pair_load((*(sp1+1))_c1__,(*(sp1+1))_c2__);
    s4 = a_chia(1,idn(ix1+1,1));
    _sse_pair_load_up((*(sp1+1))_c3__,(*(sp1+1))_c4__);
    _sse_42_gamma1_plus();
    _sse_vector_store(*s4);

/******************************* direction -1 *********************************/
     

/******************************* direction +2 *********************************/
    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    s3 = a_chia(2,idn(ix1,2));
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma2_plus();
    _sse_vector_store(*s3);

    _sse_pair_load((*(sp1+1))_c1__,(*(sp1+1))_c2__);
    s4 = a_chia(2,idn(ix1+1,2));
    _sse_pair_load_up((*(sp1+1))_c3__,(*(sp1+1))_c4__);
    _sse_42_gamma2_plus();
    _sse_vector_store(*s4);

/******************************* direction -2 *********************************/
    sp2=sp1+1;

/******************************* direction +3 *********************************/
    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    s3 = a_chia(3,idn(ix1,3));
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma3_plus();
    _sse_vector_store(*s3);
      
    _sse_pair_load((*sp2)_c1__,(*sp2)_c2__);
    s4 = a_chia(3,idn(ix1+1,3));
    _sse_pair_load_up((*sp2)_c3__,(*sp2)_c4__);
    _sse_42_gamma3_plus();
    _sse_vector_store(*s4);

/******************************* direction -3 *********************************/
    iz1=ix1+2;
    sp1=&spinor_field[iz1];

/******************************** end of loop *********************************/
  }
}


/* need gauge fields on opposite cb */
void decomp_hvv_minus(size_t lo,size_t hi, int id, const void *ptr )
{
  DECL_COMMON_ALIASES_TEMPS;
  const Arg_s *a =(Arg_s *)ptr;
  MY_SPINOR* spinor_field = a->spinfun;
  MY_SSE_VECTOR* chib = a->chifun; /* a 1-d map of a 2-d array */
  MY_GAUGE_ARRAY gauge_field = a->u;
  MY_SSE_VECTOR* s3;
  MY_SSE_VECTOR* s4;

  int cb = a->cb;
  int  low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;


  /*  printf("\nlo:%i, hi:%i, id:%i, chib:%x", lo, hi, id, chib[0][0]);*/

/************************ loop over all lattice sites *************************/
  for (ix1=low;ix1<high;ix1+=2) 
  {
/******************************* direction +0 *********************************/
    /* ...(1+gamma(0))... */
    sm1=&spinor_field[ix1];
    um1=a_gauge_field0_0(ix1); 
    um2=a_gauge_field0_1(ix1);
    
    /* printf("gauge0_0(ix1)=%f,gauge0_1(ix1)=%f\n",(*um1)[0][0][0],(*um2)[0][0][0]); */
/******************************* direction -0 *********************************/
    /* ...(1-gamma(0))... */
    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    um3=a_gauge_field0_2(ix1);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _prefetch_su3(um3);
	  
    _sse_42_gamma0_minus();
    s3 = a_chib(0,iup(ix1,0));
    _prefetch_single(s3);

    _sse_su3_inverse_multiply((*um1));
    _sse_vector_store_up(*s3);

    _sse_pair_load((*(sm1+1))_c1__,(*(sm1+1))_c2__);
    _sse_pair_load_up((*(sm1+1))_c3__,(*(sm1+1))_c4__);
	  
    _sse_42_gamma0_minus();
    s4 = a_chib(0,iup(ix1+1,0));
    _prefetch_single(s4);
    _sse_su3_inverse_multiply((*(um2)));
    _sse_vector_store_up(*s4);

    /*printf("chi(ix1)=%f,chi(ix1+1) [0][1][1]=%f\n",(*s3)[0][1][1],(*s4)[0][1][1]);
      printf("chi(ix1)=%f,chi(ix1+1) [0][2][0]=%f\n",(*s3)[0][2][0],(*s4)[0][2][0]);*/
      
/******************************* direction +1 *********************************/
    um1=um3;
    um2=a_gauge_field0_3(ix1);
   
/******************************* direction -1 *********************************/
    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
	  
    um3=a_gauge_field1_0(ix1);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _prefetch_su3(um3);
    _sse_42_gamma1_minus();
	  
    s3 = a_chib(1,iup(ix1,1));
    _prefetch_single(s3);
    _sse_su3_inverse_multiply((*um1));

    _sse_vector_store_up(*s3);

    _sse_pair_load((*(sm1+1))_c1__,(*(sm1+1))_c2__);
    _sse_pair_load_up((*(sm1+1))_c3__,(*(sm1+1))_c4__);
    _sse_42_gamma1_minus();
    s4 = a_chib(1,iup(ix1+1,1));
    _prefetch_single(s4);
    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_store_up(*s4);

/******************************* direction +2 *********************************/
    um1=um3;
    um2=a_gauge_field1_1(ix1);

/******************************* direction -2 *********************************/
    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    um3=a_gauge_field1_2(ix1);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _prefetch_su3(um3);

    _sse_42_gamma2_minus();      

    s3 = a_chib(2,iup(ix1,2));
    _prefetch_single(s3);
    _sse_su3_inverse_multiply((*um1));

    _sse_vector_store_up(*s3);
      
    _sse_pair_load((*(sm1+1))_c1__,(*(sm1+1))_c2__);
    _sse_pair_load_up((*(sm1+1))_c3__,(*(sm1+1))_c4__);
    _sse_42_gamma2_minus();      

    s4 = a_chib(2,iup(ix1+1,2));
    _prefetch_single(s4);
    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_store_up(*s4);
     
/******************************* direction +3 *********************************/
    um1=um3;
    um2=a_gauge_field1_3(ix1);
    sm2=sm1+1;

/******************************* direction -3 *********************************/
    iz1=ix1+2;
    if (iz1==high)
      iz1=0;
      
    s3 = a_chib(3,iup(ix1,3));
    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    sm1 = &spinor_field[iz1];

    _sse_42_gamma3_minus();
	  
    _prefetch_single(s3);
    _prefetch_spinor(sm1);

    _sse_su3_inverse_multiply((*um1));
    _sse_vector_store_up(*s3);
      
    um1=a_gauge_field0_0(iz1);  /* gauge packed or not this is the same */
    _prefetch_su3(um1);
      
    _sse_pair_load((*sm2)_c1__,(*sm2)_c2__);
    _sse_pair_load_up((*sm2)_c3__,(*sm2)_c4__);
    s4 = a_chib(3,iup(ix1+1,3));
    _prefetch_single(s4);

    _sse_42_gamma3_minus();
    _sse_su3_inverse_multiply((*um2));
    _sse_vector_store_up(*s4);
	  
/******************************** end of loop *********************************/
  }
}
/***************end of decomp_hvv****************/


void mvv_recons_minus(size_t lo,size_t hi, int id, const void *ptr )
{
  DECL_COMMON_ALIASES_TEMPS;
  /* if going to support unpacked gauge fields, need to treat site ix1 and site ix1+1 separately */
  /* to support unpacked gauge fields the prefetches will need to be changed */
  const Arg_s *a =(Arg_s *)ptr;
  MY_SPINOR* spinor_field = a->spinfun;
  MY_SSE_VECTOR* chia = a->chifun; /* a 1-d map of a 2-d array */
  MY_GAUGE_ARRAY gauge_field = a->u;
  MY_SSE_VECTOR* s3;
  MY_SSE_VECTOR* s4;
  int cb = a->cb;
  int  low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;

  s3 = a_chia(0,mvv_shift(low,0));
  _prefetch_single_NTA(s3);
  iy1=iup(low,0);
  sp1=&spinor_field[iy1];
  up1=a_gauge_field0_0(low);
  up2=a_gauge_field0_1(low);

/************************ loop over all lattice sites *************************/
  for (ix1=low;ix1<high;ix1+=2) 
  {
/******************************* direction +0 *********************************/
    /* ...(1-isign*gamma(0))... */
    /*printf("gaugB0_0(ix1)=%f,gauge0_1(ix1)=%f\n",(*up1)[0][0][0],(*up2)[0][0][0]);*/
    _sse_vector_load(*s3);
    s4 = a_chia(0,mvv_shift(ix1+1,0));
    /*printf("chia(ix1)=%f,chia(ix1+1) [0][1][1]=%f\n",(*s3)[0][1][1],(*s4)[0][1][1]);
      printf("chia(ix1)=%f,chia(ix1+1) [0][2][0]=%f\n",(*s3)[0][2][0],(*s4)[0][2][0]);*/
    _prefetch_single_NTA(s4);
    _prefetch_su3(up1+2); 
      
    _sse_su3_multiply((*up1));
    _sse_vector_store_up(r12_1);

    _sse_24_gamma0_plus_set(); 
    _sse_vector_store_up(r34_1);
    
    _sse_vector_load(*s4);
    s3 = a_chia(1,mvv_shift(ix1,1));
    _prefetch_single_NTA(s3);
    _sse_su3_multiply((*(up2)));
    _sse_vector_store_up(r12_2);

    _sse_24_gamma0_plus_set(); 
    _sse_vector_store_up(r34_2);
      
/******************************* direction -0 *********************************/
    /* ...(1-gamma(0))... */
    up1=a_gauge_field0_2(ix1);
    up2=a_gauge_field0_3(ix1);   
      
/******************************* direction +1 *********************************/
    _sse_vector_load(*s3);
    s4 = a_chia(1,mvv_shift(ix1+1,1));
    _prefetch_single_NTA(s4);
    _prefetch_su3(up1+2); 

    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);
      
    _sse_vector_load(r34_1);
    _sse_24_gamma1_plus();
    _sse_vector_store(r34_1);

      
    _sse_vector_load(*s4);
    s3 = a_chia(2,mvv_shift(ix1,2));
    _prefetch_single_NTA(s3);

    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);
      
    _sse_vector_load(r34_2);
    _sse_24_gamma1_plus();
    _sse_vector_store(r34_2);

/******************************* direction -1 *********************************/
    up1=a_gauge_field1_0(ix1); /* default is packed gauge fields, thats 
				* why the statements are non-intuitive */
    up2=a_gauge_field1_1(ix1);

/******************************* direction +2 *********************************/
    _sse_vector_load(*s3);
    s4 = a_chia(2,mvv_shift(ix1+1,2));
    _prefetch_single_NTA(s4);
    _prefetch_su3(up1+2); 

    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);       

    _sse_vector_load(r34_1);
    _sse_24_gamma2_plus();
    _sse_vector_store(r34_1);       

    _sse_vector_load(*s4);
    s3 = a_chia(3,mvv_shift(ix1,3));
    _prefetch_single_NTA(s3);

    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);       

    _sse_vector_load(r34_2);
    _sse_24_gamma2_plus();
    _sse_vector_store(r34_2); 

/******************************* direction -2 *********************************/
    up1=a_gauge_field1_2(ix1);
    up2=a_gauge_field1_3(ix1);

/******************************* direction +3 *********************************/
    sn1=&spinor_field[ix1];

    _sse_vector_load(*s3);
    s4 = a_chia(3,mvv_shift(ix1+1,3));
    _prefetch_single_NTA(s4);
    _prefetch_su3(up1+2); 
      
    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);
    _sse_24_gamma3_plus_rows12();
    _sse_pair_store((*sn1)_c1__,(*sn1)_c2__);

    _sse_vector_load(r34_1);
    _sse_24_gamma3_plus();
    _sse_pair_store((*sn1)_c3__,(*sn1)_c4__);   

    iz1=ix1+2;
    if (iz1==high)
      iz1=0;

    _sse_vector_load(*s4);
    s3 = a_chia(0,mvv_shift(iz1,0));
    _prefetch_single_NTA(s3);
      
    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_24_gamma3_plus_rows12();
    _sse_pair_store((*(sn1+1))_c1__,(*(sn1+1))_c2__);

    _sse_vector_load(r34_2);
    _sse_24_gamma3_plus();
    _sse_pair_store((*(sn1+1))_c3__,(*(sn1+1))_c4__); 

/******************************* direction -3 *********************************/
    up1=a_gauge_field0_0(iz1);
    up2=a_gauge_field0_1(iz1);

/******************************** end of loop *********************************/
  }
}
/******************end of mvv_recons*************************/


void recons_minus(size_t lo,size_t hi, int id, const void *ptr )	
{
  DECL_COMMON_ALIASES_TEMPS;
  const Arg_s *a = (Arg_s *)ptr;
  MY_SPINOR* spinor_field = a->spinfun;
  MY_SSE_VECTOR* hs0, *hs1, *hs2, *hs3, *hs4, *s2;
  MY_SSE_VECTOR* chib = a->chifun; /* a 1-d map of a 2-d array */
  int cb = a->cb;
  int low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;

/************************ loop over all lattice sites *************************/
  for (ix1=low;ix1<high;ix1+=2) 
  {
/******************************* direction +0 *********************************/
#define PREFDIST 4

    /* I just played around with PREFDIST to get a good value after getting a 
     * starting value from the prefetch scheduling distance function thats in
     * one of the intel manuals ...feel free to do the same if you want */
    hs0 = a_chib(0,rec_shift(ix1,0));
    sn1=&spinor_field[ix1];   
    
    _prefetch_spinor_NTA(sn1+PREFDIST);
    
    /* need to do psum[ix1][0] first for all +mu */
    /* loop over our two adjacent sites slowest inside here */

    /***** psi 1&2 site 1 ******/
    _sse_vector_load_up(*(hs0));   /* vector in xmm3-5 */
     
    /* accumulate in xmm0-2 */
    _sse_pair_load((*sn1)_c1__, (*sn1)_c2__); /*load in partial sum */
	  
    hs1 = a_chib(1,rec_shift(ix1,1));
    _sse_vector_add();
    _sse_vector_load_up(*(hs1)); /*direction +1 */ 
	 
    hs2 = a_chib(2,rec_shift(ix1,2));
    _sse_vector_add();          /* accumulating in xmm0-2 */
    _sse_vector_load_up(*(hs2));  /* direction +2 */
	
    hs3 = a_chib(3,rec_shift(ix1,3));
    _sse_vector_add();
    _sse_vector_load_up(*(hs3)); /* direction +3 */
    
    _sse_24_gamma3_minus_rows12();  /* here's that thing that doesn't 
				     * reduce to an _add in some spin basis */
    _sse_pair_store((*sn1)_c1__,(*sn1)_c2__);

    /***** psi's 3 and 4 sites 1*****/
    /*s2= a_chib(0,rec_shift(ix1+PREFDIST,0));
    _prefetch_single(s2);*/
   
    _sse_vector_load_up(*(hs0));       /* vector in xmm3-5 */
        
    /* accumulate in xmm0-2 */
    /*s2 = a_chib(1,rec_shift(ix1+PREFDIST,1));
    _prefetch_single(s2);*/

    _sse_pair_load((*sn1)_c3__,(*sn1)_c4__); /*load in partial sum */
	   
    _sse_24_gamma0_minus_add();

    /*s2 = a_chib(2,rec_shift(ix1+PREFDIST,2));
    _prefetch_single(s2);*/

    _sse_vector_load_up(*(hs1));
    _sse_24_gamma1_minus();

    /*s2 = a_chib(3,rec_shift(ix1+PREFDIST,3));
    _prefetch_single(s2);*/

    _sse_vector_load_up(*(hs2));
    _sse_24_gamma2_minus();
    hs0 = a_chib(0,rec_shift(ix1+1,0)); 
      
    _sse_vector_load_up(*(hs3));
    _sse_24_gamma3_minus();
       
    _sse_pair_store((*sn1)_c3__,(*sn1)_c4__);
    _prefetch_single(s2);

    /****** psi 1&2 site 2 ******/
    _sse_vector_load_up(*(hs0));   /* vector in xmm3-5 */
    hs1 = a_chib(1,rec_shift(ix1+1,1));
    _sse_pair_load((*(sn1+1))_c1__,(*(sn1+1))_c2__); /*load in partial sum */
    
    _sse_vector_add();
    _sse_vector_load_up(*(hs1)); /*direction +1 */
    hs2 = a_chib(2,rec_shift(ix1+1,2));
	  
    _sse_vector_add();          /* accumulating in xmm0-2 */
    _sse_vector_load_up(*(hs2)); /* direction +2 */
    hs3 = a_chib(3,rec_shift(ix1+1,3));
	   
    _sse_vector_add();
    _sse_vector_load_up(*(hs3)); /* direction +3 */
     
    _sse_24_gamma3_minus_rows12();
    _sse_pair_store((*(sn1+1))_c1__,(*(sn1+1))_c2__);
	 
    /***** psi's 3 and 4 site 2 *****/
    /*s2 = a_chib(0,rec_shift(ix1+1+PREFDIST,0));
    _prefetch_single(s2);*/

    _sse_vector_load_up(*(hs0));       /* vector in xmm3-5 */
  
    /* accumulate in xmm0-2 */
    _sse_pair_load((*(sn1+1))_c3__,(*(sn1+1))_c4__); /*load in partial sum */
    _sse_24_gamma0_minus_add();

    /*s2 = a_chib(1,rec_shift(ix1+1+PREFDIST,1));
    _prefetch_single(s2);*/
    _sse_vector_load_up(*(hs1));
  
    _sse_24_gamma1_minus();
    /*s2 = a_chib(2,rec_shift(ix1+1+PREFDIST,2));
    _prefetch_single(s2);*/
    _sse_vector_load_up(*(hs2));
   
    _sse_24_gamma2_minus();
    /*s2 = a_chib(3,rec_shift(ix1+1+PREFDIST,3));
    _prefetch_single(s2);*/
    _sse_vector_load_up(*(hs3));
    _sse_24_gamma3_minus();

    _sse_pair_store((*(sn1+1))_c3__,(*(sn1+1))_c4__); 
    /*************************end of loop ****************************/
  }
}
/*****************end of recons_minus**************/



/*****************end of isign corresponding to -1 **************************************/


/***************** start of initialization routine ***************************************/

#define Nd 4
#define Nc 3
#define Ns 4
#define Ns2 2

static MY_SSE_VECTOR* chi1;
static MY_SSE_VECTOR* chi2;
#define TAIL1(chi,mymu) (chi+subgrid_vol_cb*(1+3*mymu))
#define TAIL2(chi,mymu) (chi+subgrid_vol_cb*(2+3*mymu))

static QMP_mem_t* xchi1;
static QMP_mem_t* xchi2;


/* Nearest neighbor communication channels */
static int total_comm = 0;
static QMP_msgmem_t forw_msg[Nd][2];
static QMP_msgmem_t back_msg[Nd][2];
static QMP_msghandle_t forw_mh[Nd][2];
static QMP_msghandle_t back_mh[Nd][2];
static QMP_msghandle_t forw_all_mh;
static QMP_msghandle_t back_all_mh;


void init_sse_su3dslash(const int latt_size[])   // latt_size not used, here for scalar version
{
  const int *machine_size = QMP_get_logical_dimensions();
  const int *subgrid_cb_size = QMP_get_subgrid_dimensions();
  int bound[2][4][Nd];
  int mu, num, nsize;

  /* If we are already initialised, then increase the refcount and return */
  if (initP > 0) 
  {
    initP++;
    return;
  }


  /* Otherwise initialise */
  if (QMP_get_logical_number_of_dimensions() != Nd)
  {
    QMP_error("init_sse_su3dslash: number of logical dimensions does not match problem");
    QMP_abort(1);
  }
    

  /* Check problem size */
  for(mu=0; mu < Nd; mu++) 
    if ( latt_size[mu] == 1 ) 
    {
      QMP_error("This SSE Dslash does not support a problem size = 1. Here the lattice in dimension %d has length %d\n", mu, latt_size[mu]);
      QMP_abort(1);
    }


  num = latt_size[0] / machine_size[0];
  if ( num & 1 != 0 )
  {
    QMP_error("This SSE Dslash does not work for odd x-sublattice. Here the sublattice is odd in dimension 0 with length %d\n", num);
    QMP_abort(1);
  }


  subgrid_vol_cb = subgrid_cb_size[0];
  for(mu = 1; mu < Nd; mu++)
    subgrid_vol_cb *= subgrid_cb_size[mu];
    
  subgrid_vol = subgrid_vol_cb << 1;


  // The SSE code expects to have at least 2 sites after checkerboarding.
  if ( subgrid_vol_cb <= 1 )
  {
    QMP_error("This SSE Dslash expects there to be at least 2 subgrid sites after checkerboarding");
    QMP_abort(1);
  }


  /* Allocated space for the floating temps */
  /* Wasteful - allocate 3 times subgrid_vol_cb. Otherwise, need to pack the TAIL{1,2} offsets */
  nsize = 2*Nc*Ns2*sizeof(float)*subgrid_vol_cb*3*Nd;  /* Note 3x4 half-ferm temps */
  if ((xchi1 = QMP_allocate_aligned_memory(nsize,128,0)) == 0)
  {
    QMP_error("init_wnxtsu3dslash: could not initialize xchi1");
    QMP_abort(1);
  }
  if ((xchi2 = QMP_allocate_aligned_memory(nsize,128,0)) == 0)
  {
    QMP_error("init_wnxtsu3dslash: could not initialize xchi2");
    QMP_abort(1);
  }
    
  chi1 = (MY_SSE_VECTOR*)QMP_get_memory_pointer(xchi1);
  chi2 = (MY_SSE_VECTOR*)QMP_get_memory_pointer(xchi2); 
    
  /* Construct all the shift tables needed */
  /* Use malloc here: the aligned mem might be pinned */
  if ((newshift = (int *)malloc(Nd*subgrid_vol*4*sizeof(int))) == 0)
  {
    QMP_error("init_wnxtsu3dslash: could not initialize newshift");
    QMP_abort(1);
  }
    
  make_shift_tables(newshift, icolor_start, bound);
  icolor_end[0] = icolor_start[0] + subgrid_vol_cb;
  icolor_end[1] = icolor_start[1] + subgrid_vol_cb;
  
  /* Loop over all communicating directions and build up the two message
   * handles. If there is no communications, the message handles will not
   * be initialized 
   */
  num = 0;
    
  for(mu=0; mu < Nd; ++mu) 
  {
    if(machine_size[mu] > 1) 
    {
      if (bound[0][0][mu] == 0)
      {
	QMP_error("init_sse_dslash: type 0 message size is 0");
	QMP_abort(1);
      }

      forw_msg[num][0] = QMP_declare_msgmem(TAIL1(chi1,mu), bound[0][0][mu]*sizeof(MY_SSE_VECTOR));
      forw_msg[num][1] = QMP_declare_msgmem(TAIL2(chi1,mu), bound[0][0][mu]*sizeof(MY_SSE_VECTOR));
      forw_mh[num][0]  = QMP_declare_receive_relative(forw_msg[num][1], mu, +1, 0);
      forw_mh[num][1]  = QMP_declare_send_relative(forw_msg[num][0], mu, -1, 0);
	
      if (bound[0][1][mu] == 0)
      {
	QMP_error("init_sse_dslash: type 0 message size is 0");
	QMP_abort(1);
      }

      back_msg[num][0] = QMP_declare_msgmem(TAIL1(chi2,mu), bound[0][1][mu]*sizeof(MY_SSE_VECTOR));
      back_msg[num][1] = QMP_declare_msgmem(TAIL2(chi2,mu), bound[0][1][mu]*sizeof(MY_SSE_VECTOR));
      back_mh[num][0]  = QMP_declare_receive_relative(back_msg[num][1], mu, -1, 0);
      back_mh[num][1]  = QMP_declare_send_relative(back_msg[num][0], mu, +1, 0);
	
      num++;
    }
  }

  if (num > 0) {
    forw_all_mh = QMP_declare_multiple(&(forw_mh[0][0]), 2*num);
    back_all_mh = QMP_declare_multiple(&(back_mh[0][0]), 2*num);
  }
  
  total_comm = num;
  initP = 1;
}

void free_sse_su3dslash(void)
{
  const int *machine_size = QMP_get_logical_dimensions();
  int mu, num;

  /***************HACK********************/
  /* There is a gigE QMP_free bug  (2/6/04). For now, turn off ever free-ing !! */
  return;


  /* If we are uninitialised just return */
  if (initP == 0) {
    return;
  }

  /* Otherwise decrease the refcount */
  initP--;

  /* If the refcount has now hit 0 then free stuff */
  if( initP == 0 ) { 

    /* Free space */
    QMP_free_memory(xchi1);
    QMP_free_memory(xchi2);
    free(newshift);
    
    if (total_comm > 0) {
      
      QMP_free_msghandle(forw_all_mh);
      QMP_free_msghandle(back_all_mh);
  
      num = 0;
      
      for(mu=0; mu < Nd; ++mu) {
	
	if(machine_size[mu] > 1) {
      
	  /* QMP_free_msghandle(forw_mh[num][0]); */
	  /* QMP_free_msghandle(forw_mh[num][1]); */
	  QMP_free_msgmem(forw_msg[num][0]);
	  QMP_free_msgmem(forw_msg[num][1]);

	  /* QMP_free_msghandle(back_mh[num][0]); */
	  /* QMP_free_msghandle(back_mh[num][1]); */
	  QMP_free_msgmem(back_msg[num][0]);
	  QMP_free_msgmem(back_msg[num][1]);
	  
	  num++;
	}
      }
    }
    
    total_comm = 0;
  }
}

/***************** end of initialization routine ***************************************/

/*include passing in shift table as argument later */
/* This macro is legacy ... I am adapting it not to bother
   with the 'smp' call but it occurs too often for me to strip it completely
   for now. Is that true? */

/* #define smpscaller2(a,bleah,spinfun2,chifun2,u3,cb2,volume2) \
    a.spinfun = spinfun2;\
    a.chifun = chifun2;\
    a.u = u3;  \
    a.cb = cb2; \
    printf(""); \
    smp_scall(bleah, (size_t)(volume2), sizeof(a), &a); */

/* Stripped out the sml_scall from this version -- non threaded */
/* Hence n = 0. => lo = 0, hi = volume2, */
#ifndef DSLASH_USE_QMT_THREADS
#define smpscaller2(a,bleah,spinfun2,chifun2,u3,cb2,volume2) \
    a.spinfun = spinfun2;\
    a.chifun = chifun2;\
    a.u = u3;  \
    a.cb = cb2; \
    (*bleah)(0, volume2, 0, &a);
#else
#define smpscaller2(a,bleah,spinfun2,chifun2,u3,cb2,volume2) \
    a.spinfun = spinfun2;\
    a.chifun = chifun2;\
    a.u = u3;  \
    a.cb = cb2; \
    qmt_call((qmt_userfunc_t)bleah, volume2, &a);
#endif
void sse_su3dslash_wilson(float *u, float *psi, float *res, int isign, int cb)
{
  int mu;
  Arg_s a;

  if (initP == 0) {
    QMP_error("sse_su3dslash_wilson not initialized");
    QMP_abort(1);
  }

  if(isign==1) 
  {
//    QMP_info("SSE dslash: isign=+1, calling decomp_plus");

    smpscaller2(a,decomp_plus,
		(MY_SPINOR*)psi,
		chi1,
		(MY_GAUGE_ARRAY)u,
		cb,
		subgrid_vol_cb);

    _prefetch_su3((u+0+2*(0+3*(0+3*(0+4*(icolor_start[cb]))))));
    _prefetch_single(chi2);

//    QMP_info("SSE dslash: start forw comm");

    if (total_comm > 0)
      if (QMP_start(forw_all_mh) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	QMP_abort(1);
      }

//    QMP_info("SSE dslash: call decomp_hvv_plus");

    /*other checkerboard's worth */
    smpscaller2(a,decomp_hvv_plus,
		(MY_SPINOR*)psi,
		chi2,
		(MY_GAUGE_ARRAY)u,
		cb,
		subgrid_vol_cb);
	
//    QMP_info("SSE dslash: wait on forw comm");

    if (total_comm > 0)
      if (QMP_wait(forw_all_mh) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	QMP_abort(1);
      }

    _prefetch_su3(u+0+2*(0+3*(0+3*(0+4*(icolor_start[1-cb])))));
    _prefetch_single(chi1);

//    QMP_info("SSE dslash: start on back comm");

    if (total_comm > 0)
      if (QMP_start(back_all_mh) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in backward direction");
	QMP_abort(1);
      }

//    QMP_info("SSE dslash: call mvv_recons_plus");

    /* current cb's u */
    /* for(mu=0; mu < Nd; ++mu)	
       fprintf(stderr,"node: %i recv:first value in tail1: %f   tail2: %f\n",
       PARSMP_physical_nodeid(), (*(TAIL1(chi1,mu)))[0][0][0], (*(TAIL2(chi1,mu))[0][0][0])); */
    smpscaller2(a,mvv_recons_plus,
		(MY_SPINOR*)res,
		chi1,
		(MY_GAUGE_ARRAY)u,
		1-cb,
		subgrid_vol_cb);

//    QMP_info("SSE dslash: wait on back comm");

    if (total_comm > 0)
      if (QMP_wait(back_all_mh) != QMP_SUCCESS)
      {
	QMP_error("wnxtsu3dslash: QMP_wait failed in backward direction");
	QMP_abort(1);
      }

//    QMP_info("SSE dslash: call recons_plus");

    smpscaller2(a,recons_plus,
		(MY_SPINOR*)res, 
		chi2,
		(MY_GAUGE_ARRAY)u,	
		1-cb,
		subgrid_vol_cb);
  }		

  if(isign==-1) 
  {
    smpscaller2(a,decomp_minus,
		(MY_SPINOR*)psi,
		chi1,
		(MY_GAUGE_ARRAY)u,
		cb,
		subgrid_vol_cb);

    _prefetch_su3((u+0+2*(0+3*(0+3*(0+4*(icolor_start[cb]))))));
    _prefetch_single(chi2);

    if (total_comm > 0)
      if (QMP_start(forw_all_mh) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	QMP_abort(1);
      }

    /*other checkerboard's worth */
    smpscaller2(a,decomp_hvv_minus,
		(MY_SPINOR*)psi,
		chi2,
		(MY_GAUGE_ARRAY)u,
		cb,
		subgrid_vol_cb);
    _prefetch_single(chi2);

    if (total_comm > 0)
      if (QMP_wait(forw_all_mh) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	QMP_abort(1);
      }

    _prefetch_su3(u+0+2*(0+3*(0+3*(0+4*(icolor_start[1-cb])))));
    _prefetch_single(chi1);

    if (total_comm > 0)
      if (QMP_start(back_all_mh) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in backward direction");
	QMP_abort(1);
      }

    /*current cb's u */
    smpscaller2(a,mvv_recons_minus,
		(MY_SPINOR*)res,
		chi1,
		(MY_GAUGE_ARRAY)u,
		1-cb,
		subgrid_vol_cb);
    
    _prefetch_single(chi2);

    if (total_comm > 0)
      if (QMP_wait(back_all_mh) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_wait failed in backward direction");
	QMP_abort(1);
      }

    smpscaller2(a,recons_minus,
		(MY_SPINOR*)res, 
		chi2,
		(MY_GAUGE_ARRAY)u,	
		1-cb,
		subgrid_vol_cb);
  }		
}



#ifdef __cplusplus
}
#endif

