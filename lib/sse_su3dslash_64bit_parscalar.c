/*******************************************************************************
 * $Id: sse_su3dslash_64bit_parscalar.c,v 1.2 2007-09-13 13:52:14 bjoo Exp $
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
/*	       /    mu			  mu */
/*	       --- */
/*	       mu=0 */

/*	             Nd-1 */
/*	             --- */
/*	             \    + */
/*                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu) */
/*	             /    mu			   mu */
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

#include <sse64.h>

extern void make_shift_tables(int *shift, int icolor_start[2], int bound[2][4][4]);

#define BASE 0x3f

#include <sse_align.h>

static int subgrid_vol = 0;
static int subgrid_vol_cb = 0;
static int initP=0;

/* now overlays for spinors as arrays or structs */

static int *newshift;
static int icolor_start[2];    /* starting site for each coloring (cb) */
static int icolor_end[2];      /* end site for each coloring (cb) */

#ifdef BACKWARD 
#undef BACKWARD
#endif

#define BACKWARD 0 

#ifdef FORWARD
#undef FORWARD
#endif
#define FORWARD 1


#define iup(mysite,mymu) newshift[mymu+4*(mysite+subgrid_vol*(1))]
#define idn(mysite,mymu) newshift[mymu+4*(mysite+subgrid_vol*(0))]
#define mvv_shift(mysite,mymu) newshift[mymu+4*(mysite+subgrid_vol*(2))]
#define rec_shift(mysite,mymu) newshift[mymu+4*(mysite+subgrid_vol*(3))]


#define a_chia(mymu,mysite) (chi+mysite+3*subgrid_vol_cb*mymu)
#define a_chib(mymu,mysite) (chi+mysite+3*subgrid_vol_cb*mymu)

#define _gauge_field0_0(mysite) gauge_field[mysite][0]
#define _gauge_field0_1(mysite) gauge_field[mysite][1]
#define _gauge_field0_2(mysite) gauge_field[mysite][2]
#define _gauge_field0_3(mysite) gauge_field[mysite][3]
#define _gauge_field_opp0_0(mysite) gauge_field_opp[mysite+1][0]
#define _gauge_field_opp0_1(mysite) gauge_field_opp[mysite+1][1]
#define _gauge_field_opp0_2(mysite) gauge_field_opp[mysite+1][2]
#define _gauge_field_opp0_3(mysite) gauge_field_opp[mysite+1][3]

/* now overlays for spinors as arrays or structs */
typedef double chi_double[2] __attribute__ ((aligned (16)));
typedef chi_double chi_three[3] __attribute__ ((aligned (16)));
typedef double u_mat_array[3][3][2]  ALIGN;  /* color color re/im */ 
typedef double spinor_array[4][3][2] ALIGN; /* Nspin4 color re/im */
typedef chi_three chi_array[2]    ALIGN; /*.. Nspin2 color re/im ::note:: Nspin2 has to be slowest varying */
typedef u_mat_array (*my_mat_array)[4] ALIGN;  


#define MY_SPINOR spinor_array
#define MY_SSE_VECTOR chi_array
#define MY_SSE_HALFVECT chi_three
#define MY_SSE_DOUBLE sse_double  
#define _c1__ [0] 
#define _c2__ [1] 
#define _c3__ [2]
#define _c4__ [3]
#define rs_c1__ rs[0]
#define rs_c2__ rs[1]
#define rs_c3__ rs[2]
#define rs_c4__ rs[3]


/* now overlays for gauge matrices as arrays or structs */
#define MY_GAUGE u_mat_array
#define MY_GAUGE_ARRAY my_mat_array

#define MY_SPINOR spinor_array
#define MY_SSE_VECTOR chi_array
#define _c1__ [0]
#define _c2__ [1]
#define _c3__ [2]
#define _c4__ [3]


/* now overlays for gauge matrices as arrays or structs */
#define MY_GAUGE u_mat_array



/* end overlays */


/* macros for spin basis note: assume first two rows are linearly independent except for gamma3 */

#define _sse_load_temp_1(chi) _sse_load((chi)_c1__)
#define _sse_load_temp_2(chi) _sse_load((chi)_c2__)
#define _sse_load_up_temp_1(chi) _sse_load_up((chi)_c1__)
#define _sse_load_up_temp_2(chi) _sse_load_up((chi)_c2__)

#define _sse_store_temp_1(chi) _sse_store((chi)_c1__)
#define _sse_store_temp_2(chi) _sse_store((chi)_c2__)
#define _sse_store_up_temp_1(chi) _sse_store_up((chi)_c1__)
#define _sse_store_up_temp_2(chi) _sse_store_up((chi)_c2__)


#define _sse_psum_set(psi) \
	  rs_c1__ = (psi)_c1__; \
	  rs_c2__ = (psi)_c2__; \
	  rs_c3__ = (psi)_c3__; \
	  rs_c4__ = (psi)_c4__


/* gamma 0 */

#define _sse_42_1_gamma0_minus(sp) \
      _sse_load((sp)_c1__); \
      _sse_load_up((sp)_c4__);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

#define _sse_24_1_gamma0_minus_set() \
      _sse_store_up(rs_c1__);\
      _sse_vector_i_mul_up();\
      _sse_store_up(rs_c4__) 

#define _sse_24_1_gamma0_minus_add() \
      _sse_load(rs_c1__);\
      _sse_vector_add();\
      _sse_store(rs_c1__);\
      _sse_load(rs_c4__);\
      _sse_vector_i_mul();\
      _sse_vector_add();\
      _sse_store(rs_c4__) 
	  
#define _sse_42_2_gamma0_minus(sp) \
      _sse_load((sp)_c2__);\
      _sse_load_up((sp)_c3__);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

#define _sse_24_2_gamma0_minus_set() \
	  _sse_store_up(rs_c2__);\
      _sse_vector_i_mul_up();\
      _sse_store_up(rs_c3__)

#define _sse_24_2_gamma0_minus_add() \
	  _sse_load(rs_c2__);\
      _sse_vector_add();\
      _sse_store(rs_c2__);\
      _sse_load(rs_c3__);\
      _sse_vector_i_mul();\
      _sse_vector_add();\
      _sse_store(rs_c3__)

#define _sse_42_1_gamma0_plus(sm) \
      _sse_load((sm)_c1__);\
      _sse_load_up((sm)_c4__);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_24_1_gamma0_plus_set() \
	  _sse_store_up(rs_c1__);\
      _sse_vector_i_mul_neg_up();\
      _sse_store_up(rs_c4__)

#define _sse_24_1_gamma0_plus_add() \
	  _sse_load(rs_c1__);\
      _sse_vector_add();\
      _sse_store(rs_c1__);\
      _sse_load(rs_c4__);\
      _sse_vector_i_mul();\
      _sse_vector_sub();\
      _sse_store(rs_c4__)

#define _sse_42_2_gamma0_plus(sm) \
	  _sse_load((sm)_c2__);\
      _sse_load_up((sm)_c3__);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_24_2_gamma0_plus_set() \
      _sse_store_up(rs_c2__);\
      _sse_vector_i_mul_neg_up();  \
      _sse_store_up(rs_c3__)

#define _sse_24_2_gamma0_plus_add() \
       _sse_load(rs_c2__);\
      _sse_vector_add();\
      _sse_store(rs_c2__);\
      _sse_load(rs_c3__);\
      _sse_vector_i_mul();  \
      _sse_vector_sub();\
      _sse_store(rs_c3__)



/* gamma 1 */


#define _sse_42_1_gamma1_minus(sp) \
      _sse_load((sp)_c1__);\
      _sse_load_up((sp)_c4__);\
      _sse_vector_add()

#define _sse_24_1_gamma1_minus() \
      _sse_load(rs_c1__);\
      _sse_vector_add();\
      _sse_store(rs_c1__);\
      _sse_load(rs_c4__);\
      _sse_vector_add();\
      _sse_store(rs_c4__)
	  
#define _sse_42_2_gamma1_minus(sp) \
      _sse_load((sp)_c2__);\
      _sse_load_up((sp)_c3__);\
      _sse_vector_sub()

#define _sse_24_2_gamma1_minus() \
	  _sse_load(rs_c2__);\
      _sse_vector_add();\
      _sse_store(rs_c2__);\
      _sse_load(rs_c3__);\
      _sse_vector_sub();\
      _sse_store(rs_c3__)

#define _sse_42_1_gamma1_plus(sm) \
      _sse_load((sm)_c1__);\
      _sse_load_up((sm)_c4__);\
      _sse_vector_sub()

#define _sse_24_1_gamma1_plus() \
      _sse_load(rs_c1__);\
      _sse_vector_add();\
      _sse_store(rs_c1__);\
      _sse_load(rs_c4__);\
      _sse_vector_sub();\
      _sse_store(rs_c4__)

#define _sse_42_2_gamma1_plus(sm) \
	  _sse_load((sm)_c2__);\
      _sse_load_up((sm)_c3__);\
      _sse_vector_add()

#define _sse_24_2_gamma1_plus() \
       _sse_load(rs_c2__);\
      _sse_vector_add();\
      _sse_store(rs_c2__);\
      _sse_load(rs_c3__);\
      _sse_vector_add();\
      _sse_store(rs_c3__)





/* gamma 2 */


#define _sse_42_1_gamma2_minus(sp) \
      _sse_load((sp)_c1__);\
      _sse_load_up((sp)_c3__);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

#define _sse_24_1_gamma2_minus() \
      _sse_load(rs_c1__);\
      _sse_vector_add();\
      _sse_store(rs_c1__);\
      _sse_load(rs_c3__);\
      _sse_vector_i_mul();   \
      _sse_vector_add();\
      _sse_store(rs_c3__)
	  
#define _sse_42_2_gamma2_minus(sp) \
      _sse_load((sp)_c2__);\
      _sse_load_up((sp)_c4__);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_24_2_gamma2_minus() \
	   _sse_load(rs_c2__);\
      _sse_vector_add();\
      _sse_store(rs_c2__);\
      _sse_load(rs_c4__);\
      _sse_vector_i_mul();   \
      _sse_vector_sub();\
      _sse_store(rs_c4__)

#define _sse_42_1_gamma2_plus(sm) \
      _sse_load((sm)_c1__);\
      _sse_load_up((sm)_c3__);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_24_1_gamma2_plus() \
      _sse_load(rs_c1__);\
      _sse_vector_add();\
      _sse_store(rs_c1__);\
      _sse_load(rs_c3__);\
      _sse_vector_i_mul();   \
      _sse_vector_sub();\
      _sse_store(rs_c3__);

#define _sse_42_2_gamma2_plus(sm) \
	  _sse_load((sm)_c2__);\
      _sse_load_up((sm)_c4__);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

#define _sse_24_2_gamma2_plus() \
      _sse_load(rs_c2__);\
      _sse_vector_add();\
      _sse_store(rs_c2__);\
      _sse_load(rs_c4__);\
      _sse_vector_i_mul();     \
      _sse_vector_add();\
      _sse_store(rs_c4__)





/* gamma 3 */
#define _sse_42_1_gamma3_minus(sp) \
	  _sse_load((sp)_c1__); \
	  _sse_load_up((sp)_c3__); \
      _sse_vector_sub()

#define _sse_24_1_gamma3_minus_set() \
	  _sse_load(rs_c1__);\
      _sse_vector_add();\
       _sse_store((*rn)_c1__);\
      _sse_load(rs_c3__);\
      _sse_vector_sub();\
       _sse_store((*rn)_c3__)

#define _sse_24_1_gamma3_minus_add() \
	  _sse_load(rs_c1__);\
      _sse_vector_add();\
       _sse_store(rs_c1__);\
      _sse_load(rs_c3__);\
      _sse_vector_sub();\
       _sse_store(rs_c3__)
	  
#define _sse_42_2_gamma3_minus(sp) \
      _sse_load((sp)_c2__);\
      _sse_load_up((sp)_c4__);\
      _sse_vector_sub()

#define _sse_24_2_gamma3_minus_set() \
	  _sse_load(rs_c2__);\
      _sse_vector_add();\
       _sse_store((*rn)_c2__);\
      _sse_load(rs_c4__);\
      _sse_vector_sub();\
      _sse_store((*rn)_c4__)

#define _sse_24_2_gamma3_minus_add() \
	  _sse_load(rs_c2__);\
      _sse_vector_add();\
       _sse_store(rs_c2__);\
      _sse_load(rs_c4__);\
      _sse_vector_sub();\
      _sse_store(rs_c4__)

#define _sse_42_1_gamma3_plus(sm) \
      _sse_load((sm)_c1__);\
      _sse_load_up((sm)_c3__);\
      _sse_vector_add()

#define _sse_24_1_gamma3_plus_set() \
	  _sse_load(rs_c1__);\
      _sse_vector_add();\
       _sse_store((*rn)_c1__);\
      _sse_load(rs_c3__);\
      _sse_vector_add();\
      _sse_store((*rn)_c3__)

#define _sse_24_1_gamma3_plus_add() \
	  _sse_load(rs_c1__);\
      _sse_vector_add();\
       _sse_store(rs_c1__);\
      _sse_load(rs_c3__);\
      _sse_vector_add();\
      _sse_store(rs_c3__)

#define _sse_42_2_gamma3_plus(sm) \
	  _sse_load((sm)_c2__);\
      _sse_load_up((sm)_c4__);\
      _sse_vector_add()

#define _sse_24_2_gamma3_plus_set() \
       _sse_load(rs_c2__);\
      _sse_vector_add();\
       _sse_store((*rn)_c2__);\
      _sse_load(rs_c4__);\
      _sse_vector_add();\
       _sse_store((*rn)_c4__)

#define _sse_24_2_gamma3_plus_add() \
       _sse_load(rs_c2__);\
       _sse_vector_add();\
       _sse_store(rs_c2__);\
       _sse_load(rs_c4__);\
       _sse_vector_add();\
       _sse_store(rs_c4__);\



static int init=0;
static sse_double fact1,fact2;
static MY_SPINOR rs __attribute__ ((aligned (16)));


void D_psi_fun(size_t lo,size_t hi, int id, const void *ptr);

void D_psi_fun_plus(size_t lo,size_t hi, int id, const void *ptr);


void D_psi_fun_minus(size_t lo,size_t hi, int id, const void *ptr);


typedef struct
{
   sse_double c1,c2,c3;
} sse_half_vector __attribute__ ((aligned (16)));

typedef struct
{
   sse_half_vector c1,c2;
} sse_vector __attribute__ ((aligned (16)));


typedef struct {
  double m0;  /* the input array                          */
  int  k;  /* the phase factors                        */
  int l;
  MY_SPINOR* psi;  /* for Luescher */
  MY_SPINOR* res;
  MY_SPINOR* spinfun;
  MY_SSE_VECTOR *chifun;  /*note this must be allocated ..direction varying more slowly than vol */
  MY_GAUGE       (*u)[4];
  MY_GAUGE       (*u2)[4];
  int cb;
  /* the output array                         */
} Arg_s;





#define DECLARE_COMMON_STUFF \
	  int ix,iy,iz;\
     const Arg_s *a =(Arg_s *)ptr; \
   static int ix1,iy1,iy2,iz1;\
   const int low = lo; \
   const int high = hi;\
    MY_GAUGE (*gauge_field)[4] = a->u;\
    MY_SPINOR *psi = a->spinfun;\
	MY_SSE_VECTOR *chi = a->chifun;\
   MY_GAUGE *up ALIGN,*um ALIGN;\
   MY_SSE_VECTOR *s2 ALIGN,*s3 ALIGN, *s4 ALIGN; \
   MY_SPINOR *s ALIGN, *s1 ALIGN, *sp ALIGN,*sm ALIGN,rs ALIGN,*rn ALIGN;\
   const int cb = a->cb; \
   sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; \
   sse_double _minus_one ALIGN = {-1.0, -1.0}
 
/* watch out for sm and sp */


/* this routine is similar to wnxtsu3dslash, except instead of handling the second site's worth in the same loop, the second
spin component's worth must be handled seperately */
void decomp_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  DECLARE_COMMON_STUFF;

  for (ix1=lo;ix1<hi;ix1+=1) 
  {
#ifdef PREFDIST
#undef PREFDIST
#endif
#define PREFDIST 4
    sp=&psi[ix1];
    s1=&psi[ix1+PREFDIST];
    _prefetch_spinor(s1);
      
/******************************* direction +0 *********************************/	   
    s3 = a_chia(0,idn(ix1,0));
    /*spin decomposition of first component of halfspinor */
    _sse_42_1_gamma0_minus(*sp);
    _sse_store_temp_1(*s3);
    /*spin decomposition of first component of halfspinor */
    _sse_42_2_gamma0_minus(*sp);
    _sse_store_temp_2(*s3);

/******************************* direction +1 *********************************/
    s3 = a_chia(1,idn(ix1,1));
	   
    _sse_42_1_gamma1_minus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma1_minus(*sp);
    _sse_store_temp_2(*s3);
   
/******************************* direction +2 *********************************/
    s3 = a_chia(2,idn(ix1,2));
	   
    _sse_42_1_gamma2_minus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma2_minus(*sp);
    _sse_store_temp_2(*s3);

/******************************* direction +3 *********************************/
    s3 = a_chia(3,idn(ix1,3));
	   
    _sse_42_1_gamma3_minus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma3_minus(*sp);
    _sse_store_temp_2(*s3);
  }
}


void decomp_hvv_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  DECLARE_COMMON_STUFF;
   
  um=&_gauge_field0_0(lo);
  _prefetch_su3(um+1);
  for (ix1=lo;ix1<hi;ix1++) 
  {
    sm=&psi[ix1];
     
/******************************* direction -0 *********************************/
    s3 = a_chia(0,iup(ix1,0));
	   
    _sse_42_1_gamma0_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);

    _sse_42_2_gamma0_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);


/******************************* direction -1 *********************************/
    um++;
    _prefetch_su3(um+1);
    /* in decomp_hvv */
    s3 = a_chia(1,iup(ix1,1));
    
    _sse_42_1_gamma1_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);
    
    _sse_42_2_gamma1_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);


/******************************* direction -2 *********************************/

    um++;
    _prefetch_su3(um+1);

    s3 = a_chia(2,iup(ix1,2));
    
    _sse_42_1_gamma2_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);

    _sse_42_2_gamma2_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);

/******************************* direction -3 *********************************/

    um++;
    _prefetch_su3(um+1);

    s1=&psi[ix1+1];
    _prefetch_spinor(s1);

    s3 = a_chia(3,iup(ix1,3));
	   

    _sse_42_1_gamma3_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);

    if(ix1 < hi) _prefetch_su3(um+1);
    else _prefetch_su3(&_gauge_field0_0(lo));

    _sse_42_2_gamma3_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);

    um++;
       
	   
  }
}


void mvv_recons_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  DECLARE_COMMON_STUFF;
   	
   
 
  up=&_gauge_field0_0(lo);
  _prefetch_su3(up+1);
  s3 = a_chia(0,mvv_shift(lo,0));


  for (ix1=lo;ix1<hi;ix1++) 
  {
    rn=&psi[ix1];
	 
      
/******************************* direction +0 *********************************/	

    s4 = a_chia(1,mvv_shift(ix1,1));
    _prefetch_spinor(s4);



    _sse_load_temp_1(*s3);
    _sse_su3_multiply(*up);
    _sse_24_1_gamma0_minus_set();
	  

    _sse_load_temp_2(*s3);
    _sse_su3_multiply(*up);
    _sse_24_2_gamma0_minus_set();

    up++;
	   
/******************************* direction +1 *********************************/

     
	   

	 

    s3 = a_chia(2,mvv_shift(ix1,2));
    _prefetch_spinor(s3);
 
	   

    _sse_load_temp_1(*s4);
    _sse_su3_multiply(*up);
    _sse_24_1_gamma1_minus();
	  

    _sse_load_temp_2(*s4);
    _sse_su3_multiply(*up);
    _sse_24_2_gamma1_minus();

    up++;

/******************************* direction +2 *********************************/

	   

    _prefetch_spinor(rn);


	 

    s4 = a_chia(3,mvv_shift(ix1,3));
    _prefetch_spinor(s4);

    _sse_load_temp_1(*s3);
    _sse_su3_multiply(*up);
    _sse_24_1_gamma2_minus();
	  

    _sse_load_temp_2(*s3);
    _sse_su3_multiply(*up);
    _sse_24_2_gamma2_minus();





    up++;
    /******************************* direction +3 *********************************/     
	 

    s3 = a_chia(0,mvv_shift(ix1+1,0));
    _prefetch_spinor(s3);



    _sse_load_temp_1(*s4);
    _sse_su3_multiply(*up);
    _sse_24_1_gamma3_minus_set();
	  
    _prefetch_su3(up+1);
	  

    _sse_load_temp_2(*s4);
    _sse_su3_multiply(*up);
    _sse_24_2_gamma3_minus_set();

    up++;
  }
}

/*non-optimized recons...has too many loads and stores, but is still basis independent with the macros */
/*I'm hoping since the stuff should be in cache that the superfluous loads and stores won't be that bad */
void recons_plus_old(size_t lo, size_t hi, int id, const void *ptr)
{
  int ix,iy,iz;
  const Arg_s *a =(Arg_s *)ptr; 
	 
  static int ix1,iy1,iy2,iz1;
  const int low = lo;
  const int high = hi;
  MY_GAUGE (*gauge_field)[4] ALIGN = a->u;
  MY_SPINOR *psi = a->spinfun;
  MY_SSE_VECTOR *chi = a->chifun;
  MY_GAUGE *up,*um;
  MY_SSE_VECTOR *s2,*s3 ALIGN, *s4; 
  MY_SPINOR *s, *s1, *sp,*sm,rs ALIGN,*rn ALIGN;
  const int cb = a->cb;
  /* 	printf("\n &a:%x",&a);*/
 
  s3 = a_chia(0,rec_shift(lo,0));

#ifdef PREFDIST
#undef PREFDIST
#endif
#define PREFDIST 1
  /*printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
    printf("\n &psi[0]:%x",(&psi[0]));*/
  for (ix1=lo;ix1<hi;ix1++) 
  {
    rn=&psi[ix1];
	 
    s1=&psi[ix1+1];
    _prefetch_spinor(s1);

#define rs (*rn)
	   
    s4 = a_chia(1,rec_shift(ix1,1));
    s2 = a_chia(1,rec_shift(ix1+PREFDIST,1));
    _prefetch_single(s2);
     

    _sse_load_up_temp_1(*s3);
	   
    _sse_24_1_gamma0_plus_add();
	  

    _sse_load_up_temp_2(*s3);
	   
    _sse_24_2_gamma0_plus_add();

	   

    s3 = a_chia(2,rec_shift(ix1,2));
    s2 = a_chia(2,rec_shift(ix1+PREFDIST,2));
    _prefetch_single(s2);
 

    _sse_load_up_temp_1(*s4);
	   
    _sse_24_1_gamma1_plus();
	  

    _sse_load_up_temp_2(*s4);
	  
    _sse_24_2_gamma1_plus();

	   

    s4 = a_chia(3,rec_shift(ix1,3));
    s2 = a_chia(3,rec_shift(ix1+PREFDIST,3));
    _prefetch_single(s2);
	   

    _sse_load_up_temp_1(*s3);
	  
    _sse_24_1_gamma2_plus();
	  

    _sse_load_up_temp_2(*s3);
	   
    _sse_24_2_gamma2_plus();

	   
       	
    s3 = a_chia(0,rec_shift(ix1+1,0));
    s2 = a_chia(0,rec_shift(ix1+1+PREFDIST,0));
    _prefetch_single(s2);



    _sse_load_up_temp_1(*s4);
	   
    _sse_24_1_gamma3_plus_set();
	  


    _sse_load_up_temp_2(*s4);
	   
    _sse_24_2_gamma3_plus_set();

#undef rs
	   
    /*if(ix1==0) printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);*/
    /* could be problem because we write back to spinor field indirectly */
	   
  }
  /*printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
    printf("\n &psi[0]:%x",(&psi[0]));*/
 
}


/*optimized for SZIN spin basis */
void recons_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  int ix,iy,iz;
  const  Arg_s *a =(Arg_s *)ptr; 
	 
  static int ix1,iy1,iy2,iz1;
  const int low = lo;
  const int high = hi;
  MY_GAUGE (*gauge_field)[4] ALIGN = a->u;
  MY_SPINOR *psi = a->spinfun;
  MY_SSE_VECTOR *chi = a->chifun;
  MY_GAUGE *up,*um;
  MY_SSE_VECTOR *s2,*s3 ALIGN, *s4, *temp0,*temp1,*temp2,*temp3; 
  MY_SPINOR *s, *s1, *sp,*sm,rs ALIGN,*rn ALIGN;
  const int cb = a->cb;
  sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; 
  sse_double _minus_one ALIGN = {-1.0, -1.0};
  /* 	printf("\n &a:%x",&a);*/
 
  
#ifdef PREFDIST
#undef PREFDIST
#endif
#define PREFDIST 3
  /*printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
	printf("\n &psi[0]:%x",(&psi[0]));*/
  for (ix1=lo;ix1<hi;ix1++) 
  {
    rn=&psi[ix1];
	 
    s1=&psi[ix1+PREFDIST];
    _prefetch_nta_spinor(s1);

	  

    /* first spin component of result */
    temp0 = a_chia(0,rec_shift(ix1,0));
	  
	    

    _sse_load_up_temp_1(*(temp0));
    _sse_load((*rn)_c1__);
	   
    _sse_vector_add();

    temp1 = a_chia(1,rec_shift(ix1,1));

    _sse_load_up_temp_1(*(temp1));
    _sse_vector_add();

    temp2 = a_chia(2,rec_shift(ix1,2));

    _sse_load_up_temp_1(*(temp2));
    _sse_vector_add();

    temp3 = a_chia(3,rec_shift(ix1,3));

    _sse_load_up_temp_1(*(temp3));
    _sse_vector_add();
    _sse_store((*rn)_c1__);

/* second spin component of result */
    s2= a_chia(0,rec_shift(ix1+PREFDIST,0));
    _prefetch_single(s2);
     

    _sse_load_up_temp_2(*(temp0));
    _sse_load((*rn)_c2__);
	   
    _sse_vector_add();

    s2 = a_chia(1,rec_shift(ix1+PREFDIST,1));
    _prefetch_single(s2);

    _sse_load_up_temp_2(*(temp1));
    _sse_vector_add();

    s2 = a_chia(2,rec_shift(ix1+PREFDIST,2));
    _prefetch_single(s2);

    _sse_load_up_temp_2(*(temp2));
    _sse_vector_add();

    s2 = a_chia(3,rec_shift(ix1+PREFDIST,3));
    _prefetch_single(s2);

    _sse_load_up_temp_2(*(temp3));
    _sse_vector_add();
    _sse_store((*rn)_c2__);



		

    /* third spin component, here it gets tricky */
	  
	   
	   
     

    _sse_load_up_temp_2(*(temp0));
    _sse_load((*rn)_c3__);
	   
    _sse_vector_i_mul_up();
    _sse_vector_sub();

	  
	  

    _sse_load_up_temp_2(*(temp1));
		

    _sse_vector_add();

      
      

    _sse_load_up_temp_1(*(temp2));

    _sse_vector_i_mul();
    _sse_vector_sub();

    

    _sse_load_up_temp_1(*(temp3));

    _sse_vector_add();

    _sse_store((*rn)_c3__);

/* fourth spin component, again it gets tricky */
	  
	   
	   
     

    _sse_load_up_temp_1(*(temp0));
    _sse_load((*rn)_c4__);
	   
    _sse_vector_i_mul_up();
    _sse_vector_sub();

	  
	  

    _sse_load_up_temp_1(*(temp1));
		

    _sse_vector_sub();

      
      

    _sse_load_up_temp_2(*(temp2));

    _sse_vector_i_mul_up();
    _sse_vector_add();

    

    _sse_load_up_temp_2(*(temp3));

    _sse_vector_add();

    _sse_store((*rn)_c4__);

    /* end of loop */


    /*if(ix1==0) printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);*/
    /* could be problem because we write back to spinor field indirectly */
	   
  }
  /*printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
    printf("\n &psi[0]:%x",(&psi[0]));*/
 
}


/************now for isign = -1  **********************/


void decomp_minus(size_t lo, size_t hi, int id, const void *ptr)
{
  DECLARE_COMMON_STUFF;
   
    
 

  for (ix1=lo;ix1<hi;ix1+=1) 
  {

#ifdef PREFDIST
#undef PREFDIST
#endif
#define PREFDIST 4
    sp=&psi[ix1];
    s1=&psi[ix1+PREFDIST];
    _prefetch_spinor(s1);
  
	  
	   
    s3 = a_chia(0,idn(ix1,0));
	   
    _sse_42_1_gamma0_plus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma0_plus(*sp);
    _sse_store_temp_2(*s3);

    /*printf("chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s3)[0][0][0],(*s3)[1][0][0]);
 
      printf("chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s3)[0][2][0],(*s3)[1][2][0]); */

    s3 = a_chia(1,idn(ix1,1));
	   
    _sse_42_1_gamma1_plus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma1_plus(*sp);
    _sse_store_temp_2(*s3);
   

    /*printf("chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s3)[0][0][0],(*s3)[1][0][0]);
 
      printf("chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s3)[0][2][0],(*s3)[1][2][0]); */


    s3 = a_chia(2,idn(ix1,2));
	   
    _sse_42_1_gamma2_plus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma2_plus(*sp);
    _sse_store_temp_2(*s3);


    /*printf("chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s3)[0][0][0],(*s3)[1][0][0]);
 
      printf("chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s3)[0][2][0],(*s3)[1][2][0]); */


    s3 = a_chia(3,idn(ix1,3));
	   
    _sse_42_1_gamma3_plus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma3_plus(*sp);
    _sse_store_temp_2(*s3);

	   
    /* printf("s chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s3)[0][0][0],(*s3)[1][0][0]);
 
       printf("s chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s3)[0][2][0],(*s3)[1][2][0]); */

     
  }
}


void decomp_hvv_minus(size_t lo, size_t hi, int id, const void *ptr)
{
  DECLARE_COMMON_STUFF;
   
 
  um=&_gauge_field0_0(lo);
  _prefetch_su3(um+1);
  for (ix1=lo;ix1<hi;ix1++) 
  {
    sm=&psi[ix1];
	 
     

       	   
    s3 = a_chia(0,iup(ix1,0));
	   
    _sse_42_1_gamma0_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);

    _sse_42_2_gamma0_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);

    um++;
    _prefetch_su3(um+1);
    /* in decomp_hvv */
    s3 = a_chia(1,iup(ix1,1));
	   

    _sse_42_1_gamma1_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);

    _sse_42_2_gamma1_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);

    um++;
    _prefetch_su3(um+1);

    s3 = a_chia(2,iup(ix1,2));
	   

    _sse_42_1_gamma2_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);

    _sse_42_2_gamma2_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);

    um++;
    _prefetch_su3(um+1);

    s1=&psi[ix1+1];
    _prefetch_spinor(s1);

    s3 = a_chia(3,iup(ix1,3));
	   

    _sse_42_1_gamma3_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);

    if(ix1 < hi) _prefetch_su3(um+1);
    else _prefetch_su3(&_gauge_field0_0(lo));

    _sse_42_2_gamma3_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);

    um++;
       
	   
  }
}


void mvv_recons_minus(size_t lo, size_t hi, int id, const void *ptr)
{
  DECLARE_COMMON_STUFF;
   	
   
 
  up=&_gauge_field0_0(lo);
  _prefetch_su3(up+1);
  s3 = a_chia(0,mvv_shift(lo,0));


  for (ix1=lo;ix1<hi;ix1++) 
  {
    rn=&psi[ix1];
	 
      
    /* printf("g chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s3)[0][0][0],(*s3)[1][0][0]);
 
       printf("g chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s3)[0][2][0],(*s3)[1][2][0]); 
	   
       printf("u(ix1)[0][0][0]=%f \n", (*up)[0][0][0]);*/

    s4 = a_chia(1,mvv_shift(ix1,1));
    _prefetch_spinor(s4);



    _sse_load_temp_1(*s3);
    _sse_su3_multiply(*up);
    _sse_24_1_gamma0_plus_set();
	  

    _sse_load_temp_2(*s3);
    _sse_su3_multiply(*up);
    _sse_24_2_gamma0_plus_set();

    /* printf("rs(ix1)[0]=%f,rs(ix1)[1] [0][0]=%f\n",(rs_c1__)[0][0],(rs_c2__)[0][0]);
 
       printf("rs(ix1)[2]=%f,rs(ix1)[3] [0][0]=%f\n",(rs_c3__)[0][0],(rs_c4__)[0][0]);*/

    up++;


    /* printf("chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s4)[0][0][0],(*s4)[1][0][0]);
 
       printf("chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s4)[0][2][0],(*s4)[1][2][0]); */

    s3 = a_chia(2,mvv_shift(ix1,2));
    _prefetch_spinor(s3);
 
	   

    _sse_load_temp_1(*s4);
    _sse_su3_multiply(*up);
    _sse_24_1_gamma1_plus();
	  

    _sse_load_temp_2(*s4);
    _sse_su3_multiply(*up);
    _sse_24_2_gamma1_plus();

    up++;

    _prefetch_spinor(rn);


    /*printf("chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s3)[0][0][0],(*s3)[1][0][0]);
 
      printf("chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s3)[0][2][0],(*s3)[1][2][0]); */

    s4 = a_chia(3,mvv_shift(ix1,3));
    _prefetch_spinor(s4);

    _sse_load_temp_1(*s3);
    _sse_su3_multiply(*up);
    _sse_24_1_gamma2_plus();
	  

    _sse_load_temp_2(*s3);
    _sse_su3_multiply(*up);
    _sse_24_2_gamma2_plus();

    up++;
      
    /* printf("chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s4)[0][0][0],(*s4)[1][0][0]);
 
       printf("chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s4)[0][2][0],(*s4)[1][2][0]);  */

    s3 = a_chia(0,mvv_shift(ix1+1,0));
    _prefetch_spinor(s3);



    _sse_load_temp_1(*s4);
    _sse_su3_multiply(*up);
    _sse_24_1_gamma3_plus_set();
	  
    _prefetch_su3(up+1);
    /*else _prefetch_su3(&_gauge_field0_0(lo));*/

    _sse_load_temp_2(*s4);
    _sse_su3_multiply(*up);
    _sse_24_2_gamma3_plus_set();


    /*printf("psum(ix1)[0]=%f,psum(ix1)[1] [0][0]=%f\n",(*rn)[0][0][0],(*rn)[1][0][0]);
 
      printf("psum(ix1)[2]=%f,psum(ix1)[3] [0][0]=%f\n",(*rn)[2][0][0],(*rn)[3][0][0]);  */

    up++;

     
     
	   
  }
  /* printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
     printf("\n &psi[0]:%x",(&psi[0]));*/
}

/*non-optimized recons...has too many loads and stores, but is still basis independent with the macros */
/*I'm hoping since the stuff should be in lvl 1 cache that the superfluous loads and stores will hit in lvl1 cache */
void recons_minus_old(size_t lo, size_t hi, int id, const void *ptr)
{
  int ix,iy,iz;
  const  Arg_s *a =(Arg_s *)ptr; 
	 
  static int ix1,iy1,iy2,iz1;
  const int low = lo;
  const int high = hi;
  MY_GAUGE (*gauge_field)[4] ALIGN = a->u;
  MY_SPINOR *psi = a->spinfun;
  MY_SSE_VECTOR *chi = a->chifun;
  MY_GAUGE *up,*um;
  MY_SSE_VECTOR *s2,*s3 ALIGN, *s4; 
  MY_SPINOR *s, *s1, *sp,*sm,rs ALIGN,*rn ALIGN;
  const int cb = a->cb;
  /* 	printf("\n &a:%x",&a);*/
 
  s3 = a_chia(0,rec_shift(lo,0));

#ifdef PREFDIST
#undef PREFDIST
#endif

#define PREFDIST 1
  /*printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
    printf("\n &psi[0]:%x",(&psi[0]));*/
  for (ix1=lo;ix1<hi;ix1++) 
  {
    rn=&psi[ix1];
	 
    s1=&psi[ix1+1];
    _prefetch_spinor(s1);

#define rs (*rn)
	   
    s4 = a_chia(1,rec_shift(ix1,1));
    s2 = a_chia(1,rec_shift(ix1+PREFDIST,1));
    _prefetch_single(s2);
     

    _sse_load_up_temp_1(*s3);
	   
    _sse_24_1_gamma0_minus_add();
	  

    _sse_load_up_temp_2(*s3);
	   
    _sse_24_2_gamma0_minus_add();

	   

    s3 = a_chia(2,rec_shift(ix1,2));
    s2 = a_chia(2,rec_shift(ix1+PREFDIST,2));
    _prefetch_single(s2);
 

    _sse_load_up_temp_1(*s4);
	   
    _sse_24_1_gamma1_minus();
	  

    _sse_load_up_temp_2(*s4);
	  
    _sse_24_2_gamma1_minus();

	   

    s4 = a_chia(3,rec_shift(ix1,3));
    s2 = a_chia(3,rec_shift(ix1+PREFDIST,3));
    _prefetch_single(s2);
	   

    _sse_load_up_temp_1(*s3);
	  
    _sse_24_1_gamma2_minus();
	  

    _sse_load_up_temp_2(*s3);
	   
    _sse_24_2_gamma2_minus();

	   
       	
    s3 = a_chia(0,rec_shift(ix1+1,0));
    s2 = a_chia(0,rec_shift(ix1+1+PREFDIST,0));
    _prefetch_single(s2);



    _sse_load_up_temp_1(*s4);
	   
    _sse_24_1_gamma3_minus_set();
	  


    _sse_load_up_temp_2(*s4);
	   
    _sse_24_2_gamma3_minus_set();

#undef rs
	   
    /*if(ix1==0) printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);*/
    /* could be problem because we write back to spinor field indirectly */
	   
  }
  /*printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
    printf("\n &psi[0]:%x",(&psi[0]));*/
 
}



void recons_minus(size_t lo, size_t hi, int id, const void *ptr)
{
  int ix,iy,iz;
  const  Arg_s *a =(Arg_s *) ptr; 
	 
  static int ix1,iy1,iy2,iz1;
  const int low = lo;
  const int high = hi;
  MY_GAUGE (*gauge_field)[4] ALIGN = a->u;
  MY_SPINOR *psi = a->spinfun;
  MY_SSE_VECTOR *chi = a->chifun;
  MY_GAUGE *up,*um;
  MY_SSE_VECTOR *s2,*s3 ALIGN, *s4, *temp0,*temp1,*temp2,*temp3; 
  MY_SPINOR *s, *s1, *sp,*sm,rs ALIGN,*rn ALIGN;
  const int cb = a->cb;
  sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; 
  sse_double _minus_one ALIGN = {-1.0, -1.0};
  /* 	printf("\n &a:%x",&a);*/
 
  
#ifdef PREFDIST
#undef PREFDIST
#endif
#define PREFDIST 3
  /*printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
    printf("\n &psi[0]:%x",(&psi[0]));*/
  for (ix1=lo;ix1<hi;ix1++) 
  {
    rn=&psi[ix1];
	 
    s1=&psi[ix1+PREFDIST];
    _prefetch_nta_spinor(s1);

	  

    /* first spin component of result */
    temp0 = a_chia(0,rec_shift(ix1,0));
	  
	    

    _sse_load_up_temp_1(*(temp0));
    _sse_load((*rn)_c1__);
	   
    _sse_vector_add();

    temp1 = a_chia(1,rec_shift(ix1,1));

    _sse_load_up_temp_1(*(temp1));
    _sse_vector_add();

    temp2 = a_chia(2,rec_shift(ix1,2));

    _sse_load_up_temp_1(*(temp2));
    _sse_vector_add();

    temp3 = a_chia(3,rec_shift(ix1,3));

    _sse_load_up_temp_1(*(temp3));
    _sse_vector_add();
    _sse_store((*rn)_c1__);

/* second spin component of result */
    s2= a_chia(0,rec_shift(ix1+PREFDIST,0));
    _prefetch_single(s2);
     

    _sse_load_up_temp_2(*(temp0));
    _sse_load((*rn)_c2__);
	   
    _sse_vector_add();

    s2 = a_chia(1,rec_shift(ix1+PREFDIST,1));
    _prefetch_single(s2);

    _sse_load_up_temp_2(*(temp1));
    _sse_vector_add();

    s2 = a_chia(2,rec_shift(ix1+PREFDIST,2));
    _prefetch_single(s2);

    _sse_load_up_temp_2(*(temp2));
    _sse_vector_add();

    s2 = a_chia(3,rec_shift(ix1+PREFDIST,3));
    _prefetch_single(s2);

    _sse_load_up_temp_2(*(temp3));
    _sse_vector_add();
    _sse_store((*rn)_c2__);



		

    /* third spin component, here it gets tricky */
	  
	   
	   
     

    _sse_load_up_temp_2(*(temp0));
    _sse_load((*rn)_c3__);
	   
    _sse_vector_i_mul_up();
    _sse_vector_add();

	  
	  

    _sse_load_up_temp_2(*(temp1));
		

    _sse_vector_sub();

      
      

    _sse_load_up_temp_1(*(temp2));

    _sse_vector_i_mul();
    _sse_vector_add();

    

    _sse_load_up_temp_1(*(temp3));

    _sse_vector_sub();

    _sse_store((*rn)_c3__);

/* fourth spin component, again it gets tricky */
	  
	   
	   
     

    _sse_load_up_temp_1(*(temp0));
    _sse_load((*rn)_c4__);
	   
    _sse_vector_i_mul_up();
    _sse_vector_add();

	  
	  

    _sse_load_up_temp_1(*(temp1));
		

    _sse_vector_add();

      
      

    _sse_load_up_temp_2(*(temp2));

    _sse_vector_i_mul_up();
    _sse_vector_sub();

    

    _sse_load_up_temp_2(*(temp3));

    _sse_vector_sub();

    _sse_store((*rn)_c4__);

    /* end of loop */


    /*if(ix1==0) printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);*/
    /* could be problem because we write back to spinor field indirectly */
	   
  }
  /*printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
    printf("\n &psi[0]:%x",(&psi[0]));*/
 
}


#define Nd 4
#define Nc 3
#define Ns 4
#define Ns2 2

static MY_SSE_VECTOR* chi1;
static MY_SSE_VECTOR* chi2;
#define TAIL1(chi,mymu) (chi+subgrid_vol_cb*(1+3*mymu))
#define TAIL2(chi,mymu) (chi+0+subgrid_vol_cb*(2+3*mymu))

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
static int bound[Nd];


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
  nsize = 2*Nc*Ns2*sizeof(double)*subgrid_vol_cb*3*Nd;  /* Note 3x4 half-ferm temps */
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
/*
#define smpscaller2(a,bleah,spinfun2,chifun2,u3,cb2,volume2) \
    a.spinfun = spinfun2;\
	a.chifun = chifun2;\
	a.u = u3;  \
	a.cb = cb2; \
	printf(""); \
	smp_scall(bleah, (size_t)(volume2), sizeof(a), &a);
*/
#define smpscaller2(a,bleah,spinfun2,chifun2,u3,cb2) \
    a.spinfun = spinfun2;\
    a.chifun = chifun2;\
    a.u = u3;  \
    a.cb = cb2; \
    (*bleah)(icolor_start[cb2], icolor_end[cb2], 0, &a);


void sse_su3dslash_wilson(SSEREAL *u, SSEREAL *psi, SSEREAL *res, int isign, int cb)
{
  int count, mu;
  Arg_s a;

  if (initP == 0) {
    QMP_error("sse_su3dslash_wilson not initialised");
    QMP_abort(1);
  }

  if(isign==1) 
  {
    smpscaller2(a,decomp_plus,
		(MY_SPINOR*)psi,
		chi1,
		(MY_GAUGE_ARRAY)u,
		cb);

    _prefetch_su3((u+0+2*(0+3*(0+3*(0+4*(icolor_start[cb]))))));
    _prefetch_single(chi2);


    if (total_comm > 0)
      if (QMP_start(forw_all_mh) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	QMP_abort(1);
      }

    /*other checkerboard's worth */
    smpscaller2(a,decomp_hvv_plus,
		(MY_SPINOR*)psi,
		chi2,
		(MY_GAUGE_ARRAY)u,
		cb);
	
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

    smpscaller2(a,mvv_recons_plus,
		(MY_SPINOR*)res,
		chi1,
		(MY_GAUGE_ARRAY)u,
		1-cb);
	
    if (total_comm > 0)
      if (QMP_wait(back_all_mh) != QMP_SUCCESS)
      {
	QMP_error("wnxtsu3dslash: QMP_wait failed in backward direction");
	QMP_abort(1);
      }


    smpscaller2(a,recons_plus,
		(MY_SPINOR*)res, 
		chi2,
		(MY_GAUGE_ARRAY)u,	
		1-cb);
  }		

  if(isign==-1) 
  {
    smpscaller2(a,decomp_minus,
		(MY_SPINOR*)psi,
		chi1,
		(MY_GAUGE_ARRAY)u,
		cb);

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
		cb);
    
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
	QMP_error("wnxtsu3dslash: QMP_start failed in backward direction");
	QMP_abort(1);
      }

    /*current cb's u */
    smpscaller2(a,mvv_recons_minus,
		(MY_SPINOR*)res,
		chi1,
		(MY_GAUGE_ARRAY)u,
		1-cb);

    _prefetch_single(chi2);

    if (total_comm > 0)
      if (QMP_wait(back_all_mh) != QMP_SUCCESS)
      {
	QMP_error("wnxtsu3dslash: QMP_wait failed in backward direction"); 
	QMP_abort(1);
      }

    smpscaller2(a,recons_minus,
		(MY_SPINOR*)res, 
		chi2,
		(MY_GAUGE_ARRAY)u,	
		1-cb);
  }		
}


#ifdef __cplusplus
}
#endif
