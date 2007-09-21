/*******************************************************************************
 * $Id: sse_su3dslash_64bit_parscalar.c,v 1.3 2007-09-21 17:19:48 bjoo Exp $
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


#include <sse_config.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <qmp.h>

#include <sse64.h>
#include <sse_align.h>
#include <dispatch_parscalar.h>
#include <shift_tables_parscalar.h>

static int initP=0;

/* now overlays for spinors as arrays or structs */

static int icolor_start[2];    /* starting site for each coloring (cb) */
static int icolor_end[2];      /* end site for each coloring (cb) */

#define _gauge_field0_0(mysite) gauge_field[mysite][0]
#define _gauge_field0_1(mysite) gauge_field[mysite][1]
#define _gauge_field0_2(mysite) gauge_field[mysite][2]
#define _gauge_field0_3(mysite) gauge_field[mysite][3]
#define _gauge_field_opp0_0(mysite) gauge_field_opp[mysite+1][0]
#define _gauge_field_opp0_1(mysite) gauge_field_opp[mysite+1][1]
#define _gauge_field_opp0_2(mysite) gauge_field_opp[mysite+1][2]
#define _gauge_field_opp0_3(mysite) gauge_field_opp[mysite+1][3]

#if 0
/* now overlays for spinors as arrays or structs */
typedef double chi_double[2] __attribute__ ((aligned (16)));
typedef chi_double chi_three[3] __attribute__ ((aligned (16)));
typedef double u_mat_array[3][3][2]  ALIGN;  /* color color re/im */ 
typedef double spinor_array[4][3][2] ALIGN; /* Nspin4 color re/im */
typedef chi_three halfspinor_array[2]    ALIGN; /*.. Nspin2 color re/im ::note:: Nspin2 has to be slowest varying */
typedef u_mat_array (*my_mat_array)[4] ALIGN;  
#endif 


/* end overlays */


/* macros for spin basis note: assume first two rows are linearly independent except for gamma3 */

#define _sse_load_temp_1(chi) _sse_load((chi)[0])
#define _sse_load_temp_2(chi) _sse_load((chi)[1])
#define _sse_load_up_temp_1(chi) _sse_load_up((chi)[0])
#define _sse_load_up_temp_2(chi) _sse_load_up((chi)[1])

#define _sse_store_temp_1(chi) _sse_store((chi)[0])
#define _sse_store_temp_2(chi) _sse_store((chi)[1])
#define _sse_store_up_temp_1(chi) _sse_store_up((chi)[0])
#define _sse_store_up_temp_2(chi) _sse_store_up((chi)[1])


#define _sse_psum_set(psi) \
	  rs[0] = (psi)[0]; \
	  rs[1] = (psi)[1]; \
	  rs[2] = (psi)[2]; \
	  rs[3] = (psi)[3]


/* gamma 0 */

#define _sse_42_1_gamma0_minus(sp) \
      _sse_load((sp)[0]); \
      _sse_load_up((sp)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

#define _sse_24_1_gamma0_minus_set() \
      _sse_store_up(rs[0]);\
      _sse_vector_i_mul_up();\
      _sse_store_up(rs[3]) 

#define _sse_24_1_gamma0_minus_add() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[3]);\
      _sse_vector_i_mul();\
      _sse_vector_add();\
      _sse_store(rs[3]) 
	  
#define _sse_42_2_gamma0_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

#define _sse_24_2_gamma0_minus_set() \
	  _sse_store_up(rs[1]);\
      _sse_vector_i_mul_up();\
      _sse_store_up(rs[2])

#define _sse_24_2_gamma0_minus_add() \
	  _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store(rs[1]);\
      _sse_load(rs[2]);\
      _sse_vector_i_mul();\
      _sse_vector_add();\
      _sse_store(rs[2])

#define _sse_42_1_gamma0_plus(sm) \
      _sse_load((sm)[0]);\
      _sse_load_up((sm)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_24_1_gamma0_plus_set() \
	  _sse_store_up(rs[0]);\
      _sse_vector_i_mul_neg_up();\
      _sse_store_up(rs[3])

#define _sse_24_1_gamma0_plus_add() \
	  _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[3]);\
      _sse_vector_i_mul();\
      _sse_vector_sub();\
      _sse_store(rs[3])

#define _sse_42_2_gamma0_plus(sm) \
	  _sse_load((sm)[1]);\
      _sse_load_up((sm)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_24_2_gamma0_plus_set() \
      _sse_store_up(rs[1]);\
      _sse_vector_i_mul_neg_up();  \
      _sse_store_up(rs[2])

#define _sse_24_2_gamma0_plus_add() \
       _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store(rs[1]);\
      _sse_load(rs[2]);\
      _sse_vector_i_mul();  \
      _sse_vector_sub();\
      _sse_store(rs[2])



/* gamma 1 */


#define _sse_42_1_gamma1_minus(sp) \
      _sse_load((sp)[0]);\
      _sse_load_up((sp)[3]);\
      _sse_vector_add()

#define _sse_24_1_gamma1_minus() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[3]);\
      _sse_vector_add();\
      _sse_store(rs[3])
	  
#define _sse_42_2_gamma1_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[2]);\
      _sse_vector_sub()

#define _sse_24_2_gamma1_minus() \
	  _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store(rs[1]);\
      _sse_load(rs[2]);\
      _sse_vector_sub();\
      _sse_store(rs[2])

#define _sse_42_1_gamma1_plus(sm) \
      _sse_load((sm)[0]);\
      _sse_load_up((sm)[3]);\
      _sse_vector_sub()

#define _sse_24_1_gamma1_plus() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[3]);\
      _sse_vector_sub();\
      _sse_store(rs[3])

#define _sse_42_2_gamma1_plus(sm) \
	  _sse_load((sm)[1]);\
      _sse_load_up((sm)[2]);\
      _sse_vector_add()

#define _sse_24_2_gamma1_plus() \
       _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store(rs[1]);\
      _sse_load(rs[2]);\
      _sse_vector_add();\
      _sse_store(rs[2])





/* gamma 2 */


#define _sse_42_1_gamma2_minus(sp) \
      _sse_load((sp)[0]);\
      _sse_load_up((sp)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

#define _sse_24_1_gamma2_minus() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[2]);\
      _sse_vector_i_mul();   \
      _sse_vector_add();\
      _sse_store(rs[2])
	  
#define _sse_42_2_gamma2_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_24_2_gamma2_minus() \
	   _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store(rs[1]);\
      _sse_load(rs[3]);\
      _sse_vector_i_mul();   \
      _sse_vector_sub();\
      _sse_store(rs[3])

#define _sse_42_1_gamma2_plus(sm) \
      _sse_load((sm)[0]);\
      _sse_load_up((sm)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_add()

#define _sse_24_1_gamma2_plus() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[2]);\
      _sse_vector_i_mul();   \
      _sse_vector_sub();\
      _sse_store(rs[2]);

#define _sse_42_2_gamma2_plus(sm) \
	  _sse_load((sm)[1]);\
      _sse_load_up((sm)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

#define _sse_24_2_gamma2_plus() \
      _sse_load(rs[1]);\
      _sse_vector_add();\
      _sse_store(rs[1]);\
      _sse_load(rs[3]);\
      _sse_vector_i_mul();     \
      _sse_vector_add();\
      _sse_store(rs[3])





/* gamma 3 */
#define _sse_42_1_gamma3_minus(sp) \
	  _sse_load((sp)[0]); \
	  _sse_load_up((sp)[2]); \
      _sse_vector_sub()

#define _sse_24_1_gamma3_minus_set() \
	  _sse_load(rs[0]);\
      _sse_vector_add();\
       _sse_store((*rn)[0]);\
      _sse_load(rs[2]);\
      _sse_vector_sub();\
       _sse_store((*rn)[2])

#define _sse_24_1_gamma3_minus_add() \
	  _sse_load(rs[0]);\
      _sse_vector_add();\
       _sse_store(rs[0]);\
      _sse_load(rs[2]);\
      _sse_vector_sub();\
       _sse_store(rs[2])
	  
#define _sse_42_2_gamma3_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[3]);\
      _sse_vector_sub()

#define _sse_24_2_gamma3_minus_set() \
	  _sse_load(rs[1]);\
      _sse_vector_add();\
       _sse_store((*rn)[1]);\
      _sse_load(rs[3]);\
      _sse_vector_sub();\
      _sse_store((*rn)[3])

#define _sse_24_2_gamma3_minus_add() \
	  _sse_load(rs[1]);\
      _sse_vector_add();\
       _sse_store(rs[1]);\
      _sse_load(rs[3]);\
      _sse_vector_sub();\
      _sse_store(rs[3])

#define _sse_42_1_gamma3_plus(sm) \
      _sse_load((sm)[0]);\
      _sse_load_up((sm)[2]);\
      _sse_vector_add()

#define _sse_24_1_gamma3_plus_set() \
	  _sse_load(rs[0]);\
      _sse_vector_add();\
       _sse_store((*rn)[0]);\
      _sse_load(rs[2]);\
      _sse_vector_add();\
      _sse_store((*rn)[2])

#define _sse_24_1_gamma3_plus_add() \
	  _sse_load(rs[0]);\
      _sse_vector_add();\
       _sse_store(rs[0]);\
      _sse_load(rs[2]);\
      _sse_vector_add();\
      _sse_store(rs[2])

#define _sse_42_2_gamma3_plus(sm) \
	  _sse_load((sm)[1]);\
      _sse_load_up((sm)[3]);\
      _sse_vector_add()

#define _sse_24_2_gamma3_plus_set() \
       _sse_load(rs[1]);\
      _sse_vector_add();\
       _sse_store((*rn)[1]);\
      _sse_load(rs[3]);\
      _sse_vector_add();\
       _sse_store((*rn)[3])

#define _sse_24_2_gamma3_plus_add() \
       _sse_load(rs[1]);\
       _sse_vector_add();\
       _sse_store(rs[1]);\
       _sse_load(rs[3]);\
       _sse_vector_add();\
       _sse_store(rs[3]);\



static int init=0;
static sse_double fact1,fact2;
static spinor_array rs __attribute__ ((aligned (16)));


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






#define DECLARE_COMMON_STUFF \
 
/* watch out for sm and sp */


/* this routine is similar to wnxtsu3dslash, except instead of handling the second site's worth in the same loop, the second
spin component's worth must be handled seperately */
void decomp_plus(size_t lo, size_t hi, int id, const void *ptr)
{

  int ix,iy,iz;
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  static int ix1,iy1,iy2,iz1;
  const int low = icolor_start[cb]+lo; 
  const int high = icolor_start[cb]+hi;
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  u_mat_array *up ALIGN,*um ALIGN;
  halfspinor_array *s2 ALIGN,*s3 ALIGN, *s4 ALIGN,  *temp0,*temp1,*temp2,*temp3;
  spinor_array *s ALIGN, *s1 ALIGN, *sp ALIGN,*sm ALIGN,rs ALIGN,*rn ALIGN;
  sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; 
  sse_double _minus_one ALIGN = {-1.0, -1.0};

  for (ix1=low;ix1<high;ix1+=1) 
  {
#ifdef PREFDIST
#undef PREFDIST
#endif
#define PREFDIST 4
    sp=&psi[ix1];
    s1=&psi[ix1+PREFDIST];
    _prefetch_spinor(s1);
      
/******************************* direction +0 *********************************/	   
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,0);
    /*spin decomposition of first component of halfspinor */
    _sse_42_1_gamma0_minus(*sp);
    _sse_store_temp_1(*s3);
    /*spin decomposition of first component of halfspinor */
    _sse_42_2_gamma0_minus(*sp);
    _sse_store_temp_2(*s3);

/******************************* direction +1 *********************************/
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,1);
	   
    _sse_42_1_gamma1_minus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma1_minus(*sp);
    _sse_store_temp_2(*s3);
   
/******************************* direction +2 *********************************/
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,2);
	   
    _sse_42_1_gamma2_minus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma2_minus(*sp);
    _sse_store_temp_2(*s3);

/******************************* direction +3 *********************************/
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,3);
	   
    _sse_42_1_gamma3_minus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma3_minus(*sp);
    _sse_store_temp_2(*s3);
  }
}


void decomp_hvv_plus(size_t lo, size_t hi, int id, const void *ptr)
{

  int ix,iy,iz;
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  static int ix1,iy1,iy2,iz1;
  const int low = icolor_start[cb]+lo; 
  const int high = icolor_start[cb]+hi;
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  u_mat_array *up ALIGN,*um ALIGN;
  halfspinor_array *s2 ALIGN,*s3 ALIGN, *s4 ALIGN, *temp0,*temp1,*temp2,*temp3;
  spinor_array *s ALIGN, *s1 ALIGN, *sp ALIGN,*sm ALIGN,rs ALIGN,*rn ALIGN;
  sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; 
  sse_double _minus_one ALIGN = {-1.0, -1.0};
   
  um=&_gauge_field0_0(low);
  _prefetch_su3(um+1);
  for (ix1=low;ix1<high;ix1++) 
  {
    sm=&psi[ix1];
     
/******************************* direction -0 *********************************/
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,0);
	   
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
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,1);
    
    _sse_42_1_gamma1_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);
    
    _sse_42_2_gamma1_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);


/******************************* direction -2 *********************************/

    um++;
    _prefetch_su3(um+1);

    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,2);
    
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

    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,3);
	   

    _sse_42_1_gamma3_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);

    if(ix1 < hi) _prefetch_su3(um+1);
    else _prefetch_su3(&_gauge_field0_0(low));

    _sse_42_2_gamma3_plus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);

    um++;
       
	   
  }
}


void mvv_recons_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  int ix,iy,iz;
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  static int ix1,iy1,iy2,iz1;
  const int low = icolor_start[cb]+lo; 
  const int high = icolor_start[cb]+hi;
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  u_mat_array *up ALIGN,*um ALIGN;
  halfspinor_array *s2 ALIGN,*s3 ALIGN, *s4 ALIGN, *temp0,*temp1,*temp2,*temp3 ; 
  spinor_array *s ALIGN, *s1 ALIGN, *sp ALIGN,*sm ALIGN,rs ALIGN,*rn ALIGN;
  sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; 
  sse_double _minus_one ALIGN = {-1.0, -1.0};
   	
   
 
  up=&_gauge_field0_0(low);
  _prefetch_su3(up+1);
  s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,lo,0);


  for (ix1=low;ix1<high;ix1++) 
  {
    rn=&psi[ix1];
	 
      
/******************************* direction +0 *********************************/	

    s4 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,1);
    _prefetch_spinor(s4);



    _sse_load_temp_1(*s3);
    _sse_su3_multiply(*up);
    _sse_24_1_gamma0_minus_set();
	  

    _sse_load_temp_2(*s3);
    _sse_su3_multiply(*up);
    _sse_24_2_gamma0_minus_set();

    up++;
	   
/******************************* direction +1 *********************************/

     
	   

	 

    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,2);
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


	 

    s4 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,3);
    _prefetch_spinor(s4);

    _sse_load_temp_1(*s3);
    _sse_su3_multiply(*up);
    _sse_24_1_gamma2_minus();
	  

    _sse_load_temp_2(*s3);
    _sse_su3_multiply(*up);
    _sse_24_2_gamma2_minus();





    up++;
    /******************************* direction +3 *********************************/     
	 

    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1+1,0);
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
  /* 	printf("\n &a:%x",&a);*/
  int ix,iy,iz;
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  static int ix1,iy1,iy2,iz1;
  const int low = icolor_start[cb]+lo; 
  const int high = icolor_start[cb]+hi;
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  u_mat_array *up ALIGN,*um ALIGN;
  halfspinor_array *s2 ALIGN,*s3 ALIGN, *s4 ALIGN, *temp1,*temp2,*temp3; 
  spinor_array *s ALIGN, *s1 ALIGN, *sp ALIGN,*sm ALIGN,rs ALIGN,*rn ALIGN;
  sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; 
  sse_double _minus_one ALIGN = {-1.0, -1.0};


 
  s3 = chi + halfspinor_buffer_offset(RECONS_GATHER,lo,0);

#ifdef PREFDIST
#undef PREFDIST
#endif
#define PREFDIST 1
  /*printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
    printf("\n &psi[0]:%x",(&psi[0]));*/
  for (ix1=low;ix1<high;ix1++) 
  {
    rn=&psi[ix1];
	 
    s1=&psi[ix1+1];
    _prefetch_spinor(s1);

#define rs (*rn)
	   
    s4 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,1);
    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,1);
    _prefetch_single(s2);
     

    _sse_load_up_temp_1(*s3);
	   
    _sse_24_1_gamma0_plus_add();
	  

    _sse_load_up_temp_2(*s3);
	   
    _sse_24_2_gamma0_plus_add();

	   

    s3 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,2);
    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,2);
    _prefetch_single(s2);
 

    _sse_load_up_temp_1(*s4);
	   
    _sse_24_1_gamma1_plus();
	  

    _sse_load_up_temp_2(*s4);
	  
    _sse_24_2_gamma1_plus();

	   

    s4 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,3);
    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,3);
    _prefetch_single(s2);
	   

    _sse_load_up_temp_1(*s3);
	  
    _sse_24_1_gamma2_plus();
	  

    _sse_load_up_temp_2(*s3);
	   
    _sse_24_2_gamma2_plus();

	   
       	
    s3 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+1,0);
    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+1+PREFDIST,0);
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
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  static int ix1,iy1,iy2,iz1;
  const int low = icolor_start[cb]+lo; 
  const int high = icolor_start[cb]+hi;
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  u_mat_array *up ALIGN,*um ALIGN;
  halfspinor_array *s2 ALIGN,*s3 ALIGN, *s4 ALIGN,*temp0, *temp1,*temp2,*temp3;
  spinor_array *s ALIGN, *s1 ALIGN, *sp ALIGN,*sm ALIGN,rs ALIGN,*rn ALIGN;
  sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; 
  sse_double _minus_one ALIGN = {-1.0, -1.0};
 
  
#ifdef PREFDIST
#undef PREFDIST
#endif
#define PREFDIST 3
  /*printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
	printf("\n &psi[0]:%x",(&psi[0]));*/
  for (ix1=low;ix1<high;ix1++) 
  {
    rn=&psi[ix1];
	 
    s1=&psi[ix1+PREFDIST];
    _prefetch_nta_spinor(s1);

	  

    /* first spin component of result */
    temp0 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,0);
	  
	    

    _sse_load_up_temp_1(*(temp0));
    _sse_load((*rn)[0]);
	   
    _sse_vector_add();

    temp1 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,1);

    _sse_load_up_temp_1(*(temp1));
    _sse_vector_add();

    temp2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,2);

    _sse_load_up_temp_1(*(temp2));
    _sse_vector_add();

    temp3 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,3);

    _sse_load_up_temp_1(*(temp3));
    _sse_vector_add();
    _sse_store((*rn)[0]);

/* second spin component of result */
    s2= chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,0);
    _prefetch_single(s2);
     

    _sse_load_up_temp_2(*(temp0));
    _sse_load((*rn)[1]);
	   
    _sse_vector_add();

    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,1);
    _prefetch_single(s2);

    _sse_load_up_temp_2(*(temp1));
    _sse_vector_add();

    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,2);
    _prefetch_single(s2);

    _sse_load_up_temp_2(*(temp2));
    _sse_vector_add();

    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,3);
    _prefetch_single(s2);

    _sse_load_up_temp_2(*(temp3));
    _sse_vector_add();
    _sse_store((*rn)[1]);



		

    /* third spin component, here it gets tricky */
	  
	   
	   
     

    _sse_load_up_temp_2(*(temp0));
    _sse_load((*rn)[2]);
	   
    _sse_vector_i_mul_up();
    _sse_vector_sub();

	  
	  

    _sse_load_up_temp_2(*(temp1));
		

    _sse_vector_add();

      
      

    _sse_load_up_temp_1(*(temp2));

    _sse_vector_i_mul();
    _sse_vector_sub();

    

    _sse_load_up_temp_1(*(temp3));

    _sse_vector_add();

    _sse_store((*rn)[2]);

/* fourth spin component, again it gets tricky */
	  
	   
	   
     

    _sse_load_up_temp_1(*(temp0));
    _sse_load((*rn)[3]);
	   
    _sse_vector_i_mul_up();
    _sse_vector_sub();

	  
	  

    _sse_load_up_temp_1(*(temp1));
		

    _sse_vector_sub();

      
      

    _sse_load_up_temp_2(*(temp2));

    _sse_vector_i_mul_up();
    _sse_vector_add();

    

    _sse_load_up_temp_2(*(temp3));

    _sse_vector_add();

    _sse_store((*rn)[3]);

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
   
  int ix,iy,iz;
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  static int ix1,iy1,iy2,iz1;
  const int low = icolor_start[cb]+lo; 
  const int high = icolor_start[cb]+hi;
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  u_mat_array *up ALIGN,*um ALIGN;
  halfspinor_array *s2 ALIGN,*s3 ALIGN, *s4 ALIGN; 
  spinor_array *s ALIGN, *s1 ALIGN, *sp ALIGN,*sm ALIGN,rs ALIGN,*rn ALIGN;
  sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; 
  sse_double _minus_one ALIGN = {-1.0, -1.0};
    
 

  for (ix1=low;ix1<high;ix1+=1) 
  {

#ifdef PREFDIST
#undef PREFDIST
#endif
#define PREFDIST 4
    sp=&psi[ix1];
    s1=&psi[ix1+PREFDIST];
    _prefetch_spinor(s1);
  
	  
	   
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,0);
	   
    _sse_42_1_gamma0_plus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma0_plus(*sp);
    _sse_store_temp_2(*s3);

    /*printf("chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s3)[0][0][0],(*s3)[1][0][0]);
 
      printf("chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s3)[0][2][0],(*s3)[1][2][0]); */

    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,1);
	   
    _sse_42_1_gamma1_plus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma1_plus(*sp);
    _sse_store_temp_2(*s3);
   

    /*printf("chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s3)[0][0][0],(*s3)[1][0][0]);
 
      printf("chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s3)[0][2][0],(*s3)[1][2][0]); */


    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,2);
	   
    _sse_42_1_gamma2_plus(*sp);
    _sse_store_temp_1(*s3);

    _sse_42_2_gamma2_plus(*sp);
    _sse_store_temp_2(*s3);


    /*printf("chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s3)[0][0][0],(*s3)[1][0][0]);
 
      printf("chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s3)[0][2][0],(*s3)[1][2][0]); */


    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,3);
	   
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
  int ix,iy,iz;
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  static int ix1,iy1,iy2,iz1;
  const int low = icolor_start[cb]+lo; 
  const int high = icolor_start[cb]+hi;
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  u_mat_array *up ALIGN,*um ALIGN;
  halfspinor_array *s2 ALIGN,*s3 ALIGN, *s4 ALIGN; 
  spinor_array *s ALIGN, *s1 ALIGN, *sp ALIGN,*sm ALIGN,rs ALIGN,*rn ALIGN;
  sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; 
  sse_double _minus_one ALIGN = {-1.0, -1.0};
   
 
  um=&_gauge_field0_0(low);
  _prefetch_su3(um+1);
  for (ix1=low;ix1<high;ix1++) 
  {
    sm=&psi[ix1];
	 
     

       	   
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,0);
	   
    _sse_42_1_gamma0_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);

    _sse_42_2_gamma0_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);

    um++;
    _prefetch_su3(um+1);
    /* in decomp_hvv */
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,1);
	   

    _sse_42_1_gamma1_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);

    _sse_42_2_gamma1_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);

    um++;
    _prefetch_su3(um+1);

    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,2);
	   

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

    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,3);
	   

    _sse_42_1_gamma3_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_1(*s3);

    if(ix1 < hi) _prefetch_su3(um+1);
    else _prefetch_su3(&_gauge_field0_0(low));

    _sse_42_2_gamma3_minus(*sm);
    _sse_su3_inverse_multiply(*um);
    _sse_store_up_temp_2(*s3);

    um++;
       
	   
  }
}


void mvv_recons_minus(size_t lo, size_t hi, int id, const void *ptr)
{
   	
  int ix,iy,iz;
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  static int ix1,iy1,iy2,iz1;
  const int low = icolor_start[cb]+lo; 
  const int high = icolor_start[cb]+hi;
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  u_mat_array *up ALIGN,*um ALIGN; 
  halfspinor_array *s2 ALIGN,*s3 ALIGN, *s4 ALIGN,  *temp0,*temp1,*temp2,*temp3; 
  spinor_array *s ALIGN, *s1 ALIGN, *sp ALIGN,*sm ALIGN,rs ALIGN,*rn ALIGN;
  sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; 
  sse_double _minus_one ALIGN = {-1.0, -1.0};
   
 
  up=&_gauge_field0_0(low);
  _prefetch_su3(up+1);
  s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,low,0);


  for (ix1=low;ix1<high;ix1++) 
  {
    rn=&psi[ix1];
	 
      
    /* printf("g chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s3)[0][0][0],(*s3)[1][0][0]);
 
       printf("g chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s3)[0][2][0],(*s3)[1][2][0]); 
	   
       printf("u(ix1)[0][0][0]=%f \n", (*up)[0][0][0]);*/

    s4 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,1);
    _prefetch_spinor(s4);



    _sse_load_temp_1(*s3);
    _sse_su3_multiply(*up);
    _sse_24_1_gamma0_plus_set();
	  

    _sse_load_temp_2(*s3);
    _sse_su3_multiply(*up);
    _sse_24_2_gamma0_plus_set();

    /* printf("rs(ix1)[0]=%f,rs(ix1)[1] [0][0]=%f\n",(rs[0])[0][0],(rs[1])[0][0]);
 
       printf("rs(ix1)[2]=%f,rs(ix1)[3] [0][0]=%f\n",(rs[2])[0][0],(rs[3])[0][0]);*/

    up++;


    /* printf("chi(ix1)[0]=%f,chi(ix1)[1] [0][0]=%f\n",(*s4)[0][0][0],(*s4)[1][0][0]);
 
       printf("chi(ix1)[0]=%f,chi(ix1)[1] [2][0]=%f\n",(*s4)[0][2][0],(*s4)[1][2][0]); */

    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,2);
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

    s4 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,3);
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

    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1+1,0);
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
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  static int ix1,iy1,iy2,iz1;
  const int low = icolor_start[cb]+lo; 
  const int high = icolor_start[cb]+hi;
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  u_mat_array *up ALIGN,*um ALIGN;
  halfspinor_array *s2 ALIGN,*s3 ALIGN, *s4 ALIGN,  *temp0,*temp1,*temp2,*temp3;  
  spinor_array *s ALIGN, *s1 ALIGN, *sp ALIGN,*sm ALIGN,rs ALIGN,*rn ALIGN;
  sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; 
  sse_double _minus_one ALIGN = {-1.0, -1.0};
 
  s3 = chi + halfspinor_buffer_offset(RECONS_GATHER,lo,0);

#ifdef PREFDIST
#undef PREFDIST
#endif

#define PREFDIST 1
  /*printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
    printf("\n &psi[0]:%x",(&psi[0]));*/
  for (ix1=low;ix1<high;ix1++) 
  {
    rn=&psi[ix1];
	 
    s1=&psi[ix1+1];
    _prefetch_spinor(s1);

#define rs (*rn)
	   
    s4 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,1);
    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,1);
    _prefetch_single(s2);
     

    _sse_load_up_temp_1(*s3);
	   
    _sse_24_1_gamma0_minus_add();
	  

    _sse_load_up_temp_2(*s3);
	   
    _sse_24_2_gamma0_minus_add();

	   

    s3 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,2);
    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,2);
    _prefetch_single(s2);
 

    _sse_load_up_temp_1(*s4);
	   
    _sse_24_1_gamma1_minus();
	  

    _sse_load_up_temp_2(*s4);
	  
    _sse_24_2_gamma1_minus();

	   

    s4 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,3);
    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,3);
    _prefetch_single(s2);
	   

    _sse_load_up_temp_1(*s3);
	  
    _sse_24_1_gamma2_minus();
	  

    _sse_load_up_temp_2(*s3);
	   
    _sse_24_2_gamma2_minus();

	   
       	
    s3 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+1,0);
    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+1+PREFDIST,0);
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
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  const int cb = a->cb; 
  static int ix1,iy1,iy2,iz1;
  const int low = icolor_start[cb]+lo; 
  const int high = icolor_start[cb]+hi;
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  u_mat_array *up ALIGN,*um ALIGN;
  halfspinor_array *s2 ALIGN,*s3 ALIGN, *s4 ALIGN, *temp0,*temp1,*temp2,*temp3;   
  spinor_array *s ALIGN, *s1 ALIGN, *sp ALIGN,*sm ALIGN,rs ALIGN,*rn ALIGN;
  sse_int _sse_sgn ALIGN = {0x0,0x80000000,0x0,0x0}; 
  sse_double _minus_one ALIGN = {-1.0, -1.0};

 
  
#ifdef PREFDIST
#undef PREFDIST
#endif
#define PREFDIST 3
  /*printf("\n psi[0].c1.c1.re:%f",(psi[0]).c1.c1.re);
    printf("\n &psi[0]:%x",(&psi[0]));*/
  for (ix1=low;ix1<high;ix1++) 
  {
    rn=&psi[ix1];
	 
    s1=&psi[ix1+PREFDIST];
    _prefetch_nta_spinor(s1);

	  

    /* first spin component of result */
    temp0 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,0);
	  
	    

    _sse_load_up_temp_1(*(temp0));
    _sse_load((*rn)[0]);
	   
    _sse_vector_add();

    temp1 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,1);

    _sse_load_up_temp_1(*(temp1));
    _sse_vector_add();

    temp2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,2);

    _sse_load_up_temp_1(*(temp2));
    _sse_vector_add();

    temp3 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,3);

    _sse_load_up_temp_1(*(temp3));
    _sse_vector_add();
    _sse_store((*rn)[0]);

/* second spin component of result */
    s2= chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,0);
    _prefetch_single(s2);
     

    _sse_load_up_temp_2(*(temp0));
    _sse_load((*rn)[1]);
	   
    _sse_vector_add();

    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,1);
    _prefetch_single(s2);

    _sse_load_up_temp_2(*(temp1));
    _sse_vector_add();

    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,2);
    _prefetch_single(s2);

    _sse_load_up_temp_2(*(temp2));
    _sse_vector_add();

    s2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1+PREFDIST,3);
    _prefetch_single(s2);

    _sse_load_up_temp_2(*(temp3));
    _sse_vector_add();
    _sse_store((*rn)[1]);



		

    /* third spin component, here it gets tricky */
	  
	   
	   
     

    _sse_load_up_temp_2(*(temp0));
    _sse_load((*rn)[2]);
	   
    _sse_vector_i_mul_up();
    _sse_vector_add();

	  
	  

    _sse_load_up_temp_2(*(temp1));
		

    _sse_vector_sub();

      
      

    _sse_load_up_temp_1(*(temp2));

    _sse_vector_i_mul();
    _sse_vector_add();

    

    _sse_load_up_temp_1(*(temp3));

    _sse_vector_sub();

    _sse_store((*rn)[2]);

/* fourth spin component, again it gets tricky */
	  
	   
	   
     

    _sse_load_up_temp_1(*(temp0));
    _sse_load((*rn)[3]);
	   
    _sse_vector_i_mul_up();
    _sse_vector_add();

	  
	  

    _sse_load_up_temp_1(*(temp1));
		

    _sse_vector_add();

      
      

    _sse_load_up_temp_2(*(temp2));

    _sse_vector_i_mul_up();
    _sse_vector_sub();

    

    _sse_load_up_temp_2(*(temp3));

    _sse_vector_sub();

    _sse_store((*rn)[3]);

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

static halfspinor_array* chi1;
static halfspinor_array* chi2;
//#define TAIL1(chi,mymu) (chi+subgrid_vol_cb*(1+3*mymu))
//#define TAIL2(chi,mymu) (chi+0+subgrid_vol_cb*(2+3*mymu))

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
  int bound[2][2][Nd];
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




    
  make_shift_tables(icolor_start, bound);
  // The SSE code expects to have at least 2 sites after checkerboarding.
  if ( getSubgridVolCB() <= 1 )
  {
    QMP_error("This SSE Dslash expects there to be at least 2 subgrid sites after checkerboarding");
    QMP_abort(1);
  }

  icolor_end[0] = icolor_start[0] + getSubgridVolCB();
  icolor_end[1] = icolor_start[1] + getSubgridVolCB();

  /* Allocated space for the floating temps */
  /* Wasteful - allocate 3 times subgrid_vol_cb. Otherwise, need to pack the TAIL{1,2} halfspinor_buffer_offsets */
  nsize = 2*Nc*Ns2*sizeof(double)*getSubgridVolCB()*3*Nd;  /* Note 3x4 half-ferm temps */
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
    
  chi1 = (halfspinor_array*)QMP_get_memory_pointer(xchi1);
  chi2 = (halfspinor_array*)QMP_get_memory_pointer(xchi2);
  
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

      forw_msg[num][0] = QMP_declare_msgmem(chi1+getSubgridVolCB()*(1+3*mu), bound[0][0][mu]*sizeof(halfspinor_array));
      forw_msg[num][1] = QMP_declare_msgmem(chi1+getSubgridVolCB()*(2+3*mu), bound[0][0][mu]*sizeof(halfspinor_array));
      forw_mh[num][0]  = QMP_declare_receive_relative(forw_msg[num][1], mu, +1, 0);
      forw_mh[num][1]  = QMP_declare_send_relative(forw_msg[num][0], mu, -1, 0);
	
      if (bound[0][1][mu] == 0)
      {
	QMP_error("init_sse_dslash: type 0 message size is 0");
	QMP_abort(1);
      }

      back_msg[num][0] = QMP_declare_msgmem(chi2+getSubgridVolCB()*(1+3*mu), bound[0][1][mu]*sizeof(halfspinor_array));
      back_msg[num][1] = QMP_declare_msgmem(chi2+getSubgridVolCB()*(2+3*mu), bound[0][1][mu]*sizeof(halfspinor_array));
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
    free_shift_tables();
    
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

  if (initP == 0) {
    QMP_error("sse_su3dslash_wilson not initialised");
    QMP_abort(1);
  }

  if(isign==1) 
  {
    dispatch_to_threads(decomp_plus,
			(spinor_array*)psi,
			chi1,
			(my_mat_array)u,
			cb,
			getSubgridVolCB());

    _prefetch_su3((u+0+2*(0+3*(0+3*(0+4*(icolor_start[cb]))))));
    _prefetch_single(chi2);


    if (total_comm > 0)
      if (QMP_start(forw_all_mh) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	QMP_abort(1);
      }

    /*other checkerboard's worth */
    dispatch_to_threads(decomp_hvv_plus,
			(spinor_array*)psi,
			chi2,
			(my_mat_array)u,
			cb,
			getSubgridVolCB());
	
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

    dispatch_to_threads(mvv_recons_plus,
			(spinor_array*)res,
			chi1,
			(my_mat_array)u,
			1-cb,
			getSubgridVolCB());
	
    if (total_comm > 0)
      if (QMP_wait(back_all_mh) != QMP_SUCCESS)
      {
	QMP_error("wnxtsu3dslash: QMP_wait failed in backward direction");
	QMP_abort(1);
      }


    dispatch_to_threads(recons_plus,
			(spinor_array*)res, 
			chi2,
			(my_mat_array)u,	
			1-cb,
			getSubgridVolCB());
  }		
  
  if(isign==-1) 
  {
    dispatch_to_threads(decomp_minus,
			(spinor_array*)psi,
			chi1,
			(my_mat_array)u,
			cb,
			getSubgridVolCB());
    
    _prefetch_su3((u+0+2*(0+3*(0+3*(0+4*(icolor_start[cb]))))));
    _prefetch_single(chi2);
   
    if (total_comm > 0)
      if (QMP_start(forw_all_mh) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	QMP_abort(1);
      }
	
	
    /*other checkerboard's worth */
    dispatch_to_threads(decomp_hvv_minus,
			(spinor_array*)psi,
			chi2,
			(my_mat_array)u,
			cb,
			getSubgridVolCB());
    
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
    dispatch_to_threads(mvv_recons_minus,
			(spinor_array*)res,
			chi1,
			(my_mat_array)u,
			1-cb,
			getSubgridVolCB());
    
    _prefetch_single(chi2);

    if (total_comm > 0)
      if (QMP_wait(back_all_mh) != QMP_SUCCESS)
      {
	QMP_error("wnxtsu3dslash: QMP_wait failed in backward direction"); 
	QMP_abort(1);
      }

    dispatch_to_threads(recons_minus,
			(spinor_array*)res, 
			chi2,
			(my_mat_array)u,	
			1-cb,
			getSubgridVolCB());
  }		
}


#ifdef __cplusplus
}
#endif
