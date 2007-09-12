/*******************************************************************************
 * $Id: sse_su3dslash_64bit_scalar.c,v 1.1.1.1 2007-09-12 19:33:13 bjoo Exp $
 * 
 * Action of the 64bit single-node Wilson-Dirac operator D_w on a given spinor field
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


#include <sse_config.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef DMALLOC
#include <stdlib.h>
#else
#include <dmalloc.h>
#endif
#include <stdio.h>
#include <math.h>

/* externally callable routine: tnxtsu3dslash */
/* requires gauge fields packed by no_funnystuff_pack_gauge_field of packer64.c */
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

/*  U	      Gauge field					(Read) */
/*  Psi	      Pseudofermion field				(Read) */
/*  Res	      Pseudofermion field				(Write) */
/*		      + */
/*  ISign      D' or D'  ( +1 | -1 ) respectively		(Read) */
/*  CB	      Checkerboard of input vector			(Read) */


#if defined SSE2
#warning "using sse64 stuff"
#include <sse64.h>

#if defined P4
  #define BASE 0x3f
#else
  #define BASE 0x1f
#endif

#include <sse_align.h>

/* now overlays for spinors as arrays or structs */


#ifdef SZIN
#define SZIN_SHIFT
#define SPINORS_AS_ARRAYS
#define GAUGE_AS_ARRAYS
#define NO_U_PACK
#endif


#ifdef SZIN

static int total_vol = 1;
static int initP = 0;

extern void make_shift_tables(int *soffsets, int icolor_start[2], const int lat_size[4]);

#endif



#ifdef SZIN_SHIFT  /*look at sse.h to see what else this implies..for now defining SZIN means spinors, gauge as arrays */

static int *soffsets;
static int icolor_start[2];    /* starting site for each coloring (cb) */
static int icolor_end[2];      /* end site for each coloring (cb) */

#define iup(mysite,mymu) soffsets[mymu + Nd*(mysite + total_vol*(1))]
#define idn(mysite,mymu) soffsets[mymu + Nd*(mysite + total_vol*(0))]

/*isign = 0 => *(soffsets+mysite + global_vol_cb*2*(1-cb)+global_vol_cb*2*2*dir
  isign = 1 =>  *(soffsets+mysite + global_vol_cb*2*(1-cb)+global_vol_cb*2*2*dir+global_vol_cb) */
#else


#define iup(mysite,mymu) iup[mysite][mymu]
#define idn(mysite,mymu) idn[mysite][mymu]


#endif



#define _gauge_field0_0(mysite) gauge_field[mysite][0]
#define _gauge_field0_1(mysite) gauge_field[mysite][1]
#define _gauge_field0_2(mysite) gauge_field[mysite][2]
#define _gauge_field0_3(mysite) gauge_field[mysite][3]



 
/* now overlays for spinors as arrays or structs */
typedef double chi_double[2] __attribute__ ((aligned (16)));
typedef chi_double chi_two[2] __attribute__ ((aligned (16)));
typedef double u_mat_array[3][3][2]  ALIGN;  /* color color re/im */ 
typedef double spinor_array[4][3][2] ALIGN; /* Nspin4 color re/im */
typedef chi_two chi_array[3]    ALIGN; /*..color Nspin2 re/im ::note:: color slowest varying */
typedef u_mat_array (*my_mat_array)[4] ALIGN;  

#ifdef SPINORS_AS_ARRAYS
#define MY_SPINOR spinor_array
#define MY_SSE_VECTOR chi_array
#define MY_SSE_DOUBLE sse_double  
#define _c1__ [0] 
#define _c2__ [1] 
#define _c3__ [2]
#define _c4__ [3]
#define rs_c1__ rs[0]
#define rs_c2__ rs[1]
#define rs_c3__ rs[2]
#define rs_c4__ rs[3]

#else
#define MY_SPINOR spinor
#define MY_SSE_VECTOR sse_vector
#define MY_SSE_DOUBLE sse_double
#define _c1__ .c1
#define _c2__ .c2
#define _c3__ .c3
#define _c4__ .c4
#define rs_c1__ rs.c1
#define rs_c2__ rs.c2
#define rs_c3__ rs.c3
#define rs_c4__ rs.c4

#endif

/* now overlays for gauge matrices as arrays or structs */
#ifdef GAUGE_AS_ARRAYS


#define MY_GAUGE u_mat_array
#define MY_GAUGE_ARRAY my_mat_array
#else

#define MY_GAUGE su3

#endif




#ifdef SPINORS_AS_ARRAYS

#define MY_SPINOR spinor_array
#define MY_SSE_VECTOR chi_array
#define _c1__ [0]
#define _c2__ [1]
#define _c3__ [2]
#define _c4__ [3]


#else
#warning doing .c1 overlays
#define MY_SPINOR spinor
#define MY_SSE_VECTOR sse_vector
#define MY_SSE_DOUBLE sse_double
#define _c1__ .c1
#define _c2__ .c2
#define _c3__ .c3
#define _c4__ .c4

#endif

/* now overlays for gauge matrices as arrays or structs */
#ifdef GAUGE_AS_ARRAYS



#define MY_GAUGE u_mat_array

#else

#define MY_GAUGE su3

#endif



/* end overlays */


/* macros for spin basis note: assume first two rows are linearly independent except for gamma3 */
/*here we must handle individual spin components seperately because only one spin component fits in an xmm register...*/

#ifndef NON_SZIN_BASIS
/* use SZIN spin basis */



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




#else
/* other spin basis */

#endif

/* end spin basis section */

static int init=0;
static sse_double fact1,fact2;
static MY_SPINOR rs __attribute__ ((aligned (16)));



/***************** start of initialization routine ***************************************/

#define Nd 4

void init_sse_su3dslash(const int latt_size[])
{
  int mu, total_vol_cb;

  /* If we are already initialised, then increase the refcount and return */
  if (initP > 0) 
  {
    initP++;
    return;
  }

//  printf("init_sse_su3dslash: enter\n");

  /* Check problem and subgrid size */
  for(mu=0; mu < Nd; mu++) 
    if ( latt_size[mu] & 1 != 0 )
    {
      fprintf(stderr,"This SSE Dslash only supports even problem sizes. Here the lattice is odd in dimension %d with length %d\n", mu, latt_size[mu]);
      exit(1);
    }


  total_vol = latt_size[0];
  for(mu = 1; mu < Nd; mu++)
    total_vol *= latt_size[mu];
  total_vol_cb = total_vol >> 1;

  /* Construct all the shift tables needed */
  if ((soffsets = (int *)malloc(Nd*total_vol*2*sizeof(int))) == 0)
  {
    fprintf(stderr,"init_sse_su3dslash: could not initialize soffsets\n");
    exit(1);
  }
    
  make_shift_tables(soffsets, icolor_start, latt_size);
  icolor_end[0] = icolor_start[0] + total_vol_cb;
  icolor_end[1] = icolor_start[1] + total_vol_cb;

  initP = 1;

//  printf("init_sse_su3dslash: exit\n");
}


void free_sse_su3dslash(void)
{
  /* If we are uninitialised just return */
  if (initP == 0) {
    return;
  }

  /* Otherwise decrease the refcount */
  initP--;

  /* If the refcount has now hit 0 then free stuff */
  if( initP == 0 )
    free(soffsets);

}

/***************** end of initialization routine ***************************************/


/* prototypes */

void D_psi_fun(size_t lo,size_t hi, int id, const void *ptr);

void D_psi_fun_plus(size_t lo,size_t hi, int id, const void *ptr);


void D_psi_fun_minus(size_t lo,size_t hi, int id, const void *ptr);

typedef struct {
  MY_SPINOR *psi;  /* input spinor */
  MY_SPINOR *res;   /* output spinor */
  MY_GAUGE (*u)[4];  /* gauge field associated with the output checkerboard */
  MY_GAUGE (*u2)[4];  /* gauge field associated with the input checkerboard */
  int cb;             /* output checkerboard */
  /* the output array                         */
} Arg_s;


/* Stripped out the sml_scall from this version -- non threaded */
/* Hence n = 0. => lo = 0, hi = volume2, */

#define smpscaller2(a,bleah,spinfun2,chifun2,u3,cb2) \
    a.psi = spinfun2;\
    a.res = chifun2;\
    a.u = u3;  \
    a.cb = cb2; \
    (*bleah)(icolor_start[cb2], icolor_end[cb2], 0, &a);


void sse_su3dslash_wilson(double *u, double *psi, double *res, int isign, int cb)
{
  Arg_s a;

  if (isign == 1) 
  {
//    printf("sse_su3dslash_wilson: isign=+1, output_cb=%d\n",1-cb);

    smpscaller2(a, D_psi_fun_plus, 
		(MY_SPINOR*)psi,
		(MY_SPINOR*)res,
		(MY_GAUGE_ARRAY)u,
		1-cb);
  }

  if( isign == -1) 
  {
//    printf("sse_su3dslash_wilson: isign=-1, output_cb=%d\n",1-cb);

    smpscaller2(a, D_psi_fun_minus, 
		(MY_SPINOR*)psi,
		(MY_SPINOR*)res,
		(MY_GAUGE_ARRAY)u,
		1-cb);
  }

//  printf("sse_su3dslash_wilson: exit\n");
}

 


void D_psi_fun_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  int ix,iy,iz;
  const Arg_s *a = (const Arg_s*)ptr; 
  static int ix1,iy1,iy2,iz1;
  const int low = lo; 
  const int high = hi;
  MY_GAUGE (*gauge_field)[4] ALIGN = a->u;
  MY_SPINOR *psi = a->psi;
  MY_SPINOR *res = a->res;
  MY_GAUGE *up,*um;
  MY_SPINOR *s,*sp,*sm,*rn;            
  MY_SPINOR rs ALIGN;
  sse_int _sse_sgn ALIGN ={0x0,0x80000000,0x0,0x0};
  sse_double _minus_one ALIGN = {-1.0, -1.0};
  sse_double _conj ALIGN = {1.0, -1.0};

  /*printf("low:%i high:%i id:%i, ptr:%x\n", low, high, id, ptr);*/


  iy=iup(low,0);
  sp=&psi[iy];
  up=&_gauge_field0_0(low);
   
/************************ loop over all lattice sites *************************/
   
  for (ix=low;ix<high;ix++) 
  {
    /*s=&psi[ix];
      _prefetch_spinor(s);*/
    /*printf("low:%i high:%i ix:%i id:%i, ptr:%x\n", low, high,ix, id, ptr); */
/******************************* direction +0 *********************************/

    /*1 - gamma(0) */
    /* this routine is very similar to inxtsu3dslash, except that instead of handling the second site's worth, the routine 
       has to handle the second spin component's worth separately */

    iy=idn(ix,0);
      
    sm=&psi[iy];
    _prefetch_spinor(sm);  
    /*spin projection */
    /*output is one spin component of a halfspinor */
    _sse_42_1_gamma0_minus(*sp);
    /*multiplication with SU(3) * one component of a halfspinor*/
    _sse_su3_multiply((*up));
    /* spin reconstruction and part of the sum over directions */
    /*reconstruct a second spin component from the one you have from the multiply, do part of the sum over directions
      two more still to do for this direction */
    _sse_24_1_gamma0_minus_set();
      
     
      
    um=&_gauge_field0_0(iy);
    _prefetch_su3(um);
    /*spin project with the two other input spin components, get one spin component of a halfspinor  */
    _sse_42_2_gamma0_minus(*sp);
      
    /* multiply this one spin component times an SU(3) matrix */
    _sse_su3_multiply((*up));
    /*reconstruct a second spin component from the one you have from the multiply,do part of the sum over directions */
    _sse_24_2_gamma0_minus_set();
    /*whew, now all four spin components for this direction are accounted for. Had fun? Hope so, because the
      next directions' worth is coming at ya! */
      

/******************************* direction -0 *********************************/

    /*1 + gamma(0) */
    iy=iup(ix,1);

    sp=&psi[iy];
    _prefetch_spinor(sp);

    _sse_42_1_gamma0_plus(*sm);

      
    _sse_su3_inverse_multiply((*um));
      
    _sse_24_1_gamma0_plus_add();

      
      
    up+=1;
    _prefetch_su3(up);
      
    _sse_42_2_gamma0_plus(*sm);
     
      
    _sse_su3_inverse_multiply((*um));
      
    _sse_24_2_gamma0_plus_add();
      
/******************************* direction +1 *********************************/

    iy=idn(ix,1);
      
    sm=&psi[iy];
    _prefetch_spinor(sm);

    _sse_42_1_gamma1_minus(*sp);

    _sse_su3_multiply((*up));

    _sse_24_1_gamma1_minus();


    um=&_gauge_field0_1(iy);
    _prefetch_su3(um);

    _sse_42_2_gamma1_minus(*sp);

    _sse_su3_multiply((*up));

    _sse_24_2_gamma1_minus();
          

/******************************* direction -1 *********************************/

    iy=iup(ix,2);

    sp=&psi[iy];
    _prefetch_spinor(sp);

    _sse_42_1_gamma1_plus(*sm);
      
      
    _sse_su3_inverse_multiply((*um));
      
    _sse_24_1_gamma1_plus();

     

    up+=1;
    _prefetch_su3(up);

    _sse_42_2_gamma1_plus(*sm);
      
      
    _sse_su3_inverse_multiply((*um));
      
    _sse_24_2_gamma1_plus();
     

/******************************* direction +2 *********************************/

    iy=idn(ix,2);

    sm=&psi[iy];
    _prefetch_spinor(sm);

    _sse_42_1_gamma2_minus(*sp);


    _sse_su3_multiply((*up));

    _sse_24_1_gamma2_minus();

      
    um=&_gauge_field0_2(iy);
    _prefetch_su3(um);

    _sse_42_2_gamma2_minus(*sp);
   

    _sse_su3_multiply((*up));

    _sse_24_2_gamma2_minus();
        

/******************************* direction -2 *********************************/

    iy=iup(ix,3);

    sp=&psi[iy];
    _prefetch_spinor(sp);

    _sse_42_1_gamma2_plus(*sm);
      
      
    _sse_su3_inverse_multiply((*um));
      
    _sse_24_1_gamma2_plus();

      
    up+=1;
    _prefetch_su3(up);

    _sse_42_2_gamma2_plus(*sm);

      
      
    _sse_su3_inverse_multiply((*um));

    _sse_24_2_gamma2_plus();

        
      
/******************************* direction +3 *********************************/

    iy=idn(ix,3);

    sm=&psi[iy];
    _prefetch_spinor(sm);

    _sse_42_1_gamma3_minus(*sp);
      

    _sse_su3_multiply((*up));

    _sse_24_1_gamma3_minus_add();
      
      
    um=&_gauge_field0_3(iy);
    _prefetch_su3(um);

    _sse_42_2_gamma3_minus(*sp);
     

    _sse_su3_multiply((*up));

    _sse_24_2_gamma3_minus_add();
     
      
/******************************* direction -3 *********************************/

    iz=ix+1;
    if (iz==high)
      iz=0;

    iy=iup(iz,0);
      
    sp=&psi[iy];
    _prefetch_spinor(sp);

    _sse_42_1_gamma3_plus(*sm);
      
      
    _sse_su3_inverse_multiply((*um));

    rn=&res[ix];
      
    _sse_24_1_gamma3_plus_set();

    up=&_gauge_field0_0(iz);
    _prefetch_su3(up);

    _sse_42_2_gamma3_plus(*sm);
      
      
    _sse_su3_inverse_multiply((*um));

    _sse_24_2_gamma3_plus_set();

      
      
    /* printf("low:%i high:%i ix:%i id:%i, ptr:%x\n", low, high,ix, id, ptr);*/
/******************************** end of loop *********************************/

  }
}



void D_psi_fun_minus(size_t lo, size_t hi, int id, const void *ptr )
{
  int ix,iy,iz;
  const Arg_s *a = (const Arg_s*)ptr; 
  static int ix1,iy1,iy2,iz1;
  const int low = lo; 
  const int high = hi;
  MY_GAUGE (*gauge_field)[4] ALIGN = a->u;
  MY_SPINOR *psi = a->psi;
  MY_SPINOR *res = a->res;
  MY_GAUGE *up,*um;
  MY_SPINOR *s,*sp,*sm,*rn;            
  MY_SPINOR rs ALIGN;
  sse_int _sse_sgn ALIGN ={0x0,0x80000000,0x0,0x0};
  sse_double _minus_one ALIGN = {-1.0, -1.0};
  sse_double _conj ALIGN = {1.0, -1.0};

  /*printf("low:%i high:%i id:%i, ptr:%x\n", low, high, id, ptr);*/

  iy=iup(low,0);
  sp=&psi[iy];
  up=&_gauge_field0_0(low);
   
/************************ loop over all lattice sites *************************/
   
  for (ix=low;ix<high;ix++) 
  {
    /*s=&psi[ix];
      _prefetch_spinor(s);*/
    /*printf("low:%i high:%i ix:%i id:%i, ptr:%x\n", low, high,ix, id, ptr); */
/******************************* direction +0 *********************************/

    /*1 - gamma(0) */
	  
    iy=idn(ix,0);
      
    sm=&psi[iy];
    _prefetch_spinor(sm);  
	  
    _sse_42_1_gamma0_plus(*sp);
      
    _sse_su3_multiply((*up));

    _sse_24_1_gamma0_plus_set();
      
     
      
    um=&_gauge_field0_0(iy);
    _prefetch_su3(um);
      
    _sse_42_2_gamma0_plus(*sp);
      
      
    _sse_su3_multiply((*up));
      
    _sse_24_2_gamma0_plus_set();

    /*printf("\nrs[2][0][0]:%f\n",(rs_c2__)[0][0]);
      printf("\nrs[3][0][1]:%f\n",(rs_c3__)[0][1]);*/
	  
    /*printf("\nrs[2][0][1]:%f\n",(rs_c2__)[0][1]);
      printf("\nrs[3][0][0]:%f\n",(rs_c3__)[0][0]);*/
      

/******************************* direction -0 *********************************/

    /*1 + gamma(0) */
    iy=iup(ix,1);

    sp=&psi[iy];
    _prefetch_spinor(sp);

    _sse_42_1_gamma0_minus(*sm);

      
    _sse_su3_inverse_multiply((*um));
      
    _sse_24_1_gamma0_minus_add();

      
      
    up+=1;
    _prefetch_su3(up);
      
    _sse_42_2_gamma0_minus(*sm);
     
      
    _sse_su3_inverse_multiply((*um));
      
    _sse_24_2_gamma0_minus_add();
      
/******************************* direction +1 *********************************/

    iy=idn(ix,1);
      
    sm=&psi[iy];
    _prefetch_spinor(sm);

    _sse_42_1_gamma1_plus(*sp);

    _sse_su3_multiply((*up));

    _sse_24_1_gamma1_plus();


    um=&_gauge_field0_1(iy);
    _prefetch_su3(um);

    _sse_42_2_gamma1_plus(*sp);

    _sse_su3_multiply((*up));

    _sse_24_2_gamma1_plus();
          

/******************************* direction -1 *********************************/

    iy=iup(ix,2);

    sp=&psi[iy];
    _prefetch_spinor(sp);

    _sse_42_1_gamma1_minus(*sm);
      
      
    _sse_su3_inverse_multiply((*um));
      
    _sse_24_1_gamma1_minus();

     

    up+=1;
    _prefetch_su3(up);

    _sse_42_2_gamma1_minus(*sm);
      
      
    _sse_su3_inverse_multiply((*um));
      
    _sse_24_2_gamma1_minus();
     

/******************************* direction +2 *********************************/

    iy=idn(ix,2);

    sm=&psi[iy];
    _prefetch_spinor(sm);

    _sse_42_1_gamma2_plus(*sp);


    _sse_su3_multiply((*up));

    _sse_24_1_gamma2_plus();

      
    um=&_gauge_field0_2(iy);
    _prefetch_su3(um);

    _sse_42_2_gamma2_plus(*sp);
   

    _sse_su3_multiply((*up));

    _sse_24_2_gamma2_plus();
        

/******************************* direction -2 *********************************/

    iy=iup(ix,3);

    sp=&psi[iy];
    _prefetch_spinor(sp);

    _sse_42_1_gamma2_minus(*sm);
      
      
    _sse_su3_inverse_multiply((*um));
      
    _sse_24_1_gamma2_minus();

      
    up+=1;
    _prefetch_su3(up);

    _sse_42_2_gamma2_minus(*sm);

      
      
    _sse_su3_inverse_multiply((*um));

    _sse_24_2_gamma2_minus();

        
      
/******************************* direction +3 *********************************/

    iy=idn(ix,3);

    sm=&psi[iy];
    _prefetch_spinor(sm);

    _sse_42_1_gamma3_plus(*sp);
      

    _sse_su3_multiply((*up));

    _sse_24_1_gamma3_plus_add();
      
      
    um=&_gauge_field0_3(iy);
    _prefetch_su3(um);

    _sse_42_2_gamma3_plus(*sp);
     

    _sse_su3_multiply((*up));

    _sse_24_2_gamma3_plus_add();
     
      
/******************************* direction -3 *********************************/

    iz=ix+1;
    if (iz==high)
      iz=0;

    iy=iup(iz,0);
      
    sp=&psi[iy];
    _prefetch_spinor(sp);

    _sse_42_1_gamma3_minus(*sm);
      
      
    _sse_su3_inverse_multiply((*um));

    rn=&res[ix];
      
    _sse_24_1_gamma3_minus_set();

    up=&_gauge_field0_0(iz);
    _prefetch_su3(up);

    _sse_42_2_gamma3_minus(*sm);
      
      
    _sse_su3_inverse_multiply((*um));

    _sse_24_2_gamma3_minus_set();

      
      
    /* printf("low:%i high:%i ix:%i id:%i, ptr:%x\n", low, high,ix, id, ptr);*/
/******************************** end of loop *********************************/

  }
}




#elif defined SSE

#warning "missing that section - I think it is not SSE2"

void tnxtsu3dslash(float *u, float *psi, float *res, int isign, int cb)
{
  fprintf(stderr,"tnxtsu3dslash not implemented on this platform. Check compilation flags\n");
}


#else

#warning "missing that section - I think it is not SSE or SSE2"

void tnxtsu3dslash(float *u, float *psi, float *res, int isign, int cb)
{
  fprintf(stderr,"tnxtsu3dslash not implemented on this platform. Check compilation flags\n");
}


#endif



#ifdef __cplusplus
}
#endif



