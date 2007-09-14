/*******************************************************************************
 * $Id: sse_su3dslash_64bit_scalar.c,v 1.3 2007-09-14 19:32:11 bjoo Exp $
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


#include <stdlib.h>
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


#include <sse64.h>
#include <sse_align.h>
#include <shift_tables_scalar.h>


  
  static int initP = 0;          /* Reference count */
  static int *shift_table;       /* Shift Table */

  /** ---------------- This may need to change somehow for noncontiguous subsets ----------- **/
  static int icolor_start[2];    /* starting site for each coloring (cb) */
  static int icolor_end[2];      /* end site for each coloring (cb) */


  /* now some typedefs */
  typedef double chi_double[2] __attribute__ ((aligned (16)));   /* chi_double: Two doubles (complex?) */
  typedef chi_double chi_two[2] __attribute__ ((aligned (16)));  /* chi_two   : Four doubles */
  typedef double u_mat_array[3][3][2]  ALIGN;                    /* color color re/im */ 
  typedef double spinor_array[4][3][2] ALIGN;                    /* Nspin4 color re/im */
  typedef chi_two chi_array[3]    ALIGN;                         /* Halfspinor: re/im ::note:: color slowest varying */
  typedef u_mat_array (*my_mat_array)[4] ALIGN;  



  /* Gamma Stuff: Here we use _sse_42_1_gammaX_plusminus to apply ( 1 +/- gamma_X ) and keep the top half of 
                                                                                      halfspinor

     and                      _sse_42_2_gammaX_plusminus to apply ( 1 +/- gamma_x ) and keep the bottom half of
                                                                                      halfspinor

     Likewise the 24 variants do the reconstructs.
  */

  /* gamma 0 */

  /* Project */
#define _sse_42_1_gamma0_minus(sp) \
      _sse_load((sp)[0]); \
      _sse_load_up((sp)[3]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

  /* Recons with Store */
#define _sse_24_1_gamma0_minus_set() \
      _sse_store_up(rs[0]);\
      _sse_vector_i_mul_up();\
      _sse_store_up(rs[3]) 

  /* Recons with Accumulate */
#define _sse_24_1_gamma0_minus_add() \
      _sse_load(rs[0]);\
      _sse_vector_add();\
      _sse_store(rs[0]);\
      _sse_load(rs[3]);\
      _sse_vector_i_mul();\
      _sse_vector_add();\
      _sse_store(rs[3]) 
	
  /* Project 2 */
#define _sse_42_2_gamma0_minus(sp) \
      _sse_load((sp)[1]);\
      _sse_load_up((sp)[2]);\
      _sse_vector_i_mul();\
      _sse_vector_sub()

  /* Reconst store */
#define _sse_24_2_gamma0_minus_set() \
	  _sse_store_up(rs[1]);\
      _sse_vector_i_mul_up();\
      _sse_store_up(rs[2])

  /* Recons add ... etc */
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


  
/* end spin basis section */

static int init=0;
static sse_double fact1,fact2;
static spinor_array rs __attribute__ ((aligned (16)));



/***************** start of initialization routine ***************************************/
void init_sse_su3dslash(const int latt_size[])
{
  int mu, vol_cb;
  const int Nd=4;

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

  /* If make_shift_tables cannot allocate, it will barf */
  shift_table = make_shift_tables(icolor_start, latt_size);
  vol_cb = getTotalVolCB();

  /* Check that TotalVolCB has been set. If not it will be minus 1 */
  if ( vol_cb <= 0 ) { 
    fprintf(stderr, "Something is wrong... vol_cb =%d\n", vol_cb);
    exit(1);
  }

  /* This will need to change somehow */
  icolor_end[0] = icolor_start[0] + vol_cb;
  icolor_end[1] = icolor_start[1] + vol_cb;

  initP = 1;

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
  if( initP == 0 ) {
    free_shift_tables(&shift_table);
  }

}

/***************** end of initialization routine ***************************************/


/* prototypes - for Thread Slaves */

void D_psi_fun_plus(size_t lo,size_t hi, int id, const void *ptr);
void D_psi_fun_minus(size_t lo,size_t hi, int id, const void *ptr);

typedef struct {
  spinor_array *psi;  /* input spinor */
  spinor_array *res;   /* output spinor */
  u_mat_array (*u)[4];  /* gauge field associated with the output checkerboard */
  u_mat_array (*u2)[4];  /* gauge field associated with the input checkerboard */
  int cb;             /* output checkerboard */
  /* the output array                         */
} Arg_s;


/* Stripped out the sml_scall from this version -- non-threaded for now */
/* Hence n = 0. => lo = 0, hi = volume2, */

#define smpscaller2(a,bleah,spinfun2,chifun2,u3,cb2) \
    a.psi = spinfun2;\
    a.res = chifun2;\
    a.u = u3;  \
    a.cb = cb2; \
    (*bleah)(icolor_start[cb2], icolor_end[cb2], 0, &a);


/* External routine */
void sse_su3dslash_wilson(double *u, double *psi, double *res, int isign, int cb)
{
  Arg_s a;

  if (isign == 1)  {
    smpscaller2(a, D_psi_fun_plus, 
		(spinor_array*)psi,
		(spinor_array*)res,
		(my_mat_array)u,
		1-cb);
  }

  if( isign == -1) {
    smpscaller2(a, D_psi_fun_minus, 
		(spinor_array*)psi,
		(spinor_array*)res,
		(my_mat_array)u,
		1-cb);
  }
}

 


void D_psi_fun_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  int ix,iy,iz;                                      /* Ix is corrent site */
                                                     /* Iy is a neighbour */
                                                     /* Iz is next site in loop */

  const Arg_s *a = (const Arg_s*)ptr;                /* Downcast argument */
  const int low = lo;                                /* Start site */
  const int high = hi;                               /* End site + 1 */

  u_mat_array (*gauge_field)[4] ALIGN = a->u;        /* My gauge field */
  spinor_array *psi = a->psi;                        /* Source */
  spinor_array *res = a->res;                        /* Result */

  u_mat_array *up,*um;                               /* Pointer to FORWARD neighbour */
  spinor_array *s,*sp,*sm,*rn;                       /* Pointer to BACKWARD neighbour */


  /* This is like a prefetch 
     - we peel it off the loop */

  /* Get 4 spinor from forward direction */
  sp=&psi[forward_neighbor(shift_table,low,0) ];

  /* Get Gauge Field */
  up=&(gauge_field[low][0]);
   
  /************************ loop over all lattice sites *************************/
  for (ix=low;ix<high;ix++) 
  {

    /******************************* direction +0 *********************************/

    /* Prefetch back spinor for next dir: -0 */
    iy=backward_neighbor(shift_table,ix,0);
    sm=&psi[iy];
    _prefetch_spinor(sm);  

    /* Project Top 2 components */
    _sse_42_1_gamma0_minus(*sp);

    /* multiplication with SU(3) * one component of a halfspinor*/
    _sse_su3_multiply((*up));

    /* Spin recons Top 2 components */
    _sse_24_1_gamma0_minus_set();
      
     
    /* Prefetch gauge field for next dir: 0- */
    um=&(gauge_field[iy][0]);
    _prefetch_su3(um);


    /* Project and get Bottom 2 components */
    _sse_42_2_gamma0_minus(*sp);
      
    /* multiply */
    _sse_su3_multiply((*up));

    /* Reconstruct and set bottom 2 components */
    _sse_24_2_gamma0_minus_set();
      

    /******************************* direction -0 *********************************/
    /*1 + gamma(0) */
    /* sm and um should already be prefetched */

    /* Now prefetch forward neighbor for next dir (1+) */
    sp=&psi[ forward_neighbor(shift_table,ix,1) ];
    _prefetch_spinor(sp);

    /* Project 1 */
    _sse_42_1_gamma0_plus(*sm);

    /* Multiply */
    _sse_su3_inverse_multiply((*um));
      
    /* Recons accumulate */
    _sse_24_1_gamma0_plus_add();

      
    /* prefetch Next gauge field direction.*/
    up =&(gauge_field[ix][1]);
    _prefetch_su3(up);
      
    /* 2nd Component Project */
    _sse_42_2_gamma0_plus(*sm);
     
    /* Multiply adj */
    _sse_su3_inverse_multiply((*um));

    /* Recons/Accumulate */
    _sse_24_2_gamma0_plus_add();
      

    /******************************* direction +1 *********************************/
    /* up and sp should be prefetched */

    /* Prefetch spinor for next dir: 1- */
    iy=backward_neighbor(shift_table,ix,1);
    sm=&psi[iy];
    _prefetch_spinor(sm);

    /* 1st component: Project */
    _sse_42_1_gamma1_minus(*sp);

    /* Multiply */
    _sse_su3_multiply((*up));

    /* Recons */
    _sse_24_1_gamma1_minus();


    /* Prefetch gauge field for next dir: 1- */
    um=&(gauge_field[iy][1]);
    _prefetch_su3(um);

    /* 2nd component: Project */
    _sse_42_2_gamma1_minus(*sp);

    /* Multiply */
    _sse_su3_multiply((*up));

    /* Recon */
    _sse_24_2_gamma1_minus();


    /******************************* direction -1 *********************************/
    iy=forward_neighbor(shift_table,ix,2);
    sp=&psi[iy];
    _prefetch_spinor(sp);


    /* 1st component */
    /* Project */
    _sse_42_1_gamma1_plus(*sm);

    /* Multiply */
    _sse_su3_inverse_multiply((*um));

    /* Recon */
    _sse_24_1_gamma1_plus();

     
    /* Prefetch Next Gauge Field: 2+ */
    up = &(gauge_field[ix][2]);
    _prefetch_su3(up);


    /* 2nd component */
    /* Project */
    _sse_42_2_gamma1_plus(*sm);

    /* Multiply */
    _sse_su3_inverse_multiply((*um));

    /* Recon */
    _sse_24_2_gamma1_plus();
     

    /******************************* direction +2 *********************************/

    /* Prefetch spinor for nex direction: Dir 2- */
    iy=backward_neighbor(shift_table,ix,2);
    sm=&psi[iy];
    _prefetch_spinor(sm);

    /* 1st component: Project */
    _sse_42_1_gamma2_minus(*sp);

    /* Multiply */
    _sse_su3_multiply((*up));

    /* Recon */
    _sse_24_1_gamma2_minus();

    /* Prefetch gauge field for next direction */
    um=&(gauge_field[iy][2]);
    _prefetch_su3(um);

    /* 2nd component: Project */
    _sse_42_2_gamma2_minus(*sp);
   
    /* Multiply */
    _sse_su3_multiply((*up));

    /* Recon */
    _sse_24_2_gamma2_minus();
        

    /******************************* direction -2 *********************************/

    /* Prefetch next direction: 3+ */
    iy=forward_neighbor(shift_table,ix,3);
    sp=&psi[iy];
    _prefetch_spinor(sp);

    /* 1st component: Project */
    _sse_42_1_gamma2_plus(*sm);
      
    /* Multiply */
    _sse_su3_inverse_multiply((*um));
      
    /* Recon */
    _sse_24_1_gamma2_plus();


    /* Prefeth gauge field for next dir: 3+ */
    up = &(gauge_field[ix][3]);
    _prefetch_su3(up);


    /* 2nd component Project */
    _sse_42_2_gamma2_plus(*sm);

    /* Multiply */
    _sse_su3_inverse_multiply((*um));

    /* Recon */
    _sse_24_2_gamma2_plus();

        
      
    /******************************* direction +3 *********************************/

    /* Prefetch spinor for next direction: 3- */
    iy=backward_neighbor(shift_table,ix,3);
    sm=&psi[iy];
    _prefetch_spinor(sm);


    /* First component: Project */
    _sse_42_1_gamma3_minus(*sp);
      

    /* Multiply */
    _sse_su3_multiply((*up));

    /* Recons accumulate */
    _sse_24_1_gamma3_minus_add();
      
      
    /* Prefetch gauge field for next direction: 3- */
    um=&(gauge_field[iy][3]);
    _prefetch_su3(um);

    /* Second Component: Project */
    _sse_42_2_gamma3_minus(*sp);
     

    /* Multiply */
    _sse_su3_multiply((*up));

    /* Reconstruct add */
    _sse_24_2_gamma3_minus_add();
     
      
    /******************************* direction -3 *********************************/

    /* Next site for prefetching purposes */
    iz=ix+1;

    if (iz == high) { /* If we're on the last site, prefetch first site to avoid */
      iz=0;           /* Running beyond array bounds */
    }


    /* Prefetch next spinor */
    iy=forward_neighbor(shift_table,iz,0);
    sp=&psi[iy];
    _prefetch_spinor(sp);


    /* First component: Project */
    _sse_42_1_gamma3_plus(*sm);
      
      
    /* Multiply */
    _sse_su3_inverse_multiply((*um));

    /* Store! */
    rn=&res[ix];
    _sse_24_1_gamma3_plus_set();      

    /* Prefetch gauge field for next iteration */
    up=&(gauge_field[iz][0]);
    _prefetch_su3(up);

    /* Second component */
    /* Project */
    _sse_42_2_gamma3_plus(*sm);
      
    /* Multiply */
    _sse_su3_inverse_multiply((*um));

    /* Reconstruct */
    _sse_24_2_gamma3_plus_set();

  }
}



void D_psi_fun_minus(size_t lo, size_t hi, int id, const void *ptr )
{
  int ix,iy,iz;                          /* ix is the current site */
                                         /* iy is the neighbour for prefetching */
                                         /* iz is the prefetch site for the 
					    next loop iteration */

  const Arg_s *a = (const Arg_s*)ptr;    /* Cast the void args pointer */

  const int low = lo;                     /* First site */
  const int high = hi;                    /* Last site+1 */

  u_mat_array (*gauge_field)[4] ALIGN = a->u; /* Gauge field */
  spinor_array *psi = a->psi;                 /* Source spinor */
  spinor_array *res = a->res;                 /* Result spinor */
  u_mat_array *up,*um;                        /* us for multiply (PLUS/MINUS) */
  spinor_array *sp,*sm,*rn;                   /* spinor pointers sp sm are the 
						 neighbours, rn is the result */

  /* 'peel this off' to allow prefetching */
  iy=forward_neighbor(shift_table,low,0);
  sp=&psi[iy];
  up=&(gauge_field[low][0]);
   
/************************ loop over all lattice sites *************************/
   
  for (ix=low;ix<high;ix++) {

    /*1 - gamma(0) */

    /* Prefetch spinor for next dir: 0-  */
    iy=backward_neighbor(shift_table,ix,0);
    sm=&psi[iy];
    _prefetch_spinor(sm);  


    /* First component: Project */
    _sse_42_1_gamma0_plus(*sp);
      
    /* Multiply */
    _sse_su3_multiply((*up));

    /* Reconstruct */
    _sse_24_1_gamma0_plus_set();
      
    
    /* Prefetch gauge field for next direction: 0- */
    um=&(gauge_field[iy][0]);
    _prefetch_su3(um);
      

    /* Second component: Project */
    _sse_42_2_gamma0_plus(*sp);
      
    /* Multiply */
    _sse_su3_multiply((*up));
      
    /* Reconstruct */
    _sse_24_2_gamma0_plus_set();


    /******************************* direction -0 *********************************/

    /*1 + gamma(0) */
    
    /* Prefetch spinor for next part: 0- */
    iy=forward_neighbor(shift_table,ix,1);
    sp=&psi[iy];
    _prefetch_spinor(sp);

    /* First component: Project */
    _sse_42_1_gamma0_minus(*sm);

    /* Multiply */
    _sse_su3_inverse_multiply((*um));
      
    /* Reconstruct */
    _sse_24_1_gamma0_minus_add();

      
    /* Prefetch gauge field for next direction: 1+ */
    up = &(gauge_field[ix][1]);
    _prefetch_su3(up);
      
    /* Second component: Project */
    _sse_42_2_gamma0_minus(*sm);
     
    /* Multiply */
    _sse_su3_inverse_multiply((*um));
      
    /* Recon */
    _sse_24_2_gamma0_minus_add();
    
    /******************************* direction +1 *********************************/

    /* Prefetch spinor for next direction: 1 - */
    iy=backward_neighbor(shift_table,ix,1);
    sm=&psi[iy];
    _prefetch_spinor(sm);


    /* First Component: Project */
    _sse_42_1_gamma1_plus(*sp);

    /* Multiply */
    _sse_su3_multiply((*up));

    /* Recon */
    _sse_24_1_gamma1_plus();

    /* Prefetch gauge field for next direction: 1- */
    um=&(gauge_field[iy][1]);
    _prefetch_su3(um);

    /* Second Component: Project */
    _sse_42_2_gamma1_plus(*sp);

    /* Multiply */
    _sse_su3_multiply((*up));

    /* Reconstruct */
    _sse_24_2_gamma1_plus();
          

    /******************************* direction -1 *********************************/
    
    /* Prefetch for next direction: 2+ */
    iy=forward_neighbor(shift_table,ix,2);
    sp=&psi[iy];
    _prefetch_spinor(sp);


    /* 1st component: Project */
    _sse_42_1_gamma1_minus(*sm);

    /* Multiply */
    _sse_su3_inverse_multiply((*um));
      
    /* Reconstruct */
    _sse_24_1_gamma1_minus();

     
    /* Prefetch gauge field for next direction: 2- */
    up =&(gauge_field[ix][2]);
    _prefetch_su3(up);

    /* 2nd component: Project */
    _sse_42_2_gamma1_minus(*sm);
      
    /* Multiply */
    _sse_su3_inverse_multiply((*um));
      
    /* Reconstruct */
    _sse_24_2_gamma1_minus();
     
    
    /******************************* direction +2 *********************************/
    /* Prefetch for next direction: 2- */

    iy=backward_neighbor(shift_table,ix,2);
    sm=&psi[iy];
    _prefetch_spinor(sm);

    /* First component: Project */
    _sse_42_1_gamma2_plus(*sp);


    /* Multiply */
    _sse_su3_multiply((*up));

    /* Reconstruct */
    _sse_24_1_gamma2_plus();

      
    /* Prefetch gauge field for next direction: 2- */
    um=&(gauge_field[iy][2]);
    _prefetch_su3(um);

    /* Second component: Project */
    _sse_42_2_gamma2_plus(*sp);
   
    /* Multiply */
    _sse_su3_multiply((*up));

    /* Reconstruct */
    _sse_24_2_gamma2_plus();
        
    
    /******************************* direction -2 *********************************/
    /* Prefetch spinor for next direction: 3+ */
    iy=forward_neighbor(shift_table,ix,3);
    sp=&psi[iy];
    _prefetch_spinor(sp);


    /* First Half: Project */
    _sse_42_1_gamma2_minus(*sm);
      
    /* Multiply */
    _sse_su3_inverse_multiply((*um));
      
    /* Reconstruct */
    _sse_24_1_gamma2_minus();

    /* Prefetch gauge field for nex dir: 3+ */
    up = &(gauge_field[ix][3]);
    _prefetch_su3(up);

    /* Second Component: Project */
    _sse_42_2_gamma2_minus(*sm);

    /* Multiply */
    _sse_su3_inverse_multiply((*um));

    /* Reconstruct */
    _sse_24_2_gamma2_minus();

        
      
    /******************************* direction +3 *********************************/
    /* Prefetch spinor for the next direction: 3- */

    iy=backward_neighbor(shift_table,ix,3);
    sm=&psi[iy];
    _prefetch_spinor(sm);

    /* 1st component: Project */
    _sse_42_1_gamma3_plus(*sp);
      

    /* Multiply */
    _sse_su3_multiply((*up));

    /* Reconstruct */
    _sse_24_1_gamma3_plus_add();
      
      
    /* Prefetch gauge field for next direction: 3- */
    um=&(gauge_field[iy][3]);
    _prefetch_su3(um);

    /* 2nd component: Project */
    _sse_42_2_gamma3_plus(*sp);
     

    /* Multiply */
    _sse_su3_multiply((*up));

    /* Reconstruct add */
    _sse_24_2_gamma3_plus_add();
     
      
    /******************************* direction -3 *********************************/

    /* Next site in loop. We peeled this off the loop to start with so we can prefetch... */
    iz=ix+1; 
    if (iz==high) { /* If we are on the last site, we should prefetch the first element, to
		       avoid running past the array bounds */
      iz=0;
    }

    /* Prefetch the spinor for next site, dir 0+ */
    iy=forward_neighbor(shift_table,iz,0);
    sp=&psi[iy];
    _prefetch_spinor(sp);


    /* First component: Project */
    _sse_42_1_gamma3_minus(*sm);
      
    /* Multiply */
    _sse_su3_inverse_multiply((*um));

    /* Set pointer of the target spinor site */
    rn=&res[ix];
      
    /* Recons: Store */
    _sse_24_1_gamma3_minus_set();

    /* Prefetch the gauge field for next site dir  0+ */
    up=&(gauge_field[iz][0]);
    _prefetch_su3(up);

    /* Second Component: Project */
    _sse_42_2_gamma3_minus(*sm);
      
    /* Multiply */
    _sse_su3_inverse_multiply((*um));

    /* Reconstruct and store */
    _sse_24_2_gamma3_minus_set();
  }
}


#ifdef __cplusplus
}
#endif



