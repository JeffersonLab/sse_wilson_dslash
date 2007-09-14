/*******************************************************************************
 * $Id: sse_su3dslash_32bit_scalar.c,v 1.3 2007-09-14 19:32:11 bjoo Exp $
 * 
 * Action of the 32bit single-node Wilson-Dirac operator D_w on a given spinor field
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
 * Date: 9/11/2001
 *
 *******************************************************************************/

#include <sse_config.h>
#include <sse_align.h>
#include <shift_tables_scalar.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/* externally callable function:  inxtsu3dslash */
/* requires gauge fields packed by no_funnystuff_pack_gauge_field of packer.c */

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



  /* Volume and initialization */
  static int initP = 0;

  /* These are needed for the shift tables */
  static int *shift_table;
  static int icolor_start[2];    /* starting site for each coloring (cb) */
  static int icolor_end[2];      /* end site for each coloring (cb) */


/* now overlays for spinors as arrays or structs */
  typedef float chi_float[2];                    /* 2 floats => chi float */
  typedef chi_float chi_two[2];                  /* 2 chi_floats = 4 floats => chi_two */
  typedef float u_mat_array[3][3][2]  ALIGN;     /* color color re/im */ 
  typedef float spinor_array[4][3][2] ALIGN;     /* Nspin4 color re/im */
  typedef chi_two chi_array[3]    ALIGN; /* half vector: ..color Nspin2 re/im ::note:: color slowest varying */
  typedef u_mat_array (*my_mat_array)[4] ALIGN;  /* an array of 4 u_mat_array pointers */

  
  /* macros for spin basis note: assume first two rows are linearly independent except for gamma3 */
  /* it should be possible to change the spin basis by just correctly modifying these macros. However,
     some non-chiral spin basis may need additional modifications inside the dslash routine*/


  /* use SZIN spin basis */
  
  /* gamma 0 */

#define _sse_42_gamma0_minus()   _sse_vector_xch_i_sub()

#define _sse_42_gamma0_plus()     _sse_vector_xch_i_add()

#define _sse_24_gamma0_minus_set()  _sse_vector_xch_i_mul_up()
 
#define _sse_24_gamma0_plus_set()   _sse_vector_xch_i_mul_neg_up()

#define _sse_24_gamma0_minus_add() _sse_vector_xch_i_add()

#define _sse_24_gamma0_plus_add() _sse_vector_xch_i_sub()


/* gamma 1 */

#define _sse_42_gamma1_minus()  \
								_sse_vector_xch();\
								_sse_vector_addsub()
#define _sse_42_gamma1_plus()  \
								_sse_vector_xch();\
								_sse_vector_subadd()
#define _sse_24_gamma1_minus()  \
								_sse_vector_xch();\
								_sse_vector_subadd()
#define _sse_24_gamma1_plus()  \
								_sse_vector_xch();\
								_sse_vector_addsub()



/* gamma 2 */

#define _sse_42_gamma2_minus()   _sse_vector_i_subadd()

#define _sse_42_gamma2_plus()     _sse_vector_i_addsub()

#define _sse_24_gamma2_minus() _sse_vector_i_addsub()

#define _sse_24_gamma2_plus() _sse_vector_i_subadd()

/* gamma 3 */


#define _sse_42_gamma3_minus()   _sse_vector_sub()

#define _sse_42_gamma3_plus()     _sse_vector_add()

#define _sse_24_gamma3_minus() _sse_vector_sub()

#define _sse_24_gamma3_plus() _sse_vector_add()

#define _sse_24_gamma3_minus_rows12() _sse_vector_add()

#define _sse_24_gamma3_plus_rows12() _sse_vector_add()




/***************** start of initialization routine ***************************************/

void init_sse_su3dslash(const int latt_size[])
{
  int mu;
  int vol_cb=-1;

  /* If we are already initialised, then increase the refcount and return */
  if (initP > 0) 
  {
    initP++;
    return;
  }

  
  /* Check problem and subgrid size */
  for(mu=0; mu < 4; mu++) {
    if ( latt_size[mu] % 2 != 0 ) {
      fprintf(stderr,"This SSE Dslash only supports even problem sizes. Here the lattice is odd in dimension %d with length %d\n", mu, latt_size[mu]);
      exit(1);
    }
  }

  /* Construct all the shift tables needed */
  /* 4 dimensions * 2 directions { aka FORWARD and BACKWARD } * volume */
    /* shift_table and icolor start are set, latt_size is read */
  shift_table = make_shift_tables(icolor_start, latt_size);

  vol_cb = getTotalVolCB();
  if ( vol_cb <= 0 ) { 
    fprintf(stderr,"Something is wrong. volcb = %d\n", vol_cb);
      exit(1);
  }
  
  /* Assume cb-s are contiguous */
  icolor_end[0] = icolor_start[0] + vol_cb;
  icolor_end[1] = icolor_start[1] + vol_cb;

  /* Set flag to signify initialization */
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

  /* If the refcount has now hit 0 then free the shift tables*/
  if( initP == 0 ) {
    free_shift_tables(&shift_table);
  }
}

/***************** end of initialization routine ***************************************/


/* Apply Dslash */
void D_psi_fun_plus(size_t lo,size_t hi, int id, const void *ptr);

/* Apply Dslash Dagger */
void D_psi_fun_minus(size_t lo,size_t hi, int id, const void *ptr);

/* Thread Argument */
typedef struct {
  spinor_array *psi;    /* input spinor */
  spinor_array *res;    /* output spinor */
  u_mat_array (*u)[4];    /* gauge field on the output checkerboard */
  u_mat_array (*u2)[4];   /* gauge field on the input checkerboard */
  int cb;            /* output checkerboard  */
} Arg_s;


/* Stripped out the sml_scall from this version -- non threaded */
/* Hence n = 0. => lo = 0, hi = volume2, */

#define smpscaller2(a,bleah,spinfun2,chifun2,u3,cb2) \
    a.psi = spinfun2;\
    a.res = chifun2;\
    a.u = u3;  \
    a.cb = cb2; \
    (*bleah)(icolor_start[cb2], icolor_end[cb2], 0, &a);


/* routine sse_su3dslash_wilson
   u: base pointer to gauge field
   psi: base pointer to input spinor field on full lattice
   res: base pointer to output spinor field on full lattice
   isign: 1-->normal, -1--> swaps 1 - gamma(mu^) for 1 + gamma(mu^)
   cb: checkerboard (0/1) of input fields
*/
void sse_su3dslash_wilson(float *u, float *psi, float *res, int isign, int cb)
{
  Arg_s a;

  if (isign == 1) {
    smpscaller2(a, D_psi_fun_plus, 
		(spinor_array*)psi,
		(spinor_array*)res,
		(my_mat_array)u,
		1-cb);
  }

  if( isign == -1) 
  {
    smpscaller2(a, D_psi_fun_minus, 
		(spinor_array*)psi,
		(spinor_array*)res,
		(my_mat_array)u,
		1-cb);
  }
}

  
#include <sse32.h>


void D_psi_fun_plus(size_t lo,size_t hi, int id, const void *ptr)
{

  const Arg_s *a  = (const Arg_s*)ptr;  /* Cast the (void *) to an (Arg_s*) */
  int ix1;                              /* Index of current site */
  int iy1,iy2;                          /* Index of neighbour of current site (iy1) 
					   and current site+1 (iy2) */

  int iz1;                              /* Index of next site in loop */
  const int low  =  lo;                 /* First site for this thread */
  const int high  =  hi;                /* Last site for this thread */

  u_mat_array (*gauge_field)[4]  =  a->u; /* Packed Gauge fields */
  spinor_array *psi  =  a->psi;           /* Source spinor */
  spinor_array *res  =  a->res;           /* Result spinor */


  /* Pointers to the neighboring u-s */
  u_mat_array *up1 ALIGN;                  /* U[ x  ] */
  u_mat_array *up2 ALIGN;                  /* U[ (x+1)  ] */
  u_mat_array *um1 ALIGN;                  /* U[ x - mu ] */
  u_mat_array *um2 ALIGN;                  /* U[ (x+1)-mu ] */

  /* 4 - Spinor pointers */
  spinor_array *sp1 ALIGN;
  spinor_array *sp2 ALIGN;
  spinor_array *sm1 ALIGN;
  spinor_array *sm2 ALIGN;
  spinor_array *sn1 ALIGN;

  /* Half Vectors */
  chi_array r12_1 ALIGN; /* Site 1 upper */
  chi_array r34_1 ALIGN; /* Site 1 lower */
  chi_array r12_2 ALIGN; /* Site 2 upper */
  chi_array r34_2 ALIGN; /* Site 2 lower */

  /* note these guys need to be declared in each routine in which they are used*/
  /* Is this still true? */
  /* These are used to flip the signs of various members of a 4-wide sse register */
  static sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};      /* Negate 1st and 2nd */
  static sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};      /* Negate 1st and 3rd */
  static sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};      /* Negate 1st and 4th */
  static sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};      /* Negate 2nd and 3rd */
  static sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};      /* Negate 2nd and 4th */
  static sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};      /* Negate 3rd and 4th */
  static sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f}; /* Negate All */
  
  /* note that we want the spinors from the checkerboard opposite the one we are writing to */
  /* We are doing things in bunches of two sites */

  /* Get forward neighbour of low in the x direction */
  sp1 = &psi[ forward_neighbor(shift_table,low,0)  ];

  /* Get forward neighbour of low + 1 in the x direction. */
  sp2 = &psi[ forward_neighbor(shift_table,low+1,0)];

  /* note: the gauge fields at x - mu^ are on the checkerboard opposite the one we are writing to */
  up1 = &(gauge_field[low][0]);
  up2 = &(gauge_field[low+1][0]);

  for (ix1 = low;ix1<high;ix1+= 2) 
  {
       
    /******************************* direction +0 *********************************/

    /* ...(1-isign*gamma(0)) psi(x + \hat{0}) */

    /* Prefetch the backward neighbours for the following 1 + isign gamma(0) case */
    iy1 =  backward_neighbor(shift_table,ix1,0);
    iy2 =  backward_neighbor(shift_table,ix1+1,0); 

    /* prefetch backward neighbor spinors  */	   
    sm1 = &psi[ iy1   ];
    _prefetch_single(sm1);
    sm2 = &psi[ iy2   ];
    _prefetch_single(sm2);

    /*loads in the appropriate spinor and does the spin projection */
    /* sp denotes spinors in the x+mu^ direction, sm denotes spinors in the x - mu^ direction */
    
    /* Do it for the first site */
    _sse_pair_load((*sp1)[0],(*sp1)[1]);
    _sse_pair_load_up((*sp1)[2],(*sp1)[3]);

    /* Spin Project */
    _sse_42_gamma0_minus();
    
    /* do the SU(3) * 3x2 multiply ....3 colors...both spin components of the spin-projected halfspinor*/ 
    _sse_su3_multiply((*up1));
       
    /* Reconstruct */

    /* the output is in xmm3-5, so we store_up...this is a _set since it is the first term in the sum over mu */
    /* top two spin components, the top two rows of the 4x2 reconstruction matrix will be the identity, so just store */
    _sse_vector_store_up(r12_1);
        
    /* bottom two spin components, use the bottom two components of the 4x2 reconstruction matrix ...whatever was done to the identity
       to get that submatrix needs to be done likewise to the halfspinor to reconstruct the bottom two components */     
    _sse_24_gamma0_minus_set();
    _sse_vector_store_up(r34_1);

    /* Done */

    /* Now prefetch for the  1 + \gamma_0 U^\dagger case */

    um1 = &(gauge_field[iy1][0]);
    _prefetch_single(um1);
    um2 = &(gauge_field[iy2][0]);
    _prefetch_single(um2);

    /* Now do the spin proj, multiply, spin recon for the second site */
    /* Load spinor */
    _sse_pair_load((*(sp2))[0],(*(sp2))[1]);
    _sse_pair_load_up((*(sp2))[2],(*(sp2))[3]);

    /* Project */
    _sse_42_gamma0_minus();

    /* Multiply */
    _sse_su3_multiply((*(up2)));

    /* Reconstruct */
    _sse_vector_store_up(r12_2);
    _sse_24_gamma0_minus_set(); 
    _sse_vector_store_up(r34_2);
      
    /******************************* direction -0 *********************************/
    /* Hopefully sm1 sm2 the backward neighbours are prefetched and in cache. The shifted 
       matrices should be prefetched too */
   

      /* ...(1+isign*gamma(0))... */
    
    /* Prefetch the forward spinors for the next direction */
    sp1 = &psi[ forward_neighbor(shift_table,ix1,1) ];
    _prefetch_single(sp1);
    sp2 = &psi[ forward_neighbor(shift_table,ix1+1,1)];
    _prefetch_single(sp2);


    /* Now load backward neighbours into SSE vectors */
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);

    /* Project */
    _sse_42_gamma0_plus();
	 
    /* Multiply with inverse */
    _sse_su3_inverse_multiply((*um1));


    /* Reconstruct accumulate */

    /*ok here's where things are different..now we load the partial sum over directions of the
      output spinor, then add on to them the spin reconstructed terms (using the bottom two rows of the 4x2 reconstruction
      matrix for the spin reconstruction), and then save the temp
      and hopefully it better stay in cache */
    
    /* notation: r12_1 -- first two spin components, first of the two  adjacent sites 
       r34_1 -- last two spin components, first site
       r12_2 -- first two spin components, adjacent (ix1+1) site
       r34_2 -- last two spin components, adjacent site
       
    */

    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);

    _sse_vector_load(r34_1);
    _sse_24_gamma0_plus_add();
    _sse_vector_store(r34_1);  

    /* Prefetch gauge field for next direction */
    up1 = &(gauge_field[ix1][1]);
    _prefetch_single(up1);      
    up2 = &(gauge_field[ix1+1][1]); 
    _prefetch_single(up2); 
    

    /* Now do second site */
    /* Load */
    _sse_pair_load((*(sm2))[0],(*(sm2))[1]);
    _sse_pair_load_up((*(sm2))[2],(*(sm2))[3]);

    /* Project */
    _sse_42_gamma0_plus();

    /* Multiply */
    _sse_su3_inverse_multiply((*(um2)));

    /* reconstruct accumulate */
    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);

    _sse_vector_load(r34_2);
    _sse_24_gamma0_plus_add();
    _sse_vector_store(r34_2);  
      
    /******************************* direction +1 *********************************/
    /* OK, sp1, sp2 and up1, up2 should be prefetched */

    /* Prefetch spinors for the -1 direction */
    iy1 = backward_neighbor(shift_table,ix1,1);
    iy2 = backward_neighbor(shift_table,ix1+1,1);
    sm1 = &psi[iy1];
    _prefetch_single(sm1);
    sm2 = &psi[iy2];
    _prefetch_single(sm2);

    /* First site */
    _sse_pair_load((*sp1)[0],(*sp1)[1]);
    _sse_pair_load_up((*sp1)[2],(*sp1)[3]);

    /* Project */
    _sse_42_gamma1_minus();

    /* Multiply */
    _sse_su3_multiply((*up1));

    /* Reconstruct accumulate */
    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);
    _sse_vector_load(r34_1);
    _sse_24_gamma1_minus();
    _sse_vector_store(r34_1);

    /* Prefetch gauge field for - direction */
    um1 = &(gauge_field[iy1][1]);
    _prefetch_single(um1);      
    um2 = &(gauge_field[iy2][1]);
    _prefetch_single(um2);

    /* Second site: */
    /* Load */
    _sse_pair_load((*(sp2))[0],(*(sp2))[1]);
    _sse_pair_load_up((*(sp2))[2],(*(sp2))[3]);

    /* Project */
    _sse_42_gamma1_minus();

    /* Multiply */
    _sse_su3_multiply((*(up2)));

    /* Reconstruct accumulate */
    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2); 
      
    _sse_vector_load(r34_2);
    _sse_24_gamma1_minus();
    _sse_vector_store(r34_2);
    
    /******************************* direction -1 *********************************/

    /* Prefetch forward neighbour for direction 2+ */
    sp1 = &psi[forward_neighbor(shift_table,ix1,2)];
    _prefetch_single(sp1);
    sp2 = &psi[forward_neighbor(shift_table,ix1+1,2)];
    _prefetch_single(sp2);

    /* Site 1: */
    /* Load */
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);

    /* Project */
    _sse_42_gamma1_plus();
	  
    /* Multiply */
    _sse_su3_inverse_multiply((*um1));

    /* Reconstruct Accumulate */
    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);      
    _sse_vector_load(r34_1);
    _sse_24_gamma1_plus();
    _sse_vector_store(r34_1);      

    /* Prefetch Gauge for next case: Direction 2 + */
    up1 = &(gauge_field[ix1][2]);
    _prefetch_single(up1);      
    up2 = &(gauge_field[ix1+1][2]);
    _prefetch_single(up2);

    /* Load */
    _sse_pair_load((*(sm2))[0],(*(sm2))[1]);
    _sse_pair_load_up((*(sm2))[2],(*(sm2))[3]);

    /* Project */
    _sse_42_gamma1_plus();

    /* Multiply */
    _sse_su3_inverse_multiply((*(um2)));

    /* Recons Accumulate */
    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);      

    _sse_vector_load(r34_2);
    _sse_24_gamma1_plus();
    _sse_vector_store(r34_2); 


    /******************************* direction +2 *********************************/
    /* sp1, sp2, up1, up2 should be prefetched and in cache */


    /* Prefetch sm1 & sm2 for -ve direction */
    iy1 = backward_neighbor(shift_table,ix1,2);
    sm1 = &psi[iy1];
    _prefetch_single(sm1);
    iy2 = backward_neighbor(shift_table,ix1+1,2);
    sm2 = &psi[iy2];
    _prefetch_single(sm2);

    /* Load */
    _sse_pair_load((*sp1)[0],(*sp1)[1]);
    _sse_pair_load_up((*sp1)[2],(*sp1)[3]);

    /* Project */
    _sse_42_gamma2_minus();

    /* Multiply */
    _sse_su3_multiply((*up1));

    /* Recons Accumulate */
    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);       
    _sse_vector_load(r34_1);
    _sse_24_gamma2_minus();
    _sse_vector_store(r34_1);       

    /* Prefetch gauge field for -ve case */
    um1 = &(gauge_field[iy1][2]);
    _prefetch_single(um1);
    um2 = &(gauge_field[iy2][2]);
    _prefetch_single(um2);

    /* Second site: load */
    _sse_pair_load((*(sp2))[0],(*(sp2))[1]);
    _sse_pair_load_up((*(sp2))[2],(*(sp2))[3]);

    /* Project */
    _sse_42_gamma2_minus();

    /* Multiply */
    _sse_su3_multiply((*(up2)));

    /* Recons/Accumulate */
    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);       
    _sse_vector_load(r34_2);
    _sse_24_gamma2_minus();
    _sse_vector_store(r34_2); 

    
    /******************************* direction -2 *********************************/
    /* sm1, sm2, um1, um2, should be prefetched and in cache                      */
    /******************************************************************************/
    
    /* Prefetch spinors for direction 3+ */
    sp1 = &psi[ forward_neighbor(shift_table,ix1,3) ];
    _prefetch_single(sp1);
    sp2 = &psi[ forward_neighbor(shift_table,ix1+1,3) ];
    _prefetch_single(sp2);
    
      
    /* First site: Load */
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);

    /* Project */
    _sse_42_gamma2_plus();      

    /* Multiply */
    _sse_su3_inverse_multiply((*um1));

    /* Recons/Accumulate */
    _sse_vector_load(r12_1);
    _sse_vector_add(); 
    _sse_vector_store(r12_1);
    _sse_vector_load(r34_1);
    _sse_24_gamma2_plus();
    _sse_vector_store(r34_1);

    /* Prefetch Gauge for the 3+ case */
    up1 = &(gauge_field[ix1][3]);
    _prefetch_single(up1);
    up2 = &(gauge_field[ix1+1][3]);
    _prefetch_single(up2);

    /* Load */
    _sse_pair_load((*(sm2))[0],(*(sm2))[1]);
    _sse_pair_load_up((*(sm2))[2],(*(sm2))[3]);

    /* Project */
    _sse_42_gamma2_plus();      

    /* Multiply */
    _sse_su3_inverse_multiply((*(um2)));

    /* Recons Accumulate */
    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);
    _sse_vector_load(r34_2);
    _sse_24_gamma2_plus();
    _sse_vector_store(r34_2);


    /******************************* direction +3 *********************************/
    /* sp1, sp2, up1, up2, should be prefetched and in cache already              */

    /* Pre Fetch neigbour spinors for the backwards case                          */
    iy1 = backward_neighbor(shift_table,ix1,3);
    iy2 = backward_neighbor(shift_table,ix1+1,3);
    sm1 = &psi[iy1]; 
    _prefetch_single(iy1); 
    sm2 = &psi[iy2]; 
    _prefetch_single(iy2) ;

    /* Load */
    _sse_pair_load((*sp1)[0],(*sp1)[1]);
    _sse_pair_load_up((*sp1)[2],(*sp1)[3]);

    /* Project */
    _sse_42_gamma3_minus();
      
    /* Multiply */
    _sse_su3_multiply((*up1));

    /* Recons/Accumulate */
    _sse_vector_load(r12_1);
    _sse_24_gamma3_minus_rows12();
    _sse_vector_store(r12_1);
    _sse_vector_load(r34_1);
    _sse_24_gamma3_minus();
    _sse_vector_store(r34_1);      

     
    /* Prefetch um for - case */
    um1 = &(gauge_field[iy1][3]); 
    _prefetch_single(um1);
    um2 = &(gauge_field[iy2][3]);
    _prefetch_single(um2);

    /* Site 2: Load */
    _sse_pair_load((*sp2)[0],(*sp2)[1]);
    _sse_pair_load_up((*sp2)[2],(*sp2)[3]);

    /* Project */
    _sse_42_gamma3_minus();
      
    /* Multiply */
    _sse_su3_multiply((*(up2)));

    /* Recons / Accumulate */
    _sse_vector_load(r12_2);
    _sse_24_gamma3_minus_rows12();
    _sse_vector_store(r12_2);
    _sse_vector_load(r34_2);
    _sse_24_gamma3_minus();
    _sse_vector_store(r34_2); 

    /******************************* direction -3 *********************************/
    /* ok, here's where things change up a bit...we have to set the output instead of saving it back to temp */
    /* also, the rows12 stuff is for spin basis that have the first two rows of the 4x2 matrix anything other
       than the identity...I only allow this for direction 3 because that's the only place I've seen it in the other
       spin basis I've worked with... */

    /* Iz is the next site (in step of 2) ie it is ix + 2 */
    iz1 = ix1+2;

    if (iz1 == high) {
      iz1 = 0;       /* If we reach high then set back to 0, so we don't access data from over 
			the array bound */
    }

    /* Prefetch forward neighbour spinor */
    sp1 = &psi[ forward_neighbor(shift_table,iz1,0) ];
    _prefetch_single(sp1);
    sp2 = &psi[ forward_neighbor(shift_table,iz1+1,0) ];
    _prefetch_single(sp2);

    /* Site1 : Load */
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);
    
    /* Project */
    _sse_42_gamma3_plus();
      
    /* Multiply */
    _sse_su3_inverse_multiply((*um1));

    /* Get address of result site */
    sn1 = &res[ix1];  /*we always walk across the result lexicographically */
       
    /* Recons 1st two rows */
    _sse_vector_load(r12_1);
    _sse_24_gamma3_plus_rows12();
     
    /* Store */
    _sse_pair_store((*sn1)[0],(*sn1)[1]);

    _sse_vector_load(r34_1);
    
    /* Recons 2nd two rows */
    _sse_24_gamma3_plus();
      
    /* Store */
    _sse_pair_store((*sn1)[2],(*sn1)[3]);      


    /* Prefetch gauge field for next loop iteration (0 direction) */
    up1 = &(gauge_field[iz1][0]);
    _prefetch_single(up1);
    up2 = &(gauge_field[iz1+1][0]);
    _prefetch_single(up2);


    /* Second site */

    /* Load */
    _sse_pair_load((*sm2)[0],(*sm2)[1]);
    _sse_pair_load_up((*sm2)[2],(*sm2)[3]);

    /* Recons */
    _sse_42_gamma3_plus();
      
    /* Multiply */
    _sse_su3_inverse_multiply((*um2));

    /* Recons row 12 */
    _sse_vector_load(r12_2);
    _sse_24_gamma3_plus_rows12();
      
    /* Store */
    _sse_pair_store((*(sn1+1))[0],(*(sn1+1))[1]);

    /* Recons row 34 */
    _sse_vector_load(r34_2);
    _sse_42_gamma3_plus();
     
    /* Store */
    _sse_pair_store((*(sn1+1))[2],(*(sn1+1))[3]); 
	  

    /******************************** end of loop *********************************/
      
  }
}

/*ok, now this routine here is just like the one above except isign has a different value, which means that the
signs used for 1 +- gamma(mu^) must be swapped...*/ 
void D_psi_fun_minus(size_t lo,size_t hi, int id, const void *ptr)
{
  const Arg_s *a  = (const Arg_s*)ptr;   /* Downcast to args */
  int ix1,iy1,iy2,iz1;                   /* Coordinates ix1 - current
					    iy1 - index of neighbour
					    iy1 - index of neighbour of ix1+1 
					    iz1 - index of first of next pair (ix+2) */

  const int low  =  lo;                    /* First site */
  const int high  =  hi;                   /* Last site+1 */
  u_mat_array (*gauge_field)[4]  =  a->u;  /* Gauge field */

  spinor_array *psi  =  a->psi;            /* Source 4-spinor */
  spinor_array *res  =  a->res;            /* Result 4-spinor */

  u_mat_array *up1 ALIGN;                  /* U[ ix ] */
  u_mat_array *up2 ALIGN;                  /* U[ ix+1 ] */
  u_mat_array *um1 ALIGN;                  /* U[ ix - mu ] */
  u_mat_array *um2 ALIGN;                  /* U[ ix+1 - mu ] */

  spinor_array *sp1 ALIGN;                 /* 4 spinor psi[ix1+mu] */
  spinor_array *sp2 ALIGN;                 /* 4 spinor psi[ix1+1 + mu] */
  spinor_array *sm1 ALIGN;                 /* 4 spinor psi[ix1-mu] */
  spinor_array *sm2 ALIGN;                 /* 4 spinor psi[ix1+1 - mu */
  spinor_array *sn1 ALIGN;                 /* 4 spinor result */

  chi_array r12_1;                         /* site 1 halfspinor top half */
  chi_array r34_1;                         /* site 1 halfspinor bottom half */
  chi_array r12_2;                         /* site 2 halfspinor top half */
  chi_array r34_2;                         /* site 2 halfspinor bottom half */
  
  static sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f}; 
  static sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f}; 
  static sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f}; 
  static sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f}; 
  static sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f}; 
  static sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f}; 
  static sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f};


  /* Pull these out of the loop kind a like a prefetch */
  iy1 = forward_neighbor(shift_table,low,0); 
  iy2 = forward_neighbor(shift_table,low+1,0) ;
  sp1 = &psi[iy1]; 
  sp2 = &psi[iy2]; 
  up1 = &(gauge_field[low][0]);
  up2 = &(gauge_field[low+1][0]); 

  /************************ loop over all lattice sites *************************/

  for (ix1 = low;ix1<high;ix1+= 2) {
    /******************************* direction +0 *********************************/
    
    /* ...(1-isign*gamma(0))... */

    /* Prefetch  spinor for 0- case */
    iy1 = backward_neighbor(shift_table,ix1,0);
    iy2 = backward_neighbor(shift_table,ix1+1,0);

    sm1 = &psi[iy1];
    _prefetch_single(sm1);
    sm2 = &psi[iy2];
    _prefetch_single(sm2);

    /* Site1: Load */
    _sse_pair_load((*sp1)[0],(*sp1)[1]);
    _sse_pair_load_up((*sp1)[2],(*sp1)[3]);

    /* Project */
    _sse_42_gamma0_plus();
    
    /* Multiply */
    _sse_su3_multiply((*up1));

    /* Recons */
    _sse_vector_store_up(r12_1);
    _sse_24_gamma0_plus_set();
    _sse_vector_store_up(r34_1);

    /* Prefetch gauge field for next part */
    um1 = &(gauge_field[iy1][0]);
    _prefetch_single(um1);
    um2 = &(gauge_field[iy2][0]);
    _prefetch_single(um2);

    /* Site 2: Load */
    _sse_pair_load((*(sp2))[0],(*(sp2))[1]);
    _sse_pair_load_up((*(sp2))[2],(*(sp2))[3]);

    /* Project */
    _sse_42_gamma0_plus();

    /* Multiply */
    _sse_su3_multiply((*(up2)));

    /* Recons */
    _sse_vector_store_up(r12_2);
    _sse_24_gamma0_plus_set(); 
    _sse_vector_store_up(r34_2);
      
    /******************************* direction -0 *********************************/
    /* sm1, sm2, um1, um2 should be prefetched */

    /* ...(1+isign*gamma(0))... */

    /* Prefetch for next part (1+) */
    iy1 = forward_neighbor(shift_table,ix1,1);
    iy2 = forward_neighbor(shift_table,ix1+1,1);
      sp1 = &psi[iy1];
    _prefetch_single(sp1);
    sp2 = &psi[iy2];
    _prefetch_single(sp2);

    /* Site 1: Load */
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);

    /* Project */
    _sse_42_gamma0_minus();
	 
    /* Multiply */
    _sse_su3_inverse_multiply((*um1));

    /* Recons/Accumulate */
    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);
    _sse_vector_load(r34_1);
    _sse_24_gamma0_minus_add();
    _sse_vector_store(r34_1);  


    /* Prefetch gauge field for next part: (1+) */
    up1 = &(gauge_field[ix1][1]);
    _prefetch_single(up1);      
    up2 = &(gauge_field[ix1+1][1]);
    _prefetch_single(up2);

    /* Site 2: Load */
    _sse_pair_load((*(sm2))[0],(*(sm2))[1]);
    _sse_pair_load_up((*(sm2))[2],(*(sm2))[3]);
    
    /* Project */
    _sse_42_gamma0_minus();

    /* Multiply */
    _sse_su3_inverse_multiply((*(um2)));

    /* Recons / Accumulate */
    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);

    _sse_vector_load(r34_2);
    _sse_24_gamma0_minus_add();
    _sse_vector_store(r34_2);  
      

    /******************************* direction +1 *********************************/
    /* sp1, sp2, up1, up2 should be prefetched                                    */

    /* Prefetch Spinors for next part (1-) */
    iy1 = backward_neighbor(shift_table,ix1,1);
    iy2 = backward_neighbor(shift_table,ix1+1,1);
    sm1 = &psi[iy1];
    _prefetch_single(sm1);
    sm2 = &psi[iy2];
    _prefetch_single(sm2);

    /* Site 1: Load */
    _sse_pair_load((*sp1)[0],(*sp1)[1]);
    _sse_pair_load_up((*sp1)[2],(*sp1)[3]);

    /* Project */
    _sse_42_gamma1_plus();

    /* Multiply */
    _sse_su3_multiply((*up1));

    /* Recons / Accumulate */
    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);

    _sse_vector_load(r34_1);
    _sse_24_gamma1_plus();
    _sse_vector_store(r34_1);

    /* Prefetch gauge links for next part (1-) */
    um1 = &(gauge_field[iy1][1]);
    _prefetch_single(um1);      
    um2 = &(gauge_field[iy2][1]);
    _prefetch_single(um2); 

      
    /* Site2: Load */
    _sse_pair_load((*(sp2))[0],(*(sp2))[1]);
    _sse_pair_load_up((*(sp2))[2],(*(sp2))[3]);

    /* Project */
    _sse_42_gamma1_plus();

    /* Multiply */
    _sse_su3_multiply((*(up2)));

    /* Recons Accumulate */
    _sse_vector_load(r12_2);
    _sse_vector_add(); 
    _sse_vector_store(r12_2);  

    _sse_vector_load(r34_2); 
    _sse_24_gamma1_plus(); 
    _sse_vector_store(r34_2); 


    /******************************* direction -1 *********************************/

    /* Prefetch spinors for next part: (2+) */
    iy1 = forward_neighbor(shift_table,ix1,2);
    iy2 = forward_neighbor(shift_table,ix1+1,2);
    sp1 = &psi[iy1];
    _prefetch_single(sp1);
    sp2 = &psi[iy2];
    _prefetch_single(sp2);

    /* Site 1:  Load */
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);

    /* Project */
    _sse_42_gamma1_minus();
	  
    /* Multiply */
    _sse_su3_inverse_multiply((*um1));

    /* Recons / Accumulate */
    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);      

    _sse_vector_load(r34_1);
    _sse_24_gamma1_minus();
    _sse_vector_store(r34_1);      

    /* Prefetch gauge field for next part (2+) */
    up1 = &(gauge_field[ix1][2]);
    _prefetch_single(up1);       
    up2 = &(gauge_field[ix1+1][2]); 
    _prefetch_single(up2); 

    /* Site 2: Load */
    _sse_pair_load((*(sm2))[0],(*(sm2))[1]);
    _sse_pair_load_up((*(sm2))[2],(*(sm2))[3]);

    /* Project */
    _sse_42_gamma1_minus();

    /* Multiply */
    _sse_su3_inverse_multiply((*(um2)));

    /* Recons Accumulate */
    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);      

    _sse_vector_load(r34_2);
    _sse_24_gamma1_minus();
    _sse_vector_store(r34_2); 

    /******************************* direction +2 *********************************/
    /* sp1, sp2, up1, up2 should  be in cache                                     */

    /* Prefetch spinors for next part: 2- */
    iy1 = backward_neighbor(shift_table,ix1,2);
    sm1 = &psi[iy1];
    _prefetch_single(sm1);
     iy2 = backward_neighbor(shift_table,ix1+1,2);
    sm2 = &psi[iy2];
    _prefetch_single(sm2);

    /* Site 1: Load */
    _sse_pair_load((*sp1)[0],(*sp1)[1]);
    _sse_pair_load_up((*sp1)[2],(*sp1)[3]);

    /* Project */
    _sse_42_gamma2_plus();

    /* Multiply */
    _sse_su3_multiply((*up1));

    /* Accumulate / recons */
    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);       

    _sse_vector_load(r34_1);
    _sse_24_gamma2_plus();
    _sse_vector_store(r34_1);       

    /* Prefetch Gauge field for next part: 2- */
    um1 = &(gauge_field[iy1][2]);
    _prefetch_single(um1);
    um2 = &(gauge_field[iy2][2]);
    _prefetch_single(um2);

    /* Site 2: Load */
    _sse_pair_load((*(sp2))[0],(*(sp2))[1]);
    _sse_pair_load_up((*(sp2))[2],(*(sp2))[3]);

    /* Project */
    _sse_42_gamma2_plus();

    /* Multiply */
    _sse_su3_multiply((*(up2)));

    /* Recons / Accumulate */
    _sse_vector_load(r12_2);
    _sse_vector_add(); 
    _sse_vector_store(r12_2);        

    _sse_vector_load(r34_2);
    _sse_24_gamma2_plus();
    _sse_vector_store(r34_2); 
 
    /******************************* direction -2 *********************************/

    /* sm1, sm2, um1, um2 should be cached */

    /* Prefetch spinor for next case: 3+ */
    iy1 = forward_neighbor(shift_table,ix1,3);
    sp1 = &psi[iy1];
    _prefetch_single(sp1);
    iy2 = forward_neighbor(shift_table,ix1+1,3); 
    sp2 = &psi[iy2];
    _prefetch_single(sp2);


    /* Site 1: Load */
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);

    /* Project */
    _sse_42_gamma2_minus();      

    /* Multiply */
    _sse_su3_inverse_multiply((*um1));

    /* Recons/Accumulate*/
    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);
      
    _sse_vector_load(r34_1);
    _sse_24_gamma2_minus();
    _sse_vector_store(r34_1);

    /* Prefetch gauge for next case: 3+ */
    up1 = &(gauge_field[ix1][3]);
    _prefetch_single(up1);
    up2 = &(gauge_field[ix1+1][3]);
    _prefetch_single(up2);
     
    /* Site 2: Load */
    _sse_pair_load((*(sm2))[0],(*(sm2))[1]);
    _sse_pair_load_up((*(sm2))[2],(*(sm2))[3]);
    
    /* Project */
    _sse_42_gamma2_minus();      

    /* Multiply */
    _sse_su3_inverse_multiply((*(um2)));

    /* Recons/ Accumulate */
    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);
      
    _sse_vector_load(r34_2);
    _sse_24_gamma2_minus();
    _sse_vector_store(r34_2); 

    
    /******************************* direction +3 *********************************/
    /* sp1, sp2, up1, up2 should be in cache */

    /* Prefetch spinors for next case: dir 3- */
    iy1 = backward_neighbor(shift_table,ix1,3);
    iy2 = backward_neighbor(shift_table,ix1+1,3);
    sm1 = &psi[iy1];
    _prefetch_single(iy1);
    sm2 = &psi[iy2];
    _prefetch_single(iy2);


    /* Site 1: Load */
    _sse_pair_load((*sp1)[0],(*sp1)[1]);
    _sse_pair_load_up((*sp1)[2],(*sp1)[3]);

    /* project */
    _sse_42_gamma3_plus();
      
    /* Multiply */
    _sse_su3_multiply((*up1));

    /* Recons Accumulate */
    _sse_vector_load(r12_1);
    _sse_24_gamma3_plus_rows12();
    _sse_vector_store(r12_1);

    _sse_vector_load(r34_1);
    _sse_24_gamma3_plus();
    _sse_vector_store(r34_1);      

    /* Prefetch gauge field for next part: 3- */
    um1 = &(gauge_field[iy1][3]); 
    _prefetch_single(um1);
    um2 = &(gauge_field[iy2][3]);
    _prefetch_single(um2);

    /* Site 2: Load */
    _sse_pair_load((*sp2)[0],(*sp2)[1]);
    _sse_pair_load_up((*sp2)[2],(*sp2)[3]);

    /* Project */
    _sse_42_gamma3_plus();
      
    /* Multiply */
    _sse_su3_multiply((*(up2)));

    /* Recons Accumulate */
    _sse_vector_load(r12_2);
    _sse_24_gamma3_plus_rows12();
    _sse_vector_store(r12_2);

    _sse_vector_load(r34_2);
    _sse_24_gamma3_plus();
    _sse_vector_store(r34_2); 

    /******************************* direction -3 *********************************/
    
    /* OK this was peeled off the loop, so we only want to do it if there are sites 
       left */

    /* First of the next pair of sites - what ix1 would be after loop iteration increment */
    iz1 = ix1+2;

    /* If we have reached the last site, set this to site zero. No harm in prefetching those. */
    if (iz1 == high)
      iz1 = 0;

    /* Prefetch spinor for firt case: 0+ */
    iy1 = forward_neighbor(shift_table,iz1,0);
    sp1 = &psi[iy1];
    _prefetch_single(sp1);
    iy2 = forward_neighbor(shift_table,iz1+1,0);
    sp2 = &psi[iy2];
    _prefetch_single(sp2);

    /* Site 1: Load */
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);

    /* Project */
    _sse_42_gamma3_minus();
      
    /* Multiply */
    _sse_su3_inverse_multiply((*um1));

    /* Recons/Accumulate and set in res */
    sn1 = &res[ix1];     
   
    /* Top half spinor */
    _sse_vector_load(r12_1);
    _sse_24_gamma3_minus_rows12();
    _sse_pair_store((*sn1)[0],(*sn1)[1]);

    /* Bottom half spinor */
    _sse_vector_load(r34_1);
    _sse_24_gamma3_minus();
    _sse_pair_store((*sn1)[2],(*sn1)[3]);      

    /* Prefetch Gauge field for next case: 0+ */
    up1 = &(gauge_field[iz1][0]);
    _prefetch_single(up1);
    up2 = &(gauge_field[iz1+1][0]);
    _prefetch_single(up2);

    /* Site 2: Load */
    _sse_pair_load((*sm2)[0],(*sm2)[1]);
    _sse_pair_load_up((*sm2)[2],(*sm2)[3]);

    /* Project */
    _sse_42_gamma3_minus();
      
    /* Multiply */
    _sse_su3_inverse_multiply((*um2));

    /* Recons/Accumulate and set in res */
    /* Top Half */
    _sse_vector_load(r12_2);
    _sse_24_gamma3_minus_rows12();
    _sse_pair_store((*(sn1+1))[0],(*(sn1+1))[1]);

    /* Bottom Half */
    _sse_vector_load(r34_2);
    _sse_42_gamma3_minus();
    _sse_pair_store((*(sn1+1))[2],(*(sn1+1))[3]); 
    /******************************** end of loop *********************************/
  }

}

#ifdef __cplusplus
}
#endif
