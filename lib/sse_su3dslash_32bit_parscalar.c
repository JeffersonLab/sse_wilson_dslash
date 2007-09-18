/*******************************************************************************
 * $Id: sse_su3dslash_32bit_parscalar.c,v 1.4 2007-09-18 19:25:23 bjoo Exp $
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
#include "shift_tables_parscalar.h"
#include <sse_align.h>

static int initP=0;

static int *shift_table;
static int icolor_start[2];    /* starting site for each coloring (cb) */
static int icolor_end[2];      /* end site for each coloring (cb) */

#define buffer_address(chi,mymu,mysite) ((chi)+mysite+3*subgrid_vol_cb*mymu)

  /* The gauge field is packed so that:
       gauge_field[site][0]   <-> U(x=site   , dir = 0 ) 
       gauge_field[site][1]   <-> U(x=site+1 , dir = 0 )
       gauge_field[site][2]   <-> U(x=site   , dir = 1 )
       gauge_field[site][3]   <-> U(x=site+1 , dir = 1 )
       gauge_field[site+1][0] <-> U(x=site   , dir = 2 )
       gauge_field[site+1][1] <-> U(x=site+1 , dir = 2 )
       gauge_field[site+1][2] <-> U(x=site   , dir = 3 )
       gauge_field[site+1][3] <-> U(x=site+1 , dir = 3 )
      
       Here: gauge_field is the packed gauge field and U is the corresponding unpacked gauge field 
  */


  /* now overlays for spinors as arrays or structs. Most important are: 
     u_mat_array - is a single link matrix 
     my_mat_array - is a 4-vector of u_mat_array pointers
     half_spinor_array - is a 2 component vector of color vectors 
     spinor_array - which is a 4 component vector of color vectors */
  
  typedef float u_mat_array[3][3][2]  ALIGN;       /* color color re/im */ 
  typedef float spinor_array[4][3][2] ALIGN;       /* Nspin4 color re/im */
  typedef float halfspinor_array[3][2][2]    ALIGN;    /* Half Spinor re/im,spin2,color */
  typedef u_mat_array (*my_mat_array)[4] ALIGN;  


  /* Spin Matrices */
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


  /* Parameter struct for thread workers */
  typedef struct {
    spinor_array* spinor;           /* Spinor either read or write */
    halfspinor_array *half_spinor;  /* Half Spinor - either read or write */
    u_mat_array       (*u)[4];      /* Gauge field - packed */
    int cb;                         /* Checkerboard (source) */
  } Arg_s;
  


/****************************isign corresponding to +1  **************************/

/* the basic operation here is load up a spinor, do the spin decomposition, and store the halfspinor
to a lattice temporary */

  void decomp_plus(size_t lo,size_t hi, int id, const void *ptr) /*need to fix decomp_minus */
  {

    int ix1,iy2,iz1;
    spinor_array* s1 ALIGN;
    spinor_array* sp1 ALIGN;
    spinor_array* sp2 ALIGN;
    spinor_array* sm2 ALIGN;
    spinor_array* sm1 ALIGN;
    spinor_array* sn1 ALIGN;
    sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};
    sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};
    sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};
    sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};
    sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};
    sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};
    sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f};
    
    int subgrid_vol_cb = getSubgridVolCB();
    halfspinor_array r12_1 ALIGN,r34_1 ALIGN,r12_2 ALIGN,r34_2 ALIGN;
    
    const Arg_s *a = (Arg_s *)ptr;
    
    halfspinor_array* chia = a->half_spinor; /* needs to be changed to halfspinor_array and be an array*/
    
    halfspinor_array* s3;
    halfspinor_array* s4;
    
    int cb = a->cb;
    
    int  low = icolor_start[cb]+(int)lo;
    int high = icolor_start[cb]+(int)hi;
    
    spinor_array* spinor_field= a->spinor;
    
    sp1=&spinor_field[low];
    
    
    /************************ loop over all lattice sites *************************/
    for (ix1=low;ix1<high;ix1+=2) {
      
      s1=&spinor_field[ix1+2];
      _prefetch_spinor(s1);
      
      /* prefetched input spinor for next two sites */
      
      /******************************* direction +0 *********************************/
      
      
      /* first of two sites */
      _sse_pair_load((*sp1)[0],(*sp1)[1]);
      s3 = buffer_address(chia,0,decomp_scatter_index(shift_table,ix1,0));
      _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
      _sse_42_gamma0_minus();
      
      /* the halfspinor is now in xmm0-2 , so _store*/
      
      /* note: if the communications hardware likes its buffers out to memory instead of in cache, then non-temporal stores may be in order...check one of the #define switches above for more details */
      _sse_vector_store(*s3);
      
      /* second of two sites */
      _sse_pair_load((*(sp1+1))[0],(*(sp1+1))[1]);
      s3 = buffer_address(chia,0,decomp_scatter_index(shift_table,ix1+1,0));
      _sse_pair_load_up((*(sp1+1))[2],(*(sp1+1))[3]);
      _sse_42_gamma0_minus();
      
      _sse_vector_store(*s3);
      
      
      /******************************* direction +1 *********************************/
      _sse_pair_load((*sp1)[0],(*sp1)[1]);
      _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
      s3 = buffer_address(chia,1,decomp_scatter_index(shift_table,ix1,1));
      _sse_42_gamma1_minus();
      
      _sse_vector_store(*s3);
      
      _sse_pair_load((*(sp1+1))[0],(*(sp1+1))[1]);
      s3 = buffer_address(chia,1,decomp_scatter_index(shift_table,ix1+1,1));
      _sse_pair_load_up((*(sp1+1))[2],(*(sp1+1))[3]);
      _sse_42_gamma1_minus();
      
      _sse_vector_store(*s3);
      
      /******************************* direction +2 *********************************/
      _sse_pair_load((*sp1)[0],(*sp1)[1]);
      s3 = buffer_address(chia,2,decomp_scatter_index(shift_table,ix1,2));
      _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
      _sse_42_gamma2_minus();
      
      _sse_vector_store(*s3);
      
      _sse_pair_load((*(sp1+1))[0],(*(sp1+1))[1]);
      s3 = buffer_address(chia,2,decomp_scatter_index(shift_table,ix1+1,2));
      _sse_pair_load_up((*(sp1+1))[2],(*(sp1+1))[3]);
      _sse_42_gamma2_minus();
      
      _sse_vector_store(*s3);   	
      
      sp2=sp1+1;
      
      /******************************* direction +3 *********************************/
      _sse_pair_load((*sp1)[0],(*sp1)[1]);
      s3 = buffer_address(chia,3,decomp_scatter_index(shift_table,ix1,3));
      _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
      _sse_42_gamma3_minus();
      
      _sse_vector_store(*s3);
      
      _sse_pair_load((*sp2)[0],(*sp2)[1]);
      s3 = buffer_address(chia,3,decomp_scatter_index(shift_table,ix1+1,3));
      _sse_pair_load_up((*sp2)[2],(*sp2)[3]);
      _sse_42_gamma3_minus();
      
      _sse_vector_store(*s3);
      
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

  int ix1,iy1,iy2,iz1;
  u_mat_array* up1 ALIGN;
  u_mat_array* up2 ALIGN;
  u_mat_array* um1 ALIGN;
  u_mat_array* um2 ALIGN;
  u_mat_array* um3 ALIGN;
  spinor_array* s1 ALIGN;
  spinor_array* sp1 ALIGN;
  spinor_array* sp2 ALIGN;
  spinor_array* sm2 ALIGN;
  spinor_array* sm1 ALIGN;
  spinor_array* sn1 ALIGN;
  sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};
  sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};
  sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};
  sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};
  sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};
  sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};
  sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f};

  int subgrid_vol_cb = getSubgridVolCB();
  halfspinor_array r12_1 ALIGN,r34_1 ALIGN,r12_2 ALIGN,r34_2 ALIGN;

  const Arg_s *a = (const Arg_s *)ptr;
  spinor_array* spinor_field = a->spinor;

  halfspinor_array* chib = a->half_spinor; /* a 1-d map of a 2-d array */
  my_mat_array gauge_field = a->u;
  halfspinor_array* s3;
  halfspinor_array* s4;

  int cb = a->cb;
  int  low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;


/************************ loop over all lattice sites *************************/
  for (ix1=low;ix1<high;ix1+=2) 
  {
    /******************************* direction +0 *********************************/
    sm1=&spinor_field[ix1];
    um1=&gauge_field[ix1][0]; 
    um2=&gauge_field[ix1][1];
    
    /******************************* direction -0 *********************************/
    /* load the input spinor */  
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    um3=&gauge_field[ix1][2];
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);

    /* prefetch the next direction's worth of gauge fields */
    _prefetch_su3(um3);

    /* do the spin projection */
    _sse_42_gamma0_plus();

    s3 = buffer_address(chib,0,decomp_hvv_scatter_index(shift_table,ix1,0));
    _sse_su3_inverse_multiply((*um1));
    _sse_vector_store_up(*s3);
      
    _sse_pair_load((*(sm1+1))[0],(*(sm1+1))[1]);
    _sse_pair_load_up((*(sm1+1))[2],(*(sm1+1))[3]);
     
    _sse_42_gamma0_plus();
    s4 = buffer_address(chib,0,decomp_hvv_scatter_index(shift_table,ix1+1,0));
    /*_prefetch_single(s4);*/

    _sse_su3_inverse_multiply((*(um2)));
    _sse_vector_store_up(*s4);

    /******************************* direction +1 *********************************/
    um1=um3;
    um2=&gauge_field[ix1][3];
   
    /******************************* direction -1 *********************************/
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
	  
    um3=&gauge_field[ix1+1][0];
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);
    _prefetch_su3(um3);
    _sse_42_gamma1_plus();
	  
    s3 = buffer_address(chib,1,decomp_hvv_scatter_index(shift_table,ix1,1));

    _sse_su3_inverse_multiply((*um1));
    _sse_vector_store_up(*s3);

    _sse_pair_load((*(sm1+1))[0],(*(sm1+1))[1]);
    _sse_pair_load_up((*(sm1+1))[2],(*(sm1+1))[3]);
    _sse_42_gamma1_plus();
    s4 = buffer_address(chib,1,decomp_hvv_scatter_index(shift_table,ix1+1,1));

    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_store_up(*s4);


    /******************************* direction +2 *********************************/
    um1=um3;
    um2=&gauge_field[ix1+1][1];
    
    /******************************* direction -2 *********************************/
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    um3=&gauge_field[ix1+1][2];
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);
    _prefetch_su3(um3);

    _sse_42_gamma2_plus();      

    s3 = buffer_address(chib,2,decomp_hvv_scatter_index(shift_table,ix1,2));
    _sse_su3_inverse_multiply((*um1));

    _sse_vector_store_up(*s3);
      
    _sse_pair_load((*(sm1+1))[0],(*(sm1+1))[1]);
    _sse_pair_load_up((*(sm1+1))[2],(*(sm1+1))[3]);
    _sse_42_gamma2_plus();      

    s4 = buffer_address(chib,2,decomp_hvv_scatter_index(shift_table,ix1+1,2));
    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_store_up(*s4);
     

    /******************************* direction +3 *********************************/
    um1=um3;
    um2=&gauge_field[ix1+1][3];
    sm2=sm1+1;
	  
    /******************************* direction -3 *********************************/
    iz1=ix1+2;
      
    s3 = buffer_address(chib,3,decomp_hvv_scatter_index(shift_table,ix1,3));
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);
    sm1 = &spinor_field[iz1];

    _sse_42_gamma3_plus();
	  
    /*_prefetch_single(s3);*/
    _prefetch_spinor(sm1);
    _sse_su3_inverse_multiply((*um1));
    _sse_vector_store_up(*s3);
      
    um1=&gauge_field[iz1][0];  /* gauge packed or not this is the same */
    _prefetch_su3(um1);
      
    _sse_pair_load((*sm2)[0],(*sm2)[1]);
    _sse_pair_load_up((*sm2)[2],(*sm2)[3]);
    s4 = buffer_address(chib,3,decomp_hvv_scatter_index(shift_table,ix1+1,3));
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

  int ix1,iy1,iy2,iz1;
  u_mat_array* up1 ALIGN;
  u_mat_array* up2 ALIGN;
  u_mat_array* um1 ALIGN;
  u_mat_array* um2 ALIGN;
  u_mat_array* um3 ALIGN;
  spinor_array* s1 ALIGN;
  spinor_array* sp1 ALIGN;
  spinor_array* sp2 ALIGN;
  spinor_array* sm2 ALIGN;
  spinor_array* sm1 ALIGN;
  spinor_array* sn1 ALIGN;
  sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};
  sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};
  sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};
  sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};
  sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};
  sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};
  sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f};

  int subgrid_vol_cb = getSubgridVolCB();
  halfspinor_array r12_1 ALIGN,r34_1 ALIGN,r12_2 ALIGN,r34_2 ALIGN;


  const Arg_s *a =(Arg_s *)ptr;

  spinor_array* spinor_field = a->spinor;

  halfspinor_array* chia = a->half_spinor; /* a 1-d map of a 2-d array */
  my_mat_array gauge_field = a->u;
  halfspinor_array* s3;
  halfspinor_array* s4;
  int cb = a->cb;
  int  low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;

  s3 = buffer_address(chia,0,recons_mvv_gather_index(shift_table,low,0));
  _prefetch_single_nta(s3);
  iy1=decomp_hvv_scatter_index(shift_table,low,0);
  sp1=&spinor_field[iy1];
  up1=&gauge_field[low][0];
  up2=&gauge_field[low][1];


  /************************ loop over all lattice sites *************************/
  for (ix1=low;ix1<high;ix1+=2) 
  {
    /* s1=&spinor_field[ix1];
       _prefetch_spinor(s1);*/
      
    /******************************* direction +0 *********************************/
    /* ...(1-gamma(0))... */
    /* load from the temporary */
    _sse_vector_load(*s3);
    s4 = buffer_address(chia,0,recons_mvv_gather_index(shift_table,ix1+1,0));
    
    /*prefetch the next temp into one way of lvl2 cache, and prefetch the next gague field */
    _prefetch_single_nta(s4);
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
    s3 = buffer_address(chia,1,recons_mvv_gather_index(shift_table,ix1,1));
    _prefetch_single_nta(s3);
    _sse_su3_multiply((*(up2)));
    _sse_vector_store_up(r12_2);

    _sse_24_gamma0_minus_set(); 
    _sse_vector_store_up(r34_2);

    /***************************** direction +1 ***********************************/
    up1=&gauge_field[ix1][2];
    up2=&gauge_field[ix1][3];   
      
    /* same kind of thing again */
    _sse_vector_load(*s3);
    s4 = buffer_address(chia,1,recons_mvv_gather_index(shift_table,ix1+1,1));

    _prefetch_single_nta(s4);
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
    s3 = buffer_address(chia,2,recons_mvv_gather_index(shift_table,ix1,2));
    _prefetch_single_nta(s3);

    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);
      
    _sse_vector_load(r34_2);
    _sse_24_gamma1_minus();
    _sse_vector_store(r34_2);

    up1=&gauge_field[ix1+1][0]; /* default is packed gauge fields, thats why the 
				* statements are non-intuitive */
    up2=&gauge_field[ix1+1][1];

    /******************************* direction +2 *********************************/
    _sse_vector_load(*s3);
    s4 = buffer_address(chia,2,recons_mvv_gather_index(shift_table,ix1+1,2));
    _prefetch_single_nta(s4);
    _prefetch_su3(up1+2); 

    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);       

    _sse_vector_load(r34_1);
    _sse_24_gamma2_minus();
    _sse_vector_store(r34_1);       

    _sse_vector_load(*s4);
    s3 = buffer_address(chia,3,recons_mvv_gather_index(shift_table,ix1,3));
    _prefetch_single_nta(s3);
     
    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);       

    _sse_vector_load(r34_2);
    _sse_24_gamma2_minus();
    _sse_vector_store(r34_2); 

    up1=&gauge_field[ix1+1][2];
    up2=&gauge_field[ix1+1][3];

    /******************************* direction +3 *********************************/
    sn1=&spinor_field[ix1];

    _sse_vector_load(*s3);
    s4 = buffer_address(chia,3,recons_mvv_gather_index(shift_table,ix1+1,3));
    
    _prefetch_single_nta(s4);
    /* prefetch the gague field for direction 0, site = ix1+2 */
    _prefetch_su3(up1+2); 
      
    _sse_su3_multiply((*up1));
    
    /* ok things get tricky here again...the rows12 thing is for compatibility 
     * with other spin basis in which the first two rows of that 4x2 
     * reconstruction matrix may be something other than the 2x2 identity matrix */
    _sse_vector_load(r12_1);
    _sse_24_gamma3_minus_rows12();
    _sse_pair_store((*sn1)[0],(*sn1)[1]);

    _sse_vector_load(r34_1);
    _sse_24_gamma3_minus();
    _sse_pair_store((*sn1)[2],(*sn1)[3]);   

    iz1=ix1+2;
    if (iz1==high)
      iz1=0;

    _sse_vector_load(*s4);
    s3 = buffer_address(chia,0,recons_mvv_gather_index(shift_table,iz1,0));
    _prefetch_single_nta(s3);
      
    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_24_gamma3_minus_rows12();
    _sse_pair_store((*(sn1+1))[0],(*(sn1+1))[1]);

    _sse_vector_load(r34_2);
    _sse_24_gamma3_minus();
    _sse_pair_store((*(sn1+1))[2],(*(sn1+1))[3]); 

    up1=&gauge_field[iz1][0];
	 up2=&gauge_field[iz1][1];

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
  int ix1,iy1,iy2,iz1;
  u_mat_array* up1 ALIGN;
  u_mat_array* up2 ALIGN;
  u_mat_array* um1 ALIGN;
  u_mat_array* um2 ALIGN;
  u_mat_array* um3 ALIGN;
  spinor_array* s1 ALIGN;
  spinor_array* sp1 ALIGN;
  spinor_array* sp2 ALIGN;
  spinor_array* sm2 ALIGN;
  spinor_array* sm1 ALIGN;
  spinor_array* sn1 ALIGN;
  sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};
  sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};
  sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};
  sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};
  sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};
  sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};
  sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f};

  int subgrid_vol_cb = getSubgridVolCB();
  halfspinor_array r12_1 ALIGN,r34_1 ALIGN,r12_2 ALIGN,r34_2 ALIGN;

  const Arg_s *a = (Arg_s *)ptr;
  spinor_array* spinor_field = a->spinor;

  halfspinor_array* hs0, *hs1, *hs2,*hs3,*hs4,*hs5,*hs6,*hs7,*hs8;
  halfspinor_array* chib = a->half_spinor; /* a 1-d map of a 2-d array */
  int cb = a->cb;
  int low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;

  /* printf("\nlo:%i, hi:%i, id:%i, chib:%x", lo, hi, id, chib);*/

/************************ loop over all lattice sites *************************/
  for (ix1=low;ix1<high;ix1+=2) 
  {
#define PREFDIST 4
    hs0 = buffer_address(chib,0,recons_gather_index(shift_table,ix1,0));
    sn1=&spinor_field[ix1];   
   
    _prefetch_spinor_nta(sn1+PREFDIST);
     
    /***** psi 1&2 site 1 ******/
    _sse_vector_load_up(*(hs0));   /* vector in xmm3-5 */
	                             
    /* accumulate in xmm0-2 */
    _sse_pair_load((*sn1)[0], (*sn1)[1]); /*load in partial sum */
	  
    hs1 = buffer_address(chib,1,recons_gather_index(shift_table,ix1,1));
    _sse_vector_add();
    _sse_vector_load_up(*(hs1)); /*direction +1 */ 
	 
    hs2 = buffer_address(chib,2,recons_gather_index(shift_table,ix1,2));
    _sse_vector_add();          /* accumulating in xmm0-2 */
    _sse_vector_load_up(*(hs2));  /* direction +2 */
	
    hs3 = buffer_address(chib,3,recons_gather_index(shift_table,ix1,3));
    _sse_vector_add();
    _sse_vector_load_up(*(hs3)); /* direction +3 */
   
    _sse_24_gamma3_plus_rows12();
    _sse_pair_store((*sn1)[0],(*sn1)[1]);

    /***** psi's 3 and 4 sites 1*****/
    _sse_vector_load_up(*(hs0));       /* vector in xmm3-5 */
   
    /* accumulate in xmm0-2 */
    _sse_pair_load((*sn1)[2],(*sn1)[3]); /*load in partial sum */

    _sse_24_gamma0_plus_add();

    _sse_vector_load_up(*(hs1));
    
    _sse_24_gamma1_plus();

    _sse_vector_load_up(*(hs2));
   
    _sse_24_gamma2_plus();
      
    _sse_vector_load_up(*(hs3));
    _sse_24_gamma3_plus();
       
    _sse_pair_store((*sn1)[2],(*sn1)[3]);

    _prefetch_single(hs2);


    /****** psi 1&2 site 2 ******/
    hs0 = buffer_address(chib,0,recons_gather_index(shift_table,ix1+1,0)); 
    _sse_vector_load_up(*(hs0));   /* vector in xmm3-5 */
    hs1 = buffer_address(chib,1,recons_gather_index(shift_table,ix1+1,1));
	  
    /* accumulate in xmm0-2 */
    _sse_pair_load((*(sn1+1))[0],(*(sn1+1))[1]); /*load in partial sum */
   
    _sse_vector_add();
    _sse_vector_load_up(*(hs1)); /*direction +1 */
    hs2 = buffer_address(chib,2,recons_gather_index(shift_table,ix1+1,2));
      
    _sse_vector_add();          /* accumulating in xmm0-2 */
	  
    _sse_vector_load_up(*(hs2)); /* direction +2 */
    hs3 = buffer_address(chib,3,recons_gather_index(shift_table,ix1+1,3));

    _sse_vector_add();

    _sse_vector_load_up(*(hs3)); /* direction +3 */
     
    _sse_24_gamma3_plus_rows12();

    _sse_pair_store((*(sn1+1))[0],(*(sn1+1))[1]);
	 
      
    /***** psi's 3 and 4 site 2 *****/
    _sse_vector_load_up(*(hs0));       /* vector in xmm3-5 */
     
    _sse_pair_load((*(sn1+1))[2],(*(sn1+1))[3]); /*load in partial sum */
    _sse_24_gamma0_plus_add();

    _sse_vector_load_up(*(hs1));
    
    _sse_24_gamma1_plus();
    _sse_vector_load_up(*(hs2));
    
    _sse_24_gamma2_plus();
    _sse_vector_load_up(*(hs3));
    _sse_24_gamma3_plus();
   
    _sse_pair_store((*(sn1+1))[2],(*(sn1+1))[3]); 
  
    /*************************end of loop ****************************/
  }
}
/*****************end of recons**************/





/*************** now for isign corresponding to -1  ****************************************/

void decomp_minus(size_t lo,size_t hi, int id, const void *ptr ) /*need to fix decomp_minus */
{

  int ix1,iy1,iy2,iz1;
  u_mat_array* up1 ALIGN;
  u_mat_array* up2 ALIGN;
  u_mat_array* um1 ALIGN;
  u_mat_array* um2 ALIGN;
  u_mat_array* um3 ALIGN;
  spinor_array* s1 ALIGN;
  spinor_array* sp1 ALIGN;
  spinor_array* sp2 ALIGN;
  spinor_array* sm2 ALIGN;
  spinor_array* sm1 ALIGN;
  spinor_array* sn1 ALIGN;
  sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};
  sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};
  sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};
  sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};
  sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};
  sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};
  sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f};

  int subgrid_vol_cb = getSubgridVolCB();
  halfspinor_array r12_1 ALIGN,r34_1 ALIGN,r12_2 ALIGN,r34_2 ALIGN;

  const Arg_s *a =(Arg_s *)ptr;

  halfspinor_array* chia = a->half_spinor; /* needs to be changed to halfspinor_array and be an array*/
  halfspinor_array* s3;
  halfspinor_array* s4;
  int cb = a->cb;
  int  low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;

  spinor_array* spinor_field= a->spinor;
   

  iy1=decomp_hvv_scatter_index(shift_table,low,0);
  sp1=&spinor_field[low];
 
  /************************ loop over all lattice sites *************************/

  for (ix1=low;ix1<high;ix1+=2) 
  {
    s1=&spinor_field[ix1+2];
    _prefetch_spinor(s1);
      
    /******************************* direction +0 *********************************/
    /* ...(1-gamma(0))... */
    _sse_pair_load((*sp1)[0],(*sp1)[1]);
    s3 = buffer_address(chia,0,decomp_scatter_index(shift_table,ix1,0));
    _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
    _sse_42_gamma0_plus();
    _sse_vector_store(*s3);
      
    _sse_pair_load((*(sp1+1))[0],(*(sp1+1))[1]);
    s4 = buffer_address(chia,0,decomp_scatter_index(shift_table,ix1+1,0));
    _sse_pair_load_up((*(sp1+1))[2],(*(sp1+1))[3]);
    _sse_42_gamma0_plus();
    _sse_vector_store(*s4);

    /******************************* direction -0 *********************************/
   

    /******************************* direction +1 *********************************/
    _sse_pair_load((*sp1)[0],(*sp1)[1]);
    _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
    s3 = buffer_address(chia,1,decomp_scatter_index(shift_table,ix1,1));
    _sse_42_gamma1_plus();
    _sse_vector_store(*s3);
      
    _sse_pair_load((*(sp1+1))[0],(*(sp1+1))[1]);
    s4 = buffer_address(chia,1,decomp_scatter_index(shift_table,ix1+1,1));
    _sse_pair_load_up((*(sp1+1))[2],(*(sp1+1))[3]);
    _sse_42_gamma1_plus();
    _sse_vector_store(*s4);

    /******************************* direction -1 *********************************/
     

/******************************* direction +2 *********************************/
    _sse_pair_load((*sp1)[0],(*sp1)[1]);
    s3 = buffer_address(chia,2,decomp_scatter_index(shift_table,ix1,2));
    _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
    _sse_42_gamma2_plus();
    _sse_vector_store(*s3);

    _sse_pair_load((*(sp1+1))[0],(*(sp1+1))[1]);
    s4 = buffer_address(chia,2,decomp_scatter_index(shift_table,ix1+1,2));
    _sse_pair_load_up((*(sp1+1))[2],(*(sp1+1))[3]);
    _sse_42_gamma2_plus();
    _sse_vector_store(*s4);

    /******************************* direction -2 *********************************/
    sp2=sp1+1;

    /******************************* direction +3 *********************************/
    _sse_pair_load((*sp1)[0],(*sp1)[1]);
    s3 = buffer_address(chia,3,decomp_scatter_index(shift_table,ix1,3));
    _sse_pair_load_up((*sp1)[2],(*sp1)[3]);
    _sse_42_gamma3_plus();
    _sse_vector_store(*s3);
      
    _sse_pair_load((*sp2)[0],(*sp2)[1]);
    s4 = buffer_address(chia,3,decomp_scatter_index(shift_table,ix1+1,3));
    _sse_pair_load_up((*sp2)[2],(*sp2)[3]);
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

  int ix1,iy1,iy2,iz1;
  u_mat_array* up1 ALIGN;
  u_mat_array* up2 ALIGN;
  u_mat_array* um1 ALIGN;
  u_mat_array* um2 ALIGN;
  u_mat_array* um3 ALIGN;
  spinor_array* s1 ALIGN;
  spinor_array* sp1 ALIGN;
  spinor_array* sp2 ALIGN;
  spinor_array* sm2 ALIGN;
  spinor_array* sm1 ALIGN;
  spinor_array* sn1 ALIGN;
  sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};
  sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};
  sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};
  sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};
  sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};
  sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};
  sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f};

  int subgrid_vol_cb = getSubgridVolCB();
  halfspinor_array r12_1 ALIGN,r34_1 ALIGN,r12_2 ALIGN,r34_2 ALIGN;

  const Arg_s *a =(Arg_s *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* chib = a->half_spinor; /* a 1-d map of a 2-d array */
  my_mat_array gauge_field = a->u;
  halfspinor_array* s3;
  halfspinor_array* s4;

  int cb = a->cb;
  int  low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;



  /************************ loop over all lattice sites *************************/
  for (ix1=low;ix1<high;ix1+=2) 
  {
    /******************************* direction +0 *********************************/
    /* ...(1+gamma(0))... */
    sm1=&spinor_field[ix1];
    um1=&gauge_field[ix1][0]; 
    um2=&gauge_field[ix1][1];
    

    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    um3=&gauge_field[ix1][2];
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);
    _prefetch_su3(um3);
	  
    _sse_42_gamma0_minus();
    s3 = buffer_address(chib,0,decomp_hvv_scatter_index(shift_table,ix1,0));
    _prefetch_single(s3);

    _sse_su3_inverse_multiply((*um1));
    _sse_vector_store_up(*s3);

    _sse_pair_load((*(sm1+1))[0],(*(sm1+1))[1]);
    _sse_pair_load_up((*(sm1+1))[2],(*(sm1+1))[3]);
	  
    _sse_42_gamma0_minus();
    s4 = buffer_address(chib,0,decomp_hvv_scatter_index(shift_table,ix1+1,0));
    _prefetch_single(s4);
    _sse_su3_inverse_multiply((*(um2)));
    _sse_vector_store_up(*s4);
      
    /******************************* direction +1 *********************************/
    um1=um3;
    um2=&gauge_field[ix1][3];
    
    /******************************* direction -1 *********************************/
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
	  
    um3=&gauge_field[ix1+1][0];
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);
    _prefetch_su3(um3);
    _sse_42_gamma1_minus();
	  
    s3 = buffer_address(chib,1,decomp_hvv_scatter_index(shift_table,ix1,1));
    _prefetch_single(s3);
    _sse_su3_inverse_multiply((*um1));

    _sse_vector_store_up(*s3);

    _sse_pair_load((*(sm1+1))[0],(*(sm1+1))[1]);
    _sse_pair_load_up((*(sm1+1))[2],(*(sm1+1))[3]);
    _sse_42_gamma1_minus();
    s4 = buffer_address(chib,1,decomp_hvv_scatter_index(shift_table,ix1+1,1));
    _prefetch_single(s4);
    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_store_up(*s4);

/******************************* direction +2 *********************************/
    um1=um3;
    um2=&gauge_field[ix1+1][1];

    /******************************* direction -2 *********************************/
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    um3=&gauge_field[ix1+1][2];
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);
    _prefetch_su3(um3);

    _sse_42_gamma2_minus();      

    s3 = buffer_address(chib,2,decomp_hvv_scatter_index(shift_table,ix1,2));
    _prefetch_single(s3);
    _sse_su3_inverse_multiply((*um1));

    _sse_vector_store_up(*s3);
      
    _sse_pair_load((*(sm1+1))[0],(*(sm1+1))[1]);
    _sse_pair_load_up((*(sm1+1))[2],(*(sm1+1))[3]);
    _sse_42_gamma2_minus();      

    s4 = buffer_address(chib,2,decomp_hvv_scatter_index(shift_table,ix1+1,2));
    _prefetch_single(s4);
    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_store_up(*s4);
     
    /******************************* direction +3 *********************************/
    um1=um3;
    um2=&gauge_field[ix1+1][3];
    sm2=sm1+1;

    /******************************* direction -3 *********************************/
    iz1=ix1+2;
    if (iz1==high)
      iz1=0;
      
    s3 = buffer_address(chib,3,decomp_hvv_scatter_index(shift_table,ix1,3));
    _sse_pair_load((*sm1)[0],(*sm1)[1]);
    _sse_pair_load_up((*sm1)[2],(*sm1)[3]);
    sm1 = &spinor_field[iz1];

    _sse_42_gamma3_minus();
	  
    _prefetch_single(s3);
    _prefetch_spinor(sm1);

    _sse_su3_inverse_multiply((*um1));
    _sse_vector_store_up(*s3);
      
    um1=&gauge_field[iz1][0];  /* gauge packed or not this is the same */
    _prefetch_su3(um1);
      
    _sse_pair_load((*sm2)[0],(*sm2)[1]);
    _sse_pair_load_up((*sm2)[2],(*sm2)[3]);
    s4 = buffer_address(chib,3,decomp_hvv_scatter_index(shift_table,ix1+1,3));
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
  int ix1,iy1,iy2,iz1;
  u_mat_array* up1 ALIGN;
  u_mat_array* up2 ALIGN;
  u_mat_array* um1 ALIGN;
  u_mat_array* um2 ALIGN;
  u_mat_array* um3 ALIGN;
  spinor_array* s1 ALIGN;
  spinor_array* sp1 ALIGN;
  spinor_array* sp2 ALIGN;
  spinor_array* sm2 ALIGN;
  spinor_array* sm1 ALIGN;
  spinor_array* sn1 ALIGN;
  sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};
  sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};
  sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};
  sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};
  sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};
  sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};
  sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f};

  int subgrid_vol_cb = getSubgridVolCB();
  halfspinor_array r12_1 ALIGN,r34_1 ALIGN,r12_2 ALIGN,r34_2 ALIGN;

  /* if going to support unpacked gauge fields, need to treat site ix1 and site ix1+1 separately */
  /* to support unpacked gauge fields the prefetches will need to be changed */
  const Arg_s *a =(Arg_s *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* chia = a->half_spinor; /* a 1-d map of a 2-d array */
  my_mat_array gauge_field = a->u;
  halfspinor_array* s3;
  halfspinor_array* s4;
  int cb = a->cb;
  int  low = icolor_start[cb]+(int)lo;
  int high = icolor_start[cb]+(int)hi;

  s3 = buffer_address(chia,0,recons_mvv_gather_index(shift_table,low,0));
  _prefetch_single_nta(s3);
  iy1=decomp_hvv_scatter_index(shift_table,low,0);
  sp1=&spinor_field[iy1];
  up1=&gauge_field[low][0];
  up2=&gauge_field[low][1];

/************************ loop over all lattice sites *************************/
  for (ix1=low;ix1<high;ix1+=2) 
  {
    /******************************* direction +0 *********************************/
    _sse_vector_load(*s3);
    s4 = buffer_address(chia,0,recons_mvv_gather_index(shift_table,ix1+1,0));
    _prefetch_single_nta(s4);
    _prefetch_su3(up1+2); 
      
    _sse_su3_multiply((*up1));
    _sse_vector_store_up(r12_1);

    _sse_24_gamma0_plus_set(); 
    _sse_vector_store_up(r34_1);
    
    _sse_vector_load(*s4);
    s3 = buffer_address(chia,1,recons_mvv_gather_index(shift_table,ix1,1));
    _prefetch_single_nta(s3);
    _sse_su3_multiply((*(up2)));
    _sse_vector_store_up(r12_2);

    _sse_24_gamma0_plus_set(); 
    _sse_vector_store_up(r34_2);
      
    /******************************* direction -0 *********************************/
    /* ...(1-gamma(0))... */
    up1=&gauge_field[ix1][2];
    up2=&gauge_field[ix1][3];   
      
/******************************* direction +1 *********************************/
    _sse_vector_load(*s3);
    s4 = buffer_address(chia,1,recons_mvv_gather_index(shift_table,ix1+1,1));
    _prefetch_single_nta(s4);
    _prefetch_su3(up1+2); 

    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);
      
    _sse_vector_load(r34_1);
    _sse_24_gamma1_plus();
    _sse_vector_store(r34_1);

      
    _sse_vector_load(*s4);
    s3 = buffer_address(chia,2,recons_mvv_gather_index(shift_table,ix1,2));
    _prefetch_single_nta(s3);

    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);
      
    _sse_vector_load(r34_2);
    _sse_24_gamma1_plus();
    _sse_vector_store(r34_2);

    /******************************* direction -1 *********************************/
    up1=&gauge_field[ix1+1][0]; /* default is packed gauge fields, thats 
				* why the statements are non-intuitive */
    up2=&gauge_field[ix1+1][1];

    /******************************* direction +2 *********************************/
    _sse_vector_load(*s3);
    s4 = buffer_address(chia,2,recons_mvv_gather_index(shift_table,ix1+1,2));
    _prefetch_single_nta(s4);
    _prefetch_su3(up1+2); 

    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);       

    _sse_vector_load(r34_1);
    _sse_24_gamma2_plus();
    _sse_vector_store(r34_1);       

    _sse_vector_load(*s4);
    s3 = buffer_address(chia,3,recons_mvv_gather_index(shift_table,ix1,3));
    _prefetch_single_nta(s3);

    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);       

    _sse_vector_load(r34_2);
    _sse_24_gamma2_plus();
    _sse_vector_store(r34_2); 

    /******************************* direction -2 *********************************/
    up1=&gauge_field[ix1+1][2];
    up2=&gauge_field[ix1+1][3];

/******************************* direction +3 *********************************/
    sn1=&spinor_field[ix1];

    _sse_vector_load(*s3);
    s4 = buffer_address(chia,3,recons_mvv_gather_index(shift_table,ix1+1,3));
    _prefetch_single_nta(s4);
    _prefetch_su3(up1+2); 
      
    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);
    _sse_24_gamma3_plus_rows12();
    _sse_pair_store((*sn1)[0],(*sn1)[1]);

    _sse_vector_load(r34_1);
    _sse_24_gamma3_plus();
    _sse_pair_store((*sn1)[2],(*sn1)[3]);   

    iz1=ix1+2;
    if (iz1==high)
      iz1=0;

    _sse_vector_load(*s4);
    s3 = buffer_address(chia,0,recons_mvv_gather_index(shift_table,iz1,0));
    _prefetch_single_nta(s3);
      
    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_24_gamma3_plus_rows12();
    _sse_pair_store((*(sn1+1))[0],(*(sn1+1))[1]);

    _sse_vector_load(r34_2);
    _sse_24_gamma3_plus();
    _sse_pair_store((*(sn1+1))[2],(*(sn1+1))[3]); 

    /******************************* direction -3 *********************************/
    up1=&gauge_field[iz1][0];
    up2=&gauge_field[iz1][1];

    /******************************** end of loop *********************************/
  }
}
/******************end of mvv_recons*************************/


void recons_minus(size_t lo,size_t hi, int id, const void *ptr )	
{
  int ix1,iy1,iy2,iz1;
  u_mat_array* up1 ALIGN;
  u_mat_array* up2 ALIGN;
  u_mat_array* um1 ALIGN;
  u_mat_array* um2 ALIGN;
  u_mat_array* um3 ALIGN;
  spinor_array* s1 ALIGN;
  spinor_array* sp1 ALIGN;
  spinor_array* sp2 ALIGN;
  spinor_array* sm2 ALIGN;
  spinor_array* sm1 ALIGN;
  spinor_array* sn1 ALIGN;
  sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};
  sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};
  sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};
  sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};
  sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};
  sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};
  sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f};

  int subgrid_vol_cb = getSubgridVolCB();
  halfspinor_array r12_1 ALIGN,r34_1 ALIGN,r12_2 ALIGN,r34_2 ALIGN;

  const Arg_s *a = (Arg_s *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* hs0, *hs1, *hs2, *hs3, *hs4, *s2;
  halfspinor_array* chib = a->half_spinor; /* a 1-d map of a 2-d array */
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
    hs0 = buffer_address(chib,0,recons_gather_index(shift_table,ix1,0));
    sn1=&spinor_field[ix1];   
    
    _prefetch_spinor_nta(sn1+PREFDIST);
    
    /* need to do psum[ix1][0] first for all +mu */
    /* loop over our two adjacent sites slowest inside here */

    /***** psi 1&2 site 1 ******/
    _sse_vector_load_up(*(hs0));   /* vector in xmm3-5 */
     
    /* accumulate in xmm0-2 */
    _sse_pair_load((*sn1)[0], (*sn1)[1]); /*load in partial sum */
	  
    hs1 = buffer_address(chib,1,recons_gather_index(shift_table,ix1,1));
    _sse_vector_add();
    _sse_vector_load_up(*(hs1)); /*direction +1 */ 
	 
    hs2 = buffer_address(chib,2,recons_gather_index(shift_table,ix1,2));
    _sse_vector_add();          /* accumulating in xmm0-2 */
    _sse_vector_load_up(*(hs2));  /* direction +2 */
	
    hs3 = buffer_address(chib,3,recons_gather_index(shift_table,ix1,3));
    _sse_vector_add();
    _sse_vector_load_up(*(hs3)); /* direction +3 */
    
    _sse_24_gamma3_minus_rows12();  /* here's that thing that doesn't 
				     * reduce to an _add in some spin basis */
    _sse_pair_store((*sn1)[0],(*sn1)[1]);

    /***** psi's 3 and 4 sites 1*****/
    _sse_vector_load_up(*(hs0));       /* vector in xmm3-5 */
        
    /* accumulate in xmm0-2 */
    _sse_pair_load((*sn1)[2],(*sn1)[3]); /*load in partial sum */
	   
    _sse_24_gamma0_minus_add();

    _sse_vector_load_up(*(hs1));
    _sse_24_gamma1_minus();

    _sse_vector_load_up(*(hs2));
    _sse_24_gamma2_minus();
    hs0 = buffer_address(chib,0,recons_gather_index(shift_table,ix1+1,0)); 
      
    _sse_vector_load_up(*(hs3));
    _sse_24_gamma3_minus();
       
    _sse_pair_store((*sn1)[2],(*sn1)[3]);
    _prefetch_single(s2);

    /****** psi 1&2 site 2 ******/
    _sse_vector_load_up(*(hs0));   /* vector in xmm3-5 */
    hs1 = buffer_address(chib,1,recons_gather_index(shift_table,ix1+1,1));
    _sse_pair_load((*(sn1+1))[0],(*(sn1+1))[1]); /*load in partial sum */
    
    _sse_vector_add();
    _sse_vector_load_up(*(hs1)); /*direction +1 */
    hs2 = buffer_address(chib,2,recons_gather_index(shift_table,ix1+1,2));
	  
    _sse_vector_add();          /* accumulating in xmm0-2 */
    _sse_vector_load_up(*(hs2)); /* direction +2 */
    hs3 = buffer_address(chib,3,recons_gather_index(shift_table,ix1+1,3));
	   
    _sse_vector_add();
    _sse_vector_load_up(*(hs3)); /* direction +3 */
     
    _sse_24_gamma3_minus_rows12();
    _sse_pair_store((*(sn1+1))[0],(*(sn1+1))[1]);
	 
    /***** psi's 3 and 4 site 2 *****/
    _sse_vector_load_up(*(hs0));       /* vector in xmm3-5 */
  
    /* accumulate in xmm0-2 */
    _sse_pair_load((*(sn1+1))[2],(*(sn1+1))[3]); /*load in partial sum */
    _sse_24_gamma0_minus_add();

    _sse_vector_load_up(*(hs1));
  
    _sse_24_gamma1_minus();
    _sse_vector_load_up(*(hs2));
   
    _sse_24_gamma2_minus();
    _sse_vector_load_up(*(hs3));
    _sse_24_gamma3_minus();

    _sse_pair_store((*(sn1+1))[2],(*(sn1+1))[3]); 
    /*************************end of loop ****************************/
  }
}
/*****************end of recons_minus**************/



/*****************end of isign corresponding to -1 **************************************/


/***************** start of initialization routine ***************************************/


static QMP_mem_t* xchi1;               /* QMP Memory Structures for halfspinor arrays */
static QMP_mem_t* xchi2;               /* xchi1 => FORWARD, xchi2 => BACKWARD         */

static halfspinor_array* chi1;         /* These are the aligned pointers from the QMP Memory structures */
static halfspinor_array* chi2;         /* xchi1 <=> chi1    xchi2 <=> chi2 */



/* Nearest neighbor communication channels */
static int total_comm = 0;
static QMP_msgmem_t forw_msg[4][2];
static QMP_msgmem_t back_msg[4][2];
static QMP_msghandle_t forw_mh[4][2];
static QMP_msghandle_t back_mh[4][2];
static QMP_msghandle_t forw_all_mh;
static QMP_msghandle_t back_all_mh;

/* Initialize the Dslash */
void init_sse_su3dslash(const int latt_size[])   // latt_size not used, here for scalar version
{
  const int *machine_size = QMP_get_logical_dimensions();
  const int *subgrid_cb_size = QMP_get_subgrid_dimensions();
  int bound[2][4][4];
  int mu, num, nsize;



  /* If we are already initialised, then increase the refcount and return */
  if (initP > 0) 
  {
    initP++;
    return;
  }


  /* Otherwise initialise */
  /* Check we are in 4D */
  if (QMP_get_logical_number_of_dimensions() != 4) {
    QMP_error("init_sse_su3dslash: number of logical dimensions does not match problem");
    QMP_abort(1);
  }
    

  /* Check problem size */
  for(mu=0; mu < 4; mu++) {
    if ( latt_size[mu] == 1 ) {
      QMP_error("This SSE Dslash does not support a problem size = 1. Here the lattice in dimension %d has length %d\n", mu, latt_size[mu]);
      QMP_abort(1);
    }
  }

  num = latt_size[0] / machine_size[0];
  if ( num % 2 != 0 )
  {
    QMP_error("This SSE Dslash does not work for odd x-sublattice. Here the sublattice is odd in dimension 0 with length %d\n", num);
    QMP_abort(1);
  }


  /* Make the shift table  -- this sets the vol and vol_cb so we can call getSubgridVolCB() after it */
  shift_table = (int *)make_shift_tables(icolor_start, bound);
  
  icolor_end[0] = icolor_start[0] + getSubgridVolCB();
  icolor_end[1] = icolor_start[1] + getSubgridVolCB();

  /* Allocated space for the floating temps */
  /* Wasteful - allocate 3 times subgrid_vol_cb. Otherwise, need to pack the TAIL{1,2} offsets */
  nsize = 2*3*2*sizeof(float)*getSubgridVolCB()*3*4;  /* Note 3x4 half-fermions */

  /* xchi1 and xchi2 will hold projected spinors memory handles from QMP */
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

  /* Unwrap the half spinor pointers from the QMP Structures. BTW: This 2 step technique susks so bad! */
  chi1 = (halfspinor_array*)QMP_get_memory_pointer(xchi1);
  chi2 = (halfspinor_array*)QMP_get_memory_pointer(xchi2); 

  /* Loop over all communicating directions and build up the two message
   * handles. If there is no communications, the message handles will not
   * be initialized 
   */
  num = 0;

  /* Loop over directions */
  for(mu=0; mu < 4; ++mu) {

    if(machine_size[mu] > 1) { /* If the machine is not a scalar  in this dimensio */
    
      if (bound[0][0][mu] == 0) { /* Consistency: Check the boundary in this direction is 0 */
 	QMP_error("init_sse_dslash: type 0 message size is 0");
	QMP_abort(1);
      }

      /* Boundary is indexed as            boundary[cb][0/1][direction] 
	 and it has some symmetries: 
	 1) Scatter tble boundary on one checkerboard has to be the same as the corresponding
	    gather table on the same checkerboard so:

	         bound[cb][0][dir] = bound[cb][2][dir]  --  0 scatter index, 2 corresp. gather index 
		 bound[cb][1][dir] = bound[cb][3][dir]  --  1 scatter index, 3 corresp. gather index

	 2) A checkerboard has the same size boundary in a direction as another 
	 Follows from the checkerboards having identical geometries...

		 bound[1-cb][i][dir] = bound[cb][i][dir] for i=0,1
      */	 

      /* cb = 0 */
      forw_msg[num][0] = QMP_declare_msgmem( chi1+getSubgridVolCB()*(1+3*mu),
					     bound[0][0][mu]*sizeof(halfspinor_array));

      /* cb = 1 */
      forw_msg[num][1] = QMP_declare_msgmem( chi1+getSubgridVolCB()*(2+3*mu),
					     bound[1][0][mu]*sizeof(halfspinor_array));

      /* cb = 0: Receive from cb = 1 */
      forw_mh[num][0]  = QMP_declare_receive_relative(forw_msg[num][1], mu, +1, 0);

      /* cb = 1: send to cb = 0 */
      forw_mh[num][1]  = QMP_declare_send_relative(forw_msg[num][0], mu, -1, 0);
	
      if (bound[0][1][mu] == 0)
      {
	QMP_error("init_sse_dslash: type 0 message size is 0");
	QMP_abort(1);
      }

      /* cb = 0 */
      back_msg[num][0] = QMP_declare_msgmem( chi2+getSubgridVolCB()*(1+3*mu),
					     bound[0][1][mu]*sizeof(halfspinor_array));

      /* cb = 1 */
      back_msg[num][1] = QMP_declare_msgmem( chi2+getSubgridVolCB()*(2+3*mu), 
					     bound[1][1][mu]*sizeof(halfspinor_array));

      /* cb = 0: Receive from cb=1 */
      back_mh[num][0]  = QMP_declare_receive_relative(back_msg[num][1], mu, -1, 0);

      /* cb = 1: Receive from cb=0 */
      back_mh[num][1]  = QMP_declare_send_relative(back_msg[num][0], mu, +1, 0);
	
      num++;
    }
  }

  /* Combine the messages */
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
    free(shift_table);
    
    if (total_comm > 0) {
      
      QMP_free_msghandle(forw_all_mh);
      QMP_free_msghandle(back_all_mh);
  
      num = 0;
      
      for(mu=0; mu < 4; ++mu) {
	
	if(machine_size[mu] > 1) {
      

	  QMP_free_msgmem(forw_msg[num][0]);
	  QMP_free_msgmem(forw_msg[num][1]);


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

#ifndef DSLASH_USE_QMT_THREADS
#define smpscaller2(a,bleah,spinor2,half_spinor2,u3,cb2,volume2) \
    a.spinor = spinor2;\
    a.half_spinor = half_spinor2;\
    a.u = u3;  \
    a.cb = cb2; \
    (*bleah)(0, volume2, 0, &a);
#else
#define smpscaller2(a,bleah,spinor2,half_spinor2,u3,cb2,volume2) \
    a.spinor = spinor2;\
    a.half_spinor = half_spinor2;\
    a.u = u3;  \
    a.cb = cb2; \
    qmt_call((qmt_userfunc_t)bleah, volume2, &a);
#endif


void sse_su3dslash_wilson(float *u, float *psi, float *res, int isign, int cb)
{
  int mu;
  Arg_s a;
  int subgrid_vol_cb = getSubgridVolCB();

  if (initP == 0) {
    QMP_error("sse_su3dslash_wilson not initialized");
    QMP_abort(1);
  }

  if(isign==1) 
  {


    smpscaller2(a,decomp_plus,
		(spinor_array*)psi,
		chi1,
		(my_mat_array)u,
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




    smpscaller2(a,decomp_hvv_plus,
		(spinor_array*)psi,
		chi2,
		(my_mat_array)u,
		cb,
		subgrid_vol_cb);
	


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
		(spinor_array*)res,
		chi1,
		(my_mat_array)u,
		1-cb,
		subgrid_vol_cb);



    if (total_comm > 0)
      if (QMP_wait(back_all_mh) != QMP_SUCCESS)
      {
	QMP_error("wnxtsu3dslash: QMP_wait failed in backward direction");
	QMP_abort(1);
      }



    smpscaller2(a,recons_plus,
		(spinor_array*)res, 
		chi2,
		(my_mat_array)u,	
		1-cb,
		subgrid_vol_cb);
  }		

  if(isign==-1) 
  {
    smpscaller2(a,decomp_minus,
		(spinor_array*)psi,
		chi1,
		(my_mat_array)u,
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


    smpscaller2(a,decomp_hvv_minus,
		(spinor_array*)psi,
		chi2,
		(my_mat_array)u,
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


    smpscaller2(a,mvv_recons_minus,
		(spinor_array*)res,
		chi1,
		(my_mat_array)u,
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
		(spinor_array*)res, 
		chi2,
		(my_mat_array)u,	
		1-cb,
		subgrid_vol_cb);
  }		
}



#ifdef __cplusplus
}
#endif

