/*******************************************************************************
 * $Id: sse_su3dslash_32bit_parscalar.c,v 1.22 2008-05-12 14:50:10 bjoo Exp $
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


#include <sse_config.h>            
#include <shift_tables_parscalar.h>
#include <stdlib.h>
#include <stdio.h>
#include <qmp.h>        /* QMP for the comms */
#include <sse_align.h>  /* Alignment stuff to ensure 16 byte alignments */
#include <types32.h>    /* Types and prefetch macros */
#include <dispatch_parscalar.h>   /* Threads dispatch definition */
#include <strings.h>
#include "decomp.h"
#include "decomp_hvv.h"
#include "mvv_recons_32bit.h"
#include "recons.h"

#ifdef __cplusplus
extern "C" {
#endif

  extern int subgrid_vol;
  extern int subgrid_vol_cb;

  extern int *site_table;    /* Aligned lookup table */
  extern int *offset_table;

  static inline 
  int halfspinor_buffer_offset(HalfSpinorOffsetType type, int site, int mu)
  {
    return offset_table[mu + 4*( site + subgrid_vol*type) ];
  }

  static int initP=0;



  
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


/****************************isign corresponding to +1  **************************/

/* the basic operation here is load up a spinor, do the spin decomposition, and store the halfspinor
to a lattice temporary */

  void decomp_plus(size_t lo,size_t hi, int id, const void *ptr) /*need to fix decomp_minus */
  {
    int my_node=QMP_get_node_number();

    int ix1, ix2, iz1;                           /* Site index - iz1 used at loop end */
    spinor_array* sp1 ALIGN;                /* Spinor under consideration */
    spinor_array* sp2 ALIGN;
    const ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;          /* Cast the argument */
    halfspinor_array* chi = a->half_spinor; /* needs to be changed to halfspinor_array and be an array*/
    int cb = a->cb;
    int low = cb*subgrid_vol_cb + lo;
    int high = cb*subgrid_vol_cb + hi;


    halfspinor_array* s3 ALIGN;
    halfspinor_array* s4 ALIGN;
    halfspinor_array* s5 ALIGN;
    halfspinor_array* s6 ALIGN;

    
    spinor_array* spinor_field= a->spinor;
    int thissite;


    /************************ loop over all lattice sites *************************/
    thissite = site_table[low];
    sp1=&spinor_field[thissite]; 
    _mm_prefetch(sp1, _MM_HINT_T0);
    
 
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,low,0);
    _mm_prefetch(s3, _MM_HINT_T0);

    s4 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,low,1);
    _mm_prefetch(s4, _MM_HINT_T0);

    s5 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,low,2);
    _mm_prefetch(s5, _MM_HINT_T0);


    s6 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,low,3);
    _mm_prefetch(s6, _MM_HINT_T0);

    for (ix1=low+1; ix1<high;ix1++) {
      thissite=site_table[ix1]; // Next site  
      sp2=&spinor_field[thissite]; // For prefetching
      _mm_prefetch(sp2, _MM_HINT_T0);
  
      decomp_gamma0_minus(sp1[0], *s3);
      s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,0);
      _mm_prefetch(s3, _MM_HINT_T0);

      decomp_gamma1_minus(sp1[0], *s4);
      s4 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,1);
      _mm_prefetch(s4, _MM_HINT_T0);

      decomp_gamma2_minus(sp1[0], *s5);
      s5 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,2);
      _mm_prefetch(s5, _MM_HINT_T0);




      decomp_gamma3_minus(sp1[0], *s6);
      s6 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,3);
      _mm_prefetch(s6, _MM_HINT_T0);
      
      sp1=sp2; // For prefetching
    }

    decomp_gamma0_minus(sp1[0], *s3);
    decomp_gamma1_minus(sp1[0], *s4);
    decomp_gamma2_minus(sp1[0], *s5);
    decomp_gamma3_minus(sp1[0], *s6);


  }


/* the basic operations in this routine include loading a spinor, doing 
 * the spin projection, and multiplying the halfspinor by the appropriate 
 * gauge field, and saving the resulting halfspinor to a lattice temporary */

/* need gauge fields on opposite cb */
void decomp_hvv_plus(size_t lo,size_t hi, int id, const void *ptr)
{

  int ix1, iz1;              /* Site addresses. ix1 = current. 
				iz1 is for next loop iteration to allow some loop peeling
			        with ix1 */

  u_mat_array* um1 ALIGN;    /* Gauge pointer for 1 site */
  u_mat_array* um2 ALIGN;    /* Gauge pointer for 2nd site */
  u_mat_array* um3 ALIGN;    /* Temporary gauge pointer for prefetching */
  u_mat_array* um4 ALIGN;

  spinor_array* sm1 ALIGN;   /* spinor */
  spinor_array* sm2 ALIGN;

  const ThreadWorkerArgs *a = (const ThreadWorkerArgs *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* chi = a->half_spinor; /* a 1-d map of a 2-d array */
  my_mat_array gauge_field = a->u;

  halfspinor_array* s3 ALIGN;
  halfspinor_array* s4 ALIGN;
  halfspinor_array* s5 ALIGN;
  halfspinor_array* s6 ALIGN;

  int cb = a->cb;

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;


  /************************ loop over all lattice sites *************************/
  int thissite = site_table[low];
  s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,low,0);
  s4 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,low,1);    
  s5 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,low,2);    
  s6 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,low,3); 
  um1=&gauge_field[thissite][0];
  um2=&gauge_field[thissite][1];
  um3=&gauge_field[thissite][2];
  um4=&gauge_field[thissite][3];

  _mm_prefetch(s3, _MM_HINT_T0);
  _mm_prefetch(um1,_MM_HINT_T0);

  _mm_prefetch(s4, _MM_HINT_T0);
  _mm_prefetch(um2,_MM_HINT_T0);

  _mm_prefetch(s5, _MM_HINT_T0);
  _mm_prefetch(um3, _MM_HINT_T0);

  _mm_prefetch(s6, _MM_HINT_T0);
  _mm_prefetch(um4, _MM_HINT_T0);

  sm1=&spinor_field[thissite];

  for (ix1=low+1;ix1<high;ix1++) {
    thissite=site_table[ix1]; // Next site
    sm2=&spinor_field[thissite]; 


    /****************** direction +0 *********************************/
    decomp_hvv_gamma0_plus(*sm1,*um1,*s3);
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,0);
    um1=&gauge_field[thissite][0];
    _mm_prefetch(s3, _MM_HINT_T0);
    _mm_prefetch(um1,_MM_HINT_T0);

    decomp_hvv_gamma1_plus(*sm1,*um2,*s4);
    s4 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,1);
    um2=&gauge_field[thissite][1];
    _mm_prefetch(s4, _MM_HINT_T0);
    _mm_prefetch(um2, _MM_HINT_T0);

    _mm_prefetch(sm1,_MM_HINT_T0);

    decomp_hvv_gamma2_plus(*sm1,*um3,*s5);
    s5 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,2);
    um3=&gauge_field[thissite][2];
    _mm_prefetch(s5, _MM_HINT_T0);
    _mm_prefetch(um3, _MM_HINT_T0);


    decomp_hvv_gamma3_plus(*sm1,*um4,*s6);
    s6 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,3);
    um4=&gauge_field[thissite][3];
    _mm_prefetch(s6, _MM_HINT_T0);
    _mm_prefetch(um4, _MM_HINT_T0);
    
    sm1=sm2;
  }
  decomp_hvv_gamma0_plus(*sm1,*um1,*s3);
  decomp_hvv_gamma1_plus(*sm1,*um2,*s4);

  _mm_prefetch(sm1,_MM_HINT_T0);

  decomp_hvv_gamma2_plus(*sm1,*um3,*s5);
  decomp_hvv_gamma3_plus(*sm1,*um4,*s6);

}
/***************end of decomp_hvv****************/


/* the basic operations in this routine include loading the halfspinor 
 * from memory, multiplying it by the appropriate gauge field, doing the 
 * spin reconstruction, and summing over directions, and saving the partial 
 * sum over directions */

void mvv_recons_plus(size_t lo,size_t hi, int id, const void *ptr)
{

  int ix1, iz1;

  u_mat_array* up1 ALIGN;
  u_mat_array* up2 ALIGN;
  u_mat_array* up3 ALIGN;
  u_mat_array* up4 ALIGN;

  spinor_array* sn1 ALIGN;
  spinor_array* sn2 ALIGN;

  halfspinor_array r12_1 ALIGN, r34_1 ALIGN,r12_2 ALIGN,r34_2 ALIGN;


  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* chi = a->half_spinor; /* a 1-d map of a 2-d array */
  my_mat_array gauge_field = a->u;
  int cb = a->cb;


  halfspinor_array* s3 ALIGN;
  halfspinor_array* s4 ALIGN;
  halfspinor_array* s5 ALIGN;
  halfspinor_array* s6 ALIGN;

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;
  
  /************************ loop over all lattice sites *************************/
  int thissite = site_table[ low ];
  s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,low,0);
  s4 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,low,1);
  s5 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,low,2);
  s6 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,low,3);
  up1=&gauge_field[thissite][0]; 
  up2=&gauge_field[thissite][1];
  up3=&gauge_field[thissite][2];
  up4=&gauge_field[thissite][3];

  _mm_prefetch(s3, _MM_HINT_T0);
  _mm_prefetch(up1,_MM_HINT_T0);

  _mm_prefetch(s4, _MM_HINT_T0);
  _mm_prefetch(up2,_MM_HINT_T0);

  _mm_prefetch(s5, _MM_HINT_T0);
  _mm_prefetch(up3, _MM_HINT_T0);

  _mm_prefetch(s6, _MM_HINT_T0);
  _mm_prefetch(up4, _MM_HINT_T0);

  sn1=&spinor_field[thissite];    


  for (ix1=low+1;ix1<high;ix1++) {

    thissite=site_table[ix1];
    sn2=&spinor_field[thissite];    

    mvv_recons_gamma0_plus(*s3, *up1, r12_1, r34_1);
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,0);
    up1=&gauge_field[thissite][0]; 
    _mm_prefetch(s3, _MM_HINT_T0);
    _mm_prefetch(up1,_MM_HINT_T0);

    mvv_recons_gamma1_plus_add(*s4, *up2, r12_1, r34_1);
    s4 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,1);
    up2=&gauge_field[thissite][1]; 
    _mm_prefetch(s4, _MM_HINT_T0);
    _mm_prefetch(up2, _MM_HINT_T0);

    
    _mm_prefetch(sn1,_MM_HINT_T0);

    mvv_recons_gamma2_plus_add(*s5, *up3, r12_1, r34_1);
    s5 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,2);
    up3=&gauge_field[thissite][2];
    _mm_prefetch(s5, _MM_HINT_T0);
    _mm_prefetch(up3, _MM_HINT_T0);

    mvv_recons_gamma3_plus_add_store(*s6, *up4, r12_1, r34_1,*sn1);
    s6 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,3);
    up4=&gauge_field[thissite][3];
    _mm_prefetch(s6, _MM_HINT_T0);
    _mm_prefetch(up4, _MM_HINT_T0);
  
    sn1=sn2;
  }

  mvv_recons_gamma0_plus(*s3, *up1, r12_1, r34_1);
  mvv_recons_gamma1_plus_add(*s4, *up2, r12_1, r34_1);

  _mm_prefetch(sn1,_MM_HINT_T0);

  mvv_recons_gamma2_plus_add(*s5, *up3, r12_1, r34_1);
  mvv_recons_gamma3_plus_add_store(*s6, *up4, r12_1, r34_1,*sn1);

}


   
/* this routine takes the partial sum from mvv_recons() and loops 
 * over the output spin components, 2 at a time doing a sum over directions 
 * for each set, accumulating in xmm0-2 and loading the halfspinor 
 * temporaries into xmm3-5 */

void recons_plus(size_t lo,size_t hi, int id, const void *ptr )	
{
  int ix1;
  spinor_array* sn1 ALIGN;
  spinor_array* sn2 ALIGN;
  

  const ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* chi = a->half_spinor;
  int cb = a->cb;
 
  halfspinor_array *hs0 ALIGN;
  halfspinor_array *hs1 ALIGN;
  halfspinor_array *hs2 ALIGN;
  halfspinor_array *hs3 ALIGN;

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;
  
  const int PREFDIST=4;

  /************************ loop over all lattice sites *************************/
  int thissite = site_table[low];



  hs0 = chi + halfspinor_buffer_offset(RECONS_GATHER,low,0); 
  _mm_prefetch(hs0, _MM_HINT_NTA);

  hs1 = chi + halfspinor_buffer_offset(RECONS_GATHER,low,1); 
  _mm_prefetch(hs1, _MM_HINT_NTA);

  hs2 = chi + halfspinor_buffer_offset(RECONS_GATHER,low,2); 
  _mm_prefetch(hs2, _MM_HINT_NTA);

  hs3 = chi + halfspinor_buffer_offset(RECONS_GATHER,low,3);
  _mm_prefetch(hs3, _MM_HINT_NTA);

  sn1=&spinor_field[thissite];   
  _mm_prefetch(sn1, _MM_HINT_NTA);

  for (ix1=low+1;ix1<high;ix1++) {
    thissite = site_table[ix1];
    sn2 = &spinor_field[thissite];   
    _mm_prefetch(sn2, _MM_HINT_NTA);
   
    recons_4dir_plus(*hs0, *hs1, *hs2, *hs3, *sn1);

    hs0 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,0); 
    _mm_prefetch(hs0, _MM_HINT_NTA);

    hs1 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,1); 
    _mm_prefetch(hs1, _MM_HINT_NTA);

    hs2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,2); 
    _mm_prefetch(hs2, _MM_HINT_NTA);

    hs3 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,3); 
    _mm_prefetch(hs3, _MM_HINT_NTA);

    sn1=sn2;
    /*************************end of loop ****************************/
  }
  recons_4dir_plus(*hs0, *hs1, *hs2, *hs3, *sn1);


}
/*****************end of recons**************/





/*************** now for isign corresponding to -1  ****************************************/

void decomp_minus(size_t lo,size_t hi, int id, const void *ptr ) /*need to fix decomp_minus */
{

  int ix1,iz1;

  spinor_array* sp1 ALIGN;
  spinor_array* sp2 ALIGN;

  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;

  halfspinor_array* chi = a->half_spinor; /* needs to be changed to halfspinor_array and be an array*/

  halfspinor_array* s3 ALIGN;
  halfspinor_array* s4 ALIGN;
  halfspinor_array* s5 ALIGN;
  halfspinor_array* s6 ALIGN;

  int cb = a->cb;

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;

  spinor_array* spinor_field= a->spinor;
  int thissite;

  /************************ loop over all lattice sites *************************/

  thissite = site_table[low];
  sp1=&spinor_field[thissite]; 
  _mm_prefetch(sp1, _MM_HINT_T0);
  
  s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,low,0);
  _mm_prefetch(s3, _MM_HINT_T0);
  
  s4 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,low,1);
  _mm_prefetch(s4, _MM_HINT_T0);
  
  s5 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,low,2);
  _mm_prefetch(s5, _MM_HINT_T0);
  
  s6 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,low,3);
  _mm_prefetch(s6, _MM_HINT_T0);

  
  for (ix1=low+1;ix1<high;ix1++) {
    thissite=site_table[ix1]; // Next site  
    sp2=&spinor_field[thissite]; // For prefetching
    _mm_prefetch(sp2, _MM_HINT_T0);
  
    
    decomp_gamma0_plus(sp1[0], *s3);
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,0);
    _mm_prefetch(s3, _MM_HINT_T0);

    decomp_gamma1_plus(sp1[0], *s4);
    s4 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,1);
    _mm_prefetch(s4, _MM_HINT_T0);
      
    decomp_gamma2_plus(sp1[0], *s5);
    s5 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,2);
    _mm_prefetch(s5, _MM_HINT_T0);

    decomp_gamma3_plus(sp1[0], *s6);    
    s6 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,3);
    _mm_prefetch(s6, _MM_HINT_T0);
      
    sp1=sp2; // For prefetching
  }
  decomp_gamma0_plus(sp1[0], *s3);
  decomp_gamma1_plus(sp1[0], *s4);
  decomp_gamma2_plus(sp1[0], *s5);
  decomp_gamma3_plus(sp1[0], *s6);
  
}


/* need gauge fields on opposite cb */
void decomp_hvv_minus(size_t lo,size_t hi, int id, const void *ptr )
{

  int ix1,iz1;
  u_mat_array* um1 ALIGN;
  u_mat_array* um2 ALIGN;
  u_mat_array* um3 ALIGN;
  u_mat_array* um4 ALIGN;

  spinor_array* sm1 ALIGN;
  spinor_array* sm2 ALIGN;



  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* chi = a->half_spinor; /* a 1-d map of a 2-d array */
  my_mat_array gauge_field = a->u;
  int cb = a->cb;

  halfspinor_array* s3 ALIGN;
  halfspinor_array* s4 ALIGN;
  halfspinor_array* s5 ALIGN;
  halfspinor_array* s6 ALIGN;

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;

  /************************ loop over all lattice sites *************************/
  int thissite = site_table[low];
  s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,low,0);
  s4 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,low,1);    
  s5 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,low,2);    
  s6 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,low,3); 
  um1=&gauge_field[thissite][0];
  um2=&gauge_field[thissite][1];
  um3=&gauge_field[thissite][2];
  um4=&gauge_field[thissite][3];

  _mm_prefetch(s3, _MM_HINT_T0);
  _mm_prefetch(um1,_MM_HINT_T0);

  _mm_prefetch(s4, _MM_HINT_T0);
  _mm_prefetch(um2,_MM_HINT_T0);

  _mm_prefetch(s5, _MM_HINT_T0);
  _mm_prefetch(um3, _MM_HINT_T0);

  _mm_prefetch(s6, _MM_HINT_T0);
  _mm_prefetch(um4, _MM_HINT_T0);


  sm1 = &spinor_field[thissite];

  for (ix1=low+1;ix1<high;ix1++) {
    thissite = site_table[ix1];
    sm2=&spinor_field[thissite]; 

    /***************** direction +0 *********************************/
    decomp_hvv_gamma0_minus(*sm1, *um1, *s3);
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,0);
    um1=&gauge_field[thissite][0];
    _mm_prefetch(s3, _MM_HINT_T0);
    _mm_prefetch(um1,_MM_HINT_T0);

    decomp_hvv_gamma1_minus(*sm1, *um2, *s4);
    s4 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,1);
    um2=&gauge_field[thissite][1];
    _mm_prefetch(s4, _MM_HINT_T0);
    _mm_prefetch(um2, _MM_HINT_T0);


    decomp_hvv_gamma2_minus(*sm1, *um3, *s5);
    s5 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,2);
    um3=&gauge_field[thissite][2];
    _mm_prefetch(s5, _MM_HINT_T0);
    _mm_prefetch(um3, _MM_HINT_T0);

    _mm_prefetch(sm1,_MM_HINT_T0);
    
    decomp_hvv_gamma3_minus(*sm1, *um4, *s6);
    s6 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,3);
    um4=&gauge_field[thissite][3];

    _mm_prefetch(s6, _MM_HINT_T0);
    _mm_prefetch(um4, _MM_HINT_T0);
  
    sm1=sm2;

  }
  decomp_hvv_gamma0_minus(*sm1, *um1, *s3);
  decomp_hvv_gamma1_minus(*sm1, *um2, *s4);

  _mm_prefetch(sm1,_MM_HINT_T0);

  decomp_hvv_gamma2_minus(*sm1, *um3, *s5);
  decomp_hvv_gamma3_minus(*sm1, *um4, *s6);

}


void mvv_recons_minus(size_t lo,size_t hi, int id, const void *ptr )
{
  int ix1, iz1;
  spinor_array* sn1 ALIGN;  /* The spinor to store to */
  spinor_array* sn2 ALIGN; 
  u_mat_array* um1 ALIGN;
  u_mat_array* um2 ALIGN;
  u_mat_array* um3 ALIGN;
  u_mat_array* um4 ALIGN;

  
  /* Temporaries for the top and bottom parts of spinors. */
  halfspinor_array r12_1 ALIGN, r34_1 ALIGN, r12_2 ALIGN,r34_2 ALIGN;

  /* if going to support unpacked gauge fields, need to treat site ix1 and site ix1+1 separately */
  /* to support unpacked gauge fields the prefetches will need to be changed */
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* chi = a->half_spinor; /* a 1-d map of a 2-d array */
  my_mat_array gauge_field = a->u;
  int cb = a->cb;

  halfspinor_array* s3 ALIGN;
  halfspinor_array* s4 ALIGN;
  halfspinor_array* s5 ALIGN;
  halfspinor_array* s6 ALIGN;
  

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;

  int thissite = site_table[low];
  s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,low,0);
  s4 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,low,1);
  s5 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,low,2);
  s6 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,low,3);
  um1=&gauge_field[thissite][0]; 
  um2=&gauge_field[thissite][1];
  um3=&gauge_field[thissite][2];
  um4=&gauge_field[thissite][3];

  _mm_prefetch(s3, _MM_HINT_T0);
  _mm_prefetch(um1,_MM_HINT_T0);

  _mm_prefetch(s4, _MM_HINT_T0);
  _mm_prefetch(um2,_MM_HINT_T0);

  _mm_prefetch(s5, _MM_HINT_T0);
  _mm_prefetch(um3, _MM_HINT_T0);

  _mm_prefetch(s6, _MM_HINT_T0);
  _mm_prefetch(um4, _MM_HINT_T0);

  sn1=&spinor_field[thissite];    
/************************ loop over all lattice sites *************************/
  for (ix1=low+1;ix1<high;ix1++) {
    thissite = site_table[ix1];
    sn2 = &spinor_field[thissite];   
 
    mvv_recons_gamma0_minus(*s3, *um1, r12_1, r34_1);
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,0);
    um1=&gauge_field[thissite][0];  
    _mm_prefetch(s3, _MM_HINT_T0);
    _mm_prefetch(um1,_MM_HINT_T0);

    mvv_recons_gamma1_minus_add(*s4, *um2, r12_1, r34_1);
    s4 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,1); 
    um2= &gauge_field[thissite][1];
    _mm_prefetch(s4, _MM_HINT_T0);
    _mm_prefetch(um2, _MM_HINT_T0);


    mvv_recons_gamma2_minus_add(*s5, *um3, r12_1, r34_1);
    s5 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,2);
    um3=&gauge_field[thissite][2];  
    _mm_prefetch(s5, _MM_HINT_T0);
    _mm_prefetch(um3, _MM_HINT_T0);

    _mm_prefetch(sn1,_MM_HINT_T0);

    mvv_recons_gamma3_minus_add_store(*s6, *um4, r12_1, r34_1,*sn1);

    s6 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,3);
    um4=&gauge_field[thissite][3];  
    _mm_prefetch(s6, _MM_HINT_T0);
    _mm_prefetch(um4, _MM_HINT_T0);
  
    sn1=sn2;
    /******************************** end of loop *********************************/
  }
  mvv_recons_gamma0_minus(*s3, *um1, r12_1, r34_1);
  mvv_recons_gamma1_minus_add(*s4, *um2, r12_1, r34_1);
  mvv_recons_gamma2_minus_add(*s5, *um3, r12_1, r34_1);
  _mm_prefetch(sn1,_MM_HINT_T0);
  mvv_recons_gamma3_minus_add_store(*s6, *um4, r12_1, r34_1,*sn1);
}
/******************end of mvv_recons*************************/


void recons_minus(size_t lo,size_t hi, int id, const void *ptr )	
{
  int ix1;
  spinor_array* sn1 ALIGN;
  spinor_array* sn2 ALIGN;


  const ThreadWorkerArgs *a = (ThreadWorkerArgs *)ptr;
  spinor_array* spinor_field = a->spinor;
  halfspinor_array* chi = a->half_spinor; /* a 1-d map of a 2-d array */
  int cb = a->cb;

  halfspinor_array *hs0 ALIGN;
  halfspinor_array *hs1 ALIGN;
  halfspinor_array *hs2 ALIGN;
  halfspinor_array *hs3 ALIGN;

  int low = cb*subgrid_vol_cb + lo;
  int high = cb*subgrid_vol_cb + hi;
  int thissite = site_table[low];  
  hs0 = chi + halfspinor_buffer_offset(RECONS_GATHER,low,0); 
  _mm_prefetch(hs0, _MM_HINT_NTA);

  hs1 = chi + halfspinor_buffer_offset(RECONS_GATHER,low,1); 
  _mm_prefetch(hs1, _MM_HINT_NTA);

  hs2 = chi + halfspinor_buffer_offset(RECONS_GATHER,low,2); 
  _mm_prefetch(hs2, _MM_HINT_NTA);

  hs3 = chi + halfspinor_buffer_offset(RECONS_GATHER,low,3);
  _mm_prefetch(hs3, _MM_HINT_NTA);

  sn1=&spinor_field[thissite];   
  _mm_prefetch(sn1, _MM_HINT_NTA);

  for (ix1=low+1;ix1<high;ix1++) {
     thissite = site_table[ix1];
     sn2 = &spinor_field[thissite];   
    _mm_prefetch(sn2, _MM_HINT_NTA);

    recons_4dir_minus(*hs0, *hs1, *hs2, *hs3, *sn1);

    hs0 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,0); 
    _mm_prefetch(hs0, _MM_HINT_NTA);
    
    hs1 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,1); 
    _mm_prefetch(hs1, _MM_HINT_NTA);
    
    hs2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,2); 
    _mm_prefetch(hs2, _MM_HINT_NTA);
    
    hs3 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,3); 
    _mm_prefetch(hs3, _MM_HINT_NTA);
   
    sn1=sn2;

  }
  recons_4dir_minus(*hs0, *hs1, *hs2, *hs3, *sn1);
}
/*****************end of isign corresponding to -1 **********************/


/***************** start of initialization routine ***************************************/


static QMP_mem_t* xchi1;               /* QMP Memory Structures for halfspinor arrays */
static QMP_mem_t* xchi2;               /* xchi1 => FORWARD, xchi2 => BACKWARD         */

#if 0
static halfspinor_array* chi1;         /* These are the aligned pointers from the QMP Memory structures */
static halfspinor_array* chi2;         /* xchi1 <=> chi1    xchi2 <=> chi2 */
#else
  /* Make these visible for testing/timing */
  halfspinor_array* chi1;         /* These are the aligned pointers from the QMP Memory structures */
  halfspinor_array* chi2;         /* xchi1 <=> chi1    xchi2 <=> chi2 */
#endif


/* Nearest neighbor communication channels */
static int total_comm = 0;
static QMP_msgmem_t forw_msg[4][2];
static QMP_msgmem_t back_msg[4][2];
static QMP_msghandle_t forw_mh_send[4]; /* Send forward */
static QMP_msghandle_t forw_mh_recv[4]; /* Receive from backwards */
static QMP_msghandle_t back_mh_send[4]; /* Send backwards */
static QMP_msghandle_t back_mh_recv[4]; /* Receive from forwards */
static QMP_msghandle_t forw_all_mh_send;
static QMP_msghandle_t forw_all_mh_recv;

static QMP_msghandle_t back_all_mh_send;
static QMP_msghandle_t back_all_mh_recv;
static int recvPostedP=0;



/* Initialize the Dslash */
  void init_sse_su3dslash(const int latt_size[],
			  void (*getSiteCoords)(int coord[], int node, int linearsite),
			  int (*getLinearSiteIndex)(const int coord[]),
			  int (*nodeNumber)(const int coord[]))   // latt_size not used, here for scalar version
{
  const int *machine_size = QMP_get_logical_dimensions();
  // const int *subgrid_cb_size = QMP_get_subgrid_dimensions();
 
  /* bound[cb][forward=0/backward=1][dir] */
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
    
  /* Check problem size - 4D  */
  for(mu=0; mu < 4; mu++)  {
    if ( latt_size[mu] % 2 != 0 ) {
      fprintf(stderr,"This is a Dslash with checkerboarding in 4 dimensions. Each GLOBAL dimension must be even. In addition LOCAL dimension 0 (x) has to be even ,  Your lattice does not meet the GLOBAL requirement latt_size[%d]=%d\n", 
	      mu, latt_size[mu]);
      
      exit(1);
    }
  }
  
  num = latt_size[0] / machine_size[0];
  if ( num % 2 != 0 )
    {
      fprintf(stderr,"This is a Dslash with checkerboarding in 4 dimensions. Each GLOBAL dimension must be even. In addition LOCAL dimension 0 (x) has to be even ,  Your lattice does not meet the LOCAL requirement sublattice_size[0]=%d\n", 
	      num);
    QMP_abort(1);
  }


  /* Make the shift table  -- this sets the vol and vol_cb so we can call getSugridVolCB() after it */
  make_shift_tables(bound,
		    getSiteCoords,
		    getLinearSiteIndex,
		    nodeNumber);


  /* Allocate the halfspinor temporaries */
  /* The halfspinor is really the following:
     
        halfspinor[0][Nd][vol_cb] - half spinors from the body
	halfspinor[1][Nd][vol_cb] - half spinors from the tail to be sent forward
	halfspinor[2][Nd][vol_cb] - half spinors from the tail to be sent backward

	The Nd stands for the 4 decomposition directions (gamma_matrices)

	The body and tails are overallocated to their maximum size: vol_cb
	Hence the size is half_spinor_size * vol_cb*3*Nd
	where Nd=4, 3 is the number of buffers (body, forward tail, backward tail)
	and   half spinor size is 2spin comp * 3 color * 2 floats (=1complex)

  */
  nsize = 2*3*2*sizeof(float)*subgrid_vol_cb*3*4;  /* Note 3x4 half-fermions */

  /* xchi1 and xchi2 will hold projected spinors memory handles from QMP */
  if ((xchi1 = QMP_allocate_aligned_memory(nsize,64,0)) == 0)
  {
    QMP_error("init_wnxtsu3dslash: could not initialize xchi1");
    QMP_abort(1);
  }

  if ((xchi2 = QMP_allocate_aligned_memory(nsize,64,0)) == 0)
  {
    QMP_error("init_wnxtsu3dslash: could not initialize xchi2");
    QMP_abort(1);
  }

  /* Unwrap the half spinor pointers from the QMP Structures. BTW: This 2 step technique susks so bad! */
  chi1 = (halfspinor_array*)QMP_get_memory_pointer(xchi1);
  chi2 = (halfspinor_array*)QMP_get_memory_pointer(xchi2); 


  /* Zero these out for testing */
  bzero(chi1, nsize);
  bzero(chi2, nsize);

  /* Loop over all communicating directions and build up the two message
   * handles. If there is no communications, the message handles will not
   * be initialized 
   */
  num = 0;

  /* Loop over directions */
  for(mu=0; mu < 4; ++mu) {

    if(machine_size[mu] > 1) { /* If the machine is not a scalar  in this dimensio */
    
      if (bound[0][0][mu] == 0) { /* Consistency: Check the boundary in this direction is 0 */
 	QMP_error("init_sse_dslash: type 0 message size is 0: mu=%d", mu);
	QMP_abort(1);
      }

      /* 
	 Boundary is indexed as            boundary[cb][forward=0/back=1][direction] 
      */	 

      /* cb = 0 */
      forw_msg[num][0] = QMP_declare_msgmem( chi1+subgrid_vol_cb*(1+3*mu),
					     bound[0][0][mu]*sizeof(halfspinor_array));

      /* cb = 1 */
      forw_msg[num][1] = QMP_declare_msgmem( chi1+subgrid_vol_cb*(2+3*mu),
					     bound[1][0][mu]*sizeof(halfspinor_array));

      /* cb = 0: Receive from cb = 1 */
      forw_mh_recv[num]  = QMP_declare_receive_relative(forw_msg[num][1], mu, +1, 0);

      /* cb = 1: send to cb = 0 */
      forw_mh_send[num]  = QMP_declare_send_relative(forw_msg[num][0], mu, -1, 0);
	
      if (bound[0][1][mu] == 0)
      {
	QMP_error("init_sse_dslash: type 0 message size is 0");
	QMP_abort(1);
      }

      /* cb = 0 */
      back_msg[num][0] = QMP_declare_msgmem( chi2+subgrid_vol_cb*(1+3*mu),
					     bound[0][1][mu]*sizeof(halfspinor_array));

      /* cb = 1 */
      back_msg[num][1] = QMP_declare_msgmem( chi2+subgrid_vol_cb*(2+3*mu), 
					     bound[1][1][mu]*sizeof(halfspinor_array));

      /* cb = 0: Receive from cb=1 */
      back_mh_recv[num]  = QMP_declare_receive_relative(back_msg[num][1], mu, -1, 0);

      /* cb = 1: Receive from cb=0 */
      back_mh_send[num]  = QMP_declare_send_relative(back_msg[num][0], mu, +1, 0);
	
      num++;
    }
  }

  /* Combine the messages */
  if (num > 0) {
    forw_all_mh_send = QMP_declare_multiple(&(forw_mh_send[0]), num);
    forw_all_mh_recv = QMP_declare_multiple(&(forw_mh_recv[0]), num);

    back_all_mh_send = QMP_declare_multiple(&(back_mh_send[0]), num);
    back_all_mh_recv = QMP_declare_multiple(&(back_mh_recv[0]), num);
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

    /* Free shift table */
    free_shift_tables();
    

    if (total_comm > 0) {
      
      QMP_free_msghandle(forw_all_mh_send);
      QMP_free_msghandle(forw_all_mh_recv);

      QMP_free_msghandle(back_all_mh_send);
      QMP_free_msghandle(back_all_mh_recv);
  
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

void sse_su3dslash_prepost_receives(void) 
{
  /* Prepost all receives */
  if (total_comm > 0 && recvPostedP==0) {

    if (QMP_start(forw_all_mh_recv) != QMP_SUCCESS) {
      QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
      QMP_abort(1);
    }
    
    if (QMP_start(back_all_mh_recv) != QMP_SUCCESS) {
      QMP_error("sse_su3dslash_wilson: QMP_start failed in backward direction");
      QMP_abort(1);
    }
  
    recvPostedP=1;
  }

}

void sse_su3dslash_wilson(float *u, float *psi, float *res, int isign, int cb)
{
  

  if (initP == 0) {
    QMP_error("sse_su3dslash_wilson not initialized");
    QMP_abort(1);
  }

  if(isign==1) 
  {

    sse_su3dslash_prepost_receives();

    dispatch_to_threads(decomp_plus,
		(spinor_array*)psi,
		chi1,
		(my_mat_array)u,
		cb,
		subgrid_vol_cb);


    if(total_comm > 0) {

      if (QMP_start(forw_all_mh_send) != QMP_SUCCESS) {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	QMP_abort(1);
      }
    }

    dispatch_to_threads(decomp_hvv_plus,
		(spinor_array*)psi,
		chi2,
		(my_mat_array)u,
		cb,
		subgrid_vol_cb);
	

    /* Wait for forward comms to finish */
    if (total_comm > 0) {

      /* Finish all sends */
      if (QMP_wait(forw_all_mh_send) != QMP_SUCCESS) {
	QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	QMP_abort(1);
      }


      /* Finish all forward receives */
      if (QMP_wait(forw_all_mh_recv) != QMP_SUCCESS) {
	QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	QMP_abort(1);
      }

      /* Start back sends - receives should be preposted */

      if (QMP_start(back_all_mh_send) != QMP_SUCCESS) {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in backward direction");
	QMP_abort(1);
      }
    }

    dispatch_to_threads(mvv_recons_plus,
		(spinor_array*)res,
		chi1,
		(my_mat_array)u,
		1-cb,
		subgrid_vol_cb);


    /* Wait for back comms to complete */
    if (total_comm > 0) {

      if (QMP_wait(back_all_mh_send) != QMP_SUCCESS) {
	QMP_error("wnxtsu3dslash: QMP_wait failed in backward direction");
	QMP_abort(1);
      }


      if (QMP_wait(back_all_mh_recv) != QMP_SUCCESS) {
	QMP_error("wnxtsu3dslash: QMP_wait failed in backward direction");
	QMP_abort(1);
      }

    }




    dispatch_to_threads(recons_plus,
		(spinor_array*)res, 
		chi2,
		(my_mat_array)u,	
		1-cb,
		subgrid_vol_cb);


  }		

  if(isign==-1) 
  {

    sse_su3dslash_prepost_receives();


    dispatch_to_threads(decomp_minus,
		(spinor_array*)psi,
		chi1,
		(my_mat_array)u,
		cb,
		subgrid_vol_cb);


    /* Start sends */
    if (total_comm > 0) {
      if (QMP_start(forw_all_mh_send) != QMP_SUCCESS) {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	QMP_abort(1);
      }
    }

    dispatch_to_threads(decomp_hvv_minus,
		(spinor_array*)psi,
		chi2,
		(my_mat_array)u,
		cb,
		subgrid_vol_cb);


    /* Finish forward comms */
    if (total_comm > 0) {
      if (QMP_wait(forw_all_mh_send) != QMP_SUCCESS) {
	QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	QMP_abort(1);
      }

      if (QMP_wait(forw_all_mh_recv) != QMP_SUCCESS) {
	QMP_error("sse_su3dslash_wilson: QMP_wait failed in forward direction");
	QMP_abort(1);
      }
    }


    /* Start backward comms  - receive ought to be preposted */
    if (total_comm > 0) {
      if (QMP_start(back_all_mh_send) != QMP_SUCCESS) {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in backward direction");
	QMP_abort(1);
      }
    }


    dispatch_to_threads(mvv_recons_minus,
		(spinor_array*)res,
		chi1,
		(my_mat_array)u,
		1-cb,
		subgrid_vol_cb);
    

    /* Wait for backward comms to complete */
    if (total_comm > 0) {
      if (QMP_wait(back_all_mh_send) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_wait failed in backward direction");
	QMP_abort(1);
      }


      if (QMP_wait(back_all_mh_recv) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_wait failed in backward direction");
	QMP_abort(1);
      }

    }


    dispatch_to_threads(recons_minus,
		(spinor_array*)res, 
		chi2,
		(my_mat_array)u,	
		1-cb,
		subgrid_vol_cb);
  }		

  /* Clear this for next iteration */
  recvPostedP = 0;
}

#ifdef __cplusplus
}
#endif

