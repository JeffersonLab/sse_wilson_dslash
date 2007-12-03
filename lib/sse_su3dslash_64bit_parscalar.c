/*******************************************************************************
 * $Id: sse_su3dslash_64bit_parscalar.c,v 1.10 2007-12-03 16:36:56 bjoo Exp $
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
#include <types64.h>
#include <sse_align.h>
#include <dispatch_parscalar.h>
#include <shift_tables_parscalar.h>
#include <decomp.h>
#include <decomp_hvv.h>
#include <recons.h>
#include <mvv_recons_64bit.h>

  extern int subgrid_vol;
  extern int *offset_table;

  static inline
  int halfspinor_buffer_offset(HalfSpinorOffsetType type, int site, int mu)
  {
    return offset_table[mu + 4*( site + subgrid_vol*type) ];
  }

static int initP=0;


/* now overlays for spinors as arrays or structs */

static int icolor_start[2];    /* starting site for each coloring (cb) */
static int icolor_end[2];      /* end site for each coloring (cb) */


  /* this routine is similar to wnxtsu3dslash, except instead of handling the second site's worth in the same loop, the second spin component's worth must be handled seperately */
void decomp_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  int cb = a->cb; 
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  
  int ix1=0;
  int low = icolor_start[cb]+lo; 
  int high = icolor_start[cb]+hi;
  /*  printf("ID: %d low=%d high=%d: DecompPlus\n", id, low, high); */

  halfspinor_array *s3;
  spinor_array* sp ALIGN;

  for (ix1=low;ix1<high;ix1++) {

    sp=&psi[ix1];
    /******************************* direction +0 *********************************/	   
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,0);
    decomp_gamma0_minus(*sp, *s3);
    /******************************* direction +1 *********************************/
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,1);
    decomp_gamma1_minus(*sp, *s3);
    
    /******************************* direction +2 *********************************/
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,2);
    decomp_gamma2_minus(*sp, *s3);

    /******************************* direction +3 *********************************/
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,3);
    decomp_gamma3_minus(*sp, *s3);

  }
}


void decomp_hvv_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  int cb = a->cb; 
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  int ix1=0;
  int low = icolor_start[cb]+lo; 
  int high = icolor_start[cb]+hi;

  u_mat_array *um ALIGN;

  halfspinor_array *s3 ALIGN;

  spinor_array *sm ALIGN; 

  /*  printf("ID: %d low=%d high=%d: DecompHvvPlus\n", id, low, high);*/

  for (ix1=low;ix1<high;ix1++) 
  {
    /* Spinor to project*/
    sm=&psi[ix1];

    /******************************* direction -1 *********************************/
    um=&gauge_field[ix1][0];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,0);
    decomp_hvv_gamma0_plus(*sm, *um, *s3);
    
    /******************************* direction -1 *********************************/
    um=&gauge_field[ix1][1];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,1);
    decomp_hvv_gamma1_plus(*sm, *um, *s3);

    /******************************* direction -2 *********************************/
    um=&gauge_field[ix1][2];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,2);
    decomp_hvv_gamma2_plus(*sm, *um, *s3);

    /******************************* direction -3 *********************************/
    um=&gauge_field[ix1][3];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,3);
    decomp_hvv_gamma3_plus(*sm, *um, *s3);
  }
}


void mvv_recons_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  int cb = a->cb; 

  int ix1=0;

  int low = icolor_start[cb]+lo; 
  int high = icolor_start[cb]+hi;

  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;
  u_mat_array *up ALIGN;
  halfspinor_array *s3 ALIGN, *s4 ALIGN;
  spinor_array part_sum ALIGN, *result ALIGN;

  /* printf("ID: %d low=%d high=%d: MvvReconsPlus\n", id, low, high); */
 


  for (ix1=low;ix1<high;ix1++) {
    result=&psi[ix1];

    up=&gauge_field[ix1][0];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,0);
    mvv_recons_gamma0_plus(*s3, *up, part_sum);    

    up=&gauge_field[ix1][1];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,1);
    mvv_recons_gamma1_plus_add(*s3, *up, part_sum);    

    up=&gauge_field[ix1][2];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,2);
    mvv_recons_gamma2_plus_add(*s3, *up, part_sum);    

    up=&gauge_field[ix1][3];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,3);
    mvv_recons_gamma3_plus_add_store(*s3, *up, part_sum, *result);    

  }
}



/*optimized for SZIN spin basis */
void recons_plus(size_t lo, size_t hi, int id, const void *ptr)
{


  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  int cb = a->cb; 
  int ix1=0;

  int low = icolor_start[cb]+lo; 
  int high = icolor_start[cb]+hi;

  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;

  halfspinor_array *hs0 ALIGN , *hs1 ALIGN, *hs2 ALIGN, *hs3 ALIGN;
  spinor_array *rn ALIGN;

  /*  printf("ID: %d low=%d high=%d: ReconsPlus\n", id, low, high); */
  for (ix1=low;ix1<high;ix1++) 
  {

    rn=&psi[ix1];

    /* first spin component of result */
    hs0 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,0);
    hs1 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,1);	  
    hs2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,2);
    hs3 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,3);
    recons_4dir_plus(*hs0, *hs1, *hs2, *hs3, *rn);
  }
 
}


/************now for isign = -1  **********************/


void decomp_minus(size_t lo, size_t hi, int id, const void *ptr)
{
   

  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  int cb = a->cb; 
  int ix1=0;
  int low = icolor_start[cb]+lo; 
  int high = icolor_start[cb]+hi;

  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;

  halfspinor_array *s3 ALIGN; 
  spinor_array *sp ALIGN;
 
  /*  printf("ID: %d low=%d high=%d: DecompMinus\n", id, low, high); */
  for (ix1=low;ix1<high;ix1++) {

    sp=&psi[ix1];

    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,0);
    decomp_gamma0_plus(*sp, *s3);

    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,1);
    decomp_gamma1_plus(*sp, *s3);
    
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,2);
    decomp_gamma2_plus(*sp, *s3);
    
    s3 = chi + halfspinor_buffer_offset(DECOMP_SCATTER,ix1,3);
    decomp_gamma3_plus(*sp, *s3); 
  
  }
}


void decomp_hvv_minus(size_t lo, size_t hi, int id, const void *ptr)
{

  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  int cb = a->cb; 
  int ix1 = 0;
  int low = icolor_start[cb]+lo; 
  int high = icolor_start[cb]+hi;
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;

  u_mat_array *um ALIGN;
  halfspinor_array *s3 ALIGN;
  spinor_array  *sm ALIGN;

  /* printf("ID: %d low=%d high=%d: DecompHvvMinus\n", id, low, high); */

  for (ix1=low;ix1<high;ix1++) 
  {
    sm=&psi[ix1];

    um=&gauge_field[ix1][0];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,0);
    decomp_hvv_gamma0_minus(*sm, *um, *s3);	   

    um=&gauge_field[ix1][1];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,1);
    decomp_hvv_gamma1_minus(*sm, *um, *s3);	   

    um=&gauge_field[ix1][2];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,2);
    decomp_hvv_gamma2_minus(*sm, *um, *s3);	   
    
    um=&gauge_field[ix1][3];
    s3 = chi + halfspinor_buffer_offset(DECOMP_HVV_SCATTER,ix1,3);
    decomp_hvv_gamma3_minus(*sm, *um, *s3);	   
  }
}


void mvv_recons_minus(size_t lo, size_t hi, int id, const void *ptr)
{
   	
  
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  int cb = a->cb; 
  int ix1;
  int low = icolor_start[cb]+lo; 
  int high = icolor_start[cb]+hi;
  u_mat_array (*gauge_field)[4] = a->u;
  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;

  u_mat_array *up ALIGN;

  halfspinor_array *s3 ALIGN;

  spinor_array rs ALIGN;
  spinor_array *rn ALIGN;

  /*  printf("ID: %d low=%d high=%d: MvvReconsMinus\n", id, low, high); */

  for (ix1=low;ix1<high;ix1++) {
    rn=&psi[ix1];

    up = &gauge_field[ix1][0];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,0);
    mvv_recons_gamma0_minus(*s3, *up, rs);


    up = &gauge_field[ix1][1];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,1);
    mvv_recons_gamma1_minus_add(*s3, *up, rs);


    up = &gauge_field[ix1][2];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,2);
     mvv_recons_gamma2_minus_add(*s3, *up, rs);

    up = &gauge_field[ix1][3];
    s3 = chi + halfspinor_buffer_offset(RECONS_MVV_GATHER,ix1,3);
    mvv_recons_gamma3_minus_add_store(*s3, *up, rs,*rn);

  }

}



void recons_minus(size_t lo, size_t hi, int id, const void *ptr)
{
  
  const ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  int cb = a->cb; 
  int ix1=0;

  int low = icolor_start[cb]+lo; 
  int high = icolor_start[cb]+hi;

  spinor_array *psi = a->spinor;
  halfspinor_array *chi = a->half_spinor;

  halfspinor_array *hs0 ALIGN;
  halfspinor_array *hs1 ALIGN;
  halfspinor_array *hs2 ALIGN;
  halfspinor_array *hs3 ALIGN;   
  spinor_array  *s1 ALIGN,  *rn ALIGN;

  for (ix1=low; ix1<high; ix1++)  {
    rn=&psi[ix1];

    /* first spin component of result */
    hs0 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,0);
    hs1 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,1);
    hs2 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,2);
    hs3 = chi + halfspinor_buffer_offset(RECONS_GATHER,ix1,3);
    recons_4dir_minus(*hs0, *hs1, *hs2, *hs3, *rn);
    /* end of loop */
  }
}

static QMP_mem_t* xchi1;
static QMP_mem_t* xchi2;
static halfspinor_array* chi1;
static halfspinor_array* chi2;


/* Nearest neighbor communication channels */
static int total_comm = 0;
static QMP_msgmem_t forw_msg[4][2];
static QMP_msgmem_t back_msg[4][2];
static QMP_msghandle_t forw_mh[4][2];
static QMP_msghandle_t back_mh[4][2];
static QMP_msghandle_t forw_all_mh;
static QMP_msghandle_t back_all_mh;



void init_sse_su3dslash(const int latt_size[])   // latt_size not used, here for scalar version
{
  const int *machine_size = QMP_get_logical_dimensions();
  int bound[2][4][4];

  int mu, num, nsize;

  /* If we are already initialised, then increase the refcount and return */
  if (initP > 0) 
  {
    initP++;
    return;
  }


  /* Otherwise initialise */
  if (QMP_get_logical_number_of_dimensions() != 4)
  {
    QMP_error("init_sse_su3dslash: number of logical dimensions does not match problem");
    QMP_abort(1);
  }
    

  /* Check problem size */
  for(mu=0; mu < 4; mu++) 
    if ( latt_size[mu] == 1 ) 
    {
      QMP_error("This SSE Dslash does not support a problem size = 1. Here the lattice in dimension %d has length %d\n", mu, latt_size[mu]);
      QMP_abort(1);
    }


  num = latt_size[0] / machine_size[0];
  if ( num % 2 != 0 )
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
  nsize = 2*3*2*sizeof(double)*getSubgridVolCB()*3*4;  /* Note 3x4 half-ferm temps */
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
    
  for(mu=0; mu < 4; ++mu) 
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


void sse_su3dslash_wilson(SSEREAL *u, SSEREAL *psi, SSEREAL *res, int isign, int cb)
{

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


    if (total_comm > 0) {
      if (QMP_start(forw_all_mh) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	QMP_abort(1);
      }
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
    
      
    if (total_comm > 0) {
      if (QMP_start(forw_all_mh) != QMP_SUCCESS)
      {
	QMP_error("sse_su3dslash_wilson: QMP_start failed in forward direction");
	QMP_abort(1);
      }
    }
	
    /*other checkerboard's worth */
    dispatch_to_threads(decomp_hvv_minus,
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
