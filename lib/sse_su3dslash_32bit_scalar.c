/*******************************************************************************
 * $Id: sse_su3dslash_32bit_scalar.c,v 1.5 2007-10-02 20:40:21 bjoo Exp $
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
#include <types32.h>
#include <dispatch_scalar.h>
#include <site_dslash_32bit_scalar.h>
#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>



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

  vol_cb = getSubgridVolCB();
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


/* routine sse_su3dslash_wilson
   u: base pointer to gauge field
   psi: base pointer to input spinor field on full lattice
   res: base pointer to output spinor field on full lattice
   isign: 1-->normal, -1--> swaps 1 - gamma(mu^) for 1 + gamma(mu^)
   cb: checkerboard (0/1) of input fields
*/
void sse_su3dslash_wilson(float *u, float *psi, float *res, int isign, int cb)
{

  if (isign == 1) {
    dispatch_to_threads(D_psi_fun_plus, 
			(spinor_array*)psi,
			(spinor_array*)res,
			(my_mat_array)u,
			1-cb,
			getSubgridVolCB());
  }

  if( isign == -1) 
  {
    dispatch_to_threads(D_psi_fun_minus, 
			(spinor_array*)psi,
			(spinor_array*)res,
			(my_mat_array)u,
			1-cb,
			getSubgridVolCB());
  }
}

// #include <sse32.h>


void D_psi_fun_plus(size_t lo,size_t hi, int id, const void *ptr)
{

  const ThreadWorkerArgs *a  = (const ThreadWorkerArgs*)ptr;  /* Cast the (void *) to an (ThreadWorkerArgs*) */
  int ix1;                              /* Index of current site */
  int iy1,iy2;                          /* Index of neighbour of current site (iy1) 
					   and current site+1 (iy2) */

  int iz1;                              /* Index of next site in loop */

  u_mat_array (*gauge_field)[4]  =  a->u; /* Packed Gauge fields */
  spinor_array *psi  =  a->psi;           /* Source spinor */
  spinor_array *res  =  a->res;           /* Result spinor */

  const int cb = a->cb;

  const int low  =  icolor_start[cb]+lo;                 /* First site for this thread */
  const int high  = icolor_start[cb]+ hi;                /* Last site for this thread */

  

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
  halfspinor_array r12_1 ALIGN; /* Site 1 upper */
  halfspinor_array r34_1 ALIGN; /* Site 1 lower */
  halfspinor_array r12_2 ALIGN; /* Site 2 upper */
  halfspinor_array r34_2 ALIGN; /* Site 2 lower */


  
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
    prefetch_single(sm1);
    sm2 = &psi[ iy2   ];
    prefetch_single(sm2);

    /*loads in the appropriate spinor and does the spin projection */
    /* sp denotes spinors in the x+mu^ direction, sm denotes spinors in the x - mu^ direction */
    dslash_plus_dir0_forward(*sp1, *up1, r12_1, r34_1);

    /* Now prefetch for the  1 + \gamma_0 U^\dagger case */

    um1 = &(gauge_field[iy1][0]);
    prefetch_single(um1);
    um2 = &(gauge_field[iy2][0]);
    prefetch_single(um2);

    /* Now do the spin proj, multiply, spin recon for the second site */
    /* Load spinor */
    dslash_plus_dir0_forward(*sp2, *up2, r12_2, r34_2);

    /******************************* direction -0 *********************************/
    /* Hopefully sm1 sm2 the backward neighbours are prefetched and in cache. The shifted 
       matrices should be prefetched too */
   

      /* ...(1+isign*gamma(0))... */
    
    /* Prefetch the forward spinors for the next direction */
    sp1 = &psi[ forward_neighbor(shift_table,ix1,1) ];
    prefetch_single(sp1);
    sp2 = &psi[ forward_neighbor(shift_table,ix1+1,1)];
    prefetch_single(sp2);

    dslash_plus_dir0_backward_add(*sm1, *um1, r12_1, r34_1);

    /* Prefetch gauge field for next direction */
    up1 = &(gauge_field[ix1][1]);
    prefetch_single(up1);      
    up2 = &(gauge_field[ix1+1][1]); 
    prefetch_single(up2); 
    
    dslash_plus_dir0_backward_add(*sm2, *um2, r12_2, r34_2);
    /******************************* direction +1 *********************************/
    /* OK, sp1, sp2 and up1, up2 should be prefetched */

    /* Prefetch spinors for the -1 direction */
    iy1 = backward_neighbor(shift_table,ix1,1);
    iy2 = backward_neighbor(shift_table,ix1+1,1);
    sm1 = &psi[iy1];
    prefetch_single(sm1);
    sm2 = &psi[iy2];
    prefetch_single(sm2);

    dslash_plus_dir1_forward_add(*sp1, *up1, r12_1, r34_1);

    /* Prefetch gauge field for - direction */
    um1 = &(gauge_field[iy1][1]);
    prefetch_single(um1);      
    um2 = &(gauge_field[iy2][1]);
    prefetch_single(um2);

    dslash_plus_dir1_forward_add(*sp2, *up2, r12_2, r34_2);
    
    /******************************* direction -1 *********************************/

    /* Prefetch forward neighbour for direction 2+ */
    sp1 = &psi[forward_neighbor(shift_table,ix1,2)];
    prefetch_single(sp1);
    sp2 = &psi[forward_neighbor(shift_table,ix1+1,2)];
    prefetch_single(sp2);

    dslash_plus_dir1_backward_add(*sm1, *um1, r12_1, r34_1);

    /* Prefetch Gauge for next case: Direction 2 + */
    up1 = &(gauge_field[ix1][2]);
    prefetch_single(up1);      
    up2 = &(gauge_field[ix1+1][2]);
    prefetch_single(up2);

    dslash_plus_dir1_backward_add(*sm2, *um2, r12_2, r34_2);

    /******************************* direction +2 *********************************/
    /* sp1, sp2, up1, up2 should be prefetched and in cache */


    /* Prefetch sm1 & sm2 for -ve direction */
    iy1 = backward_neighbor(shift_table,ix1,2);
    sm1 = &psi[iy1];
    prefetch_single(sm1);
    iy2 = backward_neighbor(shift_table,ix1+1,2);
    sm2 = &psi[iy2];
    prefetch_single(sm2);

    dslash_plus_dir2_forward_add(*sp1, *up1, r12_1, r34_1);

    /* Prefetch gauge field for -ve case */
    um1 = &(gauge_field[iy1][2]);
    prefetch_single(um1);
    um2 = &(gauge_field[iy2][2]);
    prefetch_single(um2);

    dslash_plus_dir2_forward_add(*sp2, *up2, r12_2, r34_2);
    
    /******************************* direction -2 *********************************/
    /* sm1, sm2, um1, um2, should be prefetched and in cache                      */
    /******************************************************************************/
    
    /* Prefetch spinors for direction 3+ */
    sp1 = &psi[ forward_neighbor(shift_table,ix1,3) ];
    prefetch_single(sp1);
    sp2 = &psi[ forward_neighbor(shift_table,ix1+1,3) ];
    prefetch_single(sp2);

    dslash_plus_dir2_backward_add(*sm1, *um1, r12_1, r34_1);
    
      
    /* Prefetch Gauge for the 3+ case */
    up1 = &(gauge_field[ix1][3]);
    prefetch_single(up1);
    up2 = &(gauge_field[ix1+1][3]);
    prefetch_single(up2);

    dslash_plus_dir2_backward_add(*sm2, *um2, r12_2, r34_2);


    /******************************* direction +3 *********************************/
    /* sp1, sp2, up1, up2, should be prefetched and in cache already              */

    /* Pre Fetch neigbour spinors for the backwards case                          */
    iy1 = backward_neighbor(shift_table,ix1,3);
    iy2 = backward_neighbor(shift_table,ix1+1,3);
    sm1 = &psi[iy1]; 
    prefetch_single(iy1); 
    sm2 = &psi[iy2]; 
    prefetch_single(iy2) ;

    dslash_plus_dir3_forward_add(*sp1, *up1, r12_1, r34_1);

     
    /* Prefetch um for - case */
    um1 = &(gauge_field[iy1][3]); 
    prefetch_single(um1);
    um2 = &(gauge_field[iy2][3]);
    prefetch_single(um2);

    dslash_plus_dir3_forward_add(*sp2, *up2, r12_2, r34_2);

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
    prefetch_single(sp1);
    sp2 = &psi[ forward_neighbor(shift_table,iz1+1,0) ];
    prefetch_single(sp2);

    sn1 = &res[ix1];  /*we always walk across the result lexicographically */
    dslash_plus_dir3_backward_add_store(*sm1, *um1, r12_1, r34_1, *sn1);


    /* Prefetch gauge field for next loop iteration (0 direction) */
    up1 = &(gauge_field[iz1][0]);
    prefetch_single(up1);
    up2 = &(gauge_field[iz1+1][0]);
    prefetch_single(up2);

    dslash_plus_dir3_backward_add_store(*sm2, *um2, r12_2, r34_2, *(sn1+1));
	  
    /******************************** end of loop *********************************/
      
  }
}

/*ok, now this routine here is just like the one above except isign has a different value, which means that the
signs used for 1 +- gamma(mu^) must be swapped...*/ 
void D_psi_fun_minus(size_t lo,size_t hi, int id, const void *ptr)
{
  const ThreadWorkerArgs *a  = (const ThreadWorkerArgs*)ptr;   /* Downcast to args */
  int ix1,iy1,iy2,iz1;                   /* Coordinates ix1 - current
					    iy1 - index of neighbour
					    iy1 - index of neighbour of ix1+1 
					    iz1 - index of first of next pair (ix+2) */
  const int cb = a->cb;

  const int low  =  icolor_start[cb]+lo;                 /* First site for this thread */
  const int high  = icolor_start[cb]+ hi;                /* Last site for this thread */

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

  halfspinor_array r12_1;                         /* site 1 halfspinor top half */
  halfspinor_array r34_1;                         /* site 1 halfspinor bottom half */
  halfspinor_array r12_2;                         /* site 2 halfspinor top half */
  halfspinor_array r34_2;                         /* site 2 halfspinor bottom half */
  


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
    prefetch_single(sm1);
    sm2 = &psi[iy2];
    prefetch_single(sm2);

    dslash_minus_dir0_forward(*sp1, *up1, r12_1, r34_1);

    /* Prefetch gauge field for next part */
    um1 = &(gauge_field[iy1][0]);
    prefetch_single(um1);
    um2 = &(gauge_field[iy2][0]);
    prefetch_single(um2);

    dslash_minus_dir0_forward(*sp2, *up2, r12_2, r34_2);

      
    /******************************* direction -0 *********************************/
    /* sm1, sm2, um1, um2 should be prefetched */

    /* ...(1+isign*gamma(0))... */

    /* Prefetch for next part (1+) */
    iy1 = forward_neighbor(shift_table,ix1,1);
    iy2 = forward_neighbor(shift_table,ix1+1,1);
    sp1 = &psi[iy1];
    prefetch_single(sp1);
    sp2 = &psi[iy2];
    prefetch_single(sp2);

    dslash_minus_dir0_backward_add(*sm1, *um1, r12_1, r34_1);

    /* Prefetch gauge field for next part: (1+) */
    up1 = &(gauge_field[ix1][1]);
    prefetch_single(up1);      
    up2 = &(gauge_field[ix1+1][1]);
    prefetch_single(up2);

    dslash_minus_dir0_backward_add(*sm2, *um2, r12_2, r34_2);
      

    /******************************* direction +1 *********************************/
    /* sp1, sp2, up1, up2 should be prefetched                                    */

    /* Prefetch Spinors for next part (1-) */
    iy1 = backward_neighbor(shift_table,ix1,1);
    iy2 = backward_neighbor(shift_table,ix1+1,1);
    sm1 = &psi[iy1];
    prefetch_single(sm1);
    sm2 = &psi[iy2];
    prefetch_single(sm2);

    dslash_minus_dir1_forward_add(*sp1, *up1, r12_1, r34_1);

    /* Prefetch gauge links for next part (1-) */
    um1 = &(gauge_field[iy1][1]);
    prefetch_single(um1);      
    um2 = &(gauge_field[iy2][1]);
    prefetch_single(um2); 

    dslash_minus_dir1_forward_add(*sp2, *up2, r12_2, r34_2);
      


    /******************************* direction -1 *********************************/

    /* Prefetch spinors for next part: (2+) */
    iy1 = forward_neighbor(shift_table,ix1,2);
    iy2 = forward_neighbor(shift_table,ix1+1,2);
    sp1 = &psi[iy1];
    prefetch_single(sp1);
    sp2 = &psi[iy2];
    prefetch_single(sp2);

    dslash_minus_dir1_backward_add(*sm1, *um1, r12_1, r34_1);

    /* Prefetch gauge field for next part (2+) */
    up1 = &(gauge_field[ix1][2]);
    prefetch_single(up1);       
    up2 = &(gauge_field[ix1+1][2]); 
    prefetch_single(up2); 

    dslash_minus_dir1_backward_add(*sm2, *um2, r12_2, r34_2);

    /******************************* direction +2 *********************************/
    /* sp1, sp2, up1, up2 should  be in cache                                     */

    /* Prefetch spinors for next part: 2- */
    iy1 = backward_neighbor(shift_table,ix1,2);
    sm1 = &psi[iy1];
    prefetch_single(sm1);
     iy2 = backward_neighbor(shift_table,ix1+1,2);
    sm2 = &psi[iy2];
    prefetch_single(sm2);

    dslash_minus_dir2_forward_add(*sp1, *up1, r12_1, r34_1);

    /* Prefetch Gauge field for next part: 2- */
    um1 = &(gauge_field[iy1][2]);
    prefetch_single(um1);
    um2 = &(gauge_field[iy2][2]);
    prefetch_single(um2);

    dslash_minus_dir2_forward_add(*sp2, *up2, r12_2, r34_2);
 
    /******************************* direction -2 *********************************/

    /* sm1, sm2, um1, um2 should be cached */

    /* Prefetch spinor for next case: 3+ */
    iy1 = forward_neighbor(shift_table,ix1,3);
    sp1 = &psi[iy1];
    prefetch_single(sp1);
    iy2 = forward_neighbor(shift_table,ix1+1,3); 
    sp2 = &psi[iy2];
    prefetch_single(sp2);
    dslash_minus_dir2_backward_add(*sm1, *um1, r12_1, r34_1);

    /* Prefetch gauge for next case: 3+ */
    up1 = &(gauge_field[ix1][3]);
    prefetch_single(up1);
    up2 = &(gauge_field[ix1+1][3]);
    prefetch_single(up2);
     
    dslash_minus_dir2_backward_add(*sm2, *um2, r12_2, r34_2);
    
    /******************************* direction +3 *********************************/
    /* sp1, sp2, up1, up2 should be in cache */

    /* Prefetch spinors for next case: dir 3- */
    iy1 = backward_neighbor(shift_table,ix1,3);
    iy2 = backward_neighbor(shift_table,ix1+1,3);
    sm1 = &psi[iy1];
    prefetch_single(iy1);
    sm2 = &psi[iy2];
    prefetch_single(iy2);

    dslash_minus_dir3_forward_add(*sp1, *up1, r12_1, r34_1);

    /* Prefetch gauge field for next part: 3- */
    um1 = &(gauge_field[iy1][3]); 
    prefetch_single(um1);
    um2 = &(gauge_field[iy2][3]);
    prefetch_single(um2);

    dslash_minus_dir3_forward_add(*sp2, *up2, r12_2, r34_2);

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
    prefetch_single(sp1);
    iy2 = forward_neighbor(shift_table,iz1+1,0);
    sp2 = &psi[iy2];
    prefetch_single(sp2);

    sn1 = &res[ix1];     
    dslash_minus_dir3_backward_add_store(*sm1, *um1, r12_1, r34_1, *sn1);


    /* Prefetch Gauge field for next case: 0+ */
    up1 = &(gauge_field[iz1][0]);
    prefetch_single(up1);
    up2 = &(gauge_field[iz1+1][0]);
    prefetch_single(up2);

    dslash_minus_dir3_backward_add_store(*sm2, *um2, r12_2, r34_2, *(sn1+1));
    /******************************** end of loop *********************************/
  }

}

#ifdef __cplusplus
}
#endif
