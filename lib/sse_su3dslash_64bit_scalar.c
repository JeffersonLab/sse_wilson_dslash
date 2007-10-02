/*******************************************************************************
 * $Id: sse_su3dslash_64bit_scalar.c,v 1.5 2007-10-02 20:40:21 bjoo Exp $
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


#include <sse_align.h>
#include <shift_tables_scalar.h>
#include <dispatch_scalar.h>
#include <site_dslash_64bit_scalar.h>


  
  static int initP = 0;          /* Reference count */
  static int *shift_table;       /* Shift Table */

  /** ---------------- This may need to change somehow for noncontiguous subsets ----------- **/
  static int icolor_start[2];    /* starting site for each coloring (cb) */
  static int icolor_end[2];      /* end site for each coloring (cb) */



static int init=0;
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
  vol_cb = getSubgridVolCB();

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



/* External routine */
void sse_su3dslash_wilson(double *u, double *psi, double *res, int isign, int cb)
{

  if (isign == 1)  {
    dispatch_to_threads(D_psi_fun_plus, 
		(spinor_array*)psi,
		(spinor_array*)res,
		(my_mat_array)u,
		1-cb,
		getSubgridVolCB());
  }

  if( isign == -1) {
    dispatch_to_threads(D_psi_fun_minus, 
		(spinor_array*)psi,
		(spinor_array*)res,
		(my_mat_array)u,
		1-cb,
		getSubgridVolCB());
  }
}

 


void D_psi_fun_plus(size_t lo, size_t hi, int id, const void *ptr)
{
  int ix,iy,iz;                                      /* Ix is corrent site */
                                                     /* Iy is a neighbour */
                                                     /* Iz is next site in loop */

  const ThreadWorkerArgs *a = (const ThreadWorkerArgs*)ptr;                /* Downcast argument */
  int cb = a->cb;
  const int low = icolor_start[cb] + lo;                                /* Start site */
  const int high = icolor_start[cb] + hi;                               /* End site + 1 */

  u_mat_array (*gauge_field)[4] ALIGN = a->u;        /* My gauge field */
  spinor_array *psi = a->psi;                        /* Source */
  spinor_array *res = a->res;                        /* Result */

  u_mat_array *up,*um;                               /* Pointer to FORWARD neighbour */
  spinor_array *s,*sp,*sm,*rn;                       /* Pointer to BACKWARD neighbour */

  spinor_array temp;

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
    rn=&res[ix];

    /******************************   Direction 0  ********************* */
    /* Prefetch back spinor for next dir: -0 */
    iy=backward_neighbor(shift_table,ix,0);
    sm=&psi[iy];
    prefetch_spinor(sm);  

    /* Prefetch back gauge field */
    um=&(gauge_field[iy][0]);
    prefetch_su3(um);

    dslash_plus_dir0_forward(*sp,*up,*rn);
      

    /* sm and um should already be prefetched */

    /* Now prefetch forward neighbor for next dir (1+) */
    /* And gauge field */
    sp=&psi[ forward_neighbor(shift_table,ix,1) ];
    prefetch_spinor(sp);
    up =&(gauge_field[ix][1]);
    prefetch_su3(up);

    dslash_plus_dir0_backward_add(*sm,*um,*rn);


    /********************** Direction 1 ************************ */
    /* up and sp should be prefetched */

    /* Prefetch backwards spinor and gauge field */
    iy=backward_neighbor(shift_table,ix,1);
    sm=&psi[iy];
    prefetch_spinor(sm);
    um=&(gauge_field[iy][1]);
    prefetch_su3(um);


    dslash_plus_dir1_forward_add(*sp,*up,*rn);



    /* Prefetch forwards spinor and gauge field for next dir: 2+ */
    iy=forward_neighbor(shift_table,ix,2);
    sp=&psi[iy];
    prefetch_spinor(sp);
    up = &(gauge_field[ix][2]);
    prefetch_su3(up);     


    dslash_plus_dir1_backward_add(*sm,*um,*rn);



    /********************** Direction 2 ************************* */
    /* Prefetch back spinor and gauge field  */

    iy=backward_neighbor(shift_table,ix,2);
    sm=&psi[iy];
    prefetch_spinor(sm);
    um=&(gauge_field[iy][2]);
    prefetch_su3(um);


    dslash_plus_dir2_forward_add(*sp,*up,*rn);


    /* Prefetch forward spinor and gauge field for next direction: 3+ */
    iy=forward_neighbor(shift_table,ix,3);
    sp=&psi[iy];
    prefetch_spinor(sp);
    up = &(gauge_field[ix][3]);
    prefetch_su3(up);

    dslash_plus_dir2_backward_add(*sm,*um,*rn);

      
    /* ******************* Direction +3 **************************** */
    /* Prefetch back spinor and gauge field  3- */

    iy=backward_neighbor(shift_table,ix,3);
    sm=&psi[iy];
    prefetch_spinor(sm);
    um=&(gauge_field[iy][3]);
    prefetch_su3(um);

    dslash_plus_dir3_forward_add(*sp,*up,*rn);      

    /* Next site */
    iz=ix+1;

    if (iz == high) { /* If we're on the last site, prefetch first site to avoid */
      iz=0;           /* Running beyond array bounds */
    }


    /* Prefetch forward spinor and gauge field for next site, dir 0 */
    iy=forward_neighbor(shift_table,iz,0);
    sp=&psi[iy];
    prefetch_spinor(sp);
    up=&(gauge_field[iz][0]);
    prefetch_su3(up);


    dslash_plus_dir3_backward_add_store(*sm,*um,*rn);      

  }
}



void D_psi_fun_minus(size_t lo, size_t hi, int id, const void *ptr )
{
  int ix,iy,iz;                          /* ix is the current site */
                                         /* iy is the neighbour for prefetching */
                                         /* iz is the prefetch site for the 
					    next loop iteration */

  const ThreadWorkerArgs *a = (const ThreadWorkerArgs*)ptr;    /* Cast the void args pointer */
  const int cb = a->cb;

  const int low = icolor_start[cb] + lo;                     /* First site */
  const int high = icolor_start[cb] + hi;                    /* Last site+1 */

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
    rn=&res[ix]; /* Pouinter to result */


    /* Prefetch back spinor and gauge field */
    iy=backward_neighbor(shift_table,ix,0);
    sm=&psi[iy];
    prefetch_spinor(sm);  
    um=&(gauge_field[iy][0]);
    prefetch_su3(um);

    dslash_minus_dir0_forward(*sp,*up,*rn);

    
    /* Prefetch spinor and gauge field for next direction (1) */
    iy=forward_neighbor(shift_table,ix,1);
    sp=&psi[iy];
    prefetch_spinor(sp);
    up = &(gauge_field[ix][1]);
    prefetch_su3(up);

    dslash_minus_dir0_backward_add(*sm,*um,*rn);    

    /* ******************* direction +1 ******************* */
    /* Prefetch backward spinor and gauge field */

    iy=backward_neighbor(shift_table,ix,1);
    sm=&psi[iy];
    prefetch_spinor(sm);
    um=&(gauge_field[iy][1]);
    prefetch_su3(um);

    dslash_minus_dir1_forward_add(*sp,*up,*rn);

          
    /* Prefetch forward spinor and gauge field for next direction: 2 */
    iy=forward_neighbor(shift_table,ix,2);
    sp=&psi[iy];
    prefetch_spinor(sp);
    up =&(gauge_field[ix][2]);
    prefetch_su3(up);

    dslash_minus_dir1_backward_add(*sm,*um,*rn);    

    
    /* ******************* direction +2 **************************** */
    /* Prefetch backward spinor and gauge field */

    iy=backward_neighbor(shift_table,ix,2); 
    sm=&psi[iy]; 
    prefetch_spinor(sm); 
    um=&(gauge_field[iy][2]);
    prefetch_su3(um);

    dslash_minus_dir2_forward_add(*sp,*up,*rn);    
        
    
    /* Prefetch spinor and gauge field for next direction: 3+ */
    iy=forward_neighbor(shift_table,ix,3);
    sp=&psi[iy];
    prefetch_spinor(sp);

    /* Prefetch gauge field for nex dir: 3+ */
    up = &(gauge_field[ix][3]);
    prefetch_su3(up);

    dslash_minus_dir2_backward_add(*sm,*um,*rn);    

        
      
    /*******************  direction +3 ************************** */
    /* Prefetch backward spinor and gauge field */
    iy=backward_neighbor(shift_table,ix,3);
    sm=&psi[iy];
    prefetch_spinor(sm);
    um=&(gauge_field[iy][3]);
    prefetch_su3(um);

    dslash_minus_dir3_forward_add(*sp,*up,*rn);    


    /* Next site in loop. We peeled this off the loop to start with so we can prefetch... */
    iz=ix+1; 
    if (iz==high) { /* If we are on the last site, we should prefetch the first element, to
		       avoid running past the array bounds */
      iz=0;
    }

    /* Prefetch the spinor and gauge field for next site, dir 0+ */
    iy=forward_neighbor(shift_table,iz,0);
    sp=&psi[iy];
    prefetch_spinor(sp);
    /* Prefetch the gauge field for next site dir  0+ */
    up=&(gauge_field[iz][0]);
    prefetch_su3(up);

    dslash_minus_dir3_backward_add_store(*sm,*um,*rn);    

  }
}


#ifdef __cplusplus
}
#endif



