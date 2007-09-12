/*******************************************************************************
 * $Id: sse_su3dslash_32bit_scalar.c,v 1.2 2007-09-12 21:00:50 bjoo Exp $
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



/* note we will have to change the way we handle multidimensional arrays because they will be dynamic and not static */

/* note for using array indices instead of structures, will need new types, or macro overlays 
   especailly for sse_pair_load and fpr sse_pair_store */

/* iup idn shift table overlays with packed/unpacked gauge lookerupper*/

/* include stuff from the Word file */
#include <sse_align.h>

  /* Note: packed data objects must be aligned on 16-byte or more boundaries or SSE
     inline asm will segfault....I've found that I needed to redeclare constants such as
     {1,-1,1,-1} in each routine they are used in */
  
  /* note: I align on larger multiples of 16 bytes for cache line alignment */


  static int total_vol = 1;
  static int initP = 0;

  extern void make_shift_tables(int *soffsets, int icolor_start[2], const int lat_size[4]);
  
  static int *soffsets;
  static int icolor_start[2];    /* starting site for each coloring (cb) */
  static int icolor_end[2];      /* end site for each coloring (cb) */
  
#define iup(mysite,mymu) soffsets[mymu + Nd*(mysite + total_vol*(1))]
#define idn(mysite,mymu) soffsets[mymu + Nd*(mysite + total_vol*(0))]

  /*isign = 0 => *(soffsets+mysite + global_vol_cb*2*(1-cb)+global_vol_cb*2*2*dir
    isign = 1 =>  *(soffsets+mysite + global_vol_cb*2*(1-cb)+global_vol_cb*2*2*dir+global_vol_cb) */


  /* SZIN uses the NO_U_PACK clause (referring to not interleaving adjacent sites and directions' worth of gauge fields), even though the gauge fields must be packed according to pack_gauge_field,
     however the default was packed...that's why the statements look non-intuitive...packed will work if you use Luescher's lattice geometry and insert
     special switches to handle the meshing in direction 3 that he has in his code  */


#define _gauge_field0_0(mysite) gauge_field[mysite][0]
#define _gauge_field0_1(mysite) gauge_field[mysite+1][0]
#define _gauge_field0_2(mysite) gauge_field[mysite][1]
#define _gauge_field0_3(mysite) gauge_field[mysite+1][1]
#define _gauge_field1_0(mysite) gauge_field[mysite][2]
#define _gauge_field1_1(mysite) gauge_field[mysite+1][2]
#define _gauge_field1_2(mysite) gauge_field[mysite][3]
#define _gauge_field1_3(mysite) gauge_field[mysite+1][3]


/* now overlays for spinors as arrays or structs */
  typedef float chi_float[2];
  typedef chi_float chi_two[2];
  typedef float u_mat_array[3][3][2]  ALIGN;  /* color color re/im */ 
  typedef float spinor_array[4][3][2] ALIGN; /* Nspin4 color re/im */
  typedef chi_two chi_array[3]    ALIGN; /*..color Nspin2 re/im ::note:: color slowest varying */
  typedef u_mat_array (*my_mat_array)[4] ALIGN;  


#define MY_SPINOR spinor_array
#define MY_SSE_VECTOR chi_array
#define MY_SSE_FLOAT sse_float  
#define _c1__ [0] 
#define _c2__ [1] 
#define _c3__ [2]
#define _c4__ [3]

#define MY_GAUGE u_mat_array
#define MY_GAUGE_ARRAY my_mat_array


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
  MY_SPINOR *psi;    /* input spinor */
  MY_SPINOR *res;    /* output spinor */
  MY_GAUGE (*u)[4];    /* gauge field on the output checkerboard */
  MY_GAUGE (*u2)[4];   /* gauge field on the input checkerboard */
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

  
#include <sse32.h>


#define BASE 0x3f

static int init=1;
static MY_SSE_FLOAT fact1,fact2;
/*static MY_SSE_VECTOR r12_1,r34_1,r12_2,r34_2;*/
/* thse must be actually declared within the routine to be assured of
being aligned properly */




/*void D_psi_fun(size_t lo,size_t hi, int id, const void *ptr );*/




void D_psi_fun_plus(size_t lo,size_t hi, int id, const void *ptr)
{
/* declare stuff */
  const Arg_s *a  = (const Arg_s*)ptr; 
  int ix1,iy1,iy2,iz1;
  const int low  =  lo; 
  const int high  =  hi;
  MY_GAUGE (*gauge_field)[4]  =  a->u;
  MY_SPINOR *psi  =  a->psi;
  MY_SPINOR *res  =  a->res;
  /* there is nothing wrong with superfluous ALIGN's */
  MY_GAUGE *up1 ALIGN, *up2 ALIGN, *um1 ALIGN,*um2 ALIGN;   
  MY_SPINOR *s1 ALIGN,*sp1 ALIGN,*sp2 ALIGN,*sm1 ALIGN,*sm2 ALIGN,*sn1 ALIGN;
  MY_SSE_VECTOR r12_1,r34_1,r12_2,r34_2;
  /* note these guys need to be declared in each routine in which they are used*/
  static sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};
  static sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};
  static sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};
  static sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};
  static sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};
  static sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};
  static sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f};
  
//  printf("low=%i high=%i\n", lo, hi);
  
	
  /*
    if (init == 0)
    check_alignment();
  */

  /* note that we want the spinors from the checkerboard opposite the one we are writing to */
  iy1 = iup(low,0);
  iy2 = iup(low+1,0);

//  printf("iup(ix1=%d, dir=0) ->  iy1=%d iy2=%d\n", low,iy1,iy2);

  sp1 = &psi[iy1];
  sp2 = &psi[iy2];
  /* note: the gauge fields at x - mu^ are on the checkerboard opposite the one we are writing to */
  up1 = &_gauge_field0_0(low);
  up2 = &_gauge_field0_1(low);
  /************************ loop over all lattice sites *************************/

  for (ix1 = low;ix1<high;ix1+= 2) 
  {
    /*  printf("m0:%f  k:%i l:%i low:%i high:%i id:%i\n, ix1:%i\n", m0, k, l, low, high, id, ix1);*/
    /*s1 = &psi[ix1];
      _prefetch_spinor(s1);*/
      
    /******************************* direction +0 *********************************/

    /* ...(1-isign*gamma(0))... */

    /* the signs on the isign*gamma(mu^) are different than the direction in which the
       input spinors are pulled from */
      
    /* compute shift addresses by reading from table...note: it may be faster to compute
       these on the fly using MMX registers which can run concurrent with XMM */
    iy1 = idn(ix1,0);
    iy2 = idn(ix1+1,0);

//    printf("idn(ix1=%d, dir=0) ->  iy1=%d iy2=%d\n", ix1,iy1,iy2);

	  /* prefetch */	   
    sm1 = &psi[iy1];
    _prefetch_single(sm1);
    sm2 = &psi[iy2];
    _prefetch_single(sm2);

    /*loads in the appropriate spinor and does the spin projection */
    /* sp denotes spinors in the x+mu^ direction, sm denotes spinors in the x - mu^ direction */
         
    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma0_minus();
    
    /* do the SU(3) * 3x2 multiply ....3 colors...both spin components of the spin-projected halfspinor*/ 
    
    _sse_su3_multiply((*up1));
    

       
    /* the output is in xmm3-5, so we store_up...this is a _set since it is the first term in the sum over mu */

    /* top two spin components, the top two rows of the 4x2 reconstruction matrix will be the identity, so just store */
    _sse_vector_store_up(r12_1);

      
        
    /* bottom two spin components, use the bottom two components of the 4x2 reconstruction matrix ...whatever was done to the identity
       to get that submatrix needs to be done likewise to the halfspinor to reconstruct the bottom two components */     
    _sse_24_gamma0_minus_set();
    _sse_vector_store_up(r34_1);

      
    /* prefetch */

    um1 = &_gauge_field0_0(iy1);
    _prefetch_single(um1);
    um2 = &_gauge_field0_0(iy2);
    _prefetch_single(um2);

		/*same stuff again...this time, for the next site in lexical order*/

    _sse_pair_load((*(sp2))_c1__,(*(sp2))_c2__);
    _sse_pair_load_up((*(sp2))_c3__,(*(sp2))_c4__);

    _sse_42_gamma0_minus();

    _sse_su3_multiply((*(up2)));

        
    
    _sse_vector_store_up(r12_2);

     
     
    _sse_24_gamma0_minus_set(); 
    _sse_vector_store_up(r34_2);
      
    /******************************* direction -0 *********************************/


      /* ...(1+isign*gamma(0))... */

      /* same kind of stuff again */

    iy1 = iup(ix1,1);
    iy2 = iup(ix1+1,1);
  
//    printf("iup(ix1=%d, dir=1) ->  iy1=%d iy2=%d\n", ix1,iy1,iy2);

    sp1 = &psi[iy1];
    _prefetch_single(sp1);
    sp2 = &psi[iy2];
    _prefetch_single(sp2);


    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _sse_42_gamma0_plus();
	 

    _sse_su3_inverse_multiply((*um1));


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


      /* now for the other site's worth */

    up1 = &_gauge_field0_2(ix1);
    _prefetch_single(up1);      
    up2 = &_gauge_field0_3(ix1);
    _prefetch_single(up2);
    

      
    _sse_pair_load((*(sm2))_c1__,(*(sm2))_c2__);
    _sse_pair_load_up((*(sm2))_c3__,(*(sm2))_c4__);
    _sse_42_gamma0_plus();

    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);

    _sse_vector_load(r34_2);
    _sse_24_gamma0_plus_add();
    _sse_vector_store(r34_2);  
      
      /******************************* direction +1 *********************************/


      /*same kind of stuff.......different stuff comes in at direction 3 */

    iy1 = idn(ix1,1);
    iy2 = idn(ix1+1,1);
    sm1 = &psi[iy1];
    _prefetch_single(sm1);
    sm2 = &psi[iy2];
    _prefetch_single(sm2);

    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma1_minus();

    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);

  
    _sse_vector_add();
	  
    _sse_vector_store(r12_1);


    _sse_vector_load(r34_1);
    _sse_24_gamma1_minus();
    _sse_vector_store(r34_1);

    um1 = &_gauge_field0_2(iy1);
    _prefetch_single(um1);      
    um2 = &_gauge_field0_2(iy2);
    _prefetch_single(um2);

      
    _sse_pair_load((*(sp2))_c1__,(*(sp2))_c2__);
    _sse_pair_load_up((*(sp2))_c3__,(*(sp2))_c4__);
    _sse_42_gamma1_minus();

    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2); 
      
    _sse_vector_load(r34_2);
    _sse_24_gamma1_minus();
    _sse_vector_store(r34_2);

      /******************************* direction -1 *********************************/





    iy1 = iup(ix1,2);
    iy2 = iup(ix1+1,2);
    sp1 = &psi[iy1];
    _prefetch_single(sp1);
    sp2 = &psi[iy2];
    _prefetch_single(sp2);

    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _sse_42_gamma1_plus();
	  

    _sse_su3_inverse_multiply((*um1));

    _sse_vector_load(r12_1);

    _sse_vector_add();
    _sse_vector_store(r12_1);      


    _sse_vector_load(r34_1);
    _sse_24_gamma1_plus();
    _sse_vector_store(r34_1);      

    up1 = &_gauge_field1_0(ix1);
    _prefetch_single(up1);      
    up2 = &_gauge_field1_1(ix1);
    _prefetch_single(up2);

    _sse_pair_load((*(sm2))_c1__,(*(sm2))_c2__);
    _sse_pair_load_up((*(sm2))_c3__,(*(sm2))_c4__);
    _sse_42_gamma1_plus();

    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);      

    _sse_vector_load(r34_2);
    _sse_24_gamma1_plus();
    _sse_vector_store(r34_2); 

      /******************************* direction +2 *********************************/





    iy1 = idn(ix1,2);
    sm1 = &psi[iy1];
    _prefetch_single(sm1);
    iy2 = idn(ix1+1,2);
    sm2 = &psi[iy2];
    _prefetch_single(sm2);

    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma2_minus();

    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);       

    _sse_vector_load(r34_1);
    _sse_24_gamma2_minus();
    _sse_vector_store(r34_1);       

    um1 = &_gauge_field1_0(iy1);
    _prefetch_single(um1);
    um2 = &_gauge_field1_0(iy2);
    _prefetch_single(um2);

    _sse_pair_load((*(sp2))_c1__,(*(sp2))_c2__);
    _sse_pair_load_up((*(sp2))_c3__,(*(sp2))_c4__);
    _sse_42_gamma2_minus();

    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);       

    _sse_vector_load(r34_2);
    _sse_24_gamma2_minus();
    _sse_vector_store(r34_2); 

      /******************************* direction -2 *********************************/





    iy1 = iup(ix1,3);
    sp1 = &psi[iy1];
    _prefetch_single(sp1);
    iy2 = iup(ix1+1,3); 
    sp2 = &psi[iy2];
    _prefetch_single(sp2);
    
      
      
    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _sse_42_gamma2_plus();      

    _sse_su3_inverse_multiply((*um1));

    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);
      
    _sse_vector_load(r34_1);
    _sse_24_gamma2_plus();
    
    _sse_vector_store(r34_1);

    up1 = &_gauge_field1_2(ix1);
    _prefetch_single(up1);
    up2 = &_gauge_field1_3(ix1);
    _prefetch_single(up2);

    _sse_pair_load((*(sm2))_c1__,(*(sm2))_c2__);
    _sse_pair_load_up((*(sm2))_c3__,(*(sm2))_c4__);
    _sse_42_gamma2_plus();      

    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);
      
    _sse_vector_load(r34_2);
    _sse_24_gamma2_plus();
    _sse_vector_store(r34_2);

      /******************************* direction +3 *********************************/





    iy1 = idn(ix1,3);
    iy2 = idn(ix1+1,3);
    sm1 = &psi[iy1];
    _prefetch_single(iy1);
    sm2 = &psi[iy2];
    _prefetch_single(iy2);

    


    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma3_minus();
      
    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);
    _sse_24_gamma3_minus_rows12();
    _sse_vector_store(r12_1);

    _sse_vector_load(r34_1);
    _sse_24_gamma3_minus();
    _sse_vector_store(r34_1);      

     

    um1 = &_gauge_field1_2(iy1); 
    _prefetch_single(um1);
    um2 = &_gauge_field1_2(iy2);
    _prefetch_single(um2);

    _sse_pair_load((*sp2)_c1__,(*sp2)_c2__);
    _sse_pair_load_up((*sp2)_c3__,(*sp2)_c4__);
    _sse_42_gamma3_minus();
      
    _sse_su3_multiply((*(up2)));

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

    iz1 = ix1+2;
    if (iz1 == high)
      iz1 = 0;
    iy1 = iup(iz1,0);
    sp1 = &psi[iy1];
    _prefetch_single(sp1);
    iy2 = iup(iz1+1,0);
    sp2 = &psi[iy2];
    _prefetch_single(sp2);
    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _sse_42_gamma3_plus();
      
    _sse_su3_inverse_multiply((*um1));

    sn1 = &res[ix1];  /*we always walk across the result lexicographically */
     
      
    _sse_vector_load(r12_1);
    _sse_24_gamma3_plus_rows12();
      
    _sse_pair_store((*sn1)_c1__,(*sn1)_c2__);

    _sse_vector_load(r34_1);
    _sse_24_gamma3_plus();
      
    _sse_pair_store((*sn1)_c3__,(*sn1)_c4__);      

    up1 = &_gauge_field0_0(iz1);
    _prefetch_single(up1);
    up2 = &_gauge_field0_1(iz1);
    _prefetch_single(up2);

    _sse_pair_load((*sm2)_c1__,(*sm2)_c2__);
    _sse_pair_load_up((*sm2)_c3__,(*sm2)_c4__);
    _sse_42_gamma3_plus();
      
    _sse_su3_inverse_multiply((*um2));

    _sse_vector_load(r12_2);
    _sse_24_gamma3_plus_rows12();
      
    _sse_pair_store((*(sn1+1))_c1__,(*(sn1+1))_c2__);

    _sse_vector_load(r34_2);
    _sse_42_gamma3_plus();
     
    _sse_pair_store((*(sn1+1))_c3__,(*(sn1+1))_c4__); 
	  

    /******************************** end of loop *********************************/
      
  }


}

/*ok, now this routine here is just like the one above except isign has a different value, which means that the
signs used for 1 +- gamma(mu^) must be swapped...*/ 


void D_psi_fun_minus(size_t lo,size_t hi, int id, const void *ptr)
{
  const Arg_s *a  = (const Arg_s*)ptr; 
  int ix1,iy1,iy2,iz1;
  const int low  =  lo; 
  const int high  =  hi;
  MY_GAUGE (*gauge_field)[4]  =  a->u;
  MY_SPINOR *psi  =  a->psi;
  MY_SPINOR *res  =  a->res;
  MY_GAUGE *up1 ALIGN, *up2 ALIGN, *um1 ALIGN,*um2 ALIGN;   
  MY_SPINOR *s1 ALIGN,*sp1 ALIGN,*sp2 ALIGN,*sm1 ALIGN,*sm2 ALIGN,*sn1 ALIGN;
  MY_SSE_VECTOR r12_1,r34_1,r12_2,r34_2;
  
  static sse_float _sse_sgn12 ALIGN ={-1.0f,-1.0f,1.0f,1.0f};
  static sse_float _sse_sgn13 ALIGN ={-1.0f,1.0f,-1.0f,1.0f};
  static sse_float _sse_sgn14 ALIGN ={-1.0f,1.0f,1.0f,-1.0f};
  static sse_float _sse_sgn23 ALIGN ={1.0f,-1.0f,-1.0f,1.0f};
  static sse_float _sse_sgn24 ALIGN ={1.0f,-1.0f,1.0f,-1.0f};
  static sse_float _sse_sgn34 ALIGN ={1.0f,1.0f,-1.0f,-1.0f};
  static sse_float _sse_sgn1234 ALIGN = {-1.0f,-1.0f,-1.0f,-1.0f};

  /*printf("%i %i",low, high);*/


  /*
    if (init == 0)
    check_alignment();
  */


  iy1 = iup(low,0);
  iy2 = iup(low+1,0);
  sp1 = &psi[iy1];
  sp2 = &psi[iy2];
  up1 = &_gauge_field0_0(low);
  up2 = &_gauge_field0_1(low);
  /************************ loop over all lattice sites *************************/

  for (ix1 = low;ix1<high;ix1+= 2) 
  {
    /*  printf("m0:%f  k:%i l:%i low:%i high:%i id:%i\n, ix1:%i\n", m0, k, l, low, high, id, ix1);*/
    /*s1 = &psi[ix1];
	_prefetch_spinor(s1);*/
      
      /******************************* direction +0 *********************************/

      /* ...(1-isign*gamma(0))... */



    iy1 = idn(ix1,0);
    iy2 = idn(ix1+1,0);
   
    sm1 = &psi[iy1];
    _prefetch_single(sm1);
    sm2 = &psi[iy2];
    _prefetch_single(sm2);



    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma0_plus();
    

    
    _sse_su3_multiply((*up1));


       
    
     
    _sse_vector_store_up(r12_1);



     
    _sse_24_gamma0_plus_set();
    _sse_vector_store_up(r34_1);





    um1 = &_gauge_field0_0(iy1);
    _prefetch_single(um1);
    um2 = &_gauge_field0_0(iy2);
    _prefetch_single(um2);

    _sse_pair_load((*(sp2))_c1__,(*(sp2))_c2__);
    _sse_pair_load_up((*(sp2))_c3__,(*(sp2))_c4__);

    _sse_42_gamma0_plus();

    _sse_su3_multiply((*(up2)));

        
    
    _sse_vector_store_up(r12_2);

     
     
    _sse_24_gamma0_plus_set(); 
    _sse_vector_store_up(r34_2);
      
      /******************************* direction -0 *********************************/


      /* ...(1+isign*gamma(0))... */



    iy1 = iup(ix1,1);
    iy2 = iup(ix1+1,1);
  
    sp1 = &psi[iy1];
    _prefetch_single(sp1);
    sp2 = &psi[iy2];
    _prefetch_single(sp2);


    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _sse_42_gamma0_minus();
	 

    _sse_su3_inverse_multiply((*um1));




    _sse_vector_load(r12_1);

    _sse_vector_add();
    _sse_vector_store(r12_1);
    



    _sse_vector_load(r34_1);
    _sse_24_gamma0_minus_add();
    _sse_vector_store(r34_1);  




    up1 = &_gauge_field0_2(ix1);
    _prefetch_single(up1);      
    up2 = &_gauge_field0_3(ix1);
    _prefetch_single(up2);



    _sse_pair_load((*(sm2))_c1__,(*(sm2))_c2__);
    _sse_pair_load_up((*(sm2))_c3__,(*(sm2))_c4__);
    _sse_42_gamma0_minus();

    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);

    _sse_vector_load(r34_2);
    _sse_24_gamma0_minus_add();
    _sse_vector_store(r34_2);  
      
      /******************************* direction +1 *********************************/




    iy1 = idn(ix1,1);
    iy2 = idn(ix1+1,1);
    sm1 = &psi[iy1];
    _prefetch_single(sm1);
    sm2 = &psi[iy2];
    _prefetch_single(sm2);

    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma1_plus();

    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);

  
    _sse_vector_add();
	  
    _sse_vector_store(r12_1);


    _sse_vector_load(r34_1);
    _sse_24_gamma1_plus();
    _sse_vector_store(r34_1);

    um1 = &_gauge_field0_2(iy1);
    _prefetch_single(um1);      
    um2 = &_gauge_field0_2(iy2);
    _prefetch_single(um2);

      
    _sse_pair_load((*(sp2))_c1__,(*(sp2))_c2__);
    _sse_pair_load_up((*(sp2))_c3__,(*(sp2))_c4__);
    _sse_42_gamma1_plus();

    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2); 
      
    _sse_vector_load(r34_2);
    _sse_24_gamma1_plus();
    _sse_vector_store(r34_2);

      /******************************* direction -1 *********************************/





    iy1 = iup(ix1,2);
    iy2 = iup(ix1+1,2);
    sp1 = &psi[iy1];
    _prefetch_single(sp1);
    sp2 = &psi[iy2];
    _prefetch_single(sp2);

    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _sse_42_gamma1_minus();
	  

    _sse_su3_inverse_multiply((*um1));

    _sse_vector_load(r12_1);

    _sse_vector_add();
    _sse_vector_store(r12_1);      


    _sse_vector_load(r34_1);
    _sse_24_gamma1_minus();
    _sse_vector_store(r34_1);      

    up1 = &_gauge_field1_0(ix1);
    _prefetch_single(up1);      
    up2 = &_gauge_field1_1(ix1);
    _prefetch_single(up2);

    _sse_pair_load((*(sm2))_c1__,(*(sm2))_c2__);
    _sse_pair_load_up((*(sm2))_c3__,(*(sm2))_c4__);
    _sse_42_gamma1_minus();

    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);      

    _sse_vector_load(r34_2);
    _sse_24_gamma1_minus();
    _sse_vector_store(r34_2); 

      /******************************* direction +2 *********************************/





    iy1 = idn(ix1,2);
    sm1 = &psi[iy1];
    _prefetch_single(sm1);
    iy2 = idn(ix1+1,2);
    sm2 = &psi[iy2];
    _prefetch_single(sm2);

    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma2_plus();

    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);       

    _sse_vector_load(r34_1);
    _sse_24_gamma2_plus();
    _sse_vector_store(r34_1);       

    um1 = &_gauge_field1_0(iy1);
    _prefetch_single(um1);
    um2 = &_gauge_field1_0(iy2);
    _prefetch_single(um2);

    _sse_pair_load((*(sp2))_c1__,(*(sp2))_c2__);
    _sse_pair_load_up((*(sp2))_c3__,(*(sp2))_c4__);
    _sse_42_gamma2_plus();

    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);       

    _sse_vector_load(r34_2);
    _sse_24_gamma2_plus();
    _sse_vector_store(r34_2); 

      /******************************* direction -2 *********************************/





    iy1 = iup(ix1,3);
    sp1 = &psi[iy1];
    _prefetch_single(sp1);
    iy2 = iup(ix1+1,3); 
    sp2 = &psi[iy2];
    _prefetch_single(sp2);
    
      /*if (iy1<iy2)
	{
	sp2 = sp1+1;
	_prefetch_spinor(sp1);
	}
	else
	{
	sp2 = sp1-1;
	_prefetch_spinor(sp2);
	} */
      
    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _sse_42_gamma2_minus();      

    _sse_su3_inverse_multiply((*um1));

    _sse_vector_load(r12_1);
    _sse_vector_add();
    _sse_vector_store(r12_1);
      
    _sse_vector_load(r34_1);
    _sse_24_gamma2_minus();
    
    _sse_vector_store(r34_1);

    up1 = &_gauge_field1_2(ix1);
    _prefetch_single(up1);
    up2 = &_gauge_field1_3(ix1);
    _prefetch_single(up2);

    _sse_pair_load((*(sm2))_c1__,(*(sm2))_c2__);
    _sse_pair_load_up((*(sm2))_c3__,(*(sm2))_c4__);
    _sse_42_gamma2_minus();      

    _sse_su3_inverse_multiply((*(um2)));

    _sse_vector_load(r12_2);
    _sse_vector_add();
    _sse_vector_store(r12_2);
      
    _sse_vector_load(r34_2);
    _sse_24_gamma2_minus();
    _sse_vector_store(r34_2);

      /******************************* direction +3 *********************************/





    iy1 = idn(ix1,3);
    iy2 = idn(ix1+1,3);
    sm1 = &psi[iy1];
    _prefetch_single(iy1);
    sm2 = &psi[iy2];
    _prefetch_single(iy2);

    /*if (iy1<iy2)
      {
      sm2 = sm1+1;
      _prefetch_spinor(sm1);
      }
      else
      {
      sm2 = sm1-1;
      _prefetch_spinor(sm2);
      } */



    _sse_pair_load((*sp1)_c1__,(*sp1)_c2__);
    _sse_pair_load_up((*sp1)_c3__,(*sp1)_c4__);
    _sse_42_gamma3_plus();
      
    _sse_su3_multiply((*up1));

    _sse_vector_load(r12_1);
    _sse_24_gamma3_plus_rows12();
    _sse_vector_store(r12_1);

    _sse_vector_load(r34_1);
    _sse_24_gamma3_plus();
    _sse_vector_store(r34_1);      

      /*if (iy1<iy2)
	{
	um1 = &_gauge_field1_0(iy2);
	um2 = um1+1;
	_prefetch_su3(um1);
	}
	else
	{
	um2 = &_gauge_field1_0(iy1);
	um1 = um2+1;
	_prefetch_su3(um2);
	} */

    um1 = &_gauge_field1_2(iy1); 
    _prefetch_single(um1);
    um2 = &_gauge_field1_2(iy2);
    _prefetch_single(um2);

    _sse_pair_load((*sp2)_c1__,(*sp2)_c2__);
    _sse_pair_load_up((*sp2)_c3__,(*sp2)_c4__);
    _sse_42_gamma3_plus();
      
    _sse_su3_multiply((*(up2)));

    _sse_vector_load(r12_2);
    _sse_24_gamma3_plus_rows12();
    _sse_vector_store(r12_2);

    _sse_vector_load(r34_2);
    _sse_24_gamma3_plus();
    _sse_vector_store(r34_2); 

      /******************************* direction -3 *********************************/





    iz1 = ix1+2;
    if (iz1 == high)
      iz1 = 0;
    iy1 = iup(iz1,0);
    sp1 = &psi[iy1];
    _prefetch_single(sp1);
    iy2 = iup(iz1+1,0);
    sp2 = &psi[iy2];
    _prefetch_single(sp2);
    _sse_pair_load((*sm1)_c1__,(*sm1)_c2__);
    _sse_pair_load_up((*sm1)_c3__,(*sm1)_c4__);
    _sse_42_gamma3_minus();
      
    _sse_su3_inverse_multiply((*um1));

    sn1 = &res[ix1];  /*we always walk across the result lexicographically */
     
      
    _sse_vector_load(r12_1);
    _sse_24_gamma3_minus_rows12();
      
    _sse_pair_store((*sn1)_c1__,(*sn1)_c2__);

    _sse_vector_load(r34_1);
    _sse_24_gamma3_minus();
      
    _sse_pair_store((*sn1)_c3__,(*sn1)_c4__);      

    up1 = &_gauge_field0_0(iz1);
    _prefetch_single(up1);
    up2 = &_gauge_field0_1(iz1);
    _prefetch_single(up2);

    _sse_pair_load((*sm2)_c1__,(*sm2)_c2__);
    _sse_pair_load_up((*sm2)_c3__,(*sm2)_c4__);
    _sse_42_gamma3_minus();
      
    _sse_su3_inverse_multiply((*um2));

    _sse_vector_load(r12_2);
    _sse_24_gamma3_minus_rows12();
      
    _sse_pair_store((*(sn1+1))_c1__,(*(sn1+1))_c2__);

    _sse_vector_load(r34_2);
    _sse_42_gamma3_minus();
     
    _sse_pair_store((*(sn1+1))_c3__,(*(sn1+1))_c4__); 
	  

    /******************************** end of loop *********************************/
      
  }


}









































#ifdef __cplusplus
}
#endif
