/* $Id: shift_tables_scalar.c,v 1.4 2008-03-04 21:50:18 bjoo Exp $

/* Set the offset tables used by the 32-bit and 64-bit single node dslash */

/* make_shift_tables(latt_size) */
/* Get the offsets needed for neighbour comm. */
/* soffsets(position,direction,isign,cb)   */ 
/*  where  isign    = +1 : plus direction */
/*                  =  0 : negative direction */
/*         cb       =  0 : even lattice (includes origin) */
/*                  = +1 : odd lattice (does not include origin) */
/* the offsets cotain the current site, i.e the neighbour for site i  */
/* is  soffsets(i,dir,cb,mu) and NOT  i + soffset(..)    */
  
#include <shift_tables_scalar.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


  static int Nd=4;
  static int* xsite_table;
  int* site_table;

  typedef struct { 
    int cb4;
    int linearcb4;
  } InvTab4;

  /* Total problem size */
  static int tot_size[4];
  static int total_vol = -1;
  int total_vol_cb = -1;

  static void setLattSize(const int size[])
  {
    int i;
    for(i=0; i < 4; ++i) {
      tot_size[i] = size[i];
    }
    total_vol = size[0];
    for(i=1; i < 4; ++i) {
      total_vol *= size[i];
    }

    total_vol_cb = total_vol / 2;
  }


  static int* getLattSize()
  {
    return tot_size;
  }


#if 0

/* Decompose lexicographic site ipos into lattice coordinates */
static void crtesn(int coord[], int ipos, int latt_size[])
{
  int i;

  for(i=0; i < 4; ++i)
  {
    coord[i] = ipos % latt_size[i];
    ipos /= latt_size[i];
  }
}

/* Calculates the lexicographic site index from the coordinate of a site */
static int local_site(int coord[], int latt_size[])
{
  int order = 0;
  int mmu;

  for(mmu=4-1; mmu >= 1; --mmu)
    order = latt_size[mmu-1]*(coord[mmu] + order);

  order += coord[0];

  return order;
}


 
/*****************************************************************************************
/*********** These are functions taken from QDP++ and converted to C *****************/

//! Reconstruct the lattice coordinate from the node and site number
/*! 
 * This is the inverse of the nodeNumber and linearSiteIndex functions.
 * The API requires this function to be here.
 */
static void getSiteCoords(int coord[], int node, int linearsite) // ignore node
{
  int cb, cbb, m;
  int vol_cb = getSubgridVolCB();
  int cb_nrow[4];

  for(m=0; m < 4; ++m)
    cb_nrow[m] = getLattSize()[m];
  cb_nrow[0] /= 2;
  
  cb = linearsite / vol_cb;
  crtesn(coord, linearsite % vol_cb, cb_nrow);

  cbb = cb;
  for(m=1; m < 4; ++m)
    cbb += coord[m];
  cbb = (cbb % 2);

  coord[0] = 2*coord[0] + cbb;
}


//! The linearized site index for the corresponding coordinate
/*! This layout is appropriate for a 2 checkerboard (red/black) lattice */
int getLinearSiteIndex(const int coord[])
{
  int vol_cb = getSubgridVolCB();
  int cb_nrow[4];
  int cb_coord[4];
  int cb, m;

  for(m=0; m < 4; ++m)
    cb_nrow[m] = getLattSize()[m];
  cb_nrow[0] /= 2;
  
  for(m=0; m < 4; ++m)
    cb_coord[m] = coord[m];
  cb_coord[0] /= 2;    // Number of checkerboards
    
  cb = 0;
  for(m=0; m < 4; ++m)
    cb += coord[m];

  cb = ( cb % 2 ); 

  return local_site(cb_coord, cb_nrow) + cb*vol_cb;
}


/**************** END of QDP functions ****************/
  

#endif


/* Offset by 1 in direction dir */
static void offs(int temp[], const int coord[], int mu, int isign)
{
  int i;

  for(i=0; i < 4; ++i)
    temp[i] = coord[i];

  /* translate address to neighbour */
  temp[mu] = (temp[mu] + isign + 2*getLattSize()[mu]) % getLattSize()[mu];
}

/* This is a 4D parity */
static int parity(const int coord[])
{
  int m;
  int sum = 0;

  for(m=0; m < 4; ++m)
    sum += coord[m];

  return sum % 2;
}


 
  int* make_shift_tables(const int nrow[],
			 void (*getSiteCoords)(int coord[], int node, int linearsite),
			 
			 int (*getLinearSiteIndex)(const int coord[]))
  
{ 
  int dir; 
  int coord[4];
  int linear;
  int Nd = 4;
  int backward=0;
  int forward =1;
  int *shift_table;

  InvTab4 *xinvtab;
  InvTab4 *invtab;

  int x,y,z,t;
  int cboffsets[2] = {0,0};
  int cb;
  int lexico;
  int site;

  int p;

  /* Set the lattice size, get total volume and checkerboarded volume */
  setLattSize(nrow);

  /* Determine what is the starting site for each color (checkerboard) */
  // icolor_start[0] = 0;
  //icolor_start[1] = total_vol_cb;

  /* Allocate the shift table */
  if ((shift_table = (int *)malloc(4*total_vol*2*sizeof(int))) == 0) {
    fprintf(stderr,"init_sse_su3dslash: could not initialize shift_table\n");
    exit(1);
  }

  /* Allocate the site table and the shift table */
  /* Now I want to build the site table */
  /* I want it cache line aligned? */
  xsite_table = (int *)malloc(sizeof(int)*total_vol+63);
  if(xsite_table == 0x0 ) { 
    fprintf(stderr,"Couldnt allocate site table\n");
    exit(1);
  }
  site_table = (int *)((((ptrdiff_t)(xsite_table))+63L)&(-64L));

  /* Loop through sites - you can choose your path below */
  /* This is my local running */
  for(p=0; p < 2; p++) { 
    for(t=0; t < nrow[3]; t++) { 
      for(z=0; z < nrow[2]; z++) {
	for(y=0; y < nrow[1]; y++) { 
	  for(x=0; x < nrow[0]/2; x++) {
	    

	    coord[0] = 2*x+p;
	    coord[1] = y;
	    coord[2] = z; 
	    coord[3] = t;
	    
	    /* Get the site and N-parity of the chosen victim */
	    lexico = getLinearSiteIndex(coord); /* get the lexico index */
	    cb = parity(coord);
	    
	    /* Add lexico site into site_table, for current cb3 and linear */
	    /* Map (cb3, linear) -> lexico */
	    site_table[ cb*total_vol_cb + cboffsets[cb] ] = lexico;

	    /* Next linear site in this checkerboarding */
	    cboffsets[cb]++;
	  }
	}
      }
    }
  }

  /* Get the offsets needed for neighbour comm. */
  /* soffsets(position,direction,isign,cb)   */ 
  /*  where  isign    = +1 : plus direction */
  /*                  =  0 : negative direction */
  /*         cb       =  0 : even lattice (includes origin) */
  /*                  = +1 : odd lattice (does not include origin) */
  /* the offsets cotain the current site, i.e the neighbour for site i  */
  /* is  shift_table(i,dir,cb,mu) and NOT  i + soffset(..)    */
  
  /* Loop over directions and sites, building up shift tables */
  for(site = 0; site < total_vol; ++site) { 
    int fcoord[4], bcoord[4];
    int blinear, flinear;
    int ipos;

    int lexico = site_table[ site ];
    getSiteCoords(coord, 0, lexico); 
    for(dir=0; dir < 4; dir++) {

      /* Backwards displacement*/
      offs(bcoord, coord, dir, -1);
      blinear = getLinearSiteIndex(bcoord);

      /* Forward displacement */
      offs(fcoord, coord, dir, +1);
      flinear = getLinearSiteIndex(fcoord);


      /* Gather */
      shift_table[dir+Nd*site ] = blinear;
      shift_table[dir+Nd*(site+total_vol)] = flinear;
    }
  }

  return shift_table;
} 

int forward_neighbor(int *shift_table, int mysite, int mymu) {
  return shift_table[mymu + 4*(mysite + total_vol)];
}

int backward_neighbor(int *shift_table, int mysite, int mymu) {
  return shift_table[mymu + 4*mysite];
}

void free_shift_tables(int **table) {
  free(*table);
  free(xsite_table);
}

/* Total volume */
int getSubgridVol()
{
  return total_vol;
}

/* Total CB volume */
int getSubgridVolCB()
{
  return total_vol_cb;
}

#ifdef __cplusplus
}
#endif
