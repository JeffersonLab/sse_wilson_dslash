/* $Id: shift_tables_scalar.c,v 1.2 2007-09-14 19:32:11 bjoo Exp $

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

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


  /* Number of dimensions */
  static int getNumDim()
  {
    return (int)(4);
  }


/* Total problem size */
  static int tot_size[4];
  static int total_vol = -1;
  static int total_vol_cb = -1;

  static void setLattSize(const int size[])
  {
    int i;
    for(i=0; i < getNumDim(); ++i) {
      tot_size[i] = size[i];
    }
    total_vol = size[0];
    for(i=1; i < getNumDim(); ++i) {
      total_vol *= size[i];
    }

    total_vol_cb = total_vol / 2;
  }

  static int* getLattSize()
  {
    return tot_size;
  }




/* Decompose lexicographic site ipos into lattice coordinates */
static void crtesn(int coord[], int ipos, int latt_size[])
{
  int i;

  for(i=0; i < getNumDim(); ++i)
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

  for(mmu=getNumDim()-1; mmu >= 1; --mmu)
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
  int vol_cb = getTotalVolCB();
  int cb_nrow[4];

  for(m=0; m < getNumDim(); ++m)
    cb_nrow[m] = getLattSize()[m];
  cb_nrow[0] /= 2;
  
  cb = linearsite / vol_cb;
  crtesn(coord, linearsite % vol_cb, cb_nrow);

  cbb = cb;
  for(m=1; m < getNumDim(); ++m)
    cbb += coord[m];
  cbb = (cbb % 2);

  coord[0] = 2*coord[0] + cbb;
}


//! The linearized site index for the corresponding coordinate
/*! This layout is appropriate for a 2 checkerboard (red/black) lattice */
int getLinearSiteIndex(const int coord[])
{
  int vol_cb = getTotalVolCB();
  int cb_nrow[4];
  int cb_coord[4];
  int cb, m;

  for(m=0; m < getNumDim(); ++m)
    cb_nrow[m] = getLattSize()[m];
  cb_nrow[0] /= 2;
  
  for(m=0; m < getNumDim(); ++m)
    cb_coord[m] = coord[m];
  cb_coord[0] /= 2;    // Number of checkerboards
    
  cb = 0;
  for(m=0; m < getNumDim(); ++m)
    cb += coord[m];

  cb = ( cb % 2 ); 

  return local_site(cb_coord, cb_nrow) + cb*vol_cb;
}


/**************** END of QDP functions ****************/
/*****************************************************************************************



/* Offset by 1 in direction dir */
static void offs(int temp[], const int coord[], int mu, int isign)
{
  int i;

  for(i=0; i < getNumDim(); ++i)
    temp[i] = coord[i];

  /* translate address to neighbour */
  temp[mu] = (temp[mu] + isign + 2*getLattSize()[mu]) % getLattSize()[mu];
}

/* This is a 4D parity */
static int parity(const int coord[])
{
  int m;
  int sum = 0;

  for(m=0; m < getNumDim(); ++m)
    sum += coord[m];

  return sum % 2;
}


 
int* make_shift_tables(int icolor_start[2], const int nrow[])
{ 
  int dir; 
  int coord[4];
  int linear;
  int Nd = 4;
  int backward=0;
  int forward =1;
  int *shift_table;

  /* Set the lattice size, get total volume and checkerboarded volume */
  setLattSize(nrow);

  /* Determine what is the starting site for each color (checkerboard) */
  icolor_start[0] = 0;
  icolor_start[1] = total_vol_cb;

  /* Allocate the shift table */
  if ((shift_table = (int *)malloc(4*total_vol*2*sizeof(int))) == 0) {
    fprintf(stderr,"init_sse_su3dslash: could not initialize shift_table\n");
    exit(1);
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
  for(dir=0; dir < Nd; dir++) {

    /* Loop through the sites linearly */
    for(linear=0; linear < total_vol; ++linear)
    {
      int fcoord[4], bcoord[4];
      int blinear, flinear;
      int ipos;

      /* Get the global site coords from the node and linear index */
      getSiteCoords(coord, 0, linear); 

      /* Backwards displacement*/
      offs(bcoord, coord, dir, -1);
      blinear = getLinearSiteIndex(bcoord);

      /* Forward displacement */
      offs(fcoord, coord, dir, +1);
      flinear = getLinearSiteIndex(fcoord);


      /* Gather */
      shift_table[dir+Nd*linear ] = blinear;
      shift_table[dir+Nd*(linear+total_vol)] = flinear;
    }
  }

  return shift_table;
} 

int forward_neighbor(int *shift_table, int mysite, int mymu) {
  return shift_table[mymu + 4*(mysite + total_vol*(1))];
}

int backward_neighbor(int *shift_table, int mysite, int mymu) {
  return shift_table[mymu + 4*(mysite + total_vol*(0))];
}

void free_shift_tables(int **table) {
  free(*table);
}

/* Total volume */
int getTotalVol()
{
  return total_vol;
}

/* Total CB volume */
int getTotalVolCB()
{
  return total_vol_cb;
}

#ifdef __cplusplus
}
#endif
