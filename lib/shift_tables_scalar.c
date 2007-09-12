/* $Id: shift_tables_scalar.c,v 1.1.1.1 2007-09-12 19:33:13 bjoo Exp $

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
  

#ifdef __cplusplus
extern "C" {
#endif

#ifndef DMALLOC
#include <stdlib.h>
#else
#include <dmalloc.h>
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>


/* Max machine size */
#define ND  4

 

/* Number of dimensions */
static int getNumDim()
{
  return (int)(ND);
}


/* Total problem size */
static int tot_size[ND];

static void setLattSize(const int size[])
{
  int i;
  for(i=0; i < getNumDim(); ++i)
    tot_size[i] = size[i];
}

static int* getLattSize()
{
  return tot_size;
}


/* Total volume */
static int getTotalVol()
{
  static int first = 1;
  static int total_vol = 1;

  if (first == 1)
  {
    int i;
    for(i=0; i < getNumDim(); ++i)
      total_vol *= getLattSize()[i];

    first = 0;
  }
  
  return total_vol;
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
  int vol_cb = getTotalVol() >> 1;
  int cb_nrow[ND];

  for(m=0; m < getNumDim(); ++m)
    cb_nrow[m] = getLattSize()[m];
  cb_nrow[0] >>= 1;
  
  cb = linearsite / vol_cb;
  crtesn(coord, linearsite % vol_cb, cb_nrow);

  cbb = cb;
  for(m=1; m < getNumDim(); ++m)
    cbb += coord[m];
  cbb = cbb & 1;

  coord[0] = 2*coord[0] + cbb;
}


//! The linearized site index for the corresponding coordinate
/*! This layout is appropriate for a 2 checkerboard (red/black) lattice */
int getLinearSiteIndex(const int coord[])
{
  int vol_cb = getTotalVol() >> 1;
  int cb_nrow[ND];
  int cb_coord[ND];
  int cb, m;

  for(m=0; m < getNumDim(); ++m)
    cb_nrow[m] = getLattSize()[m];
  cb_nrow[0] >>= 1;
  
  for(m=0; m < getNumDim(); ++m)
    cb_coord[m] = coord[m];
  cb_coord[0] >>= 1;    // Number of checkerboards
    
  cb = 0;
  for(m=0; m < getNumDim(); ++m)
    cb += coord[m];
  cb = cb & 1;

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

static int parity(const int coord[])
{
  int m;
  int sum = 0;

  for(m=0; m < getNumDim(); ++m)
    sum += coord[m];

  return sum & 1;
}


 
void make_shift_tables(int *soffsets, int icolor_start[2], const int nrow[])
{ 
  int dir; 
  int total_vol;
  int total_vol_cb;
  int coord[ND];
  int linear;
  int Nd = getNumDim();
  
  FILE *f;

  /* Set the global value */
  setLattSize(nrow);

  total_vol = getTotalVol();
  total_vol_cb = total_vol >> 1;

  /* Determine what is the starting site for each color (checkerboard) */
  icolor_start[0] = 0;
  icolor_start[1] = total_vol_cb;

//  f = fopen("shift_table.out0", "w");  /* debugging */

#define BACKWARD 0 
#define FORWARD 1 
 
  /* Get the offsets needed for neighbour comm. */
  /* soffsets(position,direction,isign,cb)   */ 
  /*  where  isign    = +1 : plus direction */
  /*                  =  0 : negative direction */
  /*         cb       =  0 : even lattice (includes origin) */
  /*                  = +1 : odd lattice (does not include origin) */
  /* the offsets cotain the current site, i.e the neighbour for site i  */
  /* is  soffsets(i,dir,cb,mu) and NOT  i + soffset(..)    */
  
  /* Loop over directions and sites, building up shift tables */
  for(dir=0; dir < Nd; dir++) 
  { 
    for(linear=0; linear < total_vol; ++linear)
    {
      int fcoord[ND], bcoord[ND];
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

//      fprintf(f,"dir=%d crd=[%d,%d,%d,%d] cb=%d  l=%d fl=%d bl=%d\n",
//	      dir, coord[0],coord[1],coord[2],coord[3],parity(coord),
//	      linear,flinear,blinear);

      /* Gather */
      soffsets[dir+Nd*(linear+total_vol*(BACKWARD))] = blinear;
      soffsets[dir+Nd*(linear+total_vol*( FORWARD))] = flinear;
    }
  }

//  fclose(f);
} 

#ifdef __cplusplus
}
#endif
