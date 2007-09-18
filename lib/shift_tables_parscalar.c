/* $Id: shift_tables_parscalar.c,v 1.2 2007-09-18 17:20:15 bjoo Exp $ */


/* both of these must be called before the P4 dslash is called */
/* computes the scatter/gather table for the P4 32-bit parallel dslash */

/* make_shift_tables(
 * int *shift,               the external table to write the new shift table to
 * int subgrid_vol_cb         # of sites that a call to the dslash will run over per node
 * int nrow[]                original problem size (before checkerboarding)
 * int subgrid_cb_nrow[]     number of sites per subgrid per cb for each direction
 * int bound[]    array of size Nd that holds the number of boundaries per direction
 * int Nd         number of directions
 * )
 * the new shift table will be accessible with hte following index order:
 * shift(
      dir   --> direction
      site  --> site
      forward/backward --> 0 backward, 1 forward
      cb   --> cb
      type --> 0: backward=scatter, forward=gather  1: backward=gather, forward=scatter)

 * decomp uses back scatter type = 0 forw/back = 0
 * decomp_hvv uses forward scatter type = 1 forw/back = 1
 * mvv_recons uses forward gather type = 0 forw/back = 1
 * recons uses back gather type = 1 forw/back = 0

 next version needs to reorder these indices in this file and then propagate to make cb slowest varying and to ditch this forward backward and type stuff (which makes mathematical sense) and go for shift table 0, 1, 2, 3....that might help a little ..... 
*/


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

#include "qmp.h"

/* Max machine size */
#define ND  4

  static int subgrid_vol = -1;
  static int subgrid_vol_cb = -1;

/* Number of dimensions */
static int getNumDim()
{
  return (int)(QMP_get_logical_number_of_dimensions());
}


/* Subgrid lattice size */
static int* getSubgridSize()
{
  static int first = 1;
  static int subgrid_size[ND];

  if (first == 1)
  {
    int i;
    for(i=0; i < getNumDim(); ++i)
      subgrid_size[i] = QMP_get_subgrid_dimensions()[i];

    subgrid_size[0] <<= 1;

    first = 0;
  }
  
  return subgrid_size;
}


/* Total problem size */
static int* getLattSize()
{
  static int first = 1;
  static int tot_size[ND];

  if (first == 1)
  {
    const int* phys_size = QMP_get_logical_dimensions();
    int i;

    for(i=0; i < getNumDim(); ++i)
      tot_size[i] = getSubgridSize()[i]*phys_size[i];

    first = 0;
  }
  
  return tot_size;
}


/* Subgrid volume */
int getSubgridVol()
{
  return subgrid_vol;
}


/* Subgrid volume */
int getSubgridVolCB()
{
  return subgrid_vol_cb;
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

 
/*****************************************************************************************/
/*********** These are functions taken from QDP++ and converted to C *****************/


/* Returns the node number given some logical node coordinate
 * This is not meant to be speedy 
 */
static int getNodeNumberFrom(const int node_coord[]) 
{
  return QMP_get_node_number_from(node_coord);
}

/* Returns the logical node coordinates given some node number
 * This is not meant to be speedy 
 */
static void getLogicalCoordFrom(int node_coord[], int node) 
{
  int* node_crd = QMP_get_logical_coordinates_from(node);  
  int i;

  for(i=0; i < getNumDim(); ++i)
    node_coord[i] = node_crd[i];

  free(node_crd);   
}


/* The linearized site index for the corresponding coordinate
 * This layout is appropriate for a 2 checkerboard (red/black) lattice 
 */
static int getLinearSiteIndex(const int coord[])
{
  int subgrid_vol_cb = getSubgridVol() >> 1;
  int subgrid_cb_nrow[ND];
  int i, cb, m;
  int subgrid_cb_coord[ND];

  for(i=0; i < getNumDim(); ++i)
    subgrid_cb_nrow[i] = getSubgridSize()[i];
  subgrid_cb_nrow[0] >>= 1;

  cb = 0;
  for(m=0; m < getNumDim(); ++m)
    cb += coord[m];
  cb &= 1;

  subgrid_cb_coord[0] = (coord[0] >> 1) % subgrid_cb_nrow[0];
  for(i=1; i < getNumDim(); ++i)
    subgrid_cb_coord[i] = coord[i] % subgrid_cb_nrow[i];
    
  return local_site(subgrid_cb_coord, subgrid_cb_nrow) + cb*subgrid_vol_cb;
}


/* The node number for the corresponding lattice coordinate
 * 
 * This layout is appropriate for a 2 checkerboard (red/black) lattice,
 * but to find the nodeNumber this function resembles a simple lexicographic 
 * layout
 */
static int getNodeNumber(const int coord[])
{
  int tmp_coord[ND];
  int i;

  for(i=0; i < getNumDim(); ++i)
    tmp_coord[i] = coord[i] / getSubgridSize()[i];
    
  return getNodeNumberFrom(tmp_coord);
}

/* Reconstruct the lattice coordinate from the node and site number
 * 
 * This is the inverse of the nodeNumber and linearSiteIndex functions.
 * The API requires this function to be here.
 */
static void getSiteCoords(int coord[], int node, int linearsite)
{
  int subgrid_vol_cb = getSubgridVol() >> 1;
  int subgrid_cb_nrow[ND];
  int i;
  int cbb, cb, tmp_coord[ND];

  for(i=0; i < getNumDim(); ++i)
    subgrid_cb_nrow[i] = getSubgridSize()[i];
  subgrid_cb_nrow[0] >>= 1;

  /* Get the base (origins) of the absolute lattice coord */
  getLogicalCoordFrom(coord, node);
  for(i=0; i < getNumDim(); ++i)
    coord[i] *= getSubgridSize()[i];
    
  cb = linearsite / subgrid_vol_cb;
  crtesn(tmp_coord, linearsite % subgrid_vol_cb, subgrid_cb_nrow);

  /* Add on position within the node
   NOTE: the cb for the x-coord is not yet determined */
  coord[0] += 2*tmp_coord[0];
  for(i=1; i < getNumDim(); ++i)
    coord[i] += tmp_coord[i];

  /*  Determine cb including global node cb */
  cbb = cb;
  for(i=1; i < getNumDim(); ++i)
    cbb += coord[i];
  coord[0] += (cbb & 1);
}
/**************** END of QDP functions ****************/
/*****************************************************************************************/



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


int* make_shift_tables(int icolor_start[2], int bound[2][4][4])
{ 
  volatile int type, cb, dir; 
  const int my_node = QMP_get_node_number();
  int coord[4];
  int linear;
  int Nd = getNumDim();
  int nsize;
  int *shift;

  int i;

  /* Setup the subgrid volume for ever after */
  subgrid_vol = 1;
  for(i=0; i < getNumDim(); ++i) {
    subgrid_vol *= getSubgridSize()[i]; 
  }

  /* Get the checkerboard size for ever after */
  subgrid_vol_cb = subgrid_vol / 2;

  /* Allocate the shift table */
  /* The structure is as follows: There are 4 shift tables in order:

    [ Table 1 | Table 2 | Table 3 | Table 4 ]

    Table 1: decomp_hvv_scatter_index[mu][site]
    Table 2: decomp_scatter_index[mu][site]
    Table 3: recons_mvv_gather_index[mu][site]
    Table 4: recons_gather_index[mu][site]
  
  */
 
  if ((shift = (int *)malloc(Nd*subgrid_vol*4*sizeof(int))) == 0) {
    QMP_error("init_wnxtsu3dslash: could not initialize shift_table");
    QMP_abort(1);
  }


  /* Determine what is the starting site for each color (checkerboard) */
  getSiteCoords(coord, my_node, 0);   /* try linear=0 */
  icolor_start[parity(coord)] = 0;

  getSiteCoords(coord, my_node, subgrid_vol_cb);   /* try linear=subgrid_vol_cb */
  icolor_start[parity(coord)] = subgrid_vol_cb;

  /* Initialize the boundary counters */
  for(cb=0; cb < 2; cb++) {
    for(type=0; type < 4; type++) {
      for(dir=0; dir < Nd; dir++) {
	bound[cb][type][dir] = 0;	
      }
    }
  }


  /* Loop over directions and sites, building up shift tables */
  for(dir=0; dir < Nd; dir++) {
 
    for(linear=0; linear < subgrid_vol; ++linear) {
    
      int fcoord[ND], bcoord[ND];
      int fnode, bnode;
      int blinear, flinear;

      /* Get the global site coords from the node and linear index */
      getSiteCoords(coord, my_node, linear);
      cb = parity(coord);   /* global cb of the coord site */

      /* Backwards displacement*/
      offs(bcoord, coord, dir, -1);
      bnode   = getNodeNumber(bcoord);
      blinear = getLinearSiteIndex(bcoord);

      /* Forward displacement */
      offs(fcoord, coord, dir, +1);
      fnode   = getNodeNumber(fcoord);
      flinear = getLinearSiteIndex(fcoord);

      /* Scatter:  decomp_{plus,minus} */
      /* Operation: a^F(shift(x,type=0),dir) <- decomp(psi(x),dir) */ 
      /* Send backwards - also called a receive from forward */
      type = 0;
      if (bnode != my_node)
	/* append to tail 1, note in table */ 
	shift[dir+Nd*(linear+subgrid_vol*(type))] = subgrid_vol_cb + (bound[1-cb][type][dir])++;
      else
	shift[dir+Nd*(linear+subgrid_vol*(type))] = blinear % subgrid_vol_cb;
    

      /* Scatter:  decomp_hvv_{plus,minus} */
      /* Operation:  a^B(shift(x,type=1),dir) <- U^dag(x,dir)*decomp(psi(x),dir) */
      /* Send forwards - also called a receive from backward */
      type = 1;
      if (fnode != my_node)
	/* Append to tail 1 */
	shift[dir+Nd*(linear+subgrid_vol*(type))] = subgrid_vol_cb + (bound[1-cb][type][dir])++;
      else
	shift[dir+Nd*(linear+subgrid_vol*(type))] = flinear % subgrid_vol_cb;


      /* Gather:  mvv_recons_{plus,minus} */
      /* Operation:  chi(x) <-  \sum_dir U(x,dir)*a^F(shift(x,type=2),dir) */
      /* Receive from forward */
      type = 2;
      if (fnode != my_node)
	/* grab from tail 2 at the right spot */ 
	shift[dir+Nd*(linear+subgrid_vol*(type))] = 2*subgrid_vol_cb + (bound[cb][type][dir])++;
      else
	shift[dir+Nd*(linear+subgrid_vol*(type))] = linear % subgrid_vol_cb;


      /* Gather:  recons_{plus,minus} */
      /* Operation:  chi(x) +=  \sum_dir recons(a^B(shift(x,type=3),dir),dir) */
      /* Receive from backward */
      type = 3;
      if (bnode != my_node)
	/* grab from tail 2 at the right spot */ 
	shift[dir+Nd*(linear+subgrid_vol*(type))] = 2*subgrid_vol_cb + (bound[cb][type][dir])++;
      else
	shift[dir+Nd*(linear+subgrid_vol*(type))] = linear % subgrid_vol_cb;
    } 
  }

  /* Sanity check - make sure the sending and receiving counters match */
  for(cb=0; cb < 2; cb++)
    for(dir=0; dir < Nd; dir++)
    {
/*      for(type=0; type < 4; type++)
	QMP_fprintf(f,"bound[cb=%d][type=%d][dir=%d] = %d",cb,type,dir,bound[cb][type][dir]); */

      if (bound[cb][0][dir] != bound[cb][2][dir])
      {
	QMP_error("SSE Wilson dslash - make_shift_tables: type 0 and 2 send/recv counts do not match: %d %d, cb=%d",
		  bound[cb][0][dir],bound[cb][2][dir],cb);
	QMP_abort(1);
      }

      if (bound[1-cb][0][dir] != bound[cb][0][dir])
      {
	QMP_error("SSE Wilson dslash - make_shift_tables: type 0 diff. cb send/recv counts do not match: %d %d",
		  bound[1-cb][0][dir],bound[cb][0][dir]);
	QMP_abort(1);
      }

      if (bound[cb][1][dir] != bound[cb][3][dir])
      {
	QMP_error("SSE Wilson dslash - make_shift_tables: type 1 and 3 send/recv counts do not match: %d %d",
		  bound[cb][1][dir],bound[cb][3][dir]);
	QMP_abort(1);
      }

      if (bound[1-cb][1][dir] != bound[cb][1][dir])
      {
	QMP_error("SSE Wilson dslash - make_shift_tables: type 1 diff. cb send/recv counts do not match: %d %d",
		  bound[1-cb][1][dir],bound[cb][1][dir]);
	QMP_abort(1);
      }
    }
  
  return shift;

} 

/* decomp plus, decomp minus scatter using these indices */
int decomp_hvv_scatter_index(int *table, int mysite, int mymu) 
{
  return table[mymu+4*(mysite+subgrid_vol)];
}

/*  decomp_hvv_plus, decomp_hvv_minus scatters using these indices */
int decomp_scatter_index(int *table, int mysite, int mymu) 
{
  return table[mymu+4*(mysite)];
}

/*  mvv_recons_plus, gathers from this index */
int recons_mvv_gather_index(int *table, int mysite, int mymu) 
{
  return table[mymu+4*(mysite+subgrid_vol*(2))];
}

/* recons_plus, recons_minux gather from this index */
int recons_gather_index(int *table, int mysite, int mymu) 
{
  return table[mymu+4*(mysite+subgrid_vol*(3))];
}

#ifdef __cplusplus
}
#endif
