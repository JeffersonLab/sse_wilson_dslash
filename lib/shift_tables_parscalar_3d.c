

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "qmp.h"
#include "shift_tables_parscalar.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* Unaligned Offset Table */
  static int* xoffset_table_body_3d;

  /* EXPORTED: Aligned Offset Table */
  int* offset_table_body_3d;

  typedef struct { 
    int cb3;
    int linearcb3;
  } InvTab;

  /* Unaligned Site Table */
  static int* xsite_table_3d;

  /* EXPORTED: Aligned Site Table */
  int* site_table_3d;



  /* Total Subgrid Volume on Core */
  /* EXPORTED:                    */
  int subgrid_vol_3d = -1;

  /* CB Volume on Core */
  /* EXPORTED */
  int subgrid_vol_cb_3d = -1;

  /* The number of dimensions once and for all */
  /* EXPORTED */
  int Nd3 = 3;


  /* Number of dimensions */
  static int getNumDim()
  {
    return (int)(QMP_get_logical_number_of_dimensions());
  }


  /* Subgrid lattice size */
  static int* getSubgridSize()
  {
    static int first = 1;
    static int subgrid_size[4];

    if (first == 1) {
      int i;
      for(i=0; i < getNumDim(); ++i) {
	subgrid_size[i] = QMP_get_subgrid_dimensions()[i];
      }

      /* Why do we multiply by 2 after QMP_get_subgrid_dimensions() */
      subgrid_size[0]  *= 2;
      
      first = 0;
    }
  
    return subgrid_size;
  }


  /* Total problem size */
  /* File scope so no conflict with 4D */
  static int* getLattSize()
  {
    static int first = 1;
    static int tot_size[4];

    if (first == 1) {
      
      const int* phys_size = QMP_get_logical_dimensions();
      int i;
      
      for(i=0; i < getNumDim(); ++i) {
	tot_size[i] = getSubgridSize()[i]*phys_size[i];
      }
      
      first = 0;
    }
    
    return tot_size;
  }

  


  /* Decompose lexicographic site ipos into lattice coordinates */
  static void crtesn(int coord[], int ipos, int latt_size[])
  {
    int i;
    
    for(i=0; i < getNumDim(); ++i) {
      coord[i] = ipos % latt_size[i];
      ipos /= latt_size[i];
    }
  }

  /* Calculates the lexicographic site index from the coordinate of a site */
  static int local_site(int coord[], int latt_size[]) {
    int order = 0;
    int mmu;
    
    for(mmu=getNumDim()-1; mmu >= 1; --mmu) {
      order = latt_size[mmu-1]*(coord[mmu] + order);
    }

    order += coord[0];
    
    return order;
  }

 
  /*****************************************************************************************/
  /*********** These are functions taken from QDP++ and converted to C *****************/
  
  /* File scope so no link time conflict with 4D */
  /* Returns the node number given some logical node coordinate
   * This is not meant to be speedy 
   */
  static int getNodeNumberFrom(const int node_coord[]) 
  {
    return QMP_get_node_number_from(node_coord);
  }

  /* File scope so no link time conflict with 4D */
  /* Returns the logical node coordinates given some node number
   * This is not meant to be speedy 
   */
  static void getLogicalCoordFrom(int node_coord[], int node) 
  {
    int* node_crd = QMP_get_logical_coordinates_from(node);  
    int i;
    
    for(i=0; i < getNumDim(); ++i) {
      node_coord[i] = node_crd[i];
    }

    free(node_crd);   
  }

  /* File scope so no link time conflict with 4D */
  /* The linearized site index for the corresponding coordinate
   * This layout is appropriate for a 2 checkerboard (red/black) lattice 
   */
  static int getLinearSiteIndex(const int coord[])
  {
    
    int subgrid_cb_nrow[4];
    int i, cb, m;
    int subgrid_cb_coord[4];

    for(i=0; i < getNumDim(); ++i) {
      subgrid_cb_nrow[i] = getSubgridSize()[i];
    }
    subgrid_cb_nrow[0] /= 2;

    /* Compute 4D Checkerboard */
    cb = 0;
    for(m=0; m < getNumDim(); ++m) {
      cb += coord[m];
    }
    cb &= 1;
    

    /* Taking remainders here would suggest that coord can be 'off node' */
    subgrid_cb_coord[0] = (coord[0] / 2) % subgrid_cb_nrow[0];
    for(i=1; i < getNumDim(); ++i) {
      subgrid_cb_coord[i] = coord[i] % subgrid_cb_nrow[i];
    }

    return local_site(subgrid_cb_coord, subgrid_cb_nrow) + cb*subgrid_vol_cb_3d;
  }

  /* File scope so no link time conflict with 4D */
  /* The node number for the corresponding lattice coordinate
   * 
   * This layout is appropriate for a 2 checkerboard (red/black) lattice,
   * but to find the nodeNumber this function resembles a simple lexicographic 
   * layout
   */
  static int getNodeNumber(const int coord[])
  {
    int tmp_coord[4];
    int i;
    
    for(i=0; i < getNumDim(); ++i) {
      tmp_coord[i] = coord[i] / getSubgridSize()[i];
    }
    return getNodeNumberFrom(tmp_coord);
  }
  
  /* File scope so no link time conflict with 4D */
  /* Reconstruct the lattice coordinate from the node and site number
   * 
   * This is the inverse of the nodeNumber and linearSiteIndex functions.
   * The API requires this function to be here.
   */
  static void getSiteCoords(int coord[], int node, int linearsite)
  {
    int subgrid_cb_nrow[4];
    int i;
    int cbb, cb, tmp_coord[4];
    
    for(i=0; i < getNumDim(); ++i) {
      subgrid_cb_nrow[i] = getSubgridSize()[i];
    }
    subgrid_cb_nrow[0] /= 2;
    
    /* Get the base (origins) of the absolute lattice coord */
    getLogicalCoordFrom(coord, node);
    for(i=0; i < getNumDim(); ++i) {
      coord[i] *= getSubgridSize()[i];
    }
    cb = linearsite / subgrid_vol_cb_3d;
    crtesn(tmp_coord, linearsite % subgrid_vol_cb_3d, subgrid_cb_nrow);

    /* Add on position within the node
       NOTE: the cb for the x-coord is not yet determined */
    coord[0] += 2*tmp_coord[0];
    for(i=1; i < getNumDim(); ++i) {
      coord[i] += tmp_coord[i];
    }

    /*  Determine cb including global node cb - NB This is a 4D Checkerboard
        we are still using a 4D Layout */
    cbb = cb;
    for(i=1; i < getNumDim(); ++i) {
      cbb += coord[i];
    }
    coord[0] += (cbb & 1);
  }
  /**************** END of QDP functions ****************/
  /*****************************************************************************************/



  /* Offset by 1 in direction dir */
  static void offs(int temp[], const int coord[], int mu, int isign)
  {
    int i;
    
    for(i=0; i < getNumDim(); ++i) {
      temp[i] = coord[i];
    }
    /* translate address to neighbour */
    temp[mu] = (temp[mu] + isign + 2*getLattSize()[mu]) % getLattSize()[mu];
  }


  /* This is now the Nd3 parity */
  static int parity(const int coord[])
  {
    int m;
    int sum = 0;
    
    for(m=0; m < Nd3; ++m)
      sum += coord[m];
    
    return sum & 1;
  }


  /* This one does al the hard work */
  void make_shift_tables_3d(int bound[2][4][3])
  { 
    volatile int dir; 
    const int my_node = QMP_get_node_number();
    int coord[4];
    int linear;
    int lexico;
    int **shift_table;
    int x,y,z,t;
    int *subgrid_size = getSubgridSize();
    int cboffsets[2] = {0,0};
    int cb3;
    
    int i;
    int site;
    InvTab *xinvtab;
    InvTab *invtab;
        
    /* Setup the subgrid volume for ever after */
    /* Nb getNumDim() here returns the number of QMP dimensions */
    /* This is potentially undesirable - maybe I should just hardwire this to be 4 */
    /* This sets up a global, allowing subgrid_vol_cb3 to be used ever after -- Yuck! */
    subgrid_vol_3d = 1;
    for(i=0; i < getNumDim(); ++i) {
      subgrid_vol_3d *= getSubgridSize()[i]; 
    }
  
  /* Get the checkerboard size for ever after */
  /* Again this sets up a global. Yuck */
  subgrid_vol_cb_3d = subgrid_vol_3d / 2;




  /* Now I want to build the site table */
  /* I want it cache line aligned? */
  xsite_table_3d = (int *)malloc(sizeof(int)*subgrid_vol_3d+63);
  if(xsite_table_3d == 0x0 ) { 
    QMP_error("Couldnt allocate site table");
    QMP_abort(1);
  }
  site_table_3d = (int *)((((ptrdiff_t)(xsite_table_3d))+63L)&(-64L));

  xinvtab = (InvTab *)malloc(sizeof(InvTab)*subgrid_vol_3d+63);
  if(xinvtab == 0x0 ) { 
    QMP_error("Couldnt allocate site table");
    QMP_abort(1);
  }
  invtab = (InvTab *)((((ptrdiff_t)(xinvtab))+63L)&(-64L));

  /* Loop through sites - you can choose your path below */
  /* For now let's go lexicographically */  
  for(t=0; t < subgrid_size[3]; t++) { 
    for(z=0; z < subgrid_size[2]; z++) {
      for(y=0; y < subgrid_size[1]; y++) { 
	for(x=0; x < subgrid_size[0]; x++) {
	  coord[0] = x;
	  coord[1] = y;
	  coord[2] = z; 
	  coord[3] = t;

	  /* Get the site and N-parity of the chosen victim */
	  lexico = getLinearSiteIndex(coord); /* get the lexico index */
	  cb3 = parity(coord); /* Get cb3 index */

	  /* Add lexico site into site_table, for current cb3 and linear */
	  /* Map (cb3, linear) -> lexico */
	  site_table_3d[ cb3*subgrid_vol_cb_3d + cboffsets[cb3] ] = lexico;

	  /* Now the inverse map: note the (cb, linear) combination for the
	     lexico coordinates */
	  /* Map lexico -> (cb3, linear) */
	  invtab[ lexico ].cb3 = cb3;
	  invtab[ lexico ].linearcb3 = cboffsets[cb3];

	  /* Next linear site in this checkerboarding */
	  cboffsets[cb3]++;
	}
      }
    }
  }

  /* Allocate the shift table */
  /* The structure is as follows: There are 4 shift tables in order:
   
    [ Table 1 | Table 2 | Table 3 | Table 4 ]
    Table 1: decomp_scatter_index[mu][site]
    Table 2: decomp_hvv_scatter_index[mu][site]
    Table 3: recons_mvv_gather_index[mu][site]
    Table 4: recons_gather_index[mu][site]
  
    Each table is indexed as:
    table[ type ][ dir-site ]  

    where dir-site is a combination of
    direction and site indices with direction running fastest.
    
    dir_site is computed as dir + Nd3*site 

    site runs as linear+subgrid_vol_cb*cb3 - although 
    this is also rolled into a linear loop over all the sites.
    
    The values stored in the index table:
      i) If the appropriate neighbour is on-site, we give the 
      linear part of its (cb3,linear) coordinate. The CB3 is 
      implied (essentially its a neigbour so it is opposite ours)
      or it is a received buffer, and may be the same as ours.
      In any case, the stored index will only index something
      on one checkerboard and so we don't need to supply the 
      checkerboard info.

      ii) If the neighbour is off site, we index into a tail.
      
      Our temporaries will be arranged as:

      [  body half spinors ][ Tail 1 half spinors ][ Tail 2 half spinors ] 
      Tail 1 gets sent backward and tail 2 receives from forward 
      or vice versa depending on parity. In any case DECOMPs always go to
      Tail 1, and RECONS-s always go to Tail 2.
      
      We count the sizes of the tails as we go along, and fill out the 
      bound array. Typically tail 1 indices look like subgrid_vol_cb+x
      While tail 2 indices look like 2*subgrid_vol_cb + x
  */
 
  /* 4 for the four types above: */
  if ((shift_table = (int **)malloc(4*sizeof(int*))) == 0 ) {
    QMP_error("init_wnxtsu3dslash: could not initialize shift_table");
    QMP_abort(1);
    
  }
  
  /* Now the table for each of the 4 types */
  for(i=0; i < 4; i++) { 
    /* Nd3 for the  directions */
    if ((shift_table[i] = (int *)malloc(Nd3*subgrid_vol_3d*sizeof(int))) == 0) {
      QMP_error("init_wnxtsu3dslash: could not initialize shift_table");
      QMP_abort(1);
    }
  }


  /* Initialize the boundary counters */
  for(cb3=0; cb3 < 2; cb3++) {
    for(dir=0; dir < Nd3 ; dir++) {
      bound[cb3][0][dir] = 0;	
      bound[cb3][1][dir] = 0;	
      bound[cb3][2][dir] = 0;	
      bound[cb3][3][dir] = 0;	
    }
  }



  /* Loop over All sites */
  for(site=0; site < subgrid_vol_3d; ++site) {

    /* Lookup the site and get its coordinates and parity */
    linear = site_table_3d[site]; 

      /* Get the global site coords from the node and linear index */
    getSiteCoords(coord, my_node, linear);
    cb3 = parity(coord);   /* GLOBAL cb of the coord site */


    /* Loop through the directions */
    for(dir=0; dir < Nd3; dir++) {
      
      int fcoord[4], bcoord[4];
      int fnode, bnode;
      int blinear, flinear;
      
      /* Find forward and backward neighbours. Get both 
	 the node info and the linear coordinates. 
	 NB The linear coordinates are the absolute 
	 lexico ones. */
      
      /* Backwards displaced coordinate & Node */
      offs(bcoord, coord, dir, -1);
      bnode   = getNodeNumber(bcoord);
      blinear = getLinearSiteIndex(bcoord);
 
      /* Forward displaced coordinate & Node */
      offs(fcoord, coord, dir, +1);
      fnode   = getNodeNumber(fcoord);
      flinear = getLinearSiteIndex(fcoord);

      
      /* Scatter:  decomp_{plus,minus} */
      /* Operation: a^F(shift(x,type=0),dir) <- decomp(psi(x),dir) */ 
      /* Send backwards - also called a receive from forward */

      if (bnode != my_node) {      
	/* Offnode */
	/* Append to Tail 1, increase boundary count */
	shift_table[DECOMP_SCATTER][dir+Nd3*site] 
	  = subgrid_vol_cb_3d + bound[1-cb3][DECOMP_SCATTER][dir];
	
	bound[1-cb3][DECOMP_SCATTER][dir]++;
	
      }
      else {                                           
	/* On node. Note the linear part of its (cb3, linear) bit,
	   using a reverse lookup */
	shift_table[DECOMP_SCATTER][dir+Nd3*site] = 
	  invtab[blinear].linearcb3;
      }
      
      /* Scatter:  decomp_hvv_{plus,minus} */
      /* Operation:  a^B(shift(x,type=1),dir) <- U^dag(x,dir)*decomp(psi(x),dir) */
      /* Send forwards - also called a receive from backward */
      if (fnode != my_node) {
	/* Offnode */
	/* Append to Tail 1, increase boundary count */
	shift_table[DECOMP_HVV_SCATTER][dir+Nd3*site]           
	  = subgrid_vol_cb_3d + bound[1-cb3][DECOMP_HVV_SCATTER][dir];
	
	bound[1-cb3][DECOMP_HVV_SCATTER][dir]++;                  
	
      }
      else {
	/* On node. Note the linear part of its (cb3, linear) bit,
	   using a reverse lookup */
	shift_table[DECOMP_HVV_SCATTER][dir+Nd3*site]           /* Onnode */
	  = invtab[flinear].linearcb3 ;
      }
      
      /* Gather:  mvv_recons_{plus,minus} */
      /* Operation:  chi(x) <-  \sum_dir U(x,dir)*a^F(shift(x,type=2),dir) */
      /* Receive from forward */
      if (fnode != my_node) {
	/* Offnode */
	/* Append to Tail 2, increase boundary count */
	
	shift_table[RECONS_MVV_GATHER][dir+Nd3*site] =
	  2*subgrid_vol_cb_3d + (bound[cb3][RECONS_MVV_GATHER][dir]);
	
	bound[cb3][RECONS_MVV_GATHER][dir]++;
	
      }
      else {
	/* On node. Note the linear part of its (cb3, linear) bit,
	   using a reverse lookup. Note this is a recons post shift,
	   so the linear coordinate to invert is mine rather than the neighbours */
	shift_table[RECONS_MVV_GATHER][dir+Nd3*site] =
	  invtab[linear].linearcb3 ;
      }
      
      
      /* Gather:  recons_{plus,minus} */
      /* Operation:  chi(x) +=  \sum_dir recons(a^B(shift(x,type=3),dir),dir) */
      /* Receive from backward */
      if (bnode != my_node) {
	
	shift_table[RECONS_GATHER][dir+Nd3*site] = 
	  2*subgrid_vol_cb_3d + bound[cb3][RECONS_GATHER][dir];
	
	bound[cb3][RECONS_GATHER][dir]++;
      }
      else {
	/* On node. Note the linear part of its (cb3, linear) bit,
	   using a reverse lookup. Note this is a recons post shift,
	   so the linear coordinate to invert is mine rather than the neighbours */

	shift_table[RECONS_GATHER][dir+Nd3*site] = 
	  invtab[linear].linearcb3 ;
      }
      
    } 
  }
  


  /* Sanity check - make sure the sending and receiving counters match */
  for(cb3=0; cb3 < 2; cb3++) {
    for(dir=0; dir < Nd3; dir++) {
      /* Offnode */
      /* Append to Tail 2, increase boundary count */
      
      /* Sanity 1: Must have same number of boundary sites on each cb for 
	 a given operation */
      for(i = 0; i < 4; i++) { 
	if (bound[1-cb3][i][dir] != bound[cb3][i][dir]) {
	  
	  QMP_error("SSE Wilson dslash - make_shift_tables: type 0 diff. cb send/recv counts do not match: %d %d",
		    bound[1-cb3][i][dir],bound[cb3][i][dir]);
	  QMP_abort(1);
	}
      }
    }
  }
  
    
  /* Now I want to make the offset table into the half spinor temporaries */
  /* The half spinor temporaries will look like this:

    dir=0 [ Body Half Spinors ][ Tail 1 Half Spinors ][ Tail 2 Half Spinors ]
    dir=1 [ Body Half Spinors ][ Tail 1 Half Spinors ][ Tail 2 Half Spinors ]
    ...

    And each of these blocks of half spinors will be sized to vol_cb
    sites (ie half volume only).  The shift_table() for a given site and
    direction indexes into one of these lines. So the offset table essentially
    delineates which line one picks, by adding an offset of 
      3*subgrid_vol_cb*dir 
    To the shift. The result from offset table, can be used directly as a
    pointer displacement on the temporaries.

  */
  xoffset_table_body_3d = (int *)malloc(Nd3*4*subgrid_vol_3d*sizeof(int)+63);
  if( xoffset_table_body_3d == 0 ) {
    QMP_error("init_wnxtsu3dslash: could not initialize offset_table[i]");
    QMP_abort(1);
  }
  /* This is the bit what aligns straight from AMD Manual */
  offset_table_body_3d = (int *)((((ptrdiff_t)(xoffset_table_body_3d)) + 63L) & (-64L));

  for(site=0; site < subgrid_vol_3d; site++) { 
    for(dir=0; dir < Nd3; dir++) { 
      
      offset_table_body_3d[ dir + Nd3*( site + subgrid_vol_3d*DECOMP_SCATTER ) ] = 
	shift_table[DECOMP_SCATTER][dir+Nd3*site ]+3*subgrid_vol_cb_3d*dir;
      
      offset_table_body_3d[ dir + Nd3*( site + subgrid_vol_3d*DECOMP_HVV_SCATTER ) ] = 
	shift_table[DECOMP_HVV_SCATTER][dir+Nd3*site ]+3*subgrid_vol_cb_3d*dir;
      
      offset_table_body_3d[ dir + Nd3*( site + subgrid_vol_3d*RECONS_MVV_GATHER ) ] = 
	shift_table[RECONS_MVV_GATHER][dir+Nd3*site ]+3*subgrid_vol_cb_3d*dir;
      
      offset_table_body_3d[ dir + Nd3*( site + subgrid_vol_3d*RECONS_GATHER ) ] = 
	shift_table[RECONS_GATHER][dir+Nd3*site ]+3*subgrid_vol_cb_3d*dir;
      
    }
  }



  /* Free shift table - it is no longer needed. We deal solely with offsets */
  for(i=0; i < 4; i++) { 
    free( (shift_table[i]) );
  }
  free( shift_table );

  /* Free the inverse site table, it is no longer needed */
  free( xinvtab );
  
} 


void free_shift_tables_3d(void) 
{
  /* Free the offset and site tables - free their aligned variants */
  free( xoffset_table_body_3d );
  free( xsite_table_3d );

}

#ifdef __cplusplus
}
#endif
