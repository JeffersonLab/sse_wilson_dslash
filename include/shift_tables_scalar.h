#ifndef SHIFT_TABLES_SCALAR
#define SHIFT_TABLES_SCALAR

#ifdef __cplusplus 
extern "C" {
#endif

  /* These index the shift tables: */
  /* forward_neighbour(site, mu) returns site index of neighbour in fwd dir */ 
  // #define forward_neighbor(table,mysite,mymu)  table[mymu + 4*(mysite + total_vol)]

  /* backward_neighbor(site, mu) return site index of backward neighbour */
  //#define backward_neighbor(table,mysite,mymu) table[mymu + 4*mysite]

  int getTotalVol();
  int getTotalVolCB();
  int forward_neighbor(int *soffsets, int mysite, int mymu);
  int backward_neighbor(int *soffsets, int mysite, int mymu);
  int* make_shift_tables(int icolor_start[2], const int lat_size[4]);
  void free_shift_tables(int **table);
#ifdef __cplusplus
};
#endif

#endif
