#ifndef SHIFT_TABLES_PARSCALAR_H
#define SHIFT_TABLES_PARSCALAR_H

#ifdef __cplusplus
extern "C" {
#endif




  
  int decomp_scatter_index(int *table, int mysite, int mymu);
  int decomp_hvv_scatter_index(int *table, int mysite, int mymu);
  int recons_mvv_gather_index(int *table, int mysite, int mymu);
  int recons_gather_index(int *table, int mysite, int mymu);
  
  int getSubgridVol();
  int getSubgridVolCB();
  int* make_shift_tables(int icolor_start[2], int bound[2][4][4]);


#ifdef __cplusplus
};
#endif


#endif
