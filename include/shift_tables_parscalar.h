#ifndef SHIFT_TABLES_PARSCALAR_H
#define SHIFT_TABLES_PARSCALAR_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

  typedef int offset[4];
  
  #define DECOMP_SCATTER 0
  #define DECOMP_HVV_SCATTER 1
  #define RECONS_MVV_GATHER 2
  #define RECONS_GATHER 3

  int getSubgridVol();
  int getSubgridVolCB();
  void make_shift_tables(int icolor_start[2], int bound[2][2][4]);
  void free_shift_tables(void);
  
  int offset_decomp_scatter(int site, int mu);
  int offset_decomp_hvv_scatter(int mu, int site);
  int offset_recons_mvv_gather(int mu, int site);
  int offset_recons_gather(int mu, int site);

#ifdef __cplusplus
};
#endif


#endif
