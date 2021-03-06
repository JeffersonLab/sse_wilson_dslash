#ifndef SITE_DSLASH_64BIT_SCALAR
#define SITE_DSLASH_64BIT_SCALAR

#include "types64.h"
#include "sse_align.h"

#ifdef __cplusplus
extern "C" { 
#endif


  void dslash_plus_dir0_forward(spinor_array  spinor_in,
				u_mat_array  u,
				spinor_array spinor_out);

  void dslash_plus_dir0_backward_add( spinor_array  spinor_in,
				 u_mat_array  u,
				spinor_array spinor_out);

  void dslash_plus_dir1_forward_add(spinor_array  spinor_in,
				    u_mat_array  u,
				    spinor_array spinor_out);  
  
  void dslash_plus_dir1_backward_add(  spinor_array  spinor_in,
				       u_mat_array  u,
				       spinor_array spinor_out);

  void dslash_plus_dir2_forward_add(  spinor_array  spinor_in,
				      u_mat_array  u,
				      spinor_array spinor_out);  
  
  void dslash_plus_dir2_backward_add(  spinor_array  spinor_in,
				       u_mat_array  u,
				       spinor_array spinor_out);

#if 0
  void dslash_plus_dir2_backward_add_store(  spinor_array  spinor_in,
					 u_mat_array  u,
					 spinor_array spinor_out);
#endif

  void dslash_plus_dir3_forward_add(  spinor_array  spinor_in,
				      u_mat_array  u,
				      spinor_array spinor_out);
  
  void dslash_plus_dir3_backward_add_store(  spinor_array  spinor_in,
					     u_mat_array  u,
					     spinor_array spinor_out);
  
  void dslash_minus_dir0_forward(  spinor_array  spinor_in,
				   u_mat_array  u,
				   spinor_array spinor_out);
  
  void dslash_minus_dir0_backward_add(  spinor_array  spinor_in,
					u_mat_array  u,
					spinor_array spinor_out);
  
  void dslash_minus_dir1_forward_add(  spinor_array  spinor_in,
				       u_mat_array  u,
				       spinor_array spinor_out);

  
  void dslash_minus_dir1_backward_add(  spinor_array  spinor_in,
					u_mat_array  u,
					spinor_array spinor_out);
  
  void dslash_minus_dir2_forward_add(  spinor_array  spinor_in,
				       u_mat_array  u,
				       spinor_array spinor_out);
  
  void dslash_minus_dir2_backward_add(  spinor_array  spinor_in,
					u_mat_array  u,
					spinor_array spinor_out);
#if 0
  void dslash_minus_dir2_backward_add_store(  spinor_array  spinor_in,
					      u_mat_array  u,
					      spinor_array spinor_out);
#endif
  
  void dslash_minus_dir3_forward_add(  spinor_array  spinor_in,
				       u_mat_array  u,
				       spinor_array spinor_out);
  
  void dslash_minus_dir3_backward_add_store(  spinor_array  spinor_in,
					      u_mat_array  u,
					      spinor_array spinor_out);  
  

#ifdef __cplusplus 
};
#endif

#endif
