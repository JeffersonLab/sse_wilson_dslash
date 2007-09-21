#include "dispatch_parscalar.h"

#ifdef __cplusplus
extern "C" { 
#endif


#ifdef DSLASH_USE_QMT_THREADS
  /* Threaded version of the dispatch. We call the qmt_call routine
     with our func, n_sites, and argument */
#include <qmt.h>

void dispatch_to_threads(void (*func)(size_t, size_t,int, const void *),
			 spinor_array* the_spinor,
			 halfspinor_array* the_halfspinor, 
			 my_mat_array u,
			 int cb,
			 int n_sites)
{
  ThreadWorkerArgs a;

  a.spinor = the_spinor;
  a.half_spinor = the_halfspinor;
  a.u = u;
  a.cb = cb; 
  qmt_call((qmt_userfunc_t)func, n_sites, &a);
}

#else
  /* Unthreaded dispatch. We call the function directly. The 'low i
ndex'
     is the first site, and the 'thread' should do all the sites */
void dispatch_to_threads(void (*func)(size_t, size_t, int, const void*),
			 spinor_array* the_spinor,
			 halfspinor_array* the_halfspinor, 
			 my_mat_array u,
			 int cb,
			 int n_sites)
{
  ThreadWorkerArgs a;

  a.spinor = the_spinor;
  a.half_spinor = the_halfspinor;
  a.u = u;
  a.cb = cb; 
  (*func)(0, n_sites, 0, &a);
}

#endif






#ifdef __cplusplus
};
#endif
