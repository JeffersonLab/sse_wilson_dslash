#include "dispatch_scalar.h"

#ifdef __cplusplus
extern "C" { 
#endif


#ifdef DSLASH_USE_QMT_THREADS
  /* Threaded version of the dispatch. We call the qmt_call routine
     with our func, n_sites, and argument */
#include <qmt.h>
#endif

void dispatch_to_threads(void (*func)(size_t, size_t,int, const void *),
			 spinor_array* source,
			 spinor_array* result, 
			 my_mat_array u,
			 int cb,
			 int n_sites)
{
  ThreadWorkerArgs a;

  a.psi = source;
  a.res = result;
  a.u = u;
  a.cb = cb; 
#ifdef DSLASH_USE_QMT_THREADS
  qmt_call((qmt_userfunc_t)func, n_sites, &a);
#else
  (*func)(0, n_sites, 0, &a);
#endif
}



#ifdef __cplusplus
};
#endif
