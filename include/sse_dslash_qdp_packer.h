#ifndef SSE_DSLASH_QDP_PACKER_H
#define SSE_DSLASH_QDP_PACKER_H

#include "sse_config.h"

#ifndef QDP_INCLUDE
#include "qdp.h"
#endif 

namespace SSEDslash { 

  typedef PColorMatrix<RComplex<REAL>, 3> PrimitiveSU3Matrix;


  void qdp_pack_gauge(const multi1d<LatticeColorMatrix>&_u, multi1d<PrimitiveSU3Matrix>& u_tmp);

};

#endif
