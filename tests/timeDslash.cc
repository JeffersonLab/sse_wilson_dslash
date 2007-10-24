#include "unittest.h"
#include "timeDslash.h"

#include "qdp.h"
using namespace QDP;

#ifndef DSLASH_M_W_H
#include "dslashm_w.h"
#endif

#ifndef REUNIT_H
#include "reunit.h"
#endif

#include "sse_dslash.h"
#include "sse_dslash_qdp_packer.h"

using namespace Assertions;
using namespace std;

void
timeDslash::run(void) 
{
  LatticeFermion chi, chi2, psi;
  int cb;
  int isign;

  // What we consider to be small enough...
  Double small;
  if ( sizeof(SSEREAL) == 4 ) { 
    small = Double(5.0e-9);
  }
  else {
    // Adjust this...
    small = Double(1.0e-17);
  }

  // Make a random gauge field 
  multi1d<LatticeColorMatrix> u(4);

  for(int mu=0; mu < 4; mu++) { 
    gaussian(u[mu]);
    reunit(u[mu]);
  }

  // Make a random source
  gaussian(psi);

  
  // Initialize the wilson dslash
  init_sse_su3dslash(Layout::lattSize().slice());

  /// Pack the gauge fields
  multi1d<SSEDslash::PrimitiveSU3Matrix> packed_gauge;
  packed_gauge.resize( 4 * Layout::sitesOnNode() );
  SSEDslash::qdp_pack_gauge(u, packed_gauge);
 
  QDPIO::cout << endl;

  // Go through the test cases -- apply SSE dslash versus, QDP Dslash 
  for(int isign=1; isign >= -1; isign -=2) {
    for(int cb=0; cb < 2; cb++) { 
      int source_cb = 1 - cb;
      int target_cb = cb;
      chi = zero;
      chi2 = zero;

      // Apply SSE Dslash
      sse_su3dslash_wilson((SSEREAL *)&(packed_gauge[0]),
			   (SSEREAL *)&(psi.elem(0).elem(0).elem(0).real()),
			   (SSEREAL *)&(chi.elem(0).elem(0).elem(0).real()),
			   isign, source_cb);
      
      // Apply QDP Dslash
      dslash(chi2,u,psi, isign, target_cb);
      
      // Check the difference per number in chi vector
      LatticeFermion diff = chi2 -chi;

      Double diff_norm = sqrt( norm2( diff ) ) 
	/ ( Real(4*3*2*Layout::vol()) / Real(2));
	
      QDPIO::cout << "\t cb = " << source_cb << "  isign = " << isign << "  diff_norm = " << diff_norm << endl;      
      // Assert things are OK...
      assertion( toBool( diff_norm < small ) );

    }
  }

  StopWatch swatch;
  double n_secs = 10;
  int iters=1;
  double time=0;
  QDPIO::cout << endl << "\t Calibrating for " << n_secs << " seconds " << endl;
  do {
    swatch.reset();
    swatch.start();
    for(int i=0; i < iters; i++) { 
       sse_su3dslash_wilson((SSEREAL *)&(packed_gauge[0]),
			   (SSEREAL *)&(psi.elem(0).elem(0).elem(0).real()),
			   (SSEREAL *)&(chi.elem(0).elem(0).elem(0).real()),
			   1, 0);
    }
    swatch.stop();
    time=swatch.getTimeInSeconds();

    // Average time over nodes
    Internal::globalSum(time);
    time /= (double)Layout::numNodes();

    if (time < n_secs) {
      iters *=2;
      QDPIO::cout << "." << flush;
    }
  }
  while ( time < (double)n_secs );
      
  QDPIO::cout << endl;
  QDPIO::cout << "\t Timing with " << iters << " counts" << endl;

  swatch.reset();
  swatch.start();
  
  for(int i=0; i < iters; ++i) {
       sse_su3dslash_wilson((SSEREAL *)&(packed_gauge[0]),
			   (SSEREAL *)&(psi.elem(0).elem(0).elem(0).real()),
			   (SSEREAL *)&(chi.elem(0).elem(0).elem(0).real()),
			   1, 0);

  }
  swatch.stop();
  time=swatch.getTimeInSeconds();

  // Average time over nodes
  Internal::globalSum(time);
  time /= (double)Layout::numNodes();

  QDPIO::cout << "\t " << iters << " iterations in " << time << " seconds " << endl;
  double Mflops = 1390.0f*(double)(iters)*(double)(Layout::vol()/2)/1.0e6;
  double perf = Mflops/time;
  QDPIO::cout << "\t Performance is: " << perf << " Mflops in Total" << endl;
  QDPIO::cout << "\t Performance is: " << perf / (double)Layout::numNodes() << " per MPI Process" << endl;

  // Finalize the Dslash
  free_sse_su3dslash();

}
