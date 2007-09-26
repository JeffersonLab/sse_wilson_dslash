// $Id: t_dslash.cc,v 1.4 2007-09-26 20:45:10 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "unittest.h"

#include "testDslashFull.h"
#include "testDecomp.h"
#include "testDecompHvv.h"


using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize QDP++ with argc, and argv. Set Lattice Dimensions
  const int latdims[] = {4,2,6,12};

  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, latdims);
  tests.addTest(new testDslashFull(), "testDslashFull" );
  tests.addTest(new testDecomp0Minus(), "testDecomp0Minus" );
  tests.addTest(new testDecomp1Minus(), "testDecomp1Minus" );
  tests.addTest(new testDecomp2Minus(), "testDecomp2Minus" );
  tests.addTest(new testDecomp3Minus(), "testDecomp3Minus" );
  tests.addTest(new testDecomp0Plus(), "testDecomp0Plus" );
  tests.addTest(new testDecomp1Plus(), "testDecomp1Plus" );  
  tests.addTest(new testDecomp2Plus(), "testDecomp2Plus" );
  tests.addTest(new testDecomp3Plus(), "testDecomp3Plus" );
  tests.addTest(new testDecompHvv0Plus(), "testDecompHvv0Plus");
  // tests.addTest(new testDecompHvv1Plus(), "testDecompHvv1Plus");
  // tests.addTest(new testDecompHvv2Plus(), "testDecompHvv2Plus");
  // tests.addTest(new testDecompHvv3Plus(), "testDecompHvv3Plus");
  // tests.addTest(new testDecompHvv0Minus(), "testDecompHvv0Minus");
  // tests.addTest(new testDecompHvv1Minus(), "testDecompHvv1Minus");
  // tests.addTest(new testDecompHvv2Minus(), "testDecompHvv2Minus");
  // tests.addTest(new testDecompHvv3Minus(), "testDecompHvv3Minus");
  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();

}

