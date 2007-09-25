// $Id: t_dslash.cc,v 1.3 2007-09-25 20:24:35 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "unittest.h"

#include "testDslashFull.h"


using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize QDP++ with argc, and argv. Set Lattice Dimensions
  const int latdims[] = {4,2,6,12};

  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, latdims);
  tests.addTest(new testDslashFull(), "testDslashFull" );
  
  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();

}

