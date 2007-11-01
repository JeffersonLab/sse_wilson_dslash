// $Id: time_dslash_3d.cc,v 1.1 2007-11-01 15:20:18 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "unittest.h"

#include "timeDslash3D.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize QDP++ with argc, and argv. Set Lattice Dimensions
  const int latdims[] = {4,4,4,8};

  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, latdims);
  tests.addTest(new timeDslash3D(), "timeDslash3D" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();

}

