// $Id: time_dslash.cc,v 1.3 2007-10-19 13:47:07 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "unittest.h"

#include "timeDslash.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize QDP++ with argc, and argv. Set Lattice Dimensions
  const int latdims[] = {4,4,4,8};

  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, latdims);
  tests.addTest(new timeDslash(), "timeDslash" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();

}

