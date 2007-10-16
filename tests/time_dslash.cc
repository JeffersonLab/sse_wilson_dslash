// $Id: time_dslash.cc,v 1.2 2007-10-16 11:49:09 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "unittest.h"

#include "timeDslash.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize QDP++ with argc, and argv. Set Lattice Dimensions
  const int latdims[] = {6,6,6,4};

  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, latdims);
  tests.addTest(new timeDslash(), "timeDslash" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();

}

