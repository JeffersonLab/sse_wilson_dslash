// $Id: t_dslash_3d.cc,v 1.1 2007-10-31 15:11:16 bjoo Exp $

#include <iostream>
#include <cstdio>

#include "qdp.h"
#include "unittest.h"

#include "testDslash3D.h"

using namespace QDP;

int main(int argc, char **argv)
{
  // Initialize QDP++ with argc, and argv. Set Lattice Dimensions
  const int latdims[] = {4,4,4,8};

  // Initialize UnitTest jig
  TestRunner  tests(&argc, &argv, latdims);
  tests.addTest(new testDslash3D(), "testDslash3D" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();

}

