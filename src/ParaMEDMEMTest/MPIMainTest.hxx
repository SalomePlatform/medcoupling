// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#ifndef _MPIMAINTEST_HXX_
#define _MPIMAINTEST_HXX_

#include <cppunit/CompilerOutputter.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TextTestProgressListener.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestRunner.h>
#include <stdexcept>

#include <mpi.h>

#include <iostream>
#include <fstream>
#ifndef WIN32
#include <fpu_control.h>
#endif

// ============================================================================
/*!
 *  Main program source for Unit Tests with cppunit package does not depend
 *  on actual tests, so we use the same for all partial unit tests.
 */
// ============================================================================

int main(int argc, char* argv[])
{
#ifndef WIN32
  fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM);
  _FPU_SETCW(cw);
#endif
  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  // --- Create the event manager and test controller
  CPPUNIT_NS::TestResult controller;

  // ---  Add a listener that colllects test result
  CPPUNIT_NS::TestResultCollector result;
  controller.addListener( &result );        

  // ---  Add a listener that print dots as test run.
#ifdef WIN32
  CPPUNIT_NS::TextTestProgressListener progress;
#else
  CPPUNIT_NS::BriefTestProgressListener progress;
#endif
  controller.addListener( &progress );      

  // ---  Get the top level suite from the registry

  CPPUNIT_NS::Test *suite =
    CPPUNIT_NS::TestFactoryRegistry::getRegistry().makeTest();

  // ---  Adds the test to the list of test to run

  CPPUNIT_NS::TestRunner runner;
  runner.addTest( suite );
  runner.run( controller);

  // ---  Print test in a compiler compatible format.

  std::ostringstream testFileName;
  testFileName<<"UnitTestsResult"<<rank;
  std::ofstream testFile;
  testFile.open(testFileName.str().c_str(), std::ios::out |  std::ios::trunc);
  //CPPUNIT_NS::CompilerOutputter outputter( &result, std::cerr );
  CPPUNIT_NS::CompilerOutputter outputter( &result, testFile );
  outputter.write(); 

  // ---  Run the tests.

  bool wasSucessful = result.wasSuccessful();
  testFile.close();

  // ---  Return error code 1 if the one of test failed.

  MPI_Finalize();
  
  return wasSucessful ? 0 : 1;
}

#endif
