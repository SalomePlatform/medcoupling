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

#ifndef _MPIACCESSDECTEST_HXX_
#define _MPIACCESSDECTEST_HXX_

#include <cppunit/extensions/HelperMacros.h>

#include <set>
#include <string>
#include <iostream>
#include "mpi.h"

// (ABN]: too many text output in the MPIAccesTest - this renders
// the analysis complicated:
#define MPI_ACCESS_VERBOSE 0
#define debugStream \
    if (!MPI_ACCESS_VERBOSE) {} \
    else std::cout

class MPIAccessDECTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MPIAccessDECTest );
  //  CPPUNIT_TEST( test_AllToAllDECSynchronousPointToPoint ) ;
  CPPUNIT_TEST( test_AllToAllDECAsynchronousPointToPoint ) ;
  //CPPUNIT_TEST( test_AllToAllvDECSynchronousPointToPoint ) ;
  CPPUNIT_TEST( test_AllToAllvDECAsynchronousPointToPoint ) ;
  //CPPUNIT_TEST( test_AllToAllTimeDECSynchronousPointToPoint ) ;
  CPPUNIT_TEST( test_AllToAllTimeDECAsynchronousPointToPoint ) ;
  CPPUNIT_TEST( test_AllToAllvTimeDECSynchronousNative ) ;
  //CPPUNIT_TEST( test_AllToAllvTimeDECSynchronousPointToPoint ) ;
  CPPUNIT_TEST( test_AllToAllvTimeDECAsynchronousPointToPoint ) ;
  //CPPUNIT_TEST( test_AllToAllvTimeDoubleDECSynchronousPointToPoint ) ;
  CPPUNIT_TEST( test_AllToAllvTimeDoubleDECAsynchronousPointToPoint ) ;
  CPPUNIT_TEST_SUITE_END();
  

public:
 
  MPIAccessDECTest():CppUnit::TestFixture(){}
  ~MPIAccessDECTest(){}  
  void setUp(){}
  void tearDown(){}
  void test_AllToAllDECSynchronousPointToPoint() ;
  void test_AllToAllDECAsynchronousPointToPoint() ;
  void test_AllToAllvDECSynchronousPointToPoint() ;
  void test_AllToAllvDECAsynchronousPointToPoint() ;
  void test_AllToAllTimeDECSynchronousPointToPoint() ;
  void test_AllToAllTimeDECAsynchronousPointToPoint() ;
  void test_AllToAllvTimeDECSynchronousNative() ;
  void test_AllToAllvTimeDECSynchronousPointToPoint() ;
  void test_AllToAllvTimeDECAsynchronousPointToPoint() ;
  void test_AllToAllvTimeDoubleDECSynchronousPointToPoint() ;
  void test_AllToAllvTimeDoubleDECAsynchronousPointToPoint() ;

private:
  void test_AllToAllDEC( bool Asynchronous ) ;
  void test_AllToAllvDEC( bool Asynchronous ) ;
  void test_AllToAllTimeDEC( bool Asynchronous ) ;
  void test_AllToAllvTimeDEC( bool Asynchronous , bool UseMPINative ) ;
  void test_AllToAllvTimeDoubleDEC( bool Asynchronous ) ;
  };

// to automatically remove temporary files from disk
class MPIAccessDECTest_TmpFilesRemover
{
public:
  MPIAccessDECTest_TmpFilesRemover() {}
  ~MPIAccessDECTest_TmpFilesRemover();
  bool Register(const std::string theTmpFile);

private:
  std::set<std::string> myTmpFiles;
};

/*!
 *  Tool to print array to stream.
 */
template<class T>
void MPIAccessDECTest_DumpArray (std::ostream & stream, const T* array, const int length, const std::string text)
{
  stream << text << ": {";
  if (length > 0) {
    stream << array[0];
    for (int i = 1; i < length; i++) {
      stream << ", " << array[i];
    }
  }
  stream << "}" << std::endl;
};

#endif
