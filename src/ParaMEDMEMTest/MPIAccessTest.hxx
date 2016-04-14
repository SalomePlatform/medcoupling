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

#ifndef _MPIACCESSTEST_HXX_
#define _MPIACCESSTEST_HXX_

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

class MPIAccessTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MPIAccessTest );
  CPPUNIT_TEST( test_MPI_Access_Send_Recv ) ;
  CPPUNIT_TEST( test_MPI_Access_Cyclic_Send_Recv ) ;
  CPPUNIT_TEST( test_MPI_Access_SendRecv ) ;
  CPPUNIT_TEST( test_MPI_Access_ISend_IRecv ) ;
  CPPUNIT_TEST( test_MPI_Access_Cyclic_ISend_IRecv ) ;
  CPPUNIT_TEST( test_MPI_Access_ISendRecv ) ;
  CPPUNIT_TEST( test_MPI_Access_Probe ) ;
  CPPUNIT_TEST( test_MPI_Access_IProbe ) ;
  CPPUNIT_TEST( test_MPI_Access_Cancel ) ;
  CPPUNIT_TEST( test_MPI_Access_Send_Recv_Length ) ;
  CPPUNIT_TEST( test_MPI_Access_ISend_IRecv_Length ) ;
  CPPUNIT_TEST( test_MPI_Access_ISend_IRecv_Length_1 ) ;
  CPPUNIT_TEST( test_MPI_Access_Time ) ;
  CPPUNIT_TEST( test_MPI_Access_Time_0 ) ;
  CPPUNIT_TEST( test_MPI_Access_ISend_IRecv_BottleNeck ) ;
  CPPUNIT_TEST_SUITE_END();
  

public:
 
  MPIAccessTest():CppUnit::TestFixture(){}
  ~MPIAccessTest(){}  
  void setUp(){}
  void tearDown(){}
  void test_MPI_Access_Send_Recv() ;
  void test_MPI_Access_Cyclic_Send_Recv() ;
  void test_MPI_Access_SendRecv() ;
  void test_MPI_Access_ISend_IRecv() ;
  void test_MPI_Access_Cyclic_ISend_IRecv() ;
  void test_MPI_Access_ISendRecv() ;
  void test_MPI_Access_Probe() ;
  void test_MPI_Access_IProbe() ;
  void test_MPI_Access_Cancel() ;
  void test_MPI_Access_Send_Recv_Length() ;
  void test_MPI_Access_ISend_IRecv_Length() ;
  void test_MPI_Access_ISend_IRecv_Length_1() ;
  void test_MPI_Access_Time() ;
  void test_MPI_Access_Time_0() ;
  void test_MPI_Access_ISend_IRecv_BottleNeck() ;

private:
  };

// to automatically remove temporary files from disk
class MPIAccessTest_TmpFilesRemover
{
public:
  MPIAccessTest_TmpFilesRemover() {}
  ~MPIAccessTest_TmpFilesRemover();
  bool Register(const std::string theTmpFile);

private:
  std::set<std::string> myTmpFiles;
};

/*!
 *  Tool to print array to stream.
 */
template<class T>
void MPIAccessTest_DumpArray (std::ostream & stream, const T* array, const int length, const std::string text)
{
  stream << text << ": {";
  if (length > 0) {
    stream << array[0];
    for (int i = 1; i < length; i++) {
      stream << ", " << array[i];
    }
  }
  stream << "}" << std::endl;
}

#endif
