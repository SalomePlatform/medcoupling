//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef _ParaMEDMEMTEST_HXX_
#define _ParaMEDMEMTEST_HXX_

#include <cppunit/extensions/HelperMacros.h>

#include <set>
#include <string>
#include <iostream>
#include "mpi.h"


class ParaMEDMEMTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( ParaMEDMEMTest );
  CPPUNIT_TEST(testMPIProcessorGroup_constructor);
  CPPUNIT_TEST(testMPIProcessorGroup_boolean);
  CPPUNIT_TEST(testMPIProcessorGroup_rank);
  CPPUNIT_TEST(testBlockTopology_constructor);
  CPPUNIT_TEST(testBlockTopology_serialize);
  CPPUNIT_TEST(testIntersectionDEC_2D);
  CPPUNIT_TEST(testIntersectionDEC_2DP0P1);
  CPPUNIT_TEST(testIntersectionDEC_3D);

  CPPUNIT_TEST(testSynchronousEqualIntersectionWithoutInterpNativeDEC_2D);
  CPPUNIT_TEST(testSynchronousEqualIntersectionWithoutInterpDEC_2D);
  CPPUNIT_TEST(testSynchronousEqualIntersectionDEC_2D);
  CPPUNIT_TEST(testSynchronousFasterSourceIntersectionDEC_2D);
  CPPUNIT_TEST(testSynchronousSlowerSourceIntersectionDEC_2D);
  CPPUNIT_TEST(testSynchronousSlowSourceIntersectionDEC_2D);
  CPPUNIT_TEST(testSynchronousFastSourceIntersectionDEC_2D);
  CPPUNIT_TEST(testAsynchronousEqualIntersectionDEC_2D);
  CPPUNIT_TEST(testAsynchronousFasterSourceIntersectionDEC_2D);
  CPPUNIT_TEST(testAsynchronousSlowerSourceIntersectionDEC_2D);
  CPPUNIT_TEST(testAsynchronousSlowSourceIntersectionDEC_2D);
  CPPUNIT_TEST(testAsynchronousFastSourceIntersectionDEC_2D);
#ifdef MED_ENABLE_FVM
  //can be added again after FVM correction for 2D
  //  CPPUNIT_TEST(testNonCoincidentDEC_2D);
  CPPUNIT_TEST(testNonCoincidentDEC_3D); 
#endif
  CPPUNIT_TEST(testStructuredCoincidentDEC);
  CPPUNIT_TEST(testStructuredCoincidentDEC);
  CPPUNIT_TEST(testICocoTrio1);
  CPPUNIT_TEST(testMEDLoaderRead1);
  CPPUNIT_TEST(testMEDLoaderPolygonRead);
  CPPUNIT_TEST(testMEDLoaderPolyhedronRead);
  //CPPUNIT_TEST(testMEDLoaderWrite1);
  //CPPUNIT_TEST(testMEDLoaderPolygonWrite);
  CPPUNIT_TEST_SUITE_END();
  

public:
 
  ParaMEDMEMTest():CppUnit::TestFixture(){}
  ~ParaMEDMEMTest(){}  
  void setUp(){}
  void tearDown(){}
  void testMPIProcessorGroup_constructor();
  void testMPIProcessorGroup_boolean();
  void testMPIProcessorGroup_rank();
  void testBlockTopology_constructor();
  void testBlockTopology_serialize();
  void testIntersectionDEC_2D();
  void testIntersectionDEC_2DP0P1();
  void testIntersectionDEC_3D();
#ifdef MED_ENABLE_FVM
  void testNonCoincidentDEC_2D();
  void testNonCoincidentDEC_3D();
#endif
  void testStructuredCoincidentDEC();
  void testSynchronousEqualIntersectionWithoutInterpNativeDEC_2D();
  void testSynchronousEqualIntersectionWithoutInterpDEC_2D();
  void testSynchronousEqualIntersectionDEC_2D();
  void testSynchronousFasterSourceIntersectionDEC_2D();
  void testSynchronousSlowerSourceIntersectionDEC_2D();
  void testSynchronousSlowSourceIntersectionDEC_2D();
  void testSynchronousFastSourceIntersectionDEC_2D();

  void testAsynchronousEqualIntersectionDEC_2D();
  void testAsynchronousFasterSourceIntersectionDEC_2D();
  void testAsynchronousSlowerSourceIntersectionDEC_2D();
  void testAsynchronousSlowSourceIntersectionDEC_2D();
  void testAsynchronousFastSourceIntersectionDEC_2D();
  //
  void testICocoTrio1();
  //
  void testMEDLoaderRead1();
  void testMEDLoaderPolygonRead();
  void testMEDLoaderPolyhedronRead();
  void testMEDLoaderWrite1();
  void testMEDLoaderPolygonWrite();
private:
  void testNonCoincidentDEC(const std::string& filename1, 
                            const std::string& meshname1, 
                            const std::string& filename2, 
                            const std::string& meshname2,
                            int nbprocsource, double epsilon);
  void testAsynchronousIntersectionDEC_2D(double dtA, double tmaxA, 
                                          double dtB, double tmaxB,
                                          bool WithPointToPoint, bool Asynchronous, bool WithInterp, const char *srcMeth, const char *targetMeth);
  void testIntersectionDEC_2D_(const char *srcMeth, const char *targetMeth);
  void testIntersectionDEC_3D_(const char *srcMeth, const char *targetMeth);
};

// to automatically remove temporary files from disk
class ParaMEDMEMTest_TmpFilesRemover
{
public:
  ParaMEDMEMTest_TmpFilesRemover() {}
  ~ParaMEDMEMTest_TmpFilesRemover();
  bool Register(const std::string theTmpFile);

private:
  std::set<std::string> myTmpFiles;
};

/*!
 *  Tool to print array to stream.
 */
template<class T>
void ParaMEDMEMTest_DumpArray (std::ostream & stream, const T* array, const int length, const std::string text)
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
