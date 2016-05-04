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
  CPPUNIT_TEST(testMPIProcessorGroup_constructor);  // 1 and 2 procs
  CPPUNIT_TEST(testMPIProcessorGroup_boolean);      // 1 and 2 procs
  CPPUNIT_TEST(testMPIProcessorGroup_rank);         // >=2 procs
  CPPUNIT_TEST(testBlockTopology_constructor);      // >=2 procs
  CPPUNIT_TEST(testBlockTopology_serialize);        // 1 proc
  CPPUNIT_TEST(testInterpKernelDEC_1D);             // 5 procs
  CPPUNIT_TEST(testInterpKernelDEC_2DCurve);        // 5 procs
  CPPUNIT_TEST(testInterpKernelDEC_2D);             // 5 procs
  CPPUNIT_TEST(testInterpKernelDEC2_2D);            // 5 procs
  CPPUNIT_TEST(testInterpKernelDEC_2DP0P1);         // Not impl.
  CPPUNIT_TEST(testInterpKernelDEC_3D);             // 3 procs
  CPPUNIT_TEST(testInterpKernelDECNonOverlapp_2D_P0P0);     // 5 procs
  CPPUNIT_TEST(testInterpKernelDECNonOverlapp_2D_P0P1P1P0); // 5 procs
  CPPUNIT_TEST(testInterpKernelDEC2DM1D_P0P0);      // 3 procs
  CPPUNIT_TEST(testInterpKernelDECPartialProcs);    // 3 procs
  CPPUNIT_TEST(testInterpKernelDEC3DSurfEmptyBBox); // 3 procs
  CPPUNIT_TEST(testOverlapDEC1);                    // 3 procs
  CPPUNIT_TEST(testOverlapDEC1_bis);                // 3 procs
  CPPUNIT_TEST(testOverlapDEC1_ter);                // 3 procs
  CPPUNIT_TEST(testOverlapDEC2);                    // 3 procs
  CPPUNIT_TEST(testOverlapDEC2_bis);                // 3 procs
  CPPUNIT_TEST(testOverlapDEC2_ter);                // 3 procs
//  CPPUNIT_TEST(testOverlapDEC3);                    // 2 procs
//  CPPUNIT_TEST(testOverlapDEC4);                    // 2 procs

  CPPUNIT_TEST(testSynchronousEqualInterpKernelWithoutInterpNativeDEC_2D);// 5 procs
  CPPUNIT_TEST(testSynchronousEqualInterpKernelWithoutInterpDEC_2D);      // 5 procs
  CPPUNIT_TEST(testSynchronousEqualInterpKernelDEC_2D);                   // 5 procs
  CPPUNIT_TEST(testSynchronousFasterSourceInterpKernelDEC_2D);            // 5 procs
  CPPUNIT_TEST(testSynchronousSlowerSourceInterpKernelDEC_2D);            // 5 procs
  CPPUNIT_TEST(testSynchronousSlowSourceInterpKernelDEC_2D);              // 5 procs
  CPPUNIT_TEST(testSynchronousFastSourceInterpKernelDEC_2D);              // 5 procs
  CPPUNIT_TEST(testAsynchronousEqualInterpKernelDEC_2D);                  // 5 procs
  CPPUNIT_TEST(testAsynchronousFasterSourceInterpKernelDEC_2D);           // 5 procs
  CPPUNIT_TEST(testAsynchronousSlowerSourceInterpKernelDEC_2D);           // 5 procs
  CPPUNIT_TEST(testAsynchronousSlowSourceInterpKernelDEC_2D);             // 5 procs
  CPPUNIT_TEST(testAsynchronousFastSourceInterpKernelDEC_2D);             // 5 procs
#ifdef MED_ENABLE_FVM
  //can be added again after FVM correction for 2D
  //  CPPUNIT_TEST(testNonCoincidentDEC_2D);
  CPPUNIT_TEST(testNonCoincidentDEC_3D);
#endif
  CPPUNIT_TEST(testStructuredCoincidentDEC);     // 5 procs
  CPPUNIT_TEST(testICoco1);           // 2 procs
  CPPUNIT_TEST(testGauthier1);        // 4 procs
  CPPUNIT_TEST(testGauthier2);        // >= 2 procs
  CPPUNIT_TEST(testGauthier3_1);      // 4 procs
  CPPUNIT_TEST(testGauthier3_2);      // 4 procs
  CPPUNIT_TEST(testGauthier3_3);      // 5 procs
  CPPUNIT_TEST(testGauthier3_4);      // 5 procs
  CPPUNIT_TEST(testGauthier4);        // 3 procs
  CPPUNIT_TEST(testFabienAPI1);       // 3 procs
  CPPUNIT_TEST(testFabienAPI2);       // 3 procs

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
  void testInterpKernelDEC_1D();
  void testInterpKernelDEC_2DCurve();
  void testInterpKernelDEC_2D();
  void testInterpKernelDEC2_2D();
  void testInterpKernelDEC_2DP0P1();
  void testInterpKernelDEC_3D();
  void testInterpKernelDECNonOverlapp_2D_P0P0();
  void testInterpKernelDECNonOverlapp_2D_P0P1P1P0();
  void testInterpKernelDEC2DM1D_P0P0();
  void testInterpKernelDECPartialProcs();
  void testInterpKernelDEC3DSurfEmptyBBox();
  void testOverlapDEC1();
  void testOverlapDEC1_bis();
  void testOverlapDEC1_ter();
  void testOverlapDEC2();
  void testOverlapDEC2_bis();
  void testOverlapDEC2_ter();
  void testOverlapDEC3();
//  void testOverlapDEC3_bis();
  void testOverlapDEC4();
#ifdef MED_ENABLE_FVM
  void testNonCoincidentDEC_2D();
  void testNonCoincidentDEC_3D();
#endif
  void testStructuredCoincidentDEC();
  void testSynchronousEqualInterpKernelWithoutInterpNativeDEC_2D();
  void testSynchronousEqualInterpKernelWithoutInterpDEC_2D();
  void testSynchronousEqualInterpKernelDEC_2D();
  void testSynchronousFasterSourceInterpKernelDEC_2D();
  void testSynchronousSlowerSourceInterpKernelDEC_2D();
  void testSynchronousSlowSourceInterpKernelDEC_2D();
  void testSynchronousFastSourceInterpKernelDEC_2D();

  void testAsynchronousEqualInterpKernelDEC_2D();
  void testAsynchronousFasterSourceInterpKernelDEC_2D();
  void testAsynchronousSlowerSourceInterpKernelDEC_2D();
  void testAsynchronousSlowSourceInterpKernelDEC_2D();
  void testAsynchronousFastSourceInterpKernelDEC_2D();
  //
  void testICoco1();
  void testGauthier1();
  void testGauthier2();
  void testGauthier3_1();
  void testGauthier3_2();
  void testGauthier3_3();
  void testGauthier3_4();
  void testGauthier4();
  void testFabienAPI1();
  void testFabienAPI2();

  std::string getTmpDirectory();
  std::string makeTmpFile( const std::string&, const std::string& = "" );

private:
#ifdef MED_ENABLE_FVM
  void testNonCoincidentDEC(const std::string& filename1, 
                            const std::string& meshname1, 
                            const std::string& filename2, 
                            const std::string& meshname2,
                            int nbprocsource, double epsilon);
#endif
  void testAsynchronousInterpKernelDEC_2D(double dtA, double tmaxA, 
                                          double dtB, double tmaxB,
                                          bool WithPointToPoint, bool Asynchronous, bool WithInterp, const char *srcMeth, const char *targetMeth);
  void testInterpKernelDEC_2D_(const char *srcMeth, const char *targetMeth);
  void testInterpKernelDEC2_2D_(const char *srcMeth, const char *targetMeth);
  void testInterpKernelDEC_3D_(const char *srcMeth, const char *targetMeth);
  void testGauthier3_GEN(bool, int);
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
}

#endif
