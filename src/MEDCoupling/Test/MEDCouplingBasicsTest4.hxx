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
// Author : Anthony Geay (CEA/DEN)

#ifndef __MEDCOUPLINGBASICSTEST4_HXX__
#define __MEDCOUPLINGBASICSTEST4_HXX__

#include "MEDCouplingBasicsTest.hxx"

#include <map>
#include <vector>

namespace MEDCoupling
{
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingMultiFields;

  class MEDCouplingBasicsTest4 : public MEDCouplingBasicsTest
  {
    CPPUNIT_TEST_SUITE(MEDCouplingBasicsTest4);
    CPPUNIT_TEST( testDescriptionInMeshTimeUnit1 );
    CPPUNIT_TEST( testMultiFields1 );
    CPPUNIT_TEST( testFieldOverTime1 );
    CPPUNIT_TEST( testDAICheckAndPreparePermutation1 );
    CPPUNIT_TEST( testDAIChangeSurjectiveFormat1 );
    CPPUNIT_TEST( testUMeshGetCellIdsLyingOnNodes1 );
    CPPUNIT_TEST( testUMeshFindCellIdsOnBoundary1 );
    CPPUNIT_TEST( testMeshSetTime1 );
    CPPUNIT_TEST( testApplyFuncTwo1 );
    CPPUNIT_TEST( testApplyFuncThree1 );
    CPPUNIT_TEST( testFillFromAnalyticTwo1 );
    CPPUNIT_TEST( testFillFromAnalyticThree1 );
    CPPUNIT_TEST( testDAUnitVar1 );
    CPPUNIT_TEST( testGaussCoordinates1 );
    CPPUNIT_TEST( testP2Localization1 );
    CPPUNIT_TEST( testP2Localization2 );
    CPPUNIT_TEST( testGetValueOn2 );
    CPPUNIT_TEST( testDAIGetIdsNotEqual1 );
    CPPUNIT_TEST( testDAIComputeOffsets1 );
    CPPUNIT_TEST( testUMeshHexagonPrism1 );
    CPPUNIT_TEST( testDADCheckIsMonotonic );
    CPPUNIT_TEST( testCheckCoherencyDeeper1 );
    CPPUNIT_TEST( testUnPolyze2 );
    CPPUNIT_TEST( testDACpyFrom1 );
    CPPUNIT_TEST( testDAITransformWithIndArr1 );
    CPPUNIT_TEST( testDAIBuildPermArrPerLevel1 );
    CPPUNIT_TEST( testDAIOperations1 );
    CPPUNIT_TEST( testEmulateMEDMEMBDC1 );
    CPPUNIT_TEST( testGetLevArrPerCellTypes1 );
    CPPUNIT_TEST( testSortCellsInMEDFileFrmt1 );
    CPPUNIT_TEST( testBuildPartAndReduceNodes1 );
    CPPUNIT_TEST( testDAITransformWithIndArrR1 );
    CPPUNIT_TEST( testDAISplitByValueRange1 );
    CPPUNIT_TEST( testUMeshSplitProfilePerType1 );
    CPPUNIT_TEST( testDAIBuildExplicitArrByRanges1 );
    CPPUNIT_TEST( testDAIComputeOffsets2 );
    CPPUNIT_TEST( testMergeField3 );
    CPPUNIT_TEST( testGetDistributionOfTypes1 );
    CPPUNIT_TEST( testNorm2_1 );
    CPPUNIT_TEST( testNormMax1 );
    CPPUNIT_TEST( testFindAndCorrectBadOriented3DExtrudedCells1 );
    CPPUNIT_TEST( testConvertExtrudedPolyhedra1 );
    CPPUNIT_TEST( testNonRegressionCopyTinyStrings );
    CPPUNIT_TEST( testDaDSetPartOfValuesAdv1 );
    CPPUNIT_TEST( testUMeshBuildSetInstanceFromThis1 );
    CPPUNIT_TEST( testUMeshMergeMeshesCVW1 );
    CPPUNIT_TEST( testDADFindCommonTuples1 );
    CPPUNIT_TEST( testDABack1 );
    CPPUNIT_TEST( testDADGetDifferentValues1 );
    CPPUNIT_TEST( testDAIBuildOld2NewArrayFromSurjectiveFormat2 );
    CPPUNIT_TEST( testDADIReverse1 );
    CPPUNIT_TEST( testGetNodeIdsInUse1 );
    CPPUNIT_TEST( testBuildDescendingConnec2 );
    CPPUNIT_TEST( testIntersect2DMeshesTmp1 );
    CPPUNIT_TEST( testFindNodesOnLine1 );
    CPPUNIT_TEST( testIntersect2DMeshesTmp2 );
    CPPUNIT_TEST( testBuildPartOfMySelfSafe1 );
    CPPUNIT_TEST( testIntersect2DMeshesTmp3 );
    CPPUNIT_TEST_SUITE_END();
  public:
    void testDescriptionInMeshTimeUnit1();
    void testMultiFields1();
    void testFieldOverTime1();
    void testDAICheckAndPreparePermutation1();
    void testDAIChangeSurjectiveFormat1();
    void testUMeshGetCellIdsLyingOnNodes1();
    void testUMeshFindCellIdsOnBoundary1();
    void testMeshSetTime1();
    void testApplyFuncTwo1();
    void testApplyFuncThree1();
    void testFillFromAnalyticTwo1();
    void testFillFromAnalyticThree1();
    void testDAUnitVar1();
    void testGaussCoordinates1();
    void testQ1Localization1();
    void testP2Localization1();
    void testP2Localization2();
    void testGetValueOn2();
    void testDAIGetIdsNotEqual1();
    void testDAIComputeOffsets1();
    void testUMeshHexagonPrism1();
    void testDADCheckIsMonotonic();
    void testCheckCoherencyDeeper1();
    void testUnPolyze2();
    void testDACpyFrom1();
    void testDAITransformWithIndArr1();
    void testDAIBuildPermArrPerLevel1();
    void testDAIOperations1();
    void testEmulateMEDMEMBDC1();
    void testGetLevArrPerCellTypes1();
    void testSortCellsInMEDFileFrmt1();
    void testBuildPartAndReduceNodes1();
    void testDAITransformWithIndArrR1();
    void testDAISplitByValueRange1();
    void testUMeshSplitProfilePerType1();
    void testDAIBuildExplicitArrByRanges1();
    void testDAIComputeOffsets2();
    void testMergeField3();
    void testGetDistributionOfTypes1();
    void testNorm2_1();
    void testNormMax1();
    void testFindAndCorrectBadOriented3DExtrudedCells1();
    void testConvertExtrudedPolyhedra1();
    void testNonRegressionCopyTinyStrings();
    void testDaDSetPartOfValuesAdv1();
    void testUMeshBuildSetInstanceFromThis1();
    void testUMeshMergeMeshesCVW1();
    void testChangeUnderlyingMeshWithCMesh1();
    void testDADFindCommonTuples1();
    void testDABack1();
    void testDADGetDifferentValues1();
    void testDAIBuildOld2NewArrayFromSurjectiveFormat2();
    void testDADIReverse1();
    void testGetNodeIdsInUse1();
    void testBuildDescendingConnec2();
    void testIntersect2DMeshesTmp1();
    void testFindNodesOnLine1();
    void testIntersect2DMeshesTmp2();
    void testBuildPartOfMySelfSafe1();
    void testIntersect2DMeshesTmp3();
  };
}

#endif
