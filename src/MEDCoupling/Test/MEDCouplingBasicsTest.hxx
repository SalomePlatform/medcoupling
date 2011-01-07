//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#ifndef __MEDCOUPLINGBASICSTEST_HXX__
#define __MEDCOUPLINGBASICSTEST_HXX__

#include <cppunit/extensions/HelperMacros.h>

#include <map>
#include <vector>

namespace ParaMEDMEM
{
  class MEDCouplingUMesh;

  class MEDCouplingBasicsTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE(MEDCouplingBasicsTest);
    //MEDCouplingBasicsTest1.cxx
    CPPUNIT_TEST( testArray );
    CPPUNIT_TEST( testArray2 );
    CPPUNIT_TEST( testArray3 );
    CPPUNIT_TEST( testMesh );
    CPPUNIT_TEST( testMeshPointsCloud );
    CPPUNIT_TEST( testMeshM1D );
    CPPUNIT_TEST( testDeepCopy );
    CPPUNIT_TEST( testRevNodal );
    CPPUNIT_TEST( testConvertToPolyTypes );
    CPPUNIT_TEST( testDescConn2D );
    CPPUNIT_TEST( testDescConn3D );
    CPPUNIT_TEST( testFindBoundaryNodes );
    CPPUNIT_TEST( testBoundaryMesh );
    CPPUNIT_TEST( testBuildPartOfMySelf );
    CPPUNIT_TEST( testBuildPartOfMySelfNode );
    CPPUNIT_TEST( testZipCoords );
    CPPUNIT_TEST( testZipConnectivity );
    CPPUNIT_TEST( testEqualMesh );
    CPPUNIT_TEST( testEqualFieldDouble );
    CPPUNIT_TEST( testNatureChecking );
    CPPUNIT_TEST( testBuildSubMeshData );
    CPPUNIT_TEST( testExtrudedMesh1 );
    CPPUNIT_TEST( testExtrudedMesh2 );
    CPPUNIT_TEST( testExtrudedMesh3 );
    CPPUNIT_TEST( testExtrudedMesh4 );
    CPPUNIT_TEST( testFindCommonNodes );
    CPPUNIT_TEST( testCheckButterflyCells );
    CPPUNIT_TEST( testMergeMesh1 );
    CPPUNIT_TEST( testMergeMeshOnSameCoords1 );
    CPPUNIT_TEST( testMergeField1 );
    CPPUNIT_TEST( testFillFromAnalytic );
    CPPUNIT_TEST( testFillFromAnalytic2 );
    CPPUNIT_TEST( testApplyFunc );
    CPPUNIT_TEST( testApplyFunc2 );
    CPPUNIT_TEST( testOperationsOnFields );
    CPPUNIT_TEST( testOperationsOnFields2 );
    CPPUNIT_TEST( testOperationsOnFields3 );
    CPPUNIT_TEST( testOperationsOnFields4 );
    CPPUNIT_TEST( testMergeNodesOnField );
    CPPUNIT_TEST( testCheckConsecutiveCellTypes );
    CPPUNIT_TEST( testRearrange2ConsecutiveCellTypes );
    CPPUNIT_TEST( testSplitByType );
    CPPUNIT_TEST( testFuseUMeshesOnSameCoords );
    CPPUNIT_TEST( testFuseUMeshesOnSameCoords2 );
    CPPUNIT_TEST( testBuildOrthogonalField );
    CPPUNIT_TEST( testGetCellsContainingPoint );
    CPPUNIT_TEST( testGetValueOn1 );
    CPPUNIT_TEST( testCMesh0 );
    CPPUNIT_TEST( testScale );
    CPPUNIT_TEST( testTryToShareSameCoords );
    CPPUNIT_TEST( testFindNodeOnPlane );
    CPPUNIT_TEST( testRenumberCells );
    CPPUNIT_TEST( testChangeSpaceDimension );
    //MEDCouplingBasicsTest2.cxx
    CPPUNIT_TEST( testGaussPointField1 );
    CPPUNIT_TEST( testGaussPointNEField1 );
    CPPUNIT_TEST( testCellOrientation1 );
    CPPUNIT_TEST( testCellOrientation2 );
    CPPUNIT_TEST( testPolyhedronBarycenter );
    CPPUNIT_TEST( testNormL12Integ1D );
    CPPUNIT_TEST( testAreaBary2D );
    CPPUNIT_TEST( testAreaBary3D );
    CPPUNIT_TEST( testRenumberCellsForFields );
    CPPUNIT_TEST( testRenumberNodesForFields );
    CPPUNIT_TEST( testConvertQuadraticCellsToLinear );
    CPPUNIT_TEST( testCheckGeoEquivalWith );
    CPPUNIT_TEST( testCheckGeoEquivalWith2 );
    CPPUNIT_TEST( testCopyTinyStringsFromOnFields );
    CPPUNIT_TEST( testTryToShareSameCoordsPermute );
    CPPUNIT_TEST( testTryToShareSameCoordsPermute2 );
    CPPUNIT_TEST( testChangeUnderlyingMesh1 );
    CPPUNIT_TEST( testGetMaxValue1 );
    CPPUNIT_TEST( testSubstractInPlaceDM1 );
    CPPUNIT_TEST( testDotCrossProduct1 );
    CPPUNIT_TEST( testMinMaxFields1 );
    CPPUNIT_TEST( testApplyLin1 );
    CPPUNIT_TEST( testGetIdsInRange1 );
    CPPUNIT_TEST( testBuildSubPart1 );
    CPPUNIT_TEST( testDoublyContractedProduct1 );
    CPPUNIT_TEST( testDeterminant1 );
    CPPUNIT_TEST( testEigenValues1 );
    CPPUNIT_TEST( testEigenVectors1 );
    CPPUNIT_TEST( testInverse1 );
    CPPUNIT_TEST( testTrace1 );
    CPPUNIT_TEST( testDeviator1 );
    CPPUNIT_TEST( testMagnitude1 );
    CPPUNIT_TEST( testMaxPerTuple1 );
    CPPUNIT_TEST( testChangeNbOfComponents );
    CPPUNIT_TEST( testSortPerTuple1 );
    CPPUNIT_TEST( testIsEqualWithoutConsideringStr1 );
    CPPUNIT_TEST( testGetNodeIdsOfCell1 );
    CPPUNIT_TEST( testGetEdgeRatioField1 );
    CPPUNIT_TEST( testFillFromAnalytic3 );
    CPPUNIT_TEST( testFieldDoubleOpEqual1 );
    CPPUNIT_TEST( testAreaBary3D2 );
    CPPUNIT_TEST( testGetMeasureFieldCMesh1 );
    CPPUNIT_TEST( testFieldDoubleZipCoords1 );
    CPPUNIT_TEST( testFieldDoubleZipConnectivity1 );
    CPPUNIT_TEST( testDaDoubleRenumber1 );
    CPPUNIT_TEST( testDaDoubleRenumberAndReduce1 );
    CPPUNIT_TEST( testDaDoubleRenumberInPlace1 );
    CPPUNIT_TEST( testDaDoubleSelectByTupleId1 );
    CPPUNIT_TEST( testDaDoubleRenumberR1 );
    CPPUNIT_TEST( testDaDoubleRenumberInPlaceR1 );
    CPPUNIT_TEST( testDaDoubleGetMinMaxValues1 );
    CPPUNIT_TEST( testFieldDoubleGetMinMaxValues2 );
    CPPUNIT_TEST( testBuildUnstructuredCMesh1 );
    CPPUNIT_TEST( testDataArrayIntInvertO2NNO21 );
    CPPUNIT_TEST( testKeepSetSelectedComponent1 );
    CPPUNIT_TEST( testKeepSetSelectedComponent2 );
    CPPUNIT_TEST( testDAIGetIdsEqual1 );
    CPPUNIT_TEST( testDAIGetIdsEqualList1 );
    CPPUNIT_TEST( testDAFromNoInterlace1 );
    CPPUNIT_TEST( testDAToNoInterlace1 );
    CPPUNIT_TEST( testDAIsUniform1 );
    CPPUNIT_TEST( testDADFromPolarToCart1 );
    CPPUNIT_TEST( testDADFromCylToCart1 );
    CPPUNIT_TEST( testDADFromSpherToCart1 );
    CPPUNIT_TEST( testUnPolyze1 );
    CPPUNIT_TEST( testConvertDegeneratedCells1 );
    CPPUNIT_TEST( testGetNodeIdsNearPoints1 );
    CPPUNIT_TEST( testFieldCopyTinyAttrFrom1 );
    CPPUNIT_TEST( testExtrudedMesh5 );
    CPPUNIT_TEST( testExtrudedMesh6 );
    CPPUNIT_TEST( testExtrudedMesh7 );
    CPPUNIT_TEST( testSimplexize1 );
    CPPUNIT_TEST( testSimplexize2 );
    CPPUNIT_TEST( testDAMeld1 );
    CPPUNIT_TEST( testFieldMeld1 );
    CPPUNIT_TEST( testMergeNodes2 );
    CPPUNIT_TEST( testMergeField2 );
    CPPUNIT_TEST( testDAIBuildComplement1 );
    CPPUNIT_TEST( testDAIBuildUnion1 );
    CPPUNIT_TEST( testDAIBuildIntersection1 );
    CPPUNIT_TEST( testDAIDeltaShiftIndex1 );
    CPPUNIT_TEST( testDaDoubleSelectByTupleIdSafe1 );
    CPPUNIT_TEST( testAreCellsIncludedIn1 );
    CPPUNIT_TEST( testDAIBuildSubstraction1 );
    CPPUNIT_TEST( testBuildOrthogonalField2 );
    CPPUNIT_TEST( testUMInsertNextCell1 );
    CPPUNIT_TEST( testFieldOperatorDivDiffComp1 );
    CPPUNIT_TEST( testDARearrange1 );
    CPPUNIT_TEST( testGetDifferentValues1 );
    CPPUNIT_TEST( testDAIBuildPermutationArr1 );
    CPPUNIT_TEST( testAreCellsIncludedIn2 );
    CPPUNIT_TEST( testUMeshGetPartBarycenterAndOwner1 );
    CPPUNIT_TEST( testUMeshGetPartMeasureField1 );
    CPPUNIT_TEST( testUMeshBuildPartOrthogonalField1 );
    CPPUNIT_TEST( testUMeshGetTypesOfPart1 );
    CPPUNIT_TEST( testUMeshKeepCellIdsByType1 );
    //MEDCouplingBasicsTestInterp.cxx
    CPPUNIT_TEST( test2DInterpP0P0_1 );
    CPPUNIT_TEST( test2DInterpP0P0PL_1 );
    CPPUNIT_TEST( test2DInterpP0P0PL_2 );
    CPPUNIT_TEST( test2DInterpP0P0PL_3 );
    CPPUNIT_TEST( test2DInterpP0P0PL_4 );
    CPPUNIT_TEST( test2DInterpP0P1_1 );
    CPPUNIT_TEST( test2DInterpP0P1PL_1 );
    CPPUNIT_TEST( test2DInterpP0P1PL_2 );
    CPPUNIT_TEST( test2DInterpP1P0_1 );
    CPPUNIT_TEST( test2DInterpP1P0PL_1 );
    CPPUNIT_TEST( test2DInterpP1P0PL_2 );
    CPPUNIT_TEST( test2DInterpP1P1_1 );
    CPPUNIT_TEST( test2DInterpP1P1PL_1 );
    CPPUNIT_TEST( test3DSurfInterpP0P0_1 );
    CPPUNIT_TEST( test3DSurfInterpP0P0PL_1 );
    CPPUNIT_TEST( test3DSurfInterpP0P1_1 );
    CPPUNIT_TEST( test3DSurfInterpP0P1PL_1 );
    CPPUNIT_TEST( test3DSurfInterpP1P0_1 );
    CPPUNIT_TEST( test3DSurfInterpP1P0PL_1 );
    CPPUNIT_TEST( test3DSurfInterpP1P1_1 );
    CPPUNIT_TEST( test3DSurfInterpP1P1PL_1 );
    CPPUNIT_TEST( test3DSurfInterpP0P0_2 );
    CPPUNIT_TEST( test3DSurfInterpP0P0_3 );

    CPPUNIT_TEST( testInterpolationCC );
    CPPUNIT_TEST( testInterpolationCU1D );
    CPPUNIT_TEST( testInterpolationCU2D );
    CPPUNIT_TEST( testInterpolationCU3D );

    CPPUNIT_TEST( test3DInterpP0P0_1 );
    CPPUNIT_TEST( test3DInterpP0P0PL_1 );
    CPPUNIT_TEST( test3DInterpP0P0PL_2 );
    CPPUNIT_TEST( test3DInterpP0P0PL_3 );
    CPPUNIT_TEST( test3DInterpP0P0PL_4 );
    CPPUNIT_TEST( test3DInterpP0P1_1 );
    CPPUNIT_TEST( test3DInterpP0P1PL_1 );
    CPPUNIT_TEST( test3DInterpP1P0_1 );
    CPPUNIT_TEST( test3DInterpP1P0PL_1 );
    CPPUNIT_TEST( test3DInterpP1P1_1 );
    CPPUNIT_TEST( test3DInterpP1P1PL_1 );
    CPPUNIT_TEST( test3DInterpP0P0Empty );
    CPPUNIT_TEST( test2DInterpP0IntegralUniform );
    CPPUNIT_TEST( test3DSurfInterpP0IntegralUniform );
    CPPUNIT_TEST( test3DInterpP0IntegralUniform );
    CPPUNIT_TEST( test2DInterpP1IntegralUniform );
    CPPUNIT_TEST( test3DInterpP1IntegralUniform );
    CPPUNIT_TEST( test2DInterpP1P0Bary_1 );
    CPPUNIT_TEST( test3DSurfInterpP1P0Bary_1 );
    CPPUNIT_TEST( test3DInterpP1P0Bary_1 );
    CPPUNIT_TEST( test3DTo1DInterpP0P0PL_1 );

    CPPUNIT_TEST( test1DInterp_1 );
    CPPUNIT_TEST( test2DCurveInterpP0P0_1 );
    CPPUNIT_TEST( test2DCurveInterpP0P0_2 );
    CPPUNIT_TEST( test2DCurveInterpP0P1_1 );
    CPPUNIT_TEST( test2DCurveInterpP1P0_1 );
    CPPUNIT_TEST( test2DCurveInterpP1P1_1 );
    CPPUNIT_TEST_SUITE_END();
  public:
    //MEDCouplingBasicsTest1.cxx
    void testArray();
    void testArray2();
    void testArray3();
    void testMesh();
    void testMeshPointsCloud();
    void testMeshM1D();
    void testDeepCopy();
    void testRevNodal();
    void testConvertToPolyTypes();
    void testDescConn2D();
    void testDescConn3D();
    void testFindBoundaryNodes();
    void testBoundaryMesh();
    void testBuildPartOfMySelf();
    void testBuildPartOfMySelfNode();
    void testZipCoords();
    void testZipConnectivity();
    void testEqualMesh();
    void testEqualFieldDouble();
    void testNatureChecking();
    void testBuildSubMeshData();
    void testExtrudedMesh1();
    void testExtrudedMesh2();
    void testExtrudedMesh3();
    void testExtrudedMesh4();
    void testFindCommonNodes();
    void testCheckButterflyCells();
    void testMergeMesh1();
    void testMergeMeshOnSameCoords1();
    void testMergeField1();
    void testFillFromAnalytic();
    void testFillFromAnalytic2();
    void testApplyFunc();
    void testApplyFunc2();
    void testOperationsOnFields();
    void testOperationsOnFields2();
    void testOperationsOnFields3();
    void testOperationsOnFields4();
    void testMergeNodesOnField();
    void testCheckConsecutiveCellTypes();
    void testRearrange2ConsecutiveCellTypes();
    void testSplitByType();
    void testFuseUMeshesOnSameCoords();
    void testFuseUMeshesOnSameCoords2();
    void testBuildOrthogonalField();
    void testGetCellsContainingPoint();
    void testGetValueOn1();
    void testCMesh0();
    void testScale();
    void testTryToShareSameCoords();
    void testFindNodeOnPlane();
    void testRenumberCells();
    void testChangeSpaceDimension();
    //MEDCouplingBasicsTest2.cxx
    void testGaussPointField1();
    void testGaussPointNEField1();
    void testCellOrientation1();
    void testCellOrientation2();
    void testPolyhedronBarycenter();
    void testNormL12Integ1D();
    void testAreaBary2D();
    void testAreaBary3D();
    void testRenumberCellsForFields();
    void testRenumberNodesForFields();
    void testConvertQuadraticCellsToLinear();
    void testCheckGeoEquivalWith();
    void testCheckGeoEquivalWith2();
    void testCopyTinyStringsFromOnFields();
    void testTryToShareSameCoordsPermute();
    void testTryToShareSameCoordsPermute2();
    void testChangeUnderlyingMesh1();
    void testGetMaxValue1();
    void testSubstractInPlaceDM1();
    void testDotCrossProduct1();
    void testMinMaxFields1();
    void testApplyLin1();
    void testGetIdsInRange1();
    void testBuildSubPart1();
    void testDoublyContractedProduct1();
    void testDeterminant1();
    void testEigenValues1();
    void testEigenVectors1();
    void testInverse1();
    void testTrace1();
    void testDeviator1();
    void testMagnitude1();
    void testMaxPerTuple1();
    void testChangeNbOfComponents();
    void testSortPerTuple1();
    void testIsEqualWithoutConsideringStr1();
    void testGetNodeIdsOfCell1();
    void testGetEdgeRatioField1();
    void testFillFromAnalytic3();
    void testFieldDoubleOpEqual1();
    void testAreaBary3D2();
    void testGetMeasureFieldCMesh1();
    void testFieldDoubleZipCoords1();
    void testFieldDoubleZipConnectivity1();
    void testDaDoubleRenumber1();
    void testDaDoubleRenumberAndReduce1();
    void testDaDoubleRenumberInPlace1();
    void testDaDoubleSelectByTupleId1();
    void testDaDoubleRenumberR1();
    void testDaDoubleRenumberInPlaceR1();
    void testDaDoubleGetMinMaxValues1();
    void testFieldDoubleGetMinMaxValues2();
    void testBuildUnstructuredCMesh1();
    void testDataArrayIntInvertO2NNO21();
    void testKeepSetSelectedComponent1();
    void testKeepSetSelectedComponent2();
    void testDAIGetIdsEqual1();
    void testDAIGetIdsEqualList1();
    void testDAFromNoInterlace1();
    void testDAToNoInterlace1();
    void testDAIsUniform1();
    void testDADFromPolarToCart1();
    void testDADFromCylToCart1();
    void testDADFromSpherToCart1();
    void testUnPolyze1();
    void testConvertDegeneratedCells1();
    void testGetNodeIdsNearPoints1();
    void testFieldCopyTinyAttrFrom1();
    void testExtrudedMesh5();
    void testExtrudedMesh6();
    void testExtrudedMesh7();
    void testSimplexize1();
    void testSimplexize2();
    void testDAMeld1();
    void testFieldMeld1();
    void testMergeNodes2();
    void testMergeField2();
    void testDAIBuildComplement1();
    void testDAIBuildUnion1();
    void testDAIBuildIntersection1();
    void testDAIDeltaShiftIndex1();
    void testDaDoubleSelectByTupleIdSafe1();
    void testAreCellsIncludedIn1();
    void testDAIBuildSubstraction1();
    void testBuildOrthogonalField2();
    void testUMInsertNextCell1();
    void testFieldOperatorDivDiffComp1();
    void testDARearrange1();
    void testGetDifferentValues1();
    void testDAIBuildPermutationArr1();
    void testAreCellsIncludedIn2();
    void testUMeshGetPartBarycenterAndOwner1();
    void testUMeshGetPartMeasureField1();
    void testUMeshBuildPartOrthogonalField1();
    void testUMeshGetTypesOfPart1();
    void testUMeshKeepCellIdsByType1();
    //MEDCouplingBasicsTestInterp.cxx
    void test2DInterpP0P0_1();
    void test2DInterpP0P0PL_1();
    void test2DInterpP0P0PL_2();
    void test2DInterpP0P0PL_3();
    void test2DInterpP0P0PL_4();
    void test2DInterpP0P1_1();
    void test2DInterpP0P1PL_1();
    void test2DInterpP0P1PL_2();
    void test2DInterpP1P0_1();
    void test2DInterpP1P0PL_1();
    void test2DInterpP1P0PL_2();
    void test2DInterpP1P1_1();
    void test2DInterpP1P1PL_1();
    void test3DSurfInterpP0P0_1();
    void test3DSurfInterpP0P0PL_1();
    void test3DSurfInterpP0P1_1();
    void test3DSurfInterpP0P1PL_1();
    void test3DSurfInterpP1P0_1();
    void test3DSurfInterpP1P0PL_1();
    void test3DSurfInterpP1P1_1();
    void test3DSurfInterpP1P1PL_1();
    void test3DSurfInterpP0P0_2();
    void test3DSurfInterpP0P0_3();
    void test3DInterpP0P0_1();
    void test3DInterpP0P0PL_1();
    void test3DInterpP0P0PL_2();
    void test3DInterpP0P0PL_3();
    void test3DInterpP0P0PL_4();
    void test3DInterpP0P1_1();
    void test3DInterpP0P1PL_1();
    void test3DInterpP1P0_1();
    void test3DInterpP1P0PL_1();
    void test3DInterpP1P1_1();
    void test3DInterpP1P1PL_1();

    void testInterpolationCC();
    void testInterpolationCU1D();
    void testInterpolationCU2D();
    void testInterpolationCU3D();

    void test3DInterpP0P0Empty();
    void test2DInterpP0IntegralUniform();
    void test3DSurfInterpP0IntegralUniform();
    void test3DInterpP0IntegralUniform();
    void test2DInterpP1IntegralUniform();
    void test3DInterpP1IntegralUniform();
    void test2DInterpP1P0Bary_1();
    void test3DSurfInterpP1P0Bary_1();
    void test3DInterpP1P0Bary_1();
    void test3DTo1DInterpP0P0PL_1();

    void test1DInterp_1();
    void test2DCurveInterpP0P0_1();
    void test2DCurveInterpP0P0_2();
    void test2DCurveInterpP0P1_1();
    void test2DCurveInterpP1P0_1();
    void test2DCurveInterpP1P1_1();

  public:
    static MEDCouplingUMesh *build3DSourceMesh_2();
    static MEDCouplingUMesh *build3DTargetMesh_2();
    static MEDCouplingUMesh *build1DTargetMesh_1();
    static MEDCouplingUMesh *build2DSourceMesh_1();
    static MEDCouplingUMesh *build2DTargetMesh_1();
    static MEDCouplingUMesh *build2DTargetMeshPerm_1();
    static MEDCouplingUMesh *build2DTargetMesh_2();
    static MEDCouplingUMesh *buildCU1DMesh_U();
    static MEDCouplingUMesh *buildCU2DMesh_U();
    static MEDCouplingUMesh *buildCU3DMesh_U();
    static MEDCouplingUMesh *build3DSurfSourceMesh_1();
    static MEDCouplingUMesh *build3DSurfSourceMesh_2();
    static MEDCouplingUMesh *build3DSurfTargetMesh_1();
    static MEDCouplingUMesh *build3DSurfTargetMeshPerm_1();
    static MEDCouplingUMesh *build3DSurfTargetMesh_2();
    static MEDCouplingUMesh *build3DSourceMesh_1();
    static MEDCouplingUMesh *build3DTargetMesh_1();
    static MEDCouplingUMesh *build2DTargetMeshMergeNode_1();
    static MEDCouplingUMesh *build3DTargetMeshMergeNode_1();
    static MEDCouplingUMesh *build3DExtrudedUMesh_1(MEDCouplingUMesh *&mesh2D);
    static void build3DExtrudedUMesh_2(MEDCouplingUMesh *&meshN, MEDCouplingUMesh *&meshTT, MEDCouplingUMesh *&meshTF);
    static MEDCouplingUMesh *build2DTargetMeshMerged_1();
    static MEDCouplingUMesh *build2DCurveMesh(double dx, double dy);
    static MEDCouplingUMesh *build1DMesh(double dx);
    static MEDCouplingUMesh *build1DSourceMesh_2();
    static MEDCouplingUMesh *build1DTargetMesh_2();
    static MEDCouplingUMesh *build2DCurveSourceMesh_2();
    static MEDCouplingUMesh *build2DCurveTargetMesh_2();
    static MEDCouplingUMesh *build1DTargetMesh_3();
    static MEDCouplingUMesh *build2DCurveTargetMesh_3();
    static MEDCouplingUMesh *build2DTargetMesh_3();
    static MEDCouplingUMesh *build3DTargetMesh_3();
    static MEDCouplingUMesh *build2DTargetMesh_4();
    static double sumAll(const std::vector< std::map<int,double> >& matrix);
  };
}

#endif
