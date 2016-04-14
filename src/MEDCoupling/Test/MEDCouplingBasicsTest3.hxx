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

#ifndef __MEDCOUPLINGBASICSTEST3_HXX__
#define __MEDCOUPLINGBASICSTEST3_HXX__

#include "MEDCouplingBasicsTest.hxx"

#include <map>
#include <vector>

namespace MEDCoupling
{
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingMultiFields;

  class MEDCouplingBasicsTest3 : public MEDCouplingBasicsTest
  {
    CPPUNIT_TEST_SUITE(MEDCouplingBasicsTest3);
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
    CPPUNIT_TEST( testElementaryDAThrowAndSpecialCases );
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
    CPPUNIT_TEST( testDAIAggregateMulti1 );
    CPPUNIT_TEST( testMergeUMeshes2 );
    CPPUNIT_TEST( testBuild0DMeshFromCoords1 );
    CPPUNIT_TEST_SUITE_END();
  public:
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
    void testElementaryDAThrowAndSpecialCases();
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
    void testDAIAggregateMulti1();
    void testMergeUMeshes2();
    void testBuild0DMeshFromCoords1();
  };
}

#endif
