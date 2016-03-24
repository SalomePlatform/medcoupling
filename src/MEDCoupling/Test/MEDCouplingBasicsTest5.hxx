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

#ifndef __MEDCOUPLINGBASICSTEST5_HXX__
#define __MEDCOUPLINGBASICSTEST5_HXX__

#include "MEDCouplingBasicsTest.hxx"

#include <map>
#include <vector>

namespace MEDCoupling
{
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingMultiFields;

  class MEDCouplingBasicsTest5 : public MEDCouplingBasicsTest
  {
    CPPUNIT_TEST_SUITE(MEDCouplingBasicsTest5);
    CPPUNIT_TEST( testUMeshTessellate2D1 );
    CPPUNIT_TEST( testUMeshTessellate2DCurve1 );
    CPPUNIT_TEST( testIntersect2DMeshesTmp4 );
    CPPUNIT_TEST( testGetCellIdsCrossingPlane1 );
    CPPUNIT_TEST( testBuildSlice3D1 );
    CPPUNIT_TEST( testBuildSlice3DSurf1 );
    CPPUNIT_TEST( testDataArrayDoubleAdvSetting1 );
    CPPUNIT_TEST( testDataArrayIntAdvSetting1 );
    CPPUNIT_TEST( testBuildDescendingConnec2Of3DMesh1 );
    CPPUNIT_TEST( testAre2DCellsNotCorrectlyOriented1 );
    CPPUNIT_TEST( testDataArrayAbs1 );
    CPPUNIT_TEST( testGetValueOn3 );
    CPPUNIT_TEST( testGetNodeIdsOfCell2 );
    CPPUNIT_TEST( testRenumberNodesInConn1 );
    CPPUNIT_TEST( testComputeNeighborsOfCells1 );
    CPPUNIT_TEST( testCheckButterflyCellsBug1 );
    CPPUNIT_TEST( testDataArrayIntRange1 );
    CPPUNIT_TEST( testDataArrayDoubleGetMinMaxPerComponent1 );
    CPPUNIT_TEST( testDataArrayIntGetHashCode1 );
    CPPUNIT_TEST( testZipConnectivityPol1 );
    CPPUNIT_TEST( testConvexEnvelop2D1 );
    CPPUNIT_TEST( testDataArraySort1 );
    CPPUNIT_TEST( testPartitionBySpreadZone1 );
    CPPUNIT_TEST( testGiveCellsWithType1 );
    CPPUNIT_TEST( testBuildSlice3D2 );
    CPPUNIT_TEST( testComputeTupleIdsToSelectFromCellIds1 );
    CPPUNIT_TEST( testComputeSkin1 );
    CPPUNIT_TEST( testUMeshSetPartOfMySelf2 );
    CPPUNIT_TEST( testUnPolyze3 );
    CPPUNIT_TEST( testKrSpatialDiscretization1 );
    CPPUNIT_TEST( testDuplicateEachTupleNTimes1 );
    CPPUNIT_TEST( testIntersect2DMeshesTmp5 );
    CPPUNIT_TEST( testDAIBuildUnique1 );
    CPPUNIT_TEST( testDAIPartitionByDifferentValues1 );
    CPPUNIT_TEST( testDAICheckMonotonic1 );
    CPPUNIT_TEST( testIntersect2DMeshesTmp6 );
    CPPUNIT_TEST( testIntersect2DMeshesTmp7 );
    CPPUNIT_TEST( testDAIBuildSubstractionOptimized1 );
    CPPUNIT_TEST( testDAIIsStrictlyMonotonic1 );
    CPPUNIT_TEST( testSimplexize3 );
    CPPUNIT_TEST_SUITE_END();
  public:
    void testUMeshTessellate2D1();
    void testUMeshTessellate2DCurve1();
    void testIntersect2DMeshesTmp4();
    void testGetCellIdsCrossingPlane1();
    void testBuildSlice3D1();
    void testBuildSlice3DSurf1();
    void testDataArrayDoubleAdvSetting1();
    void testDataArrayIntAdvSetting1();
    void testBuildDescendingConnec2Of3DMesh1();
    void testAre2DCellsNotCorrectlyOriented1();
    void testDataArrayAbs1();
    void testGetValueOn3();
    void testGetNodeIdsOfCell2();
    void testRenumberNodesInConn1();
    void testComputeNeighborsOfCells1();
    void testCheckButterflyCellsBug1();
    void testDataArrayIntRange1();
    void testDataArrayDoubleGetMinMaxPerComponent1();
    void testDataArrayIntGetHashCode1();
    void testZipConnectivityPol1();
    void testConvexEnvelop2D1();
    void testDataArraySort1();
    void testPartitionBySpreadZone1();
    void testGiveCellsWithType1();
    void testBuildSlice3D2();
    void testComputeTupleIdsToSelectFromCellIds1();
    void testComputeSkin1();
    void testUMeshSetPartOfMySelf2();
    void testUnPolyze3();
    void testKrSpatialDiscretization1();
    void testDuplicateEachTupleNTimes1();
    void testIntersect2DMeshesTmp5();
    void testDAIBuildUnique1();
    void testDAIPartitionByDifferentValues1();
    void testDAICheckMonotonic1();
    void testIntersect2DMeshesTmp6();
    void testIntersect2DMeshesTmp7();
    void testDAIBuildSubstractionOptimized1();
    void testDAIIsStrictlyMonotonic1();
    void testSimplexize3();
  };
}

#endif
