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

#ifndef __MEDCOUPLINGBASICSTEST1_HXX__
#define __MEDCOUPLINGBASICSTEST1_HXX__

#include "MEDCouplingBasicsTest.hxx"

#include <map>
#include <vector>

namespace MEDCoupling
{
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingMultiFields;

  class MEDCouplingBasicsTest1 : public MEDCouplingBasicsTest
  {
    CPPUNIT_TEST_SUITE(MEDCouplingBasicsTest1);
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
    CPPUNIT_TEST( testCMesh1 );
    CPPUNIT_TEST( testCMesh2 );
    CPPUNIT_TEST( testScale );
    CPPUNIT_TEST( testTryToShareSameCoords );
    CPPUNIT_TEST( testFindNodeOnPlane );
    CPPUNIT_TEST( testRenumberCells );
    CPPUNIT_TEST( testChangeSpaceDimension );
    CPPUNIT_TEST( testSetConnectivity );
    CPPUNIT_TEST_SUITE_END();
  public:
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
    void testCMesh1();
    void testCMesh2();
    void testScale();
    void testTryToShareSameCoords();
    void testFindNodeOnPlane();
    void testRenumberCells();
    void testChangeSpaceDimension();
    void testSetConnectivity();
  };
}

#endif
