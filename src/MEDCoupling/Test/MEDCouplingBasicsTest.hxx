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
    CPPUNIT_TEST( testArray );
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
    CPPUNIT_TEST( testEqualMesh );
    CPPUNIT_TEST( testEqualFieldDouble );
    CPPUNIT_TEST( testNatureChecking );
    CPPUNIT_TEST( testBuildSubMeshData );
    CPPUNIT_TEST( testExtrudedMesh1 );
    CPPUNIT_TEST( testFindCommonNodes );
    CPPUNIT_TEST( testCheckButterflyCells );
    CPPUNIT_TEST( testMergeMesh1 );
    CPPUNIT_TEST( testMergeField1 );
    CPPUNIT_TEST( testFillFromAnalytic );
    CPPUNIT_TEST( testApplyFunc );
    CPPUNIT_TEST( testOperationsOnFields );
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
    CPPUNIT_TEST_SUITE_END();
  public:
    void testArray();
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
    void testEqualMesh();
    void testEqualFieldDouble();
    void testNatureChecking();
    void testBuildSubMeshData();
    void testExtrudedMesh1();
    void testFindCommonNodes();
    void testCheckButterflyCells();
    void testMergeMesh1();
    void testMergeField1();
    void testFillFromAnalytic();
    void testApplyFunc();
    void testOperationsOnFields();
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
  private:
    MEDCouplingUMesh *build3DSourceMesh_2();
    MEDCouplingUMesh *build3DTargetMesh_2();
    MEDCouplingUMesh *build1DTargetMesh_1();
    MEDCouplingUMesh *build2DSourceMesh_1();
    MEDCouplingUMesh *build2DTargetMesh_1();
    MEDCouplingUMesh *build2DTargetMeshPerm_1();
    MEDCouplingUMesh *build2DTargetMesh_2();
    MEDCouplingUMesh *buildCU1DMesh_U();
    MEDCouplingUMesh *buildCU2DMesh_U();
    MEDCouplingUMesh *buildCU3DMesh_U();
    MEDCouplingUMesh *build3DSurfSourceMesh_1();
    MEDCouplingUMesh *build3DSurfSourceMesh_2();
    MEDCouplingUMesh *build3DSurfTargetMesh_1();
    MEDCouplingUMesh *build3DSurfTargetMeshPerm_1();
    MEDCouplingUMesh *build3DSurfTargetMesh_2();
    MEDCouplingUMesh *build3DSourceMesh_1();
    MEDCouplingUMesh *build3DTargetMesh_1();
    MEDCouplingUMesh *build2DTargetMeshMergeNode_1();
    MEDCouplingUMesh *build3DTargetMeshMergeNode_1();
    MEDCouplingUMesh *build3DExtrudedUMesh_1(MEDCouplingUMesh *&mesh2D);
    MEDCouplingUMesh *build2DTargetMeshMerged_1();
    double sumAll(const std::vector< std::map<int,double> >& matrix);
  };
}

#endif
