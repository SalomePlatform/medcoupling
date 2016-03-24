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

#ifndef __MEDCOUPLINGBASICSTEST2_HXX__
#define __MEDCOUPLINGBASICSTEST2_HXX__

#include "MEDCouplingBasicsTest.hxx"

#include <map>
#include <vector>

namespace MEDCoupling
{
  class DataArrayDouble;
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingMultiFields;

  class MEDCouplingBasicsTest2 : public MEDCouplingBasicsTest
  {
    CPPUNIT_TEST_SUITE(MEDCouplingBasicsTest2);
    CPPUNIT_TEST( testGaussPointField1 );
    CPPUNIT_TEST( testGaussPointNEField1 );
    CPPUNIT_TEST( testCellOrientation1 );
    CPPUNIT_TEST( testCellOrientation2 );
    CPPUNIT_TEST( testCellOrientation3 );
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
    CPPUNIT_TEST_SUITE_END();
  public:
    void testGaussPointField1();
    void testGaussPointNEField1();
    void testCellOrientation1();
    void testCellOrientation2();
    void testCellOrientation3();
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
  };
}

#endif
