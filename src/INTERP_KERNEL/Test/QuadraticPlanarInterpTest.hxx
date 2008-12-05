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
#ifndef _QUADRATICPLANARINTERPTEST_HXX_
#define _QUADRATICPLANARINTERPTEST_HXX_

#include <cppunit/extensions/HelperMacros.h>

namespace INTERP_KERNEL
{
  class Node;
  class EdgeArcCircle;
  class QuadraticPolygon;

  class QuadraticPlanarInterpTest : public CppUnit::TestFixture
  {
    CPPUNIT_TEST_SUITE( QuadraticPlanarInterpTest );
    CPPUNIT_TEST( ReadWriteInXfigElementary );
    CPPUNIT_TEST( ReadWriteInXfigGlobal );
    CPPUNIT_TEST( BasicGeometricTools );
    CPPUNIT_TEST( IntersectionBasics );
    CPPUNIT_TEST( EdgeLinUnitary );
    CPPUNIT_TEST( IntersectionEdgeOverlapUnitarySegSeg );
    CPPUNIT_TEST( IntersectionPointOnlyUnitarySegSeg );
    CPPUNIT_TEST( IntersectArcCircleBase );
    CPPUNIT_TEST( IntersectArcCircleFull );
    CPPUNIT_TEST( IntersectArcCircleSegumentBase );
    CPPUNIT_TEST( checkInOutDetection );
    CPPUNIT_TEST( checkAssemblingBases1 );
    CPPUNIT_TEST( checkAssemblingBases2 );
    CPPUNIT_TEST( checkPolygonsIntersection1 );
    CPPUNIT_TEST( checkAreasCalculations );
    CPPUNIT_TEST( checkHighLevelFunctionTest1 );
    CPPUNIT_TEST( check1DInterpLin );
    CPPUNIT_TEST( checkNonRegression1 );
    CPPUNIT_TEST( checkNonRegression2 );
    CPPUNIT_TEST( checkNonRegression3 );
    CPPUNIT_TEST( checkNonRegression4 );
//     CPPUNIT_TEST( checkNonRegression5 );
    CPPUNIT_TEST( checkNonRegression6 );
    CPPUNIT_TEST( checkNonRegression7 );
    CPPUNIT_TEST( checkNonRegression8 );
    CPPUNIT_TEST( checkNonRegression9 );
    CPPUNIT_TEST( checkNonRegression10 );
    CPPUNIT_TEST( checkNonRegression11 );
//     CPPUNIT_TEST( checkNonRegression12 );
    CPPUNIT_TEST_SUITE_END();
  public:  
    void setUp();
    void tearDown();
    void cleanUp();
    //
    void ReadWriteInXfigElementary();
    void ReadWriteInXfigGlobal();
    void BasicGeometricTools();
    void IntersectionBasics();
    void EdgeLinUnitary();
    void IntersectionEdgeOverlapUnitarySegSeg();
    void IntersectionPointOnlyUnitarySegSeg();
    //
    void IntersectArcCircleBase();
    void IntersectArcCircleFull();
    void IntersectArcCircleSegumentBase();
    //
    void checkInOutDetection();
    //
    void checkAssemblingBases1();
    void checkAssemblingBases2();
    //
    void checkPolygonsIntersection1();
    void checkAreasCalculations();
    //
    void checkHighLevelFunctionTest1();
    //
    void check1DInterpLin();
    //
    void checkNonRegression1();
    void checkNonRegression2();
    void checkNonRegression3();
    void checkNonRegression4();
    void checkNonRegression5();
    void checkNonRegression6();
    void checkNonRegression7();
    void checkNonRegression8();
    void checkNonRegression9();
    void checkNonRegression10();
    void checkNonRegression11();
    void checkNonRegression12();
		//
		void checkNormalize();
  private:
    EdgeArcCircle *buildArcOfCircle(const double *center, double radius, double alphaStart, double alphaEnd);
    double btw2NodesAndACenter(const Node& n1, const Node& n2, const double *center);
    void checkBasicsOfPolygons(QuadraticPolygon& pol1, QuadraticPolygon& pol2, bool checkDirection);
  };
}

#endif
