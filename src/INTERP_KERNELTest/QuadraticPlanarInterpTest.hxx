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

#ifndef _QUADRATICPLANARINTERPTEST_HXX_
#define _QUADRATICPLANARINTERPTEST_HXX_

#include <cppunit/extensions/HelperMacros.h>

#include "InterpKernelTestExport.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DQuadraticPolygon.hxx"

namespace INTERP_TEST
{
  class INTERPKERNELTEST_EXPORT QuadraticPlanarInterpTest : public CppUnit::TestFixture
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
    CPPUNIT_TEST( checkPolygonsIntersection2 );
    CPPUNIT_TEST( checkAreasCalculations );
    CPPUNIT_TEST( checkBarycenterCalculations );
    CPPUNIT_TEST( checkHighLevelFunctionTest1 );
    CPPUNIT_TEST( check1DInterpLin );
    CPPUNIT_TEST( checkEpsilonCoherency1 );
    CPPUNIT_TEST( checkNonRegression1 );
    CPPUNIT_TEST( checkNonRegression2 );
    CPPUNIT_TEST( checkNonRegression3 );
    CPPUNIT_TEST( checkNonRegression4 );
    CPPUNIT_TEST( checkNonRegression5 );
    CPPUNIT_TEST( checkNonRegression6 );
    CPPUNIT_TEST( checkNonRegression7 );
    CPPUNIT_TEST( checkNonRegression8 );
    CPPUNIT_TEST( checkNonRegression9 );
    CPPUNIT_TEST( checkNonRegression10 );
    CPPUNIT_TEST( checkNonRegression11 );
    CPPUNIT_TEST( checkNonRegression12 );
    CPPUNIT_TEST ( checkNonRegression13 );
    CPPUNIT_TEST ( checkNonRegression14 );
    CPPUNIT_TEST ( checkNonRegression15 );
    CPPUNIT_TEST ( checkNonRegression16 );
    CPPUNIT_TEST ( checkNonRegression17 );
    //
    CPPUNIT_TEST ( checkNonRegressionOmar0000 );
    CPPUNIT_TEST ( checkNonRegressionOmar0001 );
    CPPUNIT_TEST ( checkNonRegressionOmar0002 );
    CPPUNIT_TEST ( checkNonRegressionOmar0003 );
    CPPUNIT_TEST ( checkNonRegressionOmar0004 );
    CPPUNIT_TEST ( checkNonRegressionOmar0005 );
    CPPUNIT_TEST ( checkNonRegressionOmar0006 );
    CPPUNIT_TEST ( checkNonRegressionOmar0007 );
    CPPUNIT_TEST ( checkNonRegressionOmar0008 );
    CPPUNIT_TEST ( checkNonRegressionOmar0009 );
    CPPUNIT_TEST ( checkNonRegressionOmar0010 );
    CPPUNIT_TEST ( checkNonRegressionOmar0011 );
    CPPUNIT_TEST ( checkNonRegressionOmar2511 );
    CPPUNIT_TEST ( checkNonRegressionOmar0012 );
    CPPUNIT_TEST ( checkNonRegressionOmar0013 );
    CPPUNIT_TEST ( checkNonRegressionOmar0014 );
    CPPUNIT_TEST ( checkNonRegressionOmar0015 );
    CPPUNIT_TEST ( checkNonRegressionOmar0016 );
    CPPUNIT_TEST ( checkNonRegressionOmar0017 );
    CPPUNIT_TEST ( checkNonRegressionOmar0018 );
    CPPUNIT_TEST ( checkNonRegressionOmar0019 );
    CPPUNIT_TEST ( checkNonRegressionOmar0020 );
    CPPUNIT_TEST ( checkNonRegressionOmar0021 );
    CPPUNIT_TEST ( checkNonRegressionOmar0022 );
    CPPUNIT_TEST ( checkNonRegressionOmar0023 );
    CPPUNIT_TEST ( checkNonRegressionOmar0024 );
    CPPUNIT_TEST ( checkNonRegressionOmar2524 );
    CPPUNIT_TEST ( checkNonRegressionOmar0025 );
    CPPUNIT_TEST ( checkNonRegressionOmar0026 );
    CPPUNIT_TEST ( checkNonRegressionOmar0027 );
    CPPUNIT_TEST ( checkNonRegressionOmar0028 );
    CPPUNIT_TEST ( checkNonRegressionOmar0029 );
    CPPUNIT_TEST ( checkNonRegressionOmar0030 );
    //
    CPPUNIT_TEST( checkNormalize );
    CPPUNIT_TEST( checkMakePartitionAbs1 );
    //
    CPPUNIT_TEST( checkIsInOrOut );
    CPPUNIT_TEST( checkGetMiddleOfPoints );
    CPPUNIT_TEST( checkGetMiddleOfPointsOriented );
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
    void checkPolygonsIntersection2();
    void checkAreasCalculations();
    void checkBarycenterCalculations();
    //
    void checkHighLevelFunctionTest1();
    //
    void check1DInterpLin();
    //
    void checkEpsilonCoherency1();
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
    void checkNonRegression13();
    void checkNonRegression14();
    void checkNonRegression15();
    void checkNonRegression16();
    void checkNonRegression17();
    //
    void checkNonRegressionOmar0000();
    void checkNonRegressionOmar0001();
    void checkNonRegressionOmar0002();
    void checkNonRegressionOmar0003();
    void checkNonRegressionOmar0004();
    void checkNonRegressionOmar0005();
    void checkNonRegressionOmar0006();
    void checkNonRegressionOmar0007();
    void checkNonRegressionOmar0008();
    void checkNonRegressionOmar0009();
    void checkNonRegressionOmar0010();
    void checkNonRegressionOmar0011();
    void checkNonRegressionOmar2511();
    void checkNonRegressionOmar0012();
    void checkNonRegressionOmar0013();
    void checkNonRegressionOmar0014();
    void checkNonRegressionOmar0015();
    void checkNonRegressionOmar0016();
    void checkNonRegressionOmar0017();
    void checkNonRegressionOmar0018();
    void checkNonRegressionOmar0019();
    void checkNonRegressionOmar0020();
    void checkNonRegressionOmar0021();
    void checkNonRegressionOmar0022();
    void checkNonRegressionOmar0023();
    void checkNonRegressionOmar0024();
    void checkNonRegressionOmar2524();
    void checkNonRegressionOmar0025();
    void checkNonRegressionOmar0026();
    void checkNonRegressionOmar0027();
    void checkNonRegressionOmar0028();
    void checkNonRegressionOmar0029();
    void checkNonRegressionOmar0030();
    //
    void checkNormalize();
    void checkMakePartitionAbs1();
    // From Adrien:
    void checkIsInOrOut();
    void checkGetMiddleOfPoints();
    void checkGetMiddleOfPointsOriented();

  private:
    INTERP_KERNEL::QuadraticPolygon *buildQuadraticPolygonCoarseInfo(const double *coords, const int *conn, int lgth);
    INTERP_KERNEL::EdgeArcCircle *buildArcOfCircle(const double *center, double radius, double alphaStart, double alphaEnd);
    double btw2NodesAndACenter(const INTERP_KERNEL::Node& n1, const INTERP_KERNEL::Node& n2, const double *center);
    void checkBasicsOfPolygons(INTERP_KERNEL::QuadraticPolygon& pol1, INTERP_KERNEL::QuadraticPolygon& pol2, bool checkDirection);
  };
}

#endif
