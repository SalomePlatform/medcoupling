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
    CPPUNIT_TEST_SUITE_END();
  public:  
    void setUp();
    void tearDown();
    void cleanUp();
    //
    void ReadWriteInXfigElementary();
    void ReadWriteInXfigGlobal();
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
  private:
    EdgeArcCircle *buildArcOfCircle(const double *center, double radius, double alphaStart, double alphaEnd);
    double btw2NodesAndACenter(const Node& n1, const Node& n2, const double *center);
    void checkBasicsOfPolygons(QuadraticPolygon& pol1, QuadraticPolygon& pol2, bool checkDirection);
  };
}

#endif
