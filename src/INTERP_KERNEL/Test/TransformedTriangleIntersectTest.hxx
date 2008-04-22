#ifndef __TU_TRANSFORMED_TRIANGLE_INTERSECT_HXX__
#define __TU_TRANSFORMED_TRIANGLE_INTERSECT_HXX__

#include <cppunit/extensions/HelperMacros.h>
#include "../TransformedTriangle.hxx"

namespace INTERP_TEST
{

  class TransformedTriangleIntersectTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( TransformedTriangleIntersectTest );

    CPPUNIT_TEST( testTriangle1 );
    CPPUNIT_TEST( testTriangle2 );
    CPPUNIT_TEST( testTriangle3 );
    CPPUNIT_TEST( testTriangle4 );
    CPPUNIT_TEST( testTriangle5 );
    CPPUNIT_TEST( testTriangle6 );
    CPPUNIT_TEST( testTriangle7 );
    CPPUNIT_TEST( testTriangle8 );
    CPPUNIT_TEST( testTriangle9 );
    CPPUNIT_TEST( testTriangle10 );
    CPPUNIT_TEST( testTriangle11 );
    CPPUNIT_TEST( testTriangle12 );
    CPPUNIT_TEST( testTriangle13 );

    CPPUNIT_TEST_SUITE_END();

    typedef INTERP_KERNEL::TransformedTriangle::TriSegment TriSegment;
    typedef INTERP_KERNEL::TransformedTriangle::DoubleProduct DoubleProduct;

  public:

    void testTriangle1();
  
    void testTriangle2();

    void testTriangle3();

    void testTriangle4();

    void testTriangle5();

    void testTriangle6();

    void testTriangle7();

    void testTriangle8();

    void testTriangle9();
  
    void testTriangle10();
  
    void testTriangle11();
  
    void testTriangle12();

    void testTriangle13();

  private:
 
  };

}






#endif
