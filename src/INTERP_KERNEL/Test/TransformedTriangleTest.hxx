#ifndef __TU_TRANSFORMED_TRIANGLE_HXX__
#define __TU_TRANSFORMED_TRIANGLE_HXX__

#include <cppunit/extensions/HelperMacros.h>
#include "../TransformedTriangle.hxx"

#define ERR_TOL 1.0e-8

using INTERP_KERNEL::TransformedTriangle;

namespace INTERP_TEST
{

  /**
   * \brief Test suite testing some of the low level methods of TransformedTriangle.
   *
   */
  class TransformedTriangleTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( TransformedTriangleTest );
    CPPUNIT_TEST( test_constructor );
    CPPUNIT_TEST( test_calcUnstableC );
    CPPUNIT_TEST( test_calcUnstableT );
    //removed because the test fails to enter the desired code branch
 //   CPPUNIT_TEST( test_calcStableC_Consistency );
    CPPUNIT_TEST_SUITE_END();

    typedef INTERP_KERNEL::TransformedTriangle::TriSegment TriSegment;
    typedef INTERP_KERNEL::TransformedTriangle::DoubleProduct DoubleProduct;

  public:
    void setUp();

    void tearDown();

    // tests
    void test_constructor();

    void test_calcUnstableC(); 

    void test_calcUnstableT();

    void test_calcStableC_Consistency();

    double p1[3], q1[3], r1[3];
    double hp1, hq1, hr1;
    double Hp1, Hq1, Hr1;

    double p2[3], q2[3], r2[3];
    double hp2, hq2, hr2;
    double Hp2, Hq2, Hr2;

    double stable_c2[24];
  
  private:
    TransformedTriangle* tri1;
    TransformedTriangle* tri2;

  };




}



#endif
