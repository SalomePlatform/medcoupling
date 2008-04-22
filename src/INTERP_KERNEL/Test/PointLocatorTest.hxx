#ifndef __TU_POINTLOCATOR_HXX__
#define __TU_POINTLOCATOR_HXX__

#include <cppunit/extensions/HelperMacros.h>
#include "../PointLocator.hxx"

namespace INTERP_TEST
{

  /**
   * \brief Test suite testing some of the low level methods of TransformedTriangle.
   *
   */
  class PointLocatorTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( PointLocatorTest );
    CPPUNIT_TEST( test_PointLocator );
    CPPUNIT_TEST_SUITE_END();

   
  public:
    void setUp();

    void tearDown();

    // tests
    void test_PointLocator();

  };




}



#endif
