#ifndef __TU_INTERPOLATIONOPTIONS_HXX__
#define __TU_INTERPOLATIONOPTIONS_HXX__

#include <cppunit/extensions/HelperMacros.h>
#include "../InterpolationOptions.hxx"
#include "MEDMEM_Field.hxx"

namespace INTERP_TEST
{

  /**
   * \brief Test suite testing some of the low level methods of TransformedTriangle.
   *
   */
  class InterpolationOptionsTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( InterpolationOptionsTest );
    CPPUNIT_TEST( test_InterpolationOptions );
    CPPUNIT_TEST_SUITE_END();

   
  public:
    void setUp();

    void tearDown();

    // tests
    void test_InterpolationOptions();

	private:
	 
  };




}



#endif
