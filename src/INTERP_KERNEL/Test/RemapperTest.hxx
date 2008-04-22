#ifndef __TU_REMAPPER_HXX__
#define __TU_REMAPPER_HXX__

#include <cppunit/extensions/HelperMacros.h>
#include "../Remapper.hxx"
#include "MEDMEM_Field.hxx"

namespace INTERP_TEST
{

  /**
   * \brief Test suite testing some of the low level methods of TransformedTriangle.
   *
   */
  class RemapperTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( RemapperTest );
    CPPUNIT_TEST( test_Remapper );
    CPPUNIT_TEST_SUITE_END();

   
  public:
    void setUp();

    void tearDown();

    // tests
    void test_Remapper();

	private:
		void absField(MEDMEM::FIELD<double>&);
  };




}



#endif
