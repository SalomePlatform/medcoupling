#ifndef __TU_BB_TREE_HXX__
#define __TU_BB_TREE_HXX__

#include <cppunit/extensions/HelperMacros.h>
#include "../BBTree.txx"

namespace INTERP_TEST
{

  /**
   * \brief Test suite testing some of the low level methods of TransformedTriangle.
   *
   */
  class BBTreeTest : public CppUnit::TestFixture
  {

    CPPUNIT_TEST_SUITE( BBTreeTest );
    CPPUNIT_TEST( test_BBTree );
    CPPUNIT_TEST_SUITE_END();

   
  public:
    void setUp();

    void tearDown();

    // tests
    void test_BBTree();

  };




}



#endif
