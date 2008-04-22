#ifndef __TU_TEST_CPPUNIT_HXX__
#define __TU_TEST_CPPUNIT_HXX__

#include <cppunit/extensions/HelperMacros.h>

/**
 * \brief Class tested by TestBogusClass : not very useful
 */
class BogusClass {
  friend class TestBogusClass;

public:
  BogusClass(double _x) : x(_x) {;} 
 
  double getX() { return x; }

private: 
  double x;
};
  
/**
 * \brief Class used to figure out CppUnit : not very useful
 *
 */
class TestBogusClass : public CppUnit::TestFixture
{

  CPPUNIT_TEST_SUITE( TestBogusClass );
  CPPUNIT_TEST( test1 );
  CPPUNIT_TEST( test2 );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {
    obj = new BogusClass(3.14);
  }

  void tearDown() {
    delete obj;
  }

  void test1() {
    // test something
    CPPUNIT_ASSERT(obj->x == 3.14);
    CPPUNIT_ASSERT(obj->getX() == obj->x);
  }

  void test2() {
    // test something else
    obj->x += 2.6;
    CPPUNIT_ASSERT(obj->getX() > 3.14);
  }

private:
  BogusClass* obj;

};








#endif
