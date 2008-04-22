#include "PointLocatorTest.hxx"
#include "PointLocator.hxx"
#include "MEDMeshMaker.hxx"

#include <iostream>
#include <list>

namespace INTERP_TEST
{


  void PointLocatorTest::setUp() 
  {
  }

 
  void PointLocatorTest::tearDown() 
  {
  }

  /**
   * Test that creates a tree in 2D and check that 
   * the results are correct in three
   * cases :
   * a non matching search
   * a standard case
   * a bbox overlapping the bboxes of the tree
   */
  void PointLocatorTest::test_PointLocator() {
		MEDMEM::MESH* mesh2D= MEDMeshMaker(2,2,MED_EN::MED_QUAD4);
		INTERP_KERNEL::PointLocator pl(*mesh2D) ;
		double x[2]={0.0,0.0};
		std::list<int> elems = pl.locate(x);
		CPPUNIT_ASSERT_EQUAL(1,(int)elems.size());
		CPPUNIT_ASSERT_EQUAL(1,(int)(*(elems.begin())));
		elems.clear();
		
		double x2[2]={0.25,0.25};
		elems = pl.locate(x2);
		CPPUNIT_ASSERT_EQUAL(1,(int)elems.size());
		CPPUNIT_ASSERT_EQUAL(1,(int)(*(elems.begin())));
		elems.clear();
		
		double x3[2]={0.5,0.5};
		elems = pl.locate(x3);
		CPPUNIT_ASSERT_EQUAL(4,(int)elems.size());
		elems.clear();

		double x4[2]={-1.0,0.0};
		elems = pl.locate(x4);
		CPPUNIT_ASSERT_EQUAL(0,(int)elems.size());
		elems.clear();

		MEDMEM::MESH* mesh3D= MEDMeshMaker(3,2,MED_EN::MED_HEXA8);
		INTERP_KERNEL::PointLocator pl3(*mesh3D);
		double xx[3]={0.0,0.0,0.0};
	  elems = pl3.locate(xx);
		CPPUNIT_ASSERT_EQUAL(1,(int)elems.size());
		CPPUNIT_ASSERT_EQUAL(1,(int)*(elems.begin()));
		elems.clear();
		
		double xx2[3]={0.25,0.25,0.25};
		elems = pl3.locate(xx2);
		CPPUNIT_ASSERT_EQUAL(1,(int)elems.size());
		CPPUNIT_ASSERT_EQUAL(1,(int)*(elems.begin()));
			elems.clear();

		double xx3[3]={0.5,0.5,0.5};
		elems = pl3.locate(xx3);
		CPPUNIT_ASSERT_EQUAL(8,(int)elems.size());
		elems.clear();
		
		double xx4[3]={-1.0,0.0,0.0};
		elems = pl3.locate(x4);
		CPPUNIT_ASSERT_EQUAL(0,(int)elems.size());
		elems.clear();


  }


}
