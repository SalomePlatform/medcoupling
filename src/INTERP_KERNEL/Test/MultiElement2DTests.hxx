#ifndef __MULTI_ELEMENT_2D_TESTS_HXX_
#define __MULTI_ELEMENT_2D_TESTS_HXX_

#include "InterpolationTestSuite.hxx"

namespace INTERP_TEST
{
  /**
   * \brief Class testing algorithm by intersecting meshes of several 
   * polygonal elements - up to a few thousand. This serves to check the 
   * filtering methods and the matrix assemblage, as well as verifying
   * that computation errors do not become unmanageable. It uses mehes of 
   * different geometries : triangle, quadrilateral.
   *
   */
  class MultiElement2DTests : public InterpolationTestSuite<2,2>
  {
    CPPUNIT_TEST_SUITE( MultiElement2DTests );
		
		CPPUNIT_TEST(SymetryTranspose2DTest);
		CPPUNIT_TEST(SelfIntersection2DTest);

    CPPUNIT_TEST_SUITE_END();

  public:
		void SymetryTranspose2DTest()
		{ 
			_testTools->_intersectionType=INTERP_KERNEL::Triangulation;
			_testTools->intersectMeshes("square1.med", "Mesh_2","square2.med","Mesh_3", 10000.);
			_testTools->_intersectionType=INTERP_KERNEL::Convex;
			_testTools->intersectMeshes("square1.med", "Mesh_2","square2.med","Mesh_3", 10000.);
		}
		void SelfIntersection2DTest()
		{ 
			IntersectionMatrix m;
			_testTools->_intersectionType=INTERP_KERNEL::Triangulation;
			_testTools->calcIntersectionMatrix("square1.med", "Mesh_2","square1.med","Mesh_2", m);
			_testTools->_intersectionType=INTERP_KERNEL::Convex;
			_testTools->calcIntersectionMatrix("square1.med", "Mesh_2","square1.med","Mesh_2", m);
		}
  };
}

#endif
