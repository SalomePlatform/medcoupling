// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

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
      //_testTools->_intersectionType=INTERP_KERNEL::Convex;// valgrind complains !
      //_testTools->calcIntersectionMatrix("square1.med", "Mesh_2","square1.med","Mesh_2", m);
    }
  };
}

#endif
