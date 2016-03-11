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

#ifndef __MULTI_ELEMENT_TETRA_TESTS_HXX_
#define __MULTI_ELEMENT_TETRA_TESTS_HXX_

#include "InterpolationTestSuite.hxx"

namespace INTERP_TEST
{
  /**
   * \brief Class testing algorithm by intersecting meshes of several 
   * elements (all tetrahedra) - up to a few thousand. This serves to check the 
   * filtering methods and the matrix assemblage, as well as verifying
   * that computation errors do not become unmanageable. It uses mehes of 
   * different geometries : tetrahedra, boxes and cylinders.
   *
   */
  class MultiElementTetraTests : public InterpolationTestSuite<3,3>
  {
    CPPUNIT_TEST_SUITE( MultiElementTetraTests );

    CPPUNIT_TEST( tetraComplexIncluded );
    CPPUNIT_TEST( dividedUnitTetraSimplerReflexive );
    CPPUNIT_TEST( dividedUnitTetraReflexive );
    CPPUNIT_TEST( nudgedDividedUnitTetraSimpler );
    CPPUNIT_TEST( nudgedDividedUnitTetra );
    CPPUNIT_TEST( dividedGenTetra );
    CPPUNIT_TEST( tinyBoxReflexive );
    CPPUNIT_TEST( moderateBoxEvenSmallerReflexive );
    CPPUNIT_TEST( moderateBoxSmallReflexive );
    CPPUNIT_TEST( boxReflexive );
    CPPUNIT_TEST( boxReflexiveModerate );
    CPPUNIT_TEST( tetraBoxes );
    CPPUNIT_TEST( moderateBoxesSmaller );
    CPPUNIT_TEST( moderateBoxes );

    CPPUNIT_TEST_SUITE_END();

  public:

    /// Tetrahedron situated totally inside another
    /// \brief Status : pass
    void tetraComplexIncluded()
    {
      _testTools->intersectMeshes("ComplexIncludedTetra", "ComplexIncludingTetra", 17.0156);
    }

    /// Unit tetrahedron divided in 4 elements intersecting itself.
    /// \brief Status : pass
    void dividedUnitTetraSimplerReflexive()
    {
      _testTools->intersectMeshes("DividedUnitTetraSimpler", "DividedUnitTetraSimpler", 0.1666667);
    }

    /// Unit tetrahedron divided in 14 elements intersecting itself.
    /// \brief Status : pass
    void dividedUnitTetraReflexive()
    {
      _testTools->intersectMeshes("DividedUnitTetra", "DividedUnitTetra", 0.1666667);
    }

    /// Unit tetrahedron divided in 4 elements intersecting slightly displaced version of itself.
    /// \brief Status : pass
    void nudgedDividedUnitTetraSimpler()
    {
      _testTools->intersectMeshes("NudgedDividedUnitTetraSimpler", "DividedUnitTetraSimpler", 0.150191);
    }

    /// Unit tetrahedron divided in 14 elements intersecting slightly displaced version of itself.
    /// \brief Status : pass
    void nudgedDividedUnitTetra()
    {
      _testTools->intersectMeshes("NudgedDividedUnitTetra", "DividedUnitTetra", 0.150191);
    }

    /// Two intersecting tetrahedra in general position, one with 23 elements, the other with 643 elements
    /// \brief Status : pass
    void dividedGenTetra()
    {
      _testTools->intersectMeshes("DividedGenTetra1",  "DividedGenTetra2", 0.546329);
    }

    /// Large box in general position with 12 elements intersecting itself
    /// \brief Status : pass
    void tinyBoxReflexive()
    {
      _testTools->intersectMeshes("TinyBox", "TinyBox", 979200);
    }

    /// Small box in general position with 33 elements intersecting itself
    /// \brief Status : pass
    void boxReflexive()
    {
      _testTools->intersectMeshes("Box3",  "Box3", 13.9954);
    }

    /// Box in general position with 67 elements intersecting itself
    /// \brief Status : pass
    void moderateBoxEvenSmallerReflexive()
    {
      _testTools->intersectMeshes("BoxEvenSmaller1", "BoxEvenSmaller1", 1.44018e6);
    }

    /// Box in general position with 544 elements intersecting itself
    /// \brief Status : pass
    void moderateBoxSmallReflexive()
    {
      _testTools->intersectMeshes("BoxModSmall1", "BoxModSmall1", 1.44018e6);
    }

    /// Large box in general position with 2943 elements intersecting itself
    /// \brief Status : pass
    void boxReflexiveModerate()
    {
      _testTools->intersectMeshes("Box1Moderate",  "Box1Moderate", 1.0e6);
    }

    /// Two intersecting boxes in general position with 12 and 18 elements
    /// \brief Status : pass
    void tetraBoxes()
    {
      _testTools->intersectMeshes("Box1", "Box2", 124.197);
    }
    
    /// Two intersecting boxes in general position with 430 and 544 elements
    /// \brief Status : pass
    void moderateBoxesSmaller()
    {
      _testTools->intersectMeshes("BoxModSmall1", "BoxModSmall2", 321853);
    }

    /// Two intersecting boxes in general position with 2943 and 3068 elements
    /// \brief Status : pass
    void moderateBoxes()
    {
      _testTools->intersectMeshes("Box1Moderate",  "Box2Moderate", 376856);
    }

  };
}

#endif
