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

#ifndef __SINGLE_ELEMENT_TETRA_TESTS_HXX_
#define __SINGLE_ELEMENT_TETRA_TESTS_HXX_

#include "InterpolationTestSuite.hxx"

namespace INTERP_TEST 
{
  /**
   * \brief Class testing algorithm by intersecting simple meshes having only one element each. This serves mainly to verify that
   * the volume calculations between elements is correct.
   *
   */
  class SingleElementTetraTests : public InterpolationTestSuite<3,3>
  {
    CPPUNIT_TEST_SUITE( SingleElementTetraTests );

    CPPUNIT_TEST( tetraReflexiveUnit );
    CPPUNIT_TEST( tetraReflexiveGeneral );
    CPPUNIT_TEST( tetraNudgedSimpler );
    CPPUNIT_TEST( tetraNudged );
    CPPUNIT_TEST( tetraCorner );
    CPPUNIT_TEST( tetraSimpleIncluded );
    CPPUNIT_TEST( tetraDegenEdge );
    CPPUNIT_TEST( tetraDegenFace );
    CPPUNIT_TEST( tetraDegenTranslatedInPlane );
    CPPUNIT_TEST( tetraHalfstripOnly );
    CPPUNIT_TEST( tetraHalfstripOnly2 );
    CPPUNIT_TEST( tetraSimpleHalfstripOnly );
    CPPUNIT_TEST( generalTetra );
    CPPUNIT_TEST( trickyTetra1 );
    //    CPPUNIT_TEST( inconsistentTetra );

    CPPUNIT_TEST_SUITE_END();

  public:

    /// Unit tetrahedron mesh intersecting itself
    /// \brief Status : pass
    void tetraReflexiveUnit()
    {
      _testTools->intersectMeshes("UnitTetra", "UnitTetra", 1.0/6.0);
    }

    /// Tetrahedron mesh with itself
    /// \brief Status : pass
    void tetraReflexiveGeneral()
    {
      _testTools->intersectMeshes("GeneralTetra", "GeneralTetra", 0.428559);
    }

    /// Unit tetrahedron mesh intersecting slightly displaced copy of itself
    /// \brief Status : pass
    void tetraNudged()
    {
      _testTools->intersectMeshes("UnitTetra", "NudgedTetra", 0.142896);
    }

    /// Single-element unit tetrahedron mesh intersecting even slightly displaced (along one axis only) copy of itself
    /// \brief Status : pass
    void tetraNudgedSimpler()
    {
      _testTools->intersectMeshes("UnitTetra", "NudgedSimpler", 0.152112);
    }

    /// Tetrahedron intersecting unit tetrahedron with in non-degenerate way around corner O
    /// \brief Status : pass
    void tetraCorner()
    {
      _testTools->intersectMeshes("UnitTetra", "CornerTetra", 0.0135435);
    }

    /// Tetrahedron situated totally inside another
    /// \brief Status : pass
    void tetraSimpleIncluded()
    {
      _testTools->intersectMeshes("SimpleIncludedTetra", "SimpleIncludingTetra", 17.0156);
    }

    /// Displaced unit tetrahedron intersecting another unit tetrahedron with which it shares an edge
    /// \brief Status : pass
    void tetraDegenEdge()
    {
      _testTools->intersectMeshes("UnitTetraDegenT", "DegenEdgeXY", 0.0);
    }

    /// Displaced unit tetrahedron intersecting another unit tetrahedron with which it shares a face
    /// \brief Status : pass
    void tetraDegenFace()
    {
      _testTools->intersectMeshes("UnitTetraDegenT", "DegenFaceXYZ", 0.0);
    }

    /// Displaced unit tetrahedron intersecting another unit tetrahedron with which it shares a part of the face XYZ
    /// \brief Status : pass
    void tetraDegenTranslatedInPlane()
    {
      _testTools->intersectMeshes("UnitTetraDegenT", "DegenTranslatedInPlane", 0.0571667);
    }

    /// Tetrahedron having only half-strip intersections with the unit tetrahedron
    /// \brief Status : pass, but does not really test what it should - does not check that the intersections are detected. No longer needed.
    void tetraHalfstripOnly()
    {
      // NB this test is not completely significant : we should also verify that 
      // there are triangles on the element that give a non-zero volume
      _testTools->intersectMeshes("HalfstripOnly", "UnitTetra", 0.0);
    }

    /// Tetrahedron having only half-strip intersections with the unit tetrahedron
    /// \brief Status : pass, but does not really test what it should - does not check that the intersections are detected. No longer needed.
    void tetraHalfstripOnly2()
    {
      // NB this test is not completely significant : we should also verify that 
      // there are triangles on the element that give a non-zero volume
      _testTools->intersectMeshes("HalfstripOnly2", "UnitTetra", 0.0);
    }
  
    /// Tetrahedron having only half-strip intersections with the unit tetrahedron
    /// \brief Status : pass, but does not really test what it should - does not check that the intersections are detected. No longer needed.
    void tetraSimpleHalfstripOnly()
    {
      // NB this test is not completely significant : we should also verify that 
      // there are triangles on the element that give a non-zero volume
      _testTools->intersectMeshes("SimpleHalfstripOnly", "UnitTetra", 0.0);
    }

    /// Two intersecting tetrahedra situated in a general position in space
    /// \brief Status : pass
    void generalTetra()
    {
      _testTools->intersectMeshes("GenTetra1", "GenTetra2", 4.91393);
    }

    /// Tetrahedron which is in a tricky position relative to unit tetrahedron.
    /// \brief Status : pass
    void trickyTetra1()
    {
      _testTools->intersectMeshes("UnitTetra", "TrickyTetra1", 0.0);
    }

    /// Two large tetrahedra which nearly share part of an edge and intersect at the origin. Created with goal of getting the as-of-yet uncovered "consistency" test
    /// part of the correction of double products covered. However, it does not succeed with this.
    /// \brief Status : fails, but is quite far-fetched as far as typical use cases are concerned
    void inconsistentTetra()
    {
      _testTools->intersectMeshes("LargeUnitTetra.med", "LargeUnitTetra", "LargeInconsistentTetra.med", "LargeInconsistent", 7.86231e7);
    }

  };
}
#endif
