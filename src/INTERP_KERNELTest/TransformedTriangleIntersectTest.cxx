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

#include "TransformedTriangleIntersectTest.hxx"
#include <iostream>

#include "Log.hxx"

/// macro to test for zero double products outside the segment-edge intersection test method
/// as is done in TransformedTriangle when OPTIMIZE is defined
#define TEST_ZERO_DP_EDGE(seg, edge) isZero[TT::NO_DP*int(seg) + int(DoubleProduct(edge))]

/// macro to test for zero double products outside the segment-corner intersection test method
/// as is done in TransformedTriangle when OPTIMIZE is defined
#define TEST_ZERO_DP_CORNER(seg, corner)                                \
  isZero[DoubleProduct(TT::NO_DP*int(seg) +  TT::EDGES_FOR_CORNER[3*corner] )] && \
  isZero[DoubleProduct(TT::NO_DP*int(seg) +  TT::EDGES_FOR_CORNER[3*corner+1] )] && \
  isZero[DoubleProduct(TT::NO_DP*int(seg) +  TT::EDGES_FOR_CORNER[3*corner+2] )]

/// macro to test for zero double products outside the segment-ray intersection test method
/// as is done in TransformedTriangle when OPTIMIZE is defined
#define TEST_ZERO_DP_RAY(seg, corner) isZero[TT::NO_DP*int(seg) + TT::DP_SEGMENT_RAY_INTERSECTION[7*(corner-1)]]

using namespace INTERP_KERNEL;

namespace INTERP_TEST
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \class TransformedTriangleIntersectTest
  /// \brief Class testing the intersection detection methods of TransformedTriangle.
  ///
  /// This class contains unit tests for the intersection methods of the TransformedTriangle class.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Each method in the class runs all the intersection tests with some triangle. The goal is to cover all
  /// the different types of intersections between a triangle and a tetrahedron. The table below gives a 
  /// a summary of what is being tested. Before each method, there is also a summary of what how the 
  /// triangle in the method intersects the unit tetrahedron.
  /// 
  /// Since performing all tests would require a large number of triangles, we have limited our coverage to 
  /// be such that each column and each row in the table below has at least one entry for each type of 
  /// intersection. The intersection forumlae are totally symmetric with respect to changing the segment
  /// (PQ, QR, or RP) of the triangle, so they only enter in a very simple way in the code. Testing 
  ///  all these cases is therefore of low priority.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Intersections tested (number indicates first triangle which contains the intersection):
  /// <PRE>
  /// -----------------------------------------------------------------------------------------------------
  /// CI  ->  P: 3      Q: 4     R: 7
  /// COH ->  P: 9      Q: 8     R: 10
  /// CAH ->  P: 4      Q: 10    R: 9
  /// -----------------------------------------------------------------------------------------------------
  /// SF  ->  (PQ, OZX) : 1   (PQ, OYZ) : 2   (PQ, OXY) : 1   (PQ, XYZ) : 3
  ///     ->  (QR, OZX) : 8   (QR, OYZ) : -   (QR, OXY) : 4   (QR, XYZ) : 7
  ///     ->  (RP, OZX) : 1   (RP, OYZ) : 3   (RP, OXY) : 7   (RP, XYZ) : 1
  /// -----------------------------------------------------------------------------------------------------
  /// SE  ->  (PQ, OX)  : 11  (PQ, OY)  : -   (PQ, OZ)  : 12  (PQ, XY)  : 2   (PQ, ZX)  : -  (PQ, YZ)  : 10
  ///     ->  (QR, OX)  : -   (QR, OY)  : -   (QR, OZ)  : -   (QR, XY)  : -   (QR, ZX)  : 9  (QR, YZ)  : -
  ///     ->  (RP, OX)  : -   (RP, OY)  : 12  (RP, OZ)  : -   (RP, XY)  : -   (RP, ZX)  : -  (RP, YZ)  : -
  /// -----------------------------------------------------------------------------------------------------
  /// SC  ->  (PQ, O)   : -   (PQ, X)   : -   (PQ, Y)   : 8   (PQ, Z)   : -
  ///     ->  (QR, O)   : -   (QR, X)   : 2   (QR, Y)   : -   (QR, Z)   : 13
  ///     ->  (RP, O)   : 11  (RP, X)   : -   (RP, Y)   : -   (RP, Z)   : -
  /// -----------------------------------------------------------------------------------------------------
  /// SHS ->  (PQ, XY)  : 3   (PQ, ZX)  : -   (PQ, YZ)  : 13
  ///     ->  (QR, XY)  : 3   (QR, ZX)  : 5   (QR, YZ)  : 3
  ///     ->  (RP, XY)  : 1   (RP, ZX)  : 4   (RP, YZ)  : -
  /// -----------------------------------------------------------------------------------------------------
  /// SR  ->  (PQ, X)   : 6   (PQ, Y)   : 5   (PQ, Z)   : -
  ///     ->  (QR, X)   : -   (QR, Y)   : -   (QR, Z)   : 6
  ///     ->  (RP, X)   : -   (RP, Y)   : -   (RP, Z)   : -
  /// -----------------------------------------------------------------------------------------------------
  /// TE  ->  OX : 4   OY : 7    OZ : 8     XY : 1     ZX : 4    YZ : 3
  /// -----------------------------------------------------------------------------------------------------
  /// TR  ->  X  : 7    Y : 6     Z : 5
  /// -----------------------------------------------------------------------------------------------------
  /// </PRE>
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Key to triangle descriptions : 
  /// CI  = Triangle corner contained in tetrahedron 
  /// COH = Triangle corner on h = 0 face of tetrahedron
  /// CAH = Triangle corner above h = 0 face of tetrahedron in z-direction
  /// SF  = Segment - facet intersection
  /// SE  = Segment - edge intersection
  /// SC  = Segment - corner intersection
  /// SHS = Segment - halfstrip intersection
  /// SR  = Segment - ray intersection
  /// TE  = Tetrahedron edge intersects triangle (surface - edge intersection)
  /// TR  = Surface - ray intersection
  ///
  /// In the descriptions for each triangle, square brackets indicate superfluous but allowed intersections
  /// that arise as by-products of for instance segment-corner intersections.
  /// E.g. A segment - corner intersection can imply three surface - edge intersections
  /// Since these "extra" intersections arise under special circumstances, they are not counted in the 
  /// table above
  ////////////////////////////////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Triangle 1 has the following intersections
  /// <PRE>
  /// CI     -
  /// COH    -
  /// CAH    -
  /// SF     (PQ, OXY), (PQ, OZX), (RP, XYZ), (RP, OZX)
  /// SE     -
  /// SC     - 
  /// SHS    (RP, XY)
  /// SR     - 
  /// TE     XY
  /// TR     - 
  /// </PRE>
  /// \brief Status : pass
  void TransformedTriangleIntersectTest::testTriangle1()
  {
    LOG(1, "+++++++ Testing triangle 1" );

    typedef TransformedTriangle TT;

    double coords[9] = 
      {
        0.4,-0.5, 0.5, // P
        0.4, 2.5,-1.0, // Q
        0.4, 2.5, 0.5  // R
      };

    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false

    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));

    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));
  
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));
  
    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT(tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

    delete tri;

  }

  /// Triangle 2 has the following intersections
  /// <PRE>
  /// CI     -
  /// COH    -
  /// CAH    -
  /// SF     (PQ, OYZ)
  /// SE     (PQ, XY)
  /// SC     (QR, X)
  /// SHS    -
  /// SR     - 
  /// TE     [OX, OZ, ZX]
  /// TR     - 
  /// </PRE>
  /// \brief Status: pass
  void TransformedTriangleIntersectTest::testTriangle2()
  {
    LOG(1, "+++++++ Testing triangle 2" );
    typedef TransformedTriangle TT;

    double coords[9] =
      {
        -0.5, 0.5, 0.25, // P
        1.5, 0.5,-0.25, // Q
        -0.5,-1.5, 0.75  // R
      };
    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false

    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));  

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));

    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));
    
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));

    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

    delete tri;
  }

  /// Triangle 3 has the following intersections
  /// <PRE>
  /// CI     P
  /// COH    -
  /// CAH    -
  /// SF     (PQ, XYZ), (RP, OYZ)
  /// SE     -
  /// SC     -
  /// SHS    (PQ, XY), (QR, YZ), (QR, XY)
  /// SR     - 
  /// TE     YZ
  /// TR     - 
  /// </PRE>
  /// \brief Status : pass
  void TransformedTriangleIntersectTest::testTriangle3()
  {
    LOG(1, "+++++++ Testing triangle 3" );
    typedef TransformedTriangle TT;

    double coords[9] =
      {
        0.35, 0.15, 0.1, // P
        0.8, 0.8, 0.8,  // Q
        -0.4, 0.3, 0.9   // R
      };

    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false

    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));  

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));
  
    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));
    
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));

    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

    delete tri;
  }

  /// Triangle 4 has the following intersections
  /// <PRE>
  /// CI     Q
  /// COH    -
  /// CAH    P
  /// SF     (PQ, XYZ), (QR, OXY)
  /// SE     -
  /// SC     -
  /// SHS    (RP, ZX)
  /// SR     - 
  /// TE     (OX, ZX)
  /// TR     - 
  /// </PRE>
  /// \brief Status : pass
  void TransformedTriangleIntersectTest::testTriangle4()
  {
    LOG(1, "+++++++ Testing triangle 4" );
    typedef TransformedTriangle TT;

    double coords[9] =
      {
        0.3, 0.3, 1.8,  // P
        0.75, 0.1, 0.1,  // Q
        0.2, -1.3, -1.4   // R
      };

    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false
  
    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));  

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));

    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));  

  
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));

    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

    delete tri;
  }

  /// Triangle 5 has the following intersections
  /// <PRE>
  /// CI     -
  /// COH    -
  /// CAH    -
  /// SF     -
  /// SE     -
  /// SC     -
  /// SHS    (QR, ZX), (QR, XY)
  /// SR     (PQ, Y)
  /// TE     -
  /// TR     Z
  /// </PRE>
  /// \brief Status : pass
  void TransformedTriangleIntersectTest::testTriangle5()
  {
    LOG(1, "+++++++ Testing triangle 5" );
  
    typedef TransformedTriangle TT;

    double coords[9] =
      {
        -0.5, 0.5, 2.3,  // P
        0.5, 1.5, 2.8,  // Q
        0.5, -2.6, 1.3   // R
      };

    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false

    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));  

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));

    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));
    
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));
  
    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceRayIntersection(TT::Z));

    delete tri;
  }

  /// Triangle 6 has the following intersections
  /// <PRE>
  /// CI     -
  /// COH    -
  /// CAH    -
  /// SF     -
  /// SE     -
  /// SC     -
  /// SHS    -
  /// SR     (PQ, X), (QR, Z) 
  /// TE     -
  /// TR     Y 
  /// </PRE>
  /// \brief Status : pass
  void TransformedTriangleIntersectTest::testTriangle6()
  {
    LOG(1, "+++++++ Testing triangle 6" );

    typedef TransformedTriangle TT;
  
    double coords[9] =
      {
        1.5, 0.5, 1.35,  // P
        0.5, -0.5, 2.1,  // Q
        -3.0, 3.0, -0.5   // R
      };
  
    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false

    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));  

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));

    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));
    
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));
  
    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

    delete tri;
  }

  /// Triangle 7 has the following intersections
  /// <PRE>
  /// CI     R
  /// COH    -
  /// CAH    -
  /// SF     (RP, OXY),(QR,XYZ)
  /// SE     -
  /// SC     -
  /// SHS    (QR, XY)
  /// SR     - 
  /// TE     OX, ZX
  /// TR     X 
  /// </PRE>
  /// \brief Status : pass
  void TransformedTriangleIntersectTest::testTriangle7()
  {

    LOG(1, "+++++++ Testing triangle 7" );
    typedef TransformedTriangle TT;

    double coords[9] =
      {
        -2.3, -1.5, -2.5,  // P
        3.1, 0.15, 0.8,  // Q
        0.3, 0.4, 0.2   // R
      };

    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false

    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));  

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));
 
    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));
 
    
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));
  
    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

    delete tri;
  }

  /// Triangle 8 has the following intersections
  /// <PRE>
  /// CI     [Q] 
  /// COH    Q
  /// CAH    -
  /// SF     (QR, OZX), [ (QR, XYZ) ]
  /// SE     -
  /// SC     (PQ,Y)
  /// SHS    -
  /// SR     - 
  /// TE     OZ, [YZ,OY,XY]
  /// TR     - 
  /// </PRE>
  /// \brief Status : pass
  void TransformedTriangleIntersectTest::testTriangle8()
  {
    LOG(1, "+++++++ Testing triangle 8" );
    typedef TransformedTriangle TT;

    double coords[9] =
      {
        -0.75, 3.25, -1.5,  // P
        0.25, 0.25, 0.5,  // Q
        -0.1, -0.4, 0.9   // R
      };

    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false

    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));  

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));

    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));
    
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));

    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));  

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

    delete tri;
  }

  /// Triangle 9 has the following intersections
  /// <PRE>
  /// CI     [P]
  /// COH    P
  /// CAH    R
  /// SF     (PQ, OZX), [(PQ, XYZ), (RP,XYZ)]
  /// SE     (QR, ZX)
  /// SC     -
  /// SHS    -
  /// SR     - 
  /// TE     [ZX]
  /// TR     - 
  /// </PRE>
  /// \brief Status : pass
  void TransformedTriangleIntersectTest::testTriangle9()
  {
    LOG(1, "+++++++ Testing triangle 9" );
    typedef TransformedTriangle TT;

    double coords[9] =
      {
        0.6, 0.2, 0.2,  // P
        0.3, -0.2, 0.8,  // Q
        0.1, 0.2, 0.8   // R
      };

    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false

    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));  

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));

    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));
    
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));
  
    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

    delete tri;
  }

  
  /// Triangle 10 has the following intersections
  /// <PRE>
  /// CI     [R]
  /// COH    R
  /// CAH    Q
  /// SF     (RP, OYZ), [ (RP, XYZ), (QR, XYZ) ]
  /// SE     (PQ, YZ)
  /// SC     -
  /// SHS    -
  /// SR     - 
  /// TE     [YZ]
  /// TR     - 
  /// </PRE>
  /// \brief Status : pass
  void TransformedTriangleIntersectTest::testTriangle10()
  {
    LOG(1, "+++++++ Testing triangle 10" );
    typedef TransformedTriangle TT;

    double coords[9] =
      {
        -0.1, 0.3, 0.6,  // P
        0.1, 0.1, 1.0,  // Q
        0.4, 0.3, 0.3   // R
      };

    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false

    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerOnXYZFacet(TT::R));  

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));

    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));
    
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));
  
    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

    delete tri;
  }

  /// Triangle 11 has the following intersections
  /// <PRE>
  /// CI     Q, R
  /// COH    -
  /// CAH    -
  /// SF     -
  /// SE     (PQ, OX)
  /// SC     (RP, O)
  /// SHS    -
  /// SR     - 
  /// TE     [OY, OZ]
  /// TR     - 
  /// </PRE>
  /// \brief Status : pass
  void TransformedTriangleIntersectTest::testTriangle11()
  {
    LOG(1, "+++++++ Testing triangle 11" );  
    typedef TransformedTriangle TT;
  
    double coords[9] =
      {
        -0.2, -0.2, -0.2,  // P
        0.2, 0.1, 0.1,  // Q
        0.3, 0.3, 0.3   // R
      };
  
    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false

    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(true , tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));  

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));

    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));
    
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));
  
    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(true, tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

    delete tri;
  }


  /// Triangle 12 has the following intersections
  /// <PRE>
  /// CI     -
  /// COH    -
  /// CAH    -
  /// SF     (QR, OXY), (QR, OZX)
  /// SE     (RP, OY), (PQ, OZ)
  /// SC     -
  /// SHS    -
  /// SR     - 
  /// TE     [OY], [OZ]
  /// TR     - 
  /// </PRE>
  /// \brief Status : pass
  void TransformedTriangleIntersectTest::testTriangle12()
  {
    LOG(1, "+++++++ Testing triangle 12" );
    typedef TransformedTriangle TT;

    double coords[9] =
      {
        -0.2, 0.2, 0.2,  // P
        0.2, -0.2, 0.3,  // Q
        0.6, 0.6, -0.6   // R
      };

    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false

    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));  

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));
    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));
    
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));
  
    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

    delete tri;
  }


  /// Triangle 13 has the following intersections
  /// <PRE>
  /// CI     -
  /// COH    -
  /// CAH    -
  /// SF     (QR, OYZ), (PQ, OXY), (PQ, XYZ)
  /// SE     -
  /// SC     (QR, Z)
  /// SHS    (PQ, YZ)
  /// SR     - 
  /// TE     [OZ, YZ, ZX]
  /// TR     - 
  /// </PRE>
  /// \brief Status : pass
  void TransformedTriangleIntersectTest::testTriangle13()
  {
    LOG(1, "+++++++ Testing triangle 13" );
    typedef TransformedTriangle TT;

    double coords[9] =
      {
        -0.2, 0.3, 5.0,  // P
        0.2, 0.1, -1.0,  // Q
        -0.2, -0.1, 3.0   // R
      };

    TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

    // run all intersection tests and ensure that the ones
    // listed with yes in the tables above return true and 
    // that the ones listed with no or not listed at all return false

    bool isZero[TT::NO_TRI_SEGMENT * TT::NO_DP];
  
    for(TriSegment seg = TT::PQ ; seg < TT::NO_TRI_SEGMENT ; seg = TT::TriSegment(seg + 1))
      {
        // check beforehand which double-products are zero
        for(DoubleProduct dp = TT::C_YZ; dp < TT::NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[TT::NO_DP*int(seg) + int(dp)] = (tri->calcStableC(seg, dp) == 0.0);
          }
      }

    // corner in tetrahedron (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::R));
  
    // corner on XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));  

    // corner above XYZ facet (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
    CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

    // segment-facet (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));

    // segment-edge (18 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::OZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::YZ) && tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::ZX) && tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::PQ, TT::XY) && tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OX) && tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OY) && tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::OZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::YZ) && tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::ZX) && tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::QR, TT::XY) && tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OX) && tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OY) && tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::OZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::YZ) && tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::ZX) && tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_EDGE(TT::RP, TT::XY) && tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

    // segment - corner (12 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::O) && tri->testSegmentCornerIntersection(TT::PQ, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::X) && tri->testSegmentCornerIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Y) && tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::PQ, TT::Z) && tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::O) && tri->testSegmentCornerIntersection(TT::QR, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::X) && tri->testSegmentCornerIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::QR, TT::Y) && tri->testSegmentCornerIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(true , TEST_ZERO_DP_CORNER(TT::QR, TT::Z) && tri->testSegmentCornerIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::O) && tri->testSegmentCornerIntersection(TT::RP, TT::O));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::X) && tri->testSegmentCornerIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Y) && tri->testSegmentCornerIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_CORNER(TT::RP, TT::Z) && tri->testSegmentCornerIntersection(TT::RP, TT::Z));
    
    // segment-halfstrip (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(true , tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));
  
    // segment-ray (9 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::X) && tri->testSegmentRayIntersection(TT::PQ, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Y) && tri->testSegmentRayIntersection(TT::PQ, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::PQ, TT::Z) && tri->testSegmentRayIntersection(TT::PQ, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::X) && tri->testSegmentRayIntersection(TT::QR, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Y) && tri->testSegmentRayIntersection(TT::QR, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::QR, TT::Z) && tri->testSegmentRayIntersection(TT::QR, TT::Z));

    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::X) && tri->testSegmentRayIntersection(TT::RP, TT::X));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Y) && tri->testSegmentRayIntersection(TT::RP, TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, TEST_ZERO_DP_RAY(TT::RP, TT::Z) && tri->testSegmentRayIntersection(TT::RP, TT::Z));

    // surface-edge (6 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OY));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::OZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::YZ));
    CPPUNIT_ASSERT_EQUAL(true , tri->testSurfaceEdgeIntersection(TT::ZX));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::XY));

    // surface-ray (3 possibilities)
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
    CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

    delete tri;
  }


} // NAMESPACE 







///// TEMPLATE ///////////////////////////////



#if 0
// Triangle x has the following intersections
// CI     -
// COH    -
// CAH    -
// SF     -
// SE     -
// SC     -
// SHS    -
// SR     - 
// TE     -
// TR     - 

void TransformedTriangleIntersectTest::testTriangleX()
{
  LOG(1, "+++++++ Testing triangle X" );
  typedef TransformedTriangle TT;

  double coords[9] =
    {
      0.0, 0.0, 0.0,  // P
      0.0, 0.0, 0.0,  // Q
      0.0, 0.0, 0.0   // R
    };

  TransformedTriangle* tri = new TransformedTriangle(&coords[0], &coords[3], &coords[6]);

  // run all intersection tests and ensure that the ones
  // listed with yes in the tables above return true and 
  // that the ones listed with no or not listed at all return false

  // corner in tetrahedron (3 possibilities)
  CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::P));
  CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::Q));
  CPPUNIT_ASSERT_EQUAL(false, tri->testCornerInTetrahedron(TT::R));
  
  // corner on XYZ facet (3 possibilities)
  CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::P));
  CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::Q));
  CPPUNIT_ASSERT_EQUAL(false, tri->testCornerOnXYZFacet(TT::R));  

  // corner above XYZ facet (3 possibilities)
  CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::P));
  CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::Q));
  CPPUNIT_ASSERT_EQUAL(false, tri->testCornerAboveXYZFacet(TT::R));  

  // segment-facet (9 possibilities)
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OZX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OYZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::OXY));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::PQ, TT::XYZ));

  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OZX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OYZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::OXY));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::QR, TT::XYZ));

  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OZX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OYZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::OXY));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentFacetIntersection(TT::RP, TT::XYZ));

  // segment-edge (18 possibilities)
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::PQ, TT::OX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::PQ, TT::OY));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::PQ, TT::OZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::PQ, TT::YZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::PQ, TT::ZX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::PQ, TT::XY));

  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::QR, TT::OX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::QR, TT::OY));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::QR, TT::OZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::QR, TT::YZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::QR, TT::ZX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::QR, TT::XY));
  
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::RP, TT::OX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::RP, TT::OY));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::RP, TT::OZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::RP, TT::YZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::RP, TT::ZX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentEdgeIntersection(TT::RP, TT::XY));

  // segment - corner (12 possibilities)
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentCornerIntersection(TT::PQ, TT::O));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentCornerIntersection(TT::PQ, TT::X));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentCornerIntersection(TT::PQ, TT::Y));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentCornerIntersection(TT::PQ, TT::Z));

  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentCornerIntersection(TT::QR, TT::O));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentCornerIntersection(TT::QR, TT::X));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentCornerIntersection(TT::QR, TT::Y));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentCornerIntersection(TT::QR, TT::Z));

  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentCornerIntersection(TT::RP, TT::O));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentCornerIntersection(TT::RP, TT::X));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentCornerIntersection(TT::RP, TT::Y));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentCornerIntersection(TT::RP, TT::Z));
    
  // segment-halfstrip (9 possibilities)
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::YZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::ZX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::PQ, TT::XY));

  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::YZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::ZX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::QR, TT::XY));
  
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::YZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::ZX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentHalfstripIntersection(TT::RP, TT::XY));
  
  // segment-ray (9 possibilities)
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentRayIntersection(TT::PQ, TT::X));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentRayIntersection(TT::PQ, TT::Y));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentRayIntersection(TT::PQ, TT::Z));

  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentRayIntersection(TT::QR, TT::X));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentRayIntersection(TT::QR, TT::Y));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentRayIntersection(TT::QR, TT::Z));

  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentRayIntersection(TT::RP, TT::X));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentRayIntersection(TT::RP, TT::Y));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSegmentRayIntersection(TT::RP, TT::Z));

  // surface-edge (6 possibilities)
  CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OY));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::OZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::YZ));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::ZX));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceEdgeIntersection(TT::XY));

  // surface-ray (3 possibilities)
  CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::X));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Y));
  CPPUNIT_ASSERT_EQUAL(false, tri->testSurfaceRayIntersection(TT::Z));

  delete tri;
}
#endif
