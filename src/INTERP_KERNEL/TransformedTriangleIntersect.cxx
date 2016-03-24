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

#include "TransformedTriangle.hxx"
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include "VectorUtils.hxx"

namespace INTERP_KERNEL
{

  // ----------------------------------------------------------------------------------
  //  Correspondance tables describing all the variations of formulas. 
  // ----------------------------------------------------------------------------------

  /// \brief Correspondance between facets and double products.
  ///
  /// This table encodes Grandy, table IV. Use 3*facet + {0,1,2} as index
  const TransformedTriangle::DoubleProduct TransformedTriangle::DP_FOR_SEG_FACET_INTERSECTION[12] = 
    {
      C_XH, C_XY, C_ZX, // OYZ
      C_YH, C_YZ, C_XY, // OZX
      C_ZH, C_ZX, C_YZ, // OXY
      C_XH, C_YH, C_ZH  // XYZ
    };

  /// \brief Signs associated with entries in DP_FOR_SEGMENT_FACET_INTERSECTION.
  /// 
  /// This table encodes Grandy, table IV. Use 3*facet + {0,1,2} as index
  const double TransformedTriangle::SIGN_FOR_SEG_FACET_INTERSECTION[12] = 
    {
      1.0, 1.0, -1.0,
      1.0, 1.0, -1.0,
      1.0, 1.0, -1.0,
      1.0, 1.0,  1.0
    };

  /// \brief Coordinates of corners of tetrahedron.
  ///
  /// Use 3*Corner + coordinate as index
  const double TransformedTriangle::COORDS_TET_CORNER[12] = 
    {
      0.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0
    };

  /// \brief Indices to use in tables DP_FOR_SEG_FACET_INTERSECTION and SIGN_FOR_SEG_FACET_INTERSECTION
  /// for the calculation of the coordinates (x,y,z) of the intersection points
  /// for Segment-Facet and Segment-Edge intersections.
  ///
  /// Use 3*facet + coordinate as index. -1 indicates that the coordinate is 0.
  const int TransformedTriangle::DP_INDEX[12] =
    {
      // x, y, z
      -1, 1, 2,  // OYZ
      5, -1, 4,  // OZX
      7, 8, -1,  // OXY
      9, 10, 11  // XYZ
    };

  /// \brief Correspondance edge - corners.
  ///
  /// Gives the two corners associated with each edge
  /// Use 2*edge + {0, 1} as index
  const TransformedTriangle::TetraCorner TransformedTriangle::CORNERS_FOR_EDGE[12] = 
    {
      O, X, // OX
      O, Y, // OY
      O, Z, // OZ
      X, Y, // XY
      Y, Z, // YZ
      Z, X  // ZX
    };

  /// \brief Correspondance edge - facets.
  ///
  /// Gives the two facets shared by and edge. Use 2*facet + {0, 1} as index
  const TransformedTriangle::TetraFacet TransformedTriangle::FACET_FOR_EDGE[12] =
    {
      OXY, OZX, // OX
      OXY, OYZ, // OY
      OZX, OYZ, // OZ
      OXY, XYZ, // XY
      OYZ, XYZ, // YZ
      OZX, XYZ  // ZX
    };

  /// \brief Correspondance corners - edges.
  ///
  /// Gives edges meeting at a given corner. Use 3*corner + {0,1,2} as index
  const TransformedTriangle::TetraEdge TransformedTriangle::EDGES_FOR_CORNER[12] =
    {
      OX, OY, OZ, // O
      OX, XY, ZX, // X
      OY, XY, YZ, // Y
      OZ, ZX, YZ  // Z
    };

  /// \brief Double products to use in halfstrip intersection tests.
  ///
  /// Use 4*(offset_edge) + {0,1,2,3} as index. offset_edge = edge - 3  (so that XY -> 0, YZ -> 1, ZX -> 2)
  /// Entries with offset 0 and 1 are for the first condition (positive product) 
  /// and those with offset 2 and 3 are for the second condition (negative product).
  const TransformedTriangle::DoubleProduct TransformedTriangle::DP_FOR_HALFSTRIP_INTERSECTION[12] =
    {
      C_10, C_01, C_ZH, C_10, // XY
      C_01, C_XY, C_XH, C_01, // YZ
      C_XY, C_10, C_YH, C_XY  // ZX
    };
  
  /// \brief Double products to use in segment-ray test.
  ///
  /// Use 7*corner_offset + {0,1,2,3,4,5,6} as index. corner_offset = corner - 1 (so that X -> 0, Y-> 1, Z->2)
  /// Entries with offset 0 are for first condition (zero double product) and the rest are for condition 3 (in the same
  /// order as in the article)
  const TransformedTriangle::DoubleProduct TransformedTriangle::DP_SEGMENT_RAY_INTERSECTION[21] = 
    {
      C_10, C_YH, C_ZH, C_01, C_XY, C_YH, C_XY, // X
      C_01, C_XH, C_ZH, C_XY, C_10, C_ZH, C_10, // Y
      C_XY, C_YH, C_XH, C_10, C_01, C_XH, C_01  // Z
    };

  /**
   * Calculates the point of intersection between the given edge of the tetrahedron and the 
   * triangle PQR. (Grandy, eq [22])
   *
   * @pre   testSurfaceEdgeIntersection(edge) returns true
   * @param edge   edge of tetrahedron
   * @param pt     array of three doubles in which to store the coordinates of the intersection point
   */
  void TransformedTriangle::calcIntersectionPtSurfaceEdge(const TetraEdge edge, double* pt) const
  {
    assert(edge < H01);

    // barycentric interpolation between points A and B 
    // : (x,y,z)* = (1-alpha)*A + alpha*B where
    // alpha = t_A / (t_A - t_B)
    
    const TetraCorner corners[2] = 
      {
        CORNERS_FOR_EDGE[2*edge],
        CORNERS_FOR_EDGE[2*edge + 1]
      };
    
    // calculate alpha
    const double tA = calcStableT(corners[0]);
    const double tB = calcStableT(corners[1]);
    const double alpha = tA / (tA - tB);

    // calculate point
    LOG(4, "corner A = " << corners[0] << " corner B = " << corners[1] );
    LOG(4, "tA = " << tA << " tB = " << tB << " alpha= " << alpha );
    for(int i = 0; i < 3; ++i)
      {

        pt[i] = (1 - alpha) * COORDS_TET_CORNER[3*corners[0] + i] + 
          alpha * COORDS_TET_CORNER[3*corners[1] + i];
#if 0
        pt[i] = (1 - alpha) * getCoordinateForTetCorner<corners[0], i>() + 
          alpha * getCoordinateForTetCorner<corners[0], i>();
#endif
        LOG(6, pt[i] );
        assert(pt[i] >= 0.0);
        assert(pt[i] <= 1.0);
      }
  }

  /**
   * Calculates the point of intersection between the given segment of the triangle
   * and the given facet of the tetrahedron. (Grandy, eq. [23])
   *
   * @pre   testSurfaceEdgeIntersection(seg, facet) returns true
   * 
   * @param seg    segment of the triangle
   * @param facet  facet of the tetrahedron
   * @param pt     array of three doubles in which to store the coordinates of the intersection point
   */
  void TransformedTriangle::calcIntersectionPtSegmentFacet(const TriSegment seg, const TetraFacet facet, double* pt) const
  {
    // calculate s
    double s = 0.0;
    for(int i = 0; i < 3; ++i)
      {
        const DoubleProduct dp = DP_FOR_SEG_FACET_INTERSECTION[3*facet + i];
        const double sign = SIGN_FOR_SEG_FACET_INTERSECTION[3*facet + i];
        s -= sign * calcStableC(seg, dp);
      }

    assert(s != 0.0);

    // calculate coordinates of intersection point
    for(int i = 0 ; i < 3; ++i)
      {
        const int dpIdx = DP_INDEX[3*facet + i];
       
        if(dpIdx < 0)
          {
            pt[i] = 0.0;
          }
        else
          {
            const DoubleProduct dp = DP_FOR_SEG_FACET_INTERSECTION[dpIdx];
            const double sign = SIGN_FOR_SEG_FACET_INTERSECTION[dpIdx];
            pt[i] = -( sign * calcStableC(seg, dp) ) / s;

            LOG(4, "SegmentFacetIntPtCalc : pt[" << i << "] = " << pt[i]  );
            LOG(4, "c(" << seg << ", " << dp << ") = " <<  sign * calcStableC(seg, dp) );
            assert(pt[i] >= 0.0); 
            assert(pt[i] <= 1.0);
          }
      }
  
  }

  /**
   * Tests if the given segment of the triangle intersects the given edge of the tetrahedron (Grandy, eq. [20]
   * If the OPTIMIZE is defined, it does not do the test the double product that should be zero.
   * @param seg    segment of the triangle
   * @param edge   edge of tetrahedron
   * @return      true if the segment intersects the edge 
   */
  bool TransformedTriangle::testSegmentEdgeIntersection(const TriSegment seg, const TetraEdge edge) const
  {
      {
        // check condition that the double products for one of the two
        // facets adjacent to the edge has a positive product
        bool isFacetCondVerified = false;
        TetraFacet facet[2];
        for(int i = 0 ; i < 2 ; ++i) 
          {
            facet[i] = FACET_FOR_EDGE[2*edge + i];
           
            // find the two c-values -> the two for the other edges of the facet
            int idx1 = 0 ; 
            int idx2 = 1;
            DoubleProduct dp1 = DP_FOR_SEG_FACET_INTERSECTION[3*facet[i] + idx1];
            DoubleProduct dp2 = DP_FOR_SEG_FACET_INTERSECTION[3*facet[i] + idx2];
           
            if(dp1 == DoubleProduct( edge ))
              {
                idx1 = 2;
                dp1 = DP_FOR_SEG_FACET_INTERSECTION[3*facet[i] + idx1];
              }
            else if(dp2 == DoubleProduct( edge ))
              {
                idx2 = 2;
                dp2 = DP_FOR_SEG_FACET_INTERSECTION[3*facet[i] + idx2];
              }
           
            const double c1 = SIGN_FOR_SEG_FACET_INTERSECTION[3*facet[i] + idx1]*calcStableC(seg, dp1);
            const double c2 = SIGN_FOR_SEG_FACET_INTERSECTION[3*facet[i] + idx2]*calcStableC(seg, dp2);

            //isFacetCondVerified = isFacetCondVerified || c1*c2 > 0.0;
            if(c1*c2 > 0.0)
              {
                isFacetCondVerified = true;
              }
          }

        if(!isFacetCondVerified)
          {
            return false;
          }
        else
          {
            return testSegmentIntersectsFacet(seg, facet[0]) || testSegmentIntersectsFacet(seg, facet[1]);
          }
      }
  }
    
  /**
   * Calculates the point of intersection between the given segment of the triangle
   * and the given edge of the tetrahedron. (Grandy, eq. [25])
   *
   * @pre   testSegmentEdgeIntersection(seg, edge) returns true
   * 
   * @param seg    segment of the triangle
   * @param edge   edge of the tetrahedron
   * @param pt     array of three doubles in which to store the coordinates of the intersection point
   */
  void TransformedTriangle::calcIntersectionPtSegmentEdge(const TriSegment seg, const TetraEdge edge, double* pt) const 
  {
    assert(edge < H01);

    // get the two facets associated with the edge
    static const TetraFacet FACETS_FOR_EDGE[12] =
      {
        OXY, OZX, // OX
        OXY, OYZ, // OY
        OZX, OYZ, // OZ
        OXY, XYZ, // XY
        OYZ, XYZ, // YZ
        OZX, XYZ  // ZX
      };

    const TetraFacet facets[2] = 
      {
        FACETS_FOR_EDGE[2*edge],
        FACETS_FOR_EDGE[2*edge + 1]
      };
    
    // calculate s for the two edges
    double s[2];
    for(int i = 0; i < 2; ++i)
      {
        s[i] = 0.0;
        for(int j = 0; j < 3; ++j)
          {
            const DoubleProduct dp = DP_FOR_SEG_FACET_INTERSECTION[3*facets[i] + j];
            const double sign = SIGN_FOR_SEG_FACET_INTERSECTION[3*facets[i] + j];
            s[i] += sign * calcStableC(seg, dp);
          }
      }

    // calculate denominator
    const double denominator = s[0]*s[0] + s[1]*s[1];

    // calculate intersection point
    for(int i = 0; i < 3; ++i)
      {
        // calculate double product values for the two faces
        double c[2];
        for(int j = 0 ; j < 2; ++j)
          {
            const int dpIdx = DP_INDEX[3*facets[j] + i];
            const DoubleProduct dp = DP_FOR_SEG_FACET_INTERSECTION[dpIdx];
            const double sign = SIGN_FOR_SEG_FACET_INTERSECTION[dpIdx];
            c[j] = dpIdx < 0.0 ? 0.0 : sign * calcStableC(seg, dp);
          }
       
        // pt[i] = (c1*s1 + c2*s2) / (s1^2 + s2^2)

        pt[i] = (c[0] * s[0] + c[1] * s[1]) / denominator;
       
        // strange bug with -O2 enabled : assertion fails when we don't have the following
        // trace - line
        //std::cout << "pt[i] = " << pt[i] << std::endl;
        //assert(pt[i] >= 0.0); // check we are in tetraeder
        //assert(pt[i] <= 1.0);
       
      }
  }

    
  /**
   * Tests if the given segment of the triangle intersects the given corner of the tetrahedron.
   * (Grandy, eq. [21]). If OPTIMIZE is defined, the double products that should be zero are not verified.
   *
   * @param seg    segment of the triangle
   * @param corner corner of the tetrahedron
   * @return      true if the segment intersects the corner
   */
  bool TransformedTriangle::testSegmentCornerIntersection(const TriSegment seg, const TetraCorner corner) const 
  {
    

    // facets meeting at a given corner
    static const TetraFacet FACETS_FOR_CORNER[12] =
      {
        OXY, OYZ, OZX, // O
        OZX, OXY, XYZ, // X
        OYZ, XYZ, OXY, // Y
        OZX, XYZ, OYZ  // Z
      };
    
    // check segment intersect a facet
    for(int i = 0 ; i < 3 ; ++i)
      {
        const TetraFacet facet = FACETS_FOR_CORNER[3*corner + i];
        if(testSegmentIntersectsFacet(seg, facet))
          {
            return true;
          }
      }
    
    return false;
  }

  /**
   * Tests if the given segment of the triangle intersects the half-strip above the 
   * given edge of the h = 0 plane. (Grandy, eq. [30])
   * 
   * @param seg    segment of the triangle
   * @param edge   edge of the h = 0 plane of the tetrahedron (XY, YZ, ZX)
   * @return      true if the upwards ray from the corner intersects the triangle. 
   */
  bool TransformedTriangle::testSegmentHalfstripIntersection(const TriSegment seg, const TetraEdge edge)
  {
    // get right index here to avoid "filling out" array
    const int edgeIndex = static_cast<int>(edge) - 3;
    
    // double products used in test
    // products 1 and 2 for each edge -> first condition in Grandy [30]
    // products 3 and 4 for each edge -> third condition
    // NB : some uncertainty whether these last are correct
    // DP_FOR_HALFSTRIP_INTERSECTION
    
    // facets to use in second condition (S_m)
    static const TetraFacet FACET_FOR_HALFSTRIP_INTERSECTION[3] = 
      {
        NO_TET_FACET, // XY -> special case : test with plane H = 0
        OYZ, // YZ
        OZX  // ZX
      };

    const double cVals[4] = 
      {
        calcStableC(seg, DP_FOR_HALFSTRIP_INTERSECTION[4*edgeIndex]),
        calcStableC(seg, DP_FOR_HALFSTRIP_INTERSECTION[4*edgeIndex + 1]),
        calcStableC(seg, DP_FOR_HALFSTRIP_INTERSECTION[4*edgeIndex + 2]),
        calcStableC(seg, DP_FOR_HALFSTRIP_INTERSECTION[4*edgeIndex + 3])
      };
    
    const TetraFacet facet =  FACET_FOR_HALFSTRIP_INTERSECTION[edgeIndex];

    
    // special case : facet H = 0
    const bool cond2 = (facet == NO_TET_FACET) ? testSegmentIntersectsHPlane(seg) : testSegmentIntersectsFacet(seg, facet);
    LOG(4, "Halfstrip tests (" << seg << ", " << edge << ") : " << (cVals[0]*cVals[1] < 0.0) << ", " << cond2 << ", " << (cVals[2]*cVals[3] > 0.0) );
    LOG(4, "c2 = " << cVals[2] << ", c3 = " << cVals[3] ); 
  
    return (cVals[0]*cVals[1] < 0.0) && cond2 && (cVals[2]*cVals[3] > 0.0);
  }

  /**
   * Calculates the point of intersection between the given segment of the triangle
   * and the halfstrip above the given edge of the tetrahedron. (Grandy, eq. [31])
   *
   * @pre   testSegmentHalfstripIntersection(seg, edge) returns true
   * 
   * @param seg    segment of the triangle
   * @param edge   edge of the tetrahedron defining the halfstrip
   * @param pt     array of three doubles in which to store the coordinates of the intersection point
   */
  void TransformedTriangle::calcIntersectionPtSegmentHalfstrip(const TriSegment seg, const TetraEdge edge, double* pt) const
  {
    assert(edge > OZ);
    assert(edge < H01);

    // get right index here to avoid "filling out" array
    const int edgeIndex = static_cast<int>(edge) - 3;
    assert(edgeIndex >= 0);
    assert(edgeIndex < 3);
    
    // Barycentric interpolation on the edge
    // for edge AB : (x,y,z)* = (1-alpha) * A + alpha * B
    // where alpha = cB / (cB - cA)

    const double cA = calcStableC(seg, DP_FOR_HALFSTRIP_INTERSECTION[4*edgeIndex]);
    const double cB = calcStableC(seg, DP_FOR_HALFSTRIP_INTERSECTION[4*edgeIndex + 1]);
    assert(cA != cB);
    
    const double alpha = cA / (cA - cB);
    
    for(int i = 0; i < 3; ++i)
      {
        const TetraCorner corners[2] = 
          {
            CORNERS_FOR_EDGE[2*edge],
            CORNERS_FOR_EDGE[2*edge + 1]
          };

        const double cornerCoords[2] = 
          {
            COORDS_TET_CORNER[3*corners[0] + i],
            COORDS_TET_CORNER[3*corners[1] + i]
          };

        pt[i] = (1 - alpha) * cornerCoords[0] + alpha * cornerCoords[1];
        LOG(6, pt[i] );
        assert(pt[i] >= 0.0);
        assert(pt[i] <= 1.0);
      }
    assert(epsilonEqualRelative(pt[0] + pt[1] + pt[2], 1.0));
  }
    
  /**
   * Tests if the given segment of triangle PQR intersects the ray pointing 
   * in the upwards z - direction from the given corner of the tetrahedron. (Grandy eq. [29])
   * If OPTIMIZE is defined, the double product that should be zero is not verified.
   * 
   * @param seg    segment of the triangle PQR
   * @param corner corner of the tetrahedron on the h = 0 facet (X, Y, or Z)
   * @return      true if the upwards ray from the corner intersects the segment. 
   */
  bool TransformedTriangle::testSegmentRayIntersection(const TriSegment seg, const TetraCorner corner) const
  {
    assert(corner == X || corner == Y || corner == Z);
    LOG(4, "Testing seg - ray intersection for seg = " << seg << ", corner = " << corner );

    // readjust index since O is not used
    const int cornerIdx = static_cast<int>(corner) - 1;

    // facets to use
    //? not sure this is correct
    static const TetraFacet FIRST_FACET_SEGMENT_RAY_INTERSECTION[3] = 
      {
        OZX, // X
        OYZ, // Y
        OZX, // Z
      };

    // cond 2
    const bool cond21 = testSegmentIntersectsFacet(seg, FIRST_FACET_SEGMENT_RAY_INTERSECTION[cornerIdx]);
    const bool cond22  = (corner == Z) ? testSegmentIntersectsFacet(seg, OYZ) : testSegmentIntersectsHPlane(seg);
    
    if(!(cond21 || cond22))
      {
        LOG(4, "SR fails at cond 2 : cond21 = " << cond21 << ", cond22 = " << cond22 );
        return false;
      }
    
    // cond 3 
    const double cVals[6] = 
      {
        calcStableC(seg, DP_SEGMENT_RAY_INTERSECTION[7*cornerIdx + 1]),
        calcStableC(seg, DP_SEGMENT_RAY_INTERSECTION[7*cornerIdx + 2]),
        calcStableC(seg, DP_SEGMENT_RAY_INTERSECTION[7*cornerIdx + 3]),
        calcStableC(seg, DP_SEGMENT_RAY_INTERSECTION[7*cornerIdx + 4]),
        calcStableC(seg, DP_SEGMENT_RAY_INTERSECTION[7*cornerIdx + 5]),
        calcStableC(seg, DP_SEGMENT_RAY_INTERSECTION[7*cornerIdx + 6]),
      };
    
    // cond. 3
    if(( (cVals[0] + cVals[1])*(cVals[2] - cVals[3]) - cVals[4]*cVals[5] ) >= 0.0)
      {
        LOG(4, "SR fails at cond 3 : " << (cVals[0] + cVals[1])*(cVals[2] - cVals[3]) - cVals[4]*cVals[5]  );
      }
    return ( (cVals[0] + cVals[1])*(cVals[2] - cVals[3]) - cVals[4]*cVals[5] ) < 0.0;
    
  } 

  // /////////////////////////////////////////////////////////////////////////////////
  //  Utility methods used in intersection tests                       ///////////////
  // /////////////////////////////////////////////////////////////////////////////////
  /**
   * Tests if the triangle PQR surrounds the axis on which the
   * given edge of the tetrahedron lies.
   *
   * @param edge   edge of tetrahedron
   * @return      true if PQR surrounds edge, false if not (see Grandy, eq. [53])
   */
  bool TransformedTriangle::testTriangleSurroundsEdge(const TetraEdge edge) const
  {
    // NB DoubleProduct enum corresponds to TetraEdge enum according to Grandy, table III
    // so we can use the edge directly
    
    const double cPQ = calcStableC(PQ, DoubleProduct(edge));
    const double cQR = calcStableC(QR, DoubleProduct(edge));
    const double cRP = calcStableC(RP, DoubleProduct(edge));

    LOG(5, "TriangleSurroundsEdge : edge = " << edge << " c = [" << cPQ << ", " << cQR << ", " << cRP << "]" );

    // if two or more c-values are zero we disallow x-edge intersection
    // Grandy, p.446
    const int numZeros = (cPQ == 0.0 ? 1 : 0) + (cQR == 0.0 ? 1 : 0) + (cRP == 0.0 ? 1 : 0);
    
    if(numZeros >= 2 ) 
      {
        LOG(5, "TriangleSurroundsEdge test fails due to too many 0 dp" ); 
      }

    return (cPQ*cQR >= 0.0) && (cQR*cRP >= 0.0) && (cRP*cPQ >= 0.0) && numZeros < 2;
  }

} // NAMESPACE
