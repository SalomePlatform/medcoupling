// Copyright (C) 2007-2024  CEA, EDF
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
#include <limits>
#include <map>
#include <utility>

#include "VectorUtils.hxx"

namespace INTERP_KERNEL
{
  
  // ----------------------------------------------------------------------------------
  //  Tables                                                       
  // ----------------------------------------------------------------------------------

  /// Table with first coordinate (a) used to calculate double product c^pq_ab = p_a * q_b - p_b * q_a (index to be used : DoubleProduct)
  const int TransformedTriangle::DP_OFFSET_1[8] = {1, 2, 0, 2, 0, 1, 4, 1};

  /// Table with second coordinate (b) used to calculate double product c^pq_ab = p_a * q_b - p_b * q_a (index to be used : DoubleProduct)
  const int TransformedTriangle::DP_OFFSET_2[8] = {2, 0, 1, 3, 3, 3, 0, 4};

  /// Coordinates used to calculate triple products by the expanding one of the three rows of the determinant (index to be used : 3*Corner + row)
  const int TransformedTriangle::COORDINATE_FOR_DETERMINANT_EXPANSION[12] =
    {
      // row 1, 2, 3
      0, 1, 2, // O
      3, 1, 2, // X
      0, 3, 2, // Y
      0, 1, 3  // Z
    };
  
  /// Double products used to calculate triple products by expanding one of the three rows of the determinant (index to be used : 3*Corner + row)
  const TransformedTriangle::DoubleProduct TransformedTriangle::DP_FOR_DETERMINANT_EXPANSION[12] = 
    {
      // row 1, 2, 3
      C_YZ, C_ZX, C_XY, // O
      C_YZ, C_ZH, C_YH, // X
      C_ZH, C_ZX, C_XH, // Y
      C_YH, C_XH, C_XY  // Z
    };
  
  /// The machine epsilon, used in precision corrections
  const double TransformedTriangle::MACH_EPS = std::numeric_limits<double>::epsilon();
  
  /// 4.0 * the machine epsilon, represents the precision of multiplication when performing corrections corrections ( f in Grandy )
  const double TransformedTriangle::MULT_PREC_F = 4.0 * TransformedTriangle::MACH_EPS;

  /// Threshold for resetting double and triple products to zero; ( F / f in Grandy )
  const double TransformedTriangle::THRESHOLD_F = 100.0;

  /// Threshold for what is considered a small enough angle to warrant correction of triple products by Grandy, [57]
  const double TransformedTriangle::TRIPLE_PRODUCT_ANGLE_THRESHOLD = 0.1;


  // Handle cases where one of the segment (or all) is (almost) in XYZ plane.
  // We follow Grandy's suggestion and perturb slightly to have exactly h=0 for the segment (Grandy p.447)
  // Note that if PQR is == to the upper facet of the unit tetra (XYZ), the tetra-corner-inclusion test should take it in,
  // thanks to Grandy [21] and the fact that S_x test is "<=0" (not <0)
  // After that, we also snap P,Q,R to the corners if they're very close.
  void TransformedTriangle::handleDegenerateCases()
  {
    static const TriCorner PT_SEG_MAP[] = {
      P, Q,
      Q, R,
      R, P
    };

    const double eps = THRESHOLD_F*TransformedTriangle::MULT_PREC_F;
    for (TriSegment seg = PQ; seg <= RP; seg = TriSegment(seg+1))
      {
        // Is h coordinate for both end of segment small enough?
        int pt1 = PT_SEG_MAP[2*seg], pt2 = PT_SEG_MAP[2*seg+1];
        if (fabs(_coords[5*pt1+3]) < eps && fabs(_coords[5*pt2+3]) < eps)
          {
            // If so, perturb x,y and z to reset h to exactly zero.
            for (auto pt: {pt1, pt2})  // thx C++17
              {
                const double correc = _coords[pt*5+3]/3.; // this should be really small!
                _coords[pt*5+0] += correc;
                _coords[pt*5+1] += correc;
                _coords[pt*5+2] += correc;
                // And then, if x,y or z very close to 0 or 1, snap exactly to tetra corner:
                for(int d=0; d < 3; d++)
                  {
                    if (fabs(_coords[5*pt+d]) < eps)    _coords[5*pt+d] = 0.0;
                    if (fabs(_coords[5*pt+d]-1) < eps)  _coords[5*pt+d] = 1.0;
                  }
                _coords[pt*5+3] = 0.0;
              }
          }
      }
  }
  
  // ----------------------------------------------------------------------------------
  //  Double and triple product calculations                           
  // ----------------------------------------------------------------------------------
  
  /**
   * Pre-calculates all double products for this triangle, and stores
   * them internally. This method makes compensation for precision errors,
   * and it is thus the "stable" double products that are stored.
   *
   */
  void TransformedTriangle::preCalculateDoubleProducts()
  {
    if(_is_double_products_calculated)
      return;

    // -- calculate all unstable double products -- store in _doubleProducts
    for(TriSegment seg = PQ ; seg <= RP ; seg = TriSegment(seg + 1))
      {
        for(DoubleProduct dp = C_YZ ; dp <= C_10 ; dp = DoubleProduct(dp + 1))
          {
            const int idx = 8*seg + dp;
            _doubleProducts[idx] = calcUnstableC(seg, dp, _deltas[idx]);
          }
      }

    std::map<double, TetraCorner> distances;

    // -- (1) for each segment : check that double products satisfy Grandy, [46]
    // -- and make corrections if not
    for(TriSegment seg = PQ ; seg <= RP ; seg = TriSegment(seg + 1))
      {
        if(!areDoubleProductsConsistent(seg))
          {
            LOG(4, "inconsistent! ");
            for(TetraCorner corner = O ; corner <= Z ; corner = TetraCorner(corner + 1))
              {
                // calculate distance corner - segment axis
                const double dist = calculateDistanceCornerSegment(corner, seg);
                distances.insert( std::make_pair( dist, corner ) );
              }

            // first element -> minimum distance (because map is sorted)
            const TetraCorner minCorner = distances.begin()->second;
            resetDoubleProducts(seg, minCorner);
            distances.clear();
          }
      }

    // -- (2) check that each double product satisfies Grandy, [47], else set to 0
    for(TriSegment seg = PQ ; seg <= RP ; seg = TriSegment(seg + 1))
      {
        for(DoubleProduct dp = C_YZ ; dp <=  C_10 ; dp = DoubleProduct(dp + 1))
          {
            const int idx = 8*seg+dp;

            if( epsilonEqual(_doubleProducts[idx], 0.0, THRESHOLD_F * MULT_PREC_F * _deltas[idx]))
              {
                // debug output
#if LOG_LEVEL >= 5
                if(_doubleProducts[8*seg + dp] != 0.0)
                  {
                    LOG(5, "Double product for (seg,dp) = (" << seg << ", " << dp << ") = " );
                    LOG(5, std::fabs(_doubleProducts[8*seg + dp]) << " is imprecise, reset to 0.0" );
                  }
#endif 

                _doubleProducts[idx] = 0.0;
              }
          }
      }

    _is_double_products_calculated = true;
  }

  /**
   * Checks if the double products for a given segment are consistent, as defined by
   * Grandy, [46]. 
   *
   * @param   seg Segment for which to check consistency of double products
   * @return  true if the double products are consistent, false if not
   */
  bool TransformedTriangle::areDoubleProductsConsistent(const TriSegment seg) const
  {
    // Careful! Here doubleProducts have not yet been corrected for roundoff errors!
    // So we need to epsilon-adjust to correctly identify zeros:
    static const DoubleProduct DP_LST[6] = {C_YZ, C_XH,
                                            C_ZX, C_YH,
                                            C_XY, C_ZH};
    double dps[6];
    for (int i = 0; i < 6; i++)
      {
        const double dp = _doubleProducts[8*seg + DP_LST[i]];
        dps[i] = dp;
      }

    const double term1 = dps[0] * dps[1];
    const double term2 = dps[2] * dps[3];
    const double term3 = dps[4] * dps[5];

    LOG(2, "for seg " << seg << " consistency " << term1 + term2 + term3 );
    LOG(2, "term1 :" << term1 << " term2 :" << term2 << " term3: " << term3 );

    // Test for "== 0.0" here is OK since doubleProduct has been fixed for rounding to zero already.
    const int num_zero = (term1 == 0.0 ? 1 : 0) + (term2 == 0.0 ? 1 : 0) + (term3 == 0.0 ? 1 : 0);
    const int num_neg = (term1 < 0.0 ? 1 : 0) + (term2 < 0.0 ? 1 : 0) + (term3 < 0.0 ? 1 : 0);
    const int num_pos = (term1 > 0.0 ? 1 : 0) + (term2 > 0.0 ? 1 : 0) + (term3 > 0.0 ? 1 : 0);

    assert( num_zero + num_neg + num_pos == 3 );

    // Calculated geometry is inconsistent if we have one of the following cases
    // * one term zero and the other two of the same sign
    // * two terms zero
    // * all terms positive
    // * all terms negative
    const bool inconsist = (num_zero == 1 && num_neg != 1) ||
                           num_zero == 2 ||
                           (num_neg == 0 && num_zero != 3) ||
                           num_neg == 3;
    if(inconsist)  {
        LOG(4, "inconsistent dp found" );
      }
    return !inconsist;
  }

  /**
   * Calculate the shortest distance between a tetrahedron corner and a triangle segment.
   * 
   * @param  corner corner of the tetrahedron
   * @param  seg    segment of the triangle
   * @return shortest distance from the corner to the segment
   */
  double TransformedTriangle::calculateDistanceCornerSegment(const TetraCorner corner, const TriSegment seg) const
  {
    // NB uses fact that TriSegment <=> TriCorner that is first point of segment (PQ <=> P)
    const TriCorner ptP_idx = TriCorner(seg);
    const TriCorner ptQ_idx = TriCorner( (seg + 1) % 3);
    
    const double ptP[3] = { _coords[5*ptP_idx], _coords[5*ptP_idx + 1], _coords[5*ptP_idx + 2]  };
    const double ptQ[3] = { _coords[5*ptQ_idx], _coords[5*ptQ_idx + 1], _coords[5*ptQ_idx + 2]  };
    
    // coordinates of corner
    const double ptTetCorner[3] = 
      { 
        COORDS_TET_CORNER[3*corner    ],
        COORDS_TET_CORNER[3*corner + 1],
        COORDS_TET_CORNER[3*corner + 2]
      };
    
    // dist^2 = ( PQ x CP )^2 / |PQ|^2 where C is the corner point
    
    // difference vectors
    const double diffPQ[3] = { ptQ[0] - ptP[0], ptQ[1] - ptP[1], ptQ[2] - ptP[2] };
    const double diffCornerP[3] = { ptP[0] - ptTetCorner[0], ptP[1] - ptTetCorner[1], ptP[2] - ptTetCorner[2] };
    
    // cross product of difference vectors
    double crossProd[3];
    cross(diffPQ, diffCornerP, crossProd);
    
    const double cross_squared = dot(crossProd, crossProd);
    const double norm_diffPQ_squared = dot(diffPQ, diffPQ);
    
    assert(norm_diffPQ_squared != 0.0);
    
    return cross_squared / norm_diffPQ_squared;
  }

  /**
   * Pre-calculates all triple products for the tetrahedron with respect to
   * this triangle, and stores them internally. This method takes into account
   * the problem of errors due to cancellation.
   *
   */
  void TransformedTriangle::preCalculateTripleProducts()
  {
    if(_is_triple_products_calculated)
      return;

    // find edge / row to use -> that whose edge makes the smallest angle to the triangle
    // use a map to find the minimum
    std::map<double, int> anglesForRows;

    LOG(4, "Precalculating triple products" );
    for(TetraCorner corner = O ; corner <= Z ; corner = TetraCorner(corner + 1))
      {
        LOG(6, "- Triple product for corner " << corner );

        for(int row = 1 ; row < 4 ; ++row) 
          {
            const DoubleProduct dp = DP_FOR_DETERMINANT_EXPANSION[3*corner + (row - 1)];

            // get edge by using correspondence between Double Product and Edge
            TetraEdge edge = TetraEdge(dp);

            // use edge only if it is surrounded by the surface
            if( _triangleSurroundsEdgeCache[edge] )
                {
                  // -- calculate angle between edge and PQR
                  const double angle = calculateAngleEdgeTriangle(edge);
                  anglesForRows.insert(std::make_pair(angle, row));
                }
          }

        if(anglesForRows.size() != 0) // we have found a good row
          {
            const double minAngle = anglesForRows.begin()->first;
            const int minRow = anglesForRows.begin()->second;

            if(minAngle < TRIPLE_PRODUCT_ANGLE_THRESHOLD)
              _tripleProducts[corner] = calcTByDevelopingRow(corner, minRow, true);
            else 
              _tripleProducts[corner] = calcTByDevelopingRow(corner, minRow, false);

            _validTP[corner] = true;
          }
        else
          {
            // this value will not be used - we set it to whatever
            LOG(6, "Triple product not calculated for corner " << corner );
            _tripleProducts[corner] = std::nan("triplep");
            _validTP[corner] = false;
          }
        anglesForRows.clear();
      }
    _is_triple_products_calculated = true;
  }

  /**
   * Calculates the angle between an edge of the tetrahedron and the triangle
   *
   * @param  edge edge of the tetrahedron
   * @return angle between triangle and edge
   */
  double TransformedTriangle::calculateAngleEdgeTriangle(const TetraEdge edge) const
  {
    // find normal to PQR - cross PQ and PR
    const double pq[3] = 
      { 
        _coords[5*Q]     - _coords[5*P], 
        _coords[5*Q + 1] - _coords[5*P + 1],
        _coords[5*Q + 2] - _coords[5*P + 2]
      };
    
    const double pr[3] = 
      { 
        _coords[5*R]     - _coords[5*P], 
        _coords[5*R + 1] - _coords[5*P + 1],
        _coords[5*R + 2] - _coords[5*P + 2]
      };
    
    double normal[3];

    cross(pq, pr, normal);
    
    static const double EDGE_VECTORS[18] =
      {
        1.0, 0.0, 0.0, // OX
        0.0, 1.0, 0.0, // OY
        0.0, 0.0, 1.0, // OZ
        -1.0, 1.0, 0.0, // XY
        0.0,-1.0, 1.0, // YZ
        1.0, 0.0,-1.0 // ZX
      };
    
    const double edgeVec[3] = { 
      EDGE_VECTORS[3*edge],
      EDGE_VECTORS[3*edge + 1],
      EDGE_VECTORS[3*edge + 2],
    };

    //return angleBetweenVectors(normal, edgeVec);

    const double lenNormal = norm(normal);
    const double lenEdgeVec = norm(edgeVec);
    const double dotProd = dot(normal, edgeVec);
    
    //? is this more stable? -> no subtraction
    //    return asin( dotProd / ( lenNormal * lenEdgeVec ) ) + 3.141592625358979 / 2.0;
    const double tmp=dotProd / ( lenNormal * lenEdgeVec );
    const double safe_tmp=std::max(std::min(tmp,1.),-1.);
    return M_PI - std::acos(safe_tmp);
  }

  /**
   * Calculates triple product associated with the given corner of tetrahedron, developing
   * the determinant by the given row. The triple product gives the signed volume of 
   * the tetrahedron between this corner and the triangle PQR. If the flag project is true, 
   * one coordinate is projected out in order to eliminate errors in the intersection point
   * calculation due to cancellation.
   * 
   * Consistency with the double product computation and potential cancellation is also done here.
   *
   *
   * @pre            double products have already been calculated
   * @param corner   corner for which the triple product is calculated
   * @param row      row (1 <= row <= 3) used to calculate the determinant
   * @param project  indicates whether or not to perform projection as inidicated in Grandy, p.446
   * @return        triple product associated with corner (see Grandy, [50]-[52])
   */
  double TransformedTriangle::calcTByDevelopingRow(const TetraCorner corner, const int row, const bool project) const
  {
    
    // OVERVIEW OF CALCULATION
    // --- sign before the determinant
    // the sign used depends on the sign in front of the triple product (Grandy, [15]),
    // and the convention used in the definition of the double products
  
    // the sign in front of the determinant gives the following schema for the three terms (I): 
    // corner/row    1    2   3
    // O (sign:+)    +    -   +
    // X (sign:-)    -    +   -
    // Y (sign:-)    -    +   -
    // Z (sign:-)    -    +   -

    // the 2x2 determinants are the following (C_AB <=> A_p*B_q - B_p*A_q, etc)
    // corner/row    1       2     3
    // O (sign:+)   C_YZ   C_XZ  C_XY
    // X (sign:-)   C_YZ   C_HZ  C_HY
    // Y (sign:-)   C_HZ   C_XZ  C_XH
    // Z (sign:-)   C_YH   C_XH  C_XY

    // these are represented in DP_FOR_DETERMINANT_EXPANSION,
    // except for the fact that certain double products are inversed (C_AB <-> C_BA)

    // comparing with the DOUBLE_PRODUCTS and using the fact that C_AB = -C_BA
    // we deduce the following schema (II) :
    // corner/row    1    2   3
    // O (sign:+)    +    -   +
    // X (sign:-)    +    -   -
    // Y (sign:-)    -    -   +
    // Z (sign:-)    +    +   +

    // comparing the two schemas (I) and (II) gives us the following matrix of signs,
    // putting 1 when the signs in (I) and (II) are equal and -1 when they are different :

    static const int SIGNS[12] = 
      {
        1, 1, 1,
        -1,-1, 1,
        1,-1,-1,
        -1, 1,-1
      };

    // find the offsets of the rows of the determinant
    const int offset = COORDINATE_FOR_DETERMINANT_EXPANSION[3 * corner + (row - 1)];
  
    const DoubleProduct dp = DP_FOR_DETERMINANT_EXPANSION[3 * corner + (row - 1)];

    const int sign = SIGNS[3 * corner + (row - 1)];

    const double cQR = calcStableC(QR, dp);
    const double cRP = calcStableC(RP, dp);
    const double cPQ = calcStableC(PQ, dp);

    double alpha = 0.0;

    // coordinate to use for projection (Grandy, [57]) with edges
    // OX, OY, OZ, XY, YZ, ZX in order : 
    // (y, z, x, h, h, h)
    // for the first three we could also use {2, 0, 1}
    static const int PROJECTION_COORDS[6] = { 1, 2, 0, 3, 3, 3 } ;

    const int coord = PROJECTION_COORDS[ dp ];

    // coordinate values for P, Q and R
    const double coordValues[3] = { _coords[5*P + coord], _coords[5*Q + coord], _coords[5*R + coord] };

    if(project)
      {
        // products coordinate values with corresponding double product
        const double coordDPProd[3] = { coordValues[0] * cQR, coordValues[1] * cRP, coordValues[2] * cPQ };

        const double sumDPProd = coordDPProd[0] + coordDPProd[1] + coordDPProd[2];
        const double sumDPProdSq = dot(coordDPProd, coordDPProd);

        //       alpha = sumDPProd / sumDPProdSq;
        alpha = (sumDPProdSq != 0.0) ? sumDPProd / sumDPProdSq : 0.0;
      }

    const double cQRbar = cQR * (1.0 - alpha * coordValues[0] * cQR);
    const double cRPbar = cRP * (1.0 - alpha * coordValues[1] * cRP);
    const double cPQbar = cPQ * (1.0 - alpha * coordValues[2] * cPQ);

    // [ABN] Triple product cancellation logic:
    // This part is not well described in Grandy (end of p.446) :
    //     "We use a method analogous to (47) to remove imprecise triple products,..."
    //
    // Our algo for cancelling a triple product:
    //   - retrieve the deltas associated with each DP involved (because a DP itself is a sum of two terms - see [42]
    //   - multiply them by the coordinate coming from the determinant expansion
    //   - and finally sum the 3 corresponding terms of the developement
    //
    // Using directly the DP (as was done before here) leads to issues, since this DP might have been cancelled
    // already earlier on, and we lost the delta information -> doing this, we risk not cancelling the triple prod
    // when we should have.
    const double cQRbar_delta = _deltas[8*QR + dp],
                 cRPbar_delta = _deltas[8*RP + dp],
                 cPQbar_delta = _deltas[8*PQ + dp];

    // check sign of barred products - should not change
    //    assert(cQRbar * cQR >= 0.0);
    //assert(cRPbar * cRP >= 0.0);
    //assert(cPQbar * cPQ >= 0.0);

    const double p_term = _coords[5*P + offset] * cQRbar;
    const double q_term = _coords[5*Q + offset] * cRPbar;
    const double r_term = _coords[5*R + offset] * cPQbar;

    const double p_delta = std::fabs(_coords[5*P + offset] * cQRbar_delta),
                 q_delta = std::fabs(_coords[5*Q + offset] * cRPbar_delta),
                 r_delta = std::fabs(_coords[5*R + offset] * cPQbar_delta);

    const double delta = p_delta + q_delta + r_delta;

    if( epsilonEqual( p_term + q_term + r_term, 0.0, THRESHOLD_F * MULT_PREC_F * delta) )
      {
        LOG(4, "Reset imprecise triple product for corner " << corner << " to zero" ); 
        return 0.0;
      }
    else
      {
        // NB : using plus also for the middle term compensates for a double product
        // which is inversely ordered
        LOG(6, "Triple product for corner " << corner << ", row " << row << " = " << sign*( p_term + q_term + r_term ) );
        return sign*( p_term + q_term + r_term );
      }

  }

}
