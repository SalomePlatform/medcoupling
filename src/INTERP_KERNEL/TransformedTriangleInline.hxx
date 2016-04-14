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

#ifndef __TRANSFORMEDTRIANGLEINLINE_HXX__
#define __TRANSFORMEDTRIANGLEINLINE_HXX__

// This file contains inline versions of some of the methods in the TransformedTriangle*.cxx files.
// It replaces those methods if OPTIMIZE is defined.
// NB : most of these methods have documentation in their corresponding .cxx - file.

// ----------------------------------------------------------------------------------
//  Optimization methods. These are only defined and used if OPTIMIZE is defined.
// -----------------------------------------------------------------------------------


inline void TransformedTriangle::preCalculateTriangleSurroundsEdge() 
{
  for(TetraEdge edge = OX ; edge <= ZX ; edge = TetraEdge(edge + 1))
    {
      _triangleSurroundsEdgeCache[edge] = testTriangleSurroundsEdge(edge);
    }
}


// ----------------------------------------------------------------------------------
//   TransformedTriangle_math.cxx                                                 
// ----------------------------------------------------------------------------------

inline void TransformedTriangle::resetDoubleProducts(const TriSegment seg, const TetraCorner corner)
{
  // set the three corresponding double products to 0.0
  static const DoubleProduct DOUBLE_PRODUCTS[12] =
    {
      C_YZ, C_ZX, C_XY, // O
      C_YZ, C_ZH, C_YH, // X
      C_ZX, C_ZH, C_XH, // Y
      C_XY, C_YH, C_XH  // Z
    };
  
  for(int i = 0 ; i < 3 ; ++i) {
    const DoubleProduct dp = DOUBLE_PRODUCTS[3*corner + i];
    
    LOG(6, std::endl << "resetting inconsistent dp :" << dp << " for corner " << corner);
    _doubleProducts[8*seg + dp] = 0.0;
  };
}

inline double TransformedTriangle::calcStableC(const TriSegment seg, const DoubleProduct dp) const
{
  return _doubleProducts[8*seg + dp];
}

inline double TransformedTriangle::calcStableT(const TetraCorner corner) const
{
  //   assert(_isTripleProductsCalculated);
  //   assert(_validTP[corner]);
  return _tripleProducts[corner];
}

inline double TransformedTriangle::calcUnstableC(const TriSegment seg, const DoubleProduct dp) const
{
  
  // find the points of the triangle
  // 0 -> P, 1 -> Q, 2 -> R 
  const int pt1 = seg;
  const int pt2 = (seg + 1) % 3;
  
  // find offsets
  const int off1 = DP_OFFSET_1[dp];
  const int off2 = DP_OFFSET_2[dp];
  
  return _coords[5*pt1 + off1] * _coords[5*pt2 + off2] - _coords[5*pt1 + off2] * _coords[5*pt2 + off1];
}

// ----------------------------------------------------------------------------------
//  TransformedTriangle_intersect.cxx                                            
// ----------------------------------------------------------------------------------
inline bool TransformedTriangle::testSurfaceEdgeIntersection(const TetraEdge edge) const 
{ 
  return _triangleSurroundsEdgeCache[edge] && testEdgeIntersectsTriangle(edge);
}

inline bool TransformedTriangle::testSegmentFacetIntersection(const TriSegment seg, const TetraFacet facet) const 
{ 
  return testFacetSurroundsSegment(seg, facet) && testSegmentIntersectsFacet(seg, facet); 
}

inline bool TransformedTriangle::testSurfaceRayIntersection(const TetraCorner corner) const
{ 
  return testTriangleSurroundsRay( corner ) && testSurfaceAboveCorner( corner ); 
}

inline bool TransformedTriangle::testCornerInTetrahedron(const TriCorner corner) const
{
  const double pt[4] = 
    {
      _coords[5*corner],     // x
      _coords[5*corner + 1], // y
      _coords[5*corner + 2], // z
      _coords[5*corner + 3]  // z
    };
  
  for(int i = 0 ; i < 4 ; ++i) 
    {
      if(pt[i] < 0.0 || pt[i] > 1.0)
        {
          return false;
        }
    }
  return true;
}

inline  bool TransformedTriangle::testCornerOnXYZFacet(const TriCorner corner) const
{
#if 0
  const double pt[4] = 
    {
      _coords[5*corner],     // x
      _coords[5*corner + 1], // y 
      _coords[5*corner + 2], // z
      _coords[5*corner + 3]  // h
    };
#endif
  const double* pt = &_coords[5*corner];
    
  if(pt[3] != 0.0) 
    {
      return false;
    }

  for(int i = 0 ; i < 3 ; ++i) 
    {
      if(pt[i] < 0.0 || pt[i] > 1.0)
        {
          return false;
        }
    }
  return true;
}

inline  bool TransformedTriangle::testCornerAboveXYZFacet(const TriCorner corner) const
{
  const double x = _coords[5*corner];
  const double y = _coords[5*corner + 1];
  const double h = _coords[5*corner + 3];
  const double H = _coords[5*corner + 4];
        
  return h < 0.0 && H >= 0.0 && x >= 0.0 && y >= 0.0;
        
}

inline bool TransformedTriangle::testEdgeIntersectsTriangle(const TetraEdge edge) const
{
  
  //  assert(edge < H01);
  
  // correspondance edge - triple products
  // for edges OX, ..., ZX (Grandy, table III)
  static const TetraCorner TRIPLE_PRODUCTS[12] = 
    {
      X, O, // OX
      Y, O, // OY
      Z, O, // OZ 
      X, Y, // XY
      Y, Z, // YZ
      Z, X, // ZX
    };

  // Grandy, [16]
  const double t1 = calcStableT(TRIPLE_PRODUCTS[2*edge]);
  const double t2 = calcStableT(TRIPLE_PRODUCTS[2*edge + 1]);

  //? should equality with zero use epsilon?
  LOG(5, "testEdgeIntersectsTriangle : t1 = " << t1 << " t2 = " << t2 );
  return (t1*t2 <= 0.0) && (t1 - t2 != 0.0);
}

inline bool TransformedTriangle::testFacetSurroundsSegment(const TriSegment seg, const TetraFacet facet) const
{
#if 0
  const double signs[3] = 
    {
      SIGN_FOR_SEG_FACET_INTERSECTION[3*facet],
      SIGN_FOR_SEG_FACET_INTERSECTION[3*facet + 1],
      SIGN_FOR_SEG_FACET_INTERSECTION[3*facet + 2]
    };
#endif

  const double* signs = &SIGN_FOR_SEG_FACET_INTERSECTION[3*facet];
  const double c1 = signs[0]*calcStableC(seg, DP_FOR_SEG_FACET_INTERSECTION[3*facet]);
  const double c2 = signs[1]*calcStableC(seg, DP_FOR_SEG_FACET_INTERSECTION[3*facet + 1]);
  const double c3 = signs[2]*calcStableC(seg, DP_FOR_SEG_FACET_INTERSECTION[3*facet + 2]);

  return (c1*c3 > 0.0) && (c2*c3 > 0.0);
}

inline bool TransformedTriangle::testSegmentIntersectsFacet(const TriSegment seg, const TetraFacet facet) const
{
  // use correspondance facet a = 0 <=> offset for coordinate a in _coords
  // and also correspondance segment AB => corner A
  const double coord1 = _coords[5*seg + facet];
  const double coord2 = _coords[5*( (seg + 1) % 3) + facet];
  
  //? should we use epsilon-equality here in second test?
  LOG(5, "coord1 : " << coord1 << " coord2 : " << coord2 );
  
  return (coord1*coord2 <= 0.0) && (coord1 != coord2);
}

inline bool TransformedTriangle::testSegmentIntersectsHPlane(const TriSegment seg) const
{
  // get the H - coordinates
  const double coord1 = _coords[5*seg + 4];
  const double coord2 = _coords[5*( (seg + 1) % 3) + 4];
  //? should we use epsilon-equality here in second test?
  LOG(5, "coord1 : " << coord1 << " coord2 : " << coord2 );
  
  return (coord1*coord2 <= 0.0) && (coord1 != coord2);
}

inline bool TransformedTriangle::testSurfaceAboveCorner(const TetraCorner corner) const
{
  // ? There seems to be an error in Grandy -> it should be C_XY instead of C_YZ in [28].
  // ? I haven't really figured out why, but it seems to work.
  const double normal = calcStableC(PQ, C_XY) + calcStableC(QR, C_XY) + calcStableC(RP, C_XY);

  LOG(6, "surface above corner " << corner << " : " << "n = " << normal << ", t = [" <<  calcTByDevelopingRow(corner, 1, false) << ", "  << calcTByDevelopingRow(corner, 2, false) << ", " << calcTByDevelopingRow(corner, 3, false) );
  LOG(6, "] - stable : " << calcStableT(corner)  );

  //? we don't care here if the triple product is "invalid", that is, the triangle does not surround one of the
  // edges going out from the corner (Grandy [53])
  if(!_validTP[corner])
    {
      return ( calcTByDevelopingRow(corner, 1, false) * normal ) >= 0.0;
    }
  else
    {
      return ( calcStableT(corner) * normal ) >= 0.0;
    }
}

inline bool TransformedTriangle::testTriangleSurroundsRay(const TetraCorner corner) const
{
  //  assert(corner == X || corner == Y || corner == Z);

  // double products to use for the possible corners
  static const DoubleProduct DP_FOR_RAY_INTERSECTION[4] = 
    {
      DoubleProduct(0),        // O - only here to fill out and make indices match
      C_10,     // X
      C_01,     // Y
      C_XY      // Z
    };

  const DoubleProduct dp = DP_FOR_RAY_INTERSECTION[corner];

  const double cPQ = calcStableC(PQ, dp);
  const double cQR = calcStableC(QR, dp);
  const double cRP = calcStableC(RP, dp);

  //? NB here we have no correction for precision - is this good?
  // Our authority Grandy says nothing
  LOG(5, "dp in triSurrRay for corner " << corner << " = [" << cPQ << ", " << cQR << ", " << cRP << "]" );

  return ( cPQ*cQR > 0.0 ) && ( cPQ*cRP > 0.0 );

}
#endif
