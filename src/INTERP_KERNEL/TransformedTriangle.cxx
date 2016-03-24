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
#include "VectorUtils.hxx"
#include "TetraAffineTransform.hxx"
#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <functional>
#include <iterator>
#include <math.h>
#include <vector>

namespace INTERP_KERNEL
{

  /**
   * \brief Class representing a circular order of a set of points around their barycenter.
   * It is used with the STL sort() algorithm to sort the point of the two polygons
   *
   */
  class ProjectedCentralCircularSortOrder
  {
  public:

    /// Enumeration of different planes to project on when calculating order
    enum CoordType { XY, XZ, YZ };
  
    /**
     * Constructor
     *
     * @param barycenter  double[3] containing the barycenter of the points to be compared
     * @param type        plane to project on when comparing. The comparison will not work if all the points are in a plane perpendicular
     *                    to the plane being projected on
     */
    ProjectedCentralCircularSortOrder(const double* barycenter, const CoordType type)
      : _aIdx((type == YZ) ? 1 : 0), 
        _bIdx((type == XY) ? 1 : 2),
        _a(barycenter[_aIdx]), 
        _b(barycenter[_bIdx])
    {
    }

    /**
     * Comparison operator.
     * Compares the relative position between two points in their ordering around the barycenter.
     *
     * @param  pt1   a double[3] representing a point
     * @param  pt2   a double[3] representing a point
     * @return       true if the angle of the difference vector between pt1 and the barycenter is greater than that 
     *               of the difference vector between pt2 and the barycenter.
     */
    bool operator()(const double* pt1, const double* pt2)
    {
      // calculate angles with the axis
      const double ang1 = atan2(pt1[_aIdx] - _a, pt1[_bIdx] - _b);
      const double ang2 = atan2(pt2[_aIdx] - _a, pt2[_bIdx] - _b);

      return ang1 > ang2;
    }

  private:
    /// index corresponding to first coordinate of plane on which points are projected
    const int _aIdx;
  
    /// index corresponding to second coordinate of plane on which points are projected
    const int _bIdx;

    /// value of first projected coordinate of the barycenter
    const double _a;
  
    /// value of second projected coordinate of the barycenter
    const double _b;
  };

  // ----------------------------------------------------------------------------------
  // TransformedTriangle PUBLIC  
  // ----------------------------------------------------------------------------------
  
  /**
   * Constructor
   *
   * The coordinates are copied to the internal member variables
   *
   * @param p   array of three doubles containing coordinates of P
   * @param q   array of three doubles containing coordinates of Q
   * @param r   array of three doubles containing coordinates of R
   */
  TransformedTriangle::TransformedTriangle(double* p, double* q, double* r)
    : _is_double_products_calculated(false),  _is_triple_products_calculated(false), _volume(0)
  {
  
    for(int i = 0 ; i < 3 ; ++i) 
      {
        // xyz coordinates
        _coords[5*P + i] = p[i];
        _coords[5*Q + i] = q[i];
        _coords[5*R + i] = r[i];
      }

    // h coordinate
    
    _coords[5*P + 3] = 1 - p[0] - p[1] - p[2];
    _coords[5*Q + 3] = 1 - q[0] - q[1] - q[2];
    _coords[5*R + 3] = 1 - r[0] - r[1] - r[2];

    // H coordinate
    _coords[5*P + 4] = 1 - p[0] - p[1];
    _coords[5*Q + 4] = 1 - q[0] - q[1];
    _coords[5*R + 4] = 1 - r[0] - r[1];

    resetNearZeroCoordinates();

    // initialise rest of data
    preCalculateDoubleProducts();

    preCalculateTriangleSurroundsEdge();

    preCalculateTripleProducts();
 
  }

  /**
   * Destructor
   *
   * Deallocates the memory used to store the points of the polygons.
   * This memory is allocated in calculateIntersectionAndProjectionPolygons().
   */
  TransformedTriangle::~TransformedTriangle()
  {
    // delete elements of polygons
    for(std::vector<double*>::iterator it = _polygonA.begin() ; it != _polygonA.end() ; ++it)
      {
        delete[] *it;
      }
    for(std::vector<double*>::iterator it = _polygonB.begin() ; it != _polygonB.end() ; ++it)
      {
        delete[] *it;
      }    
  }

  /**
   * Calculates the volume of intersection between the triangle and the 
   * unit tetrahedron.
   *
   * @return   volume of intersection of this triangle with unit tetrahedron, 
   *            as described in Grandy
   *
   */
  double TransformedTriangle::calculateIntersectionVolume()
  {
    // check first that we are not below z - plane    
    if(isTriangleBelowTetraeder())
      {
        LOG(2, " --- Triangle is below tetraeder - V = 0.0");
        return 0.0;
      }
    
    // get the sign of the volume -  equal to the sign of the z-component of the normal
    // of the triangle, u_x * v_y - u_y * v_x, where u = q - p and v = r - p
    // if it is zero, the triangle is perpendicular to the z - plane and so the volume is zero
//     const double uv_xy[4] = 
//       {
//         _coords[5*Q] - _coords[5*P], _coords[5*Q + 1] - _coords[5*P + 1], // u_x, u_y
//         _coords[5*R] - _coords[5*P], _coords[5*R + 1] - _coords[5*P + 1]  // v_x, v_y
//       };

//     double sign = uv_xy[0] * uv_xy[3] - uv_xy[1] * uv_xy[2];
    int sign = isTriangleInclinedToFacet( OXY );

    if(sign == 0 )
      {
        LOG(2, " --- Triangle is perpendicular to z-plane - V = 0.0");
        return _volume = 0.0;
      }


    // normalize sign
    //sign = sign > 0.0 ? 1.0 : -1.0;

    LOG(2, "-- Calculating intersection polygons ... ");
    calculateIntersectionAndProjectionPolygons();
    
    double barycenter[3];

    // calculate volume under A
    double volA = 0.0;
    if(_polygonA.size() > 2)
      {
        LOG(2, "---- Treating polygon A ... ");
        calculatePolygonBarycenter(A, barycenter);
        sortIntersectionPolygon(A, barycenter);
        volA = calculateVolumeUnderPolygon(A, barycenter);
        LOG(2, "Volume is " << sign * volA);
      }

    double volB = 0.0;
    // if triangle is not in h = 0 plane, calculate volume under B
    if(_polygonB.size() > 2 && !isTriangleInPlaneOfFacet(XYZ))
      {
        LOG(2, "---- Treating polygon B ... ");
       
        calculatePolygonBarycenter(B, barycenter);
        sortIntersectionPolygon(B, barycenter);
        volB = calculateVolumeUnderPolygon(B, barycenter);
        LOG(2, "Volume is " << sign * volB);
      }

    LOG(2, "volA + volB = " << sign * (volA + volB) << std::endl << "***********");
  
    return _volume = sign * (volA + volB);

  } 

  /**
   * Calculates the volume of intersection between the triangle and the
   * unit tetrahedron.
   *
   * @return   volume of intersection of this triangle with unit tetrahedron,
   *            as described in Grandy
   *
   */
  double TransformedTriangle::calculateIntersectionSurface(TetraAffineTransform* tat)
  {
    // check first that we are not below z - plane
    if(isTriangleBelowTetraeder())
      {
        LOG(2, " --- Triangle is below tetraeder - V = 0.0");
        return 0.0;
      }

    LOG(2, "-- Calculating intersection polygon ... ");
    calculateIntersectionPolygon();

    _volume = 0.;
    if(_polygonA.size() > 2) {
      double barycenter[3];
      calculatePolygonBarycenter(A, barycenter);
      sortIntersectionPolygon(A, barycenter);
      const std::size_t nbPoints = _polygonA.size();
      for(std::size_t i = 0 ; i < nbPoints ; ++i)
        tat->reverseApply(_polygonA[i], _polygonA[i]);
      _volume = calculateSurfacePolygon();
    }

    return _volume;
  }

  // ----------------------------------------------------------------------------------
  // TransformedTriangle PRIVATE
  // ----------------------------------------------------------------------------------

  /**
   * Calculates the intersection polygons A and B, performing the intersection tests
   * and storing the corresponding points in the vectors _polygonA and _polygonB.
   *
   * @post _polygonA contains the intersection polygon A and _polygonB contains the
   *       intersection polygon B.
   *
   */
  void TransformedTriangle::calculateIntersectionAndProjectionPolygons()
  {
    assert(_polygonA.size() == 0);
    assert(_polygonB.size() == 0);
    // avoid reallocations in push_back() by pre-allocating enough memory
    // we should never have more than 20 points
    _polygonA.reserve(20);
    _polygonB.reserve(20);
    // -- surface intersections
    // surface - edge
    for(TetraEdge edge = OX ; edge <= ZX ; edge = TetraEdge(edge + 1))
      {
        if(testSurfaceEdgeIntersection(edge))
          {
            double* ptA = new double[3];
            calcIntersectionPtSurfaceEdge(edge, ptA);
            _polygonA.push_back(ptA);
            LOG(3,"Surface-edge : " << vToStr(ptA) << " added to A ");
            if(edge >= XY)
              {
                double* ptB = new double[3];
                copyVector3(ptA, ptB);
                _polygonB.push_back(ptB);
                LOG(3,"Surface-edge : " << vToStr(ptB) << " added to B ");
              }
           
          }
      }

    // surface - ray
    for(TetraCorner corner = X ; corner < NO_TET_CORNER ; corner = TetraCorner(corner + 1))
      {
        if(testSurfaceRayIntersection(corner))
          {
            double* ptB = new double[3];
            copyVector3(&COORDS_TET_CORNER[3 * corner], ptB);
            _polygonB.push_back(ptB);
            LOG(3,"Surface-ray : " << vToStr(ptB) << " added to B");
          }
      }

    // -- segment intersections
    for(TriSegment seg = PQ ; seg < NO_TRI_SEGMENT ; seg = TriSegment(seg + 1))
      {

        bool isZero[NO_DP];

        // check beforehand which double-products are zero
        for(DoubleProduct dp = C_YZ; dp < NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[dp] = (calcStableC(seg, dp) == 0.0);
          }

        // segment - facet
        for(TetraFacet facet = OYZ ; facet < NO_TET_FACET ; facet = TetraFacet(facet + 1))
          {
            // is this test worth it?
            const bool doTest = 
              !isZero[DP_FOR_SEG_FACET_INTERSECTION[3*facet]] && 
              !isZero[DP_FOR_SEG_FACET_INTERSECTION[3*facet + 1]] &&
              !isZero[DP_FOR_SEG_FACET_INTERSECTION[3*facet + 2]];

            if(doTest && testSegmentFacetIntersection(seg, facet))
              {
                double* ptA = new double[3];
                calcIntersectionPtSegmentFacet(seg, facet, ptA);
                _polygonA.push_back(ptA);
                LOG(3,"Segment-facet : " << vToStr(ptA) << " added to A");
                if(facet == XYZ)
                  {
                    double* ptB = new double[3];
                    copyVector3(ptA, ptB);
                    _polygonB.push_back(ptB);
                    LOG(3,"Segment-facet : " << vToStr(ptB) << " added to B");
                  }

              }
          }

        // segment - edge
        for(TetraEdge edge = OX ; edge <= ZX ; edge = TetraEdge(edge + 1))
          {
            const DoubleProduct edge_dp = DoubleProduct(edge);

            if(isZero[edge_dp] && testSegmentEdgeIntersection(seg, edge))
              {
                double* ptA = new double[3];
                calcIntersectionPtSegmentEdge(seg, edge, ptA);
                _polygonA.push_back(ptA);
                LOG(3,"Segment-edge : " << vToStr(ptA) << " added to A");
                if(edge >= XY)
                  {
                    double* ptB = new double[3];
                    copyVector3(ptA, ptB);
                    _polygonB.push_back(ptB);
                  }
              }
          }
       
        // segment - corner
        for(TetraCorner corner = O ; corner < NO_TET_CORNER ; corner = TetraCorner(corner + 1))
          {
            const bool doTest = 
              isZero[DoubleProduct( EDGES_FOR_CORNER[3*corner] )] &&
              isZero[DoubleProduct( EDGES_FOR_CORNER[3*corner+1] )] &&
              isZero[DoubleProduct( EDGES_FOR_CORNER[3*corner+2] )];

            if(doTest && testSegmentCornerIntersection(seg, corner))
              {
                double* ptA = new double[3];
                copyVector3(&COORDS_TET_CORNER[3 * corner], ptA);
                _polygonA.push_back(ptA);
                LOG(3,"Segment-corner : " << vToStr(ptA) << " added to A");
                if(corner != O)
                  {
                    double* ptB = new double[3];
                    _polygonB.push_back(ptB);
                    copyVector3(&COORDS_TET_CORNER[3 * corner], ptB);
                    LOG(3,"Segment-corner : " << vToStr(ptB) << " added to B");
                  }
              }
          }

            // segment - ray 
            for(TetraCorner corner = X ; corner < NO_TET_CORNER ; corner = TetraCorner(corner + 1))
              {
                if(isZero[DP_SEGMENT_RAY_INTERSECTION[7*(corner-1)]] && testSegmentRayIntersection(seg, corner))
                  {
                    double* ptB = new double[3];
                    copyVector3(&COORDS_TET_CORNER[3 * corner], ptB);
                    _polygonB.push_back(ptB);
                    LOG(3,"Segment-ray : " << vToStr(ptB) << " added to B");
                  }
              }
       
            // segment - halfstrip
            for(TetraEdge edge = XY ; edge <= ZX ; edge = TetraEdge(edge + 1))
              {

#if 0
                const int edgeIdx = int(edge) - 3; // offset since we only care for edges XY - ZX
                const bool doTest = 
                  !isZero[DP_FOR_HALFSTRIP_INTERSECTION[4*edgeIdx]] &&
                  !isZero[DP_FOR_HALFSTRIP_INTERSECTION[4*edgeIdx+1]];
       

                if(doTest && testSegmentHalfstripIntersection(seg, edge))
#endif
                  if(testSegmentHalfstripIntersection(seg, edge))
                    {
                      double* ptB = new double[3];
                      calcIntersectionPtSegmentHalfstrip(seg, edge, ptB);
                      _polygonB.push_back(ptB);
                      LOG(3,"Segment-halfstrip : " << vToStr(ptB) << " added to B");
                    }
              }
      }

        // inclusion tests
        for(TriCorner corner = P ; corner < NO_TRI_CORNER ; corner = TriCorner(corner + 1))
          {
            // { XYZ - inclusion only possible if in Tetrahedron?
            // tetrahedron
            if(testCornerInTetrahedron(corner))
              {
                double* ptA = new double[3];
                copyVector3(&_coords[5*corner], ptA);
                _polygonA.push_back(ptA);
                LOG(3,"Inclusion tetrahedron : " << vToStr(ptA) << " added to A");
              }

            // XYZ - plane
            if(testCornerOnXYZFacet(corner))
              {
                double* ptB = new double[3];
                copyVector3(&_coords[5*corner], ptB);
                _polygonB.push_back(ptB);
                LOG(3,"Inclusion XYZ-plane : " << vToStr(ptB) << " added to B");
              }

            // projection on XYZ - facet
            if(testCornerAboveXYZFacet(corner))
              {
                double* ptB = new double[3];
                copyVector3(&_coords[5*corner], ptB);
                ptB[2] = 1 - ptB[0] - ptB[1];
                assert(epsilonEqual(ptB[0]+ptB[1]+ptB[2] - 1, 0.0));
                _polygonB.push_back(ptB);
                LOG(3,"Projection XYZ-plane : " << vToStr(ptB) << " added to B");
              }

          }

      }

  /**
   * Calculates the intersection polygon A, performing the intersection tests
   * and storing the corresponding point in the vector _polygonA.
   *
   * @post _polygonA contains the intersection polygon A.
   *
   */
  void TransformedTriangle::calculateIntersectionPolygon()
  {
    assert(_polygonA.size() == 0);
    // avoid reallocations in push_back() by pre-allocating enough memory
    // we should never have more than 20 points
    _polygonA.reserve(20);
    // -- surface intersections
    // surface - edge
    for(TetraEdge edge = OX ; edge <= ZX ; edge = TetraEdge(edge + 1))
      {
        if(testSurfaceEdgeIntersection(edge))
          {
            double* ptA = new double[3];
            calcIntersectionPtSurfaceEdge(edge, ptA);
            _polygonA.push_back(ptA);
            LOG(3,"Surface-edge : " << vToStr(ptA) << " added to A ");
          }
      }

    // -- segment intersections
    for(TriSegment seg = PQ ; seg < NO_TRI_SEGMENT ; seg = TriSegment(seg + 1))
      {

        bool isZero[NO_DP];

        // check beforehand which double-products are zero
        for(DoubleProduct dp = C_YZ; dp < NO_DP; dp = DoubleProduct(dp + 1))
          {
            isZero[dp] = (calcStableC(seg, dp) == 0.0);
          }

        // segment - facet
        for(TetraFacet facet = OYZ ; facet < NO_TET_FACET ; facet = TetraFacet(facet + 1))
          {
            // is this test worth it?
            const bool doTest =
              !isZero[DP_FOR_SEG_FACET_INTERSECTION[3*facet]] &&
              !isZero[DP_FOR_SEG_FACET_INTERSECTION[3*facet + 1]] &&
              !isZero[DP_FOR_SEG_FACET_INTERSECTION[3*facet + 2]];

            if(doTest && testSegmentFacetIntersection(seg, facet))
              {
                double* ptA = new double[3];
                calcIntersectionPtSegmentFacet(seg, facet, ptA);
                _polygonA.push_back(ptA);
                LOG(3,"Segment-facet : " << vToStr(ptA) << " added to A");
              }
          }

        // segment - edge
        for(TetraEdge edge = OX ; edge <= ZX ; edge = TetraEdge(edge + 1))
          {
            const DoubleProduct edge_dp = DoubleProduct(edge);

            if(isZero[edge_dp] && testSegmentEdgeIntersection(seg, edge))
              {
                double* ptA = new double[3];
                calcIntersectionPtSegmentEdge(seg, edge, ptA);
                _polygonA.push_back(ptA);
                LOG(3,"Segment-edge : " << vToStr(ptA) << " added to A");
              }
          }

        // segment - corner
        for(TetraCorner corner = O ; corner < NO_TET_CORNER ; corner = TetraCorner(corner + 1))
          {
            const bool doTest =
              isZero[DoubleProduct( EDGES_FOR_CORNER[3*corner] )] &&
              isZero[DoubleProduct( EDGES_FOR_CORNER[3*corner+1] )] &&
              isZero[DoubleProduct( EDGES_FOR_CORNER[3*corner+2] )];

            if(doTest && testSegmentCornerIntersection(seg, corner))
              {
                double* ptA = new double[3];
                copyVector3(&COORDS_TET_CORNER[3 * corner], ptA);
                _polygonA.push_back(ptA);
                LOG(3,"Segment-corner : " << vToStr(ptA) << " added to A");
              }
          }

      }

        // inclusion tests
        for(TriCorner corner = P ; corner < NO_TRI_CORNER ; corner = TriCorner(corner + 1))
          {
            // { XYZ - inclusion only possible if in Tetrahedron?
            // tetrahedron
            if(testCornerInTetrahedron(corner))
              {
                double* ptA = new double[3];
                copyVector3(&_coords[5*corner], ptA);
                _polygonA.push_back(ptA);
                LOG(3,"Inclusion tetrahedron : " << vToStr(ptA) << " added to A");
              }

          }

      }


    /**
     * Returns the surface of polygon A.
     *
     * @return the surface of polygon A.
     */
    double TransformedTriangle::calculateSurfacePolygon()
    {
      const std::size_t nbPoints = _polygonA.size();
      double pdt[3];
      double sum[3] = {0., 0., 0.};

      for(std::size_t i = 0 ; i < nbPoints ; ++i)
        {
          const double *const ptCurr = _polygonA[i];  // pt "i"
          const double *const ptNext = _polygonA[(i + 1) % nbPoints]; // pt "i+1" (pt nbPoints == pt 0)

          cross(ptCurr, ptNext, pdt);
          add(pdt, sum);
        }

      const double surface = norm(sum) * 0.5;
      LOG(2,"Surface is " << surface);
      return surface;
    }

    /**
     * Calculates the barycenters of the given intersection polygon.
     *
     * @pre  the intersection polygons have been calculated with calculateIntersectionAndProjectionPolygons()
     * 
     * @param poly        one of the two intersection polygons
     * @param barycenter  array of three doubles where barycenter is stored
     *
     */
    void TransformedTriangle::calculatePolygonBarycenter(const IntersectionPolygon poly, double* barycenter)
    {
      LOG(3,"--- Calculating polygon barycenter");

      // get the polygon points
      std::vector<double*>& polygon = (poly == A) ? _polygonA : _polygonB;

      // calculate barycenter
      const std::size_t m = polygon.size();

      for(int j = 0 ; j < 3 ; ++j)
        {
          barycenter[j] = 0.0;
        }

      if(m != 0)
        {
          for(std::size_t i = 0 ; i < m ; ++i)
            {
              const double* pt = polygon[i];
              for(int j = 0 ; j < 3 ; ++j)
                {
                  barycenter[j] += pt[j] / double(m);
                }
            }
        }
      LOG(3,"Barycenter is " << vToStr(barycenter));
    }

    /**
     * Sorts the given intersection polygon in circular order around its barycenter.
     * @pre  the intersection polygons have been calculated with calculateIntersectionAndProjectionPolygons()
     * @post the vertices in _polygonA and _polygonB are sorted in circular order around their
     *       respective barycenters
     *
     * @param poly        one of the two intersection polygons
     * @param barycenter  array of three doubles with the coordinates of the barycenter
     * 
     */
    void TransformedTriangle::sortIntersectionPolygon(const IntersectionPolygon poly, const double* barycenter)
    {
      LOG(3,"--- Sorting polygon ...");

      using INTERP_KERNEL::ProjectedCentralCircularSortOrder;
      typedef ProjectedCentralCircularSortOrder SortOrder; // change is only necessary here and in constructor
      typedef SortOrder::CoordType CoordType;

      // get the polygon points
      std::vector<double*>& polygon = (poly == A) ? _polygonA : _polygonB;

      if(polygon.size() == 0)
        return;

      // determine type of sorting
      CoordType type = SortOrder::XY;
      if(poly == A && !isTriangleInclinedToFacet( OXY )) // B is on h = 0 plane -> ok
        {
          // NB : the following test is never true if we have eliminated the
          // triangles parallel to x == 0 and y == 0 in calculateIntersectionVolume().
          // We keep the test here anyway, to avoid interdependency.

          // is triangle inclined to x == 0 ?
          if(isTriangleInclinedToFacet(OZX))
            {
              type = SortOrder::XZ;
            }
          else //if(isTriangleParallelToFacet(OYZ))
            {
              type = SortOrder::YZ;
            }
        }

      // create order object
      SortOrder order(barycenter, type);

      // sort vector with this object
      // NB : do not change place of first object, with respect to which the order
      // is defined
      sort((polygon.begin()), polygon.end(), order);
    
      LOG(3,"Sorted polygon is ");
      for(size_t i = 0 ; i < polygon.size() ; ++i)
        {
          LOG(3,vToStr(polygon[i]));
        }

    }

    /**
     * Calculates the volume between the given polygon and the z = 0 plane.
     *
     * @pre  the intersection polygones have been calculated with calculateIntersectionAndProjectionPolygons(),
     *       and they have been sorted in circular order with sortIntersectionPolygons(void)
     * 
     * @param poly        one of the two intersection polygons
     * @param barycenter  array of three doubles with the coordinates of the barycenter
     * @return           the volume between the polygon and the z = 0 plane
     *
     */
    double TransformedTriangle::calculateVolumeUnderPolygon(IntersectionPolygon poly, const double* barycenter)
    {
      LOG(2,"--- Calculating volume under polygon");

      // get the polygon points
      std::vector<double*>& polygon = (poly == A) ? _polygonA : _polygonB;

      double vol = 0.0;
      const std::size_t m = polygon.size();

      for(std::size_t i = 0 ; i < m ; ++i)
        {
          const double* ptCurr = polygon[i];  // pt "i"
          const double* ptNext = polygon[(i + 1) % m]; // pt "i+1" (pt m == pt 0)
       
          const double factor1 = ptCurr[2] + ptNext[2] + barycenter[2];
          const double factor2 = 
            ptCurr[0]*(ptNext[1] - barycenter[1]) 
            + ptNext[0]*(barycenter[1] - ptCurr[1])
            + barycenter[0]*(ptCurr[1] - ptNext[1]);
          vol += (factor1 * factor2) / 6.0;
        }

      LOG(2,"Abs. Volume is " << vol); 
      return vol;
    }


    ////////////////////////////////////////////////////////////////////////////////////
    // Detection of (very) degenerate cases                                /////////////
    ////////////////////////////////////////////////////////////////////////////////////

    /**
     * Checks if the triangle lies in the plane of a given facet
     *
     * @param facet     one of the facets of the tetrahedron
     * @return         true if PQR lies in the plane of the facet, false if not
     */
    bool TransformedTriangle::isTriangleInPlaneOfFacet(const TetraFacet facet) const
    {

      // coordinate to check
      const int coord = static_cast<int>(facet);

      for(TriCorner c = P ; c < NO_TRI_CORNER ; c = TriCorner(c + 1))
        {
          if(_coords[5*c + coord] != 0.0)
            {
              return false;
            }
        }
    
      return true;
    }

    /**
     * Checks if the triangle is parallel to the given facet
     *
     * @param facet  one of the facets of the unit tetrahedron
     * @return       true if triangle is parallel to facet, false if not
     */
    bool TransformedTriangle::isTriangleParallelToFacet(const TetraFacet facet) const
    {
      // coordinate to check
      const int coord = static_cast<int>(facet);
      return (_coords[5*P + coord] == _coords[5*Q + coord]) && (_coords[5*P + coord] == _coords[5*R + coord]);
    }

    /**
     * Checks if the triangle is not perpedicular to the given facet
     *
     * @param facet  one of the facets of the unit tetrahedron
     * @return       zero if the triangle is perpendicular to the facet,
     *               else 1 or -1 depending on the sign of cross product of facet edges
     */
    int TransformedTriangle::isTriangleInclinedToFacet(const TetraFacet facet) const
    {
      // coordinate to check
      const int coord = static_cast<int>(facet);
      const int ind1 = ( coord+1 ) % 3, ind2 = ( coord+2 ) % 3;
      const double uv_xy[4] = 
        {
          // u_x, u_y
          _coords[5*Q+ind1] - _coords[5*P+ind1], _coords[5*Q+ind2] - _coords[5*P+ind2],
          // v_x, v_y
          _coords[5*R+ind1] - _coords[5*P+ind1], _coords[5*R+ind2] - _coords[5*P+ind2]
        };

      double sign = uv_xy[0] * uv_xy[3] - uv_xy[1] * uv_xy[2];
      if(epsilonEqual(sign, 0.))
        {
          sign = 0.;
        }
      return (sign < 0.) ? -1 : (sign > 0.) ? 1 : 0;
    }

    /**
     * Determines whether the triangle is below the z-plane.
     * 
     * @return true if the z-coordinate of the three corners of the triangle are all less than 0, false otherwise.
     */
    bool TransformedTriangle::isTriangleBelowTetraeder() const
    {
      for(TriCorner c = P ; c < NO_TRI_CORNER ; c = TriCorner(c + 1))
        {
          // check z-coords for all points
          if(_coords[5*c + 2] >= 0.0)
            {
              return false;
            }
        }
      return true;
    }

    /**
     * Prints the coordinates of the triangle to std::cout
     *
     */
    void TransformedTriangle::dumpCoords() const
    {
      std::cout << "Coords : ";
      for(int i = 0 ; i < 3; ++i)
        {
          std::cout << vToStr(&_coords[5*i]) << ",";
        }
      std::cout << std::endl;
    }
    
  } // NAMESPACE 
