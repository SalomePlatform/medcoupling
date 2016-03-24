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

#ifndef __TRANSFORMED_TRIANGLE_HXX__
#define __TRANSFORMED_TRIANGLE_HXX__

#include "INTERPKERNELDefines.hxx"

#include <vector>

// Levels : 
// 1 - overview of algorithm + volume result
// 2 - algorithm detail
// 3 - intersection polygon results detail
// 4 - intersection polygon search detail
// higher -> misc. gory details of calculation

#include "Log.hxx"

#ifdef WIN32
#pragma warning(disable:4251)
#endif

namespace INTERP_TEST
{
  class TransformedTriangleTest;
  class TransformedTriangleIntersectTest;
}


namespace INTERP_KERNEL
{
  class TetraAffineTransform;

  /** \class TransformedTriangle
   * \brief Class representing one of the faces of the triangulated source polyhedron after having been transformed
   * with the affine transform that takes the target tetrahedron to the unit tetrahedron. It contains the
   * logic for calculating the volume of intersection between the triangle and the unit tetrahedron.
   * 
   * \see TransformedTriangle.hxx 
   *
   * Reference : J. Grandy, "Conservative Remapping and Region Overlays by Intersecting Arbitrary Polyhedra", 
   *             Journal of Computational Physics (1999)
   *
   */

  /** \file TransformedTriangle.hxx
   * 
   * OVERVIEW of how the class works : (details can be found in the documentation of each method)
   * 
   * Constructor : 
   * The constructor takes as arguments three pointers to double[3] vectors holding the transformed
   * coordinates of the corners of the triangle. It copies their coordinates and then proceeds to pre-calculating certain
   * entities used in the intersection calculation : the double products, triple products and the values of the function E
   * (Grandy, [53]).
   *
   * calculateIntersectionVolume() : 
   * This is the only method in the public interface. It calculates the volume under the intersection polygons
   * between the triangle and the unit tetrahedron, as described in Grandy, pp. 435-447. It does this by first calculating the
   * intersection polygons A and B, with the method calculateIntersectionPolygons(). It then calculates the barycenter of each
   * polygon in calculatePolygonBarycenter(), and sorts their points in a circular order around the barycenter in 
   * sortIntersecionPolygon(). The sorting is done with STL sort, using the order defined in the class 
   * ProjectedCentralCircularSortOrder. The volume under each polygon is then calculated with calculateVolumeUnderPolygon(), which
   * implements formula [34] in Grandy.
   *
   * calculateIntersectionPolygons() :
   * This method goes through all the possible ways in which the triangle can intersect the tetrahedron and tests for these 
   * types of intersections in accordance with the formulas described in Grandy. These tests are implemented in the test* - methods.
   *    The formulas in the article are stated for one case each only, while the calculation must take into account all cases. 
   * To this end, a number of tables, implemented as static const arrays of different types, are used. The tables 
   * mainly contain values of the different enumeration types described at the beginning of the class interface. For example, 
   * the formula Grandy gives for the segment-halfstrip intersection tests ([30]) is for use with the halfstrip above the zx edge. 
   * For the other two halfstrips (above the xy and yz edges), other double products are used, which 
   * are stored in the table DP_FOR_HALFSTRIP_INTERSECTION. This allows us to treat
   * all the edges equally, avoiding switch() - statements. It is the careful choice of order of the enumeration types that makes this
   * possible. Notably, there is a correspondance between the TetraEdge type and the DoubleProduct type (see Grandy, table III) that
   * is used throughout the code, permitting statements such as DoubleProduct(some_edge) to work.
   *    When an intersection point has been detected it is calculated with a corresponding calc* - method in the cases where it
   * is not known directly. It is then added to the polygon A and/or B as necessary.
   *
   * OPTIMIZE : 
   *    If OPTIMIZE is defined, a large number of methods will be prefixed with inline and some optimizations concerning the tests 
   * with zero double products will be used.
   */
  class INTERPKERNEL_EXPORT TransformedTriangle
  {
 

  public:

    friend class INTERP_TEST::TransformedTriangleTest;
    friend class INTERP_TEST::TransformedTriangleIntersectTest;
    /*
     * Enumerations representing the different geometric elements of the unit tetrahedron
     * and the triangle. The end element, NO_* gives the number of elements in the enumeration
     * and can be used as end element in loops.
     */

    /// Corners of tetrahedron
    enum TetraCorner { O = 0, X, Y, Z, NO_TET_CORNER };

    /// Edges of tetrahedron
    enum TetraEdge { OX = 0, OY, OZ, XY, YZ, ZX, H01, H10, NO_TET_EDGE };

    /// Facets (faces) of tetrahedron
    enum TetraFacet { OYZ = 0, OZX, OXY, XYZ, NO_TET_FACET };

    /// Corners of triangle
    enum TriCorner { P = 0, Q, R, NO_TRI_CORNER };
    
    /// Segments (edges) of triangle
    enum TriSegment { PQ = 0, QR, RP, NO_TRI_SEGMENT };
    
    /// Intersection polygons
    enum IntersectionPolygon{ A = 0, B, NO_INTERSECTION_POLYGONS };

    /// Double products
    /// NB : order corresponds to TetraEdges (Grandy, table III)
    enum DoubleProduct { C_YZ = 0, C_ZX, C_XY, C_ZH, C_XH, C_YH, C_01, C_10, NO_DP };

    TransformedTriangle(double* p, double* q, double* r); 
    ~TransformedTriangle();

    double calculateIntersectionVolume(); 
    double calculateIntersectionSurface(TetraAffineTransform* tat);

    void dumpCoords() const;

    // Queries of member values used by UnitTetraIntersectionBary

    const double* getCorner(TriCorner corner) const { return _coords + 5*corner; }

    const std::vector<double*>& getPolygonA() const { return _polygonA; }

    double getVolume() const { return _volume; }

  protected:

    TransformedTriangle() { }

    // ----------------------------------------------------------------------------------
    //  High-level methods called directly by calculateIntersectionVolume()     
    // ----------------------------------------------------------------------------------
    void calculateIntersectionAndProjectionPolygons();

    void calculatePolygonBarycenter(const IntersectionPolygon poly, double* barycenter); 

    void sortIntersectionPolygon(const IntersectionPolygon poly, const double* barycenter); 

    double calculateVolumeUnderPolygon(IntersectionPolygon poly, const double* barycenter); 

    // ----------------------------------------------------------------------------------
    //  High-level methods called directly by calculateIntersectionSurface()
    // ----------------------------------------------------------------------------------
    void calculateIntersectionPolygon();

    double calculateSurfacePolygon();

    // ----------------------------------------------------------------------------------
    //  Detection of degenerate triangles  
    // ----------------------------------------------------------------------------------

    bool isTriangleInPlaneOfFacet(const TetraFacet facet) const;
    
    bool isTriangleParallelToFacet(const TetraFacet facet) const;

    int isTriangleInclinedToFacet(const TetraFacet facet) const;

    bool isTriangleBelowTetraeder() const;

    // ----------------------------------------------------------------------------------
    //  Intersection test methods and intersection point calculations           
    // ----------------------------------------------------------------------------------
 
    inline bool testSurfaceEdgeIntersection(const TetraEdge edge) const; 

    void calcIntersectionPtSurfaceEdge(const TetraEdge edge, double* pt) const;  

    inline bool testSegmentFacetIntersection(const TriSegment seg, const TetraFacet facet) const; 

    void calcIntersectionPtSegmentFacet(const TriSegment seg, const TetraFacet facet, double* pt) const;  

    bool testSegmentEdgeIntersection(const TriSegment seg, const TetraEdge edge) const; 
 
    void calcIntersectionPtSegmentEdge(const TriSegment seg, const TetraEdge edge, double* pt) const ; 

    bool testSegmentCornerIntersection(const TriSegment seg, const TetraCorner corner) const ;

    inline bool testSurfaceRayIntersection(const TetraCorner corner) const;

    bool testSegmentHalfstripIntersection(const TriSegment seg, const TetraEdge edg);

    void calcIntersectionPtSegmentHalfstrip(const TriSegment seg, const TetraEdge edge, double* pt) const;
    
    bool testSegmentRayIntersection(const TriSegment seg, const TetraCorner corner) const;

    inline bool testCornerInTetrahedron(const TriCorner corner) const;

    inline bool testCornerOnXYZFacet(const TriCorner corner) const;

    inline bool testCornerAboveXYZFacet(const TriCorner corner) const;

    // ----------------------------------------------------------------------------------
    //  Utility methods used in intersection tests                       
    // ----------------------------------------------------------------------------------
    
    bool testTriangleSurroundsEdge(const TetraEdge edge) const;

    inline bool testEdgeIntersectsTriangle(const TetraEdge edge) const;

    inline bool testFacetSurroundsSegment(const TriSegment seg, const TetraFacet facet) const;

    inline bool testSegmentIntersectsFacet(const TriSegment seg, const TetraFacet facet) const;

    bool testSegmentIntersectsHPlane(const TriSegment seg) const;

    bool testSurfaceAboveCorner(const TetraCorner corner) const;
    
    bool testTriangleSurroundsRay(const TetraCorner corner) const;

    // ----------------------------------------------------------------------------------
    //  Double and triple product calculations                           
    // ----------------------------------------------------------------------------------
    
    void resetNearZeroCoordinates();

    bool areDoubleProductsConsistent(const TriSegment seg) const;

    void preCalculateDoubleProducts(void);

    inline void resetDoubleProducts(const TriSegment seg, const TetraCorner corner);

    double calculateDistanceCornerSegment(const TetraCorner corner, const TriSegment seg) const;
    
    void preCalculateTripleProducts(void);

    double calculateAngleEdgeTriangle(const TetraEdge edge) const;

    inline double calcStableC(const TriSegment seg, const DoubleProduct dp) const;

    inline double calcStableT(const TetraCorner corner) const;

    inline double calcUnstableC(const TriSegment seg, const DoubleProduct dp) const;

    double calcTByDevelopingRow(const TetraCorner corner, const int row = 1, const bool project = false) const;

    // ----------------------------------------------------------------------------------
    //  Member variables                                                 
    // ----------------------------------------------------------------------------------
  protected:

    /// Array holding the coordinates of the triangle's three corners
    /// order : 
    /// [ p_x, p_y, p_z, p_h, p_H, q_x, q_y, q_z, q_h, q_H, r_x, r_y, r_z, r_h, r_H ]
    double _coords[15];
    
    /// Flag showing whether the double products have been calculated yet
    bool _is_double_products_calculated;

    /// Flag showing whether the triple products have been calculated yet
    bool _is_triple_products_calculated; 

    /// Array containing the 24 double products.
    /// order : c^PQ_YZ, ... ,cPQ_10, ... c^QR_YZ, ... c^RP_YZ
    /// following order in enumeration DoubleProduct
    double _doubleProducts[24];

    /// Array containing the 4 triple products.
    /// order : t_O, t_X, t_Y, t_Z
    double _tripleProducts[4];

    /// Vector holding the points of the intersection polygon A.
    /// these points are allocated in calculateIntersectionPolygons() and liberated in the destructor
    std::vector<double*> _polygonA;
    
    /// Vector holding the points of the intersection polygon B.
    /// These points are allocated in calculateIntersectionPolygons() and liberated in the destructor
    std::vector<double*> _polygonB;
    
    /// Array holding the coordinates of the barycenter of the polygon A
    /// This point is calculated in calculatePolygonBarycenter
    double _barycenterA[3];

    /// Array holding the coordinates of the barycenter of the polygon B
    /// This point is calculated in calculatePolygonBarycenter
    //double _barycenterB[3];

    /// Array of flags indicating which of the four triple products have been correctly calculated.
    /// Used for asserts in debug mode
    bool _validTP[4];

    /// calculated volume for use of UnitTetraIntersectionBary
    double _volume;
    
    /**
     * Calls TransformedTriangle::testTriangleSurroundsEdge for edges OX to ZX and stores the result in
     * member variable array_triangleSurroundsEdgeCache. 
     *
     */
    void preCalculateTriangleSurroundsEdge();

    /// Array holding results of the test testTriangleSurroundsEdge() for all the edges. 
    /// These are calculated in preCalculateTriangleSurroundsEdge().
    bool _triangleSurroundsEdgeCache[NO_TET_EDGE];

    // ----------------------------------------------------------------------------------
    //  Constants                                                    
    // ----------------------------------------------------------------------------------

    // offsets : 0 -> x, 1 -> y, 2 -> z, 3 -> h, 4 -> H
    // corresponds to order of double products in DoubleProduct
    // so that offset[C_*] gives the right coordinate
    static const int DP_OFFSET_1[8];
    static const int DP_OFFSET_2[8];

    // the coordinates used in the expansion of triple products by a given row
    // in constellation (corner, row-1)
    // (0,1,2,3) <=> (x,y,z,h)
    static const int COORDINATE_FOR_DETERMINANT_EXPANSION[12];

    // contains the edge of the double product used when 
    // expanding the triple product determinant associated with each corner
    // by a given row
    static const DoubleProduct DP_FOR_DETERMINANT_EXPANSION[12];
    
    // values used to decide how imprecise the double products 
    // should be to set them to 0.0
    static const long double MACH_EPS;    // machine epsilon
    static const long double MULT_PREC_F; // precision of multiplications (Grandy : f)
    static const long double THRESHOLD_F; // threshold for zeroing (Grandy : F/f)

    static const double TRIPLE_PRODUCT_ANGLE_THRESHOLD;

    // correspondance facet - double product
    // Grandy, table IV
    static const DoubleProduct DP_FOR_SEG_FACET_INTERSECTION[12];

    // signs associated with entries in DP_FOR_SEGMENT_FACET_INTERSECTION
    static const double SIGN_FOR_SEG_FACET_INTERSECTION[12];
    
    // coordinates of corners of tetrahedron
    static const double COORDS_TET_CORNER[12];
    
    // indices to use in tables DP_FOR_SEG_FACET_INTERSECTION and SIGN_FOR_SEG_FACET_INTERSECTION
    // for the calculation of the coordinates (x,y,z) of the intersection points
    // for Segment-Facet and Segment-Edge intersections
    static const int DP_INDEX[12];

    // correspondance edge - corners
    static const TetraCorner CORNERS_FOR_EDGE[12];

    // correspondance edge - facets
    // facets shared by each edge
    static const TetraFacet FACET_FOR_EDGE[12];

    // correspondance edge - corners
    static const TetraEdge EDGES_FOR_CORNER[12];
   
    // double products used in segment-halfstrip test
    static const DoubleProduct DP_FOR_HALFSTRIP_INTERSECTION[12];

    // double products used in segment - ray test
    static const DoubleProduct DP_SEGMENT_RAY_INTERSECTION[21];

  };

  // include definitions of inline methods

#include "TransformedTriangleInline.hxx"
}


#endif
