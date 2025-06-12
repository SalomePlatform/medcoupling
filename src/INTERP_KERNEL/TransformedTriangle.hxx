// Copyright (C) 2007-2025  CEA, EDF
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
#include "VectorUtils.hxx"
#include "assert.h"

#include <vector>

// Levels :
// 1 - overview of algorithm + volume result
// 2 - algorithm detail
// 3 - intersection polygon results detail
// 4 - intersection polygon search detail
// higher -> misc. gory details of calculation

#include "Log.hxx"

#ifdef WIN32
#pragma warning(disable : 4251)
#endif

namespace INTERP_TEST
{
class TransformedTriangleTest;
class TransformedTriangleIntersectTest;
}  // namespace INTERP_TEST

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

/**
 * OVERVIEW of how the class works : (details can be found in the documentation of each method)
 *
 * Constructor :
 * The constructor takes as arguments three pointers to double[3] vectors holding the transformed
 * coordinates of the corners of the triangle. It copies their coordinates and then proceeds to pre-calculating certain
 * entities used in the intersection calculation : the double products, triple products and the values of the function E
 * (Grandy, [53]).
 * It is also at this point in constructor that:
 *  - the special case of PQR included in the XYZ plane is treated
 *  - the inconsistencies between double products/triple products computation is handled
 *
 * calculateIntersectionVolume() :
 * This is the only method in the public interface. It calculates the volume under the intersection polygons
 * between the triangle and the unit tetrahedron, as described in Grandy, pp. 435-447. It does this by first calculating
 * the intersection polygons A and B, with the method calculateIntersectionPolygons(). It then calculates the barycenter
 * of each polygon in calculatePolygonBarycenter(), and sorts their points in a circular order around the barycenter in
 * sortIntersecionPolygon(). The sorting is done with STL sort, using the order defined in the class
 * ProjectedCentralCircularSortOrder. The volume under each polygon is then calculated with
 * calculateVolumeUnderPolygon(), which implements formula [34] in Grandy.
 *
 * calculateIntersectionPolygons() :
 * This method goes through all the possible ways in which the triangle can intersect the tetrahedron and tests for
 * these types of intersections in accordance with the formulas described in Grandy. These tests are implemented in the
 * test* - methods. The formulas in the article are stated for one case each only, while the calculation must take into
 * account all cases. To this end, a number of tables, implemented as static const arrays of different types, are used.
 * The tables mainly contain values of the different enumeration types described at the beginning of the class
 * interface. For example, the formula Grandy gives for the segment-halfstrip intersection tests ([30]) is for use with
 * the halfstrip above the zx edge. For the other two halfstrips (above the xy and yz edges), other double products are
 * used, which are stored in the table DP_FOR_HALFSTRIP_INTERSECTION. This allows us to treat all the edges equally,
 * avoiding switch() - statements. It is the careful choice of order of the enumeration types that makes this possible.
 * Notably, there is a correspondence between the TetraEdge type and the DoubleProduct type (see Grandy, table III) that
 * is used throughout the code, permitting statements such as DoubleProduct(some_edge) to work.
 *    When an intersection point has been detected it is calculated with a corresponding calc* - method in the cases
 * where it is not known directly. It is then added to the polygon A and/or B as necessary.
 *
 */
class INTERPKERNEL_EXPORT TransformedTriangle
{
   public:
    friend class INTERP_TEST::TransformedTriangleIntersectTest;
    friend class INTERP_TEST::TransformedTriangleTest;
    /*
     * Enumerations representing the different geometric elements of the unit tetrahedron
     * and the triangle. The end element, NO_* gives the number of elements in the enumeration
     * and can be used as end element in loops.
     */

    /// Corners of tetrahedron
    enum TetraCorner
    {
        O = 0,
        X,
        Y,
        Z,
        NO_TET_CORNER
    };

    /// Edges of tetrahedron
    enum TetraEdge
    {
        OX = 0,
        OY,
        OZ,
        XY,
        YZ,
        ZX,
        H01,
        H10,
        NO_TET_EDGE
    };

    /// Facets (faces) of tetrahedron
    enum TetraFacet
    {
        OYZ = 0,
        OZX,
        OXY,
        XYZ,
        NO_TET_FACET
    };

    /// Corners of triangle
    enum TriCorner
    {
        P = 0,
        Q,
        R,
        NO_TRI_CORNER
    };

    /// Segments (edges) of triangle
    enum TriSegment
    {
        PQ = 0,
        QR,
        RP,
        NO_TRI_SEGMENT
    };

    /// Intersection polygons
    enum IntersectionPolygon
    {
        A = 0,
        B,
        NO_INTERSECTION_POLYGONS
    };

    /// Double products
    /// NB : order corresponds to TetraEdges (Grandy, table III)
    enum DoubleProduct
    {
        C_YZ = 0,
        C_ZX,
        C_XY,
        C_ZH,
        C_XH,
        C_YH,
        C_01,
        C_10,
        NO_DP
    };

    TransformedTriangle(double *p, double *q, double *r);
    ~TransformedTriangle();

    double calculateIntersectionVolume();
    double calculateIntersectionSurface(TetraAffineTransform *tat);
    void dumpCoords() const;

    // Queries of member values used by UnitTetraIntersectionBary
    const double *getCorner(TriCorner corner) const { return _coords + 5 * corner; }
    const std::vector<double *> &getPolygonA() const { return _polygonA; }
    double getVolume() const { return _volume; }

   protected:
    TransformedTriangle() {}

    // ----------------------------------------------------------------------------------
    //  High-level methods called directly by calculateIntersectionVolume()
    // ----------------------------------------------------------------------------------
    void calculateIntersectionAndProjectionPolygons();
    void calculatePolygonBarycenter(const IntersectionPolygon poly, double *barycenter);
    void sortIntersectionPolygon(const IntersectionPolygon poly, const double *barycenter);
    double calculateVolumeUnderPolygon(IntersectionPolygon poly, const double *barycenter);

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
    void calcIntersectionPtSurfaceEdge(const TetraEdge edge, double *pt) const;
    inline bool testSegmentFacetIntersection(const TriSegment seg, const TetraFacet facet) const;
    void calcIntersectionPtSegmentFacet(const TriSegment seg, const TetraFacet facet, double *pt) const;
    bool testSegmentEdgeIntersection(const TriSegment seg, const TetraEdge edge) const;
    void calcIntersectionPtSegmentEdge(const TriSegment seg, const TetraEdge edge, double *pt) const;
    bool testSegmentCornerIntersection(const TriSegment seg, const TetraCorner corner) const;
    inline bool testSurfaceRayIntersection(const TetraCorner corner) const;
    bool testSegmentHalfstripIntersection(const TriSegment seg, const TetraEdge edg);
    void calcIntersectionPtSegmentHalfstrip(const TriSegment seg, const TetraEdge edge, double *pt) const;
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
    void handleDegenerateCases();
    bool areDoubleProductsConsistent(const TriSegment seg) const;
    void preCalculateDoubleProducts();
    inline void resetDoubleProducts(const TriSegment seg, const TetraCorner corner);
    double calculateDistanceCornerSegment(const TetraCorner corner, const TriSegment seg) const;
    void preCalculateTripleProducts();
    double calculateAngleEdgeTriangle(const TetraEdge edge) const;
    inline double calcStableC(const TriSegment seg, const DoubleProduct dp) const;
    inline double calcStableT(const TetraCorner corner) const;
    inline double calcUnstableC(const TriSegment seg, const DoubleProduct dp, double &delta) const;
    double calcTByDevelopingRow(const TetraCorner corner, const int row, const bool project) const;

    // ----------------------------------------------------------------------------------
    // Debug
    // ----------------------------------------------------------------------------------
    inline const std::string &strTC(TetraCorner tc) const;
    inline const std::string &strTE(TetraEdge te) const;
    inline const std::string &strTF(TetraFacet tf) const;
    inline const std::string &strTriC(TriCorner tc) const;
    inline const std::string &strTriS(TriSegment tc) const;

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

    double _deltas[24];

    /// Array containing the 4 triple products.
    /// order : t_O, t_X, t_Y, t_Z
    /// For example t_O represent the signed volume of the tetrahedron OPQR, and is positive if PQR is oriented
    /// clockwise
    //  when seen from the vertex O.
    double _tripleProducts[4];

    /// Vector holding the points of the intersection polygon A.
    /// these points are allocated in calculateIntersectionPolygons() and liberated in the destructor
    std::vector<double *> _polygonA;

    /// Vector holding the points of the intersection polygon B.
    /// These points are allocated in calculateIntersectionPolygons() and liberated in the destructor
    std::vector<double *> _polygonB;

    /// Array holding the coordinates of the barycenter of the polygon A
    /// This point is calculated in calculatePolygonBarycenter
    double _barycenterA[3];

    /// Array holding the coordinates of the barycenter of the polygon B
    /// This point is calculated in calculatePolygonBarycenter
    // double _barycenterB[3];

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

    // values used to decide how/when imprecise the double products
    // should be to set them to 0.0
    static const double MACH_EPS;     // machine epsilon
    static const double MULT_PREC_F;  // precision of multiplications (Grandy : f)
    static const double THRESHOLD_F;  // threshold for zeroing (Grandy : F/f)

    static const double TRIPLE_PRODUCT_ANGLE_THRESHOLD;

    // correspondence facet - double product
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

    // correspondence edge - corners
    static const TetraCorner CORNERS_FOR_EDGE[12];

    // correspondence edge - facets
    // facets shared by each edge
    static const TetraFacet FACET_FOR_EDGE[12];

    // correspondence edge - corners
    static const TetraEdge EDGES_FOR_CORNER[12];

    // double products used in segment-halfstrip test
    static const DoubleProduct DP_FOR_HALFSTRIP_INTERSECTION[12];

    // double products used in segment - ray test
    static const DoubleProduct DP_SEGMENT_RAY_INTERSECTION[21];
};

inline void
TransformedTriangle::preCalculateTriangleSurroundsEdge()
{
    for (TetraEdge edge = OX; edge <= ZX; edge = TetraEdge(edge + 1))
    {
        _triangleSurroundsEdgeCache[edge] = testTriangleSurroundsEdge(edge);
    }
}

// ----------------------------------------------------------------------------------
//   TransformedTriangle_math.cxx
// ----------------------------------------------------------------------------------

inline void
TransformedTriangle::resetDoubleProducts(const TriSegment seg, const TetraCorner corner)
{
    // set the three corresponding double products to 0.0
    static const DoubleProduct DOUBLE_PRODUCTS[12] = {
        C_YZ,
        C_ZX,
        C_XY,  // O
        C_YZ,
        C_ZH,
        C_YH,  // X
        C_ZX,
        C_ZH,
        C_XH,  // Y
        C_XY,
        C_YH,
        C_XH  // Z
    };

    for (int i = 0; i < 3; ++i)
    {
        const DoubleProduct dp = DOUBLE_PRODUCTS[3 * corner + i];

        LOG(6, std::endl << "resetting inconsistent dp :" << dp << " for corner " << corner);
        _doubleProducts[8 * seg + dp] = 0.0;
    };
}

inline double
TransformedTriangle::calcStableC(const TriSegment seg, const DoubleProduct dp) const
{
    return _doubleProducts[8 * seg + dp];
}

inline double
TransformedTriangle::calcStableT(const TetraCorner corner) const
{
    assert(_validTP[corner]);
    return _tripleProducts[corner];
}

inline double
TransformedTriangle::calcUnstableC(const TriSegment seg, const DoubleProduct dp, double &delta) const
{
    // find the points of the triangle
    // 0 -> P, 1 -> Q, 2 -> R
    const int pt1 = seg;
    const int pt2 = (seg + 1) % 3;

    // find offsets
    const int off1 = DP_OFFSET_1[dp];
    const int off2 = DP_OFFSET_2[dp];

    const double prd1 = _coords[5 * pt1 + off1] * _coords[5 * pt2 + off2],
                 prd2 = _coords[5 * pt1 + off2] * _coords[5 * pt2 + off1];
    delta = std::fabs(prd1) + std::fabs(prd2);
    return prd1 - prd2;
}

// ----------------------------------------------------------------------------------
//  TransformedTriangle_intersect.cxx
// ----------------------------------------------------------------------------------
inline bool
TransformedTriangle::testSurfaceEdgeIntersection(const TetraEdge edge) const
{
    return _triangleSurroundsEdgeCache[edge] && testEdgeIntersectsTriangle(edge);
}

inline bool
TransformedTriangle::testSegmentFacetIntersection(const TriSegment seg, const TetraFacet facet) const
{
    return testFacetSurroundsSegment(seg, facet) && testSegmentIntersectsFacet(seg, facet);
}

inline bool
TransformedTriangle::testSurfaceRayIntersection(const TetraCorner corner) const
{
    return testTriangleSurroundsRay(corner) && testSurfaceAboveCorner(corner);
}

inline bool
TransformedTriangle::testCornerInTetrahedron(const TriCorner corner) const
{
    const double pt[4] = {
        _coords[5 * corner],      // x
        _coords[5 * corner + 1],  // y
        _coords[5 * corner + 2],  // z
        _coords[5 * corner + 3]   // z
    };

    for (int i = 0; i < 4; ++i)
    {
        if (pt[i] < 0.0 || pt[i] > 1.0)
        {
            return false;
        }
    }
    return true;
}

inline bool
TransformedTriangle::testCornerOnXYZFacet(const TriCorner corner) const
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
    const double *pt = &_coords[5 * corner];

    if (pt[3] != 0.0)
    {
        return false;
    }

    for (int i = 0; i < 3; ++i)
    {
        if (pt[i] < 0.0 || pt[i] > 1.0)
        {
            return false;
        }
    }
    return true;
}

inline bool
TransformedTriangle::testCornerAboveXYZFacet(const TriCorner corner) const
{
    const double x = _coords[5 * corner];
    const double y = _coords[5 * corner + 1];
    const double h = _coords[5 * corner + 3];
    const double H = _coords[5 * corner + 4];

    return h < 0.0 && H >= 0.0 && x >= 0.0 && y >= 0.0;
}

inline bool
TransformedTriangle::testEdgeIntersectsTriangle(const TetraEdge edge) const
{
    // correspondence edge - triple products for edges OX, ..., ZX (Grandy, table III)
    static const TetraCorner TRIPLE_PRODUCTS[12] = {
        X,
        O,  // OX
        Y,
        O,  // OY
        Z,
        O,  // OZ
        X,
        Y,  // XY
        Y,
        Z,  // YZ
        Z,
        X,  // ZX
    };

    // Grandy, [16]
    const double t1 = calcStableT(TRIPLE_PRODUCTS[2 * edge]);
    const double t2 = calcStableT(TRIPLE_PRODUCTS[2 * edge + 1]);

    // [ABN] Okayyy: if either t1 or t2 exactly equal zero, then it can mean two things:
    //   - either PQR is very close to the corner -> this is OK, further computation of intersection point between
    // surface and edge will produce a correct result
    //   - or, if the other triple prod is also very small, then this is a degenerate case: the edge is almost in PQR ->
    //   this is bad
    // and leads to weird intersection point computation -> we avoid this.
    // PS : here was written "// tuleap26461" -> whoo this helps :-)
    if (t1 == 0.0 || t2 == 0.0)
        if (std::fabs(t1 + t2) < THRESHOLD_F * MULT_PREC_F)
            return false;

    return (t1 * t2 <= 0.0) && !epsilonEqual(t1, t2, MULT_PREC_F);
}

inline bool
TransformedTriangle::testFacetSurroundsSegment(const TriSegment seg, const TetraFacet facet) const
{
#if 0
    const double signs[3] =
    {
      SIGN_FOR_SEG_FACET_INTERSECTION[3*facet],
      SIGN_FOR_SEG_FACET_INTERSECTION[3*facet + 1],
      SIGN_FOR_SEG_FACET_INTERSECTION[3*facet + 2]
    };
#endif

    const double *signs = &SIGN_FOR_SEG_FACET_INTERSECTION[3 * facet];
    const double c1 = signs[0] * calcStableC(seg, DP_FOR_SEG_FACET_INTERSECTION[3 * facet]);
    const double c2 = signs[1] * calcStableC(seg, DP_FOR_SEG_FACET_INTERSECTION[3 * facet + 1]);
    const double c3 = signs[2] * calcStableC(seg, DP_FOR_SEG_FACET_INTERSECTION[3 * facet + 2]);

    return (c1 * c3 > 0.0) && (c2 * c3 > 0.0);
}

inline bool
TransformedTriangle::testSegmentIntersectsFacet(const TriSegment seg, const TetraFacet facet) const
{
    // use correspondence facet a = 0 <=> offset for coordinate a in _coords
    // and also correspondence segment AB => corner A
    const double coord1 = _coords[5 * seg + facet];
    const double coord2 = _coords[5 * ((seg + 1) % 3) + facet];

    //? should we use epsilon-equality here in second test?
    LOG(5, "coord1 : " << coord1 << " coord2 : " << coord2);

    return (coord1 * coord2 <= 0.0) && (coord1 != coord2);
}

inline bool
TransformedTriangle::testSegmentIntersectsHPlane(const TriSegment seg) const
{
    // get the H - coordinates
    const double coord1 = _coords[5 * seg + 4];
    const double coord2 = _coords[5 * ((seg + 1) % 3) + 4];
    //? should we use epsilon-equality here in second test?
    LOG(5, "coord1 : " << coord1 << " coord2 : " << coord2);

    return (coord1 * coord2 <= 0.0) && (coord1 != coord2);
}

inline bool
TransformedTriangle::testSurfaceAboveCorner(const TetraCorner corner) const
{
    // There is an error in Grandy -> it should be C_XY instead of C_YZ in [28].
    //
    // Idea: the nz value (Grandy [28] corrected!) can be interpreted as a special variant of the triple product t_O
    // where the z coordinate has been set to 1. It represents the signed volume of the tet OP'Q'R' where P', Q' and R'
    // are the projection of P, Q, R on the z=1 plane. Comparing the sign of this triple product with t_X (or t_Y, ...
    // dep on the corner) indicates whether the corner is in the half-space above or below the PQR triangle, similarly
    // to what is explained in [16]. (this works even for the Z corner, since we could have chosen z=24 (instead of z=1)
    // this would not have changed the final sign test).
    const double nz = calcStableC(PQ, C_XY) + calcStableC(QR, C_XY) + calcStableC(RP, C_XY);

    // Triple product might have not been computed, but here we need one:
    const double tp = _validTP[corner] ? calcStableT(corner) : calcTByDevelopingRow(corner, 1, false);

    return tp * nz >= 0.0;
}

inline bool
TransformedTriangle::testTriangleSurroundsRay(const TetraCorner corner) const
{
    // double products to use for the possible corners
    static const DoubleProduct DP_FOR_RAY_INTERSECTION[4] = {
        DoubleProduct(0),  // O - only here to fill out and make indices match
        C_10,              // X
        C_01,              // Y
        C_XY               // Z
    };

    const DoubleProduct dp = DP_FOR_RAY_INTERSECTION[corner];

    const double cPQ = calcStableC(PQ, dp);
    const double cQR = calcStableC(QR, dp);
    const double cRP = calcStableC(RP, dp);

    return (cPQ * cQR > 0.0) && (cPQ * cRP > 0.0);
}

inline const std::string &
TransformedTriangle::strTC(TetraCorner tc) const
{
    static const std::map<TetraCorner, std::string> m = {{O, "O"}, {X, "X"}, {Y, "Y"}, {Z, "Z"}};
    return m.at(tc);
}

inline const std::string &
TransformedTriangle::strTE(TetraEdge te) const
{
    static const std::map<TetraEdge, std::string> m = {
        {OX, "OX"}, {OY, "OY"}, {OZ, "OZ"}, {XY, "XY"}, {YZ, "YZ"}, {ZX, "ZX"}, {H01, "H01"}, {H10, "H10"}
    };
    return m.at(te);
}

inline const std::string &
TransformedTriangle::strTF(TetraFacet tf) const
{
    static const std::map<TetraFacet, std::string> m = {{OYZ, "OYZ"}, {OZX, "OZX"}, {OXY, "OXY"}, {XYZ, "XYZ"}};
    return m.at(tf);
}

inline const std::string &
TransformedTriangle::strTriC(TriCorner tc) const
{
    static const std::map<TriCorner, std::string> m = {{P, "P"}, {Q, "Q"}, {R, "R"}};
    return m.at(tc);
}

inline const std::string &
TransformedTriangle::strTriS(TriSegment ts) const
{
    static const std::map<TriSegment, std::string> m = {{PQ, "PQ"}, {QR, "QR"}, {RP, "RP"}};
    return m.at(ts);
}

}  // namespace INTERP_KERNEL

#endif
