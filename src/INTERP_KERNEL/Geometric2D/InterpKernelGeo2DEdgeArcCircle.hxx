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
// Author : Anthony Geay (CEA/DEN)

#ifndef __INTERPKERNELGEO2DEDGEARCCIRCLE_HXX__
#define __INTERPKERNELGEO2DEDGEARCCIRCLE_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelGeo2DBounds.hxx"
#include "InterpKernelGeo2DEdge.hxx"
#include "InterpKernelGeo2DNode.hxx"
#include "MCIdType.hxx"
#include <list>
#include <istream>
#include <ostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <math.h>

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT ArcCArcCIntersector : public SameTypeEdgeIntersector
  {
  public:
    ArcCArcCIntersector(const EdgeArcCircle& e1, const EdgeArcCircle& e2);
    bool haveTheySameDirection() const override;
    bool areColinears() const override;
    void getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const override;
    void areOverlappedOrOnlyColinears(bool& obviousNoIntersection, bool& areOverlapped) override;
    std::list< IntersectElement > getIntersectionsCharacteristicVal() const override;
  private:
    //! return angle in ]-Pi;Pi[ - 'node' must be on curve of '_e1'
    double getAngle(Node *node) const;
    static bool internalAreColinears(const EdgeArcCircle& a1, const EdgeArcCircle& a2, double& distBetweenCenters, double& cst, double& radiusL, double centerL[2], double& raduisB, double centerB[2]);
    static bool areArcsOverlapped(const EdgeArcCircle& a1, const EdgeArcCircle& a2);
  private:
    const EdgeArcCircle& getE1() const { return (const EdgeArcCircle&)_e1; }
    const EdgeArcCircle& getE2() const { return (const EdgeArcCircle&)_e2; }
  private:
    double _dist;       // distance between the two arc centers
  };

  /**
   * Cross-type intersector: edge1 is the arc of circle, edge2 is the segment.
   */
  class INTERPKERNEL_EXPORT ArcCSegIntersector : public CrossTypeEdgeIntersector
  {
  public:
    ArcCSegIntersector(const EdgeArcCircle& e1, const EdgeLin& e2, bool reverse=true);
    //virtual overloading
    bool areColinears() const override;
    void getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const override;
    void areOverlappedOrOnlyColinears(bool& obviousNoIntersection, bool& areOverlapped) override;
    std::list< IntersectElement > getIntersectionsCharacteristicVal() const override;
  private:
    const EdgeArcCircle& getE1() const { return (const EdgeArcCircle&)_e1; }
    const EdgeLin& getE2() const { return (const EdgeLin&)_e2; }
  private:
    double _dx;           //!< X extent of the segment
    double _dy;           //!< Y extent of the segment
    double _drSq;         //!< Square of the norm of the seg
    double _cross;        //!< See areOverlappedOrOnlyColinears()
    double _deltaRoot_div_dr;    //!< See areOverlappedOrOnlyColinears()
    bool _i1S2E,_i1E2E;
  };

  class INTERPKERNEL_EXPORT EdgeArcCircle : public Edge
  {
  public:
    EdgeArcCircle(std::istream& lineInXfig);
    EdgeArcCircle(Node *start, Node *middle, Node *end, bool direction = true);
    EdgeArcCircle(double sX, double sY, double mX, double mY, double eX, double eY);
    EdgeArcCircle(Node *start, Node *end, const double *center, double radius, double angle0, double deltaAngle, bool direction=true);
    //! for tests
    void changeMiddle(Node *newMiddle);
    void dumpInXfigFile(std::ostream& stream, bool direction, int resolution, const Bounds& box) const override;
    void update(Node *m) override;
    double getAreaOfZone() const override;
    double getCurveLength() const override;
    void getBarycenter(double *bary) const override;
    void getBarycenterOfZone(double *bary) const override;
    void getMiddleOfPoints(const double *p1, const double *p2, double *mid) const override;
    void getMiddleOfPointsOriented(const double *p1, const double *p2, double *mid) const override;
    bool isIn(double characterVal) const override;
    Node *buildRepresentantOfMySelf() const override;
    bool isLower(double val1, double val2) const override;
    double getCharactValue(const Node& node) const override;
    double getCharactValueBtw0And1(const Node& node) const override;
    double getDistanceToPoint(const double *pt) const override;
    bool isNodeLyingOn(const double *coordOfNode) const override;
    TypeOfFunction getTypeOfFunc() const override { return ARC_CIRCLE; }
    void dynCastFunction(const EdgeLin * & /*seg*/,
                         const EdgeArcCircle * &arcSeg) const override { arcSeg=this; }
    const double *getCenter() const { return _center; }
    void getCenter(double *center) const { center[0]=_center[0]; center[1]=_center[1]; }
    bool doIHaveSameDirectionAs(const Edge&  /*other*/) const { return false; }
    void applySimilarity(double xBary, double yBary, double dimChar) override;
    void unApplySimilarity(double xBary, double yBary, double dimChar) override;
    double getAngle0() const { return _angle0; }
    double getRadius() const { return _radius; }
    double getAngle() const { return _angle; }
    void tesselate(const mcIdType *conn, mcIdType offset, double eps, std::vector<mcIdType>& newConn, std::vector<double>& addCoo) const;
    static EdgeArcCircle *BuildFromNodes(Node *start, Node *middle, Node *end);
    static double GetAbsoluteAngle(const double *vect, double& normVect);
    static double GetAbsoluteAngleOfNormalizedVect(double ux, double uy);
    static void GetArcOfCirclePassingThru(const double *start, const double *middle, const double *end, 
                                          double *center, double& radius, double& angleInRad, double& angleInRad0);
    //! To avoid in aggressive optimizations nan.
    static double SafeSqrt(double val) { double const ret=std::max(val,0.); return sqrt(ret); }
    static double SafeAcos(double cosAngle) { double ret=std::min(cosAngle,1.); ret=std::max(ret,-1.); return acos(ret); }
    static double SafeAsin(double sinAngle) { double ret=std::min(sinAngle,1.); ret=std::max(ret,-1.); return asin(ret); }
    //! @param start and @param angleIn in ]-Pi;Pi] and @param delta in ]-2*Pi,2*Pi[
    static bool IsIn2Pi(double start, double delta, double angleIn);
    //! 'delta' 'start' in ]-Pi;Pi[
    static bool IsAngleNotIn(double start, double delta, double angleIn);
    //! for an angle 'angle' in ]-3*Pi;3*Pi[ returns angle in ]-Pi;Pi[
    static double NormalizeAngle(double angle) { if(angle>M_PI) return angle-2.*M_PI; if(angle<-M_PI) return angle+2.*M_PI; return angle; }
  protected:
    void updateBounds();
    Edge *buildEdgeLyingOnMe(Node *start, Node *end, bool direction=true) const override;

  protected:
    //! Absolute angle where the arc starts. Value between -Pi and Pi
    double _angle0;
    //! Angular span of the arc. Value between -2Pi and 2Pi
    double _angle;
    //! Radius of the arc of circle
    double _radius;
    //! Center of the arc of circle
    double _center[2];
  };
}

#endif
