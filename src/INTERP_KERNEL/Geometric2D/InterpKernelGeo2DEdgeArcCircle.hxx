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
// Author : Anthony Geay (CEA/DEN)

#ifndef __INTERPKERNELGEO2DEDGEARCCIRCLE_HXX__
#define __INTERPKERNELGEO2DEDGEARCCIRCLE_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelGeo2DEdge.hxx"

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT ArcCArcCIntersector : public SameTypeEdgeIntersector
  {
  public:
    ArcCArcCIntersector(const EdgeArcCircle& e1, const EdgeArcCircle& e2);
    bool haveTheySameDirection() const;
    bool areColinears() const;
    void getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const;
    void areOverlappedOrOnlyColinears(const Bounds *whereToFind, bool& obviousNoIntersection, bool& areOverlapped);
    std::list< IntersectElement > getIntersectionsCharacteristicVal() const;
  private:
    //! return angle in ]-Pi;Pi[ - 'node' must be on curve of '_e1'
    double getAngle(Node *node) const;
    static bool internalAreColinears(const EdgeArcCircle& a1, const EdgeArcCircle& a2, double& distBetweenCenters, double& cst, double& radiusL, double centerL[2], double& raduisB, double centerB[2]);
    static bool areArcsOverlapped(const EdgeArcCircle& a1, const EdgeArcCircle& a2);
  private:
    const EdgeArcCircle& getE1() const { return (const EdgeArcCircle&)_e1; }
    const EdgeArcCircle& getE2() const { return (const EdgeArcCircle&)_e2; }
  private:
    double _dist;
  };

  class INTERPKERNEL_EXPORT ArcCSegIntersector : public CrossTypeEdgeIntersector
  {
  public:
    ArcCSegIntersector(const EdgeArcCircle& e1, const EdgeLin& e2, bool reverse=true);
    //virtual overloading
    bool areColinears() const;
    void getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const;
    void areOverlappedOrOnlyColinears(const Bounds *whereToFind, bool& obviousNoIntersection, bool& areOverlapped);
    std::list< IntersectElement > getIntersectionsCharacteristicVal() const;
  private:
    const EdgeArcCircle& getE1() const { return (const EdgeArcCircle&)_e1; }
    const EdgeLin& getE2() const { return (const EdgeLin&)_e2; }
  private:
    double _dx;
    double _dy;
    double _drSq;
    double _cross;
    double _determinant;
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
    void dumpInXfigFile(std::ostream& stream, bool direction, int resolution, const Bounds& box) const;
    void update(Node *m);
    double getAreaOfZone() const;
    double getCurveLength() const;
    void getBarycenter(double *bary) const;
    void getBarycenterOfZone(double *bary) const;
    void getMiddleOfPoints(const double *p1, const double *p2, double *mid) const;
    void getMiddleOfPointsOriented(const double *p1, const double *p2, double *mid) const;
    bool isIn(double characterVal) const;
    Node *buildRepresentantOfMySelf() const;
    bool isLower(double val1, double val2) const;
    double getCharactValue(const Node& node) const;
    double getCharactValueBtw0And1(const Node& node) const;
    double getDistanceToPoint(const double *pt) const;
    bool isNodeLyingOn(const double *coordOfNode) const;
    TypeOfFunction getTypeOfFunc() const { return ARC_CIRCLE; }
    void dynCastFunction(const EdgeLin * &seg,
                         const EdgeArcCircle * &arcSeg) const { arcSeg=this; }
    const double *getCenter() const { return _center; }
    void getCenter(double *center) const { center[0]=_center[0]; center[1]=_center[1]; }
    bool doIHaveSameDirectionAs(const Edge& other) const { return false; }
    void applySimilarity(double xBary, double yBary, double dimChar);
    void unApplySimilarity(double xBary, double yBary, double dimChar);
    double getAngle0() const { return _angle0; }
    double getRadius() const { return _radius; }
    double getAngle() const { return _angle; }
    void tesselate(const int *conn, int offset, double eps, std::vector<int>& newConn, std::vector<double>& addCoo) const;
    static EdgeArcCircle *BuildFromNodes(Node *start, Node *middle, Node *end);
    static double GetAbsoluteAngle(const double *vect, double& normVect);
    static double GetAbsoluteAngleOfNormalizedVect(double ux, double uy);
    static void GetArcOfCirclePassingThru(const double *start, const double *middle, const double *end, 
                                          double *center, double& radius, double& angleInRad, double& angleInRad0);
    //! To avoid in aggressive optimizations nan.
    static double SafeSqrt(double val) { double ret=std::max(val,0.); return sqrt(ret); }
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
    Edge *buildEdgeLyingOnMe(Node *start, Node *end, bool direction=true) const;
    void fillGlobalInfoAbs(bool direction, const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                           std::vector<int>& edgesThis, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int> mapAddCoo) const;
    void fillGlobalInfoAbs2(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                            std::vector<int>& edgesOther, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int>& mapAddCoo) const;
  protected:
    //!Value between -2Pi and 2Pi
    double _angle;
    //!Value between -Pi and Pi
    double _angle0;
    double _radius;
    double _center[2];
  };
}

#endif
