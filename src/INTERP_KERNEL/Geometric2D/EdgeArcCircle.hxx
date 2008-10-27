#ifndef __EDGEARCCIRCLE_HXX__
#define __EDGEARCCIRCLE_HXX__

#include "Geometric2D_defines.hxx"

#include "Edge.hxx"

namespace INTERP_KERNEL
{
  class GEOMETRIC2D_EXPORT ArcCArcCIntersector : public SameTypeIntersector
  {
    friend class EdgeArcCircle;
  public:
    ArcCArcCIntersector(const EdgeArcCircle& e1, const EdgeArcCircle& e2);
    bool haveTheySameDirection() const;
    void getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const;
    void areOverlappedOrOnlyColinears(const Bounds *whereToFind, bool& obviousNoIntersection, bool& areOverlapped);
    std::list< IntersectElement > getIntersectionsCharacteristicVal() const;
  private:
    //! return angle in ]-Pi;Pi[ - 'node' must be on curve of '_e1'
    double getAngle(Node *node) const;
    static bool areArcsOverlapped(const EdgeArcCircle& a1, const EdgeArcCircle& a2);
    static bool isIn2Pi(double start, double delta, double angleIn);
    //! 'delta' 'start' in ]-Pi;Pi[
    static bool isAngleNotIn(double start, double delta, double angleIn);
    //! for an angle 'angle' in ]-3*Pi;3*Pi[ returns angle in ]-Pi;Pi[
    static double normalizeAngle(double angle) { if(angle> M_PI) return angle-2.*M_PI; if(angle<-M_PI) return angle+2.*M_PI; return angle; }
  private:
    const EdgeArcCircle& getE1() const { return (const EdgeArcCircle&)_e1; }
    const EdgeArcCircle& getE2() const { return (const EdgeArcCircle&)_e2; }
  private:
    double _dist;
  };

  class ArcCSegIntersector : public CrossTypeIntersector
  {
  public:
    ArcCSegIntersector(const EdgeArcCircle& e1, const EdgeLin& e2, bool reverse=true);
    //virtual overloading
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
  
  class EdgeArcCircle : public Edge
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
    bool isIn(double characterVal) const;
    Node *buildRepresentantOfMySelf() const;
    bool isLower(double val1, double val2) const;
    double getCharactValue(const Node& node) const;
    double getDistanceToPoint(const double *pt) const;
    bool isNodeLyingOn(const double *coordOfNode) const;
    TypeOfFunction getTypeOfFunc() const { return ARC_CIRCLE; }
    void dynCastFunction(const EdgeLin * &seg,
                         const EdgeArcCircle * &arcSeg) const { arcSeg=this; }
    const double *getCenter() const { return _center; }
    void getCenter(double *center) const { center[0]=_center[0]; center[1]=_center[1]; }
    bool doIHaveSameDirectionAs(const Edge& other) const { return false; }
    void applySimilarity(double xBary, double yBary, double dimChar);
    double getAngle0() const { return _angle0; }
    double getRadius() const { return _radius; }
    double getAngle() const { return _angle; }
    static double getAbsoluteAngle(const double *vect, double& normVect);
    static void getArcOfCirclePassingThru(const double *start, const double *middle, const double *end, 
                                          double *center, double& radius, double& angleInRad, double& angleInRad0);
    //! To avoid in aggressive optimizations nan.
    static double safeSqrt(double val) { double ret=fmax(val,0.); return sqrt(ret); }
    static double safeAcos(double cosAngle) { double ret=fmin(cosAngle,1.); ret=fmax(ret,-1.); return acos(ret); }
    static double safeAsin(double sinAngle) { double ret=fmin(sinAngle,1.); ret=fmax(ret,-1.); return asin(ret); }
  protected:
    void updateBounds();
    Edge *buildEdgeLyingOnMe(Node *start, Node *end, bool direction=true) const;
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
