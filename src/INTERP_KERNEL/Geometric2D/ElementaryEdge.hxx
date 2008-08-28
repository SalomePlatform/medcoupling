#ifndef __ELEMENTARYEDGE_HXX__
#define __ELEMENTARYEDGE_HXX__

#include "Geometric2D_defines.hxx"

#include "AbstractEdge.hxx"
#include "InterpolationUtils.hxx"
#include "Edge.hxx"

namespace INTERP_KERNEL
{
  class GEOMETRIC2D_EXPORT ElementaryEdge : public AbstractEdge
  {
  public:
    ElementaryEdge(Edge *ptr, bool direction):_direction(direction),_ptr(ptr) { }
    ElementaryEdge(const ElementaryEdge& other);
    ~ElementaryEdge();
    bool completed() const { return false; }
    void declareOn() const { _ptr->declareOn(); }
    void declareIn() const { _ptr->declareIn(); }
    void declareOut() const { _ptr->declareOut(); }
    TypeOfEdgeLocInPolygon getLoc() const { return _ptr->getLoc(); }
    Edge *getPtr() const { return _ptr; }
    ElementaryEdge* &getLastElementary(IteratorOnComposedEdge::ItOnFixdLev &delta)  { throw Exception("Invalid call to getLastElementary"); }
    ElementaryEdge * &getFirstElementary(IteratorOnComposedEdge::ItOnFixdLev &delta)  { throw Exception("Invalid call to getFirstElementary"); }
    void reverse() { _direction=(!_direction); }
    bool isNodeIn(Node *n) const;
    double getAreaOfZone() const { return getAreaOfZoneFast(); }
    void fillBounds(Bounds& output) const;
    void getAllNodes(std::set<Node *>& output) const;
    void getBarycenter(double *bary, double& weigh) const;
    double getAreaOfZoneFast() const { double ret=_ptr->getAreaOfZone(); return _direction?ret:-ret; }
    AbstractEdge *clone() const;
    int recursiveSize() const;
    int size() const;
    TypeOfEdgeLocInPolygon locateFullyMySelfAbsolute(const ComposedEdge& pol) const;
    TypeOfEdgeLocInPolygon locateFullyMySelf(const ComposedEdge& pol, TypeOfEdgeLocInPolygon precEdgeLoc) const;
    AbstractEdge *&operator[](IteratorOnComposedEdge::ItOnFixdLev i);
    const AbstractEdge *&operator[](IteratorOnComposedEdge::ItOnFixdLev i) const;
    Node *getEndNode() const;
    Node *getStartNode() const;
    double getCurveLength() const { return _ptr->getCurveLength(); }
    bool changeEndNodeWith(Node *node) const;
    bool changeStartNodeWith(Node *node) const;
    bool intresicEqual(const AbstractEdge *other) const;
    bool intresicEqualDirSensitive(const AbstractEdge *other) const;
    void dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const;
    bool getDirection() const { return _direction; }
    bool intresincEqCoarse(const Edge *other) const;
  private:
    bool _direction;
    Edge *_ptr;
  };
}

#endif
