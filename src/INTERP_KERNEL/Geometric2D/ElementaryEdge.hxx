#ifndef __ELEMENTARYEDGE_HXX__
#define __ELEMENTARYEDGE_HXX__

#include "Geometric2D_defines.hxx"
#include "InterpolationUtils.hxx"
#include "AbstractEdge.hxx"
#include "Edge.hxx"

namespace INTERP_KERNEL
{
  class GEOMETRIC2D_EXPORT ElementaryEdge
  {
  public:
    ElementaryEdge(Edge *ptr, bool direction):_direction(direction),_ptr(ptr) { }
    ElementaryEdge(const ElementaryEdge& other);
    ~ElementaryEdge();
    bool isThereStartPoint() const { return _iterator.isValid(); }
    IteratorOnComposedEdge& getIterator() { return _iterator; }
    bool completed() const { return false; }
    void declareOn() const { _ptr->declareOn(); }
    void declareIn() const { _ptr->declareIn(); }
    void declareOut() const { _ptr->declareOut(); }
    TypeOfEdgeLocInPolygon getLoc() const { return _ptr->getLoc(); }
    Edge *getPtr() const { return _ptr; }
    void reverse() { _direction=(!_direction); }
    bool isNodeIn(Node *n) const;
    double getAreaOfZone() const { double ret=_ptr->getAreaOfZone(); return _direction?ret:-ret; }
    void fillBounds(Bounds& output) const;
    void applySimilarity(double xBary, double yBary, double dimChar) { _ptr->applySimilarity(xBary,yBary,dimChar); }
    void getAllNodes(std::set<Node *>& output) const;
    void getBarycenter(double *bary, double& weigh) const;
    ElementaryEdge *clone() const;
    void initLocations() const;
    int size() const;
    TypeOfEdgeLocInPolygon locateFullyMySelfAbsolute(const ComposedEdge& pol) const;
    TypeOfEdgeLocInPolygon locateFullyMySelf(const ComposedEdge& pol, TypeOfEdgeLocInPolygon precEdgeLoc) const;
    Node *getEndNode() const;
    Node *getStartNode() const;
    double getCurveLength() const { return _ptr->getCurveLength(); }
    bool changeEndNodeWith(Node *node) const;
    bool changeStartNodeWith(Node *node) const;
    bool intresicEqual(const ElementaryEdge *other) const;
    bool intresicEqualDirSensitive(const ElementaryEdge *other) const;
    void dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const;
    bool getDirection() const { return _direction; }
    bool intresincEqCoarse(const Edge *other) const;
  private:
    bool _direction;
    Edge *_ptr;
    IteratorOnComposedEdge _iterator;
  };
}

#endif
