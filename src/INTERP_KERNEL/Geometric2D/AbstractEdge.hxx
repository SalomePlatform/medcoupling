#ifndef __ABSTRACTEDGE_HXX__
#define __ABSTRACTEDGE_HXX__

#include "Geometric2D_defines.hxx"

#include <set>
#include <fstream>

namespace INTERP_KERNEL
{
  class Edge;
  class Node;
  class Bounds;

  class AbstractEdge;
  class ComposedEdge;
  class ElementaryEdge;

  class GEOMETRIC2D_EXPORT IteratorOnComposedEdge
  {
    friend class AbstractEdge;
    friend class ComposedEdge;
    friend class ElementaryEdge;
    friend class QuadraticPolygon;
    //! Implicitely we suppose here that at maximum we have 256 edges on a current level.
    typedef unsigned char ItOnFixdLev;
  public:
    IteratorOnComposedEdge(ComposedEdge *cont);
    void operator=(const IteratorOnComposedEdge& other);
    void first();
    void next();
    void last();
    void nextLoop();
    void previousLoop();
    bool finished() const;
    AbstractEdge *getLowestDealing() const;
    bool goToNextInOn(bool direction, int& i, int nbMax);
    ElementaryEdge * &current() { return updateNumbering(); }
    AbstractEdge *currentDirect() const;
  private:
    ElementaryEdge* &updateNumbering();
  private:
    //! this number (+1) represents the maximum intersection an edge is going to have.
    static const unsigned MAX_INTERSCT_DEPH=8;
    ComposedEdge *_container;
    mutable ItOnFixdLev _current[MAX_INTERSCT_DEPH];
  };

  class GEOMETRIC2D_EXPORT AbstractEdge
  {
  public:
    virtual ~AbstractEdge() { }
    virtual int size() const = 0;
    virtual void reverse() = 0;
    virtual bool completed() const = 0;
    virtual int recursiveSize() const = 0;
    virtual AbstractEdge *clone() const = 0;
    virtual bool isNodeIn(Node *n) const = 0;
    virtual double getAreaOfZone() const = 0;
    virtual void fillBounds(Bounds& output) const = 0;
    virtual void getAllNodes(std::set<Node *>& output) const = 0;
    virtual void getBarycenter(double *bary, double& weigh) const = 0;
    virtual ElementaryEdge* &getLastElementary(IteratorOnComposedEdge::ItOnFixdLev &delta) = 0;
    virtual ElementaryEdge* &getFirstElementary(IteratorOnComposedEdge::ItOnFixdLev &delta) = 0;
    virtual const AbstractEdge *&operator[](IteratorOnComposedEdge::ItOnFixdLev i) const = 0;
    virtual AbstractEdge *&operator[](IteratorOnComposedEdge::ItOnFixdLev i) = 0;
    virtual Node *getEndNode() const = 0;
    virtual Node *getStartNode() const = 0;
    virtual bool changeStartNodeWith(Node *node) const = 0;
    virtual bool changeEndNodeWith(Node *node) const = 0;
    virtual bool intresicEqual(const AbstractEdge *other) const = 0;
    virtual bool intresicEqualDirSensitive(const AbstractEdge *other) const = 0;
    virtual void dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const = 0;
    virtual bool intresincEqCoarse(const Edge *other) const = 0;
    virtual bool getDirection() const = 0;
  };
}

#endif
