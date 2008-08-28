#ifndef __COMPOSEDNODE_HXX__
#define __COMPOSEDNODE_HXX__

#include "Geometric2D_defines.hxx"

#include "AbstractEdge.hxx"

#include <vector>

namespace INTERP_KERNEL
{
  class GEOMETRIC2D_EXPORT ComposedEdge : public AbstractEdge
  {
  public:
    ComposedEdge() { }
    ComposedEdge(const ComposedEdge& other);
    ComposedEdge(int size):_subEdges(size) { }
    static void Delete(ComposedEdge *pt) { delete pt; }
    void reverse();
    int recursiveSize() const;
    AbstractEdge *clone() const;
    bool isNodeIn(Node *n) const;
    double getAreaOfZone() const;
    void fillBounds(Bounds& output) const;
    void getAllNodes(std::set<Node *>& output) const;
    void getBarycenter(double *bary, double& weigh) const;
    bool completed() const { return getEndNode()==getStartNode(); }
    ElementaryEdge * &getLastElementary(IteratorOnComposedEdge::ItOnFixdLev &delta);
    ElementaryEdge * &getFirstElementary(IteratorOnComposedEdge::ItOnFixdLev &delta);
    void setValueAt(int i, AbstractEdge *val) { delete _subEdges[i]; _subEdges[i]=val; }
    void setValueAt(int i, Edge *e, bool direction=true);
    void clear();
    bool empty() const { return _subEdges.empty(); }
    AbstractEdge *front() const { return _subEdges.front(); }
    AbstractEdge *back() const { return _subEdges.back(); }
    void resize(int i) { _subEdges.resize(i); }
    void pushBack(Edge *edge, bool direction=true);
    void pushBack(AbstractEdge *elem);
    int size() const { return _subEdges.size(); }
    AbstractEdge *&operator[](IteratorOnComposedEdge::ItOnFixdLev i) { return (AbstractEdge *&)_subEdges[i]; }
    const AbstractEdge *&operator[](IteratorOnComposedEdge::ItOnFixdLev i) const { return (const AbstractEdge *&)_subEdges[i]; }
    Node *getEndNode() const;
    Node *getStartNode() const;
    AbstractEdge *simplify();
    bool addEdgeIfIn(ElementaryEdge *edge);
    bool changeEndNodeWith(Node *node) const;
    bool changeStartNodeWith(Node *node) const;
    void dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const;
    bool intresicEqual(const AbstractEdge *other) const;
    bool intresicEqualDirSensitive(const AbstractEdge *other) const;
    bool isInOrOut(Node *nodeToTest) const;
    bool getDirection() const;
    bool intresincEqCoarse(const Edge *other) const;
  protected:
    ~ComposedEdge();
  private:
    void clearAll(std::vector<AbstractEdge *>::iterator startToDel);
  protected:
    std::vector<AbstractEdge *> _subEdges;
  };
}

#endif
