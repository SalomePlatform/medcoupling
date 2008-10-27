#ifndef __COMPOSEDNODE_HXX__
#define __COMPOSEDNODE_HXX__

#include "Geometric2D_defines.hxx"
#include <set>
#include <list>
#include <vector>

namespace INTERP_KERNEL
{
  class Node;
  class Edge;
  class Bounds;
  class ElementaryEdge;
  class IteratorOnComposedEdge;

  class GEOMETRIC2D_EXPORT ComposedEdge
  {
    friend class IteratorOnComposedEdge;
  public:
    ComposedEdge() { }
    ComposedEdge(const ComposedEdge& other);
    ComposedEdge(int size):_subEdges(size) { }
    static void Delete(ComposedEdge *pt) { delete pt; }
    static void SoftDelete(ComposedEdge *pt) { pt->_subEdges.clear(); delete pt; }
    void reverse();
    int recursiveSize() const { return _subEdges.size(); }
    void initLocations() const;
    ComposedEdge *clone() const;
    bool isNodeIn(Node *n) const;
    double getArea() const;
    double getPerimeter() const;
    double getHydraulicDiameter() const;
    double normalize(ComposedEdge *other);
    void fillBounds(Bounds& output) const;
    int getNbOfEdgeSonsOfSameId(int id) const;
    void applySimilarity(double xBary, double yBary, double dimChar);
    void applyGlobalSimilarity(double xBary, double yBary, double dimChar);
    void dispatchForNode(const ComposedEdge& father, std::vector<int>& nbOfCreatedNodes) const;
    void dispatchPerimeter(const std::set<int>& ids1, const std::set<int>& ids2,
                           double& part1, double& part2, double& commonPart) const;
    double dispatchPerimeterAdv(const ComposedEdge& father, std::vector<double>& result) const;
    double checkFatherHood(ElementaryEdge *edge, std::vector<double>& result) const;
    void fillAllEdgeIds(std::set<int>& ids) const;
    void getAllNodes(std::set<Node *>& output) const;
    void getBarycenter(double *bary, double& weigh) const;
    bool completed() const { return getEndNode()==getStartNode(); }
    void setValueAt(int i, Edge *e, bool direction=true);
    double getCommonLengthWith(const ComposedEdge& other) const;
    void clear();
    bool empty() const { return _subEdges.empty(); }
    ElementaryEdge *front() const { return _subEdges.front(); }
    ElementaryEdge *back() const { return _subEdges.back(); }
    void resize(int i) { _subEdges.resize(i); }
    void pushBack(Edge *edge, bool direction=true);
    void pushBack(ElementaryEdge *elem);
    void pushBack(ComposedEdge *elem);
    int size() const { return _subEdges.size(); }
    ElementaryEdge *operator[](int i) const;
    Node *getEndNode() const;
    Node *getStartNode() const;
    bool changeEndNodeWith(Node *node) const;
    bool changeStartNodeWith(Node *node) const;
    void dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const;
    bool isInOrOut(Node *nodeToTest) const;
    bool getDirection() const;
    bool intresincEqCoarse(const Edge *other) const;
  private:
    std::list<ElementaryEdge *>* getListBehind() { return &_subEdges; }
  protected:
    ~ComposedEdge();
  private:
    void clearAll(std::list<ElementaryEdge *>::iterator startToDel);
  protected:
    std::list<ElementaryEdge *> _subEdges;
  };
}

#endif
