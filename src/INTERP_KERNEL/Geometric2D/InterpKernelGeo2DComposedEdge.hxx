//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef __INTERPKERNELGEO2DCOMPOSEDNODE_HXX__
#define __INTERPKERNELGEO2DCOMPOSEDNODE_HXX__

#include "INTERPKERNELGEOMETRIC2DDefines.hxx"

#include <set>
#include <list>
#include <vector>
#include <ostream>

namespace INTERP_KERNEL
{
  class Node;
  class Edge;
  class Bounds;
  class ElementaryEdge;
  class IteratorOnComposedEdge;

  class INTERPKERNELGEOMETRIC2D_EXPORT ComposedEdge
  {
    friend class IteratorOnComposedEdge;
  public:
    ComposedEdge() { }
    ComposedEdge(const ComposedEdge& other);
    ComposedEdge(int size):_sub_edges(size) { }
    static void Delete(ComposedEdge *pt) { delete pt; }
    static void SoftDelete(ComposedEdge *pt) { pt->_sub_edges.clear(); delete pt; }
    void reverse();
    int recursiveSize() const { return _sub_edges.size(); }
    void initLocations() const;
    ComposedEdge *clone() const;
    bool isNodeIn(Node *n) const;
    double getArea() const;
    double getPerimeter() const;
    double getHydraulicDiameter() const;
    void getBarycenter(double *bary) const;
    double normalize(ComposedEdge *other, double& xBary, double& yBary);
    void fillBounds(Bounds& output) const;
    void applySimilarity(double xBary, double yBary, double dimChar);
    void applyGlobalSimilarity(double xBary, double yBary, double dimChar);
    void dispatchPerimeter(double& partConsidered) const;
    void dispatchPerimeterExcl(double& partConsidered, double& commonPart) const;
    double dispatchPerimeterAdv(const ComposedEdge& father, std::vector<double>& result) const;
    void getAllNodes(std::set<Node *>& output) const;
    void getBarycenter(double *bary, double& weigh) const;
    bool completed() const { return getEndNode()==getStartNode(); }
    void setValueAt(int i, Edge *e, bool direction=true);
    double getCommonLengthWith(const ComposedEdge& other) const;
    void clear();
    bool empty() const { return _sub_edges.empty(); }
    ElementaryEdge *front() const { return _sub_edges.front(); }
    ElementaryEdge *back() const { return _sub_edges.back(); }
    void resize(int i) { _sub_edges.resize(i); }
    void pushBack(Edge *edge, bool direction=true);
    void pushBack(ElementaryEdge *elem);
    void pushBack(ComposedEdge *elem);
    int size() const { return _sub_edges.size(); }
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
    std::list<ElementaryEdge *>* getListBehind() { return &_sub_edges; }
  protected:
    ~ComposedEdge();
  private:
    void clearAll(std::list<ElementaryEdge *>::iterator startToDel);
  protected:
    std::list<ElementaryEdge *> _sub_edges;
  };
}

#endif
