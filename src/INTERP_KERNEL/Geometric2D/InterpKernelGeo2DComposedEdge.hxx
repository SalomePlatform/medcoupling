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

#ifndef __INTERPKERNELGEO2DCOMPOSEDNODE_HXX__
#define __INTERPKERNELGEO2DCOMPOSEDNODE_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelGeo2DEdge.hxx"

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

  /**
     * A set of quadratic or linear edges, described mainly by their connectivity
     * The set is assumed to be connected, but not necessarily closed (i.e. not necessarily forming a closed polygon).
     * Some methods however requires a closed form.
  */
  class ComposedEdge
  {
    friend class IteratorOnComposedEdge;
  public:
    INTERPKERNEL_EXPORT ComposedEdge() { }
    INTERPKERNEL_EXPORT ComposedEdge(const ComposedEdge& other);
    INTERPKERNEL_EXPORT ComposedEdge(int sz):_sub_edges(sz) { }
    INTERPKERNEL_EXPORT static void Delete(ComposedEdge *pt) { delete pt; }
    INTERPKERNEL_EXPORT static void SoftDelete(ComposedEdge *pt) { pt->_sub_edges.clear(); delete pt; }
    INTERPKERNEL_EXPORT void reverse();
    INTERPKERNEL_EXPORT int recursiveSize() const { return (int)_sub_edges.size(); }
    INTERPKERNEL_EXPORT bool presenceOfOn() const;
    INTERPKERNEL_EXPORT bool presenceOfQuadraticEdge() const;
    INTERPKERNEL_EXPORT void initLocations() const;
    INTERPKERNEL_EXPORT static void InitLocationsWithOther(const ComposedEdge& first, const ComposedEdge& other);
    INTERPKERNEL_EXPORT ComposedEdge *clone() const;
    INTERPKERNEL_EXPORT bool isNodeIn(Node *n) const;
    INTERPKERNEL_EXPORT double getArea() const;
    INTERPKERNEL_EXPORT double getPerimeter() const;
    INTERPKERNEL_EXPORT double getHydraulicDiameter() const;
    INTERPKERNEL_EXPORT void getBarycenter(double *bary) const;
    INTERPKERNEL_EXPORT void getBarycenterGeneral(double *bary) const;
    INTERPKERNEL_EXPORT double normalizeMe(double& xBary, double& yBary);
    INTERPKERNEL_EXPORT double normalize(ComposedEdge *other, double& xBary, double& yBary);
    INTERPKERNEL_EXPORT double normalizeExt(ComposedEdge *other, double& xBary, double& yBary);
    INTERPKERNEL_EXPORT void unApplyGlobalSimilarityExt(ComposedEdge& other, double xBary, double yBary, double fact);
    INTERPKERNEL_EXPORT void fillBounds(Bounds& output) const;
    INTERPKERNEL_EXPORT void applySimilarity(double xBary, double yBary, double dimChar);
    INTERPKERNEL_EXPORT void applyGlobalSimilarity(double xBary, double yBary, double dimChar);
    INTERPKERNEL_EXPORT void applyGlobalSimilarity2(ComposedEdge *other, double xBary, double yBary, double dimChar);
    INTERPKERNEL_EXPORT void dispatchPerimeter(double& partConsidered) const;
    INTERPKERNEL_EXPORT void dispatchPerimeterExcl(double& partConsidered, double& commonPart) const;
    INTERPKERNEL_EXPORT double dispatchPerimeterAdv(const ComposedEdge& father, std::vector<double>& result) const;
    INTERPKERNEL_EXPORT void getAllNodes(std::set<Node *>& output) const;
    INTERPKERNEL_EXPORT void initNodeHitStatus() const;
    INTERPKERNEL_EXPORT void applySimilarityOnMyNodes(double xBary, double yBary, double dimChar) const;
    INTERPKERNEL_EXPORT void unApplySimilarityOnMyNodes(double xBary, double yBary, double dimChar) const;
    INTERPKERNEL_EXPORT void applySimilarityOnMyNodesIfNotAlreadyHit(double xBary, double yBary, double dimChar) const;
    INTERPKERNEL_EXPORT void unApplySimilarityOnMyNodesIfNotAlreadyHit(double xBary, double yBary, double dimChar) const;
    INTERPKERNEL_EXPORT void initEdgeHitStatus() const;
    INTERPKERNEL_EXPORT void applySimilarityOnMyEdges(double xBary, double yBary, double dimChar) const;
    INTERPKERNEL_EXPORT void unApplySimilarityOnMyEdges(double xBary, double yBary, double dimChar) const;
    INTERPKERNEL_EXPORT void applySimilarityOnMyEdgesIfNotAlreadyHit(double xBary, double yBary, double dimChar) const;
    INTERPKERNEL_EXPORT void unApplySimilarityOnMyEdgesIfNotAlreadyHit(double xBary, double yBary, double dimChar) const;
    INTERPKERNEL_EXPORT void getBarycenter(double *bary, double& weigh) const;
    INTERPKERNEL_EXPORT bool completed() const { return getEndNode()==getStartNode(); }
    INTERPKERNEL_EXPORT void setValueAt(int i, Edge *e, bool direction=true);
    INTERPKERNEL_EXPORT double getCommonLengthWith(const ComposedEdge& other) const;
    INTERPKERNEL_EXPORT void clear();
    INTERPKERNEL_EXPORT bool empty() const { return _sub_edges.empty(); }
    INTERPKERNEL_EXPORT ElementaryEdge *front() const { return _sub_edges.front(); }
    INTERPKERNEL_EXPORT ElementaryEdge *back() const { return _sub_edges.back(); }
    INTERPKERNEL_EXPORT void resize(int i) { _sub_edges.resize(i); }
    INTERPKERNEL_EXPORT void pushBack(Edge *edge, bool direction=true);
    INTERPKERNEL_EXPORT void pushBack(ElementaryEdge *elem);
    INTERPKERNEL_EXPORT void pushBack(ComposedEdge *elem);
    INTERPKERNEL_EXPORT int size() const { return (int)_sub_edges.size(); }
    INTERPKERNEL_EXPORT ElementaryEdge *operator[](int i) const;
    INTERPKERNEL_EXPORT Node *getEndNode() const;
    INTERPKERNEL_EXPORT Node *getStartNode() const;
    INTERPKERNEL_EXPORT bool changeEndNodeWith(Node *node) const;
    INTERPKERNEL_EXPORT bool changeStartNodeWith(Node *node) const;
    INTERPKERNEL_EXPORT void dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const;
    INTERPKERNEL_EXPORT bool isInOrOut(Node *nodeToTest) const;
    INTERPKERNEL_EXPORT bool isInOrOut2(Node *nodeToTest) const;
    INTERPKERNEL_EXPORT bool getDirection() const;
    INTERPKERNEL_EXPORT bool intresincEqCoarse(const Edge *other) const;
  private:
    std::list<ElementaryEdge *>* getListBehind() { return &_sub_edges; }
    double isInOrOutAlg(Node *nodeToTest, const std::set<Node*>& nodes, std::set< IntersectElement >& inOutSwitch) const;
  protected:
    ~ComposedEdge();
  private:
    void clearAll(std::list<ElementaryEdge *>::iterator startToDel);
  protected:
    std::list<ElementaryEdge *> _sub_edges;
  };
}

#endif
