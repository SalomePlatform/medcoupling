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

#ifndef __INTERPKERNELGEO2DELEMENTARYEDGE_HXX__
#define __INTERPKERNELGEO2DELEMENTARYEDGE_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelException.hxx"
#include "InterpKernelGeo2DAbstractEdge.hxx"
#include "InterpKernelGeo2DEdge.hxx"

namespace INTERP_KERNEL
{
  class ElementaryEdge
  {
  public:
    INTERPKERNEL_EXPORT ElementaryEdge(Edge *ptr, bool direction):_direction(direction),_ptr(ptr) { }
    INTERPKERNEL_EXPORT ElementaryEdge(const ElementaryEdge& other);
    INTERPKERNEL_EXPORT ~ElementaryEdge();
    INTERPKERNEL_EXPORT bool isThereStartPoint() const { return _iterator.isValid(); }
    INTERPKERNEL_EXPORT IteratorOnComposedEdge& getIterator() { return _iterator; }
    INTERPKERNEL_EXPORT bool completed() const { return false; }
    INTERPKERNEL_EXPORT void declareOn() const { _ptr->declareOn(); }
    INTERPKERNEL_EXPORT void declareIn() const { _ptr->declareIn(); }
    INTERPKERNEL_EXPORT void declareOut() const { _ptr->declareOut(); }
    INTERPKERNEL_EXPORT TypeOfEdgeLocInPolygon getLoc() const { return _ptr->getLoc(); }
    INTERPKERNEL_EXPORT Edge *getPtr() const { return _ptr; }
    INTERPKERNEL_EXPORT void reverse() { _direction=(!_direction); }
    INTERPKERNEL_EXPORT bool isNodeIn(Node *n) const;
    INTERPKERNEL_EXPORT double getAreaOfZone() const { double ret=_ptr->getAreaOfZone(); return _direction?ret:-ret; }
    INTERPKERNEL_EXPORT void getBarycenterOfZone(double *bary) const;
    INTERPKERNEL_EXPORT void fillBounds(Bounds& output) const;
    INTERPKERNEL_EXPORT void applySimilarity(double xBary, double yBary, double dimChar) { _ptr->applySimilarity(xBary,yBary,dimChar); }
    INTERPKERNEL_EXPORT void unApplySimilarity(double xBary, double yBary, double dimChar) { _ptr->unApplySimilarity(xBary,yBary,dimChar); }
    INTERPKERNEL_EXPORT void getAllNodes(std::set<Node *>& output) const;
    INTERPKERNEL_EXPORT void getBarycenter(double *bary, double& weigh) const;
    INTERPKERNEL_EXPORT ElementaryEdge *clone() const;
    INTERPKERNEL_EXPORT void initLocations() const;
    INTERPKERNEL_EXPORT int size() const;
    INTERPKERNEL_EXPORT TypeOfEdgeLocInPolygon locateFullyMySelfAbsolute(const ComposedEdge& pol) const;
    INTERPKERNEL_EXPORT TypeOfEdgeLocInPolygon locateFullyMySelf(const ComposedEdge& pol, TypeOfEdgeLocInPolygon precEdgeLoc) const;
    INTERPKERNEL_EXPORT Node *getEndNode() const;
    INTERPKERNEL_EXPORT Node *getStartNode() const;
    INTERPKERNEL_EXPORT double getCurveLength() const { return _ptr->getCurveLength(); }
    INTERPKERNEL_EXPORT bool changeEndNodeWith(Node *node) const;
    INTERPKERNEL_EXPORT bool changeStartNodeWith(Node *node) const;
    INTERPKERNEL_EXPORT bool intresicEqual(const ElementaryEdge *other) const;
    INTERPKERNEL_EXPORT bool intresicEqualDirSensitive(const ElementaryEdge *other) const;
    INTERPKERNEL_EXPORT void dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const;
    INTERPKERNEL_EXPORT bool getDirection() const { return _direction; }
    INTERPKERNEL_EXPORT bool intresincEqCoarse(const Edge *other) const;
    INTERPKERNEL_EXPORT bool isEqual(const ElementaryEdge& other) const;
  public:
    INTERPKERNEL_EXPORT void fillGlobalInfoAbs(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                                               std::vector<int>& edgesThis, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int> mapAddCoo) const;
    INTERPKERNEL_EXPORT void fillGlobalInfoAbs2(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                                                std::vector<int>& edgesOther, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int>& mapAddCoo) const;
    INTERPKERNEL_EXPORT static ElementaryEdge *BuildEdgeFromStartEndDir(bool direction, INTERP_KERNEL::Node *start, INTERP_KERNEL::Node *end);
  private:
    bool _direction;
    Edge *_ptr;
    IteratorOnComposedEdge _iterator;
  };
}

#endif
