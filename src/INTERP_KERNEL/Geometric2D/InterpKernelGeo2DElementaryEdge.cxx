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

#include "InterpKernelGeo2DElementaryEdge.hxx"
#include "InterpKernelException.hxx"
#include "InterpKernelGeo2DEdge.hxx"
#include "InterpKernelGeo2DComposedEdge.hxx"

using namespace INTERP_KERNEL;

ElementaryEdge::ElementaryEdge(const ElementaryEdge& other):_direction(other._direction),_ptr(other._ptr)
{
  _ptr->incrRef(); 
}

ElementaryEdge::~ElementaryEdge()
{
  if(_ptr)
    _ptr->decrRef();
}

bool ElementaryEdge::isNodeIn(Node *n) const
{
  return _ptr->getStartNode()==n || _ptr->getEndNode()==n;
}

/*!
 * \b WARNING contrary to INTERP_KERNEL::Edge::getBarycenterOfZone method called,
 * this one is cumulative.
 */
void ElementaryEdge::getBarycenterOfZone(double *bary) const
{
  double tmp[2];
  _ptr->getBarycenterOfZone(tmp);
  if(_direction)
    {
      bary[0]+=tmp[0];
      bary[1]+=tmp[1];
    }
  else
    {
      bary[0]-=tmp[0];
      bary[1]-=tmp[1];
    }
}

void ElementaryEdge::fillBounds(Bounds& output) const
{
  output.aggregate(_ptr->getBounds());
}

void ElementaryEdge::getAllNodes(std::set<Node *>& output) const
{
  output.insert(_ptr->getStartNode());
  output.insert(_ptr->getEndNode());
}

void ElementaryEdge::getBarycenter(double *bary, double& weigh) const
{
  _ptr->getBarycenter(bary);
  weigh=_ptr->getCurveLength();
}

ElementaryEdge *ElementaryEdge::clone() const
{
  return new ElementaryEdge(*this);
}

void ElementaryEdge::initLocations() const
{
  _ptr->initLocs();
}

/*!
 * WARNING use this method if and only if this is so that it is completely in/out/on of @param pol.
 */
TypeOfEdgeLocInPolygon ElementaryEdge::locateFullyMySelf(const ComposedEdge& pol, TypeOfEdgeLocInPolygon precEdgeLoc) const
{
  if(getLoc()!=FULL_UNKNOWN)
    return getLoc();
  //obvious cases
  if(precEdgeLoc==FULL_IN_1)
    {
      if(getStartNode()->getLoc()==ON_1)
        {
          declareOut();
          return getLoc();
        }
      else if(getStartNode()->getLoc()==IN_1 || getStartNode()->getLoc()==ON_TANG_1)
        {
          declareIn();
          return getLoc();
        }
    }
  if(precEdgeLoc==FULL_OUT_1)
    {
      if(getStartNode()->getLoc()==ON_1)
        {
          declareIn();
          return getLoc();
        }
      else if(getStartNode()->getLoc()==IN_1 || getStartNode()->getLoc()==ON_TANG_1)
        {
          declareOut();
          return getLoc();
        }
    }
  if(getStartNode()->getLoc()==IN_1 || getEndNode()->getLoc()==IN_1)
    {
      declareIn();
      return getLoc();
    }
  if(getStartNode()->getLoc()==OUT_1 || getEndNode()->getLoc()==OUT_1)
    {
      declareOut();
      return getLoc();
    }
  //a seek is requested
  return locateFullyMySelfAbsolute(pol);
}

TypeOfEdgeLocInPolygon ElementaryEdge::locateFullyMySelfAbsolute(const ComposedEdge& pol) const
{
  Node *node=_ptr->buildRepresentantOfMySelf(); // build barycenter used to detect if the edge is IN or OUT
  if(pol.isInOrOut(node))
    declareIn(); 
  else
    declareOut();
  node->decrRef();
  return getLoc();
}

Node *ElementaryEdge::getEndNode() const
{ 
  if(_direction)
    return _ptr->getEndNode();
  else return _ptr->getStartNode();
}

Node *ElementaryEdge::getStartNode() const
{
  if(_direction)
    return _ptr->getStartNode();
  else 
    return _ptr->getEndNode();
}

bool ElementaryEdge::changeEndNodeWith(Node *node) const
{
  if(_direction)
    return _ptr->changeEndNodeWith(node);
  else 
    return _ptr->changeStartNodeWith(node);
}

bool ElementaryEdge::changeStartNodeWith(Node *node) const
{
  if(_direction)
    return _ptr->changeStartNodeWith(node);
  else 
    return _ptr->changeEndNodeWith(node);
}

void ElementaryEdge::dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const
{
  _ptr->dumpInXfigFile(stream,_direction,resolution,box);
}

bool ElementaryEdge::intresicEqual(const ElementaryEdge *other) const
{
  return _ptr==other->_ptr;
}

bool ElementaryEdge::intresicEqualDirSensitive(const ElementaryEdge *other) const
{
  return ( _direction==other->_direction ) && (_ptr==other->_ptr);
}

bool ElementaryEdge::intresincEqCoarse(const Edge *other) const
{
  return _ptr==other;
}

bool ElementaryEdge::isEqual(const ElementaryEdge& other) const
{
  return _ptr->isEqual(*other._ptr);
}

/*!
 * Called by QuadraticPolygon::splitAbs method.
 */
void ElementaryEdge::fillGlobalInfoAbs(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                                       std::vector<int>& edgesThis, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int> mapAddCoo) const
{
  _ptr->fillGlobalInfoAbs(_direction,mapThis,mapOther,offset1,offset2,fact,baryX,baryY,edgesThis,addCoo,mapAddCoo);
}

/*!
 * Called by QuadraticPolygon::splitAbs method. Close to ElementaryEdge::fillGlobalInfoAbs method expect that here edgesOther (that replace edgesThis) is here an in/out parameter that only contains nodes
 * unsorted because the "other" mesh is not subdivided yet.
 */
void ElementaryEdge::fillGlobalInfoAbs2(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                                        std::vector<int>& edgesOther, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int>& mapAddCoo) const
{
  _ptr->fillGlobalInfoAbs2(mapThis,mapOther,offset1,offset2,fact,baryX,baryY,edgesOther,addCoo,mapAddCoo);
}

/*!
 * This method builds from a \a start node, an \a end node and a direction a new ElementaryEdge.
 */
ElementaryEdge *ElementaryEdge::BuildEdgeFromStartEndDir(bool direction, INTERP_KERNEL::Node *start, INTERP_KERNEL::Node *end)
{
  Edge *ptr=Edge::BuildEdgeFrom(start,end);
  return new ElementaryEdge(ptr,direction);
}
