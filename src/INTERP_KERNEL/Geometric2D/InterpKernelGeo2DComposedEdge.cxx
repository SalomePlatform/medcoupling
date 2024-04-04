// Copyright (C) 2007-2024  CEA, EDF
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

#include "InterpKernelGeo2DComposedEdge.hxx"
#include "InterpKernelGeo2DEdge.hxx"
#include "InterpKernelGeo2DBounds.hxx"
#include "InterpKernelGeo2DElementaryEdge.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DEdgeInfLin.hxx"
#include "InterpKernelException.hxx"
#include "InterpKernelGeo2DNode.hxx"
#include "InterpKernelGeo2DPrecision.hxx"

#include <algorithm>
#include <list>
#include <cmath>
#include <map>
#include <iostream>
#include <math.h>
#include <functional>
#include <cstddef>
#include <memory>
#include <iterator>
#include <ostream>
#include <set>
#include <vector>

using namespace INTERP_KERNEL;

ComposedEdge::ComposedEdge(const ComposedEdge& other)
{
  for(auto _sub_edge : other._sub_edges)
    _sub_edges.push_back(_sub_edge->clone());
}

ComposedEdge::~ComposedEdge()
{
  clearAll(_sub_edges.begin());
}

void ComposedEdge::setValueAt(int i, Edge *e, bool direction)
{
  auto it=_sub_edges.begin();
  for(int j=0;j<i;j++)
    it++;
  delete *it;
  *it=new ElementaryEdge(e,direction);
}

/*! \cond HIDDEN_ITEMS */
struct AbsEdgeCmp
{
  AbsEdgeCmp(ElementaryEdge *b):_b1(b) { }
  bool operator()(ElementaryEdge *a) { return a->getPtr()==_b1->getPtr();}

  ElementaryEdge *_b1;
};
/*! \endcond */

double ComposedEdge::getCommonLengthWith(const ComposedEdge& other) const
{
  double ret=0.;
  for(auto _sub_edge : _sub_edges)
    {
      if(find_if(other._sub_edges.begin(),other._sub_edges.end(),AbsEdgeCmp(_sub_edge))!=other._sub_edges.end())
        {
          const auto *tmp=static_cast<const ElementaryEdge *>(_sub_edge);
          ret+=tmp->getCurveLength();
        }
    }
  return ret;
}

void ComposedEdge::clear()
{
  clearAll(_sub_edges.begin());
  _sub_edges.clear();
}

void ComposedEdge::pushBack(Edge *edge, bool direction)
{
  _sub_edges.push_back(new ElementaryEdge(edge,direction));
}

void ComposedEdge::pushBack(ElementaryEdge *elem)
{
  _sub_edges.push_back(elem);
}

void ComposedEdge::pushBack(ComposedEdge *elem)
{
  std::list<ElementaryEdge *> *elemsOfElem=elem->getListBehind();
  _sub_edges.insert(_sub_edges.end(),elemsOfElem->begin(),elemsOfElem->end());
}

ElementaryEdge *ComposedEdge::operator[](int i) const
{
  auto iter=_sub_edges.begin();
  for(int ii=0;ii<i;ii++)
    iter++;
  return *iter;
}

void ComposedEdge::reverse()
{
  _sub_edges.reverse();
  for(auto & _sub_edge : _sub_edges)
    _sub_edge->reverse();
}

bool ComposedEdge::presenceOfOn() const
{
  bool ret=false;
  for(auto iter=_sub_edges.begin();iter!=_sub_edges.end() && !ret;iter++)
    ret=((*iter)->getLoc()==FULL_ON_1);
  return ret;
}

bool ComposedEdge::presenceOfQuadraticEdge() const
{
  bool ret=false;
  for(auto iter=_sub_edges.begin();iter!=_sub_edges.end() && !ret;iter++)
    {
      Edge *e=(*iter)->getPtr();
      if(e)
        ret=dynamic_cast<EdgeArcCircle*>(e)!=nullptr;
    }
  return ret;
}

void ComposedEdge::initLocations() const
{
  for(auto _sub_edge : _sub_edges)
    _sub_edge->initLocations();
}

/**
 * Reset the status of all edges (OUT, IN, ON) because they were potentially assigned
 * by the previous candidate processing.
 */
void ComposedEdge::InitLocationsWithOther(const ComposedEdge& first, const ComposedEdge& other)
{
  std::set<Edge *> s1,s2;
  for(auto _sub_edge : first._sub_edges)
    s1.insert(_sub_edge->getPtr());
  for(auto _sub_edge : other._sub_edges)
    s2.insert(_sub_edge->getPtr());
  first.initLocations();
  other.initLocations();
  std::vector<Edge *> s3;
  std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::back_insert_iterator< std::vector<Edge *> >(s3));
  for(auto it3 : s3)
    it3->declareOn();
}

ComposedEdge *ComposedEdge::clone() const
{
  return new ComposedEdge(*this);
}

bool ComposedEdge::isNodeIn(Node *n) const
{
  bool ret=false;
  for(auto iter=_sub_edges.begin();iter!=_sub_edges.end() && !ret;iter++)
    ret=(*iter)->isNodeIn(n);
  return ret;
}

/*!
 * This method computes the area of 'this'.
 * By definition :
 * \f[
 * Area=\int_{Polygon} dS
 * \f]
 * Thanks to Green's theorem we have.
 * \f[
 * \int_{Polygon} x \cdot dS=\sum_{0 \leq i < nb of edges} -\int_{Edge_{i}}ydx=\sum_{0 \leq i < nb of edges} AreaOfZone_{Edge_{i}}
 * \f]
 * Where \f$ AreaOfZone_{i} \f$ is computed virtually by INTERP_KERNEL::Edge::getAreaOfZone with following formula :
 * \f[
 * AreaOfZone_{i}=\int_{Edge_{i}} -ydx
 * \f]
 */
double ComposedEdge::getArea() const
{
  double ret=0.;
  for(auto _sub_edge : _sub_edges)
    ret+=_sub_edge->getAreaOfZone();
  return ret;
}

double ComposedEdge::getPerimeter() const
{
  double ret=0.;
  for(auto _sub_edge : _sub_edges)
    ret+=_sub_edge->getCurveLength();
  return ret;
}

double ComposedEdge::getHydraulicDiameter() const
{
  return 4*fabs(getArea())/getPerimeter();
}

/*!
 * This method computes barycenter of 'this' by returning xG in bary[0] and yG in bary[1].
 * By definition :
 * \f[
 * Area \cdot x_{G}=\int_{Polygon} x \cdot dS
 * \f]
 * \f[
 * Area \cdot y_{G}=\int_{Polygon} y \cdot dS
 * \f]
 * Thanks to Green's theorem we have.
 * \f[
 * \int_{Polygon} x \cdot dS=\sum_{0 \leq i < nb of edges} -\int_{Edge_{i}}yxdx
 * \f]
 * \f[
 * \int_{Polygon} y \cdot dS=\sum_{0 \leq i < nb of edges} -\int_{Edge_{i}}\frac{y^{2}}{2}dx
 * \f]
 * Area is computed using the same principle than described in INTERP_KERNEL::ComposedEdge::getArea method.
 * \f$ -\int_{Edge_{i}}yxdx \f$ and \f$ -\int_{Edge_{i}}\frac{y^{2}}{2}dx \f$ are computed virtually with INTERP_KERNEL::Edge::getBarycenterOfZone.
 */
void ComposedEdge::getBarycenter(double *bary) const
{
  bary[0]=0.;
  bary[1]=0.;
  double area=0.;
  for(auto _sub_edge : _sub_edges)
    {
      _sub_edge->getBarycenterOfZone(bary);
      area+=_sub_edge->getAreaOfZone();
    }
  bary[0]/=area;
  bary[1]/=area;
}

/*!
 * Idem ComposedEdge::getBarycenter except that the special case where _sub_edges==1 is dealt here.
 */
void ComposedEdge::getBarycenterGeneral(double *bary) const
{
  if(_sub_edges.empty())
    throw INTERP_KERNEL::Exception("ComposedEdge::getBarycenterGeneral called on an empty polygon !");
  if(_sub_edges.size()>2)
    return getBarycenter(bary);
  double w;
  _sub_edges.back()->getBarycenter(bary,w);
}

double ComposedEdge::normalizeMe(double& xBary, double& yBary)
{
  Bounds b;
  b.prepareForAggregation();
  fillBounds(b);
  double const dimChar=b.getCaracteristicDim();
  b.getBarycenter(xBary,yBary);
  applyGlobalSimilarity(xBary,yBary,dimChar);
  return dimChar;
}

double ComposedEdge::normalize(ComposedEdge *other, double& xBary, double& yBary)
{
  Bounds b;
  b.prepareForAggregation();
  fillBounds(b); 
  other->fillBounds(b);
  double const dimChar=b.getCaracteristicDim();
  b.getBarycenter(xBary,yBary);
  applyGlobalSimilarity(xBary,yBary,dimChar);
  other->applyGlobalSimilarity(xBary,yBary,dimChar);
  return dimChar;
}

/*!
 * This method operates the opposite operation than ComposedEdge::applyGlobalSimilarity.
 */
void ComposedEdge::unApplyGlobalSimilarityExt(ComposedEdge& other, double xBary, double yBary, double fact)
{
  initNodeHitStatus();
  other.initNodeHitStatus();
  unApplySimilarityOnMyNodes(xBary,yBary,fact);
  other.unApplySimilarityOnMyNodesIfNotAlreadyHit(xBary,yBary,fact);
  initEdgeHitStatus();
  other.initEdgeHitStatus();
  unApplySimilarityOnMyEdges(xBary,yBary,fact);
  other.unApplySimilarityOnMyEdgesIfNotAlreadyHit(xBary,yBary,fact);
}

double ComposedEdge::normalizeExt(ComposedEdge *other, double& xBary, double& yBary)
{
  Bounds b;
  b.prepareForAggregation();
  fillBounds(b); 
  other->fillBounds(b);
  double const dimChar=b.getCaracteristicDim();
  b.getBarycenter(xBary,yBary);
  applyGlobalSimilarity2(other,xBary,yBary,dimChar);
  return dimChar;
}

void ComposedEdge::dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const
{
  stream.precision(10);
  for(auto _sub_edge : _sub_edges)
    _sub_edge->dumpInXfigFile(stream,resolution,box);
}

void ComposedEdge::dumpToCout(const std::map<INTERP_KERNEL::Node *,int>& mapp) const
{
  int i=0;
  for(auto iter=_sub_edges.begin();iter!=_sub_edges.end();iter++, i++)
    (*iter)->dumpToCout(mapp, i);
  std::cout << std::endl;
}

Node *ComposedEdge::getEndNode() const
{
  return _sub_edges.back()->getEndNode();
}

Node *ComposedEdge::getStartNode() const
{
  return _sub_edges.front()->getStartNode();
}

bool ComposedEdge::changeEndNodeWith(Node *node) const
{
  return _sub_edges.back()->changeEndNodeWith(node);
}

bool ComposedEdge::changeStartNodeWith(Node *node) const
{
  return _sub_edges.front()->changeStartNodeWith(node);
}

void ComposedEdge::fillBounds(Bounds& output) const
{
  for(auto _sub_edge : _sub_edges)
    _sub_edge->fillBounds(output);
}

/*!
 * \b WARNING : applies similarity \b ONLY on edges without any change on Nodes. To perform a global similarity call applyGlobalSimilarity.
 */
void ComposedEdge::applySimilarity(double xBary, double yBary, double dimChar)
{
  for(auto & _sub_edge : _sub_edges)
    _sub_edge->applySimilarity(xBary,yBary,dimChar);
}

/*!
 * Perform Similarity transformation on all elements of this Nodes and Edges.
 */
void ComposedEdge::applyGlobalSimilarity(double xBary, double yBary, double dimChar)
{
  std::set<Node *> allNodes;
  getAllNodes(allNodes);
  for(auto allNode : allNodes)
    allNode->applySimilarity(xBary,yBary,dimChar);
  for(auto & _sub_edge : _sub_edges)
    _sub_edge->applySimilarity(xBary,yBary,dimChar);
}

/*!
 * Perform Similarity transformation on all elements of this Nodes and Edges on 'this' and 'other'.
 * Nodes can be shared between 'this' and 'other'.
 */
void ComposedEdge::applyGlobalSimilarity2(ComposedEdge *other, double xBary, double yBary, double dimChar)
{
  initNodeHitStatus();
  other->initNodeHitStatus();
  applySimilarityOnMyNodes(xBary,yBary,dimChar);
  other->applySimilarityOnMyNodesIfNotAlreadyHit(xBary,yBary,dimChar);
  initEdgeHitStatus();
  other->initEdgeHitStatus();
  applySimilarityOnMyEdges(xBary,yBary,dimChar);
  other->applySimilarityOnMyEdgesIfNotAlreadyHit(xBary,yBary,dimChar);
}

/*!
 * This method append to param 'partConsidered' the part of length of subedges IN or ON.
 * @param partConsidered INOUT param.
 */
void ComposedEdge::dispatchPerimeter(double& partConsidered) const
{
  for(auto _sub_edge : _sub_edges)
    {
      TypeOfEdgeLocInPolygon const loc=_sub_edge->getLoc();
      if(loc==FULL_IN_1 || loc==FULL_ON_1)
        partConsidered+=_sub_edge->getCurveLength();
    }
}

/*!
 * Idem dispatchPerimeterExcl except that when a subedge is declared as ON this subedge is counted in commonPart.
 */
void ComposedEdge::dispatchPerimeterExcl(double& partConsidered, double& commonPart) const
{
  for(auto _sub_edge : _sub_edges)
    {
      TypeOfEdgeLocInPolygon const loc=_sub_edge->getLoc();
      if(loc==FULL_IN_1)
        partConsidered+=_sub_edge->getCurveLength();
      if(loc==FULL_ON_1)
        commonPart+=_sub_edge->getCurveLength();
    }
}

void ComposedEdge::getAllNodes(std::set<Node *>& output) const
{
  auto iter=_sub_edges.begin();
  for(;iter!=_sub_edges.end();iter++)
    (*iter)->getAllNodes(output);
}

void ComposedEdge::initNodeHitStatus() const
{
  for(auto _sub_edge : _sub_edges)
    {
      _sub_edge->getStartNode()->initHitStatus();
      _sub_edge->getEndNode()->initHitStatus();
    }
}

void ComposedEdge::applySimilarityOnMyNodes(double xBary, double yBary, double dimChar) const
{
  for(auto _sub_edge : _sub_edges)
    {
      _sub_edge->getStartNode()->hitMeAlone(xBary,yBary,dimChar);
      _sub_edge->getEndNode()->hitMeAlone(xBary,yBary,dimChar);
    }
}

void ComposedEdge::unApplySimilarityOnMyNodes(double xBary, double yBary, double dimChar) const
{
  for(auto _sub_edge : _sub_edges)
    {
      _sub_edge->getStartNode()->unHitMeAlone(xBary,yBary,dimChar);
      _sub_edge->getEndNode()->unHitMeAlone(xBary,yBary,dimChar);
    }
}

void ComposedEdge::applySimilarityOnMyNodesIfNotAlreadyHit(double xBary, double yBary, double dimChar) const
{
  for(auto _sub_edge : _sub_edges)
    {
      _sub_edge->getStartNode()->hitMeAfter(xBary,yBary,dimChar);
      _sub_edge->getEndNode()->hitMeAfter(xBary,yBary,dimChar);
    }
}

void ComposedEdge::unApplySimilarityOnMyNodesIfNotAlreadyHit(double xBary, double yBary, double dimChar) const
{
  for(auto _sub_edge : _sub_edges)
    {
      _sub_edge->getStartNode()->unHitMeAfter(xBary,yBary,dimChar);
      _sub_edge->getEndNode()->unHitMeAfter(xBary,yBary,dimChar);
    }
}

void ComposedEdge::initEdgeHitStatus() const
{
  for(auto _sub_edge : _sub_edges)
    _sub_edge->getPtr()->initHitStatus();
}

void ComposedEdge::applySimilarityOnMyEdges(double xBary, double yBary, double dimChar) const
{
  for(auto _sub_edge : _sub_edges)
    _sub_edge->getPtr()->hitMeAlone(xBary,yBary,dimChar);
}

void ComposedEdge::unApplySimilarityOnMyEdges(double xBary, double yBary, double dimChar) const
{
  for(auto _sub_edge : _sub_edges)
    _sub_edge->getPtr()->unHitMeAlone(xBary,yBary,dimChar);
}

void ComposedEdge::applySimilarityOnMyEdgesIfNotAlreadyHit(double xBary, double yBary, double dimChar) const
{
  for(auto _sub_edge : _sub_edges)
    _sub_edge->getPtr()->hitMeAfter(xBary,yBary,dimChar);
}

void ComposedEdge::unApplySimilarityOnMyEdgesIfNotAlreadyHit(double xBary, double yBary, double dimChar) const
{
  for(auto _sub_edge : _sub_edges)
    _sub_edge->getPtr()->unHitMeAfter(xBary,yBary,dimChar);
}

void ComposedEdge::getBarycenter(double *bary, double& weigh) const
{
  weigh=0.; bary[0]=0.; bary[1]=0.;
  double tmp1,tmp2[2];
  for(auto _sub_edge : _sub_edges)
    {
      _sub_edge->getBarycenter(tmp2,tmp1);
      weigh+=tmp1;
      bary[0]+=tmp1*tmp2[0];
      bary[1]+=tmp1*tmp2[1];
    }
  bary[0]/=weigh;
  bary[1]/=weigh;
}

/*!
 * This method makes the hypothesis that \a nodeToTest can be either IN or OUT.
 * 
 * \sa ComposedEdge::isInOrOut2
 */
bool ComposedEdge::isInOrOut(Node *nodeToTest) const
{
  Bounds b; b.prepareForAggregation();
  fillBounds(b);
  if(b.nearlyWhere((*nodeToTest)[0],(*nodeToTest)[1])==OUT)
    return false;
  std::set< IntersectElement > inOutSwitch;
  std::set<Node *> nodes;
  getAllNodes(nodes);
  double const ref(isInOrOutAlg(nodeToTest,nodes,inOutSwitch));
  bool ret(false);
  for(const auto & iter4 : inOutSwitch)
    {
      if(iter4.getVal1()<ref)
        {
          if(iter4.getNodeOnly()->getLoc()==ON_1)
            ret=!ret;
        }
      else
        break;
    }
  return ret;
}

/*!
 * This method is close to ComposedEdge::isInOrOut behaviour except that here EPSILON is taken into account to detect if it is IN or OUT.
 * If \a nodeToTest is close to an edge in \a this, true will be returned even if it is outside informatically from \a this.
 *
 * \sa ComposedEdge::isInOrOut
 */
bool ComposedEdge::isInOrOut2(Node *nodeToTest) const
{
  std::set< IntersectElement > inOutSwitch;
  std::set<Node *> nodes;
  getAllNodes(nodes);
  for(auto node : nodes)
    if(sqrt(node->distanceWithSq(*nodeToTest))<QuadraticPlanarPrecision::getPrecision())
      return true;
  double const ref(isInOrOutAlg(nodeToTest,nodes,inOutSwitch));
  bool ret(false);
  for(const auto & iter4 : inOutSwitch)
    {
      double const val(iter4.getVal1());
      if(fabs(val-ref)>=QuadraticPlanarPrecision::getPrecision())
        {
          if(val<ref)
            {
              if(iter4.getNodeOnly()->getLoc()==ON_1)
                ret=!ret;
            }
          else
            break;
        }
      else
        return true;
    }
  return ret;
}

double ComposedEdge::isInOrOutAlg(Node *nodeToTest, const std::set<Node*>& nodes, std::set< IntersectElement >& inOutSwitch) const
{
  // searching for e1
  std::set<double> radialDistributionOfNodes;
  std::set<Node *>::const_iterator iter;
  for(iter=nodes.begin();iter!=nodes.end();iter++)
    radialDistributionOfNodes.insert(nodeToTest->getSlope(*(*iter)));
  std::vector<double> radialDistrib(radialDistributionOfNodes.begin(),radialDistributionOfNodes.end());
  radialDistributionOfNodes.clear();
  std::vector<double> radialDistrib2(radialDistrib.size());
  copy(radialDistrib.begin()+1,radialDistrib.end(),radialDistrib2.begin());
  radialDistrib2.back()=M_PI+radialDistrib.front();
  std::vector<double> radialDistrib3(radialDistrib.size());
  std::transform(radialDistrib2.begin(),radialDistrib2.end(),radialDistrib.begin(),radialDistrib3.begin(),std::minus<double>());
  auto const iter3=max_element(radialDistrib3.begin(),radialDistrib3.end());
  std::size_t const i=iter3-radialDistrib3.begin();
  // ok for e1 - Let's go.
  auto *e1=new EdgeInfLin(nodeToTest,radialDistrib[i]+radialDistrib3[i]/2.);
  double const ref=e1->getCharactValue(*nodeToTest);
  for(auto val : _sub_edges)
    {
      if(val)
        {
          Edge *e=val->getPtr();
          std::unique_ptr<EdgeIntersector> intersc(Edge::BuildIntersectorWith(e1,e));
          bool obviousNoIntersection,areOverlapped;
          intersc->areOverlappedOrOnlyColinears(obviousNoIntersection,areOverlapped);
          if(obviousNoIntersection)
            {
              continue;
            }
          if(!areOverlapped)
            {
              std::list< IntersectElement > const listOfIntesc=intersc->getIntersectionsCharacteristicVal();
              for(auto & iter2 : listOfIntesc)
                if(iter2.isIncludedByBoth())
                  inOutSwitch.insert(iter2);
            }
          //if overlapped we can forget
        }
      else
        throw Exception("Invalid use of ComposedEdge::isInOrOutAlg : only one level supported !");
    }
  e1->decrRef();
  return ref;
}

bool ComposedEdge::getDirection() const
{
  throw Exception("ComposedEdge::getDirection : no sense");
}

bool ComposedEdge::intresincEqCoarse(const Edge *other) const
{
  if(_sub_edges.size()!=1)
    return false;
  return _sub_edges.front()->intresincEqCoarse(other);
}

void ComposedEdge::clearAll(std::list<ElementaryEdge *>::iterator startToDel)
{
  for(auto iter=startToDel;iter!=_sub_edges.end();iter++)
    delete (*iter);
}
