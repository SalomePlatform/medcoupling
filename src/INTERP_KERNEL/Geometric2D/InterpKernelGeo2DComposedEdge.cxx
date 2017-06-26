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

#include "InterpKernelGeo2DComposedEdge.hxx"
#include "InterpKernelGeo2DElementaryEdge.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DEdgeInfLin.hxx"
#include "InterpKernelException.hxx"

#include <algorithm>
#include <memory>
#include <iterator>
#include <set>

using namespace INTERP_KERNEL;

ComposedEdge::ComposedEdge(const ComposedEdge& other)
{
  for(std::list<ElementaryEdge *>::const_iterator iter=other._sub_edges.begin();iter!=other._sub_edges.end();iter++)
    _sub_edges.push_back((*iter)->clone());
}

ComposedEdge::~ComposedEdge()
{
  clearAll(_sub_edges.begin());
}

void ComposedEdge::setValueAt(int i, Edge *e, bool direction)
{
  std::list<ElementaryEdge*>::iterator it=_sub_edges.begin();
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
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      if(find_if(other._sub_edges.begin(),other._sub_edges.end(),AbsEdgeCmp(*iter))!=other._sub_edges.end())
        {
          const ElementaryEdge *tmp=static_cast<const ElementaryEdge *>(*iter);
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
  std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();
  for(int ii=0;ii<i;ii++)
    iter++;
  return *iter;
}

void ComposedEdge::reverse()
{
  _sub_edges.reverse();
  for(std::list<ElementaryEdge *>::iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->reverse();
}

bool ComposedEdge::presenceOfOn() const
{
  bool ret=false;
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end() && !ret;iter++)
    ret=((*iter)->getLoc()==FULL_ON_1);
  return ret;
}

bool ComposedEdge::presenceOfQuadraticEdge() const
{
  bool ret=false;
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end() && !ret;iter++)
    {
      Edge *e=(*iter)->getPtr();
      if(e)
        ret=dynamic_cast<EdgeArcCircle*>(e)!=0;
    }
  return ret;
}

void ComposedEdge::initLocations() const
{
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->initLocations();
}

/**
 * Reset the status of all edges (OUT, IN, ON) because they were potentially assigned
 * by the previous candidate processing.
 */
void ComposedEdge::InitLocationsWithOther(const ComposedEdge& first, const ComposedEdge& other)
{
  std::set<Edge *> s1,s2;
  for(std::list<ElementaryEdge *>::const_iterator it1=first._sub_edges.begin();it1!=first._sub_edges.end();it1++)
    s1.insert((*it1)->getPtr());
  for(std::list<ElementaryEdge *>::const_iterator it2=other._sub_edges.begin();it2!=other._sub_edges.end();it2++)
    s2.insert((*it2)->getPtr());
  first.initLocations();
  other.initLocations();
  std::vector<Edge *> s3;
  std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::back_insert_iterator< std::vector<Edge *> >(s3));
  for(std::vector<Edge *>::const_iterator it3=s3.begin();it3!=s3.end();it3++)
    (*it3)->declareOn();
}

ComposedEdge *ComposedEdge::clone() const
{
  return new ComposedEdge(*this);
}

bool ComposedEdge::isNodeIn(Node *n) const
{
  bool ret=false;
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end() && !ret;iter++)
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
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    ret+=(*iter)->getAreaOfZone();
  return ret;
}

double ComposedEdge::getPerimeter() const
{
  double ret=0.;
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    ret+=(*iter)->getCurveLength();
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
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      (*iter)->getBarycenterOfZone(bary);
      area+=(*iter)->getAreaOfZone();
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
  double dimChar=b.getCaracteristicDim();
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
  double dimChar=b.getCaracteristicDim();
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
  double dimChar=b.getCaracteristicDim();
  b.getBarycenter(xBary,yBary);
  applyGlobalSimilarity2(other,xBary,yBary,dimChar);
  return dimChar;
}

void ComposedEdge::dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const
{
  stream.precision(10);
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->dumpInXfigFile(stream,resolution,box);
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
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->fillBounds(output);
}

/*!
 * \b WARNING : applies similarity \b ONLY on edges without any change on Nodes. To perform a global similarity call applyGlobalSimilarity.
 */
void ComposedEdge::applySimilarity(double xBary, double yBary, double dimChar)
{
  for(std::list<ElementaryEdge *>::iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->applySimilarity(xBary,yBary,dimChar);
}

/*!
 * Perform Similarity transformation on all elements of this Nodes and Edges.
 */
void ComposedEdge::applyGlobalSimilarity(double xBary, double yBary, double dimChar)
{
  std::set<Node *> allNodes;
  getAllNodes(allNodes);
  for(std::set<Node *>::iterator iter=allNodes.begin();iter!=allNodes.end();iter++)
    (*iter)->applySimilarity(xBary,yBary,dimChar);
  for(std::list<ElementaryEdge *>::iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->applySimilarity(xBary,yBary,dimChar);
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
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      TypeOfEdgeLocInPolygon loc=(*iter)->getLoc();
      if(loc==FULL_IN_1 || loc==FULL_ON_1)
        partConsidered+=(*iter)->getCurveLength();
    }
}

/*!
 * Idem dispatchPerimeterExcl except that when a subedge is declared as ON this subedge is counted in commonPart.
 */
void ComposedEdge::dispatchPerimeterExcl(double& partConsidered, double& commonPart) const
{
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      TypeOfEdgeLocInPolygon loc=(*iter)->getLoc();
      if(loc==FULL_IN_1)
        partConsidered+=(*iter)->getCurveLength();
      if(loc==FULL_ON_1)
        commonPart+=(*iter)->getCurveLength();
    }
}

void ComposedEdge::getAllNodes(std::set<Node *>& output) const
{
  std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();
  for(;iter!=_sub_edges.end();iter++)
    (*iter)->getAllNodes(output);
}

void ComposedEdge::initNodeHitStatus() const
{
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      (*iter)->getStartNode()->initHitStatus();
      (*iter)->getEndNode()->initHitStatus();
    }
}

void ComposedEdge::applySimilarityOnMyNodes(double xBary, double yBary, double dimChar) const
{
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      (*iter)->getStartNode()->hitMeAlone(xBary,yBary,dimChar);
      (*iter)->getEndNode()->hitMeAlone(xBary,yBary,dimChar);
    }
}

void ComposedEdge::unApplySimilarityOnMyNodes(double xBary, double yBary, double dimChar) const
{
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      (*iter)->getStartNode()->unHitMeAlone(xBary,yBary,dimChar);
      (*iter)->getEndNode()->unHitMeAlone(xBary,yBary,dimChar);
    }
}

void ComposedEdge::applySimilarityOnMyNodesIfNotAlreadyHit(double xBary, double yBary, double dimChar) const
{
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      (*iter)->getStartNode()->hitMeAfter(xBary,yBary,dimChar);
      (*iter)->getEndNode()->hitMeAfter(xBary,yBary,dimChar);
    }
}

void ComposedEdge::unApplySimilarityOnMyNodesIfNotAlreadyHit(double xBary, double yBary, double dimChar) const
{
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      (*iter)->getStartNode()->unHitMeAfter(xBary,yBary,dimChar);
      (*iter)->getEndNode()->unHitMeAfter(xBary,yBary,dimChar);
    }
}

void ComposedEdge::initEdgeHitStatus() const
{
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->getPtr()->initHitStatus();
}

void ComposedEdge::applySimilarityOnMyEdges(double xBary, double yBary, double dimChar) const
{
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->getPtr()->hitMeAlone(xBary,yBary,dimChar);
}

void ComposedEdge::unApplySimilarityOnMyEdges(double xBary, double yBary, double dimChar) const
{
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->getPtr()->unHitMeAlone(xBary,yBary,dimChar);
}

void ComposedEdge::applySimilarityOnMyEdgesIfNotAlreadyHit(double xBary, double yBary, double dimChar) const
{
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->getPtr()->hitMeAfter(xBary,yBary,dimChar);
}

void ComposedEdge::unApplySimilarityOnMyEdgesIfNotAlreadyHit(double xBary, double yBary, double dimChar) const
{
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->getPtr()->unHitMeAfter(xBary,yBary,dimChar);
}

void ComposedEdge::getBarycenter(double *bary, double& weigh) const
{
  weigh=0.; bary[0]=0.; bary[1]=0.;
  double tmp1,tmp2[2];
  for(std::list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      (*iter)->getBarycenter(tmp2,tmp1);
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
  double ref(isInOrOutAlg(nodeToTest,nodes,inOutSwitch));
  bool ret(false);
  for(std::set< IntersectElement >::iterator iter4=inOutSwitch.begin();iter4!=inOutSwitch.end();iter4++)
    {
      if((*iter4).getVal1()<ref)
        {
          if((*iter4).getNodeOnly()->getLoc()==ON_1)
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
  for(std::set<Node *>::const_iterator iter=nodes.begin();iter!=nodes.end();iter++)
    if(sqrt((*iter)->distanceWithSq(*nodeToTest))<QuadraticPlanarPrecision::getPrecision())
      return true;
  double ref(isInOrOutAlg(nodeToTest,nodes,inOutSwitch));
  bool ret(false);
  for(std::set< IntersectElement >::iterator iter4=inOutSwitch.begin();iter4!=inOutSwitch.end();iter4++)
    {
      double val((*iter4).getVal1());
      if(fabs(val-ref)>=QuadraticPlanarPrecision::getPrecision())
        {
          if(val<ref)
            {
              if((*iter4).getNodeOnly()->getLoc()==ON_1)
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
  std::vector<double>::iterator iter3=max_element(radialDistrib3.begin(),radialDistrib3.end());
  int i=iter3-radialDistrib3.begin();
  // ok for e1 - Let's go.
  EdgeInfLin *e1=new EdgeInfLin(nodeToTest,radialDistrib[i]+radialDistrib3[i]/2.);
  double ref=e1->getCharactValue(*nodeToTest);
  for(std::list<ElementaryEdge *>::const_iterator iter4=_sub_edges.begin();iter4!=_sub_edges.end();iter4++)
    {
      ElementaryEdge *val=(*iter4);
      if(val)
        {
          Edge *e=val->getPtr();
          std::auto_ptr<EdgeIntersector> intersc(Edge::BuildIntersectorWith(e1,e));
          bool obviousNoIntersection,areOverlapped;
          intersc->areOverlappedOrOnlyColinears(0,obviousNoIntersection,areOverlapped);  // first parameter never used
          if(obviousNoIntersection)
            {
              continue;
            }
          if(!areOverlapped)
            {
              std::list< IntersectElement > listOfIntesc=intersc->getIntersectionsCharacteristicVal();
              for(std::list< IntersectElement >::iterator iter2=listOfIntesc.begin();iter2!=listOfIntesc.end();iter2++)
                if((*iter2).isIncludedByBoth())
                  inOutSwitch.insert(*iter2);
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
  for(std::list<ElementaryEdge *>::iterator iter=startToDel;iter!=_sub_edges.end();iter++)
    delete (*iter);
}
