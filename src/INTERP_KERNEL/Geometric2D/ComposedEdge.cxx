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
#include "ComposedEdge.hxx"
#include "ElementaryEdge.hxx"
#include "EdgeInfLin.hxx"
#include "InterpKernelException.hxx"

#include <algorithm>
#include <iterator>
#include <set>

using namespace std;
using namespace INTERP_KERNEL;

ComposedEdge::ComposedEdge(const ComposedEdge& other)
{
  for(list<ElementaryEdge *>::const_iterator iter=other._sub_edges.begin();iter!=other._sub_edges.end();iter++)
    _sub_edges.push_back((*iter)->clone());
}

ComposedEdge::~ComposedEdge()
{
  clearAll(_sub_edges.begin());
}

void ComposedEdge::setValueAt(int i, Edge *e, bool direction)
{
  list<ElementaryEdge*>::iterator it=_sub_edges.begin();
  for(int j=0;j<i;j++)
    it++;
  delete *it;
  *it=new ElementaryEdge(e,direction);
}

struct AbsEdgeCmp
{
  AbsEdgeCmp(ElementaryEdge *b):_b1(b) { }
  bool operator()(ElementaryEdge *a) { return a->getPtr()==_b1->getPtr();}

  ElementaryEdge *_b1;
};

double ComposedEdge::getCommonLengthWith(const ComposedEdge& other) const
{
  double ret=0.;
  for(list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
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
  list<ElementaryEdge *> *elemsOfElem=elem->getListBehind();
  _sub_edges.insert(_sub_edges.end(),elemsOfElem->begin(),elemsOfElem->end());
}

ElementaryEdge *ComposedEdge::operator[](int i) const
{
  list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();
  for(int ii=0;ii<i;ii++)
    iter++;
  return *iter;
}

void ComposedEdge::reverse()
{
  _sub_edges.reverse();
  for(list<ElementaryEdge *>::iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->reverse();
}

void ComposedEdge::initLocations() const
{
  for(list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->initLocations();
}

ComposedEdge *ComposedEdge::clone() const
{
  return new ComposedEdge(*this);
}

bool ComposedEdge::isNodeIn(Node *n) const
{
  bool ret=false;
  for(list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end() && !ret;iter++)
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
  for(list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    ret+=(*iter)->getAreaOfZone();
  return ret;
}

double ComposedEdge::getPerimeter() const
{
  double ret=0.;
  for(list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
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
  for(list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      (*iter)->getBarycenterOfZone(bary);
      area+=(*iter)->getAreaOfZone();
    }
  bary[0]/=area;
  bary[1]/=area;
}

double ComposedEdge::normalize(ComposedEdge *other)
{
  Bounds b;
  b.prepareForAggregation();
  fillBounds(b); 
  other->fillBounds(b);
  double dimChar=b.getCaracteristicDim();
  double xBary,yBary;
  b.getBarycenter(xBary,yBary);
  applyGlobalSimilarity(xBary,yBary,dimChar);
  other->applyGlobalSimilarity(xBary,yBary,dimChar);
  return dimChar;
}

void ComposedEdge::dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const
{
  stream.precision(10);
  for(list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
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
  for(list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->fillBounds(output);
}

/*!
 * \b WARNING : applies similarity \b ONLY on edges without any change on Nodes. To perform a global similarity call applyGlobalSimilarity.
 */
void ComposedEdge::applySimilarity(double xBary, double yBary, double dimChar)
{
  for(list<ElementaryEdge *>::iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->applySimilarity(xBary,yBary,dimChar);
}

/*!
 * Perform Similarity transformation on all elements of this Nodes and Edges.
 */
void ComposedEdge::applyGlobalSimilarity(double xBary, double yBary, double dimChar)
{
  set<Node *> allNodes;
  getAllNodes(allNodes);
  for(set<Node *>::iterator iter=allNodes.begin();iter!=allNodes.end();iter++)
    (*iter)->applySimilarity(xBary,yBary,dimChar);
  for(list<ElementaryEdge *>::iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    (*iter)->applySimilarity(xBary,yBary,dimChar);
}

/*!
 * This method append to param 'partConsidered' the part of length of subedges IN or ON.
 * @param partConsidered INOUT param.
 */
void ComposedEdge::dispatchPerimeter(double& partConsidered) const
{
  for(list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
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
  for(list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
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
  list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();
  for(;iter!=_sub_edges.end();iter++)
    (*iter)->getAllNodes(output);
}

void ComposedEdge::getBarycenter(double *bary, double& weigh) const
{
  weigh=0.; bary[0]=0.; bary[1]=0.;
  double tmp1,tmp2[2];
  for(list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      (*iter)->getBarycenter(tmp2,tmp1);
      weigh+=tmp1;
      bary[0]+=tmp1*tmp2[0];
      bary[1]+=tmp1*tmp2[1];
    }
  bary[0]/=weigh;
  bary[1]/=weigh;
}

bool ComposedEdge::isInOrOut(Node *nodeToTest) const
{
  Bounds b; b.prepareForAggregation();
  fillBounds(b);
  if(b.nearlyWhere((*nodeToTest)[0],(*nodeToTest)[1])==OUT)
    return false;
  // searching for e1
  set<Node *> nodes;
  getAllNodes(nodes);
  set<double> radialDistributionOfNodes;
  set<Node *>::const_iterator iter;
  for(iter=nodes.begin();iter!=nodes.end();iter++)
    radialDistributionOfNodes.insert(nodeToTest->getSlope(*(*iter)));
  vector<double> radialDistrib(radialDistributionOfNodes.begin(),radialDistributionOfNodes.end());
  radialDistributionOfNodes.clear();
  vector<double> radialDistrib2(radialDistrib.size());
  copy(radialDistrib.begin()+1,radialDistrib.end(),radialDistrib2.begin());
  radialDistrib2.back()=M_PI+radialDistrib.front();
  vector<double> radialDistrib3(radialDistrib.size());
  transform(radialDistrib2.begin(),radialDistrib2.end(),radialDistrib.begin(),radialDistrib3.begin(),minus<double>());
  vector<double>::iterator iter3=max_element(radialDistrib3.begin(),radialDistrib3.end());
  int i=iter3-radialDistrib3.begin();
  // ok for e1 - Let's go.
  EdgeInfLin *e1=new EdgeInfLin(nodeToTest,radialDistrib[i]+radialDistrib3[i]/2.);
  double ref=e1->getCharactValue(*nodeToTest);
  set< IntersectElement > inOutSwitch;
  for(list<ElementaryEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      ElementaryEdge *val=(*iter);
      if(val)
        {
          Edge *e=val->getPtr();
          auto_ptr<EdgeIntersector> intersc(Edge::buildIntersectorWith(e1,e));
          bool obviousNoIntersection,areOverlapped;
          intersc->areOverlappedOrOnlyColinears(0,obviousNoIntersection,areOverlapped);
          if(obviousNoIntersection)
            {
              continue;
            }
          if(!areOverlapped)
            {
              list< IntersectElement > listOfIntesc=intersc->getIntersectionsCharacteristicVal();
              for(list< IntersectElement >::iterator iter2=listOfIntesc.begin();iter2!=listOfIntesc.end();iter2++)
                if((*iter2).isIncludedByBoth())
                  inOutSwitch.insert(*iter2);
              }
          //if overlapped we can forget
        }
      else
        throw Exception("Invalid use of ComposedEdge::isInOrOut : only one level supported !");
    }
  e1->decrRef();
  bool ret=false;
  for(set< IntersectElement >::iterator iter=inOutSwitch.begin();iter!=inOutSwitch.end();iter++)
    {
      if((*iter).getVal1()<ref)
        {
          if((*iter).getNodeOnly()->getLoc()==ON_1)
            ret=!ret;
        }
      else
        break;
    }
  return ret;
}

/*bool ComposedEdge::isInOrOut(Node *aNodeOn, Node *nodeToTest) const
{
  
  EdgeInfLin *e1=new EdgeInfLin(aNodeOn,nodeToTest);
  double ref=e1->getCharactValue(*nodeToTest);
  set< IntersectElement > inOutSwitch;
  for(vector<AbstractEdge *>::const_iterator iter=_sub_edges.begin();iter!=_sub_edges.end();iter++)
    {
      ElementaryEdge *val=dynamic_cast<ElementaryEdge *>(*iter);
      if(val)
        {
          Edge *e=val->getPtr();
          auto_ptr<Intersector> intersc(Edge::buildIntersectorWith(e1,e));
          bool obviousNoIntersection,areOverlapped;
          intersc->areOverlappedOrOnlyColinears(0,obviousNoIntersection,areOverlapped);
          if(obviousNoIntersection)
            {
              continue;
            }
          if(!areOverlapped)
            {
              list< IntersectElement > listOfIntesc=intersc->getIntersectionsCharacteristicVal();
              for(list< IntersectElement >::iterator iter2=listOfIntesc.begin();iter2!=listOfIntesc.end();iter2++)
                if((*iter2).isIncludedByBoth())
                  inOutSwitch.insert(*iter2);
              }
          //if overlapped we can forget
        }
      else
        throw Exception("Invalid use of ComposedEdge::isInOrOut : only one level supported !");
    }
  e1->decrRef();
  bool ret=false;
  for(set< IntersectElement >::iterator iter=inOutSwitch.begin();iter!=inOutSwitch.end();iter++)
    {
      if((*iter).getVal1()<ref)
        {
          if((*iter).getNodeOnly()->getLoc()==ON_1)
            ret=!ret;
        }
      else
        break;
    }
  return ret;
}*/

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

void ComposedEdge::clearAll(list<ElementaryEdge *>::iterator startToDel)
{
  for(list<ElementaryEdge *>::iterator iter=startToDel;iter!=_sub_edges.end();iter++)
    delete (*iter);
}
