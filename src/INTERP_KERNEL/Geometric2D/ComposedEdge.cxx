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
#include "InterpolationUtils.hxx"

#include <algorithm>
#include <iterator>
#include <set>

using namespace std;
using namespace INTERP_KERNEL;

ComposedEdge::ComposedEdge(const ComposedEdge& other)
{
  for(list<ElementaryEdge *>::const_iterator iter=other._subEdges.begin();iter!=other._subEdges.end();iter++)
    _subEdges.push_back((*iter)->clone());
}

ComposedEdge::~ComposedEdge()
{
  clearAll(_subEdges.begin());
}

void ComposedEdge::setValueAt(int i, Edge *e, bool direction)
{
  list<ElementaryEdge*>::iterator it=_subEdges.begin();
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
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    {
      if(find_if(other._subEdges.begin(),other._subEdges.end(),AbsEdgeCmp(*iter))!=other._subEdges.end())
        {
          const ElementaryEdge *tmp=static_cast<const ElementaryEdge *>(*iter);
          ret+=tmp->getCurveLength();
        }
    }
  return ret;
}

void ComposedEdge::clear()
{
  clearAll(_subEdges.begin());
  _subEdges.clear();
}

void ComposedEdge::pushBack(Edge *edge, bool direction)
{
  _subEdges.push_back(new ElementaryEdge(edge,direction));
}

void ComposedEdge::pushBack(ElementaryEdge *elem)
{
  _subEdges.push_back(elem);
}

void ComposedEdge::pushBack(ComposedEdge *elem)
{
  list<ElementaryEdge *> *elemsOfElem=elem->getListBehind();
  _subEdges.insert(_subEdges.end(),elemsOfElem->begin(),elemsOfElem->end());
}

ElementaryEdge *ComposedEdge::operator[](int i) const
{
  list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();
  for(int ii=0;ii<i;ii++)
    iter++;
  return *iter;
}

void ComposedEdge::reverse()
{
  _subEdges.reverse();
  for(list<ElementaryEdge *>::iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    (*iter)->reverse();
}

void ComposedEdge::initLocations() const
{
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    (*iter)->initLocations();
}

ComposedEdge *ComposedEdge::clone() const
{
  return new ComposedEdge(*this);
}

bool ComposedEdge::isNodeIn(Node *n) const
{
  bool ret=false;
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end() && !ret;iter++)
    ret=(*iter)->isNodeIn(n);
  return ret;
}

double ComposedEdge::getArea() const
{
  double ret=0.;
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    ret+=(*iter)->getAreaOfZone();
  return ret;
}

double ComposedEdge::getPerimeter() const
{
  double ret=0.;
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    ret+=(*iter)->getCurveLength();
  return ret;
}

double ComposedEdge::getHydraulicDiameter() const
{
  return 4*fabs(getArea())/getPerimeter();
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
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    (*iter)->dumpInXfigFile(stream,resolution,box);
}

Node *ComposedEdge::getEndNode() const
{
  return _subEdges.back()->getEndNode();
}

Node *ComposedEdge::getStartNode() const
{
  return _subEdges.front()->getStartNode();
}

bool ComposedEdge::changeEndNodeWith(Node *node) const
{
  return _subEdges.back()->changeEndNodeWith(node);
}

bool ComposedEdge::changeStartNodeWith(Node *node) const
{
  return _subEdges.front()->changeStartNodeWith(node);
}

void ComposedEdge::fillBounds(Bounds& output) const
{
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    (*iter)->fillBounds(output);
}

/*!
 * \b WARNING : applies similarity \b ONLY on edges without any change on Nodes. To perform a global similarity call applyGlobalSimilarity.
 */
void ComposedEdge::applySimilarity(double xBary, double yBary, double dimChar)
{
  for(list<ElementaryEdge *>::iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    (*iter)->applySimilarity(xBary,yBary,dimChar);
}

int ComposedEdge::getNbOfEdgeSonsOfSameId(int id) const
{
  int ret=0;
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    if((*iter)->getPtr()->getId()==id)
      ret++;
  return ret;
}

/*!
 * 'this' is supposed to be the father 'split'.
 * This method returns un vector 'nbOfCreatedNodes' the number of newly created nodes of each edge of this in the same order than the sub edges constituting 'this.'
 */
void ComposedEdge::dispatchForNode(const ComposedEdge& split, std::vector<int>& nbOfCreatedNodes) const
{
  int eRk=0;
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++,eRk++)
    nbOfCreatedNodes[eRk]=split.getNbOfEdgeSonsOfSameId((*iter)->getPtr()->getId())-1;
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
  for(list<ElementaryEdge *>::iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    (*iter)->applySimilarity(xBary,yBary,dimChar);
}

/*!
 * @param part1 INOUT parameter. If an edge in 'this' is found with an id in 'ids1', 'part1' is \b incremented.
 * @param part2 INOUT parameter. If an edge in 'this' is found with an id in 'ids2', 'part2' is \b incremented.
 * @param commonPart INOUT parameter. If an edge in 'this' is found with an id neither in 'ids1' nor 'ids2', 'commonPart' is \b incremented.
 */
void ComposedEdge::dispatchPerimeter(const std::set<int>& ids1, const std::set<int>& ids2, double& part1, double& part2, double& commonPart) const
{
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    {
      int curEdgeId=(*iter)->getPtr()->getId();
      double val=(*iter)->getPtr()->getCurveLength();
      if(ids1.find(curEdgeId)!=ids1.end())
        {
          if((*iter)->getLoc()!=FULL_ON_1)
            part1+=val;
          else
            commonPart+=val;
        }
      else
        if(ids2.find(curEdgeId)!=ids2.end())
          part2+=val;
        else
          commonPart+=val;
    }
}



/*!
 * For every subedges of 'this' if it comes from a part of father polygon the corresponding lgth is added in the corresponding 'result' rank.
 * If a sub-edge is found with a status equal to ON, the sub-edge is shared and the sum of all edges in the same case is returned.
 */
double ComposedEdge::dispatchPerimeterAdv(const ComposedEdge& father, std::vector<double>& result) const
{
  double ret=0;
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    ret+=father.checkFatherHood(*iter,result);
  return ret;
}

/*!
 * Givent edge 'edge' is one of subedge is father of this edge AND NOT ON on the corresponding rank in 'this' the curvelength is added.
 * If this sub edge is son and ON the corresponding length is returned by return param.
 */
double ComposedEdge::checkFatherHood(ElementaryEdge *edge, std::vector<double>& result) const
{
  double ret=0.;
  int rank=0;
  int refId=edge->getPtr()->getId();
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++,rank++)
    {
      int curId=(*iter)->getPtr()->getId();
      if(curId==refId)
        {
          if((*iter)->getLoc()!=FULL_ON_1)
            result[rank]+=edge->getCurveLength();
          else
            ret+=edge->getCurveLength();
          break;
        }
    }
  return ret;
}

void ComposedEdge::fillAllEdgeIds(std::set<int>& ids) const
{
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    ids.insert((*iter)->getPtr()->getId());
}

void ComposedEdge::getAllNodes(std::set<Node *>& output) const
{
  list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();
  for(;iter!=_subEdges.end();iter++)
    (*iter)->getAllNodes(output);
}

void ComposedEdge::getBarycenter(double *bary, double& weigh) const
{
  weigh=0.; bary[0]=0.; bary[1]=0.;
  double tmp1,tmp2[2];
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
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
  for(list<ElementaryEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    {
      ElementaryEdge *val=(*iter);
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
}

/*bool ComposedEdge::isInOrOut(Node *aNodeOn, Node *nodeToTest) const
{
  
  EdgeInfLin *e1=new EdgeInfLin(aNodeOn,nodeToTest);
  double ref=e1->getCharactValue(*nodeToTest);
  set< IntersectElement > inOutSwitch;
  for(vector<AbstractEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
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
  if(_subEdges.size()!=1)
    return false;
  return _subEdges.front()->intresincEqCoarse(other);
}

void ComposedEdge::clearAll(list<ElementaryEdge *>::iterator startToDel)
{
  for(list<ElementaryEdge *>::iterator iter=startToDel;iter!=_subEdges.end();iter++)
    delete (*iter);
}
