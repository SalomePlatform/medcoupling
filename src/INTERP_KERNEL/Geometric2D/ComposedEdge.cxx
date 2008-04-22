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
  _subEdges.resize(other._subEdges.size());
  int i=0;
  for(vector<AbstractEdge *>::const_iterator iter=other._subEdges.begin();iter!=other._subEdges.end();iter++,i++)
    _subEdges[i]=(*iter)->clone();
}

ComposedEdge::~ComposedEdge()
{
  clearAll(_subEdges.begin());
}

ElementaryEdge* &ComposedEdge::getLastElementary(IteratorOnComposedEdge::ItOnFixdLev &delta)
{
  AbstractEdge *e=_subEdges.back();
  ElementaryEdge *eCast=dynamic_cast< ElementaryEdge* >(e);
  if(eCast)
    return (ElementaryEdge* &)((AbstractEdge * &) _subEdges.back());
  delta++;
  return _subEdges.back()->getLastElementary(delta);
}

ElementaryEdge* &ComposedEdge::getFirstElementary(IteratorOnComposedEdge::ItOnFixdLev &delta)
{
  AbstractEdge *e=_subEdges.front();
  ElementaryEdge *eCast=dynamic_cast< ElementaryEdge* >(e);
  if(eCast)
    return (ElementaryEdge* &)((AbstractEdge * &) _subEdges.front());
  delta++;
  return _subEdges.front()->getFirstElementary(delta);
}

void ComposedEdge::setValueAt(int i, Edge *e, bool direction)
{
  delete _subEdges[i];
  _subEdges[i]=new ElementaryEdge(e,direction);
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

void ComposedEdge::pushBack(AbstractEdge *elem)
{
  _subEdges.push_back(elem);
}

void ComposedEdge::reverse()
{
  std::reverse(_subEdges.begin(),_subEdges.end());
  for(vector<AbstractEdge *>::iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    (*iter)->reverse();
}

int ComposedEdge::recursiveSize() const
{
  int ret=0;
  for(vector<AbstractEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    {
      ret+=(*iter)->recursiveSize();
    }
  return ret;
}

AbstractEdge *ComposedEdge::clone() const
{
  return new ComposedEdge(*this);
}

bool ComposedEdge::isNodeIn(Node *n) const
{
  bool ret=false;
  for(std::vector<AbstractEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end() && !ret;iter++)
    ret=(*iter)->isNodeIn(n);
  return ret;
}

double ComposedEdge::getAreaOfZone() const
{
  double ret=0.;
  for(std::vector<AbstractEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    ret+=(*iter)->getAreaOfZone();
  return ret;
}

void ComposedEdge::dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const
{
  stream.precision(10);
  for(std::vector<AbstractEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
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

AbstractEdge *ComposedEdge::simplify()
{
  if(size()!=1)
    return this;
  else
    {
      vector<AbstractEdge *>::iterator iter=_subEdges.begin();
      AbstractEdge *ret=*iter; iter++;
      clearAll(iter);
      _subEdges.clear();
      delete this;
      return ret;
    }
}

/*!
 * This method adds edge to this is it edge is localized as IN or ON.
 * The returned value is true if this is not empty and that edge is not IN or ON.
 */
bool ComposedEdge::addEdgeIfIn(ElementaryEdge *edge)
{
  if(edge->getLoc()!=FULL_OUT_1)
    {
      _subEdges.push_back(edge->clone());
      return false;
    }
  else if(!empty())
    {
      return true;
    }
  return false;
}

bool ComposedEdge::changeEndNodeWith(Node *node) const
{
  return _subEdges.back()->changeEndNodeWith(node);
}

bool ComposedEdge::changeStartNodeWith(Node *node) const
{
  return _subEdges.front()->changeStartNodeWith(node);
}

bool ComposedEdge::intresicEqual(const AbstractEdge *other) const
{
  throw Exception("ComposedEdge::intresicEqual : should never been called.");
}

bool ComposedEdge::intresicEqualDirSensitive(const AbstractEdge *other) const
{
  throw Exception("ComposedEdge::intresicEqualDirSensitive : should never been called.");
}

void ComposedEdge::fillBounds(Bounds& output) const
{
  vector<AbstractEdge *>::const_iterator iter=_subEdges.begin();
  for(;iter!=_subEdges.end();iter++)
    (*iter)->fillBounds(output);
}

void ComposedEdge::getAllNodes(std::set<Node *>& output) const
{
  vector<AbstractEdge *>::const_iterator iter=_subEdges.begin();
  for(;iter!=_subEdges.end();iter++)
    (*iter)->getAllNodes(output);
}

void ComposedEdge::getBarycenter(double *bary, double& weigh) const
{
  weigh=0.; bary[0]=0.; bary[1]=0.;
  double tmp1,tmp2[2];
  for(vector<AbstractEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
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

void ComposedEdge::clearAll(vector<AbstractEdge *>::iterator startToDel)
{
  for(vector<AbstractEdge *>::iterator iter=startToDel;iter!=_subEdges.end();iter++)
    delete (*iter);
}
