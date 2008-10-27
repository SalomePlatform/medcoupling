#include "Edge.hxx"
#include "EdgeLin.hxx"
#include "EdgeInfLin.hxx"
//#include "EdgeParabol.hxx"
#include "EdgeArcCircle.hxx"
#include "InterpolationUtils.hxx"

using namespace std;
using namespace INTERP_KERNEL;

int Edge::_idCounter=0;

MergePoints::MergePoints():_ass1Start1(0),_ass1End1(0),_ass1Start2(0),_ass1End2(0),
                           _ass2Start1(0),_ass2End1(0),_ass2Start2(0),_ass2End2(0)
{
}

void MergePoints::start1Replaced()
{
  unsigned nbOfAsso=getNumberOfAssociations();
  if(nbOfAsso==0)
    _ass1Start1=1;
  else
    _ass2Start1=1;
}

void MergePoints::end1Replaced()
{
  unsigned nbOfAsso=getNumberOfAssociations();
  if(nbOfAsso==0)
    _ass1End1=1;
  else
    _ass2End1=1;
}

void MergePoints::start1OnStart2()
{
  unsigned nbOfAsso=getNumberOfAssociations();
  if(nbOfAsso==0)
    {
      _ass1Start1=1;
      _ass1Start2=1;
    }
  else
    {
      _ass2Start1=1;
      _ass2Start2=1;
    }
}

void MergePoints::start1OnEnd2()
{
  unsigned nbOfAsso=getNumberOfAssociations();
  if(nbOfAsso==0)
    {
      _ass1Start1=1;
      _ass1End2=1;
    }
  else
    {
      _ass2Start1=1;
      _ass2End2=1;
    }
}

void MergePoints::end1OnStart2()
{
  unsigned nbOfAsso=getNumberOfAssociations();
  if(nbOfAsso==0)
    {
      _ass1End1=1;
      _ass1Start2=1;
    }
  else
    {
      _ass2End1=1;
      _ass2Start2=1;
    }
}

void MergePoints::end1OnEnd2()
{
  unsigned nbOfAsso=getNumberOfAssociations();
  if(nbOfAsso==0)
    {
      _ass1End1=1;
      _ass1End2=1;
    }
  else
    {
      _ass2End1=1;
      _ass2End2=1;
    }
}

bool MergePoints::isStart1(unsigned rk) const
{
  if(rk==0)
    return _ass1Start1;
  else
    return _ass2Start1;
}

bool MergePoints::isEnd1(unsigned rk) const
{
  if(rk==0)
    return _ass1End1;
  else
    return _ass2End1;
}

bool MergePoints::isStart2(unsigned rk) const
{
  if(rk==0)
    return _ass1Start2;
  else
    return _ass2Start2;
}

bool MergePoints::isEnd2(unsigned rk) const
{
  if(rk==0)
    return _ass1End2;
  else
    return _ass2End2;
}

void MergePoints::clear()
{
  _ass1Start1=0;_ass1End1=0;_ass1Start2=0;_ass1End2=0;
  _ass2Start1=0;_ass2End1=0;_ass2Start2=0;_ass2End2=0;
}

unsigned MergePoints::getNumberOfAssociations() const
{
  unsigned ret=0;
  unsigned subTot=_ass1Start1+_ass1End1+_ass1Start2+_ass1End2;
  if(subTot!=0)
    ret++;
  subTot=_ass2Start1+_ass2End1+_ass2Start2+_ass2End2;
  if(subTot!=0)
    ret++;
  return ret;
}

IntersectElement::IntersectElement(double val1, double val2, bool start1, bool end1, bool start2, bool end2, Node *node
                                   , const Edge& e1, const Edge& e2, bool keepOrder):_1S(keepOrder?start1:start2),
                                                                                     _1E(keepOrder?end1:end2),
                                                                                     _2S(keepOrder?start2:start1),
                                                                                     _2E(keepOrder?end2:end1),
                                                                                     _chararctValForE1(keepOrder?val1:val2),
                                                                                     _chararctValForE2(keepOrder?val2:val1),
                                                                                     _node(node),_e1(keepOrder?e1:e2),
                                                                                     _e2(keepOrder?e2:e1)
{
}

IntersectElement::IntersectElement(const IntersectElement& other):_1S(other._1S),_1E(other._1E),_2S(other._2S),_2E(other._2E),
                                                                  _chararctValForE1(other._chararctValForE1),
                                                                  _chararctValForE2(other._chararctValForE2),_node(other._node),
                                                                  _e1(other._e1), _e2(other._e2)
{
  if(_node)
    _node->incrRef();
}

IntersectElement& IntersectElement::operator=(const IntersectElement& other)
{
  _1S=other._1S;_1E=other._1E; _2S=other._2S; _2E=other._2E;
  _chararctValForE1=other._chararctValForE1;
  _chararctValForE2=other._chararctValForE2;
  setNode(other._node);
  return *this;
}

bool IntersectElement::operator<(const IntersectElement& other) const
{
  return _e1.isLower(_chararctValForE1,other._chararctValForE1);
}

IntersectElement::~IntersectElement()
{
  if(_node)
    _node->decrRef();
}

/*!
 * Returns 0 or 1.
 */
bool IntersectElement::isOnMergedExtremity() const
{
  if( (_1S && _2S) || (_1S && _2E) || (_1E && _2S) || (_1E && _2E) )
    return true;
  return false;
}

/*!
 * To call if isOnMergedExtremity returned true.
 */
void IntersectElement::performMerging(MergePoints& commonNode) const
{
  if(_1S && _2S)
    {
      if(_e1.changeStartNodeWith(_e2.getStartNode()))
        {
          _e2.getStartNode()->declareOnLim();
          commonNode.start1OnStart2();
        }
    }
  else if(_1S && _2E)
    {
      if(_e1.changeStartNodeWith(_e2.getEndNode()))
        {
          _e2.getEndNode()->declareOnLim();
          commonNode.start1OnEnd2();
        }
    }
  else if(_1E && _2S)
    {
      if(_e1.changeEndNodeWith(_e2.getStartNode()))
        {
          _e2.getStartNode()->declareOnLim();
          commonNode.end1OnStart2();
        }
    }
  else if(_1E && _2E)
    {
      if(_e1.changeEndNodeWith(_e2.getEndNode()))
        {
          _e2.getEndNode()->declareOnLim();
          commonNode.end1OnEnd2();
        }
    }
}

/*!
 * This methode is const because 'node' is supposed to be equal geomitrically to _node.
 */
void IntersectElement::setNode(Node *node) const
{
  if(node!=_node)
    {
      if(_node)
        ((Node *)_node)->decrRef();
      ((IntersectElement *)(this))->_node=node;
      if(_node)
        _node->incrRef();
    }
}

bool IntersectElement::isLowerOnOther(const IntersectElement& other) const
{
  return _e2.isLower(_chararctValForE2,other._chararctValForE2);
}

unsigned IntersectElement::isOnExtrForAnEdgeAndInForOtherEdge() const
{
  if(( _1S && !(_2S || _2E) ) || ( _1E && !(_2S || _2E) ))
    {
      if(_1S && !(_2S || _2E))
        setNode(_e1.getStartNode());
      else
        setNode(_e1.getEndNode());
      if(_e2.isIn(_chararctValForE2))
        return LIMIT_ON;
      return LIMIT_ALONE;
    }
  if(( _2S && !(_1S || _1E) ) || ( _2E && !(_1S || _1E)))
    {
      if(_2S && !(_1S || _1E))
        setNode(_e2.getStartNode());
      else
        setNode(_e2.getEndNode());
      if(_e1.isIn(_chararctValForE1))
        return LIMIT_ON;
      return LIMIT_ALONE;
    }
  return NO_LIMIT;
}

bool IntersectElement::isIncludedByBoth() const
{
  return _e1.isIn(_chararctValForE1) && _e2.isIn(_chararctValForE2);
}
  
bool Intersector::intersect(const Bounds *whereToFind, std::vector<Node *>& newNodes, bool& order, MergePoints& commonNode)
{
  list< IntersectElement > listOfIntesc=getIntersectionsCharacteristicVal();
  list< IntersectElement >::iterator iter;
  for(iter=listOfIntesc.begin();iter!=listOfIntesc.end();)
    {
      if((*iter).isOnMergedExtremity())
        {
          (*iter).performMerging(commonNode);
          iter=listOfIntesc.erase(iter);
          continue;
        }
      unsigned tmp=(*iter).isOnExtrForAnEdgeAndInForOtherEdge();
      if(tmp==IntersectElement::LIMIT_ALONE)
        {
          iter=listOfIntesc.erase(iter);
          continue;
        }
      else if(tmp==IntersectElement::LIMIT_ON)
        {
	  iter++;
	  continue;
	}
      if(!(*iter).isIncludedByBoth())
        {
          iter=listOfIntesc.erase(iter);
          continue;
        }
      iter++;
    }
  if(listOfIntesc.size()==0)
    return false;
  if(listOfIntesc.size()==1)
    {
      order=true;//useless
      newNodes.push_back(listOfIntesc.front().getNodeAndReleaseIt());
    }
  else
    {
      vector<IntersectElement> vecOfIntesc(listOfIntesc.begin(),listOfIntesc.end());
      listOfIntesc.clear();
      sort(vecOfIntesc.begin(),vecOfIntesc.end());
      for(vector<IntersectElement>::iterator iterV=vecOfIntesc.begin();iterV!=vecOfIntesc.end();iterV++)
        newNodes.push_back((*iterV).getNodeAndReleaseIt());
      order=vecOfIntesc.front().isLowerOnOther(vecOfIntesc.back());
    }
  return true;
}

void Intersector::obviousCaseForCurvAbscisse(Node *node, TypeOfLocInEdge& where, MergePoints& commonNode, bool& obvious) const
{
  obvious=true;
  if(node->isEqual(*_e1.getStartNode()))
    {
      where=START;
      if(_e1.changeStartNodeWith(node))
	{
	  commonNode.start1Replaced();
	  node->declareOnLim();
	}
      return ;
    }
  if(node->isEqual(*_e1.getEndNode()))
    {
      where=END;
      if(_e1.changeEndNodeWith(node))
	{
	  commonNode.end1Replaced();
	  node->declareOnLim();
	}
      return ;
    }
  obvious=false;
}

Edge::Edge(double sX, double sY, double eX, double eY):_cnt(1),_loc(FULL_UNKNOWN),_start(new Node(sX,sY)),_end(new Node(eX,eY)),_id(_idCounter++)
{
}

Edge::~Edge()
{
  _start->decrRef();
  if(_end)
    _end->decrRef();
}

bool Edge::decrRef()
{
  bool ret=(--_cnt==0);
  if(ret)
    delete this;
  return ret;
}

void Edge::declareOn() const
{
  if(_loc==FULL_UNKNOWN)
    {
      _loc=FULL_ON_1;
      _start->declareOn();
      _end->declareOn();
    }
}

void Edge::declareIn() const
{
  if(_loc==FULL_UNKNOWN)
    {
      _loc=FULL_IN_1;
      _start->declareIn();
      _end->declareIn();
    }
}

void Edge::declareOut() const
{
  if(_loc==FULL_UNKNOWN)
    {
      _loc=FULL_OUT_1;
      _start->declareOut();
      _end->declareOut();
    }
}

void Edge::fillXfigStreamForLoc(std::ostream& stream) const
{
  switch(_loc)
    {
    case FULL_IN_1:
      stream << '2';//Green
      break;
    case FULL_OUT_1:
      stream << '1';//Bleue
      break;
    case FULL_ON_1:
      stream << '4';//Red
      break;
    default:
      stream << '0';
    }
}

bool Edge::changeStartNodeWith(Node *otherStartNode) const
{
  if(_start==otherStartNode)
    return true;
  if(_start->isEqual(*otherStartNode))
    {
      (((Edge *)this)->_start)->decrRef();//un-const cast Ok thanks to 2 lines above.
      (((Edge *)this)->_start)=otherStartNode;
      _start->incrRef();
      return true;
    }
  return false;
}

bool Edge::changeStartNodeWithAndKeepTrack(Node *otherStartNode, std::vector<Node *>& track) const
{
  if(_start==otherStartNode)
    return true;
  if(_start->isEqualAndKeepTrack(*otherStartNode,track))
    {
      (((Edge *)this)->_start)->decrRef();//un-const cast Ok thanks to 2 lines above.
      (((Edge *)this)->_start)=otherStartNode;
      otherStartNode->incrRef();
      return true;
    }
  return false;
}

bool Edge::changeEndNodeWith(Node *otherEndNode) const
{
  if(_end==otherEndNode)
    return true;
  if(_end->isEqual(*otherEndNode))
    {
      (((Edge *)this)->_end)->decrRef();
      (((Edge *)this)->_end)=otherEndNode;
      _end->incrRef();
      return true;
    }
  return false;
}

bool Edge::changeEndNodeWithAndKeepTrack(Node *otherEndNode, std::vector<Node *>& track) const
{
  if(_end==otherEndNode)
    return true;
  if(_end->isEqualAndKeepTrack(*otherEndNode,track))
    {
      (((Edge *)this)->_end)->decrRef();
      (((Edge *)this)->_end)=otherEndNode;
      otherEndNode->incrRef();
      return true;
    }
  return false;
}

/*!
 * Precondition : 'start' and 'end' are lying on the same curve than 'this'.
 * Add in vec the sub edge lying on this.
 * If 'start' is equal (by pointer) to '_end' and 'end' is equal to '_end' too nothing is added.
 * If 'start' is equal (by pointer) to '_start' and 'end' is equal to '_start' too nothing is added.
 * If 'start' is equal (by pointer) to '_start' and 'end' is equal to '_end' this is added in vec.
 */
void Edge::addSubEdgeInVector(Node *start, Node *end, ComposedEdge& vec) const
{
  if((start==_start && end==_start) || (start==_end && end==_end))
    return ;
  if(start==_start && end==_end)
    {
      incrRef();
      vec.pushBack((Edge *)this);
      return ;
    }
  vec.pushBack(buildEdgeLyingOnMeWithId(start,end,true));
}

/*!
 * Retrieves a vector 'vectOutput' that is normal to 'this'. 'vectOutput' is normalized.
 */
void Edge::getNormalVector(double *vectOutput) const
{
  copy((const double *)(*_end),(const double *)(*_end)+2,vectOutput);
  transform(vectOutput,vectOutput+2,(const double *)(*_start),vectOutput,minus<double>());
  double norm=1./Node::norm(vectOutput);
  transform(vectOutput,vectOutput+2,vectOutput,bind2nd(multiplies<double>(),norm));
  double tmp=vectOutput[0];
  vectOutput[0]=vectOutput[1];
  vectOutput[1]=-tmp;
}

Edge *Edge::buildEdgeFrom(Node *start, Node *end)
{
  return new EdgeLin(start,end);
}

Edge *Edge::buildFromXfigLine(std::istream& str)
{
  unsigned char type;
  str >> type;
  if(type=='2')
    return new EdgeLin(str);
  else if(type=='5')
    return new EdgeArcCircle(str);
  else
    {
      std::cerr << "Unknown line found...";
      return 0;
    }
}

Edge *Edge::buildEdgeLyingOnMeWithId(Node *start, Node *end, bool direction) const
{
  Edge *ret=buildEdgeLyingOnMe(start,end,direction);
  ret->setId(_id);
  return ret;
}

/*!
 * \param other The Edge with which we are going to intersect.
 * \param commonNode Output. The common nodes found during operation of intersecting.
 * \param outVal1 Output filled in case true is returned. It specifies the new or not new edges by which 'this' is replaced after intersecting op.
 * \param outVal2 Output filled in case true is returned. It specifies the new or not new edges by which 'other' is replaced after intersecting op.
 * return true if the intersection between this.
 */
bool Edge::intersectWith(const Edge *other, MergePoints& commonNode,
                         ComposedEdge& outVal1, ComposedEdge& outVal2) const
{
  bool ret=true;
  Bounds *merge=_bounds.nearlyAmIIntersectingWith(other->getBounds());
  if(!merge)
    return false;
  delete merge;
  merge=0;
  Intersector *intersector=buildIntersectorWith(this,other);
  ret=intersect(this,other,intersector,merge,commonNode,outVal1,outVal2);
  delete intersector;
  return ret;
}

bool Edge::intersectOverlapped(const Edge *f1, const Edge *f2, Intersector *intersector, MergePoints& commonNode,
                               ComposedEdge& outValForF1, ComposedEdge& outValForF2)
{
  bool rev=intersector->haveTheySameDirection();
  Node *f2Start=f2->getNode(rev?START:END);
  Node *f2End=f2->getNode(rev?END:START);
  TypeOfLocInEdge place1, place2;
  intersector->getPlacements(f2Start,f2End,place1,place2,commonNode);
  int codeForIntersectionCase=combineCodes(place1,place2);
  return splitOverlappedEdges(f1,f2,f2Start,f2End,rev,codeForIntersectionCase,outValForF1,outValForF2);
}

/*!
 * Perform 1D linear interpolation. Warning distrib1 and distrib2 are expected to be in ascending mode.
 */
void Edge::interpolate1DLin(const std::vector<double>& distrib1, const std::vector<double>& distrib2, std::map<int, std::map<int,double> >& result)
{
  int nbOfV1=distrib1.size()-1;
  int nbOfV2=distrib2.size()-1;
  Node *n1=new Node(0.,0.); Node *n3=new Node(0.,0.);
  Node *n2=new Node(0.,0.); Node *n4=new Node(0.,0.);
  MergePoints commonNode;
  for(int i=0;i<nbOfV1;i++)
    {
      vector<double>::const_iterator iter=find_if(distrib2.begin()+1,distrib2.end(),bind2nd(greater_equal<double>(),distrib1[i]));
      if(iter!=distrib2.end())
        {
          for(int j=(iter-1)-distrib2.begin();j<nbOfV2;j++)
            {
              if(distrib2[j]<=distrib1[i+1])
                {
                  EdgeLin *e1=new EdgeLin(n1,n2); EdgeLin *e2=new EdgeLin(n3,n4);
                  n1->setNewCoords(distrib1[i],0.); n2->setNewCoords(distrib1[i+1],0.);
                  n3->setNewCoords(distrib2[j],0.); n4->setNewCoords(distrib2[j+1],0.);
                  ComposedEdge *f1=new ComposedEdge;
                  ComposedEdge *f2=new ComposedEdge;
                  SegSegIntersector inters(*e1,*e2);
                  bool b1,b2;
                  inters.areOverlappedOrOnlyColinears(0,b1,b2);
                  if(intersectOverlapped(e1,e2,&inters,commonNode,*f1,*f2))
                    {
                      result[i][j]=f1->getCommonLengthWith(*f2)/e1->getCurveLength();
                    }
                  ComposedEdge::Delete(f1); ComposedEdge::Delete(f2);
                  e1->decrRef(); e2->decrRef();
                }
            }
        }
    }
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef();
}

Intersector *Edge::buildIntersectorWith(const Edge *e1, const Edge *e2)
{
  Intersector *ret=0;
  const EdgeLin *tmp1=0;
  const EdgeArcCircle *tmp2=0;
  unsigned char type1=e1->getTypeOfFunc();
  e1->dynCastFunction(tmp1,tmp2);
  unsigned char type2=e2->getTypeOfFunc();
  e2->dynCastFunction(tmp1,tmp2);
  type1|=type2;
  switch(type1)
    {
    case 1:// Intersection seg/seg
      ret=new SegSegIntersector((const EdgeLin &)(*e1),(const EdgeLin &)(*e2));
      break;
    case 5:// Intersection seg/arc of circle
      ret=new ArcCSegIntersector(*tmp2,*tmp1,tmp2==e1);
      break;
    case 4:// Intersection arc/arc of circle
      ret=new ArcCArcCIntersector((const EdgeArcCircle &)(*e1),(const EdgeArcCircle &)(*e2));
      break;
    default:
      //Should never happen
      throw Exception("A non managed association of edge has been detected. Go work for intersection computation implementation.");
    }
  return ret;
}

/*!
 * See Node::applySimilarity to see signification of params.
 */
void Edge::applySimilarity(double xBary, double yBary, double dimChar)
{
  _bounds.applySimilarity(xBary,yBary,dimChar);
}

bool Edge::intersect(const Edge *f1, const Edge *f2, Intersector *intersector, const Bounds *whereToFind, MergePoints& commonNode,
                     ComposedEdge& outValForF1, ComposedEdge& outValForF2)
{
  bool obviousNoIntersection;
  bool areOverlapped;
  intersector->areOverlappedOrOnlyColinears(whereToFind,obviousNoIntersection,areOverlapped);
  if(areOverlapped)
    return intersectOverlapped(f1,f2,intersector,commonNode,outValForF1,outValForF2);
  if(obviousNoIntersection)
    return false;
  vector<Node *> newNodes;
  bool order;
  if(intersector->intersect(whereToFind,newNodes,order,commonNode))
    {
      if(newNodes.empty())
        throw Exception("Internal error occured - error in intersector implementation!");// This case should never happen
      vector<Node *>::iterator iter=newNodes.begin();
      vector<Node *>::reverse_iterator iterR=newNodes.rbegin();
      f1->addSubEdgeInVector(f1->getStartNode(),*iter,outValForF1);
      f2->addSubEdgeInVector(f2->getStartNode(),*iter,outValForF2);
      for(vector<Node *>::iterator iter=newNodes.begin();iter!=newNodes.end();iter++,iterR++)
        {
          if((iter+1)==newNodes.end())
            {
              f1->addSubEdgeInVector(*iter,f1->getEndNode(),outValForF1);
              (*iter)->decrRef();
              f2->addSubEdgeInVector(order?*iter:*iterR,f2->getEndNode(),outValForF2);
            }
          else
            {
              f1->addSubEdgeInVector(*iter,*(iter+1),outValForF1);
              (*iter)->decrRef();
              f2->addSubEdgeInVector(order?*iter:*iterR,order?*(iter+1):*(iterR+1),outValForF2);
            }
        }
      return true;
    }
  else//no intersection inside whereToFind
    return false;
}

int Edge::combineCodes(TypeOfLocInEdge code1, TypeOfLocInEdge code2)
{
  int ret=(int)code1;
  ret*=OFFSET_FOR_TYPEOFLOCINEDGE;
  ret+=(int)code2;
  return ret;
}

/*!
 * This method splits e1 and e2 into pieces as much sharable as possible. The precondition to the call of this method
 * is that e1 and e2 have been declared as overlapped by corresponding intersector built from e1 and e2 type.
 *
 * @param nS start node of e2 with the SAME DIRECTION as e1. The pointer nS should be equal to start node of e2 or to its end node.
 * @param nE end node of e2 with the SAME DIRECTION as e1. The pointer nE should be equal to start node of e2 or to its end node.
 * @param direction is param that specifies if e2 and e1 have same directions (true) or opposed (false).
 * @param code is the code returned by method Edge::combineCodes.
 */
bool Edge::splitOverlappedEdges(const Edge *e1, const Edge *e2, Node *nS, Node *nE, bool direction, int code,
                                ComposedEdge& outVal1, ComposedEdge& outVal2)
{
  Edge *tmp;
  switch(code)
    {
    case OUT_BEFORE*OFFSET_FOR_TYPEOFLOCINEDGE+START:      // OUT_BEFORE - START
    case OUT_BEFORE*OFFSET_FOR_TYPEOFLOCINEDGE+OUT_BEFORE: // OUT_BEFORE - OUT_BEFORE
    case OUT_AFTER*OFFSET_FOR_TYPEOFLOCINEDGE+OUT_AFTER:   // OUT_AFTER - OUT_AFTER
    case END*OFFSET_FOR_TYPEOFLOCINEDGE+OUT_AFTER:         // END - OUT_AFTER
      return false;
    case INSIDE*OFFSET_FOR_TYPEOFLOCINEDGE+OUT_AFTER:      // INSIDE - OUT_AFTER
      outVal1.pushBack(e1->buildEdgeLyingOnMe(e1->getStartNode(),nS,true));
      tmp=e1->buildEdgeLyingOnMe(nS,e1->getEndNode()); tmp->incrRef();
      outVal1.pushBack(tmp);
      outVal2.resize(2);
      outVal2.setValueAt(direction?0:1,tmp,direction); tmp->declareOn();
      outVal2.setValueAt(direction?1:0,e1->buildEdgeLyingOnMe(e1->getEndNode(),nE,direction));
      return true;
    case INSIDE*OFFSET_FOR_TYPEOFLOCINEDGE+INSIDE:         // INSIDE - INSIDE
      e2->incrRef(); e2->incrRef();
      outVal1.resize(3);
      outVal1.setValueAt(0,e1->buildEdgeLyingOnMe(e1->getStartNode(),nS));
      outVal1.setValueAt(1,(Edge*)e2,direction);
      outVal1.setValueAt(2,e1->buildEdgeLyingOnMe(nE,e1->getEndNode()));
      outVal2.pushBack((Edge*)e2); e2->declareOn();
      return true;
    case OUT_BEFORE*OFFSET_FOR_TYPEOFLOCINEDGE+INSIDE:     // OUT_BEFORE - INSIDE
      tmp=e1->buildEdgeLyingOnMe(e1->getStartNode(),nE); tmp->incrRef();
      outVal1.pushBack(tmp);
      outVal1.pushBack(e1->buildEdgeLyingOnMe(nE,e1->getEndNode()));
      outVal2.resize(2);
      outVal2.setValueAt(direction?0:1,e1->buildEdgeLyingOnMe(nS,e1->getStartNode(),direction));
      outVal2.setValueAt(direction?1:0,tmp,direction); tmp->declareOn();
      return true;
    case OUT_BEFORE*OFFSET_FOR_TYPEOFLOCINEDGE+OUT_AFTER:  // OUT_BEFORE - OUT_AFTER
      e1->incrRef(); e1->incrRef();
      outVal1.pushBack((Edge*)e1);
      outVal2.resize(3);
      outVal2.setValueAt(direction?0:2,e1->buildEdgeLyingOnMe(nS,e1->getStartNode(),direction));
      outVal2.setValueAt(1,(Edge*)e1,direction); e1->declareOn();
      outVal2.setValueAt(direction?2:0,e1->buildEdgeLyingOnMe(e1->getEndNode(),nE,direction));
      return true;
    case START*OFFSET_FOR_TYPEOFLOCINEDGE+END:             // START - END
      e1->incrRef(); e1->incrRef();
      outVal1.pushBack((Edge*)e1);
      outVal2.pushBack((Edge*)e1,direction); e1->declareOn();
      return true;
    case START*OFFSET_FOR_TYPEOFLOCINEDGE+OUT_AFTER:       // START - OUT_AFTER
      e1->incrRef(); e1->incrRef();
      outVal1.pushBack((Edge*)e1);
      outVal2.resize(2);
      outVal2.setValueAt(direction?0:1,(Edge*)e1,direction); e1->declareOn();
      outVal2.setValueAt(direction?1:0,e1->buildEdgeLyingOnMe(e1->getEndNode(),nE,direction));
      return true;
    case INSIDE*OFFSET_FOR_TYPEOFLOCINEDGE+END:            // INSIDE - END
      e2->incrRef(); e2->incrRef();
      outVal1.pushBack(e1->buildEdgeLyingOnMe(e1->getStartNode(),nS,true));
      outVal1.pushBack((Edge*)e2,direction);
      outVal2.pushBack((Edge*)e2); e2->declareOn();
      return true;
    case OUT_BEFORE*OFFSET_FOR_TYPEOFLOCINEDGE+END:        // OUT_BEFORE - END
      e1->incrRef(); e1->incrRef();
      outVal1.pushBack((Edge*)e1);
      outVal2.resize(2);
      outVal2.setValueAt(direction?0:1,e1->buildEdgeLyingOnMe(nS,e1->getStartNode(),direction));
      outVal2.setValueAt(direction?1:0,(Edge*)e1,direction); e1->declareOn();
      return true;
    case START*OFFSET_FOR_TYPEOFLOCINEDGE+INSIDE:          // START - INSIDE
      e2->incrRef(); e2->incrRef();
      outVal1.pushBack((Edge*)e2,direction);
      outVal1.pushBack(e1->buildEdgeLyingOnMe(nE,e1->getEndNode()));
      outVal2.pushBack((Edge*)e2); e2->declareOn();
      return true;
    default:
      throw Exception("Unexpected situation of overlapping edges : internal error occurs ! ");
    }
}
