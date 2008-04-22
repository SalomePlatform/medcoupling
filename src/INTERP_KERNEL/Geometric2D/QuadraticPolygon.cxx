#include "QuadraticPolygon.hxx"
#include "ComposedEdgeWithIt.hxx"
#include "ElementaryEdge.hxx"
#include "EdgeArcCircle.hxx"
#include "EdgeLin.hxx"
#include "Bounds.hxx"
#include "Edge.txx"

#include <fstream>

using namespace std;
using namespace INTERP_KERNEL;

namespace INTERP_KERNEL
{
  const unsigned MAX_SIZE_OF_LINE_XFIG_FILE=1024;
}

QuadraticPolygon::QuadraticPolygon(const char *file)
{
  char currentLine[MAX_SIZE_OF_LINE_XFIG_FILE];
  ifstream stream(file);
  stream.exceptions(ios_base::eofbit);
  try
    {
      do
	stream.getline(currentLine,MAX_SIZE_OF_LINE_XFIG_FILE);
      while(strcmp(currentLine,"1200 2")!=0);
      do
	{
          Edge *newEdge=Edge::buildFromXfigLine(stream);
          if(!empty())
            newEdge->changeStartNodeWith(back()->getEndNode());
	  pushBack(newEdge);
	}
      while(1);
    }
  catch(ifstream::failure& e)
    {
    }
  front()->changeStartNodeWith(back()->getEndNode());
}

QuadraticPolygon::~QuadraticPolygon()
{
}

QuadraticPolygon *QuadraticPolygon::buildLinearPolygon(std::vector<Node *>& nodes)
{
  QuadraticPolygon *ret=new QuadraticPolygon;
  int size=nodes.size();
  for(int i=0;i<size;i++)
    {
      ret->pushBack(new EdgeLin(nodes[i],nodes[(i+1)%size]));
      nodes[i]->decrRef();
    }
  return ret;
}

QuadraticPolygon *QuadraticPolygon::buildArcCirclePolygon(std::vector<Node *>& nodes)
{
  QuadraticPolygon *ret=new QuadraticPolygon;
  int size=nodes.size();
  for(int i=0;i<size/2;i++)
    {
      EdgeLin *e1,*e2;
      e1=new EdgeLin(nodes[i],nodes[i+size/2]);
      e2=new EdgeLin(nodes[i+size/2],nodes[(i+1)%(size/2)]);
      SegSegIntersector inters(*e1,*e2);
      bool colinearity=inters.areColinears();
      delete e1; delete e2;
      if(colinearity)
        ret->pushBack(new EdgeLin(nodes[i],nodes[(i+1)%(size/2)]));
      else
        ret->pushBack(new EdgeArcCircle(nodes[i],nodes[i+size/2],nodes[(i+1)%(size/2)]));
      nodes[i]->decrRef(); nodes[i+size/2]->decrRef();
    }
  return ret;
}

void QuadraticPolygon::closeMe() const
{
  if(!front()->changeStartNodeWith(back()->getEndNode()))
    throw(Exception("big error: not closed polygon..."));
}

void QuadraticPolygon::circularPermute()
{
  vector<AbstractEdge *>::iterator iter1=_subEdges.begin();
  vector<AbstractEdge *>::iterator iter2=_subEdges.begin();
  if(iter2!=_subEdges.end())
    iter2++;
  else
    return ;
  for(;iter2!=_subEdges.end();iter1++,iter2++)
    iter_swap(iter1,iter2);
}

double QuadraticPolygon::getAreaFast() const
{
  double ret=0.;
  for(vector<AbstractEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    {
      ElementaryEdge *tmp=(ElementaryEdge *)(*iter);
      ret+=tmp->getAreaOfZoneFast();
    }
  return ret;
}

double QuadraticPolygon::getPerimeterFast() const
{
  double ret=0.;
  for(vector<AbstractEdge *>::const_iterator iter=_subEdges.begin();iter!=_subEdges.end();iter++)
    {
      ElementaryEdge *tmp=(ElementaryEdge *)(*iter);
      ret+=tmp->getCurveLength();
    }
  return ret;
}

double QuadraticPolygon::getHydraulicDiameter() const
{
  return 4*fabs(getAreaFast())/getPerimeterFast();
}

void QuadraticPolygon::dumpInXfigFileWithOther(const ComposedEdge& other, const char *fileName) const
{
  ofstream file(fileName);
  const int resolution=1200;
  Bounds box;
  box.prepareForAggregation();
  fillBounds(box);
  other.fillBounds(box);
  dumpInXfigFile(file,resolution,box);
  other.ComposedEdge::dumpInXfigFile(file,resolution,box);
}

void QuadraticPolygon::dumpInXfigFile(const char *fileName) const
{
  ofstream file(fileName);
  const int resolution=1200;
  Bounds box;
  box.prepareForAggregation();
  fillBounds(box);
  dumpInXfigFile(file,resolution,box);
}

void QuadraticPolygon::dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const
{
  stream << "#FIG 3.2  Produced by xfig version 3.2.5-alpha5" << endl;
  stream << "Landscape" << endl;
  stream << "Center" << endl;
  stream << "Metric" << endl;
  stream << "Letter" << endl;
  stream << "100.00" << endl;
  stream << "Single" << endl;
  stream << "-2" << endl;
  stream << resolution << " 2" << endl;
  ComposedEdge::dumpInXfigFile(stream,resolution,box);
}

double QuadraticPolygon::intersectWith(const QuadraticPolygon& other) const
{
  double ret=0;
  vector<QuadraticPolygon *> polygs=intersectMySelfWith(other);
  for(vector<QuadraticPolygon *>::iterator iter=polygs.begin();iter!=polygs.end();iter++)
    {
      ret+=fabs((*iter)->getAreaOfZone());
      delete *iter;
    }
  return ret;
}

std::vector<QuadraticPolygon *> QuadraticPolygon::intersectMySelfWith(const QuadraticPolygon& other) const
{
  QuadraticPolygon cpyOfThis(*this);
  QuadraticPolygon cpyOfOther(other); int nbOfSplits=0;
  splitPolygonsEachOther(cpyOfThis,cpyOfOther,nbOfSplits);
  //At this point cpyOfThis and cpyOfOther have been splited at maximum edge so that in/out can been done.
  performLocatingOperation(cpyOfOther);
  return other.buildIntersectionPolygons(cpyOfThis,cpyOfOther);
}

/*!
 * This method is typically the first step of boolean operations between pol1 and pol2.
 * This method perform the minimal splitting so that at the end each edges constituting pol1 are fully either IN or OUT or ON.
 * @param pol1 IN/OUT param that is equal to 'this' when called.
 */
void QuadraticPolygon::splitPolygonsEachOther(QuadraticPolygon& pol1, QuadraticPolygon& pol2, int& nbOfSplits) const
{
  IteratorOnComposedEdge it1(&pol1),it2(&pol2);
  MergePoints merge;
  ComposedEdge *c1=new ComposedEdge;
  ComposedEdgeWithIt *c2=new ComposedEdgeWithIt;
  for(it2.first();!it2.finished();it2.next())
    {
      ComposedEdgeWithIt *dealer=dynamic_cast<ComposedEdgeWithIt *>(it2.getLowestDealing());
      if(!dealer)
        it1.first();
      else
        it1=dealer->getIterator();
      for(;!it1.finished();)
        {
          ElementaryEdge* &curE2=it2.current();
          ElementaryEdge* &curE1=it1.current();
          merge.clear(); nbOfSplits++;
          if(curE1->getPtr()->intersectWith(curE2->getPtr(),merge,*c1,*c2))
            {
              if(!curE1->getDirection()) c1->reverse();
              if(!curE2->getDirection()) c2->reverse();
              AbstractEdge *c1s=c1->simplify();
              AbstractEdge *c2s=c2->simplify();
              updateNeighbours(merge,it1,it2,c1s,c2s);
              it1.next();//to do before
              //Substitution of simple edge by sub-edges.
              AbstractEdge **tmp1=(AbstractEdge**)&curE1; delete *tmp1; // <-- destroying simple edge coming from pol1
              AbstractEdge **tmp2=(AbstractEdge**)&curE2; delete *tmp2; // <-- destroying simple edge coming from pol2
              *tmp1=c1s;
              *tmp2=c2s;
              //
              if(c2s==c2)//in this case, all elts of c2s(ComposedEdges) should start to be intersected by starting to it1.
                c2->setIterator(it1);//c2s==c2 implies that c2s is a composed edge so sub iterations requested.
              c1=new ComposedEdge;
              c2=new ComposedEdgeWithIt;
            }
          else
	    {
	      updateNeighbours(merge,it1,it2,curE1,curE2);
	      it1.next();
	    }
        }
    }
  ComposedEdge::Delete(c1); delete c2;
}

void QuadraticPolygon::performLocatingOperation(QuadraticPolygon& pol2) const
{
  IteratorOnComposedEdge it(&pol2);
  TypeOfEdgeLocInPolygon loc=FULL_ON_1;
  for(it.first();!it.finished();it.next())
    {
      ElementaryEdge *cur=it.current();
      loc=cur->locateFullyMySelf(*this,loc);
    }
}

std::vector<QuadraticPolygon *> QuadraticPolygon::buildIntersectionPolygons(const QuadraticPolygon& pol1, const QuadraticPolygon& pol2) const
{
  vector<QuadraticPolygon *> ret;
  list<QuadraticPolygon *> pol2Zip=pol2.zipConsecutiveInSegments();
  if(!pol2Zip.empty())
    closePolygons(pol2Zip,pol1,ret);
  return ret;
}

std::list<QuadraticPolygon *> QuadraticPolygon::zipConsecutiveInSegments() const
{
  list<QuadraticPolygon *> ret;
  IteratorOnComposedEdge it((ComposedEdge *)this);
  int nbOfTurns=recursiveSize();
  int i=0;
  if(!it.goToNextInOn(false,i,nbOfTurns))
    return ret;
  i=0;
  //
  while(i<nbOfTurns)
    {
      QuadraticPolygon *tmp1=new QuadraticPolygon;
      TypeOfEdgeLocInPolygon loc=it.current()->getLoc();
      while(loc!=FULL_OUT_1 && i<nbOfTurns)
        {
	  AbstractEdge *tmp3=it.current()->clone();
          tmp1->pushBack(tmp3);
          it.nextLoop(); i++;
          loc=it.current()->getLoc();
        }
      if(tmp1->empty())
        {
          delete tmp1;
          continue;
        }
      ret.push_back(tmp1);
      it.goToNextInOn(true,i,nbOfTurns);
    }
  return ret;
}

/*!
 * 'this' should be considered as pol2Simplified.
 */
void QuadraticPolygon::closePolygons(std::list<QuadraticPolygon *>& pol2Zip, const QuadraticPolygon& pol1,
                                     std::vector<QuadraticPolygon *>& results) const
{
  bool directionKnownInPol1=false;
  bool directionInPol1;
  for(list<QuadraticPolygon *>::iterator iter=pol2Zip.begin();iter!=pol2Zip.end();)
    {
      if((*iter)->completed())
        {
          results.push_back(*iter);
          directionKnownInPol1=false;
          iter=pol2Zip.erase(iter);
          continue;
        }
      if(!directionKnownInPol1)
        if(!(*iter)->amIAChanceToBeCompletedBy(pol1,*this,directionInPol1))
          { delete *iter; iter=pol2Zip.erase(iter); continue; }
	else
	  directionKnownInPol1=true;
      list<QuadraticPolygon *>::iterator iter2=iter; iter2++;
      list<QuadraticPolygon *>::iterator iter3=(*iter)->fillAsMuchAsPossibleWith(pol1,iter2,pol2Zip.end(),directionInPol1);
      if(iter3!=pol2Zip.end())
        {
          (*iter)->pushBack(*iter3);
          pol2Zip.erase(iter3);
        }
    }
}

void QuadraticPolygon::updateNeighbours(const MergePoints& merger, IteratorOnComposedEdge it1, IteratorOnComposedEdge it2,
					const AbstractEdge *e1, const AbstractEdge *e2)
{
  it1.previousLoop(); it2.previousLoop();
  ElementaryEdge *curE1=it1.current(); ElementaryEdge *curE2=it2.current();
  curE1->changeEndNodeWith(e1->getStartNode()); curE2->changeEndNodeWith(e2->getStartNode());
  it1.nextLoop(); it1.nextLoop(); it2.nextLoop(); it2.nextLoop();
  curE1->changeStartNodeWith(e1->getEndNode()); curE2->changeStartNodeWith(e2->getEndNode());
}

bool QuadraticPolygon::amIAChanceToBeCompletedBy(const QuadraticPolygon& pol1Splitted,const QuadraticPolygon& pol2NotSplitted, bool& direction)
{
  IteratorOnComposedEdge it((QuadraticPolygon *)&pol1Splitted);
  bool found=false;
  Node *n=getEndNode();
  ElementaryEdge *cur=it.current();
  for(it.first();!it.finished() && !found;)
    {
      cur=it.current();
      found=(cur->getStartNode()==n);
      if(!found)
        it.next();
    }
  if(!found)
    throw Exception("Internal error : polygons uncompatible each others. Should never happend");
  //Ok we found correspondance between this and pol1. Searching for right direction to close polygon.
  IteratorOnComposedEdge::ItOnFixdLev tmp;
  ElementaryEdge *e=getLastElementary(tmp);
  if(e->getLoc()==FULL_ON_1)
    {
      if(e->getPtr()==cur->getPtr())
        {
          direction=false;
          it.previousLoop();
          cur=it.current();
          return pol2NotSplitted.isInOrOut(cur->getStartNode());
        }
      else
        {
          direction=true;
          return pol2NotSplitted.isInOrOut(cur->getEndNode());
        }
    }
  else
    direction=cur->locateFullyMySelfAbsolute(pol2NotSplitted)==FULL_IN_1;
  return true;
}

std::list<QuadraticPolygon *>::iterator QuadraticPolygon::fillAsMuchAsPossibleWith(const QuadraticPolygon& pol1Splitted,
                                                                                   std::list<QuadraticPolygon *>::iterator iStart,
                                                                                   std::list<QuadraticPolygon *>::iterator iEnd,
                                                                                   bool direction)
{
  IteratorOnComposedEdge it((QuadraticPolygon *)&pol1Splitted);
  bool found=false;
  Node *n=getEndNode();
  ElementaryEdge *cur;
  for(it.first();!it.finished() && !found;)
    {
      cur=it.current();
      found=(cur->getStartNode()==n);
      if(!found)
        it.next();
    }
  if(!direction)
    it.previousLoop();
  Node *nodeToTest;
  std::list<QuadraticPolygon *>::iterator ret;
  do
    {
      cur=it.current();
      AbstractEdge *tmp=cur->clone();
      if(!direction)
        tmp->reverse();
      pushBack(tmp);
      nodeToTest=tmp->getEndNode();
      direction?it.nextLoop():it.previousLoop();
      ret=checkInList(nodeToTest,iStart,iEnd);
      if(completed())
        return iEnd;
    }
  while(ret==iEnd);
  return ret;
}

std::list<QuadraticPolygon *>::iterator QuadraticPolygon::checkInList(Node *n, std::list<QuadraticPolygon *>::iterator iStart,
                                                                      std::list<QuadraticPolygon *>::iterator iEnd)
{
  for(list<QuadraticPolygon *>::iterator iter=iStart;iter!=iEnd;iter++)
    if((*iter)->isNodeIn(n))
      return iter;
  return iEnd;
}
