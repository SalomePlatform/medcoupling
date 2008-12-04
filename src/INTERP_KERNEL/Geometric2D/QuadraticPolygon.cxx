//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
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
#include "QuadraticPolygon.hxx"
#include "ElementaryEdge.hxx"
#include "EdgeArcCircle.hxx"
#include "AbstractEdge.hxx"
#include "EdgeLin.hxx"
#include "Bounds.hxx"
#include "Edge.txx"

#include <fstream>
#include <iomanip>

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

void QuadraticPolygon::buildDbgFile(const std::vector<Node *>& nodes, const char *fileName)
{
  ofstream file(fileName);
  file << setprecision(16);
  file << "  double coords[]=" << endl << "    { ";
  for(vector<Node *>::const_iterator iter=nodes.begin();iter!=nodes.end();iter++)
    {
      if(iter!=nodes.begin())
        file << "," << endl << "      ";
      file << (*(*iter))[0] << ", " << (*(*iter))[1];
    }
  file << "};" << endl;
}

void QuadraticPolygon::closeMe() const
{
  if(!front()->changeStartNodeWith(back()->getEndNode()))
    throw(Exception("big error: not closed polygon..."));
}

void QuadraticPolygon::circularPermute()
{
  if(_subEdges.size()>1)
    {
      ElementaryEdge *first=_subEdges.front();
      _subEdges.pop_front();
      _subEdges.push_back(first);
    }
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

/*!
 * Warning contrary to intersectWith method this method is \b NOT const. 'this' and 'other' are modified after call of this method.
 */
double QuadraticPolygon::intersectWithAbs(QuadraticPolygon& other)
{
  double ret=0.;
  double fact=normalize(&other);
  vector<QuadraticPolygon *> polygs=intersectMySelfWith(other);
  for(vector<QuadraticPolygon *>::iterator iter=polygs.begin();iter!=polygs.end();iter++)
    {
      ret+=fabs((*iter)->getArea());
      delete *iter;
    }
  return ret*fact*fact;
}

double QuadraticPolygon::intersectWith(const QuadraticPolygon& other) const
{
  double ret=0.;
  vector<QuadraticPolygon *> polygs=intersectMySelfWith(other);
  for(vector<QuadraticPolygon *>::iterator iter=polygs.begin();iter!=polygs.end();iter++)
    {
      ret+=fabs((*iter)->getArea());
      delete *iter;
    }
  return ret;
}

void QuadraticPolygon::intersectForPerimeter(const QuadraticPolygon& other, double& perimeterThisPart, double& perimeterOtherPart, double& perimeterCommonPart, double& area) const
{
  perimeterThisPart=0.; perimeterOtherPart=0.; perimeterCommonPart=0.; area=0.;
  vector<QuadraticPolygon *> polygs=intersectMySelfWith(other);
  set<int> idsFromThis;
  set<int> idsFromOther;
  fillAllEdgeIds(idsFromThis);
  other.fillAllEdgeIds(idsFromOther);
  for(vector<QuadraticPolygon *>::iterator iter=polygs.begin();iter!=polygs.end();iter++)
    {
      area+=fabs((*iter)->getArea());
      (*iter)->dispatchPerimeter(idsFromThis,idsFromOther,perimeterThisPart,perimeterOtherPart,perimeterCommonPart);
      delete *iter;
    }
}

/*!
 * polThis.size()==this->size() and polOther.size()==other.size().
 * For each ElementaryEdge of 'this', the corresponding UNIQUE contribution in resulting polygon is in 'polThis'.
 * For each ElementaryEdge of 'other', the corresponding UNIQUE contribution in resulting polygon is in 'polOther'.
 * The returned value is the sum of common edge lengths.
 */
double QuadraticPolygon::intersectForPerimeterAdvanced(const QuadraticPolygon& other, std::vector< double >& polThis, std::vector< double >& polOther, double& area) const
{
  area=0.;
  double ret=0.;
  polThis.resize(size());
  polOther.resize(other.size());
  vector<QuadraticPolygon *> polygs=intersectMySelfWith(other);
  for(vector<QuadraticPolygon *>::iterator iter=polygs.begin();iter!=polygs.end();iter++)
    {
      area+=fabs((*iter)->getArea());
      ret+=(*iter)->dispatchPerimeterAdv(*this,polThis);
      ret+=(*iter)->dispatchPerimeterAdv(other,polOther);
      delete *iter;
    }
  return ret;
}

/*!
 * numberOfCreatedPointsPerEdge is resized to the number of edges of 'this'.
 * This method returns in ordered maner the number of newly created points per edge.
 * This method performs a split process between 'this' and 'other' that gives the result PThis.
 * Then for each edges of 'this' this method counts how many edges in Pthis have the same id.
 */
void QuadraticPolygon::intersectForPoint(const QuadraticPolygon& other, std::vector< int >& numberOfCreatedPointsPerEdge) const
{
  numberOfCreatedPointsPerEdge.resize(size());
  QuadraticPolygon cpyOfThis(*this);
  QuadraticPolygon cpyOfOther(other); int nbOfSplits=0;
  splitPolygonsEachOther(cpyOfThis,cpyOfOther,nbOfSplits);
  dispatchForNode(cpyOfThis,numberOfCreatedPointsPerEdge);
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
  ComposedEdge *c2=new ComposedEdge;
  for(it2.first();!it2.finished();it2.next())
    {
      ElementaryEdge* curE2=it2.current();
      if(!curE2->isThereStartPoint())
        it1.first();
      else
        it1=curE2->getIterator();
      for(;!it1.finished();)
        {
          
          ElementaryEdge* curE1=it1.current();
          merge.clear(); nbOfSplits++;
          if(curE1->getPtr()->intersectWith(curE2->getPtr(),merge,*c1,*c2))
            {
              if(!curE1->getDirection()) c1->reverse();
              if(!curE2->getDirection()) c2->reverse();
              updateNeighbours(merge,it1,it2,c1,c2);
              //Substitution of simple edge by sub-edges.
              delete curE1; // <-- destroying simple edge coming from pol1
              delete curE2; // <-- destroying simple edge coming from pol2
              it1.insertElemEdges(c1,true);// <-- 2nd param is true to go next.
              it2.insertElemEdges(c2,false);// <-- 2nd param is false to avoid to go next.
              curE2=it2.current();
              //
              it1.assignMySelfToAllElems(c2);//To avoid that others
              SoftDelete(c1);
              SoftDelete(c2);
              c1=new ComposedEdge;
              c2=new ComposedEdge;
            }
          else
	    {
	      updateNeighbours(merge,it1,it2,curE1,curE2);
	      it1.next();
	    }
        }
    }
  Delete(c1);
  Delete(c2);
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
	  ElementaryEdge *tmp3=it.current()->clone();
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
          SoftDelete(*iter3);
          pol2Zip.erase(iter3);
        }
    }
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
  ElementaryEdge *e=_subEdges.back();
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
	  Node *repr=cur->getPtr()->buildRepresentantOfMySelf();
          bool ret=pol2NotSplitted.isInOrOut(repr);
	  repr->decrRef();
	  return ret;
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
      ElementaryEdge *tmp=cur->clone();
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
