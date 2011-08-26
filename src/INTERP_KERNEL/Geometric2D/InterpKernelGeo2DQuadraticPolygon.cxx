// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#include "InterpKernelGeo2DQuadraticPolygon.hxx"
#include "InterpKernelGeo2DElementaryEdge.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DAbstractEdge.hxx"
#include "InterpKernelGeo2DEdgeLin.hxx"
#include "InterpKernelGeo2DBounds.hxx"
#include "InterpKernelGeo2DEdge.txx"

#include <fstream>
#include <iomanip>
#include <cstring>
#include <limits>

using namespace INTERP_KERNEL;

namespace INTERP_KERNEL
{
  const unsigned MAX_SIZE_OF_LINE_XFIG_FILE=1024;
}

QuadraticPolygon::QuadraticPolygon(const char *file)
{
  char currentLine[MAX_SIZE_OF_LINE_XFIG_FILE];
  std::ifstream stream(file);
  stream.exceptions(std::ios_base::eofbit);
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
  catch(std::ifstream::failure&)
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
  std::ofstream file(fileName);
  file << std::setprecision(16);
  file << "  double coords[]=" << std::endl << "    { ";
  for(std::vector<Node *>::const_iterator iter=nodes.begin();iter!=nodes.end();iter++)
    {
      if(iter!=nodes.begin())
        file << "," << std::endl << "      ";
      file << (*(*iter))[0] << ", " << (*(*iter))[1];
    }
  file << "};" << std::endl;
}

void QuadraticPolygon::closeMe() const
{
  if(!front()->changeStartNodeWith(back()->getEndNode()))
    throw(Exception("big error: not closed polygon..."));
}

void QuadraticPolygon::circularPermute()
{
  if(_sub_edges.size()>1)
    {
      ElementaryEdge *first=_sub_edges.front();
      _sub_edges.pop_front();
      _sub_edges.push_back(first);
    }
}

bool QuadraticPolygon::isButterfly() const
{
  for(std::list<ElementaryEdge *>::const_iterator it=_sub_edges.begin();it!=_sub_edges.end();it++)
    {
      Edge *e1=(*it)->getPtr();
      std::list<ElementaryEdge *>::const_iterator it2=it;
      it2++;
      for(;it2!=_sub_edges.end();it2++)
        {
          MergePoints commonNode;
          ComposedEdge *outVal1=new ComposedEdge;
          ComposedEdge *outVal2=new ComposedEdge;
          Edge *e2=(*it2)->getPtr();
          if(e1->intersectWith(e2,commonNode,*outVal1,*outVal2))
            {
              Delete(outVal1);
              Delete(outVal2);
              return true;
            }
          Delete(outVal1);
          Delete(outVal2);
        }
    }
  return false;
}

void QuadraticPolygon::dumpInXfigFileWithOther(const ComposedEdge& other, const char *fileName) const
{
  std::ofstream file(fileName);
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
  std::ofstream file(fileName);
  const int resolution=1200;
  Bounds box;
  box.prepareForAggregation();
  fillBounds(box);
  dumpInXfigFile(file,resolution,box);
}

void QuadraticPolygon::dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const
{
  stream << "#FIG 3.2  Produced by xfig version 3.2.5-alpha5" << std::endl;
  stream << "Landscape" << std::endl;
  stream << "Center" << std::endl;
  stream << "Metric" << std::endl;
  stream << "Letter" << std::endl;
  stream << "100.00" << std::endl;
  stream << "Single" << std::endl;
  stream << "-2" << std::endl;
  stream << resolution << " 2" << std::endl;
  ComposedEdge::dumpInXfigFile(stream,resolution,box);
}

/*!
 * Warning contrary to intersectWith method this method is \b NOT const. 'this' and 'other' are modified after call of this method.
 */
double QuadraticPolygon::intersectWithAbs(QuadraticPolygon& other)
{
  double ret=0.,xBaryBB,yBaryBB;
  double fact=normalize(&other,xBaryBB,yBaryBB);
  std::vector<QuadraticPolygon *> polygs=intersectMySelfWith(other);
  for(std::vector<QuadraticPolygon *>::iterator iter=polygs.begin();iter!=polygs.end();iter++)
    {
      ret+=fabs((*iter)->getArea());
      delete *iter;
    }
  return ret*fact*fact;
}

/*!
 * Warning This method is \b NOT const. 'this' and 'other' are modified after call of this method.
 */
double QuadraticPolygon::intersectWithAbs1D(QuadraticPolygon& other, bool& isColinear)
{
  double ret = 0., xBaryBB, yBaryBB;
  double fact = normalize(&other, xBaryBB, yBaryBB);

  QuadraticPolygon cpyOfThis(*this);
  QuadraticPolygon cpyOfOther(other);
  int nbOfSplits = 0;
  splitPolygonsEachOther(cpyOfThis, cpyOfOther, nbOfSplits);
  //At this point cpyOfThis and cpyOfOther have been splited at maximum edge so that in/out can been done.
  performLocatingOperation(cpyOfOther);

  std::list<QuadraticPolygon *> cpyOfOtherZip = cpyOfOther.zipConsecutiveInSegments();

  isColinear = false;
  if (cpyOfOtherZip.size() == 1)
    {
      QuadraticPolygon& poly = *cpyOfOtherZip.front();
      if (poly.size() == 2 && poly.front()->isEqual(*poly.back()))
        isColinear = true;
    }

  for (std::list<QuadraticPolygon *>::iterator iter = cpyOfOtherZip.begin(); iter != cpyOfOtherZip.end(); iter++)
    {
      ret += fabs((*iter)->getPerimeter());
      delete *iter;
    }
  return ret * fact / 2.;
}

/*!
 * Warning contrary to intersectWith method this method is \b NOT const. 'this' and 'other' are modified after call of this method.
 */
double QuadraticPolygon::intersectWithAbs(QuadraticPolygon& other, double* barycenter)
{
  double ret=0.,bary[2],area,xBaryBB,yBaryBB;
  barycenter[0] = barycenter[1] = 0.;
  double fact=normalize(&other,xBaryBB,yBaryBB);
  std::vector<QuadraticPolygon *> polygs=intersectMySelfWith(other);
  for(std::vector<QuadraticPolygon *>::iterator iter=polygs.begin();iter!=polygs.end();iter++)
    {
      area=fabs((*iter)->getArea());
      (*iter)->getBarycenter(bary);
      delete *iter;
      ret+=area;
      barycenter[0] += bary[0]*area;
      barycenter[1] += bary[1]*area;
    }
  if ( ret > std::numeric_limits<double>::min() )
    {
      barycenter[0]=barycenter[0]/ret*fact+xBaryBB;
      barycenter[1]=barycenter[1]/ret*fact+yBaryBB;
      
    }
  return ret*fact*fact;
}

/*!
 * \b WARNING this method is const and other is const too. \b BUT location of Edges in 'this' and 'other' are nevertheless modified.
 * This is possible because loc attribute in Edge class is mutable.
 * This implies that if 'this' or/and 'other' are reused for intersect* method initLocations has to be called on each of this/them.
 */
double QuadraticPolygon::intersectWith(const QuadraticPolygon& other) const
{
  double ret=0.;
  std::vector<QuadraticPolygon *> polygs=intersectMySelfWith(other);
  for(std::vector<QuadraticPolygon *>::iterator iter=polygs.begin();iter!=polygs.end();iter++)
    {
      ret+=fabs((*iter)->getArea());
      delete *iter;
    }
  return ret;
}

/*!
 * \b WARNING this method is const and other is const too. \b BUT location of Edges in 'this' and 'other' are nevertheless modified.
 * This is possible because loc attribute in Edge class is mutable.
 * This implies that if 'this' or/and 'other' are reused for intersect* method initLocations has to be called on each of this/them.
 */
double QuadraticPolygon::intersectWith(const QuadraticPolygon& other, double* barycenter) const
{
  double ret=0., bary[2];
  barycenter[0] = barycenter[1] = 0.;
  std::vector<QuadraticPolygon *> polygs=intersectMySelfWith(other);
  for(std::vector<QuadraticPolygon *>::iterator iter=polygs.begin();iter!=polygs.end();iter++)
    {
      double area = fabs((*iter)->getArea());
      (*iter)->getBarycenter(bary);
      delete *iter;
      ret+=area;
      barycenter[0] += bary[0]*area;
      barycenter[1] += bary[1]*area;
    }
  if ( ret > std::numeric_limits<double>::min() )
    {
      barycenter[0] /= ret;
      barycenter[1] /= ret;
    }
  return ret;
}

/*!
 * \b WARNING this method is const and other is const too. \b BUT location of Edges in 'this' and 'other' are nevertheless modified.
 * This is possible because loc attribute in Edge class is mutable.
 * This implies that if 'this' or/and 'other' are reused for intersect* method initLocations has to be called on each of this/them.
 */
void QuadraticPolygon::intersectForPerimeter(const QuadraticPolygon& other, double& perimeterThisPart, double& perimeterOtherPart, double& perimeterCommonPart) const
{
  perimeterThisPart=0.; perimeterOtherPart=0.; perimeterCommonPart=0.;
  QuadraticPolygon cpyOfThis(*this);
  QuadraticPolygon cpyOfOther(other); int nbOfSplits=0;
  splitPolygonsEachOther(cpyOfThis,cpyOfOther,nbOfSplits);
  performLocatingOperation(cpyOfOther);
  other.performLocatingOperation(cpyOfThis);
  cpyOfThis.dispatchPerimeterExcl(perimeterThisPart,perimeterCommonPart);
  cpyOfOther.dispatchPerimeterExcl(perimeterOtherPart,perimeterCommonPart);
  perimeterCommonPart/=2.;
}

/*!
 * \b WARNING this method is const and other is const too. \b BUT location of Edges in 'this' and 'other' are nevertheless modified.
 * This is possible because loc attribute in Edge class is mutable.
 * This implies that if 'this' or/and 'other' are reused for intersect* method initLocations has to be called on each of this/them.
 *
 * polThis.size()==this->size() and polOther.size()==other.size().
 * For each ElementaryEdge of 'this', the corresponding contribution in resulting polygon is in 'polThis'.
 * For each ElementaryEdge of 'other', the corresponding contribution in resulting polygon is in 'polOther'.
 * As consequence common part are counted twice (in polThis \b and in polOther).
 */
void QuadraticPolygon::intersectForPerimeterAdvanced(const QuadraticPolygon& other, std::vector< double >& polThis, std::vector< double >& polOther) const
{
  polThis.resize(size());
  polOther.resize(other.size());
  IteratorOnComposedEdge it1((QuadraticPolygon *)this);
  int edgeId=0;
  for(it1.first();!it1.finished();it1.next(),edgeId++)
    {
      ElementaryEdge* curE1=it1.current();
      QuadraticPolygon cpyOfOther(other);
      QuadraticPolygon tmp;
      tmp.pushBack(curE1->clone());
      int tmp2;
      splitPolygonsEachOther(tmp,cpyOfOther,tmp2);
      other.performLocatingOperation(tmp);
      tmp.dispatchPerimeter(polThis[edgeId]);
    }
  //
  IteratorOnComposedEdge it2((QuadraticPolygon *)&other);
  edgeId=0;
  for(it2.first();!it2.finished();it2.next(),edgeId++)
    {
      ElementaryEdge* curE2=it2.current();
      QuadraticPolygon cpyOfThis(*this);
      QuadraticPolygon tmp;
      tmp.pushBack(curE2->clone());
      int tmp2;
      splitPolygonsEachOther(tmp,cpyOfThis,tmp2);
      performLocatingOperation(tmp);
      tmp.dispatchPerimeter(polOther[edgeId]);
    }
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
  IteratorOnComposedEdge it1((QuadraticPolygon *)this);
  int edgeId=0;
  for(it1.first();!it1.finished();it1.next(),edgeId++)
    {
      ElementaryEdge* curE1=it1.current();
      QuadraticPolygon cpyOfOther(other);
      QuadraticPolygon tmp;
      tmp.pushBack(curE1->clone());
      int tmp2;
      splitPolygonsEachOther(tmp,cpyOfOther,tmp2);
      numberOfCreatedPointsPerEdge[edgeId]=tmp.recursiveSize()-1;
    }
}

/*!
 * \b WARNING this method is const and other is const too. \b BUT location of Edges in 'this' and 'other' are nevertheless modified.
 * This is possible because loc attribute in Edge class is mutable.
 * This implies that if 'this' or/and 'other' are reused for intersect* method initLocations has to be called on each of this/them.
 */
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
void QuadraticPolygon::splitPolygonsEachOther(QuadraticPolygon& pol1, QuadraticPolygon& pol2, int& nbOfSplits)
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

/*!
 * Given 2 polygons 'pol1' and 'pol2' (localized) the resulting polygons are returned.
 *
 * this : pol2 simplified.
 * @param pol1 pol1 split.
 * @param pol2 pol2 split.
 */
std::vector<QuadraticPolygon *> QuadraticPolygon::buildIntersectionPolygons(const QuadraticPolygon& pol1, const QuadraticPolygon& pol2) const
{
  std::vector<QuadraticPolygon *> ret;
  std::list<QuadraticPolygon *> pol2Zip=pol2.zipConsecutiveInSegments();
  if(!pol2Zip.empty())
    closePolygons(pol2Zip,pol1,ret);
  else
    {//borders of pol2 do not cross pol1,and pol2 borders are outside of pol1. That is to say, either pol2 and pol1
      //do not overlap or  pol1 is fully inside pol2. So in the first case no intersection, in the other case
      //the intersection is pol1.
      ElementaryEdge *e1FromPol1=pol1[0];
      TypeOfEdgeLocInPolygon loc=FULL_ON_1;
      loc=e1FromPol1->locateFullyMySelf(*this,loc);
      if(loc==FULL_IN_1)
        ret.push_back(new QuadraticPolygon(pol1));
    }
  return ret;
}

/*!
 * Returns parts of potentially non closed-polygons. Each returned polygons are not mergeable.
 * this : pol2 split and locallized.
 */
std::list<QuadraticPolygon *> QuadraticPolygon::zipConsecutiveInSegments() const
{
  std::list<QuadraticPolygon *> ret;
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
 * @param pol2zip is a list of set of edges (openned polygon) coming from split polygon 2.
 * @param pol1 is split pol1.
 * @param results the resulting \b CLOSED polygons.
 */
void QuadraticPolygon::closePolygons(std::list<QuadraticPolygon *>& pol2Zip, const QuadraticPolygon& pol1,
                                     std::vector<QuadraticPolygon *>& results) const
{
  bool directionKnownInPol1=false;
  bool directionInPol1;
  for(std::list<QuadraticPolygon *>::iterator iter=pol2Zip.begin();iter!=pol2Zip.end();)
    {
      if((*iter)->completed())
        {
          results.push_back(*iter);
          directionKnownInPol1=false;
          iter=pol2Zip.erase(iter);
          continue;
        }
      if(!directionKnownInPol1)
        {
          if(!(*iter)->amIAChanceToBeCompletedBy(pol1,*this,directionInPol1))
            { delete *iter; iter=pol2Zip.erase(iter); continue; }
          else
            directionKnownInPol1=true;
        }
      std::list<QuadraticPolygon *>::iterator iter2=iter; iter2++;
      std::list<QuadraticPolygon *>::iterator iter3=(*iter)->fillAsMuchAsPossibleWith(pol1,iter2,pol2Zip.end(),directionInPol1);
      if(iter3!=pol2Zip.end())
        {
          (*iter)->pushBack(*iter3);
          SoftDelete(*iter3);
          pol2Zip.erase(iter3);
        }
    }
}

/*!
 * 'this' is expected to be set of edges (not closed) of pol2 split.
 */
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
  ElementaryEdge *e=_sub_edges.back();
  if(e->getLoc()==FULL_ON_1)
    {
      if(e->getPtr()==cur->getPtr())
        {
          direction=false;
          it.previousLoop();
          cur=it.current();
          Node *repr=cur->getPtr()->buildRepresentantOfMySelf();
          bool ret=pol2NotSplitted.isInOrOut(repr);
          repr->decrRef();
          return ret;
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

/*!
 * This method fills as much as possible 'this' (part of pol2 split) with edges of 'pol1Splitted'.
 */
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
  for(std::list<QuadraticPolygon *>::iterator iter=iStart;iter!=iEnd;iter++)
    if((*iter)->isNodeIn(n))
      return iter;
  return iEnd;
}
