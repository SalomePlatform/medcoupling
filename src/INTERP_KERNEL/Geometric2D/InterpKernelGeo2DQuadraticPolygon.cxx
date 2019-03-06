// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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

#include "InterpKernelGeo2DQuadraticPolygon.hxx"
#include "InterpKernelGeo2DElementaryEdge.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DAbstractEdge.hxx"
#include "InterpKernelGeo2DEdgeLin.hxx"
#include "InterpKernelGeo2DBounds.hxx"
#include "InterpKernelGeo2DEdge.txx"

#include "NormalizedUnstructuredMesh.hxx"

#include <fstream>
#include <sstream>
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
          Edge *newEdge=Edge::BuildFromXfigLine(stream);
          if(!empty())
            newEdge->changeStartNodeWith(back()->getEndNode());
          pushBack(newEdge);
        }
      while(1);
  }
  catch(const std::ifstream::failure&)
  {
  }
  catch(const std::exception & ex)
  {
      // Some code before this catch throws the C++98 version of the exception (mangled
      // name is " NSt8ios_base7failureE"), but FED24 compilation of the current version of the code
      // tries to catch the C++11 version of it (mangled name "NSt8ios_base7failureB5cxx11E").
      // So we have this nasty hack to catch both versions ...

      // TODO: the below should be replaced by a better handling avoiding exception throwing.
      if (std::string(ex.what()) == "basic_ios::clear")
        {
          //std::cout << "std::ios_base::failure C++11\n";
        }
      else
        throw ex;
  }
  front()->changeStartNodeWith(back()->getEndNode());
}

QuadraticPolygon::~QuadraticPolygon()
{
}

QuadraticPolygon *QuadraticPolygon::BuildLinearPolygon(std::vector<Node *>& nodes)
{
  QuadraticPolygon *ret(new QuadraticPolygon);
  std::size_t size=nodes.size();
  for(std::size_t i=0;i<size;i++)
    {
      ret->pushBack(new EdgeLin(nodes[i],nodes[(i+1)%size]));
      nodes[i]->decrRef();
    }
  return ret;
}

QuadraticPolygon *QuadraticPolygon::BuildArcCirclePolygon(std::vector<Node *>& nodes)
{
  QuadraticPolygon *ret(new QuadraticPolygon);
  std::size_t size=nodes.size();
  for(std::size_t i=0;i<size/2;i++)

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

Edge *QuadraticPolygon::BuildLinearEdge(std::vector<Node *>& nodes)
{
  if(nodes.size()!=2)
    throw INTERP_KERNEL::Exception("QuadraticPolygon::BuildLinearEdge : input vector is expected to be of size 2 !");
  Edge *ret(new EdgeLin(nodes[0],nodes[1]));
  nodes[0]->decrRef(); nodes[1]->decrRef();
  return ret;
}

Edge *QuadraticPolygon::BuildArcCircleEdge(std::vector<Node *>& nodes)
{
  if(nodes.size()!=3)
    throw INTERP_KERNEL::Exception("QuadraticPolygon::BuildArcCircleEdge : input vector is expected to be of size 3 !");
  EdgeLin *e1(new EdgeLin(nodes[0],nodes[2])),*e2(new EdgeLin(nodes[2],nodes[1]));
  SegSegIntersector inters(*e1,*e2);
  bool colinearity=inters.areColinears();
  delete e1; delete e2;
  Edge *ret(0);
  if(colinearity)
    ret=new EdgeLin(nodes[0],nodes[1]);
  else
    ret=new EdgeArcCircle(nodes[0],nodes[2],nodes[1]);
  nodes[0]->decrRef(); nodes[1]->decrRef(); nodes[2]->decrRef();
  return ret;
}

void QuadraticPolygon::BuildDbgFile(const std::vector<Node *>& nodes, const char *fileName)
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

bool QuadraticPolygon::isButterflyAbs()
{
  INTERP_KERNEL::Bounds b;
  double xBary,yBary;
  b.prepareForAggregation();
  fillBounds(b); 
  double dimChar=b.getCaracteristicDim();
  b.getBarycenter(xBary,yBary);
  applyGlobalSimilarity(xBary,yBary,dimChar);
  //
  return isButterfly();
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
 * This method splits 'this' with 'other' into smaller pieces localizable. 'mapThis' is a map that gives the correspondence
 * between nodes contained in 'this' and node ids in a global mesh.
 * In the same way, 'mapOther' gives the correspondence between nodes contained in 'other' and node ids in a
 * global mesh from which 'other' is extracted.
 * This method has 1 out parameter : 'edgesThis', After the call of this method, it contains the nodal connectivity (including type)
 * of 'this' into globlal "this mesh".
 * This method has 2 in/out parameters : 'subDivOther' and 'addCoo'.'otherEdgeIds' is useful to put values in
 * 'edgesThis', 'subDivOther' and 'addCoo'.
 * Size of 'otherEdgeIds' has to be equal to number of ElementaryEdges in 'other'. No check of that will be done.
 * The term 'abs' in the name recalls that we normalize the mesh (spatially) so that node coordinates fit into [0;1].
 * @param offset1 is the number of nodes contained in global mesh from which 'this' is extracted.
 * @param offset2 is the sum of nodes contained in global mesh from which 'this' is extracted and 'other' is extracted.
 * @param edgesInOtherColinearWithThis will be appended at the end of the vector with colinear edge ids of other (if any)
 * @param otherEdgeIds is a vector with the same size than other before calling this method. It gives in the same order
 * the cell id in global other mesh.
 */
void QuadraticPolygon::splitAbs(QuadraticPolygon& other,
                                const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther,
                                int offset1, int offset2 ,
                                const std::vector<int>& otherEdgeIds,
                                std::vector<int>& edgesThis, int cellIdThis,
                                std::vector< std::vector<int> >& edgesInOtherColinearWithThis, std::vector< std::vector<int> >& subDivOther,
                                std::vector<double>& addCoo, std::map<int,int>& mergedNodes)
{
  double xBaryBB, yBaryBB;
  double fact=normalizeExt(&other, xBaryBB, yBaryBB);

  //
  IteratorOnComposedEdge itThis(this),itOther(&other);  // other is (part of) the tool mesh
  MergePoints merge;
  ComposedEdge *cThis=new ComposedEdge;
  ComposedEdge *cOther=new ComposedEdge;
  int i=0;
  std::map<INTERP_KERNEL::Node *,int> mapAddCoo;
  for(itOther.first();!itOther.finished();itOther.next(),i++)
    {
      // For each edge of 'other', proceed with intersections: the edge might split into sub-edges, 'otherTmp' will hold the final split result.
      // In the process of going through all edges of 'other', 'this' (which contains initially only one edge)
      // is sub-divised into several edges : each of them has to be tested when intersecting the next candidate stored in 'other'.
      QuadraticPolygon otherTmp;
      ElementaryEdge* curOther=itOther.current();
      otherTmp.pushBack(new ElementaryEdge(curOther->getPtr(),curOther->getDirection())); curOther->getPtr()->incrRef();
      IteratorOnComposedEdge itOtherTmp(&otherTmp);
      for(itOtherTmp.first();!itOtherTmp.finished();itOtherTmp.next())
        {
          ElementaryEdge* curOtherTmp=itOtherTmp.current();
          if(!curOtherTmp->isThereStartPoint())
            itThis.first();   // reset iterator on 'this'
          else
            itThis=curOtherTmp->getIterator();
          for(;!itThis.finished();)
            {
              ElementaryEdge* curThis=itThis.current();
              merge.clear();
              //
              std::map<INTERP_KERNEL::Node *,int>::const_iterator thisStart(mapThis.find(curThis->getStartNode())),thisEnd(mapThis.find(curThis->getEndNode())),
                                                                  otherStart(mapOther.find(curOtherTmp->getStartNode())),otherEnd(mapOther.find(curOtherTmp->getEndNode()));
              int thisStart2(thisStart==mapThis.end()?-1:(*thisStart).second), thisEnd2(thisEnd==mapThis.end()?-1:(*thisEnd).second),
                  otherStart2(otherStart==mapOther.end()?-1:(*otherStart).second+offset1),otherEnd2(otherEnd==mapOther.end()?-1:(*otherEnd).second+offset1);
              //
              if(curThis->getPtr()->intersectWith(curOtherTmp->getPtr(),merge,*cThis,*cOther))
                {
                  if(!curThis->getDirection()) cThis->reverse();
                  if(!curOtherTmp->getDirection()) cOther->reverse();
                  // Substitution of a single simple edge by two sub-edges resulting from the intersection
                  //    First modify the edges currently pointed by itThis and itOtherTmp so that the newly created node
                  //    becomes the end of the previous sub-edge and the beginning of the next one.
                  UpdateNeighbours(merge,itThis,itOtherTmp,cThis,cOther);
                  delete curThis;           // <-- destroying simple edge coming from pol1
                  delete curOtherTmp;       // <-- destroying simple edge coming from pol2
                  //     Then insert second part of the intersection.
                  itThis.insertElemEdges(cThis,true);       // <-- 2nd param is true to go next.
                  itOtherTmp.insertElemEdges(cOther,false); // <-- 2nd param is false to avoid to go next.
                  curOtherTmp=itOtherTmp.current();
                  //
                  itThis.assignMySelfToAllElems(cOther);
                  SoftDelete(cThis);
                  SoftDelete(cOther);
                  cThis=new ComposedEdge;
                  cOther=new ComposedEdge;
                }
              else
                {
                  UpdateNeighbours(merge,itThis,itOtherTmp,curThis,curOtherTmp);
                  itThis.next();
                }
              merge.updateMergedNodes(thisStart2,thisEnd2,otherStart2,otherEnd2,mergedNodes);
            }
        }
      // If one sub-edge of otherTmp is "ON" an edge of this, then we have colinearity (all edges in otherTmp are //)
      if(otherTmp.presenceOfOn())
        edgesInOtherColinearWithThis[otherEdgeIds[i]].push_back(cellIdThis);
      // Converting back to integer connectivity:
      if(otherTmp._sub_edges.size()>1)   // only if a new point has been added (i.e. an actual intersection was done)
        {
          int jj = 0, sz(otherTmp._sub_edges.size());
          for(std::list<ElementaryEdge *>::const_iterator it=otherTmp._sub_edges.begin();it!=otherTmp._sub_edges.end();it++, jj++)
            {
              short skipStartOrEnd = jj == 0 ? -1 : (jj == sz-1 ? 1 : 0);  // -1 means START, 1 means END, 0 other
              (*it)->fillGlobalInfoAbs2(mapThis,mapOther,offset1,offset2,
                                      fact,xBaryBB,yBaryBB, skipStartOrEnd,
                                      /*out*/ subDivOther[otherEdgeIds[i]],addCoo,mapAddCoo);
            }
        }
    }
  Delete(cThis);
  Delete(cOther);
  //
  for(std::list<ElementaryEdge *>::const_iterator it=_sub_edges.begin();it!=_sub_edges.end();it++)
    (*it)->fillGlobalInfoAbs(mapThis,mapOther,offset1,offset2,/**/fact,xBaryBB,yBaryBB,/**/edgesThis,addCoo,mapAddCoo);
  //
}

/*!
 * This method builds 'this' from its descending conn stored in crude mode (MEDCoupling).
 * Descending conn is in FORTRAN relative mode in order to give the
 * orientation of edge (see buildDescendingConnectivity2() method).
 * See appendEdgeFromCrudeDataArray() for params description.
 */
void QuadraticPolygon::buildFromCrudeDataArray(const std::map<int,INTERP_KERNEL::Node *>& mapp, bool isQuad, const int *nodalBg, const double *coords,
                                               const int *descBg, const int *descEnd, const std::vector<std::vector<int> >& intersectEdges)
{
  std::size_t nbOfSeg=std::distance(descBg,descEnd);
  for(std::size_t i=0;i<nbOfSeg;i++)
    {
      appendEdgeFromCrudeDataArray(i,mapp,isQuad,nodalBg,coords,descBg,descEnd,intersectEdges);
    }
}

void QuadraticPolygon::appendEdgeFromCrudeDataArray(std::size_t edgePos, const std::map<int,INTERP_KERNEL::Node *>& mapp, bool isQuad,
                                                    const int *nodalBg, const double *coords,
                                                    const int *descBg, const int *descEnd, const std::vector<std::vector<int> >& intersectEdges)
{
  if(!isQuad)
    {
      bool direct=descBg[edgePos]>0;
      int edgeId=abs(descBg[edgePos])-1; // back to C indexing mode
      const std::vector<int>& subEdge=intersectEdges[edgeId];
      std::size_t nbOfSubEdges=subEdge.size()/2;
      for(std::size_t j=0;j<nbOfSubEdges;j++)
        appendSubEdgeFromCrudeDataArray(0,j,direct,edgeId,subEdge,mapp);
    }
  else
    {
      std::size_t nbOfSeg=std::distance(descBg,descEnd);
      const double *st=coords+2*(nodalBg[edgePos]); 
      INTERP_KERNEL::Node *st0=new INTERP_KERNEL::Node(st[0],st[1]);
      const double *endd=coords+2*(nodalBg[(edgePos+1)%nbOfSeg]);
      INTERP_KERNEL::Node *endd0=new INTERP_KERNEL::Node(endd[0],endd[1]);
      const double *middle=coords+2*(nodalBg[edgePos+nbOfSeg]);
      INTERP_KERNEL::Node *middle0=new INTERP_KERNEL::Node(middle[0],middle[1]);
      EdgeLin *e1,*e2;
      e1=new EdgeLin(st0,middle0);
      e2=new EdgeLin(middle0,endd0);
      SegSegIntersector inters(*e1,*e2);
      bool colinearity=inters.areColinears();
      delete e1; delete e2;
      //
      bool direct=descBg[edgePos]>0;
      int edgeId=abs(descBg[edgePos])-1;
      const std::vector<int>& subEdge=intersectEdges[edgeId];
      std::size_t nbOfSubEdges=subEdge.size()/2;
      if(colinearity)
        {   
          for(std::size_t j=0;j<nbOfSubEdges;j++)
            appendSubEdgeFromCrudeDataArray(0,j,direct,edgeId,subEdge,mapp);
        }
      else
        {
          Edge *e=new EdgeArcCircle(st0,middle0,endd0,true);
          for(std::size_t j=0;j<nbOfSubEdges;j++)
            appendSubEdgeFromCrudeDataArray(e,j,direct,edgeId,subEdge,mapp);
          e->decrRef();
        }
      st0->decrRef(); endd0->decrRef(); middle0->decrRef();
    }
}

void QuadraticPolygon::appendSubEdgeFromCrudeDataArray(Edge *baseEdge, std::size_t j, bool direct, int edgeId, const std::vector<int>& subEdge, const std::map<int,INTERP_KERNEL::Node *>& mapp)
{
  std::size_t nbOfSubEdges=subEdge.size()/2;
  if(!baseEdge)
    {//it is not a quadratic subedge
      Node *start=(*mapp.find(direct?subEdge[2*j]:subEdge[2*nbOfSubEdges-2*j-1])).second;
      Node *end=(*mapp.find(direct?subEdge[2*j+1]:subEdge[2*nbOfSubEdges-2*j-2])).second;
      ElementaryEdge *e=ElementaryEdge::BuildEdgeFromStartEndDir(true,start,end);
      pushBack(e);
    }
  else
    {//it is a quadratic subedge
      Node *start=(*mapp.find(direct?subEdge[2*j]:subEdge[2*nbOfSubEdges-2*j-1])).second;
      Node *end=(*mapp.find(direct?subEdge[2*j+1]:subEdge[2*nbOfSubEdges-2*j-2])).second;
      Edge *ee=baseEdge->buildEdgeLyingOnMe(start,end);
      ElementaryEdge *eee=new ElementaryEdge(ee,true);
      pushBack(eee);
    }
}

/*!
 * This method builds from descending conn of a quadratic polygon stored in crude mode (MEDCoupling). Descending conn is in FORTRAN relative mode in order to give the
 * orientation of edge.
 */
void QuadraticPolygon::buildFromCrudeDataArray2(const std::map<int,INTERP_KERNEL::Node *>& mapp, bool isQuad, const int *nodalBg, const double *coords, const int *descBg, const int *descEnd, const std::vector<std::vector<int> >& intersectEdges2,
                                                const INTERP_KERNEL::QuadraticPolygon& pol1, const int *descBg1, const int *descEnd1, const std::vector<std::vector<int> >& intersectEdges1,
                                                const std::vector< std::vector<int> >& colinear1,
                                                std::map<int,std::vector<INTERP_KERNEL::ElementaryEdge *> >& alreadyExistingIn2)
{
  std::size_t nbOfSeg=std::distance(descBg,descEnd);
  for(std::size_t i=0;i<nbOfSeg;i++)//loop over all edges of pol2
    {
      bool direct=descBg[i]>0;
      int edgeId=abs(descBg[i])-1;//current edge id of pol2
      std::map<int,std::vector<INTERP_KERNEL::ElementaryEdge *> >::const_iterator it1=alreadyExistingIn2.find(descBg[i]),it2=alreadyExistingIn2.find(-descBg[i]);
      if(it1!=alreadyExistingIn2.end() || it2!=alreadyExistingIn2.end())
        {
          bool sameDir=(it1!=alreadyExistingIn2.end());
          const std::vector<INTERP_KERNEL::ElementaryEdge *>& edgesAlreadyBuilt=sameDir?(*it1).second:(*it2).second;
          if(sameDir)
            {
              for(std::vector<INTERP_KERNEL::ElementaryEdge *>::const_iterator it3=edgesAlreadyBuilt.begin();it3!=edgesAlreadyBuilt.end();it3++)
                {
                  Edge *ee=(*it3)->getPtr(); ee->incrRef();
                  pushBack(new ElementaryEdge(ee,(*it3)->getDirection()));
                }
            }
          else
            {
              for(std::vector<INTERP_KERNEL::ElementaryEdge *>::const_reverse_iterator it4=edgesAlreadyBuilt.rbegin();it4!=edgesAlreadyBuilt.rend();it4++)
                {
                  Edge *ee=(*it4)->getPtr(); ee->incrRef();
                  pushBack(new ElementaryEdge(ee,!(*it4)->getDirection()));
                }
            }
          continue;
        }
      bool directos=colinear1[edgeId].empty();
      std::vector<std::pair<int,std::pair<bool,int> > > idIns1;
      int offset1=0;
      if(!directos)
        {// if the current edge of pol2 has one or more colinear edges part into pol1
          const std::vector<int>& c=colinear1[edgeId];
          std::size_t nbOfEdgesIn1=std::distance(descBg1,descEnd1);
          for(std::size_t j=0;j<nbOfEdgesIn1;j++)
            {
              int edgeId1=abs(descBg1[j])-1;
              if(std::find(c.begin(),c.end(),edgeId1)!=c.end())
                {
                  idIns1.push_back(std::pair<int,std::pair<bool,int> >(edgeId1,std::pair<bool,int>(descBg1[j]>0,offset1)));// it exists an edge into pol1 given by tuple (idIn1,direct1) that is colinear at edge 'edgeId' in pol2
                  //std::pair<edgeId1); direct1=descBg1[j]>0;
                }
              offset1+=intersectEdges1[edgeId1].size()/2;//offset1 is used to find the INTERP_KERNEL::Edge * instance into pol1 that will be part of edge into pol2
            }
          directos=idIns1.empty();
        }
      if(directos)
        {//no subpart of edge 'edgeId' of pol2 is in pol1 so let's operate the same thing that QuadraticPolygon::buildFromCrudeDataArray method
          std::size_t oldSz=_sub_edges.size();
          appendEdgeFromCrudeDataArray(i,mapp,isQuad,nodalBg,coords,descBg,descEnd,intersectEdges2);
          std::size_t newSz=_sub_edges.size();
          std::size_t zeSz=newSz-oldSz;
          alreadyExistingIn2[descBg[i]].resize(zeSz);
          std::list<ElementaryEdge *>::const_reverse_iterator it5=_sub_edges.rbegin();
          for(std::size_t p=0;p<zeSz;p++,it5++)
            alreadyExistingIn2[descBg[i]][zeSz-p-1]=*it5;
        }
      else
        {//there is subpart of edge 'edgeId' of pol2 inside pol1
          const std::vector<int>& subEdge=intersectEdges2[edgeId];
          std::size_t nbOfSubEdges=subEdge.size()/2;
          for(std::size_t j=0;j<nbOfSubEdges;j++)
            {
              int idBg=direct?subEdge[2*j]:subEdge[2*nbOfSubEdges-2*j-1];
              int idEnd=direct?subEdge[2*j+1]:subEdge[2*nbOfSubEdges-2*j-2];
              bool direction11,found=false;
              bool direct1;//store if needed the direction in 1
              int offset2;
              std::size_t nbOfSubEdges1;
              for(std::vector<std::pair<int,std::pair<bool,int> > >::const_iterator it=idIns1.begin();it!=idIns1.end() && !found;it++)
                {
                  int idIn1=(*it).first;//store if needed the cell id in 1
                  direct1=(*it).second.first;
                  offset1=(*it).second.second;
                  const std::vector<int>& subEdge1PossiblyAlreadyIn1=intersectEdges1[idIn1];
                  nbOfSubEdges1=subEdge1PossiblyAlreadyIn1.size()/2;
                  offset2=0;
                  for(std::size_t k=0;k<nbOfSubEdges1 && !found;k++)
                    {//perform a loop on all subedges of pol1 that includes edge 'edgeId' of pol2. For the moment we iterate only on subedges of ['idIn1']... To improve
                      if(subEdge1PossiblyAlreadyIn1[2*k]==idBg && subEdge1PossiblyAlreadyIn1[2*k+1]==idEnd)
                        { direction11=true; found=true; }
                      else if(subEdge1PossiblyAlreadyIn1[2*k]==idEnd && subEdge1PossiblyAlreadyIn1[2*k+1]==idBg)
                        { direction11=false; found=true; }
                      else
                        offset2++;
                    }
                }
              if(!found)
                {//the current subedge of edge 'edgeId' of pol2 is not a part of the colinear edge 'idIn1' of pol1 -> build new Edge instance
                  //appendEdgeFromCrudeDataArray(j,mapp,isQuad,nodalBg,coords,descBg,descEnd,intersectEdges);
                  Node *start=(*mapp.find(idBg)).second;
                  Node *end=(*mapp.find(idEnd)).second;
                  ElementaryEdge *e=ElementaryEdge::BuildEdgeFromStartEndDir(true,start,end);
                  pushBack(e);
                  alreadyExistingIn2[descBg[i]].push_back(e);
                }
              else
                {//the current subedge of edge 'edgeId' of pol2 is part of the colinear edge 'idIn1' of pol1 -> reuse Edge instance of pol1
                  ElementaryEdge *e=pol1[offset1+(direct1?offset2:nbOfSubEdges1-offset2-1)];
                  Edge *ee=e->getPtr();
                  ee->incrRef();
                  ElementaryEdge *e2=new ElementaryEdge(ee,!(direct1^direction11));
                  pushBack(e2);
                  alreadyExistingIn2[descBg[i]].push_back(e2);
                }
            }
        }
    }
}

/*!
 * Method expected to be called on pol2. Every params not suffixed by numbered are supposed to refer to pol2 (this).
 * Method to find edges that are ON.
 */
void QuadraticPolygon::updateLocOfEdgeFromCrudeDataArray2(const int *descBg, const int *descEnd, const std::vector<std::vector<int> >& intersectEdges,
                                                          const INTERP_KERNEL::QuadraticPolygon& pol1, const int *descBg1, const int *descEnd1,
                                                          const std::vector<std::vector<int> >& intersectEdges1, const std::vector< std::vector<int> >& colinear1) const
{
  std::size_t nbOfSeg=std::distance(descBg,descEnd);
  for(std::size_t i=0;i<nbOfSeg;i++)//loop over all edges of pol2
    {
      bool direct=descBg[i]>0;
      int edgeId=abs(descBg[i])-1;//current edge id of pol2
      const std::vector<int>& c=colinear1[edgeId];
      if(c.empty())
        continue;
      const std::vector<int>& subEdge=intersectEdges[edgeId];
      std::size_t nbOfSubEdges=subEdge.size()/2;
      //
      std::size_t nbOfEdgesIn1=std::distance(descBg1,descEnd1);
      int offset1=0;
      for(std::size_t j=0;j<nbOfEdgesIn1;j++)
        {
          int edgeId1=abs(descBg1[j])-1;
          if(std::find(c.begin(),c.end(),edgeId1)!=c.end())
            {
              for(std::size_t k=0;k<nbOfSubEdges;k++)
                {
                  int idBg=direct?subEdge[2*k]:subEdge[2*nbOfSubEdges-2*k-1];
                  int idEnd=direct?subEdge[2*k+1]:subEdge[2*nbOfSubEdges-2*k-2];
                  int idIn1=edgeId1;
                  bool direct1=descBg1[j]>0;
                  const std::vector<int>& subEdge1PossiblyAlreadyIn1=intersectEdges1[idIn1];
                  std::size_t nbOfSubEdges1=subEdge1PossiblyAlreadyIn1.size()/2;
                  int offset2=0;
                  bool found=false;
                  for(std::size_t kk=0;kk<nbOfSubEdges1 && !found;kk++)
                    {
                      found=(subEdge1PossiblyAlreadyIn1[2*kk]==idBg && subEdge1PossiblyAlreadyIn1[2*kk+1]==idEnd) || (subEdge1PossiblyAlreadyIn1[2*kk]==idEnd && subEdge1PossiblyAlreadyIn1[2*kk+1]==idBg);
                      if(!found)
                        offset2++;
                    }
                  if(found)
                    {
                      ElementaryEdge *e=pol1[offset1+(direct1?offset2:nbOfSubEdges1-offset2-1)];
                      e->getPtr()->declareOn();
                    }
                }
            }
          offset1+=intersectEdges1[edgeId1].size()/2;//offset1 is used to find the INTERP_KERNEL::Edge * instance into pol1 that will be part of edge into pol2
        }
    }
}

void QuadraticPolygon::appendCrudeData(const std::map<INTERP_KERNEL::Node *,int>& mapp, double xBary, double yBary, double fact, int offset, std::vector<double>& addCoordsQuadratic, std::vector<int>& conn, std::vector<int>& connI) const
{
  int nbOfNodesInPg=0;
  bool presenceOfQuadratic=presenceOfQuadraticEdge();
  conn.push_back(presenceOfQuadratic?NORM_QPOLYG:NORM_POLYGON);
  for(std::list<ElementaryEdge *>::const_iterator it=_sub_edges.begin();it!=_sub_edges.end();it++)
    {
      Node *tmp=0;
      tmp=(*it)->getStartNode();
      std::map<INTERP_KERNEL::Node *,int>::const_iterator it1=mapp.find(tmp);
      conn.push_back((*it1).second);
      nbOfNodesInPg++;
    }
  if(presenceOfQuadratic)
    {
      int j=0;
      int off=offset+((int)addCoordsQuadratic.size())/2;
      for(std::list<ElementaryEdge *>::const_iterator it=_sub_edges.begin();it!=_sub_edges.end();it++,j++,nbOfNodesInPg++)
        {
          INTERP_KERNEL::Node *node=(*it)->getPtr()->buildRepresentantOfMySelf();
          node->unApplySimilarity(xBary,yBary,fact);
          addCoordsQuadratic.push_back((*node)[0]);
          addCoordsQuadratic.push_back((*node)[1]);
          conn.push_back(off+j);
          node->decrRef();
        }
    }
  connI.push_back(connI.back()+nbOfNodesInPg+1);
}

/*!
 * This method make the hypothesis that \a this and \a other are split at the minimum into edges that are fully IN, OUT or ON.
 * This method returns newly created polygons in \a conn and \a connI and the corresponding ids ( \a idThis, \a idOther) are stored respectively into \a nbThis and \a nbOther.
 * @param [in,out] edgesThis, parameter that keep informed the caller about the edges in this not shared by the result of intersection of \a this with \a other
 * @param [in,out] edgesBoundaryOther, parameter that stores all edges in result of intersection that are not
 */
void QuadraticPolygon::buildPartitionsAbs(QuadraticPolygon& other, std::set<INTERP_KERNEL::Edge *>& edgesThis, std::set<INTERP_KERNEL::Edge *>& edgesBoundaryOther,
                                          const std::map<INTERP_KERNEL::Node *,int>& mapp, int idThis, int idOther, int offset,
                                          std::vector<double>& addCoordsQuadratic, std::vector<int>& conn, std::vector<int>& connI,
                                          std::vector<int>& nbThis, std::vector<int>& nbOther)
{
  double xBaryBB, yBaryBB;
  double fact=normalizeExt(&other, xBaryBB, yBaryBB);
  //Locate \a this relative to \a other (edges of \a this, aka \a pol1 are marked as IN or OUT)
  other.performLocatingOperationSlow(*this);  // without any assumption
  std::vector<QuadraticPolygon *> res=buildIntersectionPolygons(*this,other);
  for(std::vector<QuadraticPolygon *>::iterator it=res.begin();it!=res.end();it++)
    {
      (*it)->appendCrudeData(mapp,xBaryBB,yBaryBB,fact,offset,addCoordsQuadratic,conn,connI);
      INTERP_KERNEL::IteratorOnComposedEdge it1(*it);
      for(it1.first();!it1.finished();it1.next())
        {
          Edge *e=it1.current()->getPtr();
          if(edgesThis.find(e)!=edgesThis.end())
            edgesThis.erase(e);
          else
            {
              if(edgesBoundaryOther.find(e)!=edgesBoundaryOther.end())
                edgesBoundaryOther.erase(e);
              else
                edgesBoundaryOther.insert(e);
            }
        }
      nbThis.push_back(idThis);
      nbOther.push_back(idOther);
      delete *it;
    }
  unApplyGlobalSimilarityExt(other,xBaryBB,yBaryBB,fact);
}

/*!
 * Remove the two 'identical' edges from the list, and drecRef the objects.
 */
void QuadraticPolygon::cleanDegeneratedConsecutiveEdges()
{
  IteratorOnComposedEdge it(this);
  ElementaryEdge * prevEdge = 0;
  if  (recursiveSize() > 2)
    for(it.first();!it.finished();it.next())
      {
        ElementaryEdge * cur = it.current();
        if (prevEdge && prevEdge->hasSameExtremities(*cur))
          {
            it.eraseCurrent();
            it.eraseCurrent();
            prevEdge = it.current();
          }
        else
            prevEdge = cur;
      }
}

/*!
 * Warning This method is \b NOT const. 'this' and 'other' are modified after call of this method.
 * 'other' is a QuadraticPolygon of \b non closed edges.
 */
double QuadraticPolygon::intersectWithAbs1D(QuadraticPolygon& other, bool& isColinear)
{
  double ret = 0., xBaryBB, yBaryBB;
  double fact = normalize(&other, xBaryBB, yBaryBB);

  QuadraticPolygon cpyOfThis(*this);
  QuadraticPolygon cpyOfOther(other);
  int nbOfSplits = 0;
  SplitPolygonsEachOther(cpyOfThis, cpyOfOther, nbOfSplits);
  //At this point cpyOfThis and cpyOfOther have been splited at maximum edge so that in/out can been done.
  performLocatingOperation(cpyOfOther);
  isColinear = false;
  for(std::list<ElementaryEdge *>::const_iterator it=cpyOfOther._sub_edges.begin();it!=cpyOfOther._sub_edges.end();it++)
    {
      switch((*it)->getLoc())
      {
        case FULL_IN_1:
          {
            ret += fabs((*it)->getPtr()->getCurveLength());
            break;
          }
        case FULL_ON_1:
          {
            isColinear=true;
            ret += fabs((*it)->getPtr()->getCurveLength());
            break;
          }
        default:
          {
          }
      }
    }
  return ret * fact;
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
  SplitPolygonsEachOther(cpyOfThis,cpyOfOther,nbOfSplits);
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
  IteratorOnComposedEdge it1(const_cast<QuadraticPolygon *>(this));
  int edgeId=0;
  for(it1.first();!it1.finished();it1.next(),edgeId++)
    {
      ElementaryEdge* curE1=it1.current();
      QuadraticPolygon cpyOfOther(other);
      QuadraticPolygon tmp;
      tmp.pushBack(curE1->clone());
      int tmp2;
      SplitPolygonsEachOther(tmp,cpyOfOther,tmp2);
      other.performLocatingOperation(tmp);
      tmp.dispatchPerimeter(polThis[edgeId]);
    }
  //
  IteratorOnComposedEdge it2(const_cast<QuadraticPolygon *>(&other));
  edgeId=0;
  for(it2.first();!it2.finished();it2.next(),edgeId++)
    {
      ElementaryEdge* curE2=it2.current();
      QuadraticPolygon cpyOfThis(*this);
      QuadraticPolygon tmp;
      tmp.pushBack(curE2->clone());
      int tmp2;
      SplitPolygonsEachOther(tmp,cpyOfThis,tmp2);
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
  IteratorOnComposedEdge it1(const_cast<QuadraticPolygon *>(this));
  int edgeId=0;
  for(it1.first();!it1.finished();it1.next(),edgeId++)
    {
      ElementaryEdge* curE1=it1.current();
      QuadraticPolygon cpyOfOther(other);
      QuadraticPolygon tmp;
      tmp.pushBack(curE1->clone());
      int tmp2;
      SplitPolygonsEachOther(tmp,cpyOfOther,tmp2);
      numberOfCreatedPointsPerEdge[edgeId]=tmp.recursiveSize()-1;
    }
}

/*!
 * \b WARNING this method is const and other is const too. \b BUT location of Edges in 'this' and 'other' are nevertheless modified.
 * This is possible because loc attribute in Edge class is mutable.
 * This implies that if 'this' or/and 'other' are reused for intersect* method initLocations has to be called on each of this/them.
 * This method is currently not used by any high level functionality.
 */
std::vector<QuadraticPolygon *> QuadraticPolygon::intersectMySelfWith(const QuadraticPolygon& other) const
{
  QuadraticPolygon cpyOfThis(*this);
  QuadraticPolygon cpyOfOther(other); int nbOfSplits=0;
  SplitPolygonsEachOther(cpyOfThis,cpyOfOther,nbOfSplits);
  //At this point cpyOfThis and cpyOfOther have been splited at maximum edge so that in/out can been done.
  performLocatingOperation(cpyOfOther);
  return other.buildIntersectionPolygons(cpyOfOther, cpyOfThis);
}

/*!
 * This method is typically the first step of boolean operations between pol1 and pol2.
 * This method perform the minimal splitting so that at the end each edges constituting pol1 are fully either IN or OUT or ON.
 * @param pol1 IN/OUT param that is equal to 'this' when called.
 */
void QuadraticPolygon::SplitPolygonsEachOther(QuadraticPolygon& pol1, QuadraticPolygon& pol2, int& nbOfSplits)
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
              UpdateNeighbours(merge,it1,it2,c1,c2);
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
              UpdateNeighbours(merge,it1,it2,curE1,curE2);
              it1.next();
            }
        }
    }
  Delete(c1);
  Delete(c2);
}

void QuadraticPolygon::performLocatingOperation(QuadraticPolygon& pol1) const
{
  IteratorOnComposedEdge it(&pol1);
  TypeOfEdgeLocInPolygon loc=FULL_ON_1;
  for(it.first();!it.finished();it.next())
    {
      ElementaryEdge *cur=it.current();
      loc=cur->locateFullyMySelf(*this,loc);//*this=pol2=other
    }
}

void QuadraticPolygon::performLocatingOperationSlow(QuadraticPolygon& pol2) const
{
  IteratorOnComposedEdge it(&pol2);
  for(it.first();!it.finished();it.next())
    {
      ElementaryEdge *cur=it.current();
      cur->locateFullyMySelfAbsolute(*this);
    }
}

/*!
 * Given 2 polygons \a pol1 and \a pol2 (localized) the resulting polygons are returned.
 *
 * this : pol1 simplified.
 * @param [in] pol1 pol1 split.
 * @param [in] pol2 pol2 split.
 */
std::vector<QuadraticPolygon *> QuadraticPolygon::buildIntersectionPolygons(const QuadraticPolygon& pol1, const QuadraticPolygon& pol2) const
{
  std::vector<QuadraticPolygon *> ret;
  // Extract from pol1, and pol1 only, all consecutive edges.
  // pol1Zip contains concatenated pieces of pol1 which are part of the resulting intersecting cell being built.
  std::list<QuadraticPolygon *> pol1Zip=pol1.zipConsecutiveInSegments();
  if(!pol1Zip.empty())
    ClosePolygons(pol1Zip,*this,pol2,ret);
  else
    {//borders of pol1 do not cross pol2,and pol1 borders are outside of pol2. That is to say, either pol1 and pol2
      //do not overlap or  pol2 is fully inside pol1. So in the first case no intersection, in the other case
      //the intersection is pol2.
      ElementaryEdge *e1FromPol2=pol2[0];
      TypeOfEdgeLocInPolygon loc=FULL_ON_1;
      loc=e1FromPol2->locateFullyMySelf(*this,loc);
      if(loc==FULL_IN_1)
        ret.push_back(new QuadraticPolygon(pol2));
    }
  return ret;
}

/*!
 * Returns parts of potentially non closed-polygons. Each returned polygons are not mergeable.
 * this : pol1 split and localized.
 */
std::list<QuadraticPolygon *> QuadraticPolygon::zipConsecutiveInSegments() const
{
  std::list<QuadraticPolygon *> ret;
  IteratorOnComposedEdge it(const_cast<QuadraticPolygon *>(this));
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
 * @param [in] pol1zip is a list of set of edges (=an opened polygon) coming from split polygon 1.
 * @param [in] pol1 should be considered as pol1Simplified.
 * @param [in] pol2 is split pol2.
 * @param [out] results the resulting \b CLOSED polygons.
 */
void QuadraticPolygon::ClosePolygons(std::list<QuadraticPolygon *>& pol1Zip, const QuadraticPolygon& pol1, const QuadraticPolygon& pol2,
                                     std::vector<QuadraticPolygon *>& results)
{
  bool directionKnownInPol2=false;
  bool directionInPol2;
  bool needCleaning = false;
  for(std::list<QuadraticPolygon *>::iterator iter=pol1Zip.begin();iter!=pol1Zip.end();)
    {
      // Build incrementally the full closed cells from the consecutive line parts already built in pol1Zip.
      // At the end of the process the item currently iterated has been totally completed (start_node=end_node)
      // This process can produce several cells.
      if((*iter)->completed())
        {
          if (needCleaning)
            (*iter)->cleanDegeneratedConsecutiveEdges();
          results.push_back(*iter);
          directionKnownInPol2=false;
          needCleaning=false;
          iter=pol1Zip.erase(iter);
          continue;
        }
      if(!directionKnownInPol2)
        {
          if(!(*iter)->haveIAChanceToBeCompletedBy(pol1,pol2,directionInPol2, needCleaning))
            { delete *iter; iter=pol1Zip.erase(iter); continue; }
          else
            directionKnownInPol2=true;
        }
      std::list<QuadraticPolygon *>::iterator iter2=iter; iter2++;
      // Fill as much as possible the current iterate (=a part of pol1) with consecutive pieces from pol2:
      std::list<QuadraticPolygon *>::iterator iter3=(*iter)->fillAsMuchAsPossibleWith(pol2,iter2,pol1Zip.end(),directionInPol2);
      // and now add a full connected piece from pol1Zip:
      if(iter3!=pol1Zip.end())
        {
          (*iter)->pushBack(*iter3);
          SoftDelete(*iter3);
          pol1Zip.erase(iter3);
        }
    }
}

/*!
 * 'this' is expected to be set of edges (not closed) of pol1 split.
 */
bool QuadraticPolygon::haveIAChanceToBeCompletedBy(const QuadraticPolygon& pol1NotSplitted, const QuadraticPolygon& pol2Splitted,
                                                   bool& direction, bool& needCleaning) const
{
  needCleaning = false;
  IteratorOnComposedEdge it2(const_cast<QuadraticPolygon *>(&pol2Splitted));
  bool found=false;
  Node *n=getEndNode();
  ElementaryEdge *cur=it2.current();
  // Find edge in pol2 whose start node is the end node of the current piece in pol1Zip (*this)
  for(it2.first();!it2.finished() && !found;)
    {
      cur=it2.current();
      found=(cur->getStartNode()==n);
      if(!found)
        it2.next();
    }
  if(!found)
    throw Exception("Internal error: polygons incompatible with each others. Should never happen!");
  //Ok we found correspondence between this and pol2. Searching for right direction to close polygon.
  ElementaryEdge *e=_sub_edges.back();
  if(e->getLoc()==FULL_ON_1)
    {
      if(e->getPtr()==cur->getPtr())
        {
          // if we have the same edge, several possibilities:
          //    - either this means that pol1 and pol2 have opposite orientation (since we matched end node with start node before)
          //    - or (more tricky, see testIntersect2DMeshes11()) that pol1 and pol2 have same orientation but 'this' turns in such
          //    a way that it attaches to pol2 on an edge in opposite orientation.
          // To sort this out, inspect localisation of next edge in pol2 wrt pol1NotSplitted.
          it2.nextLoop();
          cur=it2.current();
          Node *repr=cur->getPtr()->buildRepresentantOfMySelf();
          bool ret=pol1NotSplitted.isInOrOut(repr);
          repr->decrRef();
          direction = ret;
          needCleaning = ret; // if true we are in tricky case 2 above, we know that we will produce two consecutive overlapping edges in result
          return ret;
        }
      else  // here we don't need to go prev or next:
        {
          Node *repr=cur->getPtr()->buildRepresentantOfMySelf();
          bool ret=pol1NotSplitted.isInOrOut(repr);
          repr->decrRef();
          direction = ret;
          return ret;
        }
    }
  else
    direction=cur->locateFullyMySelfAbsolute(pol1NotSplitted)==FULL_IN_1;
  return true;
}

/*!
 * This method fills as much as possible \a this (a sub-part of pol1 split) with edges of \a pol2Splitted.
 */
std::list<QuadraticPolygon *>::iterator QuadraticPolygon::fillAsMuchAsPossibleWith(const QuadraticPolygon& pol2Splitted,
                                                                                   std::list<QuadraticPolygon *>::iterator iStart,
                                                                                   std::list<QuadraticPolygon *>::iterator iEnd,
                                                                                   bool direction)
{
  IteratorOnComposedEdge it1(const_cast<QuadraticPolygon *>(&pol2Splitted));
  bool found=false;
  Node *n=getEndNode();
  ElementaryEdge *cur1;
  for(it1.first();!it1.finished() && !found;)
    {
      cur1=it1.current();
      found=(cur1->getStartNode()==n);
      if(!found)
        it1.next();
    }
  if(!direction)
    it1.previousLoop();
  Node *nodeToTest;
  int szMax(pol2Splitted.size()+1),ii(0);   // protection against aggressive users of IntersectMeshes using invalid input meshes ...
  std::list<QuadraticPolygon *>::iterator ret;
  do
    { // Stack (consecutive) edges of pol1 into the result (no need to care about ordering, edges from pol1 are already consecutive)
      cur1=it1.current();
      ElementaryEdge *tmp=cur1->clone();
      if(!direction)
        tmp->reverse();
      pushBack(tmp);
      nodeToTest=tmp->getEndNode();
      direction?it1.nextLoop():it1.previousLoop();
      ret=CheckInList(nodeToTest,iStart,iEnd);
      if(completed())
        return iEnd;
      ii++;
    }
  while(ret==iEnd && ii<szMax);
  if(ii==szMax)// here a protection against aggressive users of IntersectMeshes of invalid input meshes
    throw INTERP_KERNEL::Exception("QuadraticPolygon::fillAsMuchAsPossibleWith : Something is invalid with input polygons !");
  return ret;
}

std::list<QuadraticPolygon *>::iterator QuadraticPolygon::CheckInList(Node *n, std::list<QuadraticPolygon *>::iterator iStart,
                                                                      std::list<QuadraticPolygon *>::iterator iEnd)
{
  for(std::list<QuadraticPolygon *>::iterator iter=iStart;iter!=iEnd;iter++)
    if((*iter)->isNodeIn(n))
      return iter;
  return iEnd;
}

/*!
* Compute the remaining parts of the intersection of mesh1 by mesh2.
* The general algorithm is to :
*   - either return full cells from pol1 that were simply not touched by mesh2
*   - or to:
*      - concatenate pieces from pol1 into consecutive pieces (a bit like zipConsecutiveSegments())
*      - complete those pieces by edges found in edgesInPol2OnBoundary, which are edges from pol2 located on the boundary of the previously built
*      intersecting cells
*/
void QuadraticPolygon::ComputeResidual(const QuadraticPolygon& pol1, const std::set<Edge *>& notUsedInPol1, const std::set<Edge *>& edgesInPol2OnBoundary,
                                       const std::map<INTERP_KERNEL::Node *,int>& mapp, int offset, int idThis,
                                       std::vector<double>& addCoordsQuadratic, std::vector<int>& conn,
                                       std::vector<int>& connI, std::vector<int>& nb1, std::vector<int>& nb2)
{
  // Initialise locations on pol1. Remember that edges found in 'notUsedInPol1' are also part of the edges forming pol1.
  pol1.initLocations();
  for(std::set<Edge *>::const_iterator it1=notUsedInPol1.begin();it1!=notUsedInPol1.end();it1++)
    { (*it1)->initLocs(); (*it1)->declareOn(); }
  for(std::set<Edge *>::const_iterator it2=edgesInPol2OnBoundary.begin();it2!=edgesInPol2OnBoundary.end();it2++)
    { (*it2)->initLocs(); (*it2)->declareIn(); }
  ////
  std::set<Edge *> notUsedInPol1L(notUsedInPol1);
  IteratorOnComposedEdge itPol1(const_cast<QuadraticPolygon *>(&pol1));
  int sz=pol1.size();
  std::list<QuadraticPolygon *> pol1Zip;
  // If none of the edges of pol1 was consumed by the rebuilding process, we can directly take pol1 as it is to form a cell:
  if(pol1.size()==(int)notUsedInPol1.size() && edgesInPol2OnBoundary.empty())
    {
      pol1.appendCrudeData(mapp,0.,0.,1.,offset,addCoordsQuadratic,conn,connI); nb1.push_back(idThis); nb2.push_back(-1);
      return ;
    }
  // Zip consecutive edges found in notUsedInPol1L and which are not overlapping with boundary edge from edgesInPol2OnBoundary:
  while(!notUsedInPol1L.empty())
    {
      // If all nodes are IN or edge is ON (i.e. as initialised at the begining of the method) then error
      for(int i=0;i<sz && (itPol1.current()->getStartNode()->getLoc()!=IN_1 || itPol1.current()->getLoc()!=FULL_ON_1);i++)
        itPol1.nextLoop();
      if(itPol1.current()->getStartNode()->getLoc()!=IN_1 || itPol1.current()->getLoc()!=FULL_ON_1)
        throw INTERP_KERNEL::Exception("Presence of a target polygon fully included in source polygon ! The partition of this leads to a non simply connex cell (with hole) ! Impossible ! Such resulting cell cannot be stored in MED cell format !");
      QuadraticPolygon *tmp1=new QuadraticPolygon;
      do
        {
          Edge *ee=itPol1.current()->getPtr();
          if(ee->getLoc()==FULL_ON_1)
            {
              ee->incrRef(); notUsedInPol1L.erase(ee);
              tmp1->pushBack(new ElementaryEdge(ee,itPol1.current()->getDirection()));
            }
          itPol1.nextLoop();
        }
      while(itPol1.current()->getStartNode()->getLoc()!=IN_1 && !notUsedInPol1L.empty());
      pol1Zip.push_back(tmp1);
    }

  ////
  std::list<QuadraticPolygon *> retPolsUnderContruction;
  std::list<Edge *> edgesInPol2OnBoundaryL(edgesInPol2OnBoundary.begin(),edgesInPol2OnBoundary.end());
  std::map<QuadraticPolygon *, std::list<QuadraticPolygon *> > pol1ZipConsumed;  // for memory management only.
  std::size_t maxNbOfTurn=edgesInPol2OnBoundaryL.size(),nbOfTurn=0,iiMNT=0;
  for(std::list<QuadraticPolygon *>::const_iterator itMNT=pol1Zip.begin();itMNT!=pol1Zip.end();itMNT++,iiMNT++)
    nbOfTurn+=(*itMNT)->size();
  maxNbOfTurn=maxNbOfTurn*nbOfTurn; maxNbOfTurn*=maxNbOfTurn;
  // [ABN] at least 3 turns for very small cases (typically one (quad) edge against one (quad or lin) edge forming a new cell)!
  maxNbOfTurn = maxNbOfTurn<3 ? 3 : maxNbOfTurn;
  nbOfTurn=0;
  while(nbOfTurn<maxNbOfTurn)  // the 'normal' way out of this loop is the break towards the end when pol1Zip is empty.
    {
      // retPolsUnderConstruction initially empty -> see if(!pol1Zip.empty()) below ...
      for(std::list<QuadraticPolygon *>::iterator itConstr=retPolsUnderContruction.begin();itConstr!=retPolsUnderContruction.end();)
        {
          if((*itConstr)->getStartNode()==(*itConstr)->getEndNode())  // reconstruction of a cell is finished
            {
              itConstr++;
              continue;
            }
          Node *curN=(*itConstr)->getEndNode();
          bool smthHappened=false;
          // Complete a partially reconstructed polygon with boundary edges by matching nodes:
          for(std::list<Edge *>::iterator it2=edgesInPol2OnBoundaryL.begin();it2!=edgesInPol2OnBoundaryL.end();)
            {
              if(curN==(*it2)->getStartNode())
                { (*it2)->incrRef(); (*itConstr)->pushBack(new ElementaryEdge(*it2,true)); curN=(*it2)->getEndNode(); smthHappened=true; it2=edgesInPol2OnBoundaryL.erase(it2); }
              else if(curN==(*it2)->getEndNode())
                { (*it2)->incrRef(); (*itConstr)->pushBack(new ElementaryEdge(*it2,false)); curN=(*it2)->getStartNode(); smthHappened=true; it2=edgesInPol2OnBoundaryL.erase(it2); }
              else
                it2++;
            }
          if(smthHappened)
            {
              for(std::list<QuadraticPolygon *>::iterator itZip=pol1Zip.begin();itZip!=pol1Zip.end();)
                {
                  if(curN==(*itZip)->getStartNode()) // we found a matching piece to append in pol1Zip. Append all of it to the current polygon being reconstr
                    {
                      for(std::list<ElementaryEdge *>::const_iterator it4=(*itZip)->_sub_edges.begin();it4!=(*itZip)->_sub_edges.end();it4++)
                        { (*it4)->getPtr()->incrRef(); bool dir=(*it4)->getDirection(); (*itConstr)->pushBack(new ElementaryEdge((*it4)->getPtr(),dir)); }
                      pol1ZipConsumed[*itConstr].push_back(*itZip);
                      curN=(*itZip)->getEndNode();
                      itZip=pol1Zip.erase(itZip);  // one zipped piece has been consumed
                      break;                       // we can stop here, pieces in pol1Zip are not connected, by definition.
                    }
                  else
                    itZip++;
                }
            }
          else
            {
              for(std::list<ElementaryEdge *>::const_iterator it5=(*itConstr)->_sub_edges.begin();it5!=(*itConstr)->_sub_edges.end();it5++)
                {
                  Edge *ee=(*it5)->getPtr();
                  if(edgesInPol2OnBoundary.find(ee)!=edgesInPol2OnBoundary.end())
                    edgesInPol2OnBoundaryL.push_back(ee);
                }
              for(std::list<QuadraticPolygon *>::iterator it6=pol1ZipConsumed[*itConstr].begin();it6!=pol1ZipConsumed[*itConstr].end();it6++)
                pol1Zip.push_front(*it6);
              pol1ZipConsumed.erase(*itConstr);
              delete *itConstr;
              itConstr=retPolsUnderContruction.erase(itConstr);
            }
        }
      if(!pol1Zip.empty())  // the filling process of retPolsUnderConstruction starts here
        {
          QuadraticPolygon *tmp=new QuadraticPolygon;
          QuadraticPolygon *first=*(pol1Zip.begin());
          for(std::list<ElementaryEdge *>::const_iterator it4=first->_sub_edges.begin();it4!=first->_sub_edges.end();it4++)
            { (*it4)->getPtr()->incrRef(); bool dir=(*it4)->getDirection(); tmp->pushBack(new ElementaryEdge((*it4)->getPtr(),dir)); }
          pol1ZipConsumed[tmp].push_back(first);
          retPolsUnderContruction.push_back(tmp);
          pol1Zip.erase(pol1Zip.begin());
        }
      else
        break;
      nbOfTurn++;
    }
  if(nbOfTurn==maxNbOfTurn)
    {
      std::ostringstream oss; oss << "Error during reconstruction of residual of cell ! It appears that either source or/and target mesh is/are not conform !";
      oss << " Number of turns is = " << nbOfTurn << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  // Convert to integer connectivity:
  for(std::list<QuadraticPolygon *>::iterator itConstr=retPolsUnderContruction.begin();itConstr!=retPolsUnderContruction.end();itConstr++)
    {
      if((*itConstr)->getStartNode()==(*itConstr)->getEndNode())  // take only fully closed reconstructed polygon
        {
          (*itConstr)->cleanDegeneratedConsecutiveEdges();
          (*itConstr)->appendCrudeData(mapp,0.,0.,1.,offset,addCoordsQuadratic,conn,connI); nb1.push_back(idThis); nb2.push_back(-1);
          for(std::list<QuadraticPolygon *>::iterator it6=pol1ZipConsumed[*itConstr].begin();it6!=pol1ZipConsumed[*itConstr].end();it6++)
            delete *it6;
          delete *itConstr;
        }
      else
        {
          std::ostringstream oss; oss << "Internal error during reconstruction of residual of cell! Non fully closed polygon built!";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}
