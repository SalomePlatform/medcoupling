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

#include "MEDCouplingUMesh.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "CellModel.hxx"
#include "VolSurfUser.txx"
#include "InterpolationUtils.hxx"
#include "PointLocatorAlgos.txx"
#include "BBTree.txx"
#include "BBTreeDst.txx"
#include "DirectedBoundingBox.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelAutoPtr.hxx"
#include "InterpKernelGeo2DNode.hxx"
#include "InterpKernelGeo2DEdgeLin.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DQuadraticPolygon.hxx"
#include "TranslationRotationMatrix.hxx"
#include "VectorUtils.hxx"
#include "MEDCouplingSkyLineArray.hxx"

#include <sstream>
#include <fstream>
#include <numeric>
#include <cstring>
#include <limits>
#include <list>
#include <assert.h>

using namespace MEDCoupling;

/// @cond INTERNAL

int InternalAddPoint(const INTERP_KERNEL::Edge *e, int id, const double *coo, int startId, int endId, DataArrayDouble& addCoo, int& nodesCnter)
{
  if(id!=-1)
    return id;
  else
    {
      int ret(nodesCnter++);
      double newPt[2];
      e->getMiddleOfPoints(coo+2*startId,coo+2*endId,newPt);
      addCoo.insertAtTheEnd(newPt,newPt+2);
      return ret;
    }
}

int InternalAddPointOriented(const INTERP_KERNEL::Edge *e, int id, const double *coo, int startId, int endId, DataArrayDouble& addCoo, int& nodesCnter)
{
  if(id!=-1)
    return id;
  else
    {
      int ret(nodesCnter++);
      double newPt[2];
      e->getMiddleOfPointsOriented(coo+2*startId,coo+2*endId,newPt);
      addCoo.insertAtTheEnd(newPt,newPt+2);
      return ret;
    }
}


void EnterTheResultOf2DCellFirst(const INTERP_KERNEL::Edge *e, int start, int stp, int nbOfEdges, bool linOrArc, const double *coords, const int *connBg, int offset, DataArrayInt *newConnOfCell, DataArrayDouble *appendedCoords, std::vector<int>& middles)
{
  int tmp[3];
  int trueStart(start>=0?start:nbOfEdges+start);
  tmp[0]=linOrArc?(int)INTERP_KERNEL::NORM_QPOLYG:(int)INTERP_KERNEL::NORM_POLYGON; tmp[1]=connBg[trueStart]; tmp[2]=connBg[stp];
  newConnOfCell->insertAtTheEnd(tmp,tmp+3);
  if(linOrArc)
    {
      if(stp-start>1)
        {
          int tmp2(0),tmp3(appendedCoords->getNumberOfTuples()/2);
          InternalAddPointOriented(e,-1,coords,tmp[1],tmp[2],*appendedCoords,tmp2);
          middles.push_back(tmp3+offset);
        }
      else
        middles.push_back(connBg[trueStart+nbOfEdges]);
    }
}

void EnterTheResultOf2DCellMiddle(const INTERP_KERNEL::Edge *e, int start, int stp, int nbOfEdges, bool linOrArc, const double *coords, const int *connBg, int offset, DataArrayInt *newConnOfCell, DataArrayDouble *appendedCoords, std::vector<int>& middles)
{
  int tmpSrt(newConnOfCell->back()),tmpEnd(connBg[stp]);
  newConnOfCell->pushBackSilent(tmpEnd);
  if(linOrArc)
    {
      if(stp-start>1)
        {
          int tmp2(0),tmp3(appendedCoords->getNumberOfTuples()/2);
          InternalAddPointOriented(e,-1,coords,tmpSrt,tmpEnd,*appendedCoords,tmp2);
          middles.push_back(tmp3+offset);
        }
      else
        middles.push_back(connBg[start+nbOfEdges]);
    }
}

void EnterTheResultOf2DCellEnd(const INTERP_KERNEL::Edge *e, int start, int stp, int nbOfEdges, bool linOrArc, const double *coords, const int *connBg, int offset, DataArrayInt *newConnOfCell, DataArrayDouble *appendedCoords, std::vector<int>& middles)
{
  // only the quadratic point to deal with:
  if(linOrArc)
    {
      if(stp-start>1)
        {
          int tmpSrt(connBg[start]),tmpEnd(connBg[stp]);
          int tmp2(0),tmp3(appendedCoords->getNumberOfTuples()/2);
          InternalAddPointOriented(e,-1,coords,tmpSrt,tmpEnd,*appendedCoords,tmp2);
          middles.push_back(tmp3+offset);
        }
      else
        middles.push_back(connBg[start+nbOfEdges]);
    }
}

void IKGeo2DInternalMapper2(INTERP_KERNEL::Node *n, const std::map<MCAuto<INTERP_KERNEL::Node>,int>& m, int forbVal0, int forbVal1, std::vector<int>& isect)
{
  MCAuto<INTERP_KERNEL::Node> nTmp(n); nTmp->incrRef();
  std::map<MCAuto<INTERP_KERNEL::Node>,int>::const_iterator it(m.find(nTmp));
  if(it==m.end())
    throw INTERP_KERNEL::Exception("Internal error in remapping !");
  int v((*it).second);
  if(v==forbVal0 || v==forbVal1)
    return ;
  if(std::find(isect.begin(),isect.end(),v)==isect.end())
    isect.push_back(v);
}

bool IKGeo2DInternalMapper(const INTERP_KERNEL::ComposedEdge& c, const std::map<MCAuto<INTERP_KERNEL::Node>,int>& m, int forbVal0, int forbVal1, std::vector<int>& isect)
{
  int sz(c.size());
  if(sz<=1)
    return false;
  bool presenceOfOn(false);
  for(int i=0;i<sz;i++)
    {
      INTERP_KERNEL::ElementaryEdge *e(c[i]);
      if(e->getLoc()!=INTERP_KERNEL::FULL_ON_1)
        continue ;
      IKGeo2DInternalMapper2(e->getStartNode(),m,forbVal0,forbVal1,isect);
      IKGeo2DInternalMapper2(e->getEndNode(),m,forbVal0,forbVal1,isect);
    }
  return presenceOfOn;
}


namespace MEDCoupling
{

  INTERP_KERNEL::Edge *MEDCouplingUMeshBuildQPFromEdge2(INTERP_KERNEL::NormalizedCellType typ, const int *bg, const double *coords2D, std::map< MCAuto<INTERP_KERNEL::Node>,int>& m)
  {
    INTERP_KERNEL::Edge *ret(0);
    MCAuto<INTERP_KERNEL::Node> n0(new INTERP_KERNEL::Node(coords2D[2*bg[0]],coords2D[2*bg[0]+1])),n1(new INTERP_KERNEL::Node(coords2D[2*bg[1]],coords2D[2*bg[1]+1]));
    m[n0]=bg[0]; m[n1]=bg[1];
    switch(typ)
    {
      case INTERP_KERNEL::NORM_SEG2:
        {
          ret=new INTERP_KERNEL::EdgeLin(n0,n1);
          break;
        }
      case INTERP_KERNEL::NORM_SEG3:
        {
          INTERP_KERNEL::Node *n2(new INTERP_KERNEL::Node(coords2D[2*bg[2]],coords2D[2*bg[2]+1])); m[n2]=bg[2];
          INTERP_KERNEL::EdgeLin *e1(new INTERP_KERNEL::EdgeLin(n0,n2)),*e2(new INTERP_KERNEL::EdgeLin(n2,n1));
          INTERP_KERNEL::SegSegIntersector inters(*e1,*e2);
          // is the SEG3 degenerated, and thus can be reduced to a SEG2?
          bool colinearity(inters.areColinears());
          delete e1; delete e2;
          if(colinearity)
            { ret=new INTERP_KERNEL::EdgeLin(n0,n1); }
          else
            { ret=new INTERP_KERNEL::EdgeArcCircle(n0,n2,n1); }
          break;
        }
      default:
        throw INTERP_KERNEL::Exception("MEDCouplingUMeshBuildQPFromEdge2 : Expecting a mesh with spaceDim==2 and meshDim==1 !");
    }
    return ret;
  }

  INTERP_KERNEL::Edge *MEDCouplingUMeshBuildQPFromEdge(INTERP_KERNEL::NormalizedCellType typ, std::map<int, std::pair<INTERP_KERNEL::Node *,bool> >& mapp2, const int *bg)
  {
    INTERP_KERNEL::Edge *ret=0;
    switch(typ)
    {
      case INTERP_KERNEL::NORM_SEG2:
        {
          ret=new INTERP_KERNEL::EdgeLin(mapp2[bg[0]].first,mapp2[bg[1]].first);
          break;
        }
      case INTERP_KERNEL::NORM_SEG3:
        {
          INTERP_KERNEL::EdgeLin *e1=new INTERP_KERNEL::EdgeLin(mapp2[bg[0]].first,mapp2[bg[2]].first);
          INTERP_KERNEL::EdgeLin *e2=new INTERP_KERNEL::EdgeLin(mapp2[bg[2]].first,mapp2[bg[1]].first);
          INTERP_KERNEL::SegSegIntersector inters(*e1,*e2);
          // is the SEG3 degenerated, and thus can be reduced to a SEG2?
          bool colinearity=inters.areColinears();
          delete e1; delete e2;
          if(colinearity)
            ret=new INTERP_KERNEL::EdgeLin(mapp2[bg[0]].first,mapp2[bg[1]].first);
          else
            ret=new INTERP_KERNEL::EdgeArcCircle(mapp2[bg[0]].first,mapp2[bg[2]].first,mapp2[bg[1]].first);
          mapp2[bg[2]].second=false;
          break;
        }
      default:
        throw INTERP_KERNEL::Exception("MEDCouplingUMeshBuildQPFromEdge : Expecting a mesh with spaceDim==2 and meshDim==1 !");
    }
    return ret;
  }

  /*!
   * This method creates a sub mesh in Geometric2D DS. The sub mesh is composed by the sub set of cells in 'candidates' taken from
   * the global mesh 'mDesc'.
   * The input mesh 'mDesc' must be so that mDim==1 and spaceDim==2.
   * 'mapp' returns a mapping between local numbering in submesh (represented by a Node*) and the global node numbering in 'mDesc'.
   */
  INTERP_KERNEL::QuadraticPolygon *MEDCouplingUMeshBuildQPFromMesh(const MEDCouplingUMesh *mDesc, const std::vector<int>& candidates,
                                                                   std::map<INTERP_KERNEL::Node *,int>& mapp)
  {
    mapp.clear();
    std::map<int, std::pair<INTERP_KERNEL::Node *,bool> > mapp2;//bool is for a flag specifying if node is boundary (true) or only a middle for SEG3.
    const double *coo=mDesc->getCoords()->getConstPointer();
    const int *c=mDesc->getNodalConnectivity()->getConstPointer();
    const int *cI=mDesc->getNodalConnectivityIndex()->getConstPointer();
    std::set<int> s;
    for(std::vector<int>::const_iterator it=candidates.begin();it!=candidates.end();it++)
      s.insert(c+cI[*it]+1,c+cI[(*it)+1]);
    for(std::set<int>::const_iterator it2=s.begin();it2!=s.end();it2++)
      {
        INTERP_KERNEL::Node *n=new INTERP_KERNEL::Node(coo[2*(*it2)],coo[2*(*it2)+1]);
        mapp2[*it2]=std::pair<INTERP_KERNEL::Node *,bool>(n,true);
      }
    INTERP_KERNEL::QuadraticPolygon *ret=new INTERP_KERNEL::QuadraticPolygon;
    for(std::vector<int>::const_iterator it=candidates.begin();it!=candidates.end();it++)
      {
        INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)c[cI[*it]];
        ret->pushBack(MEDCouplingUMeshBuildQPFromEdge(typ,mapp2,c+cI[*it]+1));
      }
    for(std::map<int, std::pair<INTERP_KERNEL::Node *,bool> >::const_iterator it2=mapp2.begin();it2!=mapp2.end();it2++)
      {
        if((*it2).second.second)
          mapp[(*it2).second.first]=(*it2).first;
        ((*it2).second.first)->decrRef();
      }
    return ret;
  }

  INTERP_KERNEL::Node *MEDCouplingUMeshBuildQPNode(int nodeId, const double *coo1, int offset1, const double *coo2, int offset2, const std::vector<double>& addCoo)
  {
    if(nodeId>=offset2)
      {
        int locId=nodeId-offset2;
        return new INTERP_KERNEL::Node(addCoo[2*locId],addCoo[2*locId+1]);
      }
    if(nodeId>=offset1)
      {
        int locId=nodeId-offset1;
        return new INTERP_KERNEL::Node(coo2[2*locId],coo2[2*locId+1]);
      }
    return new INTERP_KERNEL::Node(coo1[2*nodeId],coo1[2*nodeId+1]);
  }

  /**
   * Construct a mapping between set of Nodes and the standart MEDCoupling connectivity format (c, cI).
   */
  void MEDCouplingUMeshBuildQPFromMesh3(const double *coo1, int offset1, const double *coo2, int offset2, const std::vector<double>& addCoo,
                                        const int *desc1Bg, const int *desc1End, const std::vector<std::vector<int> >& intesctEdges1,
                                        /*output*/std::map<INTERP_KERNEL::Node *,int>& mapp, std::map<int,INTERP_KERNEL::Node *>& mappRev)
  {
    for(const int *desc1=desc1Bg;desc1!=desc1End;desc1++)
      {
        int eltId1=abs(*desc1)-1;
        for(std::vector<int>::const_iterator it1=intesctEdges1[eltId1].begin();it1!=intesctEdges1[eltId1].end();it1++)
          {
            std::map<int,INTERP_KERNEL::Node *>::const_iterator it=mappRev.find(*it1);
            if(it==mappRev.end())
              {
                INTERP_KERNEL::Node *node=MEDCouplingUMeshBuildQPNode(*it1,coo1,offset1,coo2,offset2,addCoo);
                mapp[node]=*it1;
                mappRev[*it1]=node;
              }
          }
      }
  }
}



/*!
 * Returns true if a colinearization has been found in the given cell. If false is returned the content pushed in \a newConnOfCell is equal to [ \a connBg , \a connEnd ) .
 * \a appendedCoords is a DataArrayDouble instance with number of components equal to one (even if the items are pushed by pair).
 */
bool MEDCouplingUMesh::Colinearize2DCell(const double *coords, const int *connBg, const int *connEnd, int offset, DataArrayInt *newConnOfCell, DataArrayDouble *appendedCoords)
{
  std::size_t sz(std::distance(connBg,connEnd));
  if(sz<3)//3 because 2+1(for the cell type) and 2 is the minimal number of edges of 2D cell.
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Colinearize2DCell : the input cell has invalid format !");
  sz--;
  INTERP_KERNEL::AutoPtr<int> tmpConn(new int[sz]);
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)connBg[0]));
  unsigned nbs(cm.getNumberOfSons2(connBg+1,sz));
  unsigned nbOfHit(0); // number of fusions operated
  int posBaseElt(0),posEndElt(0),nbOfTurn(0);
  const unsigned int maxNbOfHit = cm.isQuadratic() ? nbs-2 : nbs-3;  // a quad cell is authorized to end up with only two edges, a linear one has to keep 3 at least
  INTERP_KERNEL::NormalizedCellType typeOfSon;
  std::vector<int> middles;
  bool ret(false);
  for(;(nbOfTurn+nbOfHit)<nbs;nbOfTurn++)
    {
      cm.fillSonCellNodalConnectivity2(posBaseElt,connBg+1,sz,tmpConn,typeOfSon);
      std::map<MCAuto<INTERP_KERNEL::Node>,int> m;
      INTERP_KERNEL::Edge *e(MEDCouplingUMeshBuildQPFromEdge2(typeOfSon,tmpConn,coords,m));
      posEndElt = posBaseElt+1;

      // Look backward first: are the final edges of the cells colinear with the first ones?
      // This initializes posBaseElt.
      if(nbOfTurn==0)
        {
          for(unsigned i=1;i<nbs && nbOfHit<maxNbOfHit;i++) // 2nd condition is to avoid ending with a cell wih one single edge
            {
              cm.fillSonCellNodalConnectivity2(nbs-i,connBg+1,sz,tmpConn,typeOfSon);
              INTERP_KERNEL::Edge *eCand(MEDCouplingUMeshBuildQPFromEdge2(typeOfSon,tmpConn,coords,m));
              INTERP_KERNEL::EdgeIntersector *eint(INTERP_KERNEL::Edge::BuildIntersectorWith(e,eCand));
              bool isColinear=eint->areColinears();
              if(isColinear)
                {
                  nbOfHit++;
                  posBaseElt--;
                  ret=true;
                }
              delete eint;
              eCand->decrRef();
              if(!isColinear)
                break;
            }
        }
      // Now move forward:
      const unsigned fwdStart = (nbOfTurn == 0 ? 0 : posBaseElt);  // the first element to be inspected going forward
      for(unsigned j=fwdStart+1;j<nbs && nbOfHit<maxNbOfHit;j++)  // 2nd condition is to avoid ending with a cell wih one single edge
        {
          cm.fillSonCellNodalConnectivity2((int)j,connBg+1,sz,tmpConn,typeOfSon); // get edge #j's connectivity
          INTERP_KERNEL::Edge *eCand(MEDCouplingUMeshBuildQPFromEdge2(typeOfSon,tmpConn,coords,m));
          INTERP_KERNEL::EdgeIntersector *eint(INTERP_KERNEL::Edge::BuildIntersectorWith(e,eCand));
          bool isColinear(eint->areColinears());
          if(isColinear)
            {
              nbOfHit++;
              posEndElt++;
              ret=true;
            }
          delete eint;
          eCand->decrRef();
          if(!isColinear)
              break;
        }
      //push [posBaseElt,posEndElt) in newConnOfCell using e
      // The if clauses below are (volontary) not mutually exclusive: on a quad cell with 2 edges, the end of the connectivity is also its begining!
      if(nbOfTurn==0)
        // at the begining of the connectivity (insert type)
        EnterTheResultOf2DCellFirst(e,posBaseElt,posEndElt,(int)nbs,cm.isQuadratic(),coords,connBg+1,offset,newConnOfCell,appendedCoords,middles);
      else if((nbOfHit+nbOfTurn) != (nbs-1))
        // in the middle
        EnterTheResultOf2DCellMiddle(e,posBaseElt,posEndElt,(int)nbs,cm.isQuadratic(),coords,connBg+1,offset,newConnOfCell,appendedCoords,middles);
      if ((nbOfHit+nbOfTurn) == (nbs-1))
        // at the end (only quad points to deal with)
        EnterTheResultOf2DCellEnd(e,posBaseElt,posEndElt,(int)nbs,cm.isQuadratic(),coords,connBg+1,offset,newConnOfCell,appendedCoords,middles);
      posBaseElt=posEndElt;
      e->decrRef();
    }
  if(!middles.empty())
    newConnOfCell->insertAtTheEnd(middles.begin(),middles.end());
  return ret;
}



bool IsColinearOfACellOf(const std::vector< std::vector<int> >& intersectEdge1, const std::vector<int>& candidates, int start, int stop, int& retVal)
{
  if(candidates.empty())
    return false;
  for(std::vector<int>::const_iterator it=candidates.begin();it!=candidates.end();it++)
    {
      const std::vector<int>& pool(intersectEdge1[*it]);
      int tmp[2]; tmp[0]=start; tmp[1]=stop;
      if(std::search(pool.begin(),pool.end(),tmp,tmp+2)!=pool.end())
        {
          retVal=*it+1;
          return true;
        }
      tmp[0]=stop; tmp[1]=start;
      if(std::search(pool.begin(),pool.end(),tmp,tmp+2)!=pool.end())
        {
          retVal=-*it-1;
          return true;
        }
    }
  return false;
}

/*!
 * This method performs the 2nd step of Partition of 2D mesh.
 * This method has 4 inputs :
 *  - a mesh 'm1' with meshDim==1 and a SpaceDim==2
 *  - a mesh 'm2' with meshDim==1 and a SpaceDim==2
 *  - subDiv of size 'm2->getNumberOfCells()' that lists for each seg cell in 'm' the splitting node ids randomly sorted.
 * The aim of this method is to sort the splitting nodes, if any, and to put them in 'intersectEdge' output parameter based on edges of mesh 'm2'
 * Nodes end up lying consecutively on a cutted edge.
 * \param m1 is expected to be a mesh of meshDimension equal to 1 and spaceDim equal to 2. No check of that is performed by this method.
 * (Only present for its coords in case of 'subDiv' shares some nodes of 'm1')
 * \param m2 is expected to be a mesh of meshDimension equal to 1 and spaceDim equal to 2. No check of that is performed by this method.
 * \param addCoo input parameter with additional nodes linked to intersection of the 2 meshes.
 * \param[out] intersectEdge the same content as subDiv, but correclty oriented.
 */
void MEDCouplingUMesh::BuildIntersectEdges(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2,
                                           const std::vector<double>& addCoo,
                                           const std::vector< std::vector<int> >& subDiv, std::vector< std::vector<int> >& intersectEdge)
{
  int offset1=m1->getNumberOfNodes();
  int ncell=m2->getNumberOfCells();
  const int *c=m2->getNodalConnectivity()->begin();
  const int *cI=m2->getNodalConnectivityIndex()->begin();
  const double *coo=m2->getCoords()->begin();
  const double *cooBis=m1->getCoords()->begin();
  int offset2=offset1+m2->getNumberOfNodes();
  intersectEdge.resize(ncell);
  for(int i=0;i<ncell;i++,cI++)
    {
      const std::vector<int>& divs=subDiv[i];
      int nnode=cI[1]-cI[0]-1;
      std::map<int, std::pair<INTERP_KERNEL::Node *,bool> > mapp2;
      std::map<INTERP_KERNEL::Node *, int> mapp22;
      for(int j=0;j<nnode;j++)
        {
          INTERP_KERNEL::Node *nn=new INTERP_KERNEL::Node(coo[2*c[(*cI)+j+1]],coo[2*c[(*cI)+j+1]+1]);
          int nnid=c[(*cI)+j+1];
          mapp2[nnid]=std::pair<INTERP_KERNEL::Node *,bool>(nn,true);
          mapp22[nn]=nnid+offset1;
        }
      INTERP_KERNEL::Edge *e=MEDCouplingUMeshBuildQPFromEdge((INTERP_KERNEL::NormalizedCellType)c[*cI],mapp2,c+(*cI)+1);
      for(std::map<int, std::pair<INTERP_KERNEL::Node *,bool> >::const_iterator it=mapp2.begin();it!=mapp2.end();it++)
        ((*it).second.first)->decrRef();
      std::vector<INTERP_KERNEL::Node *> addNodes(divs.size());
      std::map<INTERP_KERNEL::Node *,int> mapp3;
      for(std::size_t j=0;j<divs.size();j++)
        {
          int id=divs[j];
          INTERP_KERNEL::Node *tmp=0;
          if(id<offset1)
            tmp=new INTERP_KERNEL::Node(cooBis[2*id],cooBis[2*id+1]);
          else if(id<offset2)
            tmp=new INTERP_KERNEL::Node(coo[2*(id-offset1)],coo[2*(id-offset1)+1]);//if it happens, bad news mesh 'm2' is non conform.
          else
            tmp=new INTERP_KERNEL::Node(addCoo[2*(id-offset2)],addCoo[2*(id-offset2)+1]);
          addNodes[j]=tmp;
          mapp3[tmp]=id;
        }
      e->sortIdsAbs(addNodes,mapp22,mapp3,intersectEdge[i]);
      for(std::vector<INTERP_KERNEL::Node *>::const_iterator it=addNodes.begin();it!=addNodes.end();it++)
        (*it)->decrRef();
      e->decrRef();
    }
}

MEDCouplingUMesh *BuildMesh1DCutFrom(const MEDCouplingUMesh *mesh1D, const std::vector< std::vector<int> >& intersectEdge2, const DataArrayDouble *coords1, const std::vector<double>& addCoo, const std::map<int,int>& mergedNodes, const std::vector< std::vector<int> >& colinear2, const std::vector< std::vector<int> >& intersectEdge1,
                                     MCAuto<DataArrayInt>& idsInRetColinear, MCAuto<DataArrayInt>& idsInMesh1DForIdsInRetColinear)
{
  idsInRetColinear=DataArrayInt::New(); idsInRetColinear->alloc(0,1);
  idsInMesh1DForIdsInRetColinear=DataArrayInt::New(); idsInMesh1DForIdsInRetColinear->alloc(0,1);
  int nCells(mesh1D->getNumberOfCells());
  if(nCells!=(int)intersectEdge2.size())
    throw INTERP_KERNEL::Exception("BuildMesh1DCutFrom : internal error # 1 !");
  const DataArrayDouble *coo2(mesh1D->getCoords());
  const int *c(mesh1D->getNodalConnectivity()->begin()),*ci(mesh1D->getNodalConnectivityIndex()->begin());
  const double *coo2Ptr(coo2->begin());
  int offset1(coords1->getNumberOfTuples());
  int offset2(offset1+coo2->getNumberOfTuples());
  int offset3(offset2+addCoo.size()/2);
  std::vector<double> addCooQuad;
  MCAuto<DataArrayInt> cOut(DataArrayInt::New()),ciOut(DataArrayInt::New()); cOut->alloc(0,1); ciOut->alloc(1,1); ciOut->setIJ(0,0,0);
  int tmp[4],cicnt(0),kk(0);
  for(int i=0;i<nCells;i++)
    {
      std::map<MCAuto<INTERP_KERNEL::Node>,int> m;
      INTERP_KERNEL::Edge *e(MEDCouplingUMeshBuildQPFromEdge2((INTERP_KERNEL::NormalizedCellType)c[ci[i]],c+ci[i]+1,coo2Ptr,m));
      const std::vector<int>& subEdges(intersectEdge2[i]);
      int nbSubEdge(subEdges.size()/2);
      for(int j=0;j<nbSubEdge;j++,kk++)
        {
          MCAuto<INTERP_KERNEL::Node> n1(MEDCouplingUMeshBuildQPNode(subEdges[2*j],coords1->begin(),offset1,coo2Ptr,offset2,addCoo)),n2(MEDCouplingUMeshBuildQPNode(subEdges[2*j+1],coords1->begin(),offset1,coo2Ptr,offset2,addCoo));
          MCAuto<INTERP_KERNEL::Edge> e2(e->buildEdgeLyingOnMe(n1,n2));
          INTERP_KERNEL::Edge *e2Ptr(e2);
          std::map<int,int>::const_iterator itm;
          if(dynamic_cast<INTERP_KERNEL::EdgeArcCircle *>(e2Ptr))
            {
              tmp[0]=INTERP_KERNEL::NORM_SEG3;
              itm=mergedNodes.find(subEdges[2*j]);
              tmp[1]=itm!=mergedNodes.end()?(*itm).second:subEdges[2*j];
              itm=mergedNodes.find(subEdges[2*j+1]);
              tmp[2]=itm!=mergedNodes.end()?(*itm).second:subEdges[2*j+1];
              tmp[3]=offset3+(int)addCooQuad.size()/2;
              double tmp2[2];
              e2->getBarycenter(tmp2); addCooQuad.insert(addCooQuad.end(),tmp2,tmp2+2);
              cicnt+=4;
              cOut->insertAtTheEnd(tmp,tmp+4);
              ciOut->pushBackSilent(cicnt);
            }
          else
            {
              tmp[0]=INTERP_KERNEL::NORM_SEG2;
              itm=mergedNodes.find(subEdges[2*j]);
              tmp[1]=itm!=mergedNodes.end()?(*itm).second:subEdges[2*j];
              itm=mergedNodes.find(subEdges[2*j+1]);
              tmp[2]=itm!=mergedNodes.end()?(*itm).second:subEdges[2*j+1];
              cicnt+=3;
              cOut->insertAtTheEnd(tmp,tmp+3);
              ciOut->pushBackSilent(cicnt);
            }
          int tmp00;
          if(IsColinearOfACellOf(intersectEdge1,colinear2[i],tmp[1],tmp[2],tmp00))
            {
              idsInRetColinear->pushBackSilent(kk);
              idsInMesh1DForIdsInRetColinear->pushBackSilent(tmp00);
            }
        }
      e->decrRef();
    }
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::New(mesh1D->getName(),1));
  ret->setConnectivity(cOut,ciOut,true);
  MCAuto<DataArrayDouble> arr3(DataArrayDouble::New());
  arr3->useArray(&addCoo[0],false,C_DEALLOC,(int)addCoo.size()/2,2);
  MCAuto<DataArrayDouble> arr4(DataArrayDouble::New()); arr4->useArray(&addCooQuad[0],false,C_DEALLOC,(int)addCooQuad.size()/2,2);
  std::vector<const DataArrayDouble *> coordss(4);
  coordss[0]=coords1; coordss[1]=mesh1D->getCoords(); coordss[2]=arr3; coordss[3]=arr4;
  MCAuto<DataArrayDouble> arr(DataArrayDouble::Aggregate(coordss));
  ret->setCoords(arr);
  return ret.retn();
}

MEDCouplingUMesh *BuildRefined2DCellLinear(const DataArrayDouble *coords, const int *descBg, const int *descEnd, const std::vector< std::vector<int> >& intersectEdge1)
{
  std::vector<int> allEdges;
  for(const int *it2(descBg);it2!=descEnd;it2++)
    {
      const std::vector<int>& edge1(intersectEdge1[std::abs(*it2)-1]);
      if(*it2>0)
        allEdges.insert(allEdges.end(),edge1.begin(),edge1.end());
      else
        allEdges.insert(allEdges.end(),edge1.rbegin(),edge1.rend());
    }
  std::size_t nb(allEdges.size());
  if(nb%2!=0)
    throw INTERP_KERNEL::Exception("BuildRefined2DCellLinear : internal error 1 !");
  std::size_t nbOfEdgesOf2DCellSplit(nb/2);
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::New("",2));
  ret->setCoords(coords);
  ret->allocateCells(1);
  std::vector<int> connOut(nbOfEdgesOf2DCellSplit);
  for(std::size_t kk=0;kk<nbOfEdgesOf2DCellSplit;kk++)
    connOut[kk]=allEdges[2*kk];
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYGON,connOut.size(),&connOut[0]);
  return ret.retn();
}

MEDCouplingUMesh *BuildRefined2DCellQuadratic(const DataArrayDouble *coords, const MEDCouplingUMesh *mesh2D, int cellIdInMesh2D, const int *descBg, const int *descEnd, const std::vector< std::vector<int> >& intersectEdge1)
{
  const int *c(mesh2D->getNodalConnectivity()->begin()),*ci(mesh2D->getNodalConnectivityIndex()->begin());
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)c[ci[cellIdInMesh2D]]));
  std::size_t ii(0);
  unsigned sz(cm.getNumberOfSons2(c+ci[cellIdInMesh2D]+1,ci[cellIdInMesh2D+1]-ci[cellIdInMesh2D]-1));
  if(sz!=std::distance(descBg,descEnd))
    throw INTERP_KERNEL::Exception("BuildRefined2DCellQuadratic : internal error 1 !");
  INTERP_KERNEL::AutoPtr<int> tmpPtr(new int[ci[cellIdInMesh2D+1]-ci[cellIdInMesh2D]]);
  std::vector<int> allEdges,centers;
  const double *coordsPtr(coords->begin());
  MCAuto<DataArrayDouble> addCoo(DataArrayDouble::New()); addCoo->alloc(0,1);
  int offset(coords->getNumberOfTuples());
  for(const int *it2(descBg);it2!=descEnd;it2++,ii++)
    {
      INTERP_KERNEL::NormalizedCellType typeOfSon;
      cm.fillSonCellNodalConnectivity2(ii,c+ci[cellIdInMesh2D]+1,ci[cellIdInMesh2D+1]-ci[cellIdInMesh2D]-1,tmpPtr,typeOfSon);
      const std::vector<int>& edge1(intersectEdge1[std::abs(*it2)-1]);
      if(*it2>0)
        allEdges.insert(allEdges.end(),edge1.begin(),edge1.end());
      else
        allEdges.insert(allEdges.end(),edge1.rbegin(),edge1.rend());
      if(edge1.size()==2)
        centers.push_back(tmpPtr[2]);//special case where no subsplit of edge -> reuse the original center.
      else
        {//the current edge has been subsplit -> create corresponding centers.
          std::size_t nbOfCentersToAppend(edge1.size()/2);
          std::map< MCAuto<INTERP_KERNEL::Node>,int> m;
          MCAuto<INTERP_KERNEL::Edge> ee(MEDCouplingUMeshBuildQPFromEdge2(typeOfSon,tmpPtr,coordsPtr,m));
          std::vector<int>::const_iterator it3(allEdges.end()-edge1.size());
          for(std::size_t k=0;k<nbOfCentersToAppend;k++)
            {
              double tmpp[2];
              const double *aa(coordsPtr+2*(*it3++));
              const double *bb(coordsPtr+2*(*it3++));
              ee->getMiddleOfPoints(aa,bb,tmpp);
              addCoo->insertAtTheEnd(tmpp,tmpp+2);
              centers.push_back(offset+k);
            }
        }
    }
  std::size_t nb(allEdges.size());
  if(nb%2!=0)
    throw INTERP_KERNEL::Exception("BuildRefined2DCellQuadratic : internal error 2 !");
  std::size_t nbOfEdgesOf2DCellSplit(nb/2);
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::New("",2));
  if(addCoo->empty())
    ret->setCoords(coords);
  else
    {
      addCoo->rearrange(2);
      addCoo=DataArrayDouble::Aggregate(coords,addCoo);
      ret->setCoords(addCoo);
    }
  ret->allocateCells(1);
  std::vector<int> connOut(nbOfEdgesOf2DCellSplit);
  for(std::size_t kk=0;kk<nbOfEdgesOf2DCellSplit;kk++)
    connOut[kk]=allEdges[2*kk];
  connOut.insert(connOut.end(),centers.begin(),centers.end());
  ret->insertNextCell(INTERP_KERNEL::NORM_QPOLYG,connOut.size(),&connOut[0]);
  return ret.retn();
}

/*!
 * This method creates a refinement of a cell in \a mesh2D. Those cell is defined by descending connectivity and the sorted subdivided nodal connectivity
 * of those edges.
 *
 * \param [in] mesh2D - The origin 2D mesh. \b Warning \b coords are not those of \a mesh2D. But mesh2D->getCoords()==coords[:mesh2D->getNumberOfNodes()]
 */
MEDCouplingUMesh *BuildRefined2DCell(const DataArrayDouble *coords, const MEDCouplingUMesh *mesh2D, int cellIdInMesh2D, const int *descBg, const int *descEnd, const std::vector< std::vector<int> >& intersectEdge1)
{
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(mesh2D->getTypeOfCell(cellIdInMesh2D)));
  if(!cm.isQuadratic())
    return BuildRefined2DCellLinear(coords,descBg,descEnd,intersectEdge1);
  else
    return BuildRefined2DCellQuadratic(coords,mesh2D,cellIdInMesh2D,descBg,descEnd,intersectEdge1);
}

void AddCellInMesh2D(MEDCouplingUMesh *mesh2D, const std::vector<int>& conn, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& edges)
{
  bool isQuad(false);
  for(std::vector< MCAuto<INTERP_KERNEL::Edge> >::const_iterator it=edges.begin();it!=edges.end();it++)
    {
      const INTERP_KERNEL::Edge *ee(*it);
      if(dynamic_cast<const INTERP_KERNEL::EdgeArcCircle *>(ee))
        isQuad=true;
    }
  if(!isQuad)
    mesh2D->insertNextCell(INTERP_KERNEL::NORM_POLYGON,conn.size(),&conn[0]);
  else
    {
      const double *coo(mesh2D->getCoords()->begin());
      std::size_t sz(conn.size());
      std::vector<double> addCoo;
      std::vector<int> conn2(conn);
      int offset(mesh2D->getNumberOfNodes());
      for(std::size_t i=0;i<sz;i++)
        {
          double tmp[2];
          edges[(i+1)%sz]->getMiddleOfPoints(coo+2*conn[i],coo+2*conn[(i+1)%sz],tmp);// tony a chier i+1 -> i
          addCoo.insert(addCoo.end(),tmp,tmp+2);
          conn2.push_back(offset+(int)i);
        }
      mesh2D->getCoords()->rearrange(1);
      mesh2D->getCoords()->pushBackValsSilent(&addCoo[0],&addCoo[0]+addCoo.size());
      mesh2D->getCoords()->rearrange(2);
      mesh2D->insertNextCell(INTERP_KERNEL::NORM_QPOLYG,conn2.size(),&conn2[0]);
    }
}

/*!
 * \b WARNING edges in out1 coming from \a splitMesh1D are \b NOT oriented because only used for equation of curve.
 *
 * This method cuts in 2 parts the input 2D cell given using boundaries description (\a edge1Bis and \a edge1BisPtr) using
 * a set of edges defined in \a splitMesh1D.
 */
void BuildMesh2DCutInternal2(const MEDCouplingUMesh *splitMesh1D, const std::vector<int>& edge1Bis, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& edge1BisPtr,
                             std::vector< std::vector<int> >& out0, std::vector< std::vector< MCAuto<INTERP_KERNEL::Edge> > >& out1)
{
  std::size_t nb(edge1Bis.size()/2);
  std::size_t nbOfEdgesOf2DCellSplit(nb/2);
  int iEnd(splitMesh1D->getNumberOfCells());
  if(iEnd==0)
    throw INTERP_KERNEL::Exception("BuildMesh2DCutInternal2 : internal error ! input 1D mesh must have at least one cell !");
  std::size_t ii,jj;
  const int *cSplitPtr(splitMesh1D->getNodalConnectivity()->begin()),*ciSplitPtr(splitMesh1D->getNodalConnectivityIndex()->begin());
  for(ii=0;ii<nb && edge1Bis[2*ii]!=cSplitPtr[ciSplitPtr[0]+1];ii++);
  for(jj=ii;jj<nb && edge1Bis[2*jj+1]!=cSplitPtr[ciSplitPtr[iEnd-1]+2];jj++);
  //
  if(jj==nb)
    {//the edges splitMesh1D[iStart:iEnd] does not fully cut the current 2D cell -> single output cell
      out0.resize(1); out1.resize(1);
      std::vector<int>& connOut(out0[0]);
      connOut.resize(nbOfEdgesOf2DCellSplit);
      std::vector< MCAuto<INTERP_KERNEL::Edge> >& edgesPtr(out1[0]);
      edgesPtr.resize(nbOfEdgesOf2DCellSplit);
      for(std::size_t kk=0;kk<nbOfEdgesOf2DCellSplit;kk++)
        {
          connOut[kk]=edge1Bis[2*kk];
          edgesPtr[kk]=edge1BisPtr[2*kk];
        }
    }
  else
    {
      // [i,iEnd[ contains the
      out0.resize(2); out1.resize(2);
      std::vector<int>& connOutLeft(out0[0]);
      std::vector<int>& connOutRight(out0[1]);//connOutLeft should end with edge1Bis[2*ii] and connOutRight should end with edge1Bis[2*jj+1]
      std::vector< MCAuto<INTERP_KERNEL::Edge> >& eleft(out1[0]);
      std::vector< MCAuto<INTERP_KERNEL::Edge> >& eright(out1[1]);
      for(std::size_t k=ii;k<jj+1;k++)
        { connOutLeft.push_back(edge1Bis[2*k+1]); eleft.push_back(edge1BisPtr[2*k+1]); }
      std::vector< MCAuto<INTERP_KERNEL::Edge> > ees(iEnd);
      for(int ik=0;ik<iEnd;ik++)
        {
          std::map< MCAuto<INTERP_KERNEL::Node>,int> m;
          MCAuto<INTERP_KERNEL::Edge> ee(MEDCouplingUMeshBuildQPFromEdge2((INTERP_KERNEL::NormalizedCellType)cSplitPtr[ciSplitPtr[ik]],cSplitPtr+ciSplitPtr[ik]+1,splitMesh1D->getCoords()->begin(),m));
          ees[ik]=ee;
        }
      for(int ik=iEnd-1;ik>=0;ik--)
        connOutLeft.push_back(cSplitPtr[ciSplitPtr[ik]+1]);
      for(std::size_t k=jj+1;k<nbOfEdgesOf2DCellSplit+ii;k++)
        { connOutRight.push_back(edge1Bis[2*k+1]); eright.push_back(edge1BisPtr[2*k+1]); }
      eleft.insert(eleft.end(),ees.rbegin(),ees.rend());
      for(int ik=0;ik<iEnd;ik++)
        connOutRight.push_back(cSplitPtr[ciSplitPtr[ik]+2]);
      eright.insert(eright.end(),ees.begin(),ees.end());
    }
}

struct CellInfo
{
public:
  CellInfo() { }
  CellInfo(const std::vector<int>& edges, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& edgesPtr);
public:
  std::vector<int> _edges;
  std::vector< MCAuto<INTERP_KERNEL::Edge> > _edges_ptr;
};

CellInfo::CellInfo(const std::vector<int>& edges, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& edgesPtr)
{
  std::size_t nbe(edges.size());
  std::vector<int> edges2(2*nbe); std::vector< MCAuto<INTERP_KERNEL::Edge> > edgesPtr2(2*nbe);
  for(std::size_t i=0;i<nbe;i++)
    {
      edges2[2*i]=edges[i]; edges2[2*i+1]=edges[(i+1)%nbe];
      edgesPtr2[2*i]=edgesPtr[(i+1)%nbe]; edgesPtr2[2*i+1]=edgesPtr[(i+1)%nbe];//tony a chier
    }
  _edges.resize(4*nbe); _edges_ptr.resize(4*nbe);
  std::copy(edges2.begin(),edges2.end(),_edges.begin()); std::copy(edges2.begin(),edges2.end(),_edges.begin()+2*nbe);
  std::copy(edgesPtr2.begin(),edgesPtr2.end(),_edges_ptr.begin()); std::copy(edgesPtr2.begin(),edgesPtr2.end(),_edges_ptr.begin()+2*nbe);
}

class EdgeInfo
{
public:
  EdgeInfo(int istart, int iend, const MCAuto<MEDCouplingUMesh>& mesh):_istart(istart),_iend(iend),_mesh(mesh),_left(-7),_right(-7) { }
  EdgeInfo(int istart, int iend, int pos, const MCAuto<INTERP_KERNEL::Edge>& edge):_istart(istart),_iend(iend),_edge(edge),_left(pos),_right(pos+1) { }
  bool isInMyRange(int pos) const { return pos>=_istart && pos<_iend; }
  void somethingHappendAt(int pos, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newLeft, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newRight);
  void feedEdgeInfoAt(double eps, const MEDCouplingUMesh *mesh2D, int offset, int neighbors[2]) const;
private:
  int _istart;
  int _iend;
  MCAuto<MEDCouplingUMesh> _mesh;
  MCAuto<INTERP_KERNEL::Edge> _edge;
  int _left;
  int _right;
};

void EdgeInfo::somethingHappendAt(int pos, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newLeft, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newRight)
{
  const MEDCouplingUMesh *mesh(_mesh);
  if(mesh)
    return ;
  if(_right<pos)
    return ;
  if(_left>pos)
    { _left++; _right++; return ; }
  if(_right==pos)
    {
      bool isLeft(std::find(newLeft.begin(),newLeft.end(),_edge)!=newLeft.end()),isRight(std::find(newRight.begin(),newRight.end(),_edge)!=newRight.end());
      if((isLeft && isRight) || (!isLeft && !isRight))
        throw INTERP_KERNEL::Exception("EdgeInfo::somethingHappendAt : internal error # 1 !");
      if(isLeft)
        return ;
      if(isRight)
        {
          _right++;
          return ;
        }
    }
  if(_left==pos)
    {
      bool isLeft(std::find(newLeft.begin(),newLeft.end(),_edge)!=newLeft.end()),isRight(std::find(newRight.begin(),newRight.end(),_edge)!=newRight.end());
      if((isLeft && isRight) || (!isLeft && !isRight))
        throw INTERP_KERNEL::Exception("EdgeInfo::somethingHappendAt : internal error # 2 !");
      if(isLeft)
        {
          _right++;
          return ;
        }
      if(isRight)
        {
          _left++;
          _right++;
          return ;
        }
    }
}

void EdgeInfo::feedEdgeInfoAt(double eps, const MEDCouplingUMesh *mesh2D, int offset, int neighbors[2]) const
{
  const MEDCouplingUMesh *mesh(_mesh);
  if(!mesh)
    {
      neighbors[0]=offset+_left; neighbors[1]=offset+_right;
    }
  else
    {// not fully splitting cell case
      if(mesh2D->getNumberOfCells()==1)
        {//little optimization. 1 cell no need to find in which cell mesh is !
          neighbors[0]=offset; neighbors[1]=offset;
          return;
        }
      else
        {
          MCAuto<DataArrayDouble> barys(mesh->computeCellCenterOfMass());
          int cellId(mesh2D->getCellContainingPoint(barys->begin(),eps));
          if(cellId==-1)
            throw INTERP_KERNEL::Exception("EdgeInfo::feedEdgeInfoAt : internal error !");
          neighbors[0]=offset+cellId; neighbors[1]=offset+cellId;
        }
    }
}

class VectorOfCellInfo
{
public:
  VectorOfCellInfo(const std::vector<int>& edges, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& edgesPtr);
  std::size_t size() const { return _pool.size(); }
  int getPositionOf(double eps, const MEDCouplingUMesh *mesh) const;
  void setMeshAt(std::size_t pos, const MCAuto<MEDCouplingUMesh>& mesh, int istart, int iend, const MCAuto<MEDCouplingUMesh>& mesh1DInCase, const std::vector< std::vector<int> >& edges, const std::vector< std::vector< MCAuto<INTERP_KERNEL::Edge> > >& edgePtrs);
  const std::vector<int>& getConnOf(int pos) const { return get(pos)._edges; }
  const std::vector< MCAuto<INTERP_KERNEL::Edge> >& getEdgePtrOf(int pos) const { return get(pos)._edges_ptr; }
  MCAuto<MEDCouplingUMesh> getZeMesh() const { return _ze_mesh; }
  void feedEdgeInfoAt(double eps, int pos, int offset, int neighbors[2]) const;
private:
  int getZePosOfEdgeGivenItsGlobalId(int pos) const;
  void updateEdgeInfo(int pos, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newLeft, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newRight);
  const CellInfo& get(int pos) const;
  CellInfo& get(int pos);
private:
  std::vector<CellInfo> _pool;
  MCAuto<MEDCouplingUMesh> _ze_mesh;
  std::vector<EdgeInfo> _edge_info;
};

VectorOfCellInfo::VectorOfCellInfo(const std::vector<int>& edges, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& edgesPtr):_pool(1)
{
  _pool[0]._edges=edges;
  _pool[0]._edges_ptr=edgesPtr;
}

int VectorOfCellInfo::getPositionOf(double eps, const MEDCouplingUMesh *mesh) const
{
  if(_pool.empty())
    throw INTERP_KERNEL::Exception("VectorOfCellSplitter::getPositionOf : empty !");
  if(_pool.size()==1)
    return 0;
  const MEDCouplingUMesh *zeMesh(_ze_mesh);
  if(!zeMesh)
    throw INTERP_KERNEL::Exception("VectorOfCellSplitter::getPositionOf : null aggregated mesh !");
  MCAuto<DataArrayDouble> barys(mesh->computeCellCenterOfMass());
  return zeMesh->getCellContainingPoint(barys->begin(),eps);
}

void VectorOfCellInfo::setMeshAt(std::size_t pos, const MCAuto<MEDCouplingUMesh>& mesh, int istart, int iend, const MCAuto<MEDCouplingUMesh>& mesh1DInCase, const std::vector< std::vector<int> >& edges, const std::vector< std::vector< MCAuto<INTERP_KERNEL::Edge> > >& edgePtrs)
{
  get(pos);//to check pos
  bool isFast(pos==0 && _pool.size()==1);
  std::size_t sz(edges.size());
  // dealing with edges
  if(sz==1)
    _edge_info.push_back(EdgeInfo(istart,iend,mesh1DInCase));
  else
    _edge_info.push_back(EdgeInfo(istart,iend,pos,edgePtrs[0].back()));
  //
  std::vector<CellInfo> pool(_pool.size()-1+sz);
  for(std::size_t i=0;i<pos;i++)
    pool[i]=_pool[i];
  for(std::size_t j=0;j<sz;j++)
    pool[pos+j]=CellInfo(edges[j],edgePtrs[j]);
  for(int i=pos+1;i<(int)_pool.size();i++)
    pool[i+sz-1]=_pool[i];
  _pool=pool;
  //
  if(sz==2)
    updateEdgeInfo(pos,edgePtrs[0],edgePtrs[1]);
  //
  if(isFast)
    {
      _ze_mesh=mesh;
      return ;
    }
  //
  std::vector< MCAuto<MEDCouplingUMesh> > ms;
  if(pos>0)
    {
      MCAuto<MEDCouplingUMesh> elt(static_cast<MEDCouplingUMesh *>(_ze_mesh->buildPartOfMySelfSlice(0,pos,true)));
      ms.push_back(elt);
    }
  ms.push_back(mesh);
  if(pos<_ze_mesh->getNumberOfCells()-1)
  {
    MCAuto<MEDCouplingUMesh> elt(static_cast<MEDCouplingUMesh *>(_ze_mesh->buildPartOfMySelfSlice(pos+1,_ze_mesh->getNumberOfCells(),true)));
    ms.push_back(elt);
  }
  std::vector< const MEDCouplingUMesh *> ms2(ms.size());
  for(std::size_t j=0;j<ms2.size();j++)
    ms2[j]=ms[j];
  _ze_mesh=MEDCouplingUMesh::MergeUMeshesOnSameCoords(ms2);
}

void VectorOfCellInfo::feedEdgeInfoAt(double eps, int pos, int offset, int neighbors[2]) const
{
  _edge_info[getZePosOfEdgeGivenItsGlobalId(pos)].feedEdgeInfoAt(eps,_ze_mesh,offset,neighbors);
}

int VectorOfCellInfo::getZePosOfEdgeGivenItsGlobalId(int pos) const
{
  if(pos<0)
    throw INTERP_KERNEL::Exception("VectorOfCellInfo::getZePosOfEdgeGivenItsGlobalId : invalid id ! Must be >=0 !");
  int ret(0);
  for(std::vector<EdgeInfo>::const_iterator it=_edge_info.begin();it!=_edge_info.end();it++,ret++)
    {
      if((*it).isInMyRange(pos))
        return ret;
    }
  throw INTERP_KERNEL::Exception("VectorOfCellInfo::getZePosOfEdgeGivenItsGlobalId : invalid id !");
}

void VectorOfCellInfo::updateEdgeInfo(int pos, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newLeft, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newRight)
{
  get(pos);//to check;
  if(_edge_info.empty())
    return ;
  std::size_t sz(_edge_info.size()-1);
  for(std::size_t i=0;i<sz;i++)
    _edge_info[i].somethingHappendAt(pos,newLeft,newRight);
}

const CellInfo& VectorOfCellInfo::get(int pos) const
{
  if(pos<0 || pos>=(int)_pool.size())
    throw INTERP_KERNEL::Exception("VectorOfCellSplitter::get const : invalid pos !");
  return _pool[pos];
}

CellInfo& VectorOfCellInfo::get(int pos)
{
  if(pos<0 || pos>=(int)_pool.size())
    throw INTERP_KERNEL::Exception("VectorOfCellSplitter::get : invalid pos !");
  return _pool[pos];
}

/*!
 * Given :
 * - a \b closed set of edges ( \a allEdges and \a allEdgesPtr ) that defines the split descending 2D cell.
 * - \a splitMesh1D a split 2D curve mesh contained into 2D cell defined above.
 *
 * This method returns the 2D mesh and feeds \a idsLeftRight using offset.
 *
 * Algorithm : \a splitMesh1D is cut into contiguous parts. Each contiguous parts will build incrementally the output 2D cells.
 *
 * \param [in] allEdges a list of pairs (beginNode, endNode). Linked with \a allEdgesPtr to get the equation of edge.
 */
MEDCouplingUMesh *BuildMesh2DCutInternal(double eps, const MEDCouplingUMesh *splitMesh1D, const std::vector<int>& allEdges, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& allEdgesPtr, int offset,
                                         MCAuto<DataArrayInt>& idsLeftRight)
{
  int nbCellsInSplitMesh1D(splitMesh1D->getNumberOfCells());
  if(nbCellsInSplitMesh1D==0)
    throw INTERP_KERNEL::Exception("BuildMesh2DCutInternal : internal error ! input 1D mesh must have at least one cell !");
  const int *cSplitPtr(splitMesh1D->getNodalConnectivity()->begin()),*ciSplitPtr(splitMesh1D->getNodalConnectivityIndex()->begin());
  std::size_t nb(allEdges.size()),jj;
  if(nb%2!=0)
    throw INTERP_KERNEL::Exception("BuildMesh2DCutFrom : internal error 2 !");
  std::vector<int> edge1Bis(nb*2);
  std::vector< MCAuto<INTERP_KERNEL::Edge> > edge1BisPtr(nb*2);
  std::copy(allEdges.begin(),allEdges.end(),edge1Bis.begin());
  std::copy(allEdges.begin(),allEdges.end(),edge1Bis.begin()+nb);
  std::copy(allEdgesPtr.begin(),allEdgesPtr.end(),edge1BisPtr.begin());
  std::copy(allEdgesPtr.begin(),allEdgesPtr.end(),edge1BisPtr.begin()+nb);
  //
  idsLeftRight=DataArrayInt::New(); idsLeftRight->alloc(nbCellsInSplitMesh1D*2); idsLeftRight->fillWithValue(-2); idsLeftRight->rearrange(2);
  int *idsLeftRightPtr(idsLeftRight->getPointer());
  VectorOfCellInfo pool(edge1Bis,edge1BisPtr);
  for(int iStart=0;iStart<nbCellsInSplitMesh1D;)
    {// split [0:nbCellsInSplitMesh1D) in contiguous parts [iStart:iEnd)
      int iEnd(iStart);
      for(;iEnd<nbCellsInSplitMesh1D;)
        {
          for(jj=0;jj<nb && edge1Bis[2*jj+1]!=cSplitPtr[ciSplitPtr[iEnd]+2];jj++);
          if(jj!=nb)
            break;
          else
            iEnd++;
        }
      if(iEnd<nbCellsInSplitMesh1D)
        iEnd++;
      //
      MCAuto<MEDCouplingUMesh> partOfSplitMesh1D(static_cast<MEDCouplingUMesh *>(splitMesh1D->buildPartOfMySelfSlice(iStart,iEnd,1,true)));
      int pos(pool.getPositionOf(eps,partOfSplitMesh1D));
      //
      MCAuto<MEDCouplingUMesh>retTmp(MEDCouplingUMesh::New("",2));
      retTmp->setCoords(splitMesh1D->getCoords());
      retTmp->allocateCells();

      std::vector< std::vector<int> > out0;
      std::vector< std::vector< MCAuto<INTERP_KERNEL::Edge> > > out1;

      BuildMesh2DCutInternal2(partOfSplitMesh1D,pool.getConnOf(pos),pool.getEdgePtrOf(pos),out0,out1);
      for(std::size_t cnt=0;cnt<out0.size();cnt++)
        AddCellInMesh2D(retTmp,out0[cnt],out1[cnt]);
      pool.setMeshAt(pos,retTmp,iStart,iEnd,partOfSplitMesh1D,out0,out1);
      //
      iStart=iEnd;
    }
  for(int mm=0;mm<nbCellsInSplitMesh1D;mm++)
    pool.feedEdgeInfoAt(eps,mm,offset,idsLeftRightPtr+2*mm);
  return pool.getZeMesh().retn();
}

MEDCouplingUMesh *BuildMesh2DCutFrom(double eps, int cellIdInMesh2D, const MEDCouplingUMesh *mesh2DDesc, const MEDCouplingUMesh *splitMesh1D,
                                     const int *descBg, const int *descEnd, const std::vector< std::vector<int> >& intersectEdge1, int offset,
                                     MCAuto<DataArrayInt>& idsLeftRight)
{
  const int *cdescPtr(mesh2DDesc->getNodalConnectivity()->begin()),*cidescPtr(mesh2DDesc->getNodalConnectivityIndex()->begin());
  //
  std::vector<int> allEdges;
  std::vector< MCAuto<INTERP_KERNEL::Edge> > allEdgesPtr; // for each sub edge in splitMesh2D the uncut Edge object of the original mesh2D
  for(const int *it(descBg);it!=descEnd;it++) // for all edges in the descending connectivity of the 2D mesh in relative Fortran mode
    {
      int edgeId(std::abs(*it)-1);
      std::map< MCAuto<INTERP_KERNEL::Node>,int> m;
      MCAuto<INTERP_KERNEL::Edge> ee(MEDCouplingUMeshBuildQPFromEdge2((INTERP_KERNEL::NormalizedCellType)cdescPtr[cidescPtr[edgeId]],cdescPtr+cidescPtr[edgeId]+1,mesh2DDesc->getCoords()->begin(),m));
      const std::vector<int>& edge1(intersectEdge1[edgeId]);
      if(*it>0)
        allEdges.insert(allEdges.end(),edge1.begin(),edge1.end());
      else
        allEdges.insert(allEdges.end(),edge1.rbegin(),edge1.rend());
      std::size_t sz(edge1.size());
      for(std::size_t cnt=0;cnt<sz;cnt++)
        allEdgesPtr.push_back(ee);
    }
  //
  return BuildMesh2DCutInternal(eps,splitMesh1D,allEdges,allEdgesPtr,offset,idsLeftRight);
}

bool AreEdgeEqual(const double *coo2D, const INTERP_KERNEL::CellModel& typ1, const int *conn1, const INTERP_KERNEL::CellModel& typ2, const int *conn2, double eps)
{
  if(!typ1.isQuadratic() && !typ2.isQuadratic())
    {//easy case comparison not
      return conn1[0]==conn2[0] && conn1[1]==conn2[1];
    }
  else if(typ1.isQuadratic() && typ2.isQuadratic())
    {
      bool status0(conn1[0]==conn2[0] && conn1[1]==conn2[1]);
      if(!status0)
        return false;
      if(conn1[2]==conn2[2])
        return true;
      const double *a(coo2D+2*conn1[2]),*b(coo2D+2*conn2[2]);
      double dist(sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])));
      return dist<eps;
    }
  else
    {//only one is quadratic
      bool status0(conn1[0]==conn2[0] && conn1[1]==conn2[1]);
      if(!status0)
        return false;
      const double *a(0),*bb(0),*be(0);
      if(typ1.isQuadratic())
        {
          a=coo2D+2*conn1[2]; bb=coo2D+2*conn2[0]; be=coo2D+2*conn2[1];
        }
      else
        {
          a=coo2D+2*conn2[2]; bb=coo2D+2*conn1[0]; be=coo2D+2*conn1[1];
        }
      double b[2]; b[0]=(be[0]+bb[0])/2.; b[1]=(be[1]+bb[1])/2.;
      double dist(sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1])));
      return dist<eps;
    }
}

/*!
 * This method returns among the cellIds [ \a candidatesIn2DBg , \a candidatesIn2DEnd ) in \a mesh2DSplit those exactly sharing \a cellIdInMesh1DSplitRelative in \a mesh1DSplit.
 * \a mesh2DSplit and \a mesh1DSplit are expected to share the coordinates array.
 *
 * \param [in] cellIdInMesh1DSplitRelative is in Fortran mode using sign to specify direction.
 */
int FindRightCandidateAmong(const MEDCouplingUMesh *mesh2DSplit, const int *candidatesIn2DBg, const int *candidatesIn2DEnd, const MEDCouplingUMesh *mesh1DSplit, int cellIdInMesh1DSplitRelative, double eps)
{
  if(candidatesIn2DEnd==candidatesIn2DBg)
    throw INTERP_KERNEL::Exception("FindRightCandidateAmong : internal error 1 !");
  const double *coo(mesh2DSplit->getCoords()->begin());
  if(std::distance(candidatesIn2DBg,candidatesIn2DEnd)==1)
    return *candidatesIn2DBg;
  int edgeId(std::abs(cellIdInMesh1DSplitRelative)-1);
  MCAuto<MEDCouplingUMesh> cur1D(static_cast<MEDCouplingUMesh *>(mesh1DSplit->buildPartOfMySelf(&edgeId,&edgeId+1,true)));
  if(cellIdInMesh1DSplitRelative<0)
    cur1D->changeOrientationOfCells();
  const int *c1D(cur1D->getNodalConnectivity()->begin());
  const INTERP_KERNEL::CellModel& ref1DType(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)c1D[0]));
  for(const int *it=candidatesIn2DBg;it!=candidatesIn2DEnd;it++)
    {
      MCAuto<MEDCouplingUMesh> cur2D(static_cast<MEDCouplingUMesh *>(mesh2DSplit->buildPartOfMySelf(it,it+1,true)));
      const int *c(cur2D->getNodalConnectivity()->begin()),*ci(cur2D->getNodalConnectivityIndex()->begin());
      const INTERP_KERNEL::CellModel &cm(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)c[ci[0]]));
      unsigned sz(cm.getNumberOfSons2(c+ci[0]+1,ci[1]-ci[0]-1));
      INTERP_KERNEL::AutoPtr<int> tmpPtr(new int[ci[1]-ci[0]]);
      for(unsigned it2=0;it2<sz;it2++)
        {
          INTERP_KERNEL::NormalizedCellType typeOfSon;
          cm.fillSonCellNodalConnectivity2(it2,c+ci[0]+1,ci[1]-ci[0]-1,tmpPtr,typeOfSon);
          const INTERP_KERNEL::CellModel &curCM(INTERP_KERNEL::CellModel::GetCellModel(typeOfSon));
          if(AreEdgeEqual(coo,ref1DType,c1D+1,curCM,tmpPtr,eps))
            return *it;
        }
    }
  throw INTERP_KERNEL::Exception("FindRightCandidateAmong : internal error 2 ! Unable to find the edge among split cell !");
}

/*!
 * \param [out] intersectEdge1 - for each cell in \a m1Desc returns the result of the split. The result is given using pair of int given resp start and stop.
 *                               So for all edge \a i in \a m1Desc \a  intersectEdge1[i] is of length 2*n where n is the number of sub edges.
 *                               And for each j in [1,n) intersect[i][2*(j-1)+1]==intersect[i][2*j].
 * \param [out] subDiv2 - for each cell in \a m2Desc returns nodes that split it using convention \a m1Desc first, then \a m2Desc, then addCoo
 * \param [out] colinear2 - for each cell in \a m2Desc returns the edges in \a m1Desc that are colinear to it.
 * \param [out] addCoo - nodes to be append at the end
 * \param [out] mergedNodes - gives all pair of nodes of \a m2Desc that have same location than some nodes in \a m1Desc. key is id in \a m2Desc offseted and value is id in \a m1Desc.
 */
void MEDCouplingUMesh::Intersect1DMeshes(const MEDCouplingUMesh *m1Desc, const MEDCouplingUMesh *m2Desc, double eps,
                                         std::vector< std::vector<int> >& intersectEdge1, std::vector< std::vector<int> >& colinear2, std::vector< std::vector<int> >& subDiv2, std::vector<double>& addCoo, std::map<int,int>& mergedNodes)
{
  static const int SPACEDIM=2;
  INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision arcPrec(eps);
  const int *c1(m1Desc->getNodalConnectivity()->begin()),*ci1(m1Desc->getNodalConnectivityIndex()->begin());
  // Build BB tree of all edges in the tool mesh (second mesh)
  MCAuto<DataArrayDouble> bbox1Arr(m1Desc->getBoundingBoxForBBTree(eps)),bbox2Arr(m2Desc->getBoundingBoxForBBTree(eps));
  const double *bbox1(bbox1Arr->begin()),*bbox2(bbox2Arr->begin());
  int nDescCell1(m1Desc->getNumberOfCells()),nDescCell2(m2Desc->getNumberOfCells());
  intersectEdge1.resize(nDescCell1);
  colinear2.resize(nDescCell2);
  subDiv2.resize(nDescCell2);
  BBTree<SPACEDIM,int> myTree(bbox2,0,0,m2Desc->getNumberOfCells(),-eps);

  std::vector<int> candidates1(1);
  int offset1(m1Desc->getNumberOfNodes());
  int offset2(offset1+m2Desc->getNumberOfNodes());
  for(int i=0;i<nDescCell1;i++)  // for all edges in the first mesh
    {
      std::vector<int> candidates2; // edges of mesh2 candidate for intersection
      myTree.getIntersectingElems(bbox1+i*2*SPACEDIM,candidates2);
      if(!candidates2.empty()) // candidates2 holds edges from the second mesh potentially intersecting current edge i in mesh1
        {
          std::map<INTERP_KERNEL::Node *,int> map1,map2;
          // pol2 is not necessarily a closed polygon: just a set of (quadratic) edges (same as candidates2) in the Geometric DS format
          INTERP_KERNEL::QuadraticPolygon *pol2=MEDCouplingUMeshBuildQPFromMesh(m2Desc,candidates2,map2);
          candidates1[0]=i;
          INTERP_KERNEL::QuadraticPolygon *pol1=MEDCouplingUMeshBuildQPFromMesh(m1Desc,candidates1,map1);
          // This following part is to avoid that some removed nodes (for example due to a merge between pol1 and pol2) are replaced by a newly created one
          // This trick guarantees that Node * are discriminant (i.e. form a unique identifier)
          std::set<INTERP_KERNEL::Node *> nodes;
          pol1->getAllNodes(nodes); pol2->getAllNodes(nodes);
          std::size_t szz(nodes.size());
          std::vector< MCAuto<INTERP_KERNEL::Node> > nodesSafe(szz);
          std::set<INTERP_KERNEL::Node *>::const_iterator itt(nodes.begin());
          for(std::size_t iii=0;iii<szz;iii++,itt++)
            { (*itt)->incrRef(); nodesSafe[iii]=*itt; }
          // end of protection
          // Performs egde cutting:
          pol1->splitAbs(*pol2,map1,map2,offset1,offset2,candidates2,intersectEdge1[i],i,colinear2,subDiv2,addCoo,mergedNodes);
          delete pol2;
          delete pol1;
        }
      else
        // Copy the edge (take only the two first points, ie discard quadratic point at this stage)
        intersectEdge1[i].insert(intersectEdge1[i].end(),c1+ci1[i]+1,c1+ci1[i]+3);
    }
}


/*!
 * This method is private and is the first step of Partition of 2D mesh (spaceDim==2 and meshDim==2).
 * It builds the descending connectivity of the two meshes, and then using a binary tree
 * it computes the edge intersections. This results in new points being created : they're stored in addCoo.
 * Documentation about parameters  colinear2 and subDiv2 can be found in method QuadraticPolygon::splitAbs().
 */
void MEDCouplingUMesh::IntersectDescending2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps,
                                                   std::vector< std::vector<int> >& intersectEdge1, std::vector< std::vector<int> >& colinear2, std::vector< std::vector<int> >& subDiv2,
                                                   MEDCouplingUMesh *& m1Desc, DataArrayInt *&desc1, DataArrayInt *&descIndx1, DataArrayInt *&revDesc1, DataArrayInt *&revDescIndx1,
                                                   std::vector<double>& addCoo,
                                                   MEDCouplingUMesh *& m2Desc, DataArrayInt *&desc2, DataArrayInt *&descIndx2, DataArrayInt *&revDesc2, DataArrayInt *&revDescIndx2)
{
  // Build desc connectivity
  desc1=DataArrayInt::New(); descIndx1=DataArrayInt::New(); revDesc1=DataArrayInt::New(); revDescIndx1=DataArrayInt::New();
  desc2=DataArrayInt::New();
  descIndx2=DataArrayInt::New();
  revDesc2=DataArrayInt::New();
  revDescIndx2=DataArrayInt::New();
  MCAuto<DataArrayInt> dd1(desc1),dd2(descIndx1),dd3(revDesc1),dd4(revDescIndx1);
  MCAuto<DataArrayInt> dd5(desc2),dd6(descIndx2),dd7(revDesc2),dd8(revDescIndx2);
  m1Desc=m1->buildDescendingConnectivity2(desc1,descIndx1,revDesc1,revDescIndx1);
  m2Desc=m2->buildDescendingConnectivity2(desc2,descIndx2,revDesc2,revDescIndx2);
  MCAuto<MEDCouplingUMesh> dd9(m1Desc),dd10(m2Desc);
  std::map<int,int> notUsedMap;
  Intersect1DMeshes(m1Desc,m2Desc,eps,intersectEdge1,colinear2,subDiv2,addCoo,notUsedMap);
  m1Desc->incrRef(); desc1->incrRef(); descIndx1->incrRef(); revDesc1->incrRef(); revDescIndx1->incrRef();
  m2Desc->incrRef(); desc2->incrRef(); descIndx2->incrRef(); revDesc2->incrRef(); revDescIndx2->incrRef();
}

/**
 * Private. Third step of the partitioning algorithm (Intersect2DMeshes): reconstruct full 2D cells from the
 * (newly created) nodes corresponding to the edge intersections.
 * Output params:
 * @param[out] cr, crI connectivity of the resulting mesh
 * @param[out] cNb1, cNb2 correspondance arrays giving for the merged mesh the initial cells IDs in m1 / m2
 * TODO: describe input parameters
 */
void MEDCouplingUMesh::BuildIntersecting2DCellsFromEdges(double eps, const MEDCouplingUMesh *m1, const int *desc1, const int *descIndx1,
                                                         const std::vector<std::vector<int> >& intesctEdges1, const std::vector< std::vector<int> >& colinear2,
                                                         const MEDCouplingUMesh *m2, const int *desc2, const int *descIndx2, const std::vector<std::vector<int> >& intesctEdges2,
                                                         const std::vector<double>& addCoords,
                                                         std::vector<double>& addCoordsQuadratic, std::vector<int>& cr, std::vector<int>& crI, std::vector<int>& cNb1, std::vector<int>& cNb2)
{
  static const int SPACEDIM=2;
  const double *coo1(m1->getCoords()->begin());
  const int *conn1(m1->getNodalConnectivity()->begin()),*connI1(m1->getNodalConnectivityIndex()->begin());
  int offset1(m1->getNumberOfNodes());
  const double *coo2(m2->getCoords()->begin());
  const int *conn2(m2->getNodalConnectivity()->begin()),*connI2(m2->getNodalConnectivityIndex()->begin());
  int offset2(offset1+m2->getNumberOfNodes());
  int offset3(offset2+((int)addCoords.size())/2);
  MCAuto<DataArrayDouble> bbox1Arr(m1->getBoundingBoxForBBTree(eps)),bbox2Arr(m2->getBoundingBoxForBBTree(eps));
  const double *bbox1(bbox1Arr->begin()),*bbox2(bbox2Arr->begin());
  // Here a BBTree on 2D-cells, not on segments:
  BBTree<SPACEDIM,int> myTree(bbox2,0,0,m2->getNumberOfCells(),eps);
  int ncell1(m1->getNumberOfCells());
  crI.push_back(0);
  for(int i=0;i<ncell1;i++)
    {
      std::vector<int> candidates2;
      myTree.getIntersectingElems(bbox1+i*2*SPACEDIM,candidates2);
      std::map<INTERP_KERNEL::Node *,int> mapp;
      std::map<int,INTERP_KERNEL::Node *> mappRev;
      INTERP_KERNEL::QuadraticPolygon pol1;
      INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)conn1[connI1[i]];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
      // Populate mapp and mappRev with nodes from the current cell (i) from mesh1 - this also builds the Node* objects:
      MEDCouplingUMeshBuildQPFromMesh3(coo1,offset1,coo2,offset2,addCoords,desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1,/* output */mapp,mappRev);
      // pol1 is the full cell from mesh2, in QP format, with all the additional intersecting nodes.
      pol1.buildFromCrudeDataArray(mappRev,cm.isQuadratic(),conn1+connI1[i]+1,coo1,
          desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1);
      //
      std::set<INTERP_KERNEL::Edge *> edges1;// store all edges of pol1 that are NOT consumed by intersect cells. If any after iteration over candidates2 -> a part of pol1 should appear in result
      std::set<INTERP_KERNEL::Edge *> edgesBoundary2;// store all edges that are on boundary of (pol2 intersect pol1) minus edges on pol1.
      INTERP_KERNEL::IteratorOnComposedEdge it1(&pol1);
      for(it1.first();!it1.finished();it1.next())
        edges1.insert(it1.current()->getPtr());
      //
      std::map<int,std::vector<INTERP_KERNEL::ElementaryEdge *> > edgesIn2ForShare; // common edges
      std::vector<INTERP_KERNEL::QuadraticPolygon> pol2s(candidates2.size());
      int ii=0;
      for(std::vector<int>::const_iterator it2=candidates2.begin();it2!=candidates2.end();it2++,ii++)
        {
          INTERP_KERNEL::NormalizedCellType typ2=(INTERP_KERNEL::NormalizedCellType)conn2[connI2[*it2]];
          const INTERP_KERNEL::CellModel& cm2=INTERP_KERNEL::CellModel::GetCellModel(typ2);
          // Complete mapping with elements coming from the current cell it2 in mesh2:
          MEDCouplingUMeshBuildQPFromMesh3(coo1,offset1,coo2,offset2,addCoords,desc2+descIndx2[*it2],desc2+descIndx2[*it2+1],intesctEdges2,/* output */mapp,mappRev);
          // pol2 is the new QP in the final merged result.
          pol2s[ii].buildFromCrudeDataArray2(mappRev,cm2.isQuadratic(),conn2+connI2[*it2]+1,coo2,desc2+descIndx2[*it2],desc2+descIndx2[*it2+1],intesctEdges2,
              pol1,desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1,colinear2, /* output */ edgesIn2ForShare);
        }
      ii=0;
      for(std::vector<int>::const_iterator it2=candidates2.begin();it2!=candidates2.end();it2++,ii++)
        {
          INTERP_KERNEL::ComposedEdge::InitLocationsWithOther(pol1,pol2s[ii]);
          pol2s[ii].updateLocOfEdgeFromCrudeDataArray2(desc2+descIndx2[*it2],desc2+descIndx2[*it2+1],intesctEdges2,pol1,desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1,colinear2);
          //MEDCouplingUMeshAssignOnLoc(pol1,pol2,desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1,desc2+descIndx2[*it2],desc2+descIndx2[*it2+1],intesctEdges2,colinear2);
          pol1.buildPartitionsAbs(pol2s[ii],edges1,edgesBoundary2,mapp,i,*it2,offset3,addCoordsQuadratic,cr,crI,cNb1,cNb2);
        }
      // Deals with remaining (non-consumed) edges from m1: these are the edges that were never touched
      // by m2 but that we still want to keep in the final result.
      if(!edges1.empty())
        {
          try
          {
              INTERP_KERNEL::QuadraticPolygon::ComputeResidual(pol1,edges1,edgesBoundary2,mapp,offset3,i,addCoordsQuadratic,cr,crI,cNb1,cNb2);
          }
          catch(INTERP_KERNEL::Exception& e)
          {
              std::ostringstream oss; oss << "Error when computing residual of cell #" << i << " in source/m1 mesh ! Maybe the neighbours of this cell in mesh are not well connected !\n" << "The deep reason is the following : " << e.what();
              throw INTERP_KERNEL::Exception(oss.str());
          }
        }
      for(std::map<int,INTERP_KERNEL::Node *>::const_iterator it=mappRev.begin();it!=mappRev.end();it++)
        (*it).second->decrRef();
    }
}

void InsertNodeInConnIfNecessary(int nodeIdToInsert, std::vector<int>& conn, const double *coords, double eps)
{
  std::vector<int>::iterator it(std::find(conn.begin(),conn.end(),nodeIdToInsert));
  if(it!=conn.end())
    return ;
  std::size_t sz(conn.size());
  std::size_t found(std::numeric_limits<std::size_t>::max());
  for(std::size_t i=0;i<sz;i++)
    {
      int pt0(conn[i]),pt1(conn[(i+1)%sz]);
      double v1[3]={coords[3*pt1+0]-coords[3*pt0+0],coords[3*pt1+1]-coords[3*pt0+1],coords[3*pt1+2]-coords[3*pt0+2]},v2[3]={coords[3*nodeIdToInsert+0]-coords[3*pt0+0],coords[3*nodeIdToInsert+1]-coords[3*pt0+1],coords[3*nodeIdToInsert+2]-coords[3*pt0+2]};
      double normm(sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]));
      std::transform(v1,v1+3,v1,std::bind2nd(std::multiplies<double>(),1./normm));
      std::transform(v2,v2+3,v2,std::bind2nd(std::multiplies<double>(),1./normm));
      double v3[3];
      v3[0]=v1[1]*v2[2]-v1[2]*v2[1]; v3[1]=v1[2]*v2[0]-v1[0]*v2[2]; v3[2]=v1[0]*v2[1]-v1[1]*v2[0];
      double normm2(sqrt(v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2])),dotTest(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
      if(normm2<eps)
        if(dotTest>eps && dotTest<1.-eps)
          {
            found=i;
            break;
          }
    }
  if(found==std::numeric_limits<std::size_t>::max())
    throw INTERP_KERNEL::Exception("InsertNodeInConnIfNecessary : not found point !");
  conn.insert(conn.begin()+(found+1)%sz,nodeIdToInsert);
}

void SplitIntoToPart(const std::vector<int>& conn, int pt0, int pt1, std::vector<int>& part0, std::vector<int>& part1)
{
  std::size_t sz(conn.size());
  std::vector<int> *curPart(&part0);
  for(std::size_t i=0;i<sz;i++)
    {
      int nextt(conn[(i+1)%sz]);
      (*curPart).push_back(nextt);
      if(nextt==pt0 || nextt==pt1)
        {
          if(curPart==&part0)
            curPart=&part1;
          else
            curPart=&part0;
          (*curPart).push_back(nextt);
        }
    }
}

/*!
 * this method method splits cur cells 3D Surf in sub cells 3DSurf using the previous subsplit. This method is the last one used to clip.
 */
void MEDCouplingUMesh::buildSubCellsFromCut(const std::vector< std::pair<int,int> >& cut3DSurf,
                                            const int *desc, const int *descIndx, const double *coords, double eps,
                                            std::vector<std::vector<int> >& res) const
{
  checkFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSubCellsFromCut works on umeshes with meshdim equal to 3 and spaceDim equal to 3 too!");
  const int *nodal3D(_nodal_connec->begin()),*nodalIndx3D(_nodal_connec_index->begin());
  int nbOfCells(getNumberOfCells());
  if(nbOfCells!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSubCellsFromCut works only with single cell presently !");
  for(int i=0;i<nbOfCells;i++)
    {
      int offset(descIndx[i]),nbOfFaces(descIndx[i+1]-offset);
      for(int j=0;j<nbOfFaces;j++)
        {
          const std::pair<int,int>& p=cut3DSurf[desc[offset+j]];
          const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)nodal3D[nodalIndx3D[i]]));
          int sz=nodalIndx3D[i+1]-nodalIndx3D[i]-1;
          INTERP_KERNEL::AutoPtr<int> tmp(new int[sz]);
          INTERP_KERNEL::NormalizedCellType cmsId;
          unsigned nbOfNodesSon(cm.fillSonCellNodalConnectivity2(j,nodal3D+nodalIndx3D[i]+1,sz,tmp,cmsId));
          std::vector<int> elt((int *)tmp,(int *)tmp+nbOfNodesSon);
          if(p.first!=-1 && p.second!=-1)
            {
              if(p.first!=-2)
                {
                  InsertNodeInConnIfNecessary(p.first,elt,coords,eps);
                  InsertNodeInConnIfNecessary(p.second,elt,coords,eps);
                  std::vector<int> elt1,elt2;
                  SplitIntoToPart(elt,p.first,p.second,elt1,elt2);
                  res.push_back(elt1);
                  res.push_back(elt2);
                }
              else
                res.push_back(elt);
            }
          else
            res.push_back(elt);
        }
    }
}

/*!
 * It is the linear part of MEDCouplingUMesh::split2DCells. Here no additionnal nodes will be added in \b this. So coordinates pointer remain unchanged (is not even touch).
 *
 * \sa MEDCouplingUMesh::split2DCells
 */
void MEDCouplingUMesh::split2DCellsLinear(const DataArrayInt *desc, const DataArrayInt *descI, const DataArrayInt *subNodesInSeg, const DataArrayInt *subNodesInSegI)
{
  checkConnectivityFullyDefined();
  int ncells(getNumberOfCells()),lgthToReach(getNodalConnectivityArrayLen()+subNodesInSeg->getNumberOfTuples());
  MCAuto<DataArrayInt> c(DataArrayInt::New()); c->alloc((std::size_t)lgthToReach);
  const int *subPtr(subNodesInSeg->begin()),*subIPtr(subNodesInSegI->begin()),*descPtr(desc->begin()),*descIPtr(descI->begin()),*oldConn(getNodalConnectivity()->begin());
  int *cPtr(c->getPointer()),*ciPtr(getNodalConnectivityIndex()->getPointer());
  int prevPosOfCi(ciPtr[0]);
  for(int i=0;i<ncells;i++,ciPtr++,descIPtr++)
    {
      int offset(descIPtr[0]),sz(descIPtr[1]-descIPtr[0]),deltaSz(0);
      *cPtr++=(int)INTERP_KERNEL::NORM_POLYGON; *cPtr++=oldConn[prevPosOfCi+1];
      for(int j=0;j<sz;j++)
        {
          int offset2(subIPtr[descPtr[offset+j]]),sz2(subIPtr[descPtr[offset+j]+1]-subIPtr[descPtr[offset+j]]);
          for(int k=0;k<sz2;k++)
            *cPtr++=subPtr[offset2+k];
          if(j!=sz-1)
            *cPtr++=oldConn[prevPosOfCi+j+2];
          deltaSz+=sz2;
        }
      prevPosOfCi=ciPtr[1];
      ciPtr[1]=ciPtr[0]+1+sz+deltaSz;//sz==old nb of nodes because (nb of subedges=nb of nodes for polygons)
    }
  if(c->end()!=cPtr)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::split2DCellsLinear : Some of edges to be split are orphan !");
  _nodal_connec->decrRef();
  _nodal_connec=c.retn(); _types.clear(); _types.insert(INTERP_KERNEL::NORM_POLYGON);
}


/*!
 * It is the quadratic part of MEDCouplingUMesh::split2DCells. Here some additionnal nodes can be added at the end of coordinates array object.
 *
 * \return  int - the number of new nodes created.
 * \sa MEDCouplingUMesh::split2DCells
 */
int MEDCouplingUMesh::split2DCellsQuadratic(const DataArrayInt *desc, const DataArrayInt *descI, const DataArrayInt *subNodesInSeg, const DataArrayInt *subNodesInSegI, const DataArrayInt *mid, const DataArrayInt *midI)
{
  checkConsistencyLight();
  int ncells(getNumberOfCells()),lgthToReach(getNodalConnectivityArrayLen()+2*subNodesInSeg->getNumberOfTuples()),nodesCnt(getNumberOfNodes());
  MCAuto<DataArrayInt> c(DataArrayInt::New()); c->alloc((std::size_t)lgthToReach);
  MCAuto<DataArrayDouble> addCoo(DataArrayDouble::New()); addCoo->alloc(0,1);
  const int *subPtr(subNodesInSeg->begin()),*subIPtr(subNodesInSegI->begin()),*descPtr(desc->begin()),*descIPtr(descI->begin()),*oldConn(getNodalConnectivity()->begin());
  const int *midPtr(mid->begin()),*midIPtr(midI->begin());
  const double *oldCoordsPtr(getCoords()->begin());
  int *cPtr(c->getPointer()),*ciPtr(getNodalConnectivityIndex()->getPointer());
  int prevPosOfCi(ciPtr[0]);
  for(int i=0;i<ncells;i++,ciPtr++,descIPtr++)
    {
      int offset(descIPtr[0]),sz(descIPtr[1]-descIPtr[0]),deltaSz(sz);
      for(int j=0;j<sz;j++)
        { int sz2(subIPtr[descPtr[offset+j]+1]-subIPtr[descPtr[offset+j]]); deltaSz+=sz2; }
      *cPtr++=(int)INTERP_KERNEL::NORM_QPOLYG; cPtr[0]=oldConn[prevPosOfCi+1];
      for(int j=0;j<sz;j++)//loop over subedges of oldConn
        {
          int offset2(subIPtr[descPtr[offset+j]]),sz2(subIPtr[descPtr[offset+j]+1]-subIPtr[descPtr[offset+j]]),offset3(midIPtr[descPtr[offset+j]]);
          if(sz2==0)
            {
              if(j<sz-1)
                cPtr[1]=oldConn[prevPosOfCi+2+j];
              cPtr[deltaSz]=oldConn[prevPosOfCi+1+j+sz]; cPtr++;
              continue;
            }
          std::vector<INTERP_KERNEL::Node *> ns(3);
          ns[0]=new INTERP_KERNEL::Node(oldCoordsPtr[2*oldConn[prevPosOfCi+1+j]],oldCoordsPtr[2*oldConn[prevPosOfCi+1+j]+1]);
          ns[1]=new INTERP_KERNEL::Node(oldCoordsPtr[2*oldConn[prevPosOfCi+1+(1+j)%sz]],oldCoordsPtr[2*oldConn[prevPosOfCi+1+(1+j)%sz]+1]);
          ns[2]=new INTERP_KERNEL::Node(oldCoordsPtr[2*oldConn[prevPosOfCi+1+sz+j]],oldCoordsPtr[2*oldConn[prevPosOfCi+1+sz+j]+1]);
          MCAuto<INTERP_KERNEL::Edge> e(INTERP_KERNEL::QuadraticPolygon::BuildArcCircleEdge(ns));
          for(int k=0;k<sz2;k++)//loop over subsplit of current subedge
            {
              cPtr[1]=subPtr[offset2+k];
              cPtr[deltaSz]=InternalAddPoint(e,midPtr[offset3+k],oldCoordsPtr,cPtr[0],cPtr[1],*addCoo,nodesCnt); cPtr++;
            }
          int tmpEnd(oldConn[prevPosOfCi+1+(j+1)%sz]);
          if(j!=sz-1)
            { cPtr[1]=tmpEnd; }
          cPtr[deltaSz]=InternalAddPoint(e,midPtr[offset3+sz2],oldCoordsPtr,cPtr[0],tmpEnd,*addCoo,nodesCnt); cPtr++;
        }
      prevPosOfCi=ciPtr[1]; cPtr+=deltaSz;
      ciPtr[1]=ciPtr[0]+1+2*deltaSz;//sz==old nb of nodes because (nb of subedges=nb of nodes for polygons)
    }
  if(c->end()!=cPtr)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::split2DCellsQuadratic : Some of edges to be split are orphan !");
  _nodal_connec->decrRef();
  _nodal_connec=c.retn(); _types.clear(); _types.insert(INTERP_KERNEL::NORM_QPOLYG);
  addCoo->rearrange(2);
  MCAuto<DataArrayDouble> coo(DataArrayDouble::Aggregate(getCoords(),addCoo));//info are copied from getCoords() by using Aggregate
  setCoords(coo);
  return addCoo->getNumberOfTuples();
}


/// @endcond

/*!
 * Partitions the first given 2D mesh using the second given 2D mesh as a tool, and
 * returns a result mesh constituted by polygons.
 * Thus the final result contains all nodes from m1 plus new nodes. However it doesn't necessarily contains
 * all nodes from m2.
 * The meshes should be in 2D space. In
 * addition, returns two arrays mapping cells of the result mesh to cells of the input
 * meshes.
 *  \param [in] m1 - the first input mesh which is a partitioned object. The mesh must be so that each point in the space covered by \a m1
 *                      must be covered exactly by one entity, \b no \b more. If it is not the case, some tools are available to heal the mesh (conformize2D, mergeNodes)
 *  \param [in] m2 - the second input mesh which is a partition tool. The mesh must be so that each point in the space covered by \a m2
 *                      must be covered exactly by one entity, \b no \b more. If it is not the case, some tools are available to heal the mesh (conformize2D, mergeNodes)
 *  \param [in] eps - precision used to detect coincident mesh entities.
 *  \param [out] cellNb1 - a new instance of DataArrayInt holding for each result
 *         cell an id of the cell of \a m1 it comes from. The caller is to delete
 *         this array using decrRef() as it is no more needed.
 *  \param [out] cellNb2 - a new instance of DataArrayInt holding for each result
 *         cell an id of the cell of \a m2 it comes from. -1 value means that a
 *         result cell comes from a cell (or part of cell) of \a m1 not overlapped by
 *         any cell of \a m2. The caller is to delete this array using decrRef() as
 *         it is no more needed.
 *  \return MEDCouplingUMesh * - the result 2D mesh which is a new instance of
 *         MEDCouplingUMesh. The caller is to delete this mesh using decrRef() as it
 *         is no more needed.
 *  \throw If the coordinates array is not set in any of the meshes.
 *  \throw If the nodal connectivity of cells is not defined in any of the meshes.
 *  \throw If any of the meshes is not a 2D mesh in 2D space.
 *
 *  \sa conformize2D, mergeNodes
 */
MEDCouplingUMesh *MEDCouplingUMesh::Intersect2DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2,
                                                      double eps, DataArrayInt *&cellNb1, DataArrayInt *&cellNb2)
{
  if(!m1 || !m2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Intersect2DMeshes : input meshes must be not NULL !");
  m1->checkFullyDefined();
  m2->checkFullyDefined();
  INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision arcPrec(eps);
  if(m1->getMeshDimension()!=2 || m1->getSpaceDimension()!=2 || m2->getMeshDimension()!=2 || m2->getSpaceDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Intersect2DMeshes works on umeshes m1 AND m2  with meshdim equal to 2 and spaceDim equal to 2 too!");

  // Step 1: compute all edge intersections (new nodes)
  std::vector< std::vector<int> > intersectEdge1, colinear2, subDiv2;
  MEDCouplingUMesh *m1Desc=0,*m2Desc=0; // descending connec. meshes
  DataArrayInt *desc1=0,*descIndx1=0,*revDesc1=0,*revDescIndx1=0,*desc2=0,*descIndx2=0,*revDesc2=0,*revDescIndx2=0;
  std::vector<double> addCoo,addCoordsQuadratic;  // coordinates of newly created nodes
  IntersectDescending2DMeshes(m1,m2,eps,intersectEdge1,colinear2, subDiv2,
                              m1Desc,desc1,descIndx1,revDesc1,revDescIndx1,
                              addCoo, m2Desc,desc2,descIndx2,revDesc2,revDescIndx2);
  revDesc1->decrRef(); revDescIndx1->decrRef(); revDesc2->decrRef(); revDescIndx2->decrRef();
  MCAuto<DataArrayInt> dd1(desc1),dd2(descIndx1),dd3(desc2),dd4(descIndx2);
  MCAuto<MEDCouplingUMesh> dd5(m1Desc),dd6(m2Desc);

  // Step 2: re-order newly created nodes according to the ordering found in m2
  std::vector< std::vector<int> > intersectEdge2;
  BuildIntersectEdges(m1Desc,m2Desc,addCoo,subDiv2,intersectEdge2);
  subDiv2.clear(); dd5=0; dd6=0;

  // Step 3:
  std::vector<int> cr,crI; //no DataArrayInt because interface with Geometric2D
  std::vector<int> cNb1,cNb2; //no DataArrayInt because interface with Geometric2D
  BuildIntersecting2DCellsFromEdges(eps,m1,desc1->begin(),descIndx1->begin(),intersectEdge1,colinear2,m2,desc2->begin(),descIndx2->begin(),intersectEdge2,addCoo,
                                    /* outputs -> */addCoordsQuadratic,cr,crI,cNb1,cNb2);

  // Step 4: Prepare final result:
  MCAuto<DataArrayDouble> addCooDa(DataArrayDouble::New());
  addCooDa->alloc((int)(addCoo.size())/2,2);
  std::copy(addCoo.begin(),addCoo.end(),addCooDa->getPointer());
  MCAuto<DataArrayDouble> addCoordsQuadraticDa(DataArrayDouble::New());
  addCoordsQuadraticDa->alloc((int)(addCoordsQuadratic.size())/2,2);
  std::copy(addCoordsQuadratic.begin(),addCoordsQuadratic.end(),addCoordsQuadraticDa->getPointer());
  std::vector<const DataArrayDouble *> coordss(4);
  coordss[0]=m1->getCoords(); coordss[1]=m2->getCoords(); coordss[2]=addCooDa; coordss[3]=addCoordsQuadraticDa;
  MCAuto<DataArrayDouble> coo(DataArrayDouble::Aggregate(coordss));
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::New("Intersect2D",2));
  MCAuto<DataArrayInt> conn(DataArrayInt::New()); conn->alloc((int)cr.size(),1); std::copy(cr.begin(),cr.end(),conn->getPointer());
  MCAuto<DataArrayInt> connI(DataArrayInt::New()); connI->alloc((int)crI.size(),1); std::copy(crI.begin(),crI.end(),connI->getPointer());
  MCAuto<DataArrayInt> c1(DataArrayInt::New()); c1->alloc((int)cNb1.size(),1); std::copy(cNb1.begin(),cNb1.end(),c1->getPointer());
  MCAuto<DataArrayInt> c2(DataArrayInt::New()); c2->alloc((int)cNb2.size(),1); std::copy(cNb2.begin(),cNb2.end(),c2->getPointer());
  ret->setConnectivity(conn,connI,true);
  ret->setCoords(coo);
  cellNb1=c1.retn(); cellNb2=c2.retn();
  return ret.retn();
}

/*!
 * Partitions the first given 2D mesh using the second given 1D mesh as a tool.
 * Thus the final result contains the aggregation of nodes of \a mesh2D, then nodes of \a mesh1D, then new nodes that are the result of the intersection
 * and finaly, in case of quadratic polygon the centers of edges new nodes.
 * The meshes should be in 2D space. In addition, returns two arrays mapping cells of the resulting mesh to cells of the input.
 *
 * \param [in] mesh2D - the 2D mesh (spacedim=meshdim=2) to be intersected using \a mesh1D tool. The mesh must be so that each point in the space covered by \a mesh2D
 *                      must be covered exactly by one entity, \b no \b more. If it is not the case, some tools are available to heal the mesh (conformize2D, mergeNodes)
 * \param [in] mesh1D - the 1D mesh (spacedim=2 meshdim=1) the is the tool that will be used to intersect \a mesh2D. \a mesh1D must be ordered consecutively. If it is not the case
 *                      you can invoke orderConsecutiveCells1D on \a mesh1D.
 * \param [in] eps - precision used to perform intersections and localization operations.
 * \param [out] splitMesh2D - the result of the split of \a mesh2D mesh.
 * \param [out] splitMesh1D - the result of the split of \a mesh1D mesh.
 * \param [out] cellIdInMesh2D - the array that gives for each cell id \a i in \a splitMesh2D the id in \a mesh2D it comes from.
 *                               So this array has a number of tuples equal to the number of cells of \a splitMesh2D and a number of component equal to 1.
 * \param [out] cellIdInMesh1D - the array of pair that gives for each cell id \a i in \a splitMesh1D the cell in \a splitMesh2D on the left for the 1st component
 *                               and the cell in \a splitMesh2D on the right for the 2nt component. -1 means no cell.
 *                               So this array has a number of tuples equal to the number of cells of \a splitMesh1D and a number of components equal to 2.
 *
 * \sa Intersect2DMeshes, orderConsecutiveCells1D, conformize2D, mergeNodes
 */
void MEDCouplingUMesh::Intersect2DMeshWith1DLine(const MEDCouplingUMesh *mesh2D, const MEDCouplingUMesh *mesh1D, double eps, MEDCouplingUMesh *&splitMesh2D, MEDCouplingUMesh *&splitMesh1D, DataArrayInt *&cellIdInMesh2D, DataArrayInt *&cellIdInMesh1D)
{
  if(!mesh2D || !mesh1D)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Intersect2DMeshWith1DLine : input meshes must be not NULL !");
  mesh2D->checkFullyDefined();
  mesh1D->checkFullyDefined();
  const std::vector<std::string>& compNames(mesh2D->getCoords()->getInfoOnComponents());
  if(mesh2D->getMeshDimension()!=2 || mesh2D->getSpaceDimension()!=2 || mesh1D->getMeshDimension()!=1 || mesh1D->getSpaceDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Intersect2DMeshWith1DLine works with mesh2D with spacedim=meshdim=2 and mesh1D with meshdim=1 spaceDim=2 !");
  // Step 1: compute all edge intersections (new nodes)
  std::vector< std::vector<int> > intersectEdge1, colinear2, subDiv2;
  std::vector<double> addCoo,addCoordsQuadratic;  // coordinates of newly created nodes
  INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision arcPrec(eps);
  //
  // Build desc connectivity
  DataArrayInt *desc1(DataArrayInt::New()),*descIndx1(DataArrayInt::New()),*revDesc1(DataArrayInt::New()),*revDescIndx1(DataArrayInt::New());
  MCAuto<DataArrayInt> dd1(desc1),dd2(descIndx1),dd3(revDesc1),dd4(revDescIndx1);
  MCAuto<MEDCouplingUMesh> m1Desc(mesh2D->buildDescendingConnectivity2(desc1,descIndx1,revDesc1,revDescIndx1));
  std::map<int,int> mergedNodes;
  Intersect1DMeshes(m1Desc,mesh1D,eps,intersectEdge1,colinear2,subDiv2,addCoo,mergedNodes);
  // use mergeNodes to fix intersectEdge1
  for(std::vector< std::vector<int> >::iterator it0=intersectEdge1.begin();it0!=intersectEdge1.end();it0++)
    {
      std::size_t n((*it0).size()/2);
      int eltStart((*it0)[0]),eltEnd((*it0)[2*n-1]);
      std::map<int,int>::const_iterator it1;
      it1=mergedNodes.find(eltStart);
      if(it1!=mergedNodes.end())
        (*it0)[0]=(*it1).second;
      it1=mergedNodes.find(eltEnd);
      if(it1!=mergedNodes.end())
        (*it0)[2*n-1]=(*it1).second;
    }
  //
  MCAuto<DataArrayDouble> addCooDa(DataArrayDouble::New());
  addCooDa->useArray(&addCoo[0],false,C_DEALLOC,(int)addCoo.size()/2,2);
  // Step 2: re-order newly created nodes according to the ordering found in m2
  std::vector< std::vector<int> > intersectEdge2;
  BuildIntersectEdges(m1Desc,mesh1D,addCoo,subDiv2,intersectEdge2);
  subDiv2.clear();
  // Step 3: compute splitMesh1D
  MCAuto<DataArrayInt> idsInRet1Colinear,idsInDescMesh2DForIdsInRetColinear;
  MCAuto<DataArrayInt> ret2(DataArrayInt::New()); ret2->alloc(0,1);
  MCAuto<MEDCouplingUMesh> ret1(BuildMesh1DCutFrom(mesh1D,intersectEdge2,mesh2D->getCoords(),addCoo,mergedNodes,colinear2,intersectEdge1,
      idsInRet1Colinear,idsInDescMesh2DForIdsInRetColinear));
  MCAuto<DataArrayInt> ret3(DataArrayInt::New()); ret3->alloc(ret1->getNumberOfCells()*2,1); ret3->fillWithValue(std::numeric_limits<int>::max()); ret3->rearrange(2);
  MCAuto<DataArrayInt> idsInRet1NotColinear(idsInRet1Colinear->buildComplement(ret1->getNumberOfCells()));
  // deal with cells in mesh2D that are not cut but only some of their edges are
  MCAuto<DataArrayInt> idsInDesc2DToBeRefined(idsInDescMesh2DForIdsInRetColinear->deepCopy());
  idsInDesc2DToBeRefined->abs(); idsInDesc2DToBeRefined->applyLin(1,-1);
  idsInDesc2DToBeRefined=idsInDesc2DToBeRefined->buildUnique();
  MCAuto<DataArrayInt> out0s;//ids in mesh2D that are impacted by the fact that some edges of \a mesh1D are part of the edges of those cells
  if(!idsInDesc2DToBeRefined->empty())
    {
      DataArrayInt *out0(0),*outi0(0);
      MEDCouplingUMesh::ExtractFromIndexedArrays(idsInDesc2DToBeRefined->begin(),idsInDesc2DToBeRefined->end(),dd3,dd4,out0,outi0);
      MCAuto<DataArrayInt> outi0s(outi0);
      out0s=out0;
      out0s=out0s->buildUnique();
      out0s->sort(true);
    }
  //
  MCAuto<MEDCouplingUMesh> ret1NonCol(static_cast<MEDCouplingUMesh *>(ret1->buildPartOfMySelf(idsInRet1NotColinear->begin(),idsInRet1NotColinear->end())));
  MCAuto<DataArrayDouble> baryRet1(ret1NonCol->computeCellCenterOfMass());
  MCAuto<DataArrayInt> elts,eltsIndex;
  mesh2D->getCellsContainingPoints(baryRet1->begin(),baryRet1->getNumberOfTuples(),eps,elts,eltsIndex);
  MCAuto<DataArrayInt> eltsIndex2(DataArrayInt::New()); eltsIndex2->alloc(0,1);
  if (eltsIndex->getNumberOfTuples() > 1)
    eltsIndex2 = eltsIndex->deltaShiftIndex();
  MCAuto<DataArrayInt> eltsIndex3(eltsIndex2->findIdsEqual(1));
  if(eltsIndex2->count(0)+eltsIndex3->getNumberOfTuples()!=ret1NonCol->getNumberOfCells())
    throw INTERP_KERNEL::Exception("Intersect2DMeshWith1DLine : internal error 1 !");
  MCAuto<DataArrayInt> cellsToBeModified(elts->buildUnique());
  MCAuto<DataArrayInt> untouchedCells(cellsToBeModified->buildComplement(mesh2D->getNumberOfCells()));
  if((DataArrayInt *)out0s)
    untouchedCells=untouchedCells->buildSubstraction(out0s);//if some edges in ret1 are colinear to descending mesh of mesh2D remove cells from untouched one
  std::vector< MCAuto<MEDCouplingUMesh> > outMesh2DSplit;
  // OK all is ready to insert in ret2 mesh
  if(!untouchedCells->empty())
    {// the most easy part, cells in mesh2D not impacted at all
      outMesh2DSplit.push_back(static_cast<MEDCouplingUMesh *>(mesh2D->buildPartOfMySelf(untouchedCells->begin(),untouchedCells->end())));
      outMesh2DSplit.back()->setCoords(ret1->getCoords());
      ret2->pushBackValsSilent(untouchedCells->begin(),untouchedCells->end());
    }
  if((DataArrayInt *)out0s)
    {// here dealing with cells in out0s but not in cellsToBeModified
      MCAuto<DataArrayInt> fewModifiedCells(out0s->buildSubstraction(cellsToBeModified));
      const int *rdptr(dd3->begin()),*rdiptr(dd4->begin()),*dptr(dd1->begin()),*diptr(dd2->begin());
      for(const int *it=fewModifiedCells->begin();it!=fewModifiedCells->end();it++)
        {
          outMesh2DSplit.push_back(BuildRefined2DCell(ret1->getCoords(),mesh2D,*it,dptr+diptr[*it],dptr+diptr[*it+1],intersectEdge1));
          ret1->setCoords(outMesh2DSplit.back()->getCoords());
        }
      int offset(ret2->getNumberOfTuples());
      ret2->pushBackValsSilent(fewModifiedCells->begin(),fewModifiedCells->end());
      MCAuto<DataArrayInt> partOfRet3(DataArrayInt::New()); partOfRet3->alloc(2*idsInRet1Colinear->getNumberOfTuples(),1);
      partOfRet3->fillWithValue(std::numeric_limits<int>::max()); partOfRet3->rearrange(2);
      int kk(0),*ret3ptr(partOfRet3->getPointer());
      for(const int *it=idsInDescMesh2DForIdsInRetColinear->begin();it!=idsInDescMesh2DForIdsInRetColinear->end();it++,kk++)
        {
          int faceId(std::abs(*it)-1);
          for(const int *it2=rdptr+rdiptr[faceId];it2!=rdptr+rdiptr[faceId+1];it2++)
            {
              int tmp(fewModifiedCells->findIdFirstEqual(*it2));
              if(tmp!=-1)
                {
                  if(std::find(dptr+diptr[*it2],dptr+diptr[*it2+1],-(*it))!=dptr+diptr[*it2+1])
                    ret3ptr[2*kk]=tmp+offset;
                  if(std::find(dptr+diptr[*it2],dptr+diptr[*it2+1],(*it))!=dptr+diptr[*it2+1])
                    ret3ptr[2*kk+1]=tmp+offset;
                }
              else
                {//the current edge is shared by a 2D cell that will be split just after
                  if(std::find(dptr+diptr[*it2],dptr+diptr[*it2+1],-(*it))!=dptr+diptr[*it2+1])
                    ret3ptr[2*kk]=-(*it2+1);
                  if(std::find(dptr+diptr[*it2],dptr+diptr[*it2+1],(*it))!=dptr+diptr[*it2+1])
                    ret3ptr[2*kk+1]=-(*it2+1);
                }
            }
        }
      m1Desc->setCoords(ret1->getCoords());
      ret1NonCol->setCoords(ret1->getCoords());
      ret3->setPartOfValues3(partOfRet3,idsInRet1Colinear->begin(),idsInRet1Colinear->end(),0,2,1,true);
      if(!outMesh2DSplit.empty())
        {
          DataArrayDouble *da(outMesh2DSplit.back()->getCoords());
          for(std::vector< MCAuto<MEDCouplingUMesh> >::iterator itt=outMesh2DSplit.begin();itt!=outMesh2DSplit.end();itt++)
            (*itt)->setCoords(da);
        }
    }
  cellsToBeModified=cellsToBeModified->buildUniqueNotSorted();
  for(const int *it=cellsToBeModified->begin();it!=cellsToBeModified->end();it++)
    {
      MCAuto<DataArrayInt> idsNonColPerCell(elts->findIdsEqual(*it));
      idsNonColPerCell->transformWithIndArr(eltsIndex3->begin(),eltsIndex3->end());
      MCAuto<DataArrayInt> idsNonColPerCell2(idsInRet1NotColinear->selectByTupleId(idsNonColPerCell->begin(),idsNonColPerCell->end()));
      MCAuto<MEDCouplingUMesh> partOfMesh1CuttingCur2DCell(static_cast<MEDCouplingUMesh *>(ret1NonCol->buildPartOfMySelf(idsNonColPerCell->begin(),idsNonColPerCell->end())));
      MCAuto<DataArrayInt> partOfRet3;
      MCAuto<MEDCouplingUMesh> splitOfOneCell(BuildMesh2DCutFrom(eps,*it,m1Desc,partOfMesh1CuttingCur2DCell,dd1->begin()+dd2->getIJ(*it,0),dd1->begin()+dd2->getIJ((*it)+1,0),intersectEdge1,ret2->getNumberOfTuples(),partOfRet3));
      ret3->setPartOfValues3(partOfRet3,idsNonColPerCell2->begin(),idsNonColPerCell2->end(),0,2,1,true);
      outMesh2DSplit.push_back(splitOfOneCell);
      for(std::size_t i=0;i<splitOfOneCell->getNumberOfCells();i++)
        ret2->pushBackSilent(*it);
    }
  //
  std::size_t nbOfMeshes(outMesh2DSplit.size());
  std::vector<const MEDCouplingUMesh *> tmp(nbOfMeshes);
  for(std::size_t i=0;i<nbOfMeshes;i++)
    tmp[i]=outMesh2DSplit[i];
  //
  ret1->getCoords()->setInfoOnComponents(compNames);
  MCAuto<MEDCouplingUMesh> ret2D(MEDCouplingUMesh::MergeUMeshesOnSameCoords(tmp));
  // To finish - filter ret3 - std::numeric_limits<int>::max() -> -1 - negate values must be resolved.
  ret3->rearrange(1);
  MCAuto<DataArrayInt> edgesToDealWith(ret3->findIdsStrictlyNegative());
  for(const int *it=edgesToDealWith->begin();it!=edgesToDealWith->end();it++)
    {
      int old2DCellId(-ret3->getIJ(*it,0)-1);
      MCAuto<DataArrayInt> candidates(ret2->findIdsEqual(old2DCellId));
      ret3->setIJ(*it,0,FindRightCandidateAmong(ret2D,candidates->begin(),candidates->end(),ret1,*it%2==0?-((*it)/2+1):(*it)/2+1,eps));// div by 2 because 2 components natively in ret3
    }
  ret3->changeValue(std::numeric_limits<int>::max(),-1);
  ret3->rearrange(2);
  //
  splitMesh1D=ret1.retn();
  splitMesh2D=ret2D.retn();
  cellIdInMesh2D=ret2.retn();
  cellIdInMesh1D=ret3.retn();
}

/*!
 * \b WARNING this method is \b potentially \b non \b const (if returned array is empty).
 * \b WARNING this method lead to have a non geometric type sorted mesh (for MED file users) !
 * This method performs a conformization of \b this. So if a edge in \a this can be split into entire edges in \a this this method
 * will suppress such edges to use sub edges in \a this. So this method does not add nodes in \a this if merged edges are both linear (INTERP_KERNEL::NORM_SEG2).
 * In the other cases new nodes can be created. If any are created, they will be appended at the end of the coordinates object before the invokation of this method.
 *
 * Whatever the returned value, this method does not alter the order of cells in \a this neither the orientation of cells.
 * The modified cells, if any, are systematically declared as NORM_POLYGON or NORM_QPOLYG depending on the initial quadraticness of geometric type.
 *
 * This method expects that \b this has a meshDim equal 2 and spaceDim equal to 2 too.
 * This method expects that all nodes in \a this are not closer than \a eps.
 * If it is not the case you can invoke MEDCouplingUMesh::mergeNodes before calling this method.
 *
 * \param [in] eps the relative error to detect merged edges.
 * \return DataArrayInt  * - The list of cellIds in \a this that have been subdivided. If empty, nothing changed in \a this (as if it were a const method). The array is a newly allocated array
 *                           that the user is expected to deal with.
 *
 * \throw If \a this is not coherent.
 * \throw If \a this has not spaceDim equal to 2.
 * \throw If \a this has not meshDim equal to 2.
 * \sa MEDCouplingUMesh::mergeNodes, MEDCouplingUMesh::split2DCells
 */
DataArrayInt *MEDCouplingUMesh::conformize2D(double eps)
{
  static const int SPACEDIM=2;
  checkConsistencyLight();
  if(getSpaceDimension()!=2 || getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::conformize2D : This method only works for meshes with spaceDim=2 and meshDim=2 !");
  MCAuto<DataArrayInt> desc1(DataArrayInt::New()),descIndx1(DataArrayInt::New()),revDesc1(DataArrayInt::New()),revDescIndx1(DataArrayInt::New());
  MCAuto<MEDCouplingUMesh> mDesc(buildDescendingConnectivity(desc1,descIndx1,revDesc1,revDescIndx1));
  const int *c(mDesc->getNodalConnectivity()->begin()),*ci(mDesc->getNodalConnectivityIndex()->begin()),*rd(revDesc1->begin()),*rdi(revDescIndx1->begin());
  MCAuto<DataArrayDouble> bboxArr(mDesc->getBoundingBoxForBBTree(eps));
  const double *bbox(bboxArr->begin()),*coords(getCoords()->begin());
  int nCell(getNumberOfCells()),nDescCell(mDesc->getNumberOfCells());
  std::vector< std::vector<int> > intersectEdge(nDescCell),overlapEdge(nDescCell);
  std::vector<double> addCoo;
  BBTree<SPACEDIM,int> myTree(bbox,0,0,nDescCell,-eps);
  INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision arcPrec(eps);
  for(int i=0;i<nDescCell;i++)
    {
      std::vector<int> candidates;
      myTree.getIntersectingElems(bbox+i*2*SPACEDIM,candidates);
      for(std::vector<int>::const_iterator it=candidates.begin();it!=candidates.end();it++)
        if(*it>i)
          {
            std::map<MCAuto<INTERP_KERNEL::Node>,int> m;
            INTERP_KERNEL::Edge *e1(MEDCouplingUMeshBuildQPFromEdge2((INTERP_KERNEL::NormalizedCellType)c[ci[i]],c+ci[i]+1,coords,m)),
                *e2(MEDCouplingUMeshBuildQPFromEdge2((INTERP_KERNEL::NormalizedCellType)c[ci[*it]],c+ci[*it]+1,coords,m));
            INTERP_KERNEL::MergePoints merge;
            INTERP_KERNEL::QuadraticPolygon c1,c2;
            e1->intersectWith(e2,merge,c1,c2);
            e1->decrRef(); e2->decrRef();
            if(IKGeo2DInternalMapper(c1,m,c[ci[i]+1],c[ci[i]+2],intersectEdge[i]))
              overlapEdge[i].push_back(*it);
            if(IKGeo2DInternalMapper(c2,m,c[ci[*it]+1],c[ci[*it]+2],intersectEdge[*it]))
              overlapEdge[*it].push_back(i);
          }
    }
  // splitting done. sort intersect point in intersectEdge.
  std::vector< std::vector<int> > middle(nDescCell);
  int nbOf2DCellsToBeSplit(0);
  bool middleNeedsToBeUsed(false);
  std::vector<bool> cells2DToTreat(nDescCell,false);
  for(int i=0;i<nDescCell;i++)
    {
      std::vector<int>& isect(intersectEdge[i]);
      int sz((int)isect.size());
      if(sz>1)
        {
          std::map<MCAuto<INTERP_KERNEL::Node>,int> m;
          INTERP_KERNEL::Edge *e(MEDCouplingUMeshBuildQPFromEdge2((INTERP_KERNEL::NormalizedCellType)c[ci[i]],c+ci[i]+1,coords,m));
          e->sortSubNodesAbs(coords,isect);
          e->decrRef();
        }
      if(sz!=0)
        {
          int idx0(rdi[i]),idx1(rdi[i+1]);
          if(idx1-idx0!=1)
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::conformize2D : internal error #0 !");
          if(!cells2DToTreat[rd[idx0]])
            {
              cells2DToTreat[rd[idx0]]=true;
              nbOf2DCellsToBeSplit++;
            }
          // try to reuse at most eventual 'middle' of SEG3
          std::vector<int>& mid(middle[i]);
          mid.resize(sz+1,-1);
          if((INTERP_KERNEL::NormalizedCellType)c[ci[i]]==INTERP_KERNEL::NORM_SEG3)
            {
              middleNeedsToBeUsed=true;
              const std::vector<int>& candidates(overlapEdge[i]);
              std::vector<int> trueCandidates;
              for(std::vector<int>::const_iterator itc=candidates.begin();itc!=candidates.end();itc++)
                if((INTERP_KERNEL::NormalizedCellType)c[ci[*itc]]==INTERP_KERNEL::NORM_SEG3)
                  trueCandidates.push_back(*itc);
              int stNode(c[ci[i]+1]),endNode(isect[0]);
              for(int j=0;j<sz+1;j++)
                {
                  for(std::vector<int>::const_iterator itc=trueCandidates.begin();itc!=trueCandidates.end();itc++)
                    {
                      int tmpSt(c[ci[*itc]+1]),tmpEnd(c[ci[*itc]+2]);
                      if((tmpSt==stNode && tmpEnd==endNode) || (tmpSt==endNode && tmpEnd==stNode))
                        { mid[j]=*itc; break; }
                    }
                  stNode=endNode;
                  endNode=j<sz-1?isect[j+1]:c[ci[i]+2];
                }
            }
        }
    }
  MCAuto<DataArrayInt> ret(DataArrayInt::New()),notRet(DataArrayInt::New()); ret->alloc(nbOf2DCellsToBeSplit,1);
  if(nbOf2DCellsToBeSplit==0)
    return ret.retn();
  //
  int *retPtr(ret->getPointer());
  for(int i=0;i<nCell;i++)
    if(cells2DToTreat[i])
      *retPtr++=i;
  //
  MCAuto<DataArrayInt> mSafe,nSafe,oSafe,pSafe,qSafe,rSafe;
  DataArrayInt *m(0),*n(0),*o(0),*p(0),*q(0),*r(0);
  MEDCouplingUMesh::ExtractFromIndexedArrays(ret->begin(),ret->end(),desc1,descIndx1,m,n); mSafe=m; nSafe=n;
  DataArrayInt::PutIntoToSkylineFrmt(intersectEdge,o,p); oSafe=o; pSafe=p;
  if(middleNeedsToBeUsed)
    { DataArrayInt::PutIntoToSkylineFrmt(middle,q,r); qSafe=q; rSafe=r; }
  MCAuto<MEDCouplingUMesh> modif(static_cast<MEDCouplingUMesh *>(buildPartOfMySelf(ret->begin(),ret->end(),true)));
  int nbOfNodesCreated(modif->split2DCells(mSafe,nSafe,oSafe,pSafe,qSafe,rSafe));
  setCoords(modif->getCoords());//if nbOfNodesCreated==0 modif and this have the same coordinates pointer so this line has no effect. But for quadratic cases this line is important.
  setPartOfMySelf(ret->begin(),ret->end(),*modif);
  {
    bool areNodesMerged; int newNbOfNodes;
    if(nbOfNodesCreated!=0)
      MCAuto<DataArrayInt> tmp(mergeNodes(eps,areNodesMerged,newNbOfNodes));
  }
  return ret.retn();
}

/*!
 * This non const method works on 2D mesh. This method scans every cell in \a this and look if each edge constituting this cell is not mergeable with neighbors edges of that cell.
 * If yes, the cell is "repaired" to minimize at most its number of edges. So this method do not change the overall shape of cells in \a this (with eps precision).
 * This method do not take care of shared edges between cells, so this method can lead to a non conform mesh (\a this). If a conform mesh is required you're expected
 * to invoke MEDCouplingUMesh::mergeNodes and MEDCouplingUMesh::conformize2D right after this call.
 * This method works on any 2D geometric types of cell (even static one). If a cell is touched its type becomes dynamic automaticaly. For 2D "repaired" quadratic cells
 * new nodes for center of merged edges is are systematically created and appended at the end of the previously existing nodes.
 *
 * If the returned array is empty it means that nothing has changed in \a this (as if it were a const method). If the array is not empty the connectivity of \a this is modified
 * using new instance, idem for coordinates.
 *
 * If \a this is constituted by only linear 2D cells, this method is close to the computation of the convex hull of each cells in \a this.
 *
 * \return DataArrayInt  * - The list of cellIds in \a this that have at least one edge colinearized.
 *
 * \throw If \a this is not coherent.
 * \throw If \a this has not spaceDim equal to 2.
 * \throw If \a this has not meshDim equal to 2.
 *
 * \sa MEDCouplingUMesh::conformize2D, MEDCouplingUMesh::mergeNodes, MEDCouplingUMesh::convexEnvelop2D.
 */
DataArrayInt *MEDCouplingUMesh::colinearize2D(double eps)
{
  MCAuto<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
  checkConsistencyLight();
  if(getSpaceDimension()!=2 || getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::colinearize2D : This method only works for meshes with spaceDim=2 and meshDim=2 !");
  INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision arcPrec(eps);
  int nbOfCells(getNumberOfCells()),nbOfNodes(getNumberOfNodes());
  const int *cptr(_nodal_connec->begin()),*ciptr(_nodal_connec_index->begin());
  MCAuto<DataArrayInt> newc(DataArrayInt::New()),newci(DataArrayInt::New()); newci->alloc(nbOfCells+1,1); newc->alloc(0,1); newci->setIJ(0,0,0);
  MCAuto<DataArrayDouble> appendedCoords(DataArrayDouble::New()); appendedCoords->alloc(0,1);//1 not 2 it is not a bug.
  const double *coords(_coords->begin());
  int *newciptr(newci->getPointer());
  for(int i=0;i<nbOfCells;i++,newciptr++,ciptr++)
    {
      if(Colinearize2DCell(coords,cptr+ciptr[0],cptr+ciptr[1],nbOfNodes,newc,appendedCoords))
        ret->pushBackSilent(i);
      newciptr[1]=newc->getNumberOfTuples();
    }
  //
  if(ret->empty())
    return ret.retn();
  if(!appendedCoords->empty())
    {
      appendedCoords->rearrange(2);
      MCAuto<DataArrayDouble> newCoords(DataArrayDouble::Aggregate(getCoords(),appendedCoords));//treat info on components
      //non const part
      setCoords(newCoords);
    }
  //non const part
  setConnectivity(newc,newci,true);
  return ret.retn();
}

///@cond INTERNAL
/**
 * c, cI describe a wire mesh in 3D space, in local numbering
 * startNode, endNode in global numbering
 *\return true if the segment is indeed split
 */
bool MEDCouplingUMesh::OrderPointsAlongLine(const double * coo, int startNode, int endNode,
                                            const int * c, const int * cI, const int *idsBg, const int *endBg,
                                            std::vector<int> & pointIds, std::vector<int> & hitSegs)
{
  using namespace std;

  const int SPACEDIM=3;
  typedef pair<double, int> PairDI;
  set< PairDI > x;
  for (const int * it = idsBg; it != endBg; ++it)
    {
      assert(c[cI[*it]] == INTERP_KERNEL::NORM_SEG2);
      int start = c[cI[*it]+1], end = c[cI[*it]+2];
      x.insert(make_pair(coo[start*SPACEDIM], start));  // take only X coords
      x.insert(make_pair(coo[end*SPACEDIM], end));
    }

  vector<PairDI> xx(x.begin(), x.end());
  sort(xx.begin(),xx.end());
  pointIds.reserve(xx.size());

  // Keep what is inside [startNode, endNode]:
  int go = 0;
  for (vector<PairDI>::const_iterator it=xx.begin(); it != xx.end(); ++it)
    {
      const int idx = (*it).second;
      if (!go)
        {
          if (idx == startNode)   go = 1;
          if (idx == endNode)     go = 2;
          if (go)                 pointIds.push_back(idx);
          continue;
        }
      pointIds.push_back(idx);
      if (idx == endNode || idx == startNode)
        break;
    }

//  vector<int> pointIds2(pointIds.size()+2);
//  copy(pointIds.begin(), pointIds.end(), pointIds2.data()+1);
//  pointIds2[0] = startNode;
//  pointIds2[pointIds2.size()-1] = endNode;

  if (go == 2)
    reverse(pointIds.begin(), pointIds.end());

  // Now identify smaller segments that are not sub-divided - those won't need any further treatment:
  for (const int * it = idsBg; it != endBg; ++it)
    {
      int start = c[cI[*it]+1], end = c[cI[*it]+2];
      vector<int>::const_iterator itStart = find(pointIds.begin(), pointIds.end(), start);
      if (itStart == pointIds.end()) continue;
      vector<int>::const_iterator itEnd = find(pointIds.begin(), pointIds.end(), end);
      if (itEnd == pointIds.end())               continue;
      if (abs(distance(itEnd, itStart)) != 1)    continue;
      hitSegs.push_back(*it);   // segment is undivided.
    }

  return (pointIds.size() > 2); // something else apart start and end node
}

void MEDCouplingUMesh::ReplaceEdgeInFace(const int * sIdxConn, const int * sIdxConnE, int startNode, int endNode,
                                          const std::vector<int>& insidePoints, std::vector<int>& modifiedFace)
{
  using namespace std;
  int dst = distance(sIdxConn, sIdxConnE);
  modifiedFace.reserve(dst + insidePoints.size()-2);
  modifiedFace.resize(dst);
  copy(sIdxConn, sIdxConnE, modifiedFace.data());

  vector<int>::iterator shortEnd = modifiedFace.begin()+dst;
  vector<int>::iterator startPos = find(modifiedFace.begin(), shortEnd , startNode);
  if (startPos == shortEnd)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ReplaceEdgeInFace: internal error, should never happen!");
  vector<int>::iterator endPos = find(modifiedFace.begin(),shortEnd, endNode);
  if (endPos == shortEnd)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ReplaceEdgeInFace: internal error, should never happen!");
  int d = distance(startPos, endPos);
  if (d == 1 || d == (1-dst)) // don't use modulo, for neg numbers, result is implementation defined ...
    modifiedFace.insert(++startPos, ++insidePoints.begin(), --insidePoints.end());  // insidePoints also contains start and end node. Those dont need to be inserted.
  else
    modifiedFace.insert(++endPos, ++insidePoints.rbegin(), --insidePoints.rend());
}

///@endcond


/*!
 * \b WARNING this method is \b potentially \b non \b const (if returned array is not empty).
 * \b WARNING this method lead to have a non geometric type sorted mesh (for MED file users) !
 * This method performs a conformization of \b this.
 *
 * Only polyhedron cells are supported. You can call convertAllToPoly()
 *
 * This method expects that \b this has a meshDim equal 3 and spaceDim equal to 3 too.
 * This method expects that all nodes in \a this are not closer than \a eps.
 * If it is not the case you can invoke MEDCouplingUMesh::mergeNodes before calling this method.
 *
 * \param [in] eps the relative error to detect merged edges.
 * \return DataArrayInt  * - The list of cellIds in \a this that have been subdivided. If empty, nothing changed in \a this (as if it were a const method). The array is a newly allocated array
 *                           that the user is expected to deal with.
 *
 * \throw If \a this is not coherent.
 * \throw If \a this has not spaceDim equal to 3.
 * \throw If \a this has not meshDim equal to 3.
 * \sa MEDCouplingUMesh::mergeNodes, MEDCouplingUMesh::conformize2D, MEDCouplingUMesh::convertAllToPoly()
 */
DataArrayInt *MEDCouplingUMesh::conformize3D(double eps)
{
  using namespace std;

  static const int SPACEDIM=3;
  checkConsistencyLight();
  if(getSpaceDimension()!=3 || getMeshDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::conformize3D : This method only works for meshes with spaceDim=3 and meshDim=3!");
  if(_types.size() != 1 || *(_types.begin()) != INTERP_KERNEL::NORM_POLYHED)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::conformize3D : This method only works for polyhedrons! Call convertAllToPoly first.");

  MCAuto<MEDCouplingSkyLineArray> connSla(MEDCouplingSkyLineArray::BuildFromPolyhedronConn(getNodalConnectivity(), getNodalConnectivityIndex()));
  const double * coo(_coords->begin());
  MCAuto<DataArrayInt> ret(DataArrayInt::New());

  {
    /*************************
     *  STEP 1  -- faces
     *************************/
    MCAuto<DataArrayInt> descDNU(DataArrayInt::New()),descIDNU(DataArrayInt::New()),revDesc(DataArrayInt::New()),revDescI(DataArrayInt::New());
    MCAuto<MEDCouplingUMesh> mDesc(buildDescendingConnectivity(descDNU,descIDNU,revDesc,revDescI));
    const int *revDescIP(revDescI->getConstPointer()), *revDescP(revDesc->getConstPointer());
    const int *cDesc(mDesc->getNodalConnectivity()->begin()),*cIDesc(mDesc->getNodalConnectivityIndex()->begin());
    MCAuto<MEDCouplingSkyLineArray> connSlaDesc(MEDCouplingSkyLineArray::New(mDesc->getNodalConnectivityIndex(), mDesc->getNodalConnectivity()));

    // Build BBTree
    MCAuto<DataArrayDouble> bboxArr(mDesc->getBoundingBoxForBBTree(eps));
    const double *bbox(bboxArr->begin()); getCoords()->begin();
    int nDescCell(mDesc->getNumberOfCells());
    BBTree<SPACEDIM,int> myTree(bbox,0,0,nDescCell,-eps);
    // Surfaces - handle biggest first
    MCAuto<MEDCouplingFieldDouble> surfF = mDesc->getMeasureField(true);
    DataArrayDouble * surfs = surfF->getArray();
    // Normal field
    MCAuto<MEDCouplingFieldDouble> normalsF = mDesc->buildOrthogonalField();
    DataArrayDouble * normals = normalsF->getArray();
    const double * normalsP = normals->getConstPointer();

    // Sort faces by decreasing surface:
    vector< pair<double,int> > S;
    for(std::size_t i=0;i < surfs->getNumberOfTuples();i++)
      {
        pair<double,int> p = make_pair(surfs->begin()[i], i);
        S.push_back(p);
      }
    sort(S.rbegin(),S.rend()); // reverse sort
    vector<bool> hit(nDescCell);
    fill(hit.begin(), hit.end(), false);
    vector<int> hitPoly; // the final result: which 3D cells have been modified.

    for( vector<pair<double,int> >::const_iterator it = S.begin(); it != S.end(); it++)
      {
        int faceIdx = (*it).second;
        if (hit[faceIdx]) continue;

        vector<int> candidates, cands2;
        myTree.getIntersectingElems(bbox+faceIdx*2*SPACEDIM,candidates);
        // Keep only candidates whose normal matches the normal of current face
        for(vector<int>::const_iterator it2=candidates.begin();it2!=candidates.end();it2++)
          {
            bool col = INTERP_KERNEL::isColinear3D(normalsP + faceIdx*SPACEDIM, normalsP + *(it2)*SPACEDIM, eps);
            if (*it2 != faceIdx && col)
              cands2.push_back(*it2);
          }
        if (!cands2.size())
          continue;

        // Now rotate, and match barycenters -- this is where we will bring Intersect2DMeshes later
        INTERP_KERNEL::TranslationRotationMatrix rotation;
        INTERP_KERNEL::TranslationRotationMatrix::Rotate3DTriangle(coo+SPACEDIM*(cDesc[cIDesc[faceIdx]+1]),
                                                                   coo+SPACEDIM*(cDesc[cIDesc[faceIdx]+2]),
                                                                   coo+SPACEDIM*(cDesc[cIDesc[faceIdx]+3]), rotation);

        MCAuto<MEDCouplingUMesh> mPartRef(mDesc->buildPartOfMySelfSlice(faceIdx, faceIdx+1,1,false));  // false=zipCoords is called
        MCAuto<MEDCouplingUMesh> mPartCand(mDesc->buildPartOfMySelf(&cands2[0], &cands2[0]+cands2.size(), false)); // false=zipCoords is called
        double * cooPartRef(mPartRef->_coords->getPointer());
        double * cooPartCand(mPartCand->_coords->getPointer());
        for (std::size_t ii = 0; ii < mPartRef->_coords->getNumberOfTuples(); ii++)
          rotation.transform_vector(cooPartRef+SPACEDIM*ii);
        for (std::size_t ii = 0; ii < mPartCand->_coords->getNumberOfTuples(); ii++)
          rotation.transform_vector(cooPartCand+SPACEDIM*ii);

        // Localize faces in 2D thanks to barycenters
        MCAuto<DataArrayDouble> baryPart = mPartCand->computeCellCenterOfMass();
        vector<int> compo; compo.push_back(2);
        MCAuto<DataArrayDouble> baryPartZ = baryPart->keepSelectedComponents(compo);
        MCAuto<DataArrayInt> idsGoodPlane = baryPartZ->findIdsInRange(-eps, +eps);
        if (!idsGoodPlane->getNumberOfTuples())
          continue;

        baryPart = baryPart->selectByTupleId(*idsGoodPlane);

        compo[0] = 0; compo.push_back(1);
        MCAuto<DataArrayDouble> baryPartXY = baryPart->keepSelectedComponents(compo);
        mPartRef->changeSpaceDimension(2,0.0);
        MCAuto<DataArrayInt> cc(DataArrayInt::New()), ccI(DataArrayInt::New());
        mPartRef->getCellsContainingPoints(baryPartXY->begin(), baryPartXY->getNumberOfTuples(), eps, cc, ccI);

        if (!cc->getNumberOfTuples())
          continue;
        MCAuto<DataArrayInt> dsi = ccI->deltaShiftIndex();

        {
          MCAuto<DataArrayInt> tmp = dsi->findIdsInRange(0, 2);
          if (tmp->getNumberOfTuples() != dsi->getNumberOfTuples())
            {
              ostringstream oss;
              oss << "MEDCouplingUMesh::conformize3D: Non expected non-conformity. Only simple (=partition-like) non-conformities are handled. Face #" << faceIdx << " violates this condition!";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }

        MCAuto<DataArrayInt> ids = dsi->findIdsEqual(1);
        // Boundary face:
        if (!ids->getNumberOfTuples())
          continue;

        double checkSurf=0.0;
        const int * idsGoodPlaneP(idsGoodPlane->begin());
        for (const int * ii = ids->begin(); ii != ids->end(); ii++)
          {
            int faceIdx2 = cands2[idsGoodPlaneP[*ii]];
            hit[faceIdx2] = true;
            checkSurf += surfs->begin()[faceIdx2];
          }
        if (fabs(checkSurf - surfs->begin()[faceIdx]) > eps)
          {
            ostringstream oss;
            oss << "MEDCouplingUMesh::conformize3D: Non expected non-conformity. Only simple (=partition-like) non-conformities are handled. Face #" << faceIdx << " violates this condition!";
            throw INTERP_KERNEL::Exception(oss.str());
          }

        // For all polyhedrons using this face, replace connectivity:
        vector<int> polyIndices, packsIds, facePack;
        for (int ii=revDescIP[faceIdx]; ii < revDescIP[faceIdx+1]; ii++)
          polyIndices.push_back(revDescP[ii]);
        ret->pushBackValsSilent(polyIndices.data(),polyIndices.data()+polyIndices.size());

        // Current face connectivity
        const int * sIdxConn = cDesc + cIDesc[faceIdx] + 1;
        const int * sIdxConnE = cDesc + cIDesc[faceIdx+1];
        connSla->findPackIds(polyIndices, sIdxConn, sIdxConnE, packsIds);
        // Deletion of old faces
        int jj=0;
        for (vector<int>::const_iterator it2=polyIndices.begin(); it2!=polyIndices.end(); ++it2, ++jj)
          {
            if (packsIds[jj] == -1)
              // The below should never happen - if a face is used several times, with a different layout of the nodes
              // it means that is is already conform, so it is *not* hit by the algorithm. The algorithm only hits
              // faces which are actually used only once, by a single cell. This is different for edges below.
              throw INTERP_KERNEL::Exception("MEDCouplingUMesh::conformize3D: Could not find face in connectivity! Internal error.");
            else
              connSla->deletePack(*it2, packsIds[jj]);
          }
        // Insertion of new faces:
        for (const int * ii = ids->begin(); ii != ids->end(); ii++)
          {
            // Build pack from the face to insert:
            int faceIdx2 = cands2[idsGoodPlane->getIJ(*ii,0)];
            int facePack2Sz;
            const int * facePack2 = connSlaDesc->getSimplePackSafePtr(faceIdx2, facePack2Sz); // contains the type!
            // Insert it in all hit polyhedrons:
            for (vector<int>::const_iterator it2=polyIndices.begin(); it2!=polyIndices.end(); ++it2)
              connSla->pushBackPack(*it2, facePack2+1, facePack2+facePack2Sz);  // without the type
          }
      }
  }  // end step1

  // Set back modified connectivity
  MCAuto<DataArrayInt> cAuto; cAuto.takeRef(_nodal_connec);
  MCAuto<DataArrayInt> cIAuto; cIAuto.takeRef(_nodal_connec_index);
  connSla->convertToPolyhedronConn(cAuto, cIAuto);

  {
    /************************
     *  STEP 2 -- edges
     ************************/
    // Now we have a face-conform mesh.

    // Recompute descending
    MCAuto<DataArrayInt> desc(DataArrayInt::New()),descI(DataArrayInt::New()),revDesc(DataArrayInt::New()),revDescI(DataArrayInt::New());
    // Rebuild desc connectivity with orientation this time!!
    MCAuto<MEDCouplingUMesh> mDesc(buildDescendingConnectivity2(desc,descI,revDesc,revDescI));
    const int *revDescIP(revDescI->getConstPointer()), *revDescP(revDesc->getConstPointer());
    const int *descIP(descI->getConstPointer()), *descP(desc->getConstPointer());
    const int *cDesc(mDesc->getNodalConnectivity()->begin()),*cIDesc(mDesc->getNodalConnectivityIndex()->begin());
    MCAuto<DataArrayInt> ciDeepC(mDesc->getNodalConnectivityIndex()->deepCopy());
    MCAuto<DataArrayInt> cDeepC(mDesc->getNodalConnectivity()->deepCopy());
    MCAuto<MEDCouplingSkyLineArray> connSlaDesc(MEDCouplingSkyLineArray::New(ciDeepC, cDeepC));
    MCAuto<DataArrayInt> desc2(DataArrayInt::New()),descI2(DataArrayInt::New()),revDesc2(DataArrayInt::New()),revDescI2(DataArrayInt::New());
    MCAuto<MEDCouplingUMesh> mDesc2 = mDesc->buildDescendingConnectivity(desc2,descI2,revDesc2,revDescI2);
//    std::cout << "writing!\n";
//    mDesc->writeVTK("/tmp/toto_desc_confInter.vtu");
//    mDesc2->writeVTK("/tmp/toto_desc2_confInter.vtu");
    const int *revDescIP2(revDescI2->getConstPointer()), *revDescP2(revDesc2->getConstPointer());
    const int *cDesc2(mDesc2->getNodalConnectivity()->begin()),*cIDesc2(mDesc2->getNodalConnectivityIndex()->begin());
    MCAuto<DataArrayDouble> bboxArr(mDesc2->getBoundingBoxForBBTree(eps));
    const double *bbox2(bboxArr->begin());
    int nDesc2Cell=mDesc2->getNumberOfCells();
    BBTree<SPACEDIM,int> myTree2(bbox2,0,0,nDesc2Cell,-eps);

    // Edges - handle longest first
    MCAuto<MEDCouplingFieldDouble> lenF = mDesc2->getMeasureField(true);
    DataArrayDouble * lens = lenF->getArray();

    // Sort edges by decreasing length:
    vector<pair<double,int> > S;
    for(std::size_t i=0;i < lens->getNumberOfTuples();i++)
      {
        pair<double,int> p = make_pair(lens->getIJ(i, 0), i);
        S.push_back(p);
      }
    sort(S.rbegin(),S.rend()); // reverse sort

    vector<bool> hit(nDesc2Cell);
    fill(hit.begin(), hit.end(), false);

    for( vector<pair<double,int> >::const_iterator it = S.begin(); it != S.end(); it++)
      {
        int eIdx = (*it).second;
        if (hit[eIdx])
          continue;

        vector<int> candidates, cands2;
        myTree2.getIntersectingElems(bbox2+eIdx*2*SPACEDIM,candidates);
        // Keep only candidates colinear with current edge
        double vCurr[3];
        unsigned start = cDesc2[cIDesc2[eIdx]+1], end = cDesc2[cIDesc2[eIdx]+2];
        for (int i3=0; i3 < 3; i3++)  // TODO: use fillSonCellNodalConnectivity2 or similar?
          vCurr[i3] = coo[start*SPACEDIM+i3] - coo[end*SPACEDIM+i3];
        for(vector<int>::const_iterator it2=candidates.begin();it2!=candidates.end();it2++)
          {
            double vOther[3];
            unsigned start2 = cDesc2[cIDesc2[*it2]+1], end2 = cDesc2[cIDesc2[*it2]+2];
            for (int i3=0; i3 < 3; i3++)
              vOther[i3] = coo[start2*SPACEDIM+i3] - coo[end2*SPACEDIM+i3];
            bool col = INTERP_KERNEL::isColinear3D(vCurr, vOther, eps);
            // Warning: different from faces: we need to keep eIdx in the final list of candidates because we need
            // to have its nodes inside the sub mesh mPartCand below (needed in OrderPointsAlongLine())
            if (col)
              cands2.push_back(*it2);
          }
        if (cands2.size() == 1 && cands2[0] == eIdx)  // see warning above
          continue;

        // Now rotate edges to bring them on Ox
        int startNode = cDesc2[cIDesc2[eIdx]+1];
        int endNode = cDesc2[cIDesc2[eIdx]+2];
        INTERP_KERNEL::TranslationRotationMatrix rotation;
        INTERP_KERNEL::TranslationRotationMatrix::Rotate3DBipoint(coo+SPACEDIM*startNode, coo+SPACEDIM*endNode, rotation);
        MCAuto<MEDCouplingUMesh> mPartRef(mDesc2->buildPartOfMySelfSlice(eIdx, eIdx+1,1,false));  // false=zipCoords is called
        MCAuto<MEDCouplingUMesh> mPartCand(mDesc2->buildPartOfMySelf(&cands2[0], &cands2[0]+cands2.size(), true)); // true=zipCoords is called
        MCAuto<DataArrayInt> nodeMap(mPartCand->zipCoordsTraducer());
        int nbElemsNotM1;
        {
          MCAuto<DataArrayInt> tmp(nodeMap->findIdsNotEqual(-1));
          nbElemsNotM1 = tmp->getNbOfElems();
        }
        MCAuto<DataArrayInt>  nodeMapInv = nodeMap->invertArrayO2N2N2O(nbElemsNotM1);
        double * cooPartRef(mPartRef->_coords->getPointer());
        double * cooPartCand(mPartCand->_coords->getPointer());
        for (std::size_t ii = 0; ii < mPartRef->_coords->getNumberOfTuples(); ii++)
          rotation.transform_vector(cooPartRef+SPACEDIM*ii);
        for (std::size_t ii = 0; ii < mPartCand->_coords->getNumberOfTuples(); ii++)
          rotation.transform_vector(cooPartCand+SPACEDIM*ii);


        // Eliminate all edges for which y or z is not null
        MCAuto<DataArrayDouble> baryPart = mPartCand->computeCellCenterOfMass();
        vector<int> compo; compo.push_back(1);
        MCAuto<DataArrayDouble> baryPartY = baryPart->keepSelectedComponents(compo);
        compo[0] = 2;
        MCAuto<DataArrayDouble> baryPartZ = baryPart->keepSelectedComponents(compo);
        MCAuto<DataArrayInt> idsGoodLine1 = baryPartY->findIdsInRange(-eps, +eps);
        MCAuto<DataArrayInt> idsGoodLine2 = baryPartZ->findIdsInRange(-eps, +eps);
        MCAuto<DataArrayInt> idsGoodLine = idsGoodLine1->buildIntersection(idsGoodLine2);
        if (!idsGoodLine->getNumberOfTuples())
          continue;

        // Now the ordering along the Ox axis:
        std::vector<int> insidePoints, hitSegs;
        bool isSplit = OrderPointsAlongLine(mPartCand->_coords->getConstPointer(), nodeMap->begin()[startNode], nodeMap->begin()[endNode],
            mPartCand->getNodalConnectivity()->begin(), mPartCand->getNodalConnectivityIndex()->begin(),
            idsGoodLine->begin(), idsGoodLine->end(),
            /*out*/insidePoints, hitSegs);
        // Optim: smaller segments completly included in eIdx and not split won't need any further treatment:
        for (vector<int>::const_iterator its=hitSegs.begin(); its != hitSegs.end(); ++its)
          hit[cands2[*its]] = true;

        if (!isSplit)  // current segment remains in one piece
          continue;

        // Get original node IDs in global coords array
        for (std::vector<int>::iterator iit = insidePoints.begin(); iit!=insidePoints.end(); ++iit)
          *iit = nodeMapInv->begin()[*iit];

        vector<int> polyIndices, packsIds, facePack;
        // For each face implying this edge
        for (int ii=revDescIP2[eIdx]; ii < revDescIP2[eIdx+1]; ii++)
          {
            int faceIdx = revDescP2[ii];
            // each cell where this face is involved connectivity will be modified:
            ret->pushBackValsSilent(revDescP + revDescIP[faceIdx], revDescP + revDescIP[faceIdx+1]);

            // Current face connectivity
            const int * sIdxConn = cDesc + cIDesc[faceIdx] + 1;
            const int * sIdxConnE = cDesc + cIDesc[faceIdx+1];

            vector<int> modifiedFace;
            ReplaceEdgeInFace(sIdxConn, sIdxConnE, startNode, endNode, insidePoints, /*out*/modifiedFace);
            modifiedFace.insert(modifiedFace.begin(), INTERP_KERNEL::NORM_POLYGON);
            connSlaDesc->replaceSimplePack(faceIdx, modifiedFace.data(), modifiedFace.data()+modifiedFace.size());
          }
      }

    // Rebuild 3D connectivity from descending:
    MCAuto<MEDCouplingSkyLineArray> newConn(MEDCouplingSkyLineArray::New());
    MCAuto<DataArrayInt> superIdx(DataArrayInt::New());  superIdx->alloc(getNumberOfCells()+1); superIdx->fillWithValue(0);
    MCAuto<DataArrayInt> idx(DataArrayInt::New());       idx->alloc(1);                         idx->fillWithValue(0);
    MCAuto<DataArrayInt> vals(DataArrayInt::New());      vals->alloc(0);
    newConn->set3(superIdx, idx, vals);
    for(std::size_t ii = 0; ii < getNumberOfCells(); ii++)
      for (int jj=descIP[ii]; jj < descIP[ii+1]; jj++)
        {
          int sz, faceIdx = abs(descP[jj])-1;
          bool orient = descP[jj]>0;
          const int * p = connSlaDesc->getSimplePackSafePtr(faceIdx, sz);
          if (orient)
            newConn->pushBackPack(ii, p+1, p+sz);  // +1 to skip type
          else
            {
              vector<int> rev(sz-1);
              for (int kk=0; kk<sz-1; kk++) rev[kk]=*(p+sz-kk-1);
              newConn->pushBackPack(ii, rev.data(), rev.data()+sz-1);
            }
        }
    // And finally:
    newConn->convertToPolyhedronConn(cAuto, cIAuto);
  } // end step2

  ret = ret->buildUniqueNotSorted();
  return ret.retn();
}


