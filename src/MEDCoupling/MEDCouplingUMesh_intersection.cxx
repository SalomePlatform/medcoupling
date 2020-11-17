// Copyright (C) 2007-2020  CEA/DEN, EDF R&D
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

mcIdType InternalAddPoint(const INTERP_KERNEL::Edge *e, mcIdType id, const double *coo, mcIdType startId, mcIdType endId, DataArrayDouble& addCoo, mcIdType& nodesCnter)
{
  if(id!=-1)
    return id;
  else
    {
      mcIdType ret(nodesCnter++);
      double newPt[2];
      e->getMiddleOfPoints(coo+2*startId,coo+2*endId,newPt);
      addCoo.insertAtTheEnd(newPt,newPt+2);
      return ret;
    }
}

mcIdType InternalAddPointOriented(const INTERP_KERNEL::Edge *e, mcIdType id, const double *coo, mcIdType startId, mcIdType endId, DataArrayDouble& addCoo, mcIdType& nodesCnter)
{
  if(id!=-1)
    return id;
  else
    {
      mcIdType ret(nodesCnter++);
      double newPt[2];
      e->getMiddleOfPointsOriented(coo+2*startId,coo+2*endId,newPt);
      addCoo.insertAtTheEnd(newPt,newPt+2);
      return ret;
    }
}


void EnterTheResultOf2DCellFirst(const INTERP_KERNEL::Edge *e, int start, int stp, int nbOfEdges, bool linOrArc, const double *coords, const mcIdType *connBg, mcIdType offset, DataArrayIdType *newConnOfCell, DataArrayDouble *appendedCoords, std::vector<mcIdType>& middles)
{
  mcIdType tmp[3];
  int trueStart(start>=0?start:nbOfEdges+start);
  tmp[0]=ToIdType(linOrArc?INTERP_KERNEL::NORM_QPOLYG:INTERP_KERNEL::NORM_POLYGON); tmp[1]=connBg[trueStart]; tmp[2]=connBg[stp];
  newConnOfCell->insertAtTheEnd(tmp,tmp+3);
  if(linOrArc)
    {
      if(stp-start>1)
        {
          mcIdType tmp2(0),tmp3(appendedCoords->getNumberOfTuples()/2);
          InternalAddPointOriented(e,-1,coords,tmp[1],tmp[2],*appendedCoords,tmp2);
          middles.push_back(tmp3+offset);
        }
      else
        middles.push_back(connBg[trueStart+nbOfEdges]);
    }
}

void EnterTheResultOf2DCellMiddle(const INTERP_KERNEL::Edge *e, int start, int stp, int nbOfEdges, bool linOrArc, const double *coords, const mcIdType *connBg, mcIdType offset, DataArrayIdType *newConnOfCell, DataArrayDouble *appendedCoords, std::vector<mcIdType>& middles)
{
  mcIdType tmpSrt(newConnOfCell->back()),tmpEnd(connBg[stp]);
  newConnOfCell->pushBackSilent(tmpEnd);
  if(linOrArc)
    {
      if(stp-start>1)
        {
          mcIdType tmp2(0),tmp3(appendedCoords->getNumberOfTuples()/2);
          InternalAddPointOriented(e,-1,coords,tmpSrt,tmpEnd,*appendedCoords,tmp2);
          middles.push_back(tmp3+offset);
        }
      else
        middles.push_back(connBg[start+nbOfEdges]);
    }
}

void EnterTheResultOf2DCellEnd(const INTERP_KERNEL::Edge *e, int start, int stp, int nbOfEdges, bool linOrArc, const double *coords, const mcIdType *connBg, mcIdType offset, DataArrayIdType *newConnOfCell, DataArrayDouble *appendedCoords, std::vector<mcIdType>& middles)
{
  // only the quadratic point to deal with:
  if(linOrArc)
    {
      if(stp-start>1)  // if we are covering more than one segment we need to create a new mid point
        {
          mcIdType tmpSrt(connBg[start]),tmpEnd(connBg[stp % nbOfEdges]);  // % to handle last seg.
          mcIdType tmp2(0),tmp3(appendedCoords->getNumberOfTuples()/2);
          InternalAddPointOriented(e,-1,coords,tmpSrt,tmpEnd,*appendedCoords,tmp2);
          middles.push_back(tmp3+offset);
        }
      else
        middles.push_back(connBg[start+nbOfEdges]);
    }
}

void IKGeo2DInternalMapper2(INTERP_KERNEL::Node *n, const std::map<MCAuto<INTERP_KERNEL::Node>,mcIdType>& m, mcIdType forbVal0, mcIdType forbVal1, std::vector<mcIdType>& isect)
{
  MCAuto<INTERP_KERNEL::Node> nTmp(n); nTmp->incrRef();
  std::map<MCAuto<INTERP_KERNEL::Node>,mcIdType>::const_iterator it(m.find(nTmp));
  if(it==m.end())
    throw INTERP_KERNEL::Exception("Internal error in remapping !");
  mcIdType v((*it).second);
  if(v==forbVal0 || v==forbVal1)
    return ;
  if(std::find(isect.begin(),isect.end(),v)==isect.end())
    isect.push_back(v);
}

bool IKGeo2DInternalMapper(const INTERP_KERNEL::ComposedEdge& c, const std::map<MCAuto<INTERP_KERNEL::Node>,mcIdType>& m, mcIdType forbVal0, mcIdType forbVal1, std::vector<mcIdType>& isect)
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

  INTERP_KERNEL::Edge *MEDCouplingUMeshBuildQPFromEdge2(INTERP_KERNEL::NormalizedCellType typ, const mcIdType *bg, const double *coords2D, std::map< MCAuto<INTERP_KERNEL::Node>,mcIdType>& m)
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

  INTERP_KERNEL::Edge *MEDCouplingUMeshBuildQPFromEdge(INTERP_KERNEL::NormalizedCellType typ, std::map<mcIdType, INTERP_KERNEL::NodeWithUsage >& mapp2, const mcIdType *bg)
  {
    INTERP_KERNEL::Edge *ret=0;

    mapp2[bg[0]].second = INTERP_KERNEL::USAGE_LINEAR;
    mapp2[bg[1]].second = INTERP_KERNEL::USAGE_LINEAR;

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
          if (mapp2[bg[2]].second != INTERP_KERNEL::USAGE_LINEAR) // switch the node usage to quadratic only if it is not used as an extreme point for another edge
            mapp2[bg[2]].second = INTERP_KERNEL::USAGE_QUADRATIC_ONLY;
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
  INTERP_KERNEL::QuadraticPolygon *MEDCouplingUMeshBuildQPFromMesh(const MEDCouplingUMesh *mDesc, const std::vector<mcIdType>& candidates,
                                                                   std::map<INTERP_KERNEL::Node *,mcIdType>& mapp)
  {
    mapp.clear();
    std::map<mcIdType, INTERP_KERNEL::NodeWithUsage > mapp2;  // the last var is a flag specifying if node is an extreme node of the seg (LINEAR) or only a middle for SEG3 (QUADRATIC_ONLY).
    const double *coo=mDesc->getCoords()->getConstPointer();
    const mcIdType *c=mDesc->getNodalConnectivity()->getConstPointer();
    const mcIdType *cI=mDesc->getNodalConnectivityIndex()->getConstPointer();
    std::set<mcIdType> s;
    for(std::vector<mcIdType>::const_iterator it=candidates.begin();it!=candidates.end();it++)
      s.insert(c+cI[*it]+1,c+cI[(*it)+1]);
    for(std::set<mcIdType>::const_iterator it2=s.begin();it2!=s.end();it2++)
      {
        INTERP_KERNEL::Node *n=new INTERP_KERNEL::Node(coo[2*(*it2)],coo[2*(*it2)+1]);
        mapp2[*it2]=INTERP_KERNEL::NodeWithUsage(n,INTERP_KERNEL::USAGE_UNKNOWN);
      }
    INTERP_KERNEL::QuadraticPolygon *ret=new INTERP_KERNEL::QuadraticPolygon;
    for(std::vector<mcIdType>::const_iterator it=candidates.begin();it!=candidates.end();it++)
      {
        INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)c[cI[*it]];
        ret->pushBack(MEDCouplingUMeshBuildQPFromEdge(typ,mapp2,c+cI[*it]+1));
      }
    for(std::map<mcIdType, INTERP_KERNEL::NodeWithUsage >::const_iterator it2=mapp2.begin();it2!=mapp2.end();it2++)
      {
        if((*it2).second.second == INTERP_KERNEL::USAGE_LINEAR)
          mapp[(*it2).second.first]=(*it2).first;
        ((*it2).second.first)->decrRef();
      }
    return ret;
  }

  INTERP_KERNEL::QuadraticPolygon *MEDCouplingUMeshBuildQPFromMeshWithTree(const MEDCouplingUMesh *mDesc, const std::vector<mcIdType>& candidates,
                                                                   std::map<INTERP_KERNEL::Node *,mcIdType>& mapp,
                                                                   const BBTreePts<2,mcIdType> & nodeTree,
                                                                   const std::map<mcIdType, INTERP_KERNEL::Node *>& mapRev)
  {
    mapp.clear();
    std::map<mcIdType, INTERP_KERNEL::NodeWithUsage > mapp2;  // the last var is a flag specifying if node is an extreme node of the seg (LINEAR) or only a middle for SEG3 (QUADRATIC_ONLY).
    const double *coo=mDesc->getCoords()->getConstPointer();
    const mcIdType *c=mDesc->getNodalConnectivity()->getConstPointer();
    const mcIdType *cI=mDesc->getNodalConnectivityIndex()->getConstPointer();
    std::set<mcIdType> s;
    for(std::vector<mcIdType>::const_iterator it=candidates.begin();it!=candidates.end();it++)
      s.insert(c+cI[*it]+1,c+cI[(*it)+1]);
    for(std::set<mcIdType>::const_iterator it2=s.begin();it2!=s.end();it2++)
      {
        INTERP_KERNEL::Node *n;
        // Look for a potential node to merge
        std::vector<mcIdType> candNode;
        nodeTree.getElementsAroundPoint(coo+2*(*it2), candNode);
        if (candNode.size() > 2)
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::MEDCouplingUMeshBuildQPFromMeshWithTree(): some nodes are not properly merged (within eps) in input mesh!");
        bool node_created=false;
        if  (candNode.size())
          {
            auto itt=mapRev.find(candNode[0]);
            if (itt != mapRev.end())  // we might hit a node which is in the coords array but not used in the connectivity in which case it won't be in the revMap
              {
                node_created=true;
                n = (*itt).second;
                n->incrRef();
              }
          }
        if(!node_created)
          n = new INTERP_KERNEL::Node(coo[2*(*it2)],coo[2*(*it2)+1]);
        mapp2[*it2]=INTERP_KERNEL::NodeWithUsage(n,INTERP_KERNEL::USAGE_UNKNOWN);
      }
    INTERP_KERNEL::QuadraticPolygon *ret=new INTERP_KERNEL::QuadraticPolygon;
    for(std::vector<mcIdType>::const_iterator it=candidates.begin();it!=candidates.end();it++)
      {
        INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)c[cI[*it]];
        ret->pushBack(MEDCouplingUMeshBuildQPFromEdge(typ,mapp2,c+cI[*it]+1));  // this call will set quad points to false in the map
      }
    for(std::map<mcIdType, INTERP_KERNEL::NodeWithUsage >::const_iterator it2=mapp2.begin();it2!=mapp2.end();it2++)
      {
        if((*it2).second.second == INTERP_KERNEL::USAGE_LINEAR)
          mapp[(*it2).second.first]=(*it2).first;
        ((*it2).second.first)->decrRef();
      }
    return ret;
  }

  INTERP_KERNEL::Node *MEDCouplingUMeshBuildQPNode(mcIdType nodeId, const double *coo1, mcIdType offset1, const double *coo2, mcIdType offset2, const std::vector<double>& addCoo)
  {
    if(nodeId>=offset2)
      {
        mcIdType locId=nodeId-offset2;
        return new INTERP_KERNEL::Node(addCoo[2*locId],addCoo[2*locId+1]);
      }
    if(nodeId>=offset1)
      {
        mcIdType locId=nodeId-offset1;
        return new INTERP_KERNEL::Node(coo2[2*locId],coo2[2*locId+1]);
      }
    return new INTERP_KERNEL::Node(coo1[2*nodeId],coo1[2*nodeId+1]);
  }

  /**
   * Construct a mapping between set of Nodes and the standard MEDCoupling connectivity format (c, cI).
   */
  void MEDCouplingUMeshBuildQPFromMesh3(const double *coo1, mcIdType offset1, const double *coo2, mcIdType offset2, const std::vector<double>& addCoo,
                                        const mcIdType *desc1Bg, const mcIdType *desc1End, const std::vector<std::vector<mcIdType> >& intesctEdges1,
                                        /*output*/std::map<INTERP_KERNEL::Node *,mcIdType>& mapp, std::map<mcIdType,INTERP_KERNEL::Node *>& mappRev)
  {
    for(const mcIdType *desc1=desc1Bg;desc1!=desc1End;desc1++)
      {
        mcIdType eltId1=std::abs(*desc1)-1;
        for(std::vector<mcIdType>::const_iterator it1=intesctEdges1[eltId1].begin();it1!=intesctEdges1[eltId1].end();it1++)
          {
            std::map<mcIdType,INTERP_KERNEL::Node *>::const_iterator it=mappRev.find(*it1);
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
 * \param forbiddenPoints the list of points that should not be removed in the process
 */
bool MEDCouplingUMesh::Colinearize2DCell(const double *coords, const mcIdType *connBg, const mcIdType *connEnd, mcIdType offset,
                                         const std::map<mcIdType, bool>& forbiddenPoints,
                                         DataArrayIdType *newConnOfCell, DataArrayDouble *appendedCoords)
{
  std::size_t sz(std::distance(connBg,connEnd));
  if(sz<3)//3 because 2+1(for the cell type) and 2 is the minimal number of edges of 2D cell.
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Colinearize2DCell : the input cell has invalid format !");
  sz--;
  INTERP_KERNEL::AutoPtr<mcIdType> tmpConn(new mcIdType[sz]);
  INTERP_KERNEL::AutoPtr<mcIdType> tmpConn2(new mcIdType[sz]);
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)connBg[0]));
  unsigned nbs(cm.getNumberOfSons2(connBg+1,ToIdType(sz)));
  unsigned nbOfHit(0); // number of fusions operated
  int posBaseElt(0),posEndElt(0),nbOfTurn(0);
  const std::size_t maxNbOfHit = cm.isQuadratic() ? nbs-2 : nbs-3;  // a quad cell is authorized to end up with only two edges, a linear one has to keep 3 at least
  INTERP_KERNEL::NormalizedCellType typeOfSon;
  std::vector<mcIdType> middles;
  bool ret(false);
  for(;(nbOfTurn+nbOfHit)<nbs;nbOfTurn++)
    {
      cm.fillSonCellNodalConnectivity2(posBaseElt,connBg+1,ToIdType(sz),tmpConn,typeOfSon);
      std::map<MCAuto<INTERP_KERNEL::Node>,mcIdType> m;
      INTERP_KERNEL::Edge *e(MEDCouplingUMeshBuildQPFromEdge2(typeOfSon,tmpConn,coords,m));
      posEndElt = posBaseElt+1;

      // Look backward first: are the final edges of the cells colinear with the first ones?
      // This initializes posBaseElt.
      if(nbOfTurn==0)
        {
          for(unsigned i=1;i<nbs && nbOfHit<maxNbOfHit;i++) // 2nd condition is to avoid ending with a cell with one single edge
            {
              cm.fillSonCellNodalConnectivity2(nbs-i,connBg+1,ToIdType(sz),tmpConn2,typeOfSon);
              // Identify common point:
              mcIdType commPoint = std::find((mcIdType *)tmpConn, tmpConn+2, tmpConn2[0]) != tmpConn+2 ? tmpConn2[0] : tmpConn2[1];
              auto itE(forbiddenPoints.end());
              if (forbiddenPoints.find(commPoint) != itE) // is the junction point in the list of points we can not remove?
                break;
              INTERP_KERNEL::Edge *eCand(MEDCouplingUMeshBuildQPFromEdge2(typeOfSon,tmpConn2,coords,m));
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
              // Update last connectivity
              std::copy((mcIdType *)tmpConn2, tmpConn2+sz, (mcIdType *)tmpConn);
            }
        }
      // Now move forward:
      const unsigned fwdStart = (nbOfTurn == 0 ? 0 : posBaseElt);  // the first element to be inspected going forward
      for(unsigned j=fwdStart+1;j<nbs && nbOfHit<maxNbOfHit;j++)  // 2nd condition is to avoid ending with a cell with one single edge
        {
          cm.fillSonCellNodalConnectivity2(j,connBg+1,ToIdType(sz),tmpConn2,typeOfSon); // get edge #j's connectivity
          // Identify common point:
          mcIdType commPoint = std::find((mcIdType *)tmpConn, tmpConn+2, tmpConn2[0]) != tmpConn+2 ? tmpConn2[0] : tmpConn2[1];
          auto itE(forbiddenPoints.end());
          if (forbiddenPoints.find(commPoint) != itE) // is the junction point in the list of points we can not remove?
            break;
          INTERP_KERNEL::Edge *eCand(MEDCouplingUMeshBuildQPFromEdge2(typeOfSon,tmpConn2,coords,m));
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
          // Update last connectivity
          std::copy((mcIdType *)tmpConn2, tmpConn2+sz, (mcIdType *)tmpConn);
        }
      //push [posBaseElt,posEndElt) in newConnOfCell using e
      // The if clauses below are (voluntary) not mutually exclusive: on a quad cell with 2 edges, the end of the connectivity is also its beginning!
      if(nbOfTurn==0)
        // at the beginning of the connectivity (insert type)
        EnterTheResultOf2DCellFirst(e,posBaseElt,posEndElt,nbs,cm.isQuadratic(),coords,connBg+1,offset,newConnOfCell,appendedCoords,middles);
      else if((nbOfHit+nbOfTurn) != (nbs-1))
        // in the middle
        EnterTheResultOf2DCellMiddle(e,posBaseElt,posEndElt,nbs,cm.isQuadratic(),coords,connBg+1,offset,newConnOfCell,appendedCoords,middles);
      if ((nbOfHit+nbOfTurn) == (nbs-1))
        // at the end (only quad points to deal with)
        EnterTheResultOf2DCellEnd(e,posBaseElt,posEndElt,nbs,cm.isQuadratic(),coords,connBg+1,offset,newConnOfCell,appendedCoords,middles);
      posBaseElt=posEndElt;
      e->decrRef();
    }
  if(!middles.empty())
    newConnOfCell->insertAtTheEnd(middles.begin(),middles.end());
  return ret;
}



bool IsColinearOfACellOf(const std::vector< std::vector<mcIdType> >& intersectEdge1, const std::vector<mcIdType>& candidates, mcIdType start, mcIdType stop, mcIdType& retVal)
{
  if(candidates.empty())
    return false;
  for(std::vector<mcIdType>::const_iterator it=candidates.begin();it!=candidates.end();it++)
    {
      const std::vector<mcIdType>& pool(intersectEdge1[*it]);
      mcIdType tmp[2]; tmp[0]=start; tmp[1]=stop;
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
 *  - subDiv of size 'm2->getNumberOfCells()' that lists for each seg cell in 'm2' the splitting node ids randomly sorted.
 * The aim of this method is to sort the splitting nodes, if any, and to put them in 'intersectEdge' output parameter based on edges of mesh 'm2'
 * Nodes end up lying consecutively on a cutted edge.
 * \param m1 is expected to be a mesh of meshDimension equal to 1 and spaceDim equal to 2. No check of that is performed by this method.
 * (Only present for its coords in case of 'subDiv' shares some nodes of 'm1')
 * \param m2 is expected to be a mesh of meshDimension equal to 1 and spaceDim equal to 2. No check of that is performed by this method.
 * \param addCoo input parameter with additional nodes linked to intersection of the 2 meshes.
 * \param[out] intersectEdge the same content as subDiv, but correctly oriented.
 */
void MEDCouplingUMesh::BuildIntersectEdges(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2,
                                           const std::vector<double>& addCoo,
                                           const std::vector< std::vector<mcIdType> >& subDiv, std::vector< std::vector<mcIdType> >& intersectEdge)
{
  mcIdType offset1=m1->getNumberOfNodes();
  mcIdType ncell2=m2->getNumberOfCells();
  const mcIdType *c=m2->getNodalConnectivity()->begin();
  const mcIdType *cI=m2->getNodalConnectivityIndex()->begin();
  const double *coo=m2->getCoords()->begin();
  const double *cooBis=m1->getCoords()->begin();
  mcIdType offset2=offset1+m2->getNumberOfNodes();
  intersectEdge.resize(ncell2);
  for(mcIdType i=0;i<ncell2;i++,cI++)
    {
      const std::vector<mcIdType>& divs=subDiv[i];
      mcIdType nnode=cI[1]-cI[0]-1;
      std::map<mcIdType, INTERP_KERNEL::NodeWithUsage > mapp2;
      std::map<INTERP_KERNEL::Node *, mcIdType> mapp22;
      for(mcIdType j=0;j<nnode;j++)
        {
          INTERP_KERNEL::Node *nn=new INTERP_KERNEL::Node(coo[2*c[(*cI)+j+1]],coo[2*c[(*cI)+j+1]+1]);
          mcIdType nnid=c[(*cI)+j+1];
          mapp2[nnid]=INTERP_KERNEL::NodeWithUsage(nn,INTERP_KERNEL::USAGE_UNKNOWN);
          mapp22[nn]=nnid+offset1;
        }
      INTERP_KERNEL::Edge *e=MEDCouplingUMeshBuildQPFromEdge((INTERP_KERNEL::NormalizedCellType)c[*cI],mapp2,c+(*cI)+1);
      for(std::map<mcIdType, INTERP_KERNEL::NodeWithUsage >::const_iterator it=mapp2.begin();it!=mapp2.end();it++)
        ((*it).second.first)->decrRef();
      std::vector<INTERP_KERNEL::Node *> addNodes(divs.size());
      std::map<INTERP_KERNEL::Node *,mcIdType> mapp3;
      for(std::size_t j=0;j<divs.size();j++)
        {
          mcIdType id=divs[j];
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

/*
 *  Build the final 1D mesh resulting from the newly created points after intersection with the segments of the descending 2D mesh.
 *  @param[out] idsInRetColinear IDs of edges in the result (ret) that are colinears to one of the segment of the descending 2D mesh. Indexing scheme
 *  is the one of the ret 1D mesh.
 *  @param[out] idsInMesh1DForIdsInRetColinear same IDs as above in the descending 2D mesh
 */
MEDCouplingUMesh *BuildMesh1DCutFrom(const MEDCouplingUMesh *mesh1D, const std::vector< std::vector<mcIdType> >& intersectEdge2,
                                     const DataArrayDouble *coords1, const std::vector<double>& addCoo, const std::map<mcIdType,mcIdType>& mergedNodes,
                                     const std::vector< std::vector<mcIdType> >& colinear2, const std::vector< std::vector<mcIdType> >& intersectEdge1,
                                     MCAuto<DataArrayIdType>& idsInRetColinear, MCAuto<DataArrayIdType>& idsInMesh1DForIdsInRetColinear)
{
  idsInRetColinear=DataArrayIdType::New(); idsInRetColinear->alloc(0,1);
  idsInMesh1DForIdsInRetColinear=DataArrayIdType::New(); idsInMesh1DForIdsInRetColinear->alloc(0,1);
  mcIdType nCells=mesh1D->getNumberOfCells();
  if(nCells!=ToIdType(intersectEdge2.size()))
    throw INTERP_KERNEL::Exception("BuildMesh1DCutFrom : internal error # 1 !");
  const DataArrayDouble *coo2(mesh1D->getCoords());
  const mcIdType *c(mesh1D->getNodalConnectivity()->begin()),*ci(mesh1D->getNodalConnectivityIndex()->begin());
  const double *coo2Ptr(coo2->begin());
  mcIdType offset1(coords1->getNumberOfTuples());
  mcIdType offset2(offset1+coo2->getNumberOfTuples());
  mcIdType offset3(offset2+ToIdType(addCoo.size())/2);
  std::vector<double> addCooQuad;
  MCAuto<DataArrayIdType> cOut(DataArrayIdType::New()),ciOut(DataArrayIdType::New()); cOut->alloc(0,1); ciOut->alloc(1,1); ciOut->setIJ(0,0,0);
  mcIdType tmp[4],cicnt(0),kk(0);
  for(mcIdType i=0;i<nCells;i++)
    {
      std::map<MCAuto<INTERP_KERNEL::Node>,mcIdType> m;
      INTERP_KERNEL::Edge *e(MEDCouplingUMeshBuildQPFromEdge2((INTERP_KERNEL::NormalizedCellType)c[ci[i]],c+ci[i]+1,coo2Ptr,m));
      const std::vector<mcIdType>& subEdges(intersectEdge2[i]);
      mcIdType nbSubEdge=ToIdType(subEdges.size()/2);
      for(mcIdType j=0;j<nbSubEdge;j++,kk++)
        {
          MCAuto<INTERP_KERNEL::Node> n1(MEDCouplingUMeshBuildQPNode(subEdges[2*j],coords1->begin(),offset1,coo2Ptr,offset2,addCoo)),
                                      n2(MEDCouplingUMeshBuildQPNode(subEdges[2*j+1],coords1->begin(),offset1,coo2Ptr,offset2,addCoo));
          MCAuto<INTERP_KERNEL::Edge> e2(e->buildEdgeLyingOnMe(n1,n2));
          INTERP_KERNEL::Edge *e2Ptr(e2);
          std::map<mcIdType,mcIdType>::const_iterator itm;
          if(dynamic_cast<INTERP_KERNEL::EdgeArcCircle *>(e2Ptr))
            {
              tmp[0]=INTERP_KERNEL::NORM_SEG3;
              itm=mergedNodes.find(subEdges[2*j]);
              tmp[1]=itm!=mergedNodes.end()?(*itm).second:subEdges[2*j];
              itm=mergedNodes.find(subEdges[2*j+1]);
              tmp[2]=itm!=mergedNodes.end()?(*itm).second:subEdges[2*j+1];
              tmp[3]=offset3+ToIdType(addCooQuad.size()/2);
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
          mcIdType tmp00;
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
  arr3->useArray(&addCoo[0],false,DeallocType::C_DEALLOC,addCoo.size()/2,2);
  MCAuto<DataArrayDouble> arr4(DataArrayDouble::New()); arr4->useArray(&addCooQuad[0],false,DeallocType::C_DEALLOC,addCooQuad.size()/2,2);
  std::vector<const DataArrayDouble *> coordss(4);
  coordss[0]=coords1; coordss[1]=mesh1D->getCoords(); coordss[2]=arr3; coordss[3]=arr4;
  MCAuto<DataArrayDouble> arr(DataArrayDouble::Aggregate(coordss));
  ret->setCoords(arr);
  return ret.retn();
}

MEDCouplingUMesh *BuildRefined2DCellLinear(const DataArrayDouble *coords, const mcIdType *descBg, const mcIdType *descEnd, const std::vector< std::vector<mcIdType> >& intersectEdge1)
{
  std::vector<mcIdType> allEdges;
  for(const mcIdType *it2(descBg);it2!=descEnd;it2++)
    {
      const std::vector<mcIdType>& edge1(intersectEdge1[std::abs(*it2)-1]);
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
  std::vector<mcIdType> connOut(nbOfEdgesOf2DCellSplit);
  for(std::size_t kk=0;kk<nbOfEdgesOf2DCellSplit;kk++)
    connOut[kk]=allEdges[2*kk];
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYGON,ToIdType(connOut.size()),&connOut[0]);
  return ret.retn();
}

MEDCouplingUMesh *BuildRefined2DCellQuadratic(const DataArrayDouble *coords, const MEDCouplingUMesh *mesh2D, mcIdType cellIdInMesh2D, const mcIdType *descBg, const mcIdType *descEnd, const std::vector< std::vector<mcIdType> >& intersectEdge1)
{
  const mcIdType *c(mesh2D->getNodalConnectivity()->begin()),*ci(mesh2D->getNodalConnectivityIndex()->begin());
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)c[ci[cellIdInMesh2D]]));
  int ii(0);
  unsigned sz(cm.getNumberOfSons2(c+ci[cellIdInMesh2D]+1,ci[cellIdInMesh2D+1]-ci[cellIdInMesh2D]-1));
  if(sz!=std::distance(descBg,descEnd))
    throw INTERP_KERNEL::Exception("BuildRefined2DCellQuadratic : internal error 1 !");
  INTERP_KERNEL::AutoPtr<mcIdType> tmpPtr(new mcIdType[ci[cellIdInMesh2D+1]-ci[cellIdInMesh2D]]);
  std::vector<mcIdType> allEdges,centers;
  const double *coordsPtr(coords->begin());
  MCAuto<DataArrayDouble> addCoo(DataArrayDouble::New()); addCoo->alloc(0,1);
  mcIdType offset(coords->getNumberOfTuples());
  for(const mcIdType *it2(descBg);it2!=descEnd;it2++,ii++)
    {
      INTERP_KERNEL::NormalizedCellType typeOfSon;
      cm.fillSonCellNodalConnectivity2(ii,c+ci[cellIdInMesh2D]+1,ci[cellIdInMesh2D+1]-ci[cellIdInMesh2D]-1,tmpPtr,typeOfSon);
      const std::vector<mcIdType>& edge1(intersectEdge1[std::abs(*it2)-1]);
      if(*it2>0)
        allEdges.insert(allEdges.end(),edge1.begin(),edge1.end());
      else
        allEdges.insert(allEdges.end(),edge1.rbegin(),edge1.rend());
      if(edge1.size()==2)
        centers.push_back(tmpPtr[2]);//special case where no subsplit of edge -> reuse the original center.
      else
        {//the current edge has been subsplit -> create corresponding centers.
          mcIdType nbOfCentersToAppend=ToIdType(edge1.size()/2);
          std::map< MCAuto<INTERP_KERNEL::Node>,mcIdType> m;
          MCAuto<INTERP_KERNEL::Edge> ee(MEDCouplingUMeshBuildQPFromEdge2(typeOfSon,tmpPtr,coordsPtr,m));
          std::vector<mcIdType>::const_iterator it3(allEdges.end()-edge1.size());
          for(mcIdType k=0;k<nbOfCentersToAppend;k++)
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
  std::vector<mcIdType> connOut(nbOfEdgesOf2DCellSplit);
  for(std::size_t kk=0;kk<nbOfEdgesOf2DCellSplit;kk++)
    connOut[kk]=allEdges[2*kk];
  connOut.insert(connOut.end(),centers.begin(),centers.end());
  ret->insertNextCell(INTERP_KERNEL::NORM_QPOLYG,ToIdType(connOut.size()),&connOut[0]);
  return ret.retn();
}

/*!
 * This method creates a refinement of a cell in \a mesh2D. Those cell is defined by descending connectivity and the sorted subdivided nodal connectivity
 * of those edges.
 *
 * \param [in] mesh2D - The origin 2D mesh. \b Warning \b coords are not those of \a mesh2D. But mesh2D->getCoords()==coords[:mesh2D->getNumberOfNodes()]
 */
MEDCouplingUMesh *BuildRefined2DCell(const DataArrayDouble *coords, const MEDCouplingUMesh *mesh2D, mcIdType cellIdInMesh2D, const mcIdType *descBg, const mcIdType *descEnd, const std::vector< std::vector<mcIdType> >& intersectEdge1)
{
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(mesh2D->getTypeOfCell(cellIdInMesh2D)));
  if(!cm.isQuadratic())
    return BuildRefined2DCellLinear(coords,descBg,descEnd,intersectEdge1);
  else
    return BuildRefined2DCellQuadratic(coords,mesh2D,cellIdInMesh2D,descBg,descEnd,intersectEdge1);
}

void AddCellInMesh2D(MEDCouplingUMesh *mesh2D, const std::vector<mcIdType>& conn, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& edges)
{
  bool isQuad(false);
  for(std::vector< MCAuto<INTERP_KERNEL::Edge> >::const_iterator it=edges.begin();it!=edges.end();it++)
    {
      const INTERP_KERNEL::Edge *ee(*it);
      if(dynamic_cast<const INTERP_KERNEL::EdgeArcCircle *>(ee))
        isQuad=true;
    }
  if(!isQuad)
    mesh2D->insertNextCell(INTERP_KERNEL::NORM_POLYGON,ToIdType(conn.size()),&conn[0]);
  else
    {
      const double *coo(mesh2D->getCoords()->begin());
      std::size_t sz(conn.size());
      std::vector<double> addCoo;
      std::vector<mcIdType> conn2(conn);
      mcIdType offset(mesh2D->getNumberOfNodes());
      for(std::size_t i=0;i<sz;i++)
        {
          double tmp[2];
          edges[(i+1)%sz]->getMiddleOfPoints(coo+2*conn[i],coo+2*conn[(i+1)%sz],tmp);// tony a chier i+1 -> i
          addCoo.insert(addCoo.end(),tmp,tmp+2);
          conn2.push_back(offset+ToIdType(i));
        }
      mesh2D->getCoords()->rearrange(1);
      mesh2D->getCoords()->pushBackValsSilent(&addCoo[0],&addCoo[0]+addCoo.size());
      mesh2D->getCoords()->rearrange(2);
      mesh2D->insertNextCell(INTERP_KERNEL::NORM_QPOLYG,ToIdType(conn2.size()),&conn2[0]);
    }
}

/*!
 * \b WARNING edges in out1 coming from \a splitMesh1D are \b NOT oriented because only used for equation of curve.
 *
 * This method cuts in 2 parts the input 2D cell given using boundaries description (\a edge1Bis and \a edge1BisPtr) using
 * a set of edges defined in \a splitMesh1D.
 */
void BuildMesh2DCutInternal2(const MEDCouplingUMesh *splitMesh1D, const std::vector<mcIdType>& edge1Bis, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& edge1BisPtr,
                             std::vector< std::vector<mcIdType> >& out0, std::vector< std::vector< MCAuto<INTERP_KERNEL::Edge> > >& out1)
{
  std::size_t nb(edge1Bis.size()/2);
  std::size_t nbOfEdgesOf2DCellSplit(nb/2);
  mcIdType iEnd=splitMesh1D->getNumberOfCells();
  if(iEnd==0)
    throw INTERP_KERNEL::Exception("BuildMesh2DCutInternal2 : internal error ! input 1D mesh must have at least one cell !");
  std::size_t ii,jj;
  const mcIdType *cSplitPtr(splitMesh1D->getNodalConnectivity()->begin()),*ciSplitPtr(splitMesh1D->getNodalConnectivityIndex()->begin());
  for(ii=0;ii<nb && edge1Bis[2*ii]!=cSplitPtr[ciSplitPtr[0]+1];ii++);
  for(jj=ii;jj<nb && edge1Bis[2*jj+1]!=cSplitPtr[ciSplitPtr[iEnd-1]+2];jj++);
  //
  if(jj==nb)
    {//the edges splitMesh1D[iStart:iEnd] does not fully cut the current 2D cell -> single output cell
      out0.resize(1); out1.resize(1);
      std::vector<mcIdType>& connOut(out0[0]);
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
      std::vector<mcIdType>& connOutLeft(out0[0]);
      std::vector<mcIdType>& connOutRight(out0[1]);//connOutLeft should end with edge1Bis[2*ii] and connOutRight should end with edge1Bis[2*jj+1]
      std::vector< MCAuto<INTERP_KERNEL::Edge> >& eleft(out1[0]);
      std::vector< MCAuto<INTERP_KERNEL::Edge> >& eright(out1[1]);
      for(std::size_t k=ii;k<jj+1;k++)
        { connOutLeft.push_back(edge1Bis[2*k+1]); eleft.push_back(edge1BisPtr[2*k+1]); }
      std::vector< MCAuto<INTERP_KERNEL::Edge> > ees(iEnd);
      for(mcIdType ik=0;ik<iEnd;ik++)
        {
          std::map< MCAuto<INTERP_KERNEL::Node>,mcIdType> m;
          MCAuto<INTERP_KERNEL::Edge> ee(MEDCouplingUMeshBuildQPFromEdge2((INTERP_KERNEL::NormalizedCellType)cSplitPtr[ciSplitPtr[ik]],cSplitPtr+ciSplitPtr[ik]+1,splitMesh1D->getCoords()->begin(),m));
          ees[ik]=ee;
        }
      for(mcIdType ik=iEnd-1;ik>=0;ik--)
        connOutLeft.push_back(cSplitPtr[ciSplitPtr[ik]+1]);
      for(std::size_t k=jj+1;k<nbOfEdgesOf2DCellSplit+ii;k++)
        { connOutRight.push_back(edge1Bis[2*k+1]); eright.push_back(edge1BisPtr[2*k+1]); }
      eleft.insert(eleft.end(),ees.rbegin(),ees.rend());
      for(mcIdType ik=0;ik<iEnd;ik++)
        connOutRight.push_back(cSplitPtr[ciSplitPtr[ik]+2]);
      eright.insert(eright.end(),ees.begin(),ees.end());
    }
}

struct CellInfo
{
public:
  CellInfo() { }
  CellInfo(const std::vector<mcIdType>& edges, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& edgesPtr);
public:
  std::vector<mcIdType> _edges;
  std::vector< MCAuto<INTERP_KERNEL::Edge> > _edges_ptr;
};

CellInfo::CellInfo(const std::vector<mcIdType>& edges, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& edgesPtr)
{
  std::size_t nbe(edges.size());
  std::vector<mcIdType> edges2(2*nbe); std::vector< MCAuto<INTERP_KERNEL::Edge> > edgesPtr2(2*nbe);
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
  EdgeInfo(mcIdType istart, mcIdType iend, const MCAuto<MEDCouplingUMesh>& mesh):_istart(istart),_iend(iend),_mesh(mesh),_left(-7),_right(-7) { }
  EdgeInfo(mcIdType istart, mcIdType iend, mcIdType pos, const MCAuto<INTERP_KERNEL::Edge>& edge):_istart(istart),_iend(iend),_edge(edge),_left(pos),_right(pos+1) { }
  bool isInMyRange(mcIdType pos) const { return pos>=_istart && pos<_iend; }
  void somethingHappendAt(mcIdType pos, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newLeft, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newRight);
  void feedEdgeInfoAt(double eps, const MEDCouplingUMesh *mesh2D, mcIdType offset, mcIdType neighbors[2]) const;
private:
  mcIdType _istart;
  mcIdType _iend;
  MCAuto<MEDCouplingUMesh> _mesh;
  MCAuto<INTERP_KERNEL::Edge> _edge;
  mcIdType _left;    // index (local numbering) of the left 2D cell bordering the edge '_edge'
  mcIdType _right;   // same as above, right side.
};

/*
 * Update indices of left and right 2D cell bordering the current edge.
 */
void EdgeInfo::somethingHappendAt(mcIdType pos, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newLeft, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newRight)
{
  const MEDCouplingUMesh *mesh(_mesh);
  if(mesh)
    return ;
  if(_right<pos)
    return ;
  if(_left>pos)
    { _left++; _right++; return ; }
  if (_right > pos && _left != pos)
    { _right++; return ; }
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

void EdgeInfo::feedEdgeInfoAt(double eps, const MEDCouplingUMesh *mesh2D, mcIdType offset, mcIdType neighbors[2]) const
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
          mcIdType cellId(mesh2D->getCellContainingPoint(barys->begin(),eps));
          if(cellId==-1)
            throw INTERP_KERNEL::Exception("EdgeInfo::feedEdgeInfoAt : internal error !");
          neighbors[0]=offset+cellId; neighbors[1]=offset+cellId;
        }
    }
}

class VectorOfCellInfo
{
public:
  VectorOfCellInfo(const std::vector<mcIdType>& edges, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& edgesPtr);
  std::size_t size() const { return _pool.size(); }
  mcIdType getPositionOf(double eps, const MEDCouplingUMesh *mesh) const;
  void setMeshAt(mcIdType pos, const MCAuto<MEDCouplingUMesh>& mesh, mcIdType istart, mcIdType iend, const MCAuto<MEDCouplingUMesh>& mesh1DInCase, const std::vector< std::vector<mcIdType> >& edges, const std::vector< std::vector< MCAuto<INTERP_KERNEL::Edge> > >& edgePtrs);
  const std::vector<mcIdType>& getConnOf(mcIdType pos) const { return get(pos)._edges; }
  const std::vector< MCAuto<INTERP_KERNEL::Edge> >& getEdgePtrOf(mcIdType pos) const { return get(pos)._edges_ptr; }
  MCAuto<MEDCouplingUMesh> getZeMesh() const { return _ze_mesh; }
  void feedEdgeInfoAt(double eps, mcIdType pos, mcIdType offset, mcIdType neighbors[2]) const;
private:
  mcIdType getZePosOfEdgeGivenItsGlobalId(mcIdType pos) const;
  void updateEdgeInfo(mcIdType pos, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newLeft, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newRight);
  const CellInfo& get(mcIdType pos) const;
  CellInfo& get(mcIdType pos);
private:
  std::vector<CellInfo> _pool;         // for a newly created 2D cell, the list of edges ToIdType( and edges ptr constiuing it
  MCAuto<MEDCouplingUMesh> _ze_mesh;   // the aggregated mesh
  std::vector<EdgeInfo> _edge_info;    // for each new edge added when cuting the 2D cell, the information on left and right bordering 2D cell
};

VectorOfCellInfo::VectorOfCellInfo(const std::vector<mcIdType>& edges, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& edgesPtr):_pool(1)
{
  _pool[0]._edges=edges;
  _pool[0]._edges_ptr=edgesPtr;
}

mcIdType VectorOfCellInfo::getPositionOf(double eps, const MEDCouplingUMesh *mesh) const
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

void VectorOfCellInfo::setMeshAt(mcIdType pos, const MCAuto<MEDCouplingUMesh>& mesh, mcIdType istart, mcIdType iend,
                                 const MCAuto<MEDCouplingUMesh>& mesh1DInCase, const std::vector< std::vector<mcIdType> >& edges,
                                 const std::vector< std::vector< MCAuto<INTERP_KERNEL::Edge> > >& edgePtrs)
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
  for(mcIdType i=0;i<pos;i++)
    pool[i]=_pool[i];
  for(std::size_t j=0;j<sz;j++)
    pool[pos+j]=CellInfo(edges[j],edgePtrs[j]);
  for(std::size_t i=pos+1;i<_pool.size();i++)
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

void VectorOfCellInfo::feedEdgeInfoAt(double eps, mcIdType pos, mcIdType offset, mcIdType neighbors[2]) const
{
  _edge_info[getZePosOfEdgeGivenItsGlobalId(pos)].feedEdgeInfoAt(eps,_ze_mesh,offset,neighbors);
}

mcIdType VectorOfCellInfo::getZePosOfEdgeGivenItsGlobalId(mcIdType pos) const
{
  if(pos<0)
    throw INTERP_KERNEL::Exception("VectorOfCellInfo::getZePosOfEdgeGivenItsGlobalId : invalid id ! Must be >=0 !");
  mcIdType ret(0);
  for(std::vector<EdgeInfo>::const_iterator it=_edge_info.begin();it!=_edge_info.end();it++,ret++)
    {
      if((*it).isInMyRange(pos))
        return ret;
    }
  throw INTERP_KERNEL::Exception("VectorOfCellInfo::getZePosOfEdgeGivenItsGlobalId : invalid id !");
}

void VectorOfCellInfo::updateEdgeInfo(mcIdType pos, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newLeft, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& newRight)
{
  get(pos);//to perform the sanity check;
  if(_edge_info.empty())
    return ;
  std::size_t sz(_edge_info.size()-1);
  for(std::size_t i=0;i<sz;i++)
    _edge_info[i].somethingHappendAt(pos,newLeft,newRight);
}

const CellInfo& VectorOfCellInfo::get(mcIdType pos) const
{
  if(pos<0 || pos>=ToIdType(_pool.size()))
    throw INTERP_KERNEL::Exception("VectorOfCellSplitter::get const : invalid pos !");
  return _pool[pos];
}

CellInfo& VectorOfCellInfo::get(mcIdType pos)
{
  if(pos<0 || pos>=ToIdType(_pool.size()))
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
 * \param [in] allEdges a list of pairs (beginNode, endNode). Represents all edges (already cut) in the single 2D cell being handled here. Linked with \a allEdgesPtr to get the equation of edge.
 */
MEDCouplingUMesh *BuildMesh2DCutInternal(double eps, MEDCouplingUMesh *splitMesh1D, const std::vector<mcIdType>& allEdges, const std::vector< MCAuto<INTERP_KERNEL::Edge> >& allEdgesPtr, mcIdType offset,
                                         MCAuto<DataArrayIdType>& idsLeftRight)
{
  mcIdType nbCellsInSplitMesh1D=splitMesh1D->getNumberOfCells();
  if(nbCellsInSplitMesh1D==0)
    throw INTERP_KERNEL::Exception("BuildMesh2DCutInternal : internal error ! input 1D mesh must have at least one cell !");
  const mcIdType *cSplitPtr(splitMesh1D->getNodalConnectivity()->begin()),*ciSplitPtr(splitMesh1D->getNodalConnectivityIndex()->begin());
  std::size_t nb(allEdges.size()),jj;
  if(nb%2!=0)
    throw INTERP_KERNEL::Exception("BuildMesh2DCutFrom : internal error 2 !");
  std::vector<mcIdType> edge1Bis(nb*2);
  std::vector< MCAuto<INTERP_KERNEL::Edge> > edge1BisPtr(nb*2);
  std::copy(allEdges.begin(),allEdges.end(),edge1Bis.begin());
  std::copy(allEdges.begin(),allEdges.end(),edge1Bis.begin()+nb);
  std::copy(allEdgesPtr.begin(),allEdgesPtr.end(),edge1BisPtr.begin());
  std::copy(allEdgesPtr.begin(),allEdgesPtr.end(),edge1BisPtr.begin()+nb);
  //
  idsLeftRight=DataArrayIdType::New(); idsLeftRight->alloc(nbCellsInSplitMesh1D*2); idsLeftRight->fillWithValue(-2); idsLeftRight->rearrange(2);
  mcIdType *idsLeftRightPtr(idsLeftRight->getPointer());
  VectorOfCellInfo pool(edge1Bis,edge1BisPtr);

  // Compute contiguous parts of splitMesh1D. We can not make the full assumption that segments are consecutive in the connectivity
  // (even if the user correctly called orderConsecutiveCells1D()). Indeed the tool might be a closed line whose junction point is in
  // splitMesh1D. There can be only one such a point, and if this happens this is necessarily at the start
  // of the connectivity.
  MCAuto <DataArrayIdType> renumb(DataArrayIdType::New());
  renumb->alloc(nbCellsInSplitMesh1D,1);
  const mcIdType * renumbP(renumb->begin());

  mcIdType i, first=cSplitPtr[1];
  // Follow 1D line backward as long as it is connected:
  for (i=nbCellsInSplitMesh1D-1; cSplitPtr[ciSplitPtr[i]+2] == first; i--)
    first=cSplitPtr[ciSplitPtr[i]+1];
  if (i < nbCellsInSplitMesh1D-1)
    {
      // Build circular permutation to shift consecutive edges together
      renumb->iota(i+1);
      renumb->applyModulus(nbCellsInSplitMesh1D);
      splitMesh1D->renumberCells(renumbP, false);
      cSplitPtr = splitMesh1D->getNodalConnectivity()->begin();
      ciSplitPtr = splitMesh1D->getNodalConnectivityIndex()->begin();
    }
  else
    renumb->iota();
  //
  // The 1D first piece is used to intersect the 2D cell resulting in max two 2D cells.
  // The next 1D piece is localised (getPositionOf()) into this previous cut.
  // The result of the next intersection replaces the former single 2D cell that has been cut in the
  // pool. The neighbourhood information detained by pool._edge_info is also updated so that left and right
  // adjacent 2D cell of a 1D piece is kept up to date.
  // And so on and so forth.
  for(mcIdType iStart=0;iStart<nbCellsInSplitMesh1D;)
    {// split [0:nbCellsInSplitMesh1D) in contiguous parts [iStart:iEnd)
      mcIdType iEnd(iStart);
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

      MCAuto<MEDCouplingUMesh> partOfSplitMesh1D(static_cast<MEDCouplingUMesh *>(splitMesh1D->buildPartOfMySelfSlice(iStart,iEnd,1,true)));
      mcIdType pos(pool.getPositionOf(eps,partOfSplitMesh1D));
      //
      MCAuto<MEDCouplingUMesh>retTmp(MEDCouplingUMesh::New("",2));
      retTmp->setCoords(splitMesh1D->getCoords());
      retTmp->allocateCells();

      std::vector< std::vector<mcIdType> > out0;
      std::vector< std::vector< MCAuto<INTERP_KERNEL::Edge> > > out1;

      BuildMesh2DCutInternal2(partOfSplitMesh1D,pool.getConnOf(pos),pool.getEdgePtrOf(pos),out0,out1);
      for(std::size_t cnt=0;cnt<out0.size();cnt++)
        AddCellInMesh2D(retTmp,out0[cnt],out1[cnt]);
      pool.setMeshAt(pos,retTmp,iStart,iEnd,partOfSplitMesh1D,out0,out1);
      //
      iStart=iEnd;
    }
  for(mcIdType mm=0;mm<nbCellsInSplitMesh1D;mm++)
    pool.feedEdgeInfoAt(eps,renumbP[mm],offset,idsLeftRightPtr+2*mm);

  return pool.getZeMesh().retn();
}

/*
 * splitMesh1D is an input parameter but might have its cells renumbered.
 */
MEDCouplingUMesh *BuildMesh2DCutFrom(double eps, mcIdType cellIdInMesh2D, const MEDCouplingUMesh *mesh2DDesc, MEDCouplingUMesh *splitMesh1D,
                                     const mcIdType *descBg, const mcIdType *descEnd, const std::vector< std::vector<mcIdType> >& intersectEdge1, mcIdType offset,
                                     MCAuto<DataArrayIdType>& idsLeftRight)
{
  const mcIdType *cdescPtr(mesh2DDesc->getNodalConnectivity()->begin()),*cidescPtr(mesh2DDesc->getNodalConnectivityIndex()->begin());
  //
  std::vector<mcIdType> allEdges;
  std::vector< MCAuto<INTERP_KERNEL::Edge> > allEdgesPtr; // for each sub edge in splitMesh2D the uncut Edge object of the original mesh2D
  for(const mcIdType *it(descBg);it!=descEnd;it++) // for all edges in the descending connectivity of the 2D mesh in relative Fortran mode
    {
      mcIdType edgeId(std::abs(*it)-1);
      std::map< MCAuto<INTERP_KERNEL::Node>,mcIdType> m;
      MCAuto<INTERP_KERNEL::Edge> ee(MEDCouplingUMeshBuildQPFromEdge2((INTERP_KERNEL::NormalizedCellType)cdescPtr[cidescPtr[edgeId]],cdescPtr+cidescPtr[edgeId]+1,mesh2DDesc->getCoords()->begin(),m));
      const std::vector<mcIdType>& edge1(intersectEdge1[edgeId]);
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

bool AreEdgeEqual(const double *coo2D, const INTERP_KERNEL::CellModel& typ1, const mcIdType *conn1, const INTERP_KERNEL::CellModel& typ2, const mcIdType *conn2, double eps)
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
mcIdType FindRightCandidateAmong(const MEDCouplingUMesh *mesh2DSplit, const mcIdType *candidatesIn2DBg, const mcIdType *candidatesIn2DEnd, const MEDCouplingUMesh *mesh1DSplit, mcIdType cellIdInMesh1DSplitRelative, double eps)
{
  if(candidatesIn2DEnd==candidatesIn2DBg)
    throw INTERP_KERNEL::Exception("FindRightCandidateAmong : internal error 1 !");
  const double *coo(mesh2DSplit->getCoords()->begin());
  if(std::distance(candidatesIn2DBg,candidatesIn2DEnd)==1)
    return *candidatesIn2DBg;
  mcIdType edgeId(std::abs(cellIdInMesh1DSplitRelative)-1);
  MCAuto<MEDCouplingUMesh> cur1D(static_cast<MEDCouplingUMesh *>(mesh1DSplit->buildPartOfMySelf(&edgeId,&edgeId+1,true)));
  if(cellIdInMesh1DSplitRelative<0)
    cur1D->changeOrientationOfCells();
  const mcIdType *c1D(cur1D->getNodalConnectivity()->begin());
  const INTERP_KERNEL::CellModel& ref1DType(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)c1D[0]));
  for(const mcIdType *it=candidatesIn2DBg;it!=candidatesIn2DEnd;it++)
    {
      MCAuto<MEDCouplingUMesh> cur2D(static_cast<MEDCouplingUMesh *>(mesh2DSplit->buildPartOfMySelf(it,it+1,true)));
      const mcIdType *c(cur2D->getNodalConnectivity()->begin()),*ci(cur2D->getNodalConnectivityIndex()->begin());
      const INTERP_KERNEL::CellModel &cm(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)c[ci[0]]));
      unsigned sz(cm.getNumberOfSons2(c+ci[0]+1,ci[1]-ci[0]-1));
      INTERP_KERNEL::AutoPtr<mcIdType> tmpPtr(new mcIdType[ci[1]-ci[0]]);
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
 * \param [out] intersectEdge1 - for each cell in \a m1Desc returns the result of the split. The result is given using pair of mcIdType given resp start and stop.
 *                               So for all edge \a i in \a m1Desc \a  intersectEdge1[i] is of length 2*n where n is the number of sub edges.
 *                               And for each j in [1,n) intersect[i][2*(j-1)+1]==intersect[i][2*j].
 * \param [out] subDiv2 - for each cell in \a m2Desc returns nodes that split it using convention \a m1Desc first, then \a m2Desc, then addCoo
 * \param [out] colinear2 - for each cell in \a m2Desc returns the edges in \a m1Desc that are colinear to it.
 * \param [out] addCoo - nodes to be append at the end
 * \param [out] mergedNodes - gives all pair of nodes of \a m2Desc that have same location than some nodes in \a m1Desc. key is id in \a m2Desc offsetted and value is id in \a m1Desc.
 */
void MEDCouplingUMesh::Intersect1DMeshes(const MEDCouplingUMesh *m1Desc, const MEDCouplingUMesh *m2Desc, double eps,
                                         std::vector< std::vector<mcIdType> >& intersectEdge1, std::vector< std::vector<mcIdType> >& colinear2, std::vector< std::vector<mcIdType> >& subDiv2, std::vector<double>& addCoo, std::map<mcIdType,mcIdType>& mergedNodes)
{
  static const int SPACEDIM=2;
  INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);
  const mcIdType *c1(m1Desc->getNodalConnectivity()->begin()),*ci1(m1Desc->getNodalConnectivityIndex()->begin());
  // Build BB tree of all edges in the tool mesh (second mesh)
  MCAuto<DataArrayDouble> bbox1Arr(m1Desc->getBoundingBoxForBBTree(eps)),bbox2Arr(m2Desc->getBoundingBoxForBBTree(eps));
  const double *bbox1(bbox1Arr->begin()),*bbox2(bbox2Arr->begin());
  mcIdType nDescCell1=m1Desc->getNumberOfCells(),nDescCell2=m2Desc->getNumberOfCells();
  intersectEdge1.resize(nDescCell1);
  colinear2.resize(nDescCell2);
  subDiv2.resize(nDescCell2);
  BBTree<SPACEDIM,mcIdType> myTree(bbox2,0,0,m2Desc->getNumberOfCells(),-eps);
  BBTreePts<SPACEDIM,mcIdType> treeNodes2(m2Desc->getCoords()->begin(),0,0,m2Desc->getCoords()->getNumberOfTuples(),eps);

  std::vector<mcIdType> candidates1(1);
  mcIdType offset1(m1Desc->getNumberOfNodes());
  mcIdType offset2(offset1+m2Desc->getNumberOfNodes());
  for(mcIdType i=0;i<nDescCell1;i++)  // for all edges in the first mesh
    {
      std::vector<mcIdType> candidates2; // edges of mesh2 candidate for intersection
      myTree.getIntersectingElems(bbox1+i*2*SPACEDIM,candidates2);
      if(!candidates2.empty()) // candidates2 holds edges from the second mesh potentially intersecting current edge i in mesh1
        {
          std::map<INTERP_KERNEL::Node *,mcIdType> map1,map2;
          std::map<mcIdType, INTERP_KERNEL::Node *> revMap2;
          // pol2 is not necessarily a closed polygon: just a set of (quadratic) edges (same as candidates2) in the Geometric DS format
          INTERP_KERNEL::QuadraticPolygon *pol2=MEDCouplingUMeshBuildQPFromMesh(m2Desc,candidates2,map2);
          // Build revMap2
          for (auto& kv : map2)
            revMap2[kv.second] = kv.first;
          candidates1[0]=i;
          // In the construction of pol1 we might reuse nodes from pol2, that we have identified as to be merged.
          INTERP_KERNEL::QuadraticPolygon *pol1=MEDCouplingUMeshBuildQPFromMeshWithTree(m1Desc,candidates1,map1,treeNodes2, revMap2);
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
          // Performs edge cutting:
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
                                                   std::vector< std::vector<mcIdType> >& intersectEdge1, std::vector< std::vector<mcIdType> >& colinear2, std::vector< std::vector<mcIdType> >& subDiv2,
                                                   MEDCouplingUMesh *& m1Desc, DataArrayIdType *&desc1, DataArrayIdType *&descIndx1, DataArrayIdType *&revDesc1, DataArrayIdType *&revDescIndx1,
                                                   std::vector<double>& addCoo,
                                                   MEDCouplingUMesh *& m2Desc, DataArrayIdType *&desc2, DataArrayIdType *&descIndx2, DataArrayIdType *&revDesc2, DataArrayIdType *&revDescIndx2)
{
  // Build desc connectivity
  desc1=DataArrayIdType::New(); descIndx1=DataArrayIdType::New(); revDesc1=DataArrayIdType::New(); revDescIndx1=DataArrayIdType::New();
  desc2=DataArrayIdType::New();
  descIndx2=DataArrayIdType::New();
  revDesc2=DataArrayIdType::New();
  revDescIndx2=DataArrayIdType::New();
  MCAuto<DataArrayIdType> dd1(desc1),dd2(descIndx1),dd3(revDesc1),dd4(revDescIndx1);
  MCAuto<DataArrayIdType> dd5(desc2),dd6(descIndx2),dd7(revDesc2),dd8(revDescIndx2);
  m1Desc=m1->buildDescendingConnectivity2(desc1,descIndx1,revDesc1,revDescIndx1);
  m2Desc=m2->buildDescendingConnectivity2(desc2,descIndx2,revDesc2,revDescIndx2);
  MCAuto<MEDCouplingUMesh> dd9(m1Desc),dd10(m2Desc);
  std::map<mcIdType,mcIdType> notUsedMap;
  Intersect1DMeshes(m1Desc,m2Desc,eps,intersectEdge1,colinear2,subDiv2,addCoo,notUsedMap);
  m1Desc->incrRef(); desc1->incrRef(); descIndx1->incrRef(); revDesc1->incrRef(); revDescIndx1->incrRef();
  m2Desc->incrRef(); desc2->incrRef(); descIndx2->incrRef(); revDesc2->incrRef(); revDescIndx2->incrRef();
}

/**
 * Private. Third step of the partitioning algorithm (Intersect2DMeshes): reconstruct full 2D cells from the
 * (newly created) nodes corresponding to the edge intersections.
 * Output params:
 * @param[out] cr connectivity of the resulting mesh
 * @param[out] crI connectivity of the resulting mesh
 * @param[out] cNb1 correspondence arrays giving for the merged mesh the initial cells IDs in m1 / m2
 * @param[out] cNb2 correspondence arrays giving for the merged mesh the initial cells IDs in m1 / m2
 * TODO: describe input parameters
 */
void MEDCouplingUMesh::BuildIntersecting2DCellsFromEdges(double eps, const MEDCouplingUMesh *m1, const mcIdType *desc1, const mcIdType *descIndx1,
                                                         const std::vector<std::vector<mcIdType> >& intesctEdges1, const std::vector< std::vector<mcIdType> >& colinear2,
                                                         const MEDCouplingUMesh *m2, const mcIdType *desc2, const mcIdType *descIndx2, const std::vector<std::vector<mcIdType> >& intesctEdges2,
                                                         const std::vector<double>& addCoords,
                                                         std::vector<double>& addCoordsQuadratic, std::vector<mcIdType>& cr, std::vector<mcIdType>& crI, std::vector<mcIdType>& cNb1, std::vector<mcIdType>& cNb2)
{
  static const int SPACEDIM=2;
  const double *coo1(m1->getCoords()->begin());
  const mcIdType *conn1(m1->getNodalConnectivity()->begin()),*connI1(m1->getNodalConnectivityIndex()->begin());
  mcIdType offset1(m1->getNumberOfNodes());
  const double *coo2(m2->getCoords()->begin());
  const mcIdType *conn2(m2->getNodalConnectivity()->begin()),*connI2(m2->getNodalConnectivityIndex()->begin());
  mcIdType offset2(offset1+m2->getNumberOfNodes());
  mcIdType offset3(offset2+ToIdType(addCoords.size())/2);
  MCAuto<DataArrayDouble> bbox1Arr(m1->getBoundingBoxForBBTree(eps)),bbox2Arr(m2->getBoundingBoxForBBTree(eps));
  const double *bbox1(bbox1Arr->begin()),*bbox2(bbox2Arr->begin());
  // Here a BBTree on 2D-cells, not on segments:
  BBTree<SPACEDIM,mcIdType> myTree(bbox2,0,0,m2->getNumberOfCells(),eps);
  mcIdType ncell1=m1->getNumberOfCells();
  crI.push_back(0);
  for(mcIdType i=0;i<ncell1;i++)
    {
      std::vector<mcIdType> candidates2;
      myTree.getIntersectingElems(bbox1+i*2*SPACEDIM,candidates2);
      std::map<INTERP_KERNEL::Node *,mcIdType> mapp;
      std::map<mcIdType,INTERP_KERNEL::Node *> mappRev;
      INTERP_KERNEL::QuadraticPolygon pol1;
      INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)conn1[connI1[i]];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
      // Populate mapp and mappRev with nodes from the current cell (i) from mesh1 - this also builds the Node* objects:
      MEDCouplingUMeshBuildQPFromMesh3(coo1,offset1,coo2,offset2,addCoords,desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1,/* output */mapp,mappRev);
      // pol1 is the full cell from mesh1, in QP format, with all the additional intersecting nodes.
      pol1.buildFromCrudeDataArray(mappRev,cm.isQuadratic(),conn1+connI1[i]+1,coo1,
          desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1);
      //
      std::set<INTERP_KERNEL::Edge *> edges1;// store all edges of pol1 that are NOT consumed by intersect cells. If any after iteration over candidates2 -> a part of pol1 should appear in result
      std::set<INTERP_KERNEL::Edge *> edgesBoundary2;// store all edges that are on boundary of (pol2 intersect pol1) minus edges on pol1.
      INTERP_KERNEL::IteratorOnComposedEdge it1(&pol1);
      for(it1.first();!it1.finished();it1.next())
        edges1.insert(it1.current()->getPtr());
      //
      std::map<mcIdType,std::vector<INTERP_KERNEL::ElementaryEdge *> > edgesIn2ForShare; // common edges
      std::vector<INTERP_KERNEL::QuadraticPolygon> pol2s(candidates2.size());
      mcIdType ii=0;
      // Build, for each intersecting cell candidate from mesh2, the corresponding QP.
      // Again all the additional intersecting nodes are there.
      for(std::vector<mcIdType>::const_iterator it2=candidates2.begin();it2!=candidates2.end();it2++,ii++)
        {
          INTERP_KERNEL::NormalizedCellType typ2=(INTERP_KERNEL::NormalizedCellType)conn2[connI2[*it2]];
          const INTERP_KERNEL::CellModel& cm2=INTERP_KERNEL::CellModel::GetCellModel(typ2);
          // Complete mapping with elements coming from the current cell it2 in mesh2:
          MEDCouplingUMeshBuildQPFromMesh3(coo1,offset1,coo2,offset2,addCoords,desc2+descIndx2[*it2],desc2+descIndx2[*it2+1],intesctEdges2,/* output */mapp,mappRev);
          // pol2 is the new QP in the final merged result.
          pol2s[ii].buildFromCrudeDataArray2(mappRev,cm2.isQuadratic(),conn2+connI2[*it2]+1,coo2,desc2+descIndx2[*it2],desc2+descIndx2[*it2+1],intesctEdges2,
              pol1,desc1+descIndx1[i],desc1+descIndx1[i+1],intesctEdges1,colinear2, /* output */ edgesIn2ForShare);
        }
      // The cleaning below must be done after the full construction of all pol2s to correctly deal with shared edges:
      for (auto &p: pol2s)
        p.cleanDegeneratedConsecutiveEdges();
      edgesIn2ForShare.clear();  // removing temptation to use it further since it might now contain invalid edges.
      ///
      ii=0;
      // Now rebuild intersected cells from all this:
      for(std::vector<mcIdType>::const_iterator it2=candidates2.begin();it2!=candidates2.end();it2++,ii++)
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
      for(std::map<mcIdType,INTERP_KERNEL::Node *>::const_iterator it=mappRev.begin();it!=mappRev.end();it++)
        (*it).second->decrRef();
    }
}

void InsertNodeInConnIfNecessary(mcIdType nodeIdToInsert, std::vector<mcIdType>& conn, const double *coords, double eps)
{
  std::vector<mcIdType>::iterator it(std::find(conn.begin(),conn.end(),nodeIdToInsert));
  if(it!=conn.end())
    return ;
  std::size_t sz(conn.size());
  std::size_t found(std::numeric_limits<std::size_t>::max());
  for(std::size_t i=0;i<sz;i++)
    {
      mcIdType pt0(conn[i]),pt1(conn[(i+1)%sz]);
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

void SplitIntoToPart(const std::vector<mcIdType>& conn, mcIdType pt0, mcIdType pt1, std::vector<mcIdType>& part0, std::vector<mcIdType>& part1)
{
  std::size_t sz(conn.size());
  std::vector<mcIdType> *curPart(&part0);
  for(std::size_t i=0;i<sz;i++)
    {
      mcIdType nextt(conn[(i+1)%sz]);
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
void MEDCouplingUMesh::buildSubCellsFromCut(const std::vector< std::pair<mcIdType,mcIdType> >& cut3DSurf,
                                            const mcIdType *desc, const mcIdType *descIndx, const double *coords, double eps,
                                            std::vector<std::vector<mcIdType> >& res) const
{
  checkFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSubCellsFromCut works on umeshes with meshdim equal to 3 and spaceDim equal to 3 too!");
  const mcIdType *nodal3D(_nodal_connec->begin()),*nodalIndx3D(_nodal_connec_index->begin());
  mcIdType nbOfCells=getNumberOfCells();
  if(nbOfCells!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildSubCellsFromCut works only with single cell presently !");
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      mcIdType offset(descIndx[i]),nbOfFaces(descIndx[i+1]-offset);
      for(int j=0;j<nbOfFaces;j++)
        {
          const std::pair<mcIdType,mcIdType>& p=cut3DSurf[desc[offset+j]];
          const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)nodal3D[nodalIndx3D[i]]));
          mcIdType sz=nodalIndx3D[i+1]-nodalIndx3D[i]-1;
          INTERP_KERNEL::AutoPtr<mcIdType> tmp(new mcIdType[sz]);
          INTERP_KERNEL::NormalizedCellType cmsId;
          unsigned nbOfNodesSon(cm.fillSonCellNodalConnectivity2(j,nodal3D+nodalIndx3D[i]+1,sz,tmp,cmsId));
          std::vector<mcIdType> elt((mcIdType *)tmp,(mcIdType *)tmp+nbOfNodesSon);
          if(p.first!=-1 && p.second!=-1)
            {
              if(p.first!=-2)
                {
                  InsertNodeInConnIfNecessary(p.first,elt,coords,eps);
                  InsertNodeInConnIfNecessary(p.second,elt,coords,eps);
                  std::vector<mcIdType> elt1,elt2;
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
 * It is the linear part of MEDCouplingUMesh::split2DCells. Here no additional nodes will be added in \b this. So coordinates pointer remain unchanged (is not even touch).
 *
 * \sa MEDCouplingUMesh::split2DCells
 */
void MEDCouplingUMesh::split2DCellsLinear(const DataArrayIdType *desc, const DataArrayIdType *descI, const DataArrayIdType *subNodesInSeg, const DataArrayIdType *subNodesInSegI)
{
  checkConnectivityFullyDefined();
  mcIdType ncells=getNumberOfCells();
  mcIdType lgthToReach(getNodalConnectivityArrayLen()+subNodesInSeg->getNumberOfTuples());
  MCAuto<DataArrayIdType> c(DataArrayIdType::New()); c->alloc((std::size_t)lgthToReach);
  const mcIdType *subPtr(subNodesInSeg->begin()),*subIPtr(subNodesInSegI->begin()),*descPtr(desc->begin()),*descIPtr(descI->begin()),*oldConn(getNodalConnectivity()->begin());
  mcIdType *cPtr(c->getPointer()),*ciPtr(getNodalConnectivityIndex()->getPointer());
  mcIdType prevPosOfCi(ciPtr[0]);
  for(mcIdType i=0;i<ncells;i++,ciPtr++,descIPtr++)
    {
      mcIdType offset(descIPtr[0]),sz(descIPtr[1]-descIPtr[0]),deltaSz(0);
      *cPtr++=ToIdType(INTERP_KERNEL::NORM_POLYGON); *cPtr++=oldConn[prevPosOfCi+1];
      for(mcIdType j=0;j<sz;j++)
        {
          mcIdType offset2(subIPtr[descPtr[offset+j]]),sz2(subIPtr[descPtr[offset+j]+1]-subIPtr[descPtr[offset+j]]);
          for(mcIdType k=0;k<sz2;k++)
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
 * It is the quadratic part of MEDCouplingUMesh::split2DCells. Here some additional nodes can be added at the end of coordinates array object.
 *
 * \return  mcIdType - the number of new nodes created.
 * \sa MEDCouplingUMesh::split2DCells
 */
mcIdType MEDCouplingUMesh::split2DCellsQuadratic(const DataArrayIdType *desc, const DataArrayIdType *descI, const DataArrayIdType *subNodesInSeg, const DataArrayIdType *subNodesInSegI, const DataArrayIdType *mid, const DataArrayIdType *midI)
{
  checkConsistencyLight();
  mcIdType ncells=getNumberOfCells();
  mcIdType lgthToReach(getNodalConnectivityArrayLen()+2*subNodesInSeg->getNumberOfTuples());
  mcIdType nodesCnt(getNumberOfNodes());
  MCAuto<DataArrayIdType> c(DataArrayIdType::New()); c->alloc((std::size_t)lgthToReach);
  MCAuto<DataArrayDouble> addCoo(DataArrayDouble::New()); addCoo->alloc(0,1);
  const mcIdType *subPtr(subNodesInSeg->begin()),*subIPtr(subNodesInSegI->begin()),*descPtr(desc->begin()),*descIPtr(descI->begin()),*oldConn(getNodalConnectivity()->begin());
  const mcIdType *midPtr(mid->begin()),*midIPtr(midI->begin());
  const double *oldCoordsPtr(getCoords()->begin());
  mcIdType *cPtr(c->getPointer()),*ciPtr(getNodalConnectivityIndex()->getPointer());
  mcIdType prevPosOfCi(ciPtr[0]);
  for(mcIdType i=0;i<ncells;i++,ciPtr++,descIPtr++)
    {
      mcIdType offset(descIPtr[0]),sz(descIPtr[1]-descIPtr[0]),deltaSz(sz);
      for(mcIdType j=0;j<sz;j++)
        { mcIdType sz2(subIPtr[descPtr[offset+j]+1]-subIPtr[descPtr[offset+j]]); deltaSz+=sz2; }
      *cPtr++=ToIdType(INTERP_KERNEL::NORM_QPOLYG); cPtr[0]=oldConn[prevPosOfCi+1];
      for(mcIdType j=0;j<sz;j++)//loop over subedges of oldConn
        {
          mcIdType offset2(subIPtr[descPtr[offset+j]]),sz2(subIPtr[descPtr[offset+j]+1]-subIPtr[descPtr[offset+j]]),offset3(midIPtr[descPtr[offset+j]]);
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
          for(mcIdType k=0;k<sz2;k++)//loop over subsplit of current subedge
            {
              cPtr[1]=subPtr[offset2+k];
              cPtr[deltaSz]=InternalAddPoint(e,midPtr[offset3+k],oldCoordsPtr,cPtr[0],cPtr[1],*addCoo,nodesCnt); cPtr++;
            }
          mcIdType tmpEnd(oldConn[prevPosOfCi+1+(j+1)%sz]);
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
 * \b WARNING: the two meshes should be correctly oriented for this method to work properly. Methods changeSpaceDimension() and
 * orientCorrectly2DCells() can be used for this.
 * \b WARNING: the two meshes should be "clean" (no un-merged nodes, no non-conformal cells)
 *  \param [in] m1 - the first input mesh which is a partitioned object. The mesh must be so that each point in the space covered by \a m1
 *                      must be covered exactly by one entity, \b no \b more. If it is not the case, some tools are available to heal the mesh (conformize2D, mergeNodes)
 *  \param [in] m2 - the second input mesh which is a partition tool. The mesh must be so that each point in the space covered by \a m2
 *                      must be covered exactly by one entity, \b no \b more. If it is not the case, some tools are available to heal the mesh (conformize2D, mergeNodes)
 *  \param [in] eps - precision used to detect coincident mesh entities.
 *  \param [out] cellNb1 - a new instance of DataArrayIdType holding for each result
 *         cell an id of the cell of \a m1 it comes from. The caller is to delete
 *         this array using decrRef() as it is no more needed.
 *  \param [out] cellNb2 - a new instance of DataArrayIdType holding for each result
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
                                                      double eps, DataArrayIdType *&cellNb1, DataArrayIdType *&cellNb2)
{
  if(!m1 || !m2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Intersect2DMeshes : input meshes must be not NULL !");
  m1->checkFullyDefined();
  m2->checkFullyDefined();
  INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);
  if(m1->getMeshDimension()!=2 || m1->getSpaceDimension()!=2 || m2->getMeshDimension()!=2 || m2->getSpaceDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Intersect2DMeshes works on umeshes m1 AND m2  with meshdim equal to 2 and spaceDim equal to 2 too!");

  // Step 1: compute all edge intersections (new nodes)
  std::vector< std::vector<mcIdType> > intersectEdge1, colinear2, subDiv2;
  MEDCouplingUMesh *m1Desc=0,*m2Desc=0; // descending connec. meshes
  DataArrayIdType *desc1=0,*descIndx1=0,*revDesc1=0,*revDescIndx1=0,*desc2=0,*descIndx2=0,*revDesc2=0,*revDescIndx2=0;
  std::vector<double> addCoo,addCoordsQuadratic;  // coordinates of newly created nodes
  IntersectDescending2DMeshes(m1,m2,eps,intersectEdge1,colinear2, subDiv2,
                              m1Desc,desc1,descIndx1,revDesc1,revDescIndx1,
                              addCoo, m2Desc,desc2,descIndx2,revDesc2,revDescIndx2);
  revDesc1->decrRef(); revDescIndx1->decrRef(); revDesc2->decrRef(); revDescIndx2->decrRef();
  MCAuto<DataArrayIdType> dd1(desc1),dd2(descIndx1),dd3(desc2),dd4(descIndx2);
  MCAuto<MEDCouplingUMesh> dd5(m1Desc),dd6(m2Desc);

  // Step 2: re-order newly created nodes according to the ordering found in m2
  std::vector< std::vector<mcIdType> > intersectEdge2;
  BuildIntersectEdges(m1Desc,m2Desc,addCoo,subDiv2,intersectEdge2);
  subDiv2.clear(); dd5=0; dd6=0;

  // Step 3:
  std::vector<mcIdType> cr,crI; //no DataArrayIdType because interface with Geometric2D
  std::vector<mcIdType> cNb1,cNb2; //no DataArrayIdType because interface with Geometric2D
  BuildIntersecting2DCellsFromEdges(eps,m1,desc1->begin(),descIndx1->begin(),intersectEdge1,colinear2,m2,desc2->begin(),descIndx2->begin(),intersectEdge2,addCoo,
                                    /* outputs -> */addCoordsQuadratic,cr,crI,cNb1,cNb2);

  // Step 4: Prepare final result:
  MCAuto<DataArrayDouble> addCooDa(DataArrayDouble::New());
  addCooDa->alloc(addCoo.size()/2,2);
  std::copy(addCoo.begin(),addCoo.end(),addCooDa->getPointer());
  MCAuto<DataArrayDouble> addCoordsQuadraticDa(DataArrayDouble::New());
  addCoordsQuadraticDa->alloc(addCoordsQuadratic.size()/2,2);
  std::copy(addCoordsQuadratic.begin(),addCoordsQuadratic.end(),addCoordsQuadraticDa->getPointer());
  std::vector<const DataArrayDouble *> coordss(4);
  coordss[0]=m1->getCoords(); coordss[1]=m2->getCoords(); coordss[2]=addCooDa; coordss[3]=addCoordsQuadraticDa;
  MCAuto<DataArrayDouble> coo(DataArrayDouble::Aggregate(coordss));
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::New("Intersect2D",2));
  MCAuto<DataArrayIdType> conn(DataArrayIdType::New()); conn->alloc(cr.size(),1); std::copy(cr.begin(),cr.end(),conn->getPointer());
  MCAuto<DataArrayIdType> connI(DataArrayIdType::New()); connI->alloc(crI.size(),1); std::copy(crI.begin(),crI.end(),connI->getPointer());
  MCAuto<DataArrayIdType> c1(DataArrayIdType::New()); c1->alloc(cNb1.size(),1); std::copy(cNb1.begin(),cNb1.end(),c1->getPointer());
  MCAuto<DataArrayIdType> c2(DataArrayIdType::New()); c2->alloc(cNb2.size(),1); std::copy(cNb2.begin(),cNb2.end(),c2->getPointer());
  ret->setConnectivity(conn,connI,true);
  ret->setCoords(coo);
  cellNb1=c1.retn(); cellNb2=c2.retn();
  return ret.retn();
}

/*!
 * Partitions the first given 2D mesh using the second given 1D mesh as a tool.
 * Thus the final result contains the aggregation of nodes of \a mesh2D, then nodes of \a mesh1D, then new nodes that are the result of the intersection
 * and finally, in case of quadratic polygon the centers of edges new nodes.
 * The meshes should be in 2D space. In addition, returns two arrays mapping cells of the resulting mesh to cells of the input.
 *
 * \b WARNING: the 2D mesh should be correctly oriented for this method to work properly. Methods changeSpaceDimension() and
 * orientCorrectly2DCells() can be used for this.
 * \b WARNING: the two meshes should be "clean" (no un-merged nodes, no non-conformal cells)
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
void MEDCouplingUMesh::Intersect2DMeshWith1DLine(const MEDCouplingUMesh *mesh2D, const MEDCouplingUMesh *mesh1D, double eps, MEDCouplingUMesh *&splitMesh2D, MEDCouplingUMesh *&splitMesh1D, DataArrayIdType *&cellIdInMesh2D, DataArrayIdType *&cellIdInMesh1D)
{
  if(!mesh2D || !mesh1D)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Intersect2DMeshWith1DLine : input meshes must be not NULL !");
  mesh2D->checkFullyDefined();
  mesh1D->checkFullyDefined();
  const std::vector<std::string>& compNames(mesh2D->getCoords()->getInfoOnComponents());
  if(mesh2D->getMeshDimension()!=2 || mesh2D->getSpaceDimension()!=2 || mesh1D->getMeshDimension()!=1 || mesh1D->getSpaceDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::Intersect2DMeshWith1DLine works with mesh2D with spacedim=meshdim=2 and mesh1D with meshdim=1 spaceDim=2 !");
  // Step 1: compute all edge intersections (new nodes)
  std::vector< std::vector<mcIdType> > intersectEdge1, colinear2, subDiv2;
  std::vector<double> addCoo,addCoordsQuadratic;  // coordinates of newly created nodes
  INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);
  //
  // Build desc connectivity
  DataArrayIdType *desc1(DataArrayIdType::New()),*descIndx1(DataArrayIdType::New()),*revDesc1(DataArrayIdType::New()),*revDescIndx1(DataArrayIdType::New());
  MCAuto<DataArrayIdType> dd1(desc1),dd2(descIndx1),dd3(revDesc1),dd4(revDescIndx1);
  MCAuto<MEDCouplingUMesh> m1Desc(mesh2D->buildDescendingConnectivity2(desc1,descIndx1,revDesc1,revDescIndx1));
  std::map<mcIdType,mcIdType> mergedNodes;
  Intersect1DMeshes(m1Desc,mesh1D,eps,intersectEdge1,colinear2,subDiv2,addCoo,mergedNodes);
  // use mergeNodes to fix intersectEdge1
  for(std::vector< std::vector<mcIdType> >::iterator it0=intersectEdge1.begin();it0!=intersectEdge1.end();it0++)
    {
      std::size_t n((*it0).size()/2);
      mcIdType eltStart((*it0)[0]),eltEnd((*it0)[2*n-1]);
      std::map<mcIdType,mcIdType>::const_iterator it1;
      it1=mergedNodes.find(eltStart);
      if(it1!=mergedNodes.end())
        (*it0)[0]=(*it1).second;
      it1=mergedNodes.find(eltEnd);
      if(it1!=mergedNodes.end())
        (*it0)[2*n-1]=(*it1).second;
    }
  //
  MCAuto<DataArrayDouble> addCooDa(DataArrayDouble::New());
  addCooDa->useArray(&addCoo[0],false,DeallocType::C_DEALLOC,addCoo.size()/2,2);
  // Step 2: re-order newly created nodes according to the ordering found in m2
  std::vector< std::vector<mcIdType> > intersectEdge2;
  BuildIntersectEdges(m1Desc,mesh1D,addCoo,subDiv2,intersectEdge2);
  subDiv2.clear();
  // Step 3: compute splitMesh1D
  MCAuto<DataArrayIdType> idsInRet1Colinear,idsInDescMesh2DForIdsInRetColinear;
  MCAuto<DataArrayIdType> ret2(DataArrayIdType::New()); ret2->alloc(0,1);
  MCAuto<MEDCouplingUMesh> ret1(BuildMesh1DCutFrom(mesh1D,intersectEdge2,mesh2D->getCoords(),addCoo,mergedNodes,colinear2,intersectEdge1,
                                                   idsInRet1Colinear,idsInDescMesh2DForIdsInRetColinear));

  // ### Colinearity fix :
  // if a node in ret1 has been merged with an already existing node in mesh2D,
  // we might end up with edges in ret1 that are colinear with some edges in mesh2D
  // even if none of the edges of the two original meshes were colinear.
  // this procedure detects such edges and adds them in idsInRet1Colinear/idsInDescMesh2DForIdsInRetColinear
  // a- find edges in ret1 that are in contact with 2 cells
  MCAuto<DataArrayDouble> centerOfMassRet1(ret1->computeCellCenterOfMass());
  MCAuto<DataArrayIdType> cells,cellsIndex;
  mesh2D->getCellsContainingPoints(centerOfMassRet1->begin(),centerOfMassRet1->getNumberOfTuples(),eps,cells,cellsIndex);
  MCAuto<DataArrayIdType> cellsIndex2(DataArrayIdType::New()); cellsIndex2->alloc(0,1);
  if (cellsIndex->getNumberOfTuples() > 1)
    cellsIndex2 = cellsIndex->deltaShiftIndex();
  MCAuto<DataArrayIdType> idsInRet1With2Contacts(cellsIndex2->findIdsEqual(2));

  MCAuto<DataArrayIdType> realIdsInDesc2D(desc1->deepCopy());
  realIdsInDesc2D->abs(); realIdsInDesc2D->applyLin(1,-1);
  const mcIdType *cRet1(ret1->getNodalConnectivity()->begin()),*ciRet1(ret1->getNodalConnectivityIndex()->begin());
  for(const mcIdType *it=idsInRet1With2Contacts->begin();it!=idsInRet1With2Contacts->end();it++)
    {
      // b- find the edge that the 2 cells in m1Desc have in common:
      // this is the edge which is colinear with the one in ret1
      const mcIdType* cellId1 = cells->begin() + cellsIndex->begin()[*it];
      const mcIdType* cellId2 = cells->begin() + cellsIndex->begin()[*it]+1;

      std::set<mcIdType> s1(realIdsInDesc2D->begin()+dd2->begin()[*cellId1], realIdsInDesc2D->begin()+dd2->begin()[*cellId1+1]);
      std::set<mcIdType> s2(realIdsInDesc2D->begin()+dd2->begin()[*cellId2], realIdsInDesc2D->begin()+dd2->begin()[*cellId2+1]);

      std::vector<mcIdType> commonEdgeId;
      std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(), std::back_inserter(commonEdgeId));

      // c- find correct orientation for commonEdgeId
      const mcIdType* firstNodeInColinearEdgeRet1 = cRet1 + ciRet1[*it]+1;
      const mcIdType* secondNodeInColinearEdgeRet1 =  cRet1 + ciRet1[*it]+2;
      mcIdType commonEdgeIdCorrectlyOriented;
      if(IsColinearOfACellOf(intersectEdge1, commonEdgeId, *firstNodeInColinearEdgeRet1,*secondNodeInColinearEdgeRet1, commonEdgeIdCorrectlyOriented))
        {
          idsInRet1Colinear->pushBackSilent(*it);
          idsInDescMesh2DForIdsInRetColinear->pushBackSilent(commonEdgeIdCorrectlyOriented);
        }
    }
  // ### End colinearity fix

  MCAuto<DataArrayIdType> ret3(DataArrayIdType::New()); ret3->alloc(ret1->getNumberOfCells()*2,1); ret3->fillWithValue(std::numeric_limits<mcIdType>::max()); ret3->rearrange(2);
  MCAuto<DataArrayIdType> idsInRet1NotColinear(idsInRet1Colinear->buildComplement(ret1->getNumberOfCells()));
  // deal with cells in mesh2D that are not cut but only some of their edges are
  MCAuto<DataArrayIdType> idsInDesc2DToBeRefined(idsInDescMesh2DForIdsInRetColinear->deepCopy());
  idsInDesc2DToBeRefined->abs(); idsInDesc2DToBeRefined->applyLin(1,-1);
  idsInDesc2DToBeRefined=idsInDesc2DToBeRefined->buildUnique();

  MCAuto<DataArrayIdType> out0s;//ids in mesh2D that are impacted by the fact that some edges of \a mesh1D are part of the edges of those cells
  if(!idsInDesc2DToBeRefined->empty())
    {
      DataArrayIdType *out0(0),*outi0(0);
      DataArrayIdType::ExtractFromIndexedArrays(idsInDesc2DToBeRefined->begin(),idsInDesc2DToBeRefined->end(),dd3,dd4,out0,outi0);
      MCAuto<DataArrayIdType> outi0s(outi0);
      out0s=out0;
      out0s=out0s->buildUnique();
      out0s->sort(true);
    }
  //
  MCAuto<MEDCouplingUMesh> ret1NonCol(static_cast<MEDCouplingUMesh *>(ret1->buildPartOfMySelf(idsInRet1NotColinear->begin(),idsInRet1NotColinear->end())));
  MCAuto<DataArrayDouble> baryRet1(centerOfMassRet1->selectByTupleId(idsInRet1NotColinear->begin(), idsInRet1NotColinear->end()));
  DataArrayIdType *out(0),*outi(0);
  DataArrayIdType::ExtractFromIndexedArrays(idsInRet1NotColinear->begin(),idsInRet1NotColinear->end(),cells,cellsIndex,out,outi);
  MCAuto<DataArrayIdType> elts(out),eltsIndex(outi);

  MCAuto<DataArrayIdType> eltsIndex2(DataArrayIdType::New()); eltsIndex2->alloc(0,1);
  if (eltsIndex->getNumberOfTuples() > 1)
    eltsIndex2 = eltsIndex->deltaShiftIndex();
  MCAuto<DataArrayIdType> eltsIndex3(eltsIndex2->findIdsEqual(1));
  if(eltsIndex2->count(0)+eltsIndex3->getNumberOfTuples()!=ret1NonCol->getNumberOfCells())
    throw INTERP_KERNEL::Exception("Intersect2DMeshWith1DLine : internal error 1 !");
  MCAuto<DataArrayIdType> cellsToBeModified(elts->buildUnique());
  MCAuto<DataArrayIdType> untouchedCells(cellsToBeModified->buildComplement(mesh2D->getNumberOfCells()));
  if((DataArrayIdType *)out0s)
    untouchedCells=untouchedCells->buildSubstraction(out0s);//if some edges in ret1 are colinear to descending mesh of mesh2D remove cells from untouched one
  std::vector< MCAuto<MEDCouplingUMesh> > outMesh2DSplit;
  // OK all is ready to insert in ret2 mesh
  if(!untouchedCells->empty())
    {// the most easy part, cells in mesh2D not impacted at all
      outMesh2DSplit.push_back(static_cast<MEDCouplingUMesh *>(mesh2D->buildPartOfMySelf(untouchedCells->begin(),untouchedCells->end())));
      outMesh2DSplit.back()->setCoords(ret1->getCoords());
      ret2->pushBackValsSilent(untouchedCells->begin(),untouchedCells->end());
    }
  if((DataArrayIdType *)out0s)
    {// here dealing with cells in out0s but not in cellsToBeModified
      MCAuto<DataArrayIdType> fewModifiedCells(out0s->buildSubstraction(cellsToBeModified));
      const mcIdType *rdptr(dd3->begin()),*rdiptr(dd4->begin()),*dptr(dd1->begin()),*diptr(dd2->begin());
      for(const mcIdType *it=fewModifiedCells->begin();it!=fewModifiedCells->end();it++)
        {
          outMesh2DSplit.push_back(BuildRefined2DCell(ret1->getCoords(),mesh2D,*it,dptr+diptr[*it],dptr+diptr[*it+1],intersectEdge1));
          ret1->setCoords(outMesh2DSplit.back()->getCoords());
        }
      mcIdType offset(ret2->getNumberOfTuples());
      ret2->pushBackValsSilent(fewModifiedCells->begin(),fewModifiedCells->end());
      MCAuto<DataArrayIdType> partOfRet3(DataArrayIdType::New()); partOfRet3->alloc(2*idsInRet1Colinear->getNumberOfTuples(),1);
      partOfRet3->fillWithValue(std::numeric_limits<mcIdType>::max()); partOfRet3->rearrange(2);
      mcIdType kk(0),*ret3ptr(partOfRet3->getPointer());
      for(const mcIdType *it=idsInDescMesh2DForIdsInRetColinear->begin();it!=idsInDescMesh2DForIdsInRetColinear->end();it++,kk++)
        {
          mcIdType faceId(std::abs(*it)-1);
          for(const mcIdType *it2=rdptr+rdiptr[faceId];it2!=rdptr+rdiptr[faceId+1];it2++)
            {
              mcIdType tmp(fewModifiedCells->findIdFirstEqual(*it2));
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
  for(const mcIdType *it=cellsToBeModified->begin();it!=cellsToBeModified->end();it++)
    {
      MCAuto<DataArrayIdType> idsNonColPerCell(elts->findIdsEqual(*it));
      idsNonColPerCell->transformWithIndArr(eltsIndex3->begin(),eltsIndex3->end());
      MCAuto<DataArrayIdType> idsNonColPerCell2(idsInRet1NotColinear->selectByTupleId(idsNonColPerCell->begin(),idsNonColPerCell->end()));
      MCAuto<MEDCouplingUMesh> partOfMesh1CuttingCur2DCell(static_cast<MEDCouplingUMesh *>(ret1NonCol->buildPartOfMySelf(idsNonColPerCell->begin(),idsNonColPerCell->end())));
      MCAuto<DataArrayIdType> partOfRet3;
      MCAuto<MEDCouplingUMesh> splitOfOneCell(BuildMesh2DCutFrom(eps,*it,m1Desc,partOfMesh1CuttingCur2DCell,dd1->begin()+dd2->getIJ(*it,0),dd1->begin()+dd2->getIJ((*it)+1,0),intersectEdge1,ret2->getNumberOfTuples(),partOfRet3));
      ret3->setPartOfValues3(partOfRet3,idsNonColPerCell2->begin(),idsNonColPerCell2->end(),0,2,1,true);
      outMesh2DSplit.push_back(splitOfOneCell);
      for(mcIdType i=0;i<splitOfOneCell->getNumberOfCells();i++)
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
  // To finish - filter ret3 - std::numeric_limits<mcIdType>::max() -> -1 - negate values must be resolved.
  ret3->rearrange(1);
  MCAuto<DataArrayIdType> edgesToDealWith(ret3->findIdsStrictlyNegative());
  for(const mcIdType *it=edgesToDealWith->begin();it!=edgesToDealWith->end();it++)
    {
      mcIdType old2DCellId(-ret3->getIJ(*it,0)-1);
      MCAuto<DataArrayIdType> candidates(ret2->findIdsEqual(old2DCellId));
      ret3->setIJ(*it,0,FindRightCandidateAmong(ret2D,candidates->begin(),candidates->end(),ret1,*it%2==0?-((*it)/2+1):(*it)/2+1,eps));// div by 2 because 2 components natively in ret3
    }
  ret3->changeValue(std::numeric_limits<mcIdType>::max(),-1);
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
 * This method performs a conformization of \b this. So if a edge in \a this can be split into entire edges in \a this method
 * will suppress such edges to use sub edges in \a this. So this method does not add nodes in \a this if merged edges are both linear (INTERP_KERNEL::NORM_SEG2).
 * In the other cases new nodes can be created. If any are created, they will be appended at the end of the coordinates object before the invocation of this method.
 *
 * Whatever the returned value, this method does not alter the order of cells in \a this neither the orientation of cells.
 * The modified cells, if any, are systematically declared as NORM_POLYGON or NORM_QPOLYG depending on the initial quadraticness of geometric type.
 *
 * This method expects that \b this has a meshDim equal 2 and spaceDim equal to 2 too.
 * This method expects that all nodes in \a this are not closer than \a eps.
 * If it is not the case you can invoke MEDCouplingUMesh::mergeNodes before calling this method.
 *
 * \param [in] eps the relative error to detect merged edges.
 * \return DataArrayIdType  * - The list of cellIds in \a this that have been subdivided. If empty, nothing changed in \a this (as if it were a const method). The array is a newly allocated array
 *                           that the user is expected to deal with.
 *
 * \throw If \a this is not coherent.
 * \throw If \a this has not spaceDim equal to 2.
 * \throw If \a this has not meshDim equal to 2.
 * \sa MEDCouplingUMesh::mergeNodes, MEDCouplingUMesh::split2DCells
 */
DataArrayIdType *MEDCouplingUMesh::conformize2D(double eps)
{
  static const int SPACEDIM=2;
  checkConsistencyLight();
  if(getSpaceDimension()!=2 || getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::conformize2D : This method only works for meshes with spaceDim=2 and meshDim=2 !");
  MCAuto<DataArrayIdType> desc1(DataArrayIdType::New()),descIndx1(DataArrayIdType::New()),revDesc1(DataArrayIdType::New()),revDescIndx1(DataArrayIdType::New());
  MCAuto<MEDCouplingUMesh> mDesc(buildDescendingConnectivity(desc1,descIndx1,revDesc1,revDescIndx1));
  const mcIdType *c(mDesc->getNodalConnectivity()->begin()),*ci(mDesc->getNodalConnectivityIndex()->begin()),*rd(revDesc1->begin()),*rdi(revDescIndx1->begin());
  MCAuto<DataArrayDouble> bboxArr(mDesc->getBoundingBoxForBBTree(eps));
  const double *bbox(bboxArr->begin()),*coords(getCoords()->begin());
  mcIdType nCell=getNumberOfCells(),nDescCell=mDesc->getNumberOfCells();
  std::vector< std::vector<mcIdType> > intersectEdge(nDescCell),overlapEdge(nDescCell);
  std::vector<double> addCoo;
  BBTree<SPACEDIM,mcIdType> myTree(bbox,0,0,nDescCell,-eps);
  INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);
  for(mcIdType i=0;i<nDescCell;i++)
    {
      std::vector<mcIdType> candidates;
      myTree.getIntersectingElems(bbox+i*2*SPACEDIM,candidates);
      for(std::vector<mcIdType>::const_iterator it=candidates.begin();it!=candidates.end();it++)
        if(*it>i)  // we're dealing with pair of edges, no need to treat the same pair twice
          {
            std::map<MCAuto<INTERP_KERNEL::Node>,mcIdType> m;
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
  std::vector< std::vector<mcIdType> > middle(nDescCell);
  mcIdType nbOf2DCellsToBeSplit(0);
  bool middleNeedsToBeUsed(false);
  std::vector<bool> cells2DToTreat(nDescCell,false);
  for(mcIdType i=0;i<nDescCell;i++)
    {
      std::vector<mcIdType>& isect(intersectEdge[i]);
      std::size_t sz(isect.size());
      if(sz>1)
        {
          std::map<MCAuto<INTERP_KERNEL::Node>,mcIdType> m;
          INTERP_KERNEL::Edge *e(MEDCouplingUMeshBuildQPFromEdge2((INTERP_KERNEL::NormalizedCellType)c[ci[i]],c+ci[i]+1,coords,m));
          e->sortSubNodesAbs(coords,isect);
          e->decrRef();
        }
      if(sz!=0)
        {
          mcIdType idx0(rdi[i]),idx1(rdi[i+1]);
          if(idx1-idx0!=1)
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::conformize2D : internal error #0 !");
          if(!cells2DToTreat[rd[idx0]])
            {
              cells2DToTreat[rd[idx0]]=true;
              nbOf2DCellsToBeSplit++;
            }
          // try to reuse at most eventual 'middle' of SEG3
          std::vector<mcIdType>& mid(middle[i]);
          mid.resize(sz+1,-1);
          if((INTERP_KERNEL::NormalizedCellType)c[ci[i]]==INTERP_KERNEL::NORM_SEG3)
            {
              middleNeedsToBeUsed=true;
              const std::vector<mcIdType>& candidates(overlapEdge[i]);
              std::vector<mcIdType> trueCandidates;
              for(std::vector<mcIdType>::const_iterator itc=candidates.begin();itc!=candidates.end();itc++)
                if((INTERP_KERNEL::NormalizedCellType)c[ci[*itc]]==INTERP_KERNEL::NORM_SEG3)
                  trueCandidates.push_back(*itc);
              mcIdType stNode(c[ci[i]+1]),endNode(isect[0]);
              for(std::size_t j=0;j<sz+1;j++)
                {
                  for(std::vector<mcIdType>::const_iterator itc=trueCandidates.begin();itc!=trueCandidates.end();itc++)
                    {
                      mcIdType tmpSt(c[ci[*itc]+1]),tmpEnd(c[ci[*itc]+2]);
                      if((tmpSt==stNode && tmpEnd==endNode) || (tmpSt==endNode && tmpEnd==stNode))
                        { mid[j]=*itc; break; }
                    }
                  stNode=endNode;
                  endNode=j<sz-1?isect[j+1]:c[ci[i]+2];
                }
            }
        }
    }
  MCAuto<DataArrayIdType> ret(DataArrayIdType::New()),notRet(DataArrayIdType::New()); ret->alloc(nbOf2DCellsToBeSplit,1);
  if(nbOf2DCellsToBeSplit==0)
    return ret.retn();
  //
  mcIdType *retPtr(ret->getPointer());
  for(mcIdType i=0;i<nCell;i++)
    if(cells2DToTreat[i])
      *retPtr++=i;
  //
  MCAuto<DataArrayIdType> mSafe,nSafe,oSafe,pSafe,qSafe,rSafe;
  DataArrayIdType *m(0),*n(0),*o(0),*p(0),*q(0),*r(0);
  DataArrayIdType::ExtractFromIndexedArrays(ret->begin(),ret->end(),desc1,descIndx1,m,n); mSafe=m; nSafe=n;
  DataArrayIdType::PutIntoToSkylineFrmt(intersectEdge,o,p); oSafe=o; pSafe=p;
  if(middleNeedsToBeUsed)
    { DataArrayIdType::PutIntoToSkylineFrmt(middle,q,r); qSafe=q; rSafe=r; }
  MCAuto<MEDCouplingUMesh> modif(static_cast<MEDCouplingUMesh *>(buildPartOfMySelf(ret->begin(),ret->end(),true)));
  mcIdType nbOfNodesCreated(modif->split2DCells(mSafe,nSafe,oSafe,pSafe,qSafe,rSafe));
  setCoords(modif->getCoords());//if nbOfNodesCreated==0 modif and this have the same coordinates pointer so this line has no effect. But for quadratic cases this line is important.
  setPartOfMySelf(ret->begin(),ret->end(),*modif);
  {
    bool areNodesMerged; mcIdType newNbOfNodes;
    if(nbOfNodesCreated!=0)
      MCAuto<DataArrayIdType> tmp(mergeNodes(eps,areNodesMerged,newNbOfNodes));
  }
  return ret.retn();
}

/*!
 * This non const method works on 2D mesh. This method scans every cell in \a this and look if each edge constituting this cell is not mergeable with neighbors edges of that cell.
 * If yes, the cell is "repaired" to minimize at most its number of edges. So this method do not change the overall shape of cells in \a this (with eps precision).
 * This method do not take care of shared edges between cells, so this method can lead to a non conform mesh (\a this). If a conform mesh is required you're expected
 * to invoke MEDCouplingUMesh::mergeNodes and MEDCouplingUMesh::conformize2D right after this call.
 * This method works on any 2D geometric types of cell (even static one). If a cell is touched its type becomes dynamic automatically. For 2D "repaired" quadratic cells
 * new nodes for center of merged edges is are systematically created and appended at the end of the previously existing nodes.
 *
 * If the returned array is empty it means that nothing has changed in \a this (as if it were a const method). If the array is not empty the connectivity of \a this is modified
 * using new instance, idem for coordinates.
 *
 * If \a this is constituted by only linear 2D cells, this method is close to the computation of the convex hull of each cells in \a this.
 *
 * \return DataArrayIdType  * - The list of cellIds in \a this that have at least one edge colinearized.
 *
 * \throw If \a this is not coherent.
 * \throw If \a this has not spaceDim equal to 2.
 * \throw If \a this has not meshDim equal to 2.
 *
 * \sa MEDCouplingUMesh::conformize2D, MEDCouplingUMesh::mergeNodes, MEDCouplingUMesh::convexEnvelop2D.
 */
DataArrayIdType *MEDCouplingUMesh::colinearize2D(double eps)
{
  return internalColinearize2D(eps, false);
}

/*!
 * Performs exactly the same job as colinearize2D, except that this function does not create new non-conformal cells.
 * In a given 2D cell, if two edges are colinear and the junction point is used by a third edge, the two edges will not be
 * merged, contrary to colinearize2D().
 *
 * \sa MEDCouplingUMesh::colinearize2D
 */
DataArrayIdType *MEDCouplingUMesh::colinearizeKeepingConform2D(double eps)
{
  return internalColinearize2D(eps, true);
}


/*!
 * \param stayConform is set to True, will not fuse two edges sharing a node that has (strictly) more than 2 egdes connected to it
 */
DataArrayIdType *MEDCouplingUMesh::internalColinearize2D(double eps, bool stayConform)
{
  MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(0,1);
  checkConsistencyLight();
  if(getSpaceDimension()!=2 || getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::colinearize2D : This method only works for meshes with spaceDim=2 and meshDim=2 !");
  INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);
  mcIdType nbOfCells(getNumberOfCells()),nbOfNodes(getNumberOfNodes());
  const mcIdType *cptr(_nodal_connec->begin()),*ciptr(_nodal_connec_index->begin());
  MCAuto<DataArrayIdType> newc(DataArrayIdType::New()),newci(DataArrayIdType::New()); newci->alloc(nbOfCells+1,1); newc->alloc(0,1); newci->setIJ(0,0,0);
  std::map<mcIdType, bool> forbiddenPoints;  // list of points that can not be removed (or it will break conformity)
  if(stayConform) //
    {
      // A point that is used by more than 2 edges can not be removed without breaking conformity:
      MCAuto<DataArrayIdType> desc1(DataArrayIdType::New()),descI1(DataArrayIdType::New()),revDesc1(DataArrayIdType::New()),revDescI1(DataArrayIdType::New());
      MCAuto<MEDCouplingUMesh> mDesc1D(buildDescendingConnectivity(desc1,descI1,revDesc1,revDescI1));
      MCAuto<DataArrayIdType> desc2(DataArrayIdType::New()),descI2(DataArrayIdType::New()),revDesc2(DataArrayIdType::New()),revDescI2(DataArrayIdType::New());
      MCAuto<MEDCouplingUMesh> mDesc0D(mDesc1D->buildDescendingConnectivity(desc2,descI2,revDesc2,revDescI2));
      MCAuto<DataArrayIdType> dsi(revDescI2->deltaShiftIndex());
      MCAuto<DataArrayIdType> ids(dsi->findIdsGreaterThan(2));
      const mcIdType * cPtr(mDesc0D->getNodalConnectivity()->begin());
      for(auto it = ids->begin(); it != ids->end(); it++)
         forbiddenPoints[cPtr[2*(*it)+1]] = true;  // we know that a 0D mesh has a connectivity of the form [NORM_POINT1, i1, NORM_POINT1, i2, ...]
    }

  MCAuto<DataArrayDouble> appendedCoords(DataArrayDouble::New()); appendedCoords->alloc(0,1);//1 not 2 it is not a bug.
  const double *coords(_coords->begin());
  mcIdType *newciptr(newci->getPointer());
  for(mcIdType i=0;i<nbOfCells;i++,newciptr++,ciptr++)
    {
      if(Colinearize2DCell(coords,cptr+ciptr[0],cptr+ciptr[1],nbOfNodes,forbiddenPoints, /*out*/ newc,appendedCoords))
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
bool MEDCouplingUMesh::OrderPointsAlongLine(const double * coo, mcIdType startNode, mcIdType endNode,
                                            const mcIdType * c, const mcIdType * cI, const mcIdType *idsBg, const mcIdType *endBg,
                                            std::vector<mcIdType> & pointIds, std::vector<mcIdType> & hitSegs)
{
  using namespace std;

  const int SPACEDIM=3;
  typedef pair<double, mcIdType> PairDI;
  set< PairDI > x;
  for (const mcIdType * it = idsBg; it != endBg; ++it)
    {
      assert(c[cI[*it]] == INTERP_KERNEL::NORM_SEG2);
      mcIdType start = c[cI[*it]+1], end = c[cI[*it]+2];
      x.insert(make_pair(coo[start*SPACEDIM], start));  // take only X coords
      x.insert(make_pair(coo[end*SPACEDIM], end));
    }

  vector<PairDI> xx(x.begin(), x.end());
  sort(xx.begin(),xx.end());
  pointIds.reserve(xx.size());

  // Keep what is inside [startNode, endNode]:
  mcIdType go = 0;
  for (vector<PairDI>::const_iterator it=xx.begin(); it != xx.end(); ++it)
    {
      const mcIdType idx = (*it).second;
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

//  vector<mcIdType> pointIds2(pointIds.size()+2);
//  copy(pointIds.begin(), pointIds.end(), pointIds2.data()+1);
//  pointIds2[0] = startNode;
//  pointIds2[pointIds2.size()-1] = endNode;

  if (go == 2)
    reverse(pointIds.begin(), pointIds.end());

  // Now identify smaller segments that are not sub-divided - those won't need any further treatment:
  for (const mcIdType * it = idsBg; it != endBg; ++it)
    {
      mcIdType start = c[cI[*it]+1], end = c[cI[*it]+2];
      vector<mcIdType>::const_iterator itStart = find(pointIds.begin(), pointIds.end(), start);
      if (itStart == pointIds.end()) continue;
      vector<mcIdType>::const_iterator itEnd = find(pointIds.begin(), pointIds.end(), end);
      if (itEnd == pointIds.end())               continue;
      if (abs(distance(itEnd, itStart)) != 1)    continue;
      hitSegs.push_back(*it);   // segment is undivided.
    }

  return (pointIds.size() > 2); // something else apart start and end node
}

void MEDCouplingUMesh::ReplaceEdgeInFace(const mcIdType * sIdxConn, const mcIdType * sIdxConnE, mcIdType startNode, mcIdType endNode,
                                          const std::vector<mcIdType>& insidePoints, std::vector<mcIdType>& modifiedFace)
{
  using namespace std;
  size_t dst = distance(sIdxConn, sIdxConnE);
  modifiedFace.reserve(dst + insidePoints.size()-2);
  modifiedFace.resize(dst);
  copy(sIdxConn, sIdxConnE, modifiedFace.data());

  vector<mcIdType>::iterator shortEnd = modifiedFace.begin()+dst;
  vector<mcIdType>::iterator startPos = find(modifiedFace.begin(), shortEnd , startNode);
  if (startPos == shortEnd)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ReplaceEdgeInFace: internal error, should never happen!");
  vector<mcIdType>::iterator endPos = find(modifiedFace.begin(),shortEnd, endNode);
  if (endPos == shortEnd)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ReplaceEdgeInFace: internal error, should never happen!");
  size_t d = distance(startPos, endPos);
  if (d == 1 || d == (1-dst)) // don't use modulo, for neg numbers, result is implementation defined ...
    modifiedFace.insert(++startPos, ++insidePoints.begin(), --insidePoints.end());  // insidePoints also contains start and end node. Those don't need to be inserted.
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
 * \return DataArrayIdType  * - The list of cellIds in \a this that have been subdivided. If empty, nothing changed in \a this (as if it were a const method). The array is a newly allocated array
 *                           that the user is expected to deal with.
 *
 * \throw If \a this is not coherent.
 * \throw If \a this has not spaceDim equal to 3.
 * \throw If \a this has not meshDim equal to 3.
 * \sa MEDCouplingUMesh::mergeNodes, MEDCouplingUMesh::conformize2D, MEDCouplingUMesh::convertAllToPoly()
 */
DataArrayIdType *MEDCouplingUMesh::conformize3D(double eps)
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
  MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(0,1);

  {
    /*************************
     *  STEP 1  -- faces
     *************************/
    MCAuto<DataArrayIdType> descDNU(DataArrayIdType::New()),descIDNU(DataArrayIdType::New()),revDesc(DataArrayIdType::New()),revDescI(DataArrayIdType::New());
    MCAuto<MEDCouplingUMesh> mDesc(buildDescendingConnectivity(descDNU,descIDNU,revDesc,revDescI));
    const mcIdType *revDescIP(revDescI->getConstPointer()), *revDescP(revDesc->getConstPointer());
    const mcIdType *cDesc(mDesc->getNodalConnectivity()->begin()),*cIDesc(mDesc->getNodalConnectivityIndex()->begin());
    MCAuto<MEDCouplingSkyLineArray> connSlaDesc(MEDCouplingSkyLineArray::New(mDesc->getNodalConnectivityIndex(), mDesc->getNodalConnectivity()));

    // Build BBTree
    MCAuto<DataArrayDouble> bboxArr(mDesc->getBoundingBoxForBBTree(eps));
    const double *bbox(bboxArr->begin()); getCoords()->begin();
    mcIdType nDescCell=mDesc->getNumberOfCells();
    BBTree<SPACEDIM,mcIdType> myTree(bbox,0,0,nDescCell,-eps);
    // Surfaces - handle biggest first
    MCAuto<MEDCouplingFieldDouble> surfF = mDesc->getMeasureField(true);
    DataArrayDouble * surfs = surfF->getArray();
    // Normal field
    MCAuto<MEDCouplingFieldDouble> normalsF = mDesc->buildOrthogonalField();
    DataArrayDouble * normals = normalsF->getArray();
    const double * normalsP = normals->getConstPointer();

    // Sort faces by decreasing surface:
    vector< pair<double,mcIdType> > S;
    for(mcIdType i=0;i < surfs->getNumberOfTuples();i++)
      {
        pair<double,mcIdType> p = make_pair(surfs->begin()[i], i);
        S.push_back(p);
      }
    sort(S.rbegin(),S.rend()); // reverse sort
    vector<bool> hit(nDescCell);
    fill(hit.begin(), hit.end(), false);
    vector<mcIdType> hitPoly; // the final result: which 3D cells have been modified.

    for( vector<pair<double,mcIdType> >::const_iterator it = S.begin(); it != S.end(); it++)
      {
        mcIdType faceIdx = (*it).second;
        if (hit[faceIdx]) continue;

        vector<mcIdType> candidates, cands2;
        myTree.getIntersectingElems(bbox+faceIdx*2*SPACEDIM,candidates);
        // Keep only candidates whose normal matches the normal of current face
        for(vector<mcIdType>::const_iterator it2=candidates.begin();it2!=candidates.end();it2++)
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
        for (mcIdType ii = 0; ii < mPartRef->_coords->getNumberOfTuples(); ii++)
          rotation.transform_vector(cooPartRef+SPACEDIM*ii);
        for (mcIdType ii = 0; ii < mPartCand->_coords->getNumberOfTuples(); ii++)
          rotation.transform_vector(cooPartCand+SPACEDIM*ii);

        // Localize faces in 2D thanks to barycenters
        MCAuto<DataArrayDouble> baryPart = mPartCand->computeCellCenterOfMass();
        vector<std::size_t> compo; compo.push_back(2);
        MCAuto<DataArrayDouble> baryPartZ = baryPart->keepSelectedComponents(compo);
        MCAuto<DataArrayIdType> idsGoodPlane = baryPartZ->findIdsInRange(-eps, +eps);
        if (!idsGoodPlane->getNumberOfTuples())
          continue;

        baryPart = baryPart->selectByTupleId(*idsGoodPlane);

        compo[0] = 0; compo.push_back(1);
        MCAuto<DataArrayDouble> baryPartXY = baryPart->keepSelectedComponents(compo);
        mPartRef->changeSpaceDimension(2,0.0);
        MCAuto<DataArrayIdType> cc(DataArrayIdType::New()), ccI(DataArrayIdType::New());
        mPartRef->getCellsContainingPoints(baryPartXY->begin(), baryPartXY->getNumberOfTuples(), eps, cc, ccI);

        if (!cc->getNumberOfTuples())
          continue;
        MCAuto<DataArrayIdType> dsi = ccI->deltaShiftIndex();

        {
          MCAuto<DataArrayIdType> tmp = dsi->findIdsInRange(0, 2);
          if (tmp->getNumberOfTuples() != dsi->getNumberOfTuples())
            {
              ostringstream oss;
              oss << "MEDCouplingUMesh::conformize3D: Non expected non-conformity. Only simple (=partition-like) non-conformities are handled. Face #" << faceIdx << " violates this condition!";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }

        MCAuto<DataArrayIdType> ids = dsi->findIdsEqual(1);
        // Boundary face:
        if (!ids->getNumberOfTuples())
          continue;

        double checkSurf=0.0;
        const mcIdType * idsGoodPlaneP(idsGoodPlane->begin());
        for (const mcIdType * ii = ids->begin(); ii != ids->end(); ii++)
          {
            mcIdType faceIdx2 = cands2[idsGoodPlaneP[*ii]];
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
        vector<mcIdType> polyIndices, packsIds, facePack;
        for (mcIdType ii=revDescIP[faceIdx]; ii < revDescIP[faceIdx+1]; ii++)
          polyIndices.push_back(revDescP[ii]);
        ret->pushBackValsSilent(polyIndices.data(),polyIndices.data()+polyIndices.size());

        // Current face connectivity
        const mcIdType * sIdxConn = cDesc + cIDesc[faceIdx] + 1;
        const mcIdType * sIdxConnE = cDesc + cIDesc[faceIdx+1];
        connSla->findPackIds(polyIndices, sIdxConn, sIdxConnE, packsIds);
        // Deletion of old faces
        mcIdType jj=0;
        for (vector<mcIdType>::const_iterator it2=polyIndices.begin(); it2!=polyIndices.end(); ++it2, ++jj)
          {
            if (packsIds[jj] == -1)
              // The below should never happen - if a face is used several times, with a different layout of the nodes
              // it means that it is already conform, so it is *not* hit by the algorithm. The algorithm only hits
              // faces which are actually used only once, by a single cell. This is different for edges below.
              throw INTERP_KERNEL::Exception("MEDCouplingUMesh::conformize3D: Could not find face in connectivity! Internal error.");
            else
              connSla->deletePack(*it2, packsIds[jj]);
          }
        // Insertion of new faces:
        for (const mcIdType * ii = ids->begin(); ii != ids->end(); ii++)
          {
            // Build pack from the face to insert:
            mcIdType faceIdx2 = cands2[idsGoodPlane->getIJ(*ii,0)];
            mcIdType facePack2Sz;
            const mcIdType * facePack2 = connSlaDesc->getSimplePackSafePtr(faceIdx2, facePack2Sz); // contains the type!
            // Insert it in all hit polyhedrons:
            for (vector<mcIdType>::const_iterator it2=polyIndices.begin(); it2!=polyIndices.end(); ++it2)
              connSla->pushBackPack(*it2, facePack2+1, facePack2+facePack2Sz);  // without the type
          }
      }
  }  // end step1

  // Set back modified connectivity
  MCAuto<DataArrayIdType> cAuto; cAuto.takeRef(_nodal_connec);
  MCAuto<DataArrayIdType> cIAuto; cIAuto.takeRef(_nodal_connec_index);
  connSla->convertToPolyhedronConn(cAuto, cIAuto);

  {
    /************************
     *  STEP 2 -- edges
     ************************/
    // Now we have a face-conform mesh.

    // Recompute descending
    MCAuto<DataArrayIdType> desc(DataArrayIdType::New()),descI(DataArrayIdType::New()),revDesc(DataArrayIdType::New()),revDescI(DataArrayIdType::New());
    // Rebuild desc connectivity with orientation this time!!
    MCAuto<MEDCouplingUMesh> mDesc(buildDescendingConnectivity2(desc,descI,revDesc,revDescI));
    const mcIdType *revDescIP(revDescI->getConstPointer()), *revDescP(revDesc->getConstPointer());
    const mcIdType *descIP(descI->getConstPointer()), *descP(desc->getConstPointer());
    const mcIdType *cDesc(mDesc->getNodalConnectivity()->begin()),*cIDesc(mDesc->getNodalConnectivityIndex()->begin());
    MCAuto<DataArrayIdType> ciDeepC(mDesc->getNodalConnectivityIndex()->deepCopy());
    MCAuto<DataArrayIdType> cDeepC(mDesc->getNodalConnectivity()->deepCopy());
    MCAuto<MEDCouplingSkyLineArray> connSlaDesc(MEDCouplingSkyLineArray::New(ciDeepC, cDeepC));
    MCAuto<DataArrayIdType> desc2(DataArrayIdType::New()),descI2(DataArrayIdType::New()),revDesc2(DataArrayIdType::New()),revDescI2(DataArrayIdType::New());
    MCAuto<MEDCouplingUMesh> mDesc2 = mDesc->buildDescendingConnectivity(desc2,descI2,revDesc2,revDescI2);
//    std::cout << "writing!\n";
//    mDesc->writeVTK("/tmp/toto_desc_confInter.vtu");
//    mDesc2->writeVTK("/tmp/toto_desc2_confInter.vtu");
    const mcIdType *revDescIP2(revDescI2->getConstPointer()), *revDescP2(revDesc2->getConstPointer());
    const mcIdType *cDesc2(mDesc2->getNodalConnectivity()->begin()),*cIDesc2(mDesc2->getNodalConnectivityIndex()->begin());
    MCAuto<DataArrayDouble> bboxArr(mDesc2->getBoundingBoxForBBTree(eps));
    const double *bbox2(bboxArr->begin());
    mcIdType nDesc2Cell=mDesc2->getNumberOfCells();
    BBTree<SPACEDIM,mcIdType> myTree2(bbox2,0,0,nDesc2Cell,-eps);

    // Edges - handle longest first
    MCAuto<MEDCouplingFieldDouble> lenF = mDesc2->getMeasureField(true);
    DataArrayDouble * lens = lenF->getArray();

    // Sort edges by decreasing length:
    vector<pair<double,mcIdType> > S;
    for(mcIdType i=0;i < lens->getNumberOfTuples();i++)
      {
        pair<double,mcIdType> p = make_pair(lens->getIJ(i, 0), i);
        S.push_back(p);
      }
    sort(S.rbegin(),S.rend()); // reverse sort

    vector<bool> hit(nDesc2Cell);
    fill(hit.begin(), hit.end(), false);

    for( vector<pair<double,mcIdType> >::const_iterator it = S.begin(); it != S.end(); it++)
      {
        mcIdType eIdx = (*it).second;
        if (hit[eIdx])
          continue;

        vector<mcIdType> candidates, cands2;
        myTree2.getIntersectingElems(bbox2+eIdx*2*SPACEDIM,candidates);
        // Keep only candidates colinear with current edge
        double vCurr[3];
        mcIdType start = cDesc2[cIDesc2[eIdx]+1], end = cDesc2[cIDesc2[eIdx]+2];
        for (mcIdType i3=0; i3 < 3; i3++)  // TODO: use fillSonCellNodalConnectivity2 or similar?
          vCurr[i3] = coo[start*SPACEDIM+i3] - coo[end*SPACEDIM+i3];
        for(vector<mcIdType>::const_iterator it2=candidates.begin();it2!=candidates.end();it2++)
          {
            double vOther[3];
            mcIdType start2 = cDesc2[cIDesc2[*it2]+1], end2 = cDesc2[cIDesc2[*it2]+2];
            for (mcIdType i3=0; i3 < 3; i3++)
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
        mcIdType startNode = cDesc2[cIDesc2[eIdx]+1];
        mcIdType endNode = cDesc2[cIDesc2[eIdx]+2];
        INTERP_KERNEL::TranslationRotationMatrix rotation;
        INTERP_KERNEL::TranslationRotationMatrix::Rotate3DBipoint(coo+SPACEDIM*startNode, coo+SPACEDIM*endNode, rotation);
        MCAuto<MEDCouplingUMesh> mPartRef(mDesc2->buildPartOfMySelfSlice(eIdx, eIdx+1,1,false));  // false=zipCoords is called
        MCAuto<MEDCouplingUMesh> mPartCand(mDesc2->buildPartOfMySelf(&cands2[0], &cands2[0]+cands2.size(), true)); // true=zipCoords is called
        MCAuto<DataArrayIdType> nodeMap(mPartCand->zipCoordsTraducer());
        mcIdType nbElemsNotM1;
        {
          MCAuto<DataArrayIdType> tmp(nodeMap->findIdsNotEqual(-1));
          nbElemsNotM1 = tmp->getNbOfElems();
        }
        MCAuto<DataArrayIdType>  nodeMapInv = nodeMap->invertArrayO2N2N2O(nbElemsNotM1);
        double * cooPartRef(mPartRef->_coords->getPointer());
        double * cooPartCand(mPartCand->_coords->getPointer());
        for (mcIdType ii = 0; ii < mPartRef->_coords->getNumberOfTuples(); ii++)
          rotation.transform_vector(cooPartRef+SPACEDIM*ii);
        for (mcIdType ii = 0; ii < mPartCand->_coords->getNumberOfTuples(); ii++)
          rotation.transform_vector(cooPartCand+SPACEDIM*ii);


        // Eliminate all edges for which y or z is not null
        MCAuto<DataArrayDouble> baryPart = mPartCand->computeCellCenterOfMass();
        vector<std::size_t> compo; compo.push_back(1);
        MCAuto<DataArrayDouble> baryPartY = baryPart->keepSelectedComponents(compo);
        compo[0] = 2;
        MCAuto<DataArrayDouble> baryPartZ = baryPart->keepSelectedComponents(compo);
        MCAuto<DataArrayIdType> idsGoodLine1 = baryPartY->findIdsInRange(-eps, +eps);
        MCAuto<DataArrayIdType> idsGoodLine2 = baryPartZ->findIdsInRange(-eps, +eps);
        MCAuto<DataArrayIdType> idsGoodLine = idsGoodLine1->buildIntersection(idsGoodLine2);
        if (!idsGoodLine->getNumberOfTuples())
          continue;

        // Now the ordering along the Ox axis:
        std::vector<mcIdType> insidePoints, hitSegs;
        bool isSplit = OrderPointsAlongLine(mPartCand->_coords->getConstPointer(), nodeMap->begin()[startNode], nodeMap->begin()[endNode],
            mPartCand->getNodalConnectivity()->begin(), mPartCand->getNodalConnectivityIndex()->begin(),
            idsGoodLine->begin(), idsGoodLine->end(),
            /*out*/insidePoints, hitSegs);
        // Optim: smaller segments completely included in eIdx and not split won't need any further treatment:
        for (vector<mcIdType>::const_iterator its=hitSegs.begin(); its != hitSegs.end(); ++its)
          hit[cands2[*its]] = true;

        if (!isSplit)  // current segment remains in one piece
          continue;

        // Get original node IDs in global coords array
        for (std::vector<mcIdType>::iterator iit = insidePoints.begin(); iit!=insidePoints.end(); ++iit)
          *iit = nodeMapInv->begin()[*iit];

        vector<mcIdType> polyIndices, packsIds, facePack;
        // For each face implying this edge
        for (mcIdType ii=revDescIP2[eIdx]; ii < revDescIP2[eIdx+1]; ii++)
          {
            mcIdType faceIdx = revDescP2[ii];
            // each cell where this face is involved connectivity will be modified:
            ret->pushBackValsSilent(revDescP + revDescIP[faceIdx], revDescP + revDescIP[faceIdx+1]);

            // Current face connectivity
            const mcIdType * sIdxConn = cDesc + cIDesc[faceIdx] + 1;
            const mcIdType * sIdxConnE = cDesc + cIDesc[faceIdx+1];

            vector<mcIdType> modifiedFace;
            ReplaceEdgeInFace(sIdxConn, sIdxConnE, startNode, endNode, insidePoints, /*out*/modifiedFace);
            modifiedFace.insert(modifiedFace.begin(), INTERP_KERNEL::NORM_POLYGON);
            connSlaDesc->replaceSimplePack(faceIdx, modifiedFace.data(), modifiedFace.data()+modifiedFace.size());
          }
      }

    // Rebuild 3D connectivity from descending:
    MCAuto<MEDCouplingSkyLineArray> newConn(MEDCouplingSkyLineArray::New());
    MCAuto<DataArrayIdType> superIdx(DataArrayIdType::New());  superIdx->alloc(getNumberOfCells()+1); superIdx->fillWithValue(0);
    MCAuto<DataArrayIdType> idx(DataArrayIdType::New());       idx->alloc(1);                         idx->fillWithValue(0);
    MCAuto<DataArrayIdType> vals(DataArrayIdType::New());      vals->alloc(0);
    newConn->set3(superIdx, idx, vals);
    mcIdType nbCells=getNumberOfCells();
    for(mcIdType ii = 0; ii < nbCells; ii++)
      for (mcIdType jj=descIP[ii]; jj < descIP[ii+1]; jj++)
        {
          mcIdType sz, faceIdx = abs(descP[jj])-1;
          bool orient = descP[jj]>0;
          const mcIdType * p = connSlaDesc->getSimplePackSafePtr(faceIdx, sz);
          if (orient)
            newConn->pushBackPack(ii, p+1, p+sz);  // +1 to skip type
          else
            {
              vector<mcIdType> rev(sz-1);
              for (mcIdType kk=0; kk<sz-1; kk++) rev[kk]=*(p+sz-kk-1);
              newConn->pushBackPack(ii, rev.data(), rev.data()+sz-1);
            }
        }
    // And finally:
    newConn->convertToPolyhedronConn(cAuto, cIAuto);
  } // end step2

  ret = ret->buildUniqueNotSorted();
  return ret.retn();
}


