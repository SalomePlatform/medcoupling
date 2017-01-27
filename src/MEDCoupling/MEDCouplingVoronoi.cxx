// Copyright (C) 2007-2017  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

#include "MEDCouplingVoronoi.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MCAuto.txx"

using namespace MEDCoupling;

MCAuto<MEDCouplingUMesh> ComputeBigCellFrom(const double pt1[2], const double pt2[2], const std::vector<double>& bbox, double eps)
{
  static const double FACT=1.2;
  MCAuto<MEDCouplingCMesh> m(MEDCouplingCMesh::New());
  MCAuto<DataArrayDouble> arr1(DataArrayDouble::New()); arr1->alloc(2,1) ; arr1->setIJ(0,0,bbox[0]); arr1->setIJ(1,0,bbox[1]);
  MCAuto<DataArrayDouble> arr2(DataArrayDouble::New()); arr2->alloc(2,1) ; arr2->setIJ(0,0,bbox[2]); arr2->setIJ(1,0,bbox[3]);
  m->setCoords(arr1,arr2);
  static const double PT[2]={0.,0.};
  m->scale(PT,FACT);
  MCAuto<MEDCouplingUMesh> mu(m->buildUnstructured());
  double l(std::max(bbox[1]-bbox[0],bbox[3]-bbox[2]));
  double middle[2]={(pt1[0]+pt2[0])/2.,(pt1[1]+pt2[1])/2.};
  double v[2]={pt1[0],pt1[1]};
  DataArrayDouble::Rotate2DAlg(middle,M_PI/2,1,v,v);
  v[0]=middle[0]-v[0]; v[1]=middle[1]-v[1];
  {
    double nor(sqrt(v[0]*v[0]+v[1]*v[1]));
    v[0]/=nor; v[1]/=nor;
  }
  MCAuto<MEDCouplingUMesh> line(MEDCouplingUMesh::New("line",1));
  {
    MCAuto<DataArrayDouble> coo(DataArrayDouble::New()); coo->alloc(2,2);
    coo->setIJ(0,0,middle[0]-2.*l*v[0]); coo->setIJ(0,1,middle[1]-2.*l*v[1]); coo->setIJ(1,0,middle[0]+2.*l*v[0]); coo->setIJ(1,1,middle[1]+2.*l*v[1]);
    line->setCoords(coo);
  }
  line->allocateCells();
  static const int CONN[2]={0,1};
  line->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,CONN);
  MCAuto<MEDCouplingUMesh> sp2,sp1;
  {
    DataArrayInt *cellNb1(0),*cellNb2(0);
    MEDCouplingUMesh *sp2Pt(0),*sp1Pt(0);
    MEDCouplingUMesh::Intersect2DMeshWith1DLine(mu,line,eps,sp2Pt,sp1Pt,cellNb1,cellNb2);
    sp1=sp1Pt; sp2=sp2Pt;
    MCAuto<DataArrayInt> cellNb10(cellNb1),cellNb20(cellNb2);
  }
  std::vector<int> ccp;
  sp2->getCellsContainingPoint(pt1,eps,ccp);
  if(ccp.size()!=1)
    throw INTERP_KERNEL::Exception("ComputeBigCellFrom : expected single element !");
  MCAuto<MEDCouplingUMesh> ret(sp2->buildPartOfMySelfSlice(ccp[0],ccp[0]+1,1,true));
  ret->zipCoords();
  return ret;
}

MCAuto<MEDCouplingUMesh> MergeVorCells(const std::vector< MCAuto<MEDCouplingUMesh> >& vcs, double eps)
{
  std::size_t sz(vcs.size());
  if(sz<1)
    throw INTERP_KERNEL::Exception("MergeVorCells : len of input vec expected to be >= 1 !");
  if(sz==1)
    return vcs[0];
  MCAuto<MEDCouplingUMesh> p;
  {
    std::vector< const MEDCouplingUMesh * > vcsBis(VecAutoToVecOfCstPt(vcs));
    p=MEDCouplingUMesh::MergeUMeshes(vcsBis);
  }
  p->zipCoords();
  {
    bool dummy; int dummy2;
    MCAuto<DataArrayInt> dummy3(p->mergeNodes(eps,dummy,dummy2));
  }
  MCAuto<DataArrayInt> edgeToKeep;
  MCAuto<MEDCouplingUMesh> p0;
  {
    MCAuto<DataArrayInt> d(DataArrayInt::New()),di(DataArrayInt::New()),rd(DataArrayInt::New()),rdi(DataArrayInt::New());
    p0=p->buildDescendingConnectivity(d,di,rd,rdi);
    MCAuto<DataArrayInt> dsi(rdi->deltaShiftIndex());
    edgeToKeep=dsi->findIdsEqual(1);
  }
  MCAuto<MEDCouplingUMesh> skinOfRes(p0->buildPartOfMySelf(edgeToKeep->begin(),edgeToKeep->end()));
  skinOfRes->zipCoords();
  if(skinOfRes->getNumberOfCells()!=skinOfRes->getNumberOfNodes())
    throw INTERP_KERNEL::Exception("MergeVorCells : result of merge looks bad !");
  MCAuto<DataArrayInt> d(skinOfRes->orderConsecutiveCells1D());
  MCAuto<MEDCoupling1SGTUMesh> skinOfRes2;
  {
    MCAuto<MEDCouplingUMesh> part(skinOfRes->buildPartOfMySelf(d->begin(),d->end()));
    skinOfRes2=MEDCoupling1SGTUMesh::New(part);
  }
  MCAuto<DataArrayInt> c(skinOfRes2->getNodalConnectivity()->deepCopy());
  c->circularPermutation(1);
  c->rearrange(2);
  std::vector< MCAuto<DataArrayInt> > vdi(c->explodeComponents());
  if(!vdi[0]->isEqual(*vdi[1]))
    throw INTERP_KERNEL::Exception("MergeVorCells : internal error !");
  MCAuto<MEDCouplingUMesh> m(MEDCouplingUMesh::New("",2));
  m->setCoords(skinOfRes2->getCoords());
  m->allocateCells();
  m->insertNextCell(INTERP_KERNEL::NORM_POLYGON,vdi[0]->getNumberOfTuples(),vdi[0]->begin());
  return m;
}

MCAuto<MEDCouplingUMesh> MEDCoupling::Voronoize2D(const MEDCouplingUMesh *m, const DataArrayDouble *points, double eps)
{
  if(!m || !points)
    throw INTERP_KERNEL::Exception("Voronoize2D : null pointer !");
  m->checkConsistencyLight();
  points->checkAllocated();
  if(m->getMeshDimension()!=2 || m->getSpaceDimension()!=2 || points->getNumberOfComponents()!=2)
    throw INTERP_KERNEL::Exception("Voronoize2D : spacedim must be equal to 2 and meshdim also equal to 2 !");
  if(m->getNumberOfCells()!=1)
    throw INTERP_KERNEL::Exception("Voronoize2D : mesh is expected to have only one cell !");
  int nbPts(points->getNumberOfTuples());
  if(nbPts<1)
    throw INTERP_KERNEL::Exception("Voronoize2D : at least one point expected !");
  std::vector<double> bbox(4);
  m->getBoundingBox(&bbox[0]);
  std::vector< MCAuto<MEDCouplingUMesh> > l0(1,MCAuto<MEDCouplingUMesh>(m->deepCopy()));
  const double *pts(points->begin());
  for(int i=1;i<nbPts;i++)
    {
      MCAuto<MEDCouplingUMesh> vorTess;
      {
        std::vector< const MEDCouplingUMesh * > l0Bis(VecAutoToVecOfCstPt(l0));
        vorTess=MEDCouplingUMesh::MergeUMeshes(l0Bis);
      }
      {
        bool dummy;
        int newNbNodes;
        MCAuto<DataArrayInt> dummy3(vorTess->mergeNodes(eps,dummy,newNbNodes));
      }
      std::vector<int> polygsToIterOn;
      const double *pt(pts+i*2);
      vorTess->getCellsContainingPoint(pt,eps,polygsToIterOn);
      if(polygsToIterOn.size()<1)
        throw INTERP_KERNEL::Exception("Voronoize2D : presence of a point outside the given cell !");
      std::set<int> elemsToDo,elemsDone; elemsToDo.insert(polygsToIterOn[0]);
      std::vector< MCAuto<MEDCouplingUMesh> > newVorCells;
      while(!elemsToDo.empty())
        {
          int poly(*elemsToDo.begin()); elemsToDo.erase(elemsToDo.begin()); elemsDone.insert(poly);
          const double *seed(pts+2*poly);
          MCAuto<MEDCouplingUMesh> cell(ComputeBigCellFrom(pt,seed,bbox,eps));
          MCAuto<MEDCouplingUMesh> tile(l0[poly]);
          tile->zipCoords();
          MCAuto<MEDCouplingUMesh> a;
          MCAuto<DataArrayInt> b,c;
          {
            DataArrayInt *bPtr(0),*cPtr(0);
            a=MEDCouplingUMesh::Intersect2DMeshes(tile,cell,eps,bPtr,cPtr);
            b=bPtr; c=cPtr;
          }
          MCAuto<DataArrayInt> part(c->findIdsEqual(-1));
          if(part->getNumberOfTuples()!=1)
            throw INTERP_KERNEL::Exception("Voronoize2D : internal error");
          MCAuto<MEDCouplingUMesh> newVorCell;
          {
            MCAuto<DataArrayInt> tmp(part->buildComplement(a->getNumberOfCells()));
            newVorCell=a->buildPartOfMySelf(tmp->begin(),tmp->end());
          }
          newVorCell->zipCoords();
          MCAuto<MEDCouplingUMesh> modifiedCell(a->buildPartOfMySelf(part->begin(),part->end()));
          modifiedCell->zipCoords();
          l0[poly]=modifiedCell;
          //
          MCAuto<DataArrayInt> ids;
          {
            DataArrayInt *tmp(0);
            bool sta(a->getCoords()->areIncludedInMe(cell->getCoords(),eps,tmp));
            ids=tmp;
            if(!sta)
              throw INTERP_KERNEL::Exception("Voronoize2D : internal error 2 !");
          }
          MCAuto<DataArrayDouble> newCoords;
          {
            MCAuto<DataArrayInt> tmp(ids->buildComplement(a->getNumberOfNodes()));
            newCoords=a->getCoords()->selectByTupleId(tmp->begin(),tmp->end());
          }
          const double *cPtr(newCoords->begin());
          for(int i=0;i<newCoords->getNumberOfTuples();i++,cPtr+=2)
            {
              std::set<int> zeCandidates;
              {
                std::vector<int> zeCandidatesTmp;
                vorTess->getCellsContainingPoint(cPtr,eps,zeCandidatesTmp);
                zeCandidates.insert(zeCandidatesTmp.begin(),zeCandidatesTmp.end());
              }
              std::set<int> tmp,newElementsToDo;
              std::set_difference(zeCandidates.begin(),zeCandidates.end(),elemsDone.begin(),elemsDone.end(),std::inserter(tmp,tmp.begin()));
              std::set_union(elemsToDo.begin(),elemsToDo.end(),tmp.begin(),tmp.end(),std::inserter(newElementsToDo,newElementsToDo.begin()));
              elemsToDo=newElementsToDo;
            }
          newVorCells.push_back(newVorCell);
        }
      l0.push_back(MergeVorCells(newVorCells,eps));
    }
  std::vector< const MEDCouplingUMesh * > l0Bis(VecAutoToVecOfCstPt(l0));
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::MergeUMeshes(l0Bis));
  return ret;
}
