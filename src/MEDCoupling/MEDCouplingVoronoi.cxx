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

#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "Interpolation2D.txx"
#include "Interpolation3DSurf.hxx"

using namespace MEDCoupling;

Voronizer::~Voronizer()
{
}

int Voronizer1D::getDimension() const
{
  return 1;
}

int Voronizer2D::getDimension() const
{
  return 2;
}

int Voronizer3D::getDimension() const
{
  return 3;
}

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


MCAuto<MEDCouplingUMesh> MergeVorCells2D(MEDCouplingUMesh *p, double eps, bool isZip)
{
  MCAuto<DataArrayInt> edgeToKeep;
  MCAuto<MEDCouplingUMesh> p0;
  {
    MCAuto<DataArrayInt> d(DataArrayInt::New()),di(DataArrayInt::New()),rd(DataArrayInt::New()),rdi(DataArrayInt::New());
    p0=p->buildDescendingConnectivity(d,di,rd,rdi);
    MCAuto<DataArrayInt> dsi(rdi->deltaShiftIndex());
    edgeToKeep=dsi->findIdsEqual(1);
  }
  MCAuto<MEDCouplingUMesh> skinOfRes(p0->buildPartOfMySelf(edgeToKeep->begin(),edgeToKeep->end()));
  if(isZip)
    {
      skinOfRes->zipCoords();
      if(skinOfRes->getNumberOfCells()!=skinOfRes->getNumberOfNodes())
        throw INTERP_KERNEL::Exception("MergeVorCells : result of merge looks bad !");
    }
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
  return MergeVorCells2D(p,eps,true);
}

/*!
 * suppress additional sub points on edges
 */
MCAuto<MEDCouplingUMesh> SimplifyPolygon(const MEDCouplingUMesh *m, double eps)
{
  if(m->getNumberOfCells()!=1)
    throw INTERP_KERNEL::Exception("SimplifyPolygon : internal error !");
  const int *conn(m->getNodalConnectivity()->begin()),*conni(m->getNodalConnectivityIndex()->begin());
  int nbPtsInPolygon(conni[1]-conni[0]-1);
  const double *coo(m->getCoords()->begin());
  std::vector<int> resConn;
  for(int i=0;i<nbPtsInPolygon;i++)
    {
      int prev(conn[(i+nbPtsInPolygon-1)%nbPtsInPolygon+1]),current(conn[i%nbPtsInPolygon+1]),zeNext(conn[(i+1)%nbPtsInPolygon+1]);
      double a[3]={
        coo[3*prev+0]-coo[3*current+0],
        coo[3*prev+1]-coo[3*current+1],
        coo[3*prev+2]-coo[3*current+2],
      },b[3]={
        coo[3*current+0]-coo[3*zeNext+0],
        coo[3*current+1]-coo[3*zeNext+1],
        coo[3*current+2]-coo[3*zeNext+2],
      };
      double c[3]={a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
      if(sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2])>eps)
        resConn.push_back(current);
    }
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::New("",2));
  ret->setCoords(m->getCoords());
  ret->allocateCells();
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYGON,resConn.size(),&resConn[0]);
  return ret;
}

MCAuto<MEDCouplingUMesh> MergeVorCells3D(const std::vector< MCAuto<MEDCouplingUMesh> >& vcs, double eps)
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
  MCAuto<DataArrayDouble> eqn(skinOfRes->computePlaneEquationOf3DFaces());
  MCAuto<DataArrayInt> comm,commI;
  {
    DataArrayInt *a(0),*b(0);
    eqn->findCommonTuples(eps,0,a,b);
    comm=a; commI=b;
    //comm=DataArrayInt::New(); comm->alloc(0,1); commI=DataArrayInt::New(); commI->alloc(1,1); commI->setIJ(0,0,0);
  }
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::New("",3));
  ret->setCoords(skinOfRes->getCoords());
  ret->allocateCells();
  std::vector<int> conn;
  int jj(0);
  for(int i=0;i<commI->getNumberOfTuples()-1;i++,jj++)
    {
      if(jj!=0)
        conn.push_back(-1);
      MCAuto<MEDCouplingUMesh> tmp(skinOfRes->buildPartOfMySelf(comm->begin()+commI->getIJ(i,0),comm->begin()+commI->getIJ(i+1,0),true));
      MCAuto<MEDCouplingUMesh> tmp2;
      if(commI->getIJ(i+1,0)-commI->getIJ(i,0)==1)
        tmp2=tmp;
      else
        tmp2=MergeVorCells2D(tmp,eps,false);
      tmp2=SimplifyPolygon(tmp2,eps);
      const int *cPtr(tmp2->getNodalConnectivity()->begin()),*ciPtr(tmp2->getNodalConnectivityIndex()->begin());
      conn.insert(conn.end(),cPtr+1,cPtr+ciPtr[1]);
    }
  MCAuto<DataArrayInt> remain(comm->buildComplement(skinOfRes->getNumberOfCells()));
  {
    MCAuto<MEDCouplingUMesh> tmp(skinOfRes->buildPartOfMySelf(remain->begin(),remain->end(),true));
    const int *cPtr(tmp->getNodalConnectivity()->begin()),*ciPtr(tmp->getNodalConnectivityIndex()->begin());
    for(int i=0;i<remain->getNumberOfTuples();i++,jj++)
      {
        if(jj!=0)
          conn.push_back(-1);
        conn.insert(conn.end(),cPtr+ciPtr[i]+1,cPtr+ciPtr[i+1]);
      }
  }
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,conn.size(),&conn[0]);
  return ret;
}

MCAuto<MEDCouplingUMesh> MergeVorCells1D(const std::vector< MCAuto<MEDCouplingUMesh> >& vcs, double eps)
{
  static const int CONN_SEG2_DFT[2]={0,1};
  if(vcs.empty())
    throw INTERP_KERNEL::Exception("MergeVorCells1D : internal error 1 !");
  if(vcs.size()==1)
    return vcs[0];
  if(vcs.size()>2)
    throw INTERP_KERNEL::Exception("MergeVorCells1D : internal error 2 !");
  double a0,b0,a1,b1;
  {
    const int *connPtr(vcs[0]->getNodalConnectivity()->begin());
    const double *coordPtr(vcs[0]->getCoords()->begin());
    a0=coordPtr[connPtr[1]]; b0=coordPtr[connPtr[2]];
  }
  {
    const int *connPtr(vcs[1]->getNodalConnectivity()->begin());
    const double *coordPtr(vcs[1]->getCoords()->begin());
    a1=coordPtr[connPtr[1]]; b1=coordPtr[connPtr[2]];
  }
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::New("",1)); ret->allocateCells(); ret->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,CONN_SEG2_DFT);
  MCAuto<DataArrayDouble> coo(DataArrayDouble::New()); coo->alloc(2,1); ret->setCoords(coo);
  if(fabs(b0-a1)<eps)
    { coo->setIJ(0,0,a0); coo->setIJ(1,0,b1); }
  else if(fabs(b1-a0)<eps)
    { coo->setIJ(0,0,b0); coo->setIJ(1,0,a1); }
  return ret;
}

MCAuto<MEDCouplingUMesh> MEDCoupling::Voronizer1D::doIt(const MEDCouplingUMesh *m, const DataArrayDouble *points, double eps) const
{
  static const int CONN_SEG2_DFT[2]={0,1};
  if(!m || !points)
    throw INTERP_KERNEL::Exception("Voronoize1D : null pointer !");
  m->checkConsistencyLight();
  points->checkAllocated();
  if(m->getMeshDimension()!=1 || m->getSpaceDimension()!=1 || points->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("Voronoize1D : spacedim must be equal to 1 and meshdim also equal to 1 !");
  if(m->getNumberOfCells()!=1)
    throw INTERP_KERNEL::Exception("Voronoize1D : mesh is expected to have only one cell !");
  int nbPts(points->getNumberOfTuples());
  if(nbPts<1)
    throw INTERP_KERNEL::Exception("Voronoize1D : at least one point expected !");
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
      const double *pt(pts+i);
      vorTess->getCellsContainingPoint(pt,eps,polygsToIterOn);
      if(polygsToIterOn.empty())
        throw INTERP_KERNEL::Exception("Voronoize1D : a point is outside domain !");
      if(polygsToIterOn.size()>2)
        throw INTERP_KERNEL::Exception("Voronoize1D : overlap of points !");
      std::vector< MCAuto<MEDCouplingUMesh> > newVorCells;
      for(std::vector<int>::const_iterator it=polygsToIterOn.begin();it!=polygsToIterOn.end();it++)
        {
          int poly(*it);
          //
          double seed(pts[poly]),zept(*pt);
          double mid((seed+zept)/2.);
          //
          MCAuto<MEDCouplingUMesh> tile(l0[poly]);
          tile->zipCoords();
          double a,b;
          {
            const int *connPtr(tile->getNodalConnectivity()->begin());
            const double *coordPtr(tile->getCoords()->begin());
            a=coordPtr[connPtr[1]]; b=coordPtr[connPtr[2]];
          }
          double pol0[2],pol1[2];
          MCAuto<DataArrayDouble> t0(DataArrayDouble::New()); t0->alloc(3,1); t0->setIJ(0,0,zept); t0->setIJ(1,0,mid); t0->setIJ(2,0,seed);
          t0->applyLin(1.,-a);
          if(t0->isMonotonic(true,eps))
            { pol0[0]=a; pol0[1]=mid; pol1[0]=mid; pol1[1]=b; }
          else
            { pol1[0]=a; pol1[1]=mid; pol0[0]=mid; pol0[1]=b; }
          MCAuto<MEDCouplingUMesh> modifiedCell(MEDCouplingUMesh::New("",1)); modifiedCell->allocateCells();
          MCAuto<DataArrayDouble> coo1(DataArrayDouble::New()); coo1->alloc(2,1); coo1->setIJ(0,0,pol1[0]); coo1->setIJ(1,0,pol1[1]);
          modifiedCell->setCoords(coo1); modifiedCell->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,CONN_SEG2_DFT);
          //
          MCAuto<MEDCouplingUMesh> newVorCell(MEDCouplingUMesh::New("",1)); newVorCell->allocateCells();
          MCAuto<DataArrayDouble> coo2(DataArrayDouble::New()); coo2->alloc(2,1); coo2->setIJ(0,0,pol0[0]); coo2->setIJ(1,0,pol0[1]);
          newVorCell->setCoords(coo2); newVorCell->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,CONN_SEG2_DFT);
          //
          l0[poly]=modifiedCell;
          newVorCells.push_back(newVorCell);
        }
      l0.push_back(MergeVorCells1D(newVorCells,eps));
    }
  std::vector< const MEDCouplingUMesh * > l0Bis(VecAutoToVecOfCstPt(l0));
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::MergeUMeshes(l0Bis));
  {
    bool dummy; int dummy2;
    MCAuto<DataArrayInt> dummy3(ret->mergeNodes(eps,dummy,dummy2));
  }
  return ret;
}

MCAuto<MEDCouplingUMesh> MEDCoupling::Voronizer2D::doIt(const MEDCouplingUMesh *m, const DataArrayDouble *points, double eps) const
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
          for(int j=0;j<newCoords->getNumberOfTuples();j++,cPtr+=2)
            {
              std::set<int> zeCandidates;
              {
                std::vector<int> zeCandidatesTmp;
                vorTess->getCellsContainingPoint(cPtr,eps,zeCandidatesTmp);
                zeCandidates.insert(zeCandidatesTmp.begin(),zeCandidatesTmp.end());
              }
              std::set<int> tmp2,newElementsToDo;
              std::set_difference(zeCandidates.begin(),zeCandidates.end(),elemsDone.begin(),elemsDone.end(),std::inserter(tmp2,tmp2.begin()));
              std::set_union(elemsToDo.begin(),elemsToDo.end(),tmp2.begin(),tmp2.end(),std::inserter(newElementsToDo,newElementsToDo.begin()));
              elemsToDo=newElementsToDo;
            }
          newVorCells.push_back(newVorCell);
        }
      l0.push_back(MergeVorCells(newVorCells,eps));
    }
  std::vector< const MEDCouplingUMesh * > l0Bis(VecAutoToVecOfCstPt(l0));
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::MergeUMeshes(l0Bis));
  {
    bool dummy; int dummy2;
    MCAuto<DataArrayInt> dummy3(ret->mergeNodes(eps,dummy,dummy2));
  }
  return ret;
}

MCAuto<MEDCouplingUMesh> Split3DCellInParts(const MEDCouplingUMesh *m, const double pt[3], const double seed[3], double eps, int tmp[2])
{
  if(m->getMeshDimension()!=3 || m->getSpaceDimension()!=3 || m->getNumberOfCells()!=1)
    throw INTERP_KERNEL::Exception("Split3DCellInParts : expecting a 3D with exactly one cell !");
  double middle[3]={(pt[0]+seed[0])/2.,(pt[1]+seed[1])/2.,(pt[2]+seed[2])/2.};
  double vec[3]={pt[0]-seed[0],pt[1]-seed[1],pt[2]-seed[2]};
  MCAuto<MEDCouplingUMesh> res(m->clipSingle3DCellByPlane(middle,vec,eps));
  return res;
}

MCAuto<MEDCouplingUMesh> MEDCoupling::Voronizer3D::doIt(const MEDCouplingUMesh *m, const DataArrayDouble *points, double eps) const
{
  double eps2(1.-sqrt(eps));// 2nd eps for interpolation. Here the eps is computed to feet cos(eps) ~ 1-eps^2
  if(!m || !points)
    throw INTERP_KERNEL::Exception("Voronoize3D : null pointer !");
  m->checkConsistencyLight();
  points->checkAllocated();
  if(m->getMeshDimension()!=3 || m->getSpaceDimension()!=3 || points->getNumberOfComponents()!=3)
    throw INTERP_KERNEL::Exception("Voronoize3D : spacedim must be equal to 3 and meshdim also equal to 3 !");
  if(m->getNumberOfCells()!=1)
    throw INTERP_KERNEL::Exception("Voronoize3D : mesh is expected to have only one cell !");
  int nbPts(points->getNumberOfTuples());
  if(nbPts<1)
    throw INTERP_KERNEL::Exception("Voronoize3D : at least one point expected !");
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
      const double *pt(pts+i*3);
      vorTess->getCellsContainingPoint(pt,eps,polygsToIterOn);
      if(polygsToIterOn.size()<1)
        throw INTERP_KERNEL::Exception("Voronoize3D : presence of a point outside the given cell !");
      std::vector< MCAuto<MEDCouplingUMesh> > newVorCells;
      for(int poly=0;poly<vorTess->getNumberOfCells();poly++)
        {
          const double *seed(pts+3*poly);
          MCAuto<MEDCouplingUMesh> tile(l0[poly]);
          tile->zipCoords();
          int tmp[2];
          MCAuto<MEDCouplingUMesh> cells;
          try
            {
              cells=Split3DCellInParts(tile,pt,seed,eps,tmp);
            }
          catch(INTERP_KERNEL::Exception& e)
            {
              continue;
            }
          MCAuto<MEDCouplingUMesh> newVorCell(cells->buildPartOfMySelfSlice(1,2,1,true));
          newVorCell->zipCoords();
          MCAuto<MEDCouplingUMesh> modifiedCell(cells->buildPartOfMySelfSlice(0,1,1,true));
          modifiedCell->zipCoords();
          newVorCells.push_back(newVorCell);
          l0[poly]=modifiedCell;
        }
      l0.push_back(MergeVorCells3D(newVorCells,eps));
    }
  std::vector< const MEDCouplingUMesh * > l0Bis(VecAutoToVecOfCstPt(l0));
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::MergeUMeshes(l0Bis));
  {
    bool dummy; int dummy2;
    MCAuto<DataArrayInt> dummy3(ret->mergeNodes(eps,dummy,dummy2));
  }
  return ret;
}
