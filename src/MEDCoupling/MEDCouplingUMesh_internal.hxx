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
// Author : Anthony Geay (EdF)

using namespace MEDCoupling;

class MinusOneSonsGenerator
{
public:
  MinusOneSonsGenerator(const INTERP_KERNEL::CellModel& cm):_cm(cm) { }
  unsigned getNumberOfSons2(const int *conn, int lgth) const { return _cm.getNumberOfSons2(conn,lgth); }
  unsigned fillSonCellNodalConnectivity2(int sonId, const int *nodalConn, int lgth, int *sonNodalConn, INTERP_KERNEL::NormalizedCellType& typeOfSon) const { return _cm.fillSonCellNodalConnectivity2(sonId,nodalConn,lgth,sonNodalConn,typeOfSon); }
  static const int DELTA=1;
private:
  const INTERP_KERNEL::CellModel& _cm;
};

class MinusOneSonsGeneratorBiQuadratic
{
public:
  MinusOneSonsGeneratorBiQuadratic(const INTERP_KERNEL::CellModel& cm):_cm(cm) { }
  unsigned getNumberOfSons2(const int *conn, int lgth) const { return _cm.getNumberOfSons2(conn,lgth); }
  unsigned fillSonCellNodalConnectivity2(int sonId, const int *nodalConn, int lgth, int *sonNodalConn, INTERP_KERNEL::NormalizedCellType& typeOfSon) const { return _cm.fillSonCellNodalConnectivity4(sonId,nodalConn,lgth,sonNodalConn,typeOfSon); }
  static const int DELTA=1;
private:
  const INTERP_KERNEL::CellModel& _cm;
};

class MinusTwoSonsGenerator
{
public:
  MinusTwoSonsGenerator(const INTERP_KERNEL::CellModel& cm):_cm(cm) { }
  unsigned getNumberOfSons2(const int *conn, int lgth) const { return _cm.getNumberOfEdgesIn3D(conn,lgth); }
  unsigned fillSonCellNodalConnectivity2(int sonId, const int *nodalConn, int lgth, int *sonNodalConn, INTERP_KERNEL::NormalizedCellType& typeOfSon) const { return _cm.fillSonEdgesNodalConnectivity3D(sonId,nodalConn,lgth,sonNodalConn,typeOfSon); }
  static const int DELTA=2;
private:
  const INTERP_KERNEL::CellModel& _cm;
};

class MicroEdgesGenerator2D
{
public:
  MicroEdgesGenerator2D(const INTERP_KERNEL::CellModel& cm):_cm(cm) { }
  unsigned getNumberOfSons2(const int *conn, int lgth) const { return _cm.getNumberOfMicroEdges(); }
  unsigned fillSonCellNodalConnectivity2(int sonId, const int *nodalConn, int lgth, int *sonNodalConn, INTERP_KERNEL::NormalizedCellType& typeOfSon) const { return _cm.fillMicroEdgeNodalConnectivity(sonId,nodalConn,sonNodalConn,typeOfSon); }
  static const int DELTA=1;
private:
  const INTERP_KERNEL::CellModel& _cm;
};

class MicroEdgesGenerator3D
{
public:
  MicroEdgesGenerator3D(const INTERP_KERNEL::CellModel& cm):_cm(cm) { }
  unsigned getNumberOfSons2(const int *conn, int lgth) const { return _cm.getNumberOfMicroEdges(); }
  unsigned fillSonCellNodalConnectivity2(int sonId, const int *nodalConn, int lgth, int *sonNodalConn, INTERP_KERNEL::NormalizedCellType& typeOfSon) const { return _cm.fillMicroEdgeNodalConnectivity(sonId,nodalConn,sonNodalConn,typeOfSon); }
  static const int DELTA=2;
private:
  const INTERP_KERNEL::CellModel& _cm;
};

int MEDCouplingFastNbrer(int id, unsigned nb, const INTERP_KERNEL::CellModel& cm, bool compute, const int *conn1, const int *conn2);
int MEDCouplingOrientationSensitiveNbrer(int id, unsigned nb, const INTERP_KERNEL::CellModel& cm, bool compute, const int *conn1, const int *conn2);

namespace MEDCoupling
{
  template<const int SPACEDIMM>
  class DummyClsMCUG
  {
  public:
    static const int MY_SPACEDIM=SPACEDIMM;
    static const int MY_MESHDIM=8;
    typedef int MyConnType;
    static const INTERP_KERNEL::NumberingPolicy My_numPol=INTERP_KERNEL::ALL_C_MODE;
    // begin
    // useless, but for windows compilation ...
    const double* getCoordinatesPtr() const { return 0; }
    const int* getConnectivityPtr() const { return 0; }
    const int* getConnectivityIndexPtr() const { return 0; }
    INTERP_KERNEL::NormalizedCellType getTypeOfElement(int) const { return (INTERP_KERNEL::NormalizedCellType)0; }
    // end
  };
}

template<int SPACEDIM>
void MEDCouplingUMesh::getCellsContainingPointsAlg(const double *coords, const double *pos, int nbOfPoints,
                                                   double eps, MCAuto<DataArrayInt>& elts, MCAuto<DataArrayInt>& eltsIndex) const
{
  // Override precision for this method only:
  INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision arcPrec(eps);

  elts=DataArrayInt::New(); eltsIndex=DataArrayInt::New(); eltsIndex->alloc(nbOfPoints+1,1); eltsIndex->setIJ(0,0,0); elts->alloc(0,1);
  int *eltsIndexPtr(eltsIndex->getPointer());
  MCAuto<DataArrayDouble> bboxArr(getBoundingBoxForBBTree(eps));
  const double *bbox(bboxArr->begin());
  int nbOfCells=getNumberOfCells();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  double bb[2*SPACEDIM];
  BBTree<SPACEDIM,int> myTree(&bbox[0],0,0,nbOfCells,-eps);
  for(int i=0;i<nbOfPoints;i++)
    {
      eltsIndexPtr[i+1]=eltsIndexPtr[i];
      for(int j=0;j<SPACEDIM;j++)
        {
          bb[2*j]=pos[SPACEDIM*i+j];
          bb[2*j+1]=pos[SPACEDIM*i+j];
        }
      std::vector<int> candidates;
      myTree.getIntersectingElems(bb,candidates);
      for(std::vector<int>::const_iterator iter=candidates.begin();iter!=candidates.end();iter++)
        {
          int sz(connI[(*iter)+1]-connI[*iter]-1);
          INTERP_KERNEL::NormalizedCellType ct((INTERP_KERNEL::NormalizedCellType)conn[connI[*iter]]);
          bool status(false);
          if(ct!=INTERP_KERNEL::NORM_POLYGON && ct!=INTERP_KERNEL::NORM_QPOLYG)
            status=INTERP_KERNEL::PointLocatorAlgos<DummyClsMCUG<SPACEDIM> >::isElementContainsPoint(pos+i*SPACEDIM,ct,coords,conn+connI[*iter]+1,sz,eps);
          else
            {
              if(SPACEDIM!=2)
                throw INTERP_KERNEL::Exception("MEDCouplingUMesh::getCellsContainingPointsAlg : not implemented yet for POLYGON and QPOLYGON in spaceDim 3 !");
              std::vector<INTERP_KERNEL::Node *> nodes(sz);
              INTERP_KERNEL::QuadraticPolygon *pol(0);
              for(int j=0;j<sz;j++)
                {
                  int nodeId(conn[connI[*iter]+1+j]);
                  nodes[j]=new INTERP_KERNEL::Node(coords[nodeId*SPACEDIM],coords[nodeId*SPACEDIM+1]);
                }
              if(!INTERP_KERNEL::CellModel::GetCellModel(ct).isQuadratic())
                pol=INTERP_KERNEL::QuadraticPolygon::BuildLinearPolygon(nodes);
              else
                pol=INTERP_KERNEL::QuadraticPolygon::BuildArcCirclePolygon(nodes);
              INTERP_KERNEL::Node *n(new INTERP_KERNEL::Node(pos[i*SPACEDIM],pos[i*SPACEDIM+1]));
              double a(0.),b(0.),c(0.);
              a=pol->normalizeMe(b,c); n->applySimilarity(b,c,a);
              status=pol->isInOrOut2(n);
              delete pol; n->decrRef();
            }
          if(status)
            {
              eltsIndexPtr[i+1]++;
              elts->pushBackSilent(*iter);
            }
        }
    }
}

/*!
 * \b WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done.
 */
template<class SonsGenerator>
MEDCouplingUMesh *MEDCouplingUMesh::buildDescendingConnectivityGen(DataArrayInt *desc, DataArrayInt *descIndx, DataArrayInt *revDesc, DataArrayInt *revDescIndx, DimM1DescNbrer nbrer) const
{
  if(!desc || !descIndx || !revDesc || !revDescIndx)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildDescendingConnectivityGen : present of a null pointer in input !");
  checkConnectivityFullyDefined();
  int nbOfCells=getNumberOfCells();
  int nbOfNodes=getNumberOfNodes();
  MCAuto<DataArrayInt> revNodalIndx=DataArrayInt::New(); revNodalIndx->alloc(nbOfNodes+1,1); revNodalIndx->fillWithZero();
  int *revNodalIndxPtr=revNodalIndx->getPointer();
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  std::string name="Mesh constituent of "; name+=getName();
  MCAuto<MEDCouplingUMesh> ret=MEDCouplingUMesh::New(name,getMeshDimension()-SonsGenerator::DELTA);
  ret->setCoords(getCoords());
  ret->allocateCells(2*nbOfCells);
  descIndx->alloc(nbOfCells+1,1);
  MCAuto<DataArrayInt> revDesc2(DataArrayInt::New()); revDesc2->reserve(2*nbOfCells);
  int *descIndxPtr=descIndx->getPointer(); *descIndxPtr++=0;
  for(int eltId=0;eltId<nbOfCells;eltId++,descIndxPtr++)
    {
      int pos=connIndex[eltId];
      int posP1=connIndex[eltId+1];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[pos]);
      SonsGenerator sg(cm);
      unsigned nbOfSons=sg.getNumberOfSons2(conn+pos+1,posP1-pos-1);
      INTERP_KERNEL::AutoPtr<int> tmp=new int[posP1-pos];
      for(unsigned i=0;i<nbOfSons;i++)
        {
          INTERP_KERNEL::NormalizedCellType cmsId;
          unsigned nbOfNodesSon=sg.fillSonCellNodalConnectivity2(i,conn+pos+1,posP1-pos-1,tmp,cmsId);
          for(unsigned k=0;k<nbOfNodesSon;k++)
            if(tmp[k]>=0)
              revNodalIndxPtr[tmp[k]+1]++;
          ret->insertNextCell(cmsId,nbOfNodesSon,tmp);
          revDesc2->pushBackSilent(eltId);
        }
      descIndxPtr[0]=descIndxPtr[-1]+(int)nbOfSons;
    }
  int nbOfCellsM1=ret->getNumberOfCells();
  std::transform(revNodalIndxPtr+1,revNodalIndxPtr+nbOfNodes+1,revNodalIndxPtr,revNodalIndxPtr+1,std::plus<int>());
  MCAuto<DataArrayInt> revNodal=DataArrayInt::New(); revNodal->alloc(revNodalIndx->back(),1);
  std::fill(revNodal->getPointer(),revNodal->getPointer()+revNodalIndx->back(),-1);
  int *revNodalPtr=revNodal->getPointer();
  const int *connM1=ret->getNodalConnectivity()->getConstPointer();
  const int *connIndexM1=ret->getNodalConnectivityIndex()->getConstPointer();
  for(int eltId=0;eltId<nbOfCellsM1;eltId++)
    {
      const int *strtNdlConnOfCurCell=connM1+connIndexM1[eltId]+1;
      const int *endNdlConnOfCurCell=connM1+connIndexM1[eltId+1];
      for(const int *iter=strtNdlConnOfCurCell;iter!=endNdlConnOfCurCell;iter++)
        if(*iter>=0)//for polyhedrons
          *std::find_if(revNodalPtr+revNodalIndxPtr[*iter],revNodalPtr+revNodalIndxPtr[*iter+1],std::bind2nd(std::equal_to<int>(),-1))=eltId;
    }
  //
  DataArrayInt *commonCells=0,*commonCellsI=0;
  FindCommonCellsAlg(3,0,ret->getNodalConnectivity(),ret->getNodalConnectivityIndex(),revNodal,revNodalIndx,commonCells,commonCellsI);
  MCAuto<DataArrayInt> commonCellsTmp(commonCells),commonCellsITmp(commonCellsI);
  const int *commonCellsPtr(commonCells->getConstPointer()),*commonCellsIPtr(commonCellsI->getConstPointer());
  int newNbOfCellsM1=-1;
  MCAuto<DataArrayInt> o2nM1=DataArrayInt::ConvertIndexArrayToO2N(nbOfCellsM1,commonCells->begin(),
                                                                                                            commonCellsI->begin(),commonCellsI->end(),newNbOfCellsM1);
  std::vector<bool> isImpacted(nbOfCellsM1,false);
  for(const int *work=commonCellsI->begin();work!=commonCellsI->end()-1;work++)
    for(int work2=work[0];work2!=work[1];work2++)
      isImpacted[commonCellsPtr[work2]]=true;
  const int *o2nM1Ptr=o2nM1->getConstPointer();
  MCAuto<DataArrayInt> n2oM1=o2nM1->invertArrayO2N2N2OBis(newNbOfCellsM1);
  const int *n2oM1Ptr=n2oM1->getConstPointer();
  MCAuto<MEDCouplingUMesh> ret2=static_cast<MEDCouplingUMesh *>(ret->buildPartOfMySelf(n2oM1->begin(),n2oM1->end(),true));
  ret2->copyTinyInfoFrom(this);
  desc->alloc(descIndx->back(),1);
  int *descPtr=desc->getPointer();
  const INTERP_KERNEL::CellModel& cmsDft=INTERP_KERNEL::CellModel::GetCellModel(INTERP_KERNEL::NORM_POINT1);
  for(int i=0;i<nbOfCellsM1;i++,descPtr++)
    {
      if(!isImpacted[i])
        *descPtr=nbrer(o2nM1Ptr[i],0,cmsDft,false,0,0);
      else
        {
          if(i!=n2oM1Ptr[o2nM1Ptr[i]])
            {
              const INTERP_KERNEL::CellModel& cms=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)connM1[connIndexM1[i]]);
              *descPtr=nbrer(o2nM1Ptr[i],connIndexM1[i+1]-connIndexM1[i]-1,cms,true,connM1+connIndexM1[n2oM1Ptr[o2nM1Ptr[i]]]+1,connM1+connIndexM1[i]+1);
            }
          else
            *descPtr=nbrer(o2nM1Ptr[i],0,cmsDft,false,0,0);
        }
    }
  revDesc->reserve(newNbOfCellsM1);
  revDescIndx->alloc(newNbOfCellsM1+1,1);
  int *revDescIndxPtr=revDescIndx->getPointer(); *revDescIndxPtr++=0;
  const int *revDesc2Ptr=revDesc2->getConstPointer();
  for(int i=0;i<newNbOfCellsM1;i++,revDescIndxPtr++)
    {
      int oldCellIdM1=n2oM1Ptr[i];
      if(!isImpacted[oldCellIdM1])
        {
          revDesc->pushBackSilent(revDesc2Ptr[oldCellIdM1]);
          revDescIndxPtr[0]=revDescIndxPtr[-1]+1;
        }
      else
        {
          for(int j=commonCellsIPtr[0];j<commonCellsIPtr[1];j++)
            revDesc->pushBackSilent(revDesc2Ptr[commonCellsPtr[j]]);
          revDescIndxPtr[0]=revDescIndxPtr[-1]+commonCellsIPtr[1]-commonCellsIPtr[0];
          commonCellsIPtr++;
        }
    }
  //
  return ret2.retn();
}


