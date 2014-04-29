// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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

#include "MEDCouplingStructuredMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingUMesh.hxx"

#include <numeric>

using namespace ParaMEDMEM;

MEDCouplingStructuredMesh::MEDCouplingStructuredMesh()
{
}

MEDCouplingStructuredMesh::MEDCouplingStructuredMesh(const MEDCouplingStructuredMesh& other, bool deepCopy):MEDCouplingMesh(other)
{
}

MEDCouplingStructuredMesh::~MEDCouplingStructuredMesh()
{
}

std::size_t MEDCouplingStructuredMesh::getHeapMemorySizeWithoutChildren() const
{
  return MEDCouplingMesh::getHeapMemorySizeWithoutChildren();
}

void MEDCouplingStructuredMesh::copyTinyStringsFrom(const MEDCouplingMesh *other)
{
  MEDCouplingMesh::copyTinyStringsFrom(other);
}

bool MEDCouplingStructuredMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const
{
  return MEDCouplingMesh::isEqualIfNotWhy(other,prec,reason);
}

INTERP_KERNEL::NormalizedCellType MEDCouplingStructuredMesh::getTypeOfCell(int cellId) const
{
  return GetGeoTypeGivenMeshDimension(getMeshDimension());
}

INTERP_KERNEL::NormalizedCellType MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension(int meshDim)
{
  switch(meshDim)
  {
    case 3:
      return INTERP_KERNEL::NORM_HEXA8;
    case 2:
      return INTERP_KERNEL::NORM_QUAD4;
    case 1:
      return INTERP_KERNEL::NORM_SEG2;
    case 0:
      return INTERP_KERNEL::NORM_POINT1;
    default:
      throw INTERP_KERNEL::Exception("Unexpected dimension for MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension !");
  }
}

std::set<INTERP_KERNEL::NormalizedCellType> MEDCouplingStructuredMesh::getAllGeoTypes() const
{
  std::set<INTERP_KERNEL::NormalizedCellType> ret2;
  ret2.insert(getTypeOfCell(0));
  return ret2;
}

int MEDCouplingStructuredMesh::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  int ret=getNumberOfCells();
  if(type==getTypeOfCell(0))
    return ret;
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getTypeOfCell(0));
  std::ostringstream oss; oss << "MEDCouplingStructuredMesh::getNumberOfCellsWithType : no specified type ! Type available is " << cm.getRepr() << " !";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

DataArrayInt *MEDCouplingStructuredMesh::giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  if(getTypeOfCell(0)==type)
    {
      ret->alloc(getNumberOfCells(),1);
      ret->iota(0);
    }
  else
    ret->alloc(0,1);
  return ret.retn();
}

DataArrayInt *MEDCouplingStructuredMesh::computeNbOfNodesPerCell() const
{
  int nbCells=getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbCells,1);
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getTypeOfCell(0));
  ret->fillWithValue((int)cm.getNumberOfNodes());
  return ret.retn();
}

DataArrayInt *MEDCouplingStructuredMesh::computeNbOfFacesPerCell() const
{
  int nbCells=getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbCells,1);
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getTypeOfCell(0));
  ret->fillWithValue((int)cm.getNumberOfSons());
  return ret.retn();
}

/*!
 * This method computes effective number of nodes per cell. That is to say nodes appearing several times in nodal connectivity of a cell,
 * will be counted only once here whereas it will be counted several times in MEDCouplingMesh::computeNbOfNodesPerCell method.
 * Here for structured mesh it returns exactly as MEDCouplingStructuredMesh::computeNbOfNodesPerCell does.
 *
 * \return DataArrayInt * - new object to be deallocated by the caller.
 */
DataArrayInt *MEDCouplingStructuredMesh::computeEffectiveNbOfNodesPerCell() const
{
  return computeNbOfNodesPerCell();
}

void MEDCouplingStructuredMesh::getNodeIdsOfCell(int cellId, std::vector<int>& conn) const
{
  int meshDim=getMeshDimension();
  int tmpCell[3],tmpNode[3];
  getSplitCellValues(tmpCell);
  getSplitNodeValues(tmpNode);
  int tmp2[3];
  GetPosFromId(cellId,meshDim,tmpCell,tmp2);
  switch(meshDim)
  {
    case 1:
      conn.push_back(tmp2[0]); conn.push_back(tmp2[0]+1);
      break;
    case 2:
      conn.push_back(tmp2[1]*tmpNode[1]+tmp2[0]); conn.push_back(tmp2[1]*tmpNode[1]+tmp2[0]+1);
      conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+1); conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]);
      break;
    case 3:
      conn.push_back(tmp2[1]*tmpNode[1]+tmp2[0]+tmp2[2]*tmpNode[2]); conn.push_back(tmp2[1]*tmpNode[1]+tmp2[0]+1+tmp2[2]*tmpNode[2]);
      conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+1+tmp2[2]*tmpNode[2]); conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+tmp2[2]*tmpNode[2]);
      conn.push_back(tmp2[1]*tmpNode[1]+tmp2[0]+(tmp2[2]+1)*tmpNode[2]); conn.push_back(tmp2[1]*tmpNode[1]+tmp2[0]+1+(tmp2[2]+1)*tmpNode[2]);
      conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+1+(tmp2[2]+1)*tmpNode[2]); conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+(tmp2[2]+1)*tmpNode[2]);
      break;
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::getNodeIdsOfCell : big problem spacedim must be in 1,2 or 3 !");
  };
}

/*!
 * This method returns the mesh dimension of \a this. It can be different from space dimension in case of a not null dimension contains only one node.
 */
int MEDCouplingStructuredMesh::getMeshDimension() const
{
  std::vector<int> ngs(getNodeGridStructure());
  int ret(0),pos(0);
  for(std::vector<int>::const_iterator it=ngs.begin();it!=ngs.end();it++,pos++)
    {
      if(*it<=0)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::getMeshDimension : At pos #" << pos << " number of nodes is " << *it << " ! Must be > 0 !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(*it>1)
        ret++;
    }
  return ret;
}

/*!
 * This method returns the space dimension by only considering the node grid structure.
 * For cartesian mesh the returned value is equal to those returned by getSpaceDimension.
 * But for curvelinear is could be different !
 */
int MEDCouplingStructuredMesh::getSpaceDimensionOnNodeStruct() const
{
  std::vector<int> nodeStr(getNodeGridStructure());
  int spd1(0),pos(0);
  for(std::vector<int>::const_iterator it=nodeStr.begin();it!=nodeStr.end();it++,pos++)
    {
      int elt(*it);
      if(elt<=0)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::getSpaceDimensionOnNodeStruct : At pos #" << pos << " value of node grid structure is " << *it << " ! must be >=1 !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      spd1++;
    }
  return spd1;
}

/*!
 * This method returns the number of cells of unstructured sub level mesh, without building it.
 */
int MEDCouplingStructuredMesh::getNumberOfCellsOfSubLevelMesh() const
{
  std::vector<int> cgs(getCellGridStructure());
  return GetNumberOfCellsOfSubLevelMesh(cgs,getMeshDimension());
}

/*!
 * See MEDCouplingUMesh::getDistributionOfTypes for more information
 */
std::vector<int> MEDCouplingStructuredMesh::getDistributionOfTypes() const
{
  //only one type of cell
  std::vector<int> ret(3);
  ret[0]=getTypeOfCell(0);
  ret[1]=getNumberOfCells();
  ret[2]=-1; //ret[3*k+2]==-1 because it has no sense here
  return ret;
}

/*!
 * This method tries to minimize at most the number of deep copy.
 * So if \a idsPerType is not empty it can be returned directly (without copy, but with ref count incremented) in return.
 * 
 * See MEDCouplingUMesh::checkTypeConsistencyAndContig for more information
 */
DataArrayInt *MEDCouplingStructuredMesh::checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const
{
  int nbOfCells=getNumberOfCells();
  if(code.size()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : invalid input code should be exactly of size 3 !");
  if(code[0]!=(int)getTypeOfCell(0))
    {
      std::ostringstream oss; oss << "MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : Mismatch of geometric type ! Asking for " << code[0] << " whereas the geometric type is \a this is " << getTypeOfCell(0) << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(code[2]==-1)
    {
      if(code[1]==nbOfCells)
        return 0;
      else
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : mismatch between the number of cells in this (" << nbOfCells << ") and the number of non profile (" << code[1] << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(code[2]!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : single geo type mesh ! 0 or -1 is expected at pos #2 of input code !");
  if(idsPerType.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : input code points to DataArrayInt #0 whereas the size of idsPerType is not equal to 1 !");
  const DataArrayInt *pfl=idsPerType[0];
  if(!pfl)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : the input code points to a NULL DataArrayInt at rank 0 !");
  if(pfl->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : input profile should have exactly one component !");
  pfl->checkAllIdsInRange(0,nbOfCells);
  pfl->incrRef();
  return const_cast<DataArrayInt *>(pfl);
}

/*!
 * This method is the opposite of MEDCouplingUMesh::checkTypeConsistencyAndContig method. Given a list of cells in \a profile it returns a list of sub-profiles sorted by geo type.
 * The result is put in the array \a idsPerType. In the returned parameter \a code, foreach i \a code[3*i+2] refers (if different from -1) to a location into the \a idsPerType.
 * This method has 1 input \a profile and 3 outputs \a code \a idsInPflPerType and \a idsPerType.
 * 
 * \param [out] code is a vector of size 3*n where n is the number of different geometric type in \a this \b reduced to the profile \a profile. \a code has exactly the same semantic than in MEDCouplingUMesh::checkTypeConsistencyAndContig method.
 * \param [out] idsInPflPerType is a vector of size of different geometric type in the subpart defined by \a profile of \a this ( equal to \a code.size()/3). For each i,
 *              \a idsInPflPerType[i] stores the tuple ids in \a profile that correspond to the geometric type code[3*i+0]
 * \param [out] idsPerType is a vector of size of different sub profiles needed to be defined to represent the profile \a profile for a given geometric type.
 *              This vector can be empty in case of all geometric type cells are fully covered in ascending in the given input \a profile.
 * 
 * \warning for performance reasons no deep copy will be performed, if \a profile can been used as this in output parameters \a idsInPflPerType and \a idsPerType.
 *
 * \throw if \a profile has not exactly one component. It throws too, if \a profile contains some values not in [0,getNumberOfCells()) or if \a this is not fully defined
 *
 *  \b Example1: <br>
 *          - Before \a this has 3 cells \a profile contains [0,1,2]
 *          - After \a code contains [NORM_...,nbCells,-1], \a idsInPflPerType [[0,1,2]] and \a idsPerType is empty <br>
 * 
 *  \b Example2: <br>
 *          - Before \a this has 3 cells \a profile contains [1,2]
 *          - After \a code contains [NORM_...,nbCells,0], \a idsInPflPerType [[0,1]] and \a idsPerType is [[1,2]] <br>

 */
void MEDCouplingStructuredMesh::splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const
{
  if(!profile || !profile->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::splitProfilePerType : input profile is NULL or not allocated !");
  if(profile->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::splitProfilePerType : input profile should have exactly one component !");
  int nbTuples=profile->getNumberOfTuples();
  int nbOfCells=getNumberOfCells();
  code.resize(3); idsInPflPerType.resize(1);
  code[0]=(int)getTypeOfCell(0); code[1]=nbOfCells;
  idsInPflPerType.resize(1);
  if(profile->isIdentity() && nbTuples==nbOfCells)
    {
      code[2]=-1;
      idsInPflPerType[0]=0;
      idsPerType.clear();
      return ;
    }
  code[1]=profile->getNumberOfTuples();
  code[2]=0;
  profile->checkAllIdsInRange(0,nbOfCells);
  idsPerType.resize(1);
  idsPerType[0]=profile->deepCpy();
  idsInPflPerType[0]=DataArrayInt::Range(0,nbTuples,1);
}

/*!
 * Creates a new unstructured mesh (MEDCoupling1SGTUMesh) from \a this structured one.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is to
 * delete this array using decrRef() as it is no more needed. 
 *  \throw If \a this->getMeshDimension() is not among [1,2,3].
 */
MEDCoupling1SGTUMesh *MEDCouplingStructuredMesh::build1SGTUnstructured() const
{
  int meshDim(getMeshDimension()),spaceDim(getSpaceDimensionOnNodeStruct());
  if((meshDim<0 || meshDim>3) || (spaceDim<0 || spaceDim>3))
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::build1SGTUnstructured : meshdim and spacedim must be in [1,2,3] !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coords(getCoordinatesAndOwner());
  int ns[3];
  getNodeGridStructure(ns);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn(Build1GTNodalConnectivity(ns,ns+spaceDim));
  MEDCouplingAutoRefCountObjectPtr<MEDCoupling1SGTUMesh> ret(MEDCoupling1SGTUMesh::New(getName(),GetGeoTypeGivenMeshDimension(meshDim)));
  ret->setNodalConnectivity(conn); ret->setCoords(coords);
  return ret.retn();
}

/*!
 * This method returns the unstructured mesh (having single geometric type) of the sub level mesh of \a this.
 * This method is equivalent to computing MEDCouplingUMesh::buildDescendingConnectivity on the unstructurized \a this mesh.
 * 
 * The caller is to delete the returned mesh using decrRef() as it is no more needed. 
 */
MEDCoupling1SGTUMesh *MEDCouplingStructuredMesh::build1SGTSubLevelMesh() const
{
  int meshDim(getMeshDimension());
  if(meshDim<1 || meshDim>3)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::build1SGTSubLevelMesh : meshdim must be in [2,3] !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coords(getCoordinatesAndOwner());
  int ns[3];
  getNodeGridStructure(ns);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn(Build1GTNodalConnectivityOfSubLevelMesh(ns,ns+meshDim));
  MEDCouplingAutoRefCountObjectPtr<MEDCoupling1SGTUMesh> ret(MEDCoupling1SGTUMesh::New(getName(),GetGeoTypeGivenMeshDimension(meshDim-1)));
  ret->setNodalConnectivity(conn); ret->setCoords(coords);
  return ret.retn();
}

/*!
 * Creates a new unstructured mesh (MEDCouplingUMesh) from \a this structured one.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is to
 * delete this array using decrRef() as it is no more needed. 
 *  \throw If \a this->getMeshDimension() is not among [1,2,3].
 */
MEDCouplingUMesh *MEDCouplingStructuredMesh::buildUnstructured() const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCoupling1SGTUMesh> ret0(build1SGTUnstructured());
  return ret0->buildUnstructured();
}

/*!
 * Creates a new MEDCouplingUMesh containing a part of cells of \a this mesh.
 * The cells to include to the
 * result mesh are specified by an array of cell ids.
 *  \param [in] start - an array of cell ids to include to the result mesh.
 *  \param [in] end - specifies the end of the array \a start, so that
 *              the last value of \a start is \a end[ -1 ].
 *  \return MEDCouplingMesh * - a new instance of MEDCouplingUMesh. The caller is to
 *         delete this mesh using decrRef() as it is no more needed. 
 */
MEDCouplingMesh *MEDCouplingStructuredMesh::buildPart(const int *start, const int *end) const
{
  MEDCouplingUMesh *um=buildUnstructured();
  MEDCouplingMesh *ret=um->buildPart(start,end);
  um->decrRef();
  return ret;
}

MEDCouplingMesh *MEDCouplingStructuredMesh::buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const
{
  std::vector<int> cgs(getCellGridStructure());
  std::vector< std::pair<int,int> > cellPartFormat,nodePartFormat;
  if(IsPartStructured(start,end,cgs,cellPartFormat))
    {
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingStructuredMesh> ret(buildStructuredSubPart(cellPartFormat));
      nodePartFormat=cellPartFormat;
      for(std::vector< std::pair<int,int> >::iterator it=nodePartFormat.begin();it!=nodePartFormat.end();it++)
        (*it).second++;
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp1(BuildExplicitIdsFrom(getNodeGridStructure(),nodePartFormat));
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp2(DataArrayInt::New()); tmp2->alloc(getNumberOfNodes(),1);
      tmp2->fillWithValue(-1);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp3(DataArrayInt::New()); tmp3->alloc(tmp1->getNumberOfTuples(),1); tmp3->iota(0);
      tmp2->setPartOfValues3(tmp3,tmp1->begin(),tmp1->end(),0,1,1);
      arr=tmp2.retn();
      return ret.retn();
    }
  else
    {
      MEDCouplingUMesh *um=buildUnstructured();
      MEDCouplingMesh *ret=um->buildPartAndReduceNodes(start,end,arr);
      um->decrRef();
      return ret;
    }
}

DataArrayInt *MEDCouplingStructuredMesh::simplexize(int policy)
{
  throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::simplexize : not available for Cartesian mesh !");
}

/*!
 * Returns a new MEDCouplingFieldDouble holding normal vectors to cells of \a this
 * 2D mesh. The computed vectors have 3 components and are normalized.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on
 *          cells and one time. The caller is to delete this field using decrRef() as
 *          it is no more needed.
 *  \throw If \a this->getMeshDimension() != 2.
 */
MEDCouplingFieldDouble *MEDCouplingStructuredMesh::buildOrthogonalField() const
{
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("Expected a MEDCouplingStructuredMesh with meshDim == 2 !");
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  DataArrayDouble *array=DataArrayDouble::New();
  int nbOfCells=getNumberOfCells();
  array->alloc(nbOfCells,3);
  double *vals=array->getPointer();
  for(int i=0;i<nbOfCells;i++)
    { vals[3*i]=0.; vals[3*i+1]=0.; vals[3*i+2]=1.; }
  ret->setArray(array);
  array->decrRef();
  ret->setMesh(this);
  return ret;
}

void MEDCouplingStructuredMesh::getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const
{
  std::vector<int> ngs(getNodeGridStructure());
  int dim(getSpaceDimension());
  switch(dim)
  {
    case 1:
      return GetReverseNodalConnectivity1(ngs,revNodal,revNodalIndx);
    case 2:
      return GetReverseNodalConnectivity2(ngs,revNodal,revNodalIndx);
    case 3:
      return GetReverseNodalConnectivity3(ngs,revNodal,revNodalIndx);
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::getReverseNodalConnectivity : only dimensions 1, 2 and 3 are supported !");
  }
}

void MEDCouplingStructuredMesh::GetReverseNodalConnectivity1(const std::vector<int>& ngs, DataArrayInt *revNodal, DataArrayInt *revNodalIndx)
{
  int nbNodes(ngs[0]);
  revNodalIndx->alloc(nbNodes+1,1);
  if(nbNodes==0)
    { revNodal->alloc(0,1); revNodalIndx->setIJ(0,0,0); return ; }
  if(nbNodes==1)
    { revNodal->alloc(1,1); revNodal->setIJ(0,0,0); revNodalIndx->setIJ(0,0,0); revNodalIndx->setIJ(1,0,1); return ; }
  revNodal->alloc(2*(nbNodes-1),1);
  int *rn(revNodal->getPointer()),*rni(revNodalIndx->getPointer());
  *rni++=0; *rni=1; *rn++=0;
  for(int i=1;i<nbNodes-1;i++,rni++)
    {
      rn[0]=i-1; rn[1]=i;
      rni[1]=rni[0]+2;
      rn+=2;
    }
  rn[0]=nbNodes-2; rni[1]=rni[0]+1;
}

void MEDCouplingStructuredMesh::GetReverseNodalConnectivity2(const std::vector<int>& ngs, DataArrayInt *revNodal, DataArrayInt *revNodalIndx)
{
  int nbNodesX(ngs[0]),nbNodesY(ngs[1]);
  int nbNodes(nbNodesX*nbNodesY);
  if(nbNodesX==0 || nbNodesY==0)
    { revNodal->alloc(0,1); revNodalIndx->setIJ(0,0,0); return ; }
  if(nbNodesX==1 || nbNodesY==1)
    { std::vector<int> ngs2(1); ngs2[0]=std::max(nbNodesX,nbNodesY); return GetReverseNodalConnectivity1(ngs2,revNodal,revNodalIndx); }
  revNodalIndx->alloc(nbNodes+1,1);
  int nbCellsX(nbNodesX-1),nbCellsY(nbNodesY-1);
  revNodal->alloc(4*(nbNodesX-2)*(nbNodesY-2)+2*2*(nbNodesX-2)+2*2*(nbNodesY-2)+4,1);
  int *rn(revNodal->getPointer()),*rni(revNodalIndx->getPointer());
  *rni++=0; *rni=1; *rn++=0;
  for(int i=1;i<nbNodesX-1;i++,rni++,rn+=2)
    {
      rn[0]=i-1; rn[1]=i;
      rni[1]=rni[0]+2;
    }
  rni[1]=rni[0]+1; *rn++=nbCellsX-1;
  rni++;
  for(int j=1;j<nbNodesY-1;j++)
    {
      int off(nbCellsX*(j-1)),off2(nbCellsX*j);
      rni[1]=rni[0]+2; rn[0]=off; rn[1]=off2;
      rni++; rn+=2;
      for(int i=1;i<nbNodesX-1;i++,rni++,rn+=4)
        {
          rn[0]=i-1+off; rn[1]=i+off; rn[2]=i-1+off2; rn[3]=i+off2;
          rni[1]=rni[0]+4;
        }
      rni[1]=rni[0]+2; rn[0]=off+nbCellsX-1; rn[1]=off2+nbCellsX-1;
      rni++; rn+=2;
    }
  int off3(nbCellsX*(nbCellsY-1));
  rni[1]=rni[0]+1;
  rni++; *rn++=off3;
  for(int i=1;i<nbNodesX-1;i++,rni++,rn+=2)
    {
      rn[0]=i-1+off3; rn[1]=i+off3;
      rni[1]=rni[0]+2;
    }
  rni[1]=rni[0]+1; rn[0]=nbCellsX*nbCellsY-1;
}

void MEDCouplingStructuredMesh::GetReverseNodalConnectivity3(const std::vector<int>& ngs, DataArrayInt *revNodal, DataArrayInt *revNodalIndx)
{
  int nbNodesX(ngs[0]),nbNodesY(ngs[1]),nbNodesZ(ngs[2]);
  int nbNodes(nbNodesX*nbNodesY*nbNodesZ);
  if(nbNodesX==0 || nbNodesY==0 || nbNodesZ==0)
    { revNodal->alloc(0,1); revNodalIndx->setIJ(0,0,0); return ; }
  if(nbNodesX==1 || nbNodesY==1 || nbNodesZ==1)
    {
      std::vector<int> ngs2(2);
      int pos(0);
      bool pass(false);
      for(int i=0;i<3;i++)
        {
          if(pass)
            { ngs2[pos++]=ngs[i]; }
          else
            {
              pass=ngs[i]==1;
              if(!pass)
                { ngs2[pos++]=ngs[i]; }
            }
        }
      return GetReverseNodalConnectivity2(ngs2,revNodal,revNodalIndx);
    }
  revNodalIndx->alloc(nbNodes+1,1);
  int nbCellsX(nbNodesX-1),nbCellsY(nbNodesY-1),nbCellsZ(nbNodesZ-1);
  revNodal->alloc(8*(nbNodesX-2)*(nbNodesY-2)*(nbNodesZ-2)+4*(2*(nbNodesX-2)*(nbNodesY-2)+2*(nbNodesX-2)*(nbNodesZ-2)+2*(nbNodesY-2)*(nbNodesZ-2))+2*4*(nbNodesX-2)+2*4*(nbNodesY-2)+2*4*(nbNodesZ-2)+8,1);
  int *rn(revNodal->getPointer()),*rni(revNodalIndx->getPointer());
  *rni=0;
  for(int k=0;k<nbNodesZ;k++)
    {
      bool factZ(k!=0 && k!=nbNodesZ-1);
      int offZ0((k-1)*nbCellsX*nbCellsY),offZ1(k*nbCellsX*nbCellsY);
      for(int j=0;j<nbNodesY;j++)
        {
          bool factYZ(factZ && (j!=0 && j!=nbNodesY-1));
          int off00((j-1)*nbCellsX+offZ0),off01(j*nbCellsX+offZ0),off10((j-1)*nbCellsX+offZ1),off11(j*nbCellsX+offZ1);
          for(int i=0;i<nbNodesX;i++,rni++)
            {
              int fact(factYZ && (i!=0 && i!=nbNodesX-1));
              if(fact)
                {//most of points fall in this part of code
                  rn[0]=off00+i-1; rn[1]=off00+i; rn[2]=off01+i-1; rn[3]=off01+i;
                  rn[4]=off10+i-1; rn[5]=off10+i; rn[6]=off11+i-1; rn[7]=off11+i;
                  rni[1]=rni[0]+8;
                  rn+=8;
                }
              else
                {
                  int *rnRef(rn);
                  if(k>=1 && j>=1 && i>=1)
                    *rn++=off00+i-1;
                  if(k>=1 && j>=1 && i<nbCellsX)
                    *rn++=off00+i;
                  if(k>=1 && j<nbCellsY && i>=1)
                    *rn++=off01+i-1;
                  if(k>=1 && j<nbCellsY && i<nbCellsX)
                    *rn++=off01+i;
                  //
                  if(k<nbCellsZ && j>=1 && i>=1)
                    *rn++=off10+i-1;
                  if(k<nbCellsZ && j>=1 && i<nbCellsX)
                    *rn++=off10+i;
                  if(k<nbCellsZ && j<nbCellsY && i>=1)
                    *rn++=off11+i-1;
                  if(k<nbCellsZ && j<nbCellsY && i<nbCellsX)
                    *rn++=off11+i;
                  rni[1]=rni[0]+(int)(std::distance(rnRef,rn));
                }
            }
        }
    }
}

/*!
 * \return DataArrayInt * - newly allocated instance of nodal connectivity compatible for MEDCoupling1SGTMesh instance
 */
DataArrayInt *MEDCouplingStructuredMesh::Build1GTNodalConnectivity(const int *nodeStBg, const int *nodeStEnd)
{
  int zippedNodeSt[3];
  int dim(ZipNodeStructure(nodeStBg,nodeStEnd,zippedNodeSt));
  switch(dim)
  {
    case 0:
      {
        MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn(DataArrayInt::New());
        conn->alloc(1,1); conn->setIJ(0,0,0);
        return conn.retn();
      }
    case 1:
      return Build1GTNodalConnectivity1D(zippedNodeSt);
    case 2:
      return Build1GTNodalConnectivity2D(zippedNodeSt);
    case 3:
      return Build1GTNodalConnectivity3D(zippedNodeSt);
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::Build1GTNodalConnectivity : only dimension in [0,1,2,3] supported !");
  }
}

DataArrayInt *MEDCouplingStructuredMesh::Build1GTNodalConnectivityOfSubLevelMesh(const int *nodeStBg, const int *nodeStEnd)
{
  std::size_t dim(std::distance(nodeStBg,nodeStEnd));
  switch(dim)
  {
    case 3:
      return Build1GTNodalConnectivityOfSubLevelMesh3D(nodeStBg);
    case 2:
      return Build1GTNodalConnectivityOfSubLevelMesh2D(nodeStBg);
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::Build1GTNodalConnectivityOfSubLevelMesh: only dimension in [2,3] supported !");
  }
}

/*!
 * This method retrieves the number of entities (it can be cells or nodes) given a range in compact standard format
 * used in methods like BuildExplicitIdsFrom,IsPartStructured.
 *
 * \sa BuildExplicitIdsFrom,IsPartStructured
 */
int MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt(const std::vector< std::pair<int,int> >& partCompactFormat)
{
  int ret(1);
  bool isFetched(false);
  std::size_t ii(0);
  for(std::vector< std::pair<int,int> >::const_iterator it=partCompactFormat.begin();it!=partCompactFormat.end();it++,ii++)
    {
      int a((*it).first),b((*it).second);
      if(a<0 || b<0 || b-a<0)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt : invalid input at dimension " << ii << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(b-a>0)
        {
          isFetched=true;
          ret*=(b-a);
        }
    }
  return isFetched?ret:0;
}

DataArrayInt *MEDCouplingStructuredMesh::Build1GTNodalConnectivity1D(const int *nodeStBg)
{
  int nbOfCells(*nodeStBg-1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn(DataArrayInt::New());
  conn->alloc(2*nbOfCells,1);
  int *cp=conn->getPointer();
  for(int i=0;i<nbOfCells;i++)
    {
      cp[2*i+0]=i;
      cp[2*i+1]=i+1;
    }
  return conn.retn();
}

DataArrayInt *MEDCouplingStructuredMesh::Build1GTNodalConnectivity2D(const int *nodeStBg)
{
  int n1=nodeStBg[0]-1;
  int n2=nodeStBg[1]-1;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn(DataArrayInt::New());
  conn->alloc(4*n1*n2,1);
  int *cp=conn->getPointer();
  int pos=0;
  for(int j=0;j<n2;j++)
    for(int i=0;i<n1;i++,pos++)
      {
        cp[4*pos+0]=i+1+j*(n1+1);
        cp[4*pos+1]=i+j*(n1+1);
        cp[4*pos+2]=i+(j+1)*(n1+1);
        cp[4*pos+3]=i+1+(j+1)*(n1+1);
      }
  return conn.retn();
}

DataArrayInt *MEDCouplingStructuredMesh::Build1GTNodalConnectivity3D(const int *nodeStBg)
{
  int n1=nodeStBg[0]-1;
  int n2=nodeStBg[1]-1;
  int n3=nodeStBg[2]-1;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn(DataArrayInt::New());
  conn->alloc(8*n1*n2*n3,1);
  int *cp=conn->getPointer();
  int pos=0;
  for(int k=0;k<n3;k++)
    for(int j=0;j<n2;j++)
      for(int i=0;i<n1;i++,pos++)
        {
          int tmp=(n1+1)*(n2+1);
          cp[8*pos+0]=i+1+j*(n1+1)+k*tmp;
          cp[8*pos+1]=i+j*(n1+1)+k*tmp;
          cp[8*pos+2]=i+(j+1)*(n1+1)+k*tmp;
          cp[8*pos+3]=i+1+(j+1)*(n1+1)+k*tmp;
          cp[8*pos+4]=i+1+j*(n1+1)+(k+1)*tmp;
          cp[8*pos+5]=i+j*(n1+1)+(k+1)*tmp;
          cp[8*pos+6]=i+(j+1)*(n1+1)+(k+1)*tmp;
          cp[8*pos+7]=i+1+(j+1)*(n1+1)+(k+1)*tmp;
        }
  return conn.retn();
}

DataArrayInt *MEDCouplingStructuredMesh::Build1GTNodalConnectivityOfSubLevelMesh3D(const int *nodeStBg)
{
  std::vector<int> ngs(3);
  int n0(nodeStBg[0]-1),n1(nodeStBg[1]-1),n2(nodeStBg[2]-1); ngs[0]=n0; ngs[1]=n1; ngs[2]=n2;
  int off0(nodeStBg[0]),off1(nodeStBg[0]*nodeStBg[1]);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn(DataArrayInt::New());
  conn->alloc(4*GetNumberOfCellsOfSubLevelMesh(ngs,3));
  int *cp(conn->getPointer());
  //X
  for(int i=0;i<nodeStBg[0];i++)
    for(int j=0;j<n1;j++)
      for(int k=0;k<n2;k++,cp+=4)
        { cp[0]=k*off1+j*off0+i; cp[1]=(k+1)*off1+j*off0+i; cp[2]=(k+1)*off1+(j+1)*off0+i; cp[3]=k*off1+(j+1)*off0+i; }
  //Y
  for(int j=0;j<nodeStBg[1];j++)
    for(int i=0;i<n0;i++)
      for(int k=0;k<n2;k++,cp+=4)
        { cp[0]=k*off1+j*off0+i; cp[1]=(k+1)*off1+j*off0+i; cp[2]=(k+1)*off1+j*off0+(i+1); cp[3]=k*off1+j*off0+(i+1); }
  //Z
  for(int k=0;k<nodeStBg[2];k++)
    for(int i=0;i<n0;i++)
      for(int j=0;j<n1;j++,cp+=4)
        { cp[0]=k*off1+j*off0+i; cp[1]=k*off1+j*off0+(i+1); cp[2]=k*off1+(j+1)*off0+(i+1); cp[3]=k*off1+(j+1)*off0+i; }
  return conn.retn();
}

/*!
 * This method computes given the nodal structure defined by [ \a nodeStBg , \a nodeStEnd ) the zipped form.
 * std::distance( \a nodeStBg, \a nodeStEnd ) is equal to the space dimension. The returned value is equal to
 * the meshDimension (or the zipped spaceDimension).
 *
 * \param [out] zipNodeSt - The zipped node strucutre
 * \return int - the
 */
int MEDCouplingStructuredMesh::ZipNodeStructure(const int *nodeStBg, const int *nodeStEnd, int zipNodeSt[3])
{
  int spaceDim((int)std::distance(nodeStBg,nodeStEnd));
  if(spaceDim>3 || spaceDim<1)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ZipNodeStructure : spaceDim must in [1,2,3] !");
  zipNodeSt[0]=0; zipNodeSt[1]=0; zipNodeSt[2]=0;
  int zippedI(0);
  for(int i=0;i<spaceDim;i++)
    {
      int elt(nodeStBg[i]);
      if(elt<1)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::ZipNodeStructure : the input nodal structure at pos#" << i << "(" << nodeStBg[i] << ") is invalid !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(elt>=2)
        zipNodeSt[zippedI++]=elt;
    }
  return zippedI;
}

DataArrayInt *MEDCouplingStructuredMesh::Build1GTNodalConnectivityOfSubLevelMesh2D(const int *nodeStBg)
{
  std::vector<int> ngs(2);
  int n0(nodeStBg[0]-1),n1(nodeStBg[1]-1); ngs[0]=n0; ngs[1]=n1;
  int off0(nodeStBg[0]);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn(DataArrayInt::New());
  conn->alloc(2*GetNumberOfCellsOfSubLevelMesh(ngs,2));
  int *cp(conn->getPointer());
  //X
  for(int i=0;i<nodeStBg[0];i++)
    for(int j=0;j<n1;j++,cp+=2)
      { cp[0]=j*off0+i; cp[1]=(j+1)*off0+i; }
  //Y
  for(int j=0;j<nodeStBg[1];j++)
    for(int i=0;i<n0;i++,cp+=2)
      { cp[0]=j*off0+i; cp[1]=j*off0+(i+1); }
  return conn.retn();
}

/*!
 * Returns a cell id by its (i,j,k) index. The cell is located between the i-th and
 * ( i + 1 )-th nodes along X axis etc.
 *  \param [in] i - a index of node coordinates array along X axis.
 *  \param [in] j - a index of node coordinates array along Y axis.
 *  \param [in] k - a index of node coordinates array along Z axis.
 *  \return int - a cell id in \a this mesh.
 */
int MEDCouplingStructuredMesh::getCellIdFromPos(int i, int j, int k) const
{
  int tmp[3]={i,j,k};
  int tmp2[3];
  int meshDim(getMeshDimension());
  getSplitCellValues(tmp2);
  std::transform(tmp,tmp+meshDim,tmp2,tmp,std::multiplies<int>());
  return std::accumulate(tmp,tmp+meshDim,0);
}

/*!
 * Returns a node id by its (i,j,k) index.
 *  \param [in] i - a index of node coordinates array along X axis.
 *  \param [in] j - a index of node coordinates array along Y axis.
 *  \param [in] k - a index of node coordinates array along Z axis.
 *  \return int - a node id in \a this mesh.
 */
int MEDCouplingStructuredMesh::getNodeIdFromPos(int i, int j, int k) const
{
  int tmp[3]={i,j,k};
  int tmp2[3];
  int spaceDim(getSpaceDimension());
  getSplitNodeValues(tmp2);
  std::transform(tmp,tmp+spaceDim,tmp2,tmp,std::multiplies<int>());
  return std::accumulate(tmp,tmp+spaceDim,0);
}


int MEDCouplingStructuredMesh::getNumberOfCells() const
{
  std::vector<int> ngs(getNodeGridStructure());
  int ret(1);
  bool isCatched(false);
  std::size_t ii(0);
  for(std::vector<int>::const_iterator it=ngs.begin();it!=ngs.end();it++,ii++)
    {
      int elt(*it);
      if(elt<=0)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::getNumberOfCells : at pos #" << ii << " the number of nodes in nodeStructure is " << *it << " ! Must be > 0 !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(elt>1)
        {
          ret*=elt-1;
          isCatched=true;
        }
    }
  return isCatched?ret:0;
}

int MEDCouplingStructuredMesh::getNumberOfNodes() const
{
  std::vector<int> ngs(getNodeGridStructure());
  int ret(1);
  for(std::vector<int>::const_iterator it=ngs.begin();it!=ngs.end();it++)
    ret*=*it;
  return ret;
}

void MEDCouplingStructuredMesh::GetPosFromId(int nodeId, int meshDim, const int *split, int *res)
{
  int work=nodeId;
  for(int i=meshDim-1;i>=0;i--)
    {
      int pos=work/split[i];
      work=work%split[i];
      res[i]=pos;
    }
}

std::vector<int> MEDCouplingStructuredMesh::getCellGridStructure() const
{
  std::vector<int> ret(getNodeGridStructure());
  std::transform(ret.begin(),ret.end(),ret.begin(),std::bind2nd(std::plus<int>(),-1));
  return ret;
}

/*!
 * This method states if given part ids [ \a startIds, \a stopIds) and a structure \a st returns if it can be considered as a structured dataset.
 * If true is returned \a partCompactFormat will contain the information to build the corresponding part.
 *
 * \sa MEDCouplingStructuredMesh::BuildExplicitIdsFrom, MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt
 */
bool MEDCouplingStructuredMesh::IsPartStructured(const int *startIds, const int *stopIds, const std::vector<int>& st, std::vector< std::pair<int,int> >& partCompactFormat)
{
  int dim((int)st.size());
  partCompactFormat.resize(dim);
  if(dim<1 || dim>3)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::isPartStructured : input structure must be of dimension in [1,2,3] !");
  std::vector<int> tmp2(dim),tmp(dim),tmp3(dim),tmp4(dim); tmp2[0]=1;
  for(int i=1;i<dim;i++)
    tmp2[i]=tmp2[i-1]*st[i-1];
  std::size_t sz(std::distance(startIds,stopIds));
  if(sz==0)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::IsPartStructured : empty input !");
  GetPosFromId(*startIds,dim,&tmp2[0],&tmp[0]);
  partCompactFormat.resize(dim);
  for(int i=0;i<dim;i++)
    partCompactFormat[i].first=tmp[i];
  if(tmp[dim-1]<0 || tmp[dim-1]>=st[dim-1])
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::IsPartStructured : first id in input is not in valid range !");
  if(sz==1)
    {
      for(int i=0;i<dim;i++)
        partCompactFormat[i].second=tmp[i]+1;
      return true;
    }
  GetPosFromId(startIds[sz-1],dim,&tmp2[0],&tmp3[0]);
  int szExp(1);
  for(int i=0;i<dim;i++)
    {
      if(tmp3[i]<0 || tmp3[i]>=st[i])
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::IsPartStructured : last id in input is not in valid range !");
      partCompactFormat[i].second=tmp3[i]+1;
      tmp4[i]=partCompactFormat[i].second-partCompactFormat[i].first;
      if(tmp4[i]<=0)
        return false;
      szExp*=tmp4[i];
    }
  if(szExp!=(int)sz)
    return false;
  const int *w(startIds);
  switch(dim)
  {
    case 3:
      {
        for(int i=0;i<tmp4[2];i++)
          {
            int a=tmp2[2]*(partCompactFormat[2].first+i);
            for(int j=0;j<tmp4[1];j++)
              {
                int b=tmp2[1]*(partCompactFormat[1].first+j);
                for(int k=0;k<tmp4[0];k++,w++)
                  {
                    if(partCompactFormat[0].first+k+b+a!=*w)
                      return false;
                  }
              }
          }
        return true;
      }
    case 2:
      {
        for(int j=0;j<tmp4[1];j++)
          {
            int b=tmp2[1]*(partCompactFormat[1].first+j);
            for(int k=0;k<tmp4[0];k++,w++)
              {
                if(partCompactFormat[0].first+k+b!=*w)
                  return false;
              }
          }
        return true;
      }
    case 1:
      {
        for(int k=0;k<tmp4[0];k++,w++)
          {
            if(partCompactFormat[0].first+k!=*w)
              return false;
          }
        return true;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::IsPartStructured : internal error !");
  }
}

/*!
 * This method builds the explicit entity array from the structure in \a st and the range in \a partCompactFormat.
 * If the range contains invalid values regarding sructure an exception will be thrown.
 *
 * \return DataArrayInt * - a new object.
 * \sa MEDCouplingStructuredMesh::IsPartStructured, MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt
 */
DataArrayInt *MEDCouplingStructuredMesh::BuildExplicitIdsFrom(const std::vector<int>& st, const std::vector< std::pair<int,int> >& partCompactFormat)
{
  if(st.size()!=partCompactFormat.size())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : input arrays must have the same size !");
  int nbOfItems(1);
  std::vector<int> dims(st.size());
  for(std::size_t i=0;i<st.size();i++)
    {
      if(partCompactFormat[i].first<0 || partCompactFormat[i].first>st[i])
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : invalid input range 1 !");
      if(partCompactFormat[i].second<0 || partCompactFormat[i].second>st[i])
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : invalid input range 2 !");
      if(partCompactFormat[i].second<=partCompactFormat[i].first)
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : invalid input range 3 !");
      dims[i]=partCompactFormat[i].second-partCompactFormat[i].first;
      nbOfItems*=dims[i];
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New());
  ret->alloc(nbOfItems,1);
  int *pt(ret->getPointer());
  switch(st.size())
  {
    case 3:
      {
        for(int i=0;i<dims[2];i++)
          {
            int a=(partCompactFormat[2].first+i)*st[0]*st[1];
            for(int j=0;j<dims[1];j++)
              {
                int b=(partCompactFormat[1].first+j)*st[0];
                for(int k=0;k<dims[0];k++,pt++)
                  *pt=partCompactFormat[0].first+k+b+a;
              }
          }
        break;
      }
    case 2:
      {
        for(int j=0;j<dims[1];j++)
          {
            int b=(partCompactFormat[1].first+j)*st[0];
            for(int k=0;k<dims[0];k++,pt++)
              *pt=partCompactFormat[0].first+k+b;
          }
        break;
      }
    case 1:
      {
        for(int k=0;k<dims[0];k++,pt++)
          *pt=partCompactFormat[0].first+k;
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : Dimension supported are 1,2 or 3 !");
  }
  return ret.retn();
}

int MEDCouplingStructuredMesh::GetNumberOfCellsOfSubLevelMesh(const std::vector<int>& cgs, int mdim)
{
  int ret(0);
  for(int i=0;i<mdim;i++)
    {
      int locRet(1);
      for(int j=0;j<mdim;j++)
        if(j!=i)
          locRet*=cgs[j];
        else
          locRet*=cgs[j]+1;
      ret+=locRet;
    }
  return ret;
}
