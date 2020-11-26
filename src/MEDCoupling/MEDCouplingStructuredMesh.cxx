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

#include "MEDCouplingStructuredMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingIMesh.hxx"//tony to throw when optimization will be performed in AssignPartOfFieldOfDoubleUsing

#include <numeric>

using namespace MEDCoupling;

MEDCouplingStructuredMesh::MEDCouplingStructuredMesh()
{
}

MEDCouplingStructuredMesh::MEDCouplingStructuredMesh(const MEDCouplingStructuredMesh& other, bool deepCpy):MEDCouplingMesh(other)
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

INTERP_KERNEL::NormalizedCellType MEDCouplingStructuredMesh::getTypeOfCell(mcIdType cellId) const
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

mcIdType MEDCouplingStructuredMesh::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  mcIdType ret(getNumberOfCells());
  if(type==getTypeOfCell(0))
    return ret;
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(getTypeOfCell(0)));
  std::ostringstream oss; oss << "MEDCouplingStructuredMesh::getNumberOfCellsWithType : no specified type ! Type available is " << cm.getRepr() << " !";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

DataArrayIdType *MEDCouplingStructuredMesh::giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  if(getTypeOfCell(0)==type)
    {
      ret->alloc(getNumberOfCells(),1);
      ret->iota(0);
    }
  else
    ret->alloc(0,1);
  return ret.retn();
}

DataArrayIdType *MEDCouplingStructuredMesh::computeNbOfNodesPerCell() const
{
  std::size_t nbCells=getNumberOfCells();
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  ret->alloc(nbCells,1);
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getTypeOfCell(0));
  ret->fillWithValue(ToIdType(cm.getNumberOfNodes()));
  return ret.retn();
}

DataArrayIdType *MEDCouplingStructuredMesh::computeNbOfFacesPerCell() const
{
  std::size_t nbCells=getNumberOfCells();
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  ret->alloc(nbCells,1);
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getTypeOfCell(0));
  ret->fillWithValue(ToIdType(cm.getNumberOfSons()));
  return ret.retn();
}

/*!
 * This method computes effective number of nodes per cell. That is to say nodes appearing several times in nodal connectivity of a cell,
 * will be counted only once here whereas it will be counted several times in MEDCouplingMesh::computeNbOfNodesPerCell method.
 * Here for structured mesh it returns exactly as MEDCouplingStructuredMesh::computeNbOfNodesPerCell does.
 *
 * \return DataArrayIdType * - new object to be deallocated by the caller.
 */
DataArrayIdType *MEDCouplingStructuredMesh::computeEffectiveNbOfNodesPerCell() const
{
  return computeNbOfNodesPerCell();
}

void MEDCouplingStructuredMesh::getNodeIdsOfCell(mcIdType cellId, std::vector<mcIdType>& conn) const
{
  int meshDim=getMeshDimension();
  mcIdType tmpCell[3],tmpNode[3];
  getSplitCellValues(tmpCell);
  getSplitNodeValues(tmpNode);
  mcIdType tmp2[3];
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
  std::vector<mcIdType> ngs(getNodeGridStructure());
  int ret(0),pos(0);
  for(std::vector<mcIdType>::const_iterator it=ngs.begin();it!=ngs.end();it++,pos++)
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
  std::vector<mcIdType> nodeStr(getNodeGridStructure());
  int spd1(0),pos(0);
  for(std::vector<mcIdType>::const_iterator it=nodeStr.begin();it!=nodeStr.end();it++,pos++)
    {
      mcIdType elt(*it);
      if(elt<=0)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::getSpaceDimensionOnNodeStruct : At pos #" << pos << " value of node grid structure is " << *it << " ! must be >=1 !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      spd1++;
    }
  return spd1;
}

void MEDCouplingStructuredMesh::getSplitCellValues(mcIdType *res) const
{
  std::vector<mcIdType> strct(getCellGridStructure());
  std::vector<mcIdType> ret(MEDCouplingStructuredMesh::GetSplitVectFromStruct(strct));
  std::copy(ret.begin(),ret.end(),res);
}

void MEDCouplingStructuredMesh::getSplitNodeValues(mcIdType *res) const
{
  std::vector<mcIdType> strct(getNodeGridStructure());
  std::vector<mcIdType> ret(MEDCouplingStructuredMesh::GetSplitVectFromStruct(strct));
  std::copy(ret.begin(),ret.end(),res);
}

/*!
 * This method returns the number of cells of unstructured sub level mesh, without building it.
 */
mcIdType MEDCouplingStructuredMesh::getNumberOfCellsOfSubLevelMesh() const
{
  std::vector<mcIdType> cgs(getCellGridStructure());
  return GetNumberOfCellsOfSubLevelMesh(cgs,getMeshDimension());
}

/*!
 * See MEDCouplingUMesh::getDistributionOfTypes for more information
 */
std::vector<mcIdType> MEDCouplingStructuredMesh::getDistributionOfTypes() const
{
  //only one type of cell
  std::vector<mcIdType> ret(3);
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
DataArrayIdType *MEDCouplingStructuredMesh::checkTypeConsistencyAndContig(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const
{
  mcIdType nbOfCells=getNumberOfCells();
  if(code.size()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : invalid input code should be exactly of size 3 !");
  if(code[0]!=ToIdType(getTypeOfCell(0)))
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
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : input code points to DataArrayIdType #0 whereas the size of idsPerType is not equal to 1 !");
  const DataArrayIdType *pfl=idsPerType[0];
  if(!pfl)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : the input code points to a NULL DataArrayIdType at rank 0 !");
  if(pfl->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : input profile should have exactly one component !");
  pfl->checkAllIdsInRange(0,nbOfCells);
  pfl->incrRef();
  return const_cast<DataArrayIdType *>(pfl);
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
void MEDCouplingStructuredMesh::splitProfilePerType(const DataArrayIdType *profile, std::vector<mcIdType>& code, std::vector<DataArrayIdType *>& idsInPflPerType, std::vector<DataArrayIdType *>& idsPerType, bool smartPflKiller) const
{
  if(!profile || !profile->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::splitProfilePerType : input profile is NULL or not allocated !");
  if(profile->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::splitProfilePerType : input profile should have exactly one component !");
  mcIdType nbTuples(profile->getNumberOfTuples());
  mcIdType nbOfCells=getNumberOfCells();
  code.resize(3); idsInPflPerType.resize(1);
  code[0]=ToIdType(getTypeOfCell(0)); code[1]=nbOfCells;
  idsInPflPerType.resize(1);
  if(smartPflKiller && profile->isIota(nbOfCells))
    {
      code[2]=-1;
      idsInPflPerType[0]=profile->deepCopy();
      idsPerType.clear();
      return ;
    }
  code[1]=profile->getNumberOfTuples();
  code[2]=0;
  profile->checkAllIdsInRange(0,nbOfCells);
  idsPerType.resize(1);
  idsPerType[0]=profile->deepCopy();
  idsInPflPerType[0]=DataArrayIdType::Range(0,nbTuples,1);
}

/*!
 * Creates a new unstructured mesh (MEDCoupling1SGTUMesh) from \a this structured one.
 *
 * In the returned mesh, the nodes are ordered with the first axis varying first: (X0,Y0), (X1,Y0),  ... (X0,Y1), (X1,Y1), ...
 * and the cells are ordered with the same logic, i.e. in (i,j) notation: (0,0), (1,0), (2,0), ... (0,1), (1,1), ...
 *
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is to
 * delete this array using decrRef() as it is no more needed.
 *  \throw If \a this->getMeshDimension() is not among [1,2,3].
 */
MEDCoupling1SGTUMesh *MEDCouplingStructuredMesh::build1SGTUnstructured() const
{
  int meshDim(getMeshDimension()),spaceDim(getSpaceDimensionOnNodeStruct());
  if((meshDim<0 || meshDim>3) || (spaceDim<0 || spaceDim>3))
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::build1SGTUnstructured : meshdim and spacedim must be in [1,2,3] !");
  MCAuto<DataArrayDouble> coords(getCoordinatesAndOwner());
  mcIdType ns[3];
  getNodeGridStructure(ns);
  MCAuto<DataArrayIdType> conn(Build1GTNodalConnectivity(ns,ns+spaceDim));
  MCAuto<MEDCoupling1SGTUMesh> ret(MEDCoupling1SGTUMesh::New(getName(),GetGeoTypeGivenMeshDimension(meshDim)));
  ret->setNodalConnectivity(conn); ret->setCoords(coords);
  try
    { ret->copyTinyInfoFrom(this); }
  catch(INTERP_KERNEL::Exception&) { }
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
  MCAuto<DataArrayDouble> coords(getCoordinatesAndOwner());
  mcIdType ns[3];
  getNodeGridStructure(ns);
  MCAuto<DataArrayIdType> conn(Build1GTNodalConnectivityOfSubLevelMesh(ns,ns+meshDim));
  MCAuto<MEDCoupling1SGTUMesh> ret(MEDCoupling1SGTUMesh::New(getName(),GetGeoTypeGivenMeshDimension(meshDim-1)));
  ret->setNodalConnectivity(conn); ret->setCoords(coords);
  return ret.retn();
}

/*!
 * Creates a new unstructured mesh (MEDCouplingUMesh) from \a this structured one.
 *
 * In the returned mesh, the nodes are ordered with the first axis varying first: (X0,Y0), (X1,Y0),  ... (X0,Y1), (X1,Y1), ...
 * and the cells are ordered with the same logic, i.e. in (i,j) notation: (0,0), (1,0), (2,0), ... (0,1), (1,1), ...
 *
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is to
 * delete this array using decrRef() as it is no more needed.
 *  \throw If \a this->getMeshDimension() is not among [1,2,3].
 */
MEDCouplingUMesh *MEDCouplingStructuredMesh::buildUnstructured() const
{
  MCAuto<MEDCoupling1SGTUMesh> ret0(build1SGTUnstructured());
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
MEDCouplingMesh *MEDCouplingStructuredMesh::buildPart(const mcIdType *start, const mcIdType *end) const
{
  MCAuto<MEDCouplingUMesh> um(buildUnstructured());
  return um->buildPart(start,end);
}

MEDCouplingMesh *MEDCouplingStructuredMesh::buildPartAndReduceNodes(const mcIdType *start, const mcIdType *end, DataArrayIdType*& arr) const
{
  std::vector<mcIdType> cgs(getCellGridStructure());
  std::vector< std::pair<mcIdType,mcIdType> > cellPartFormat,nodePartFormat;
  if(IsPartStructured(start,end,cgs,cellPartFormat))
    {
      MCAuto<MEDCouplingStructuredMesh> ret(buildStructuredSubPart(cellPartFormat));
      nodePartFormat=cellPartFormat;
      for(std::vector< std::pair<mcIdType,mcIdType> >::iterator it=nodePartFormat.begin();it!=nodePartFormat.end();it++)
        (*it).second++;
      MCAuto<DataArrayIdType> tmp1(BuildExplicitIdsFrom(getNodeGridStructure(),nodePartFormat));
      MCAuto<DataArrayIdType> tmp2(DataArrayIdType::New()); tmp2->alloc(getNumberOfNodes(),1);
      tmp2->fillWithValue(-1);
      MCAuto<DataArrayIdType> tmp3(DataArrayIdType::New()); tmp3->alloc(tmp1->getNumberOfTuples(),1); tmp3->iota(0);
      tmp2->setPartOfValues3(tmp3,tmp1->begin(),tmp1->end(),0,1,1);
      arr=tmp2.retn();
      return ret.retn();
    }
  else
    {
      MCAuto<MEDCouplingUMesh> um(buildUnstructured());
      return um->buildPartAndReduceNodes(start,end,arr);
    }
}

DataArrayIdType *MEDCouplingStructuredMesh::simplexize(int policy)
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
  MCAuto<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME));
  MCAuto<DataArrayDouble> array(DataArrayDouble::New());
  mcIdType nbOfCells=getNumberOfCells();
  array->alloc(nbOfCells,3);
  double *vals(array->getPointer());
  for(mcIdType i=0;i<nbOfCells;i++)
    { vals[3*i]=0.; vals[3*i+1]=0.; vals[3*i+2]=1.; }
  ret->setArray(array);
  ret->setMesh(this);
  return ret.retn();
}

void MEDCouplingStructuredMesh::getReverseNodalConnectivity(DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx) const
{
  std::vector<mcIdType> ngs(getNodeGridStructure());
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

void MEDCouplingStructuredMesh::GetReverseNodalConnectivity1(const std::vector<mcIdType>& ngs, DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx)
{
  mcIdType nbNodes(ngs[0]);
  revNodalIndx->alloc(nbNodes+1,1);
  if(nbNodes==0)
    { revNodal->alloc(0,1); revNodalIndx->setIJ(0,0,0); return ; }
  if(nbNodes==1)
    { revNodal->alloc(1,1); revNodal->setIJ(0,0,0); revNodalIndx->setIJ(0,0,0); revNodalIndx->setIJ(1,0,1); return ; }
  revNodal->alloc(2*(nbNodes-1),1);
  mcIdType *rn(revNodal->getPointer()),*rni(revNodalIndx->getPointer());
  *rni++=0; *rni=1; *rn++=0;
  for(mcIdType i=1;i<nbNodes-1;i++,rni++)
    {
      rn[0]=i-1; rn[1]=i;
      rni[1]=rni[0]+2;
      rn+=2;
    }
  rn[0]=nbNodes-2; rni[1]=rni[0]+1;
}

void MEDCouplingStructuredMesh::GetReverseNodalConnectivity2(const std::vector<mcIdType>& ngs, DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx)
{
  mcIdType nbNodesX(ngs[0]),nbNodesY(ngs[1]);
  mcIdType nbNodes(nbNodesX*nbNodesY);
  if(nbNodesX==0 || nbNodesY==0)
    { revNodal->alloc(0,1); revNodalIndx->setIJ(0,0,0); return ; }
  if(nbNodesX==1 || nbNodesY==1)
    { std::vector<mcIdType> ngs2(1); ngs2[0]=std::max(nbNodesX,nbNodesY); return GetReverseNodalConnectivity1(ngs2,revNodal,revNodalIndx); }
  revNodalIndx->alloc(nbNodes+1,1);
  mcIdType nbCellsX(nbNodesX-1),nbCellsY(nbNodesY-1);
  revNodal->alloc(4*(nbNodesX-2)*(nbNodesY-2)+2*2*(nbNodesX-2)+2*2*(nbNodesY-2)+4,1);
  mcIdType *rn(revNodal->getPointer()),*rni(revNodalIndx->getPointer());
  *rni++=0; *rni=1; *rn++=0;
  for(mcIdType i=1;i<nbNodesX-1;i++,rni++,rn+=2)
    {
      rn[0]=i-1; rn[1]=i;
      rni[1]=rni[0]+2;
    }
  rni[1]=rni[0]+1; *rn++=nbCellsX-1;
  rni++;
  for(mcIdType j=1;j<nbNodesY-1;j++)
    {
      mcIdType off(nbCellsX*(j-1)),off2(nbCellsX*j);
      rni[1]=rni[0]+2; rn[0]=off; rn[1]=off2;
      rni++; rn+=2;
      for(mcIdType i=1;i<nbNodesX-1;i++,rni++,rn+=4)
        {
          rn[0]=i-1+off; rn[1]=i+off; rn[2]=i-1+off2; rn[3]=i+off2;
          rni[1]=rni[0]+4;
        }
      rni[1]=rni[0]+2; rn[0]=off+nbCellsX-1; rn[1]=off2+nbCellsX-1;
      rni++; rn+=2;
    }
  mcIdType off3(nbCellsX*(nbCellsY-1));
  rni[1]=rni[0]+1;
  rni++; *rn++=off3;
  for(mcIdType i=1;i<nbNodesX-1;i++,rni++,rn+=2)
    {
      rn[0]=i-1+off3; rn[1]=i+off3;
      rni[1]=rni[0]+2;
    }
  rni[1]=rni[0]+1; rn[0]=nbCellsX*nbCellsY-1;
}

void MEDCouplingStructuredMesh::GetReverseNodalConnectivity3(const std::vector<mcIdType>& ngs, DataArrayIdType *revNodal, DataArrayIdType *revNodalIndx)
{
  mcIdType nbNodesX(ngs[0]),nbNodesY(ngs[1]),nbNodesZ(ngs[2]);
  mcIdType nbNodes(nbNodesX*nbNodesY*nbNodesZ);
  if(nbNodesX==0 || nbNodesY==0 || nbNodesZ==0)
    { revNodal->alloc(0,1); revNodalIndx->setIJ(0,0,0); return ; }
  if(nbNodesX==1 || nbNodesY==1 || nbNodesZ==1)
    {
      std::vector<mcIdType> ngs2(2);
      mcIdType pos(0);
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
  mcIdType nbCellsX(nbNodesX-1),nbCellsY(nbNodesY-1),nbCellsZ(nbNodesZ-1);
  revNodal->alloc(8*(nbNodesX-2)*(nbNodesY-2)*(nbNodesZ-2)+4*(2*(nbNodesX-2)*(nbNodesY-2)+2*(nbNodesX-2)*(nbNodesZ-2)+2*(nbNodesY-2)*(nbNodesZ-2))+2*4*(nbNodesX-2)+2*4*(nbNodesY-2)+2*4*(nbNodesZ-2)+8,1);
  mcIdType *rn(revNodal->getPointer()),*rni(revNodalIndx->getPointer());
  *rni=0;
  for(mcIdType k=0;k<nbNodesZ;k++)
    {
      bool factZ(k!=0 && k!=nbNodesZ-1);
      mcIdType offZ0((k-1)*nbCellsX*nbCellsY),offZ1(k*nbCellsX*nbCellsY);
      for(mcIdType j=0;j<nbNodesY;j++)
        {
          bool factYZ(factZ && (j!=0 && j!=nbNodesY-1));
          mcIdType off00((j-1)*nbCellsX+offZ0),off01(j*nbCellsX+offZ0),off10((j-1)*nbCellsX+offZ1),off11(j*nbCellsX+offZ1);
          for(mcIdType i=0;i<nbNodesX;i++,rni++)
            {
              mcIdType fact(factYZ && (i!=0 && i!=nbNodesX-1));
              if(fact)
                {//most of points fall in this part of code
                  rn[0]=off00+i-1; rn[1]=off00+i; rn[2]=off01+i-1; rn[3]=off01+i;
                  rn[4]=off10+i-1; rn[5]=off10+i; rn[6]=off11+i-1; rn[7]=off11+i;
                  rni[1]=rni[0]+8;
                  rn+=8;
                }
              else
                {
                  mcIdType *rnRef(rn);
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
                  rni[1]=rni[0]+ToIdType(std::distance(rnRef,rn));
                }
            }
        }
    }
}

/*!
 * \return DataArrayIdType * - newly allocated instance of nodal connectivity compatible for MEDCoupling1SGTMesh instance
 */
DataArrayIdType *MEDCouplingStructuredMesh::Build1GTNodalConnectivity(const mcIdType *nodeStBg, const mcIdType *nodeStEnd)
{
  mcIdType zippedNodeSt[3];
  mcIdType dim(ZipNodeStructure(nodeStBg,nodeStEnd,zippedNodeSt));
  switch(dim)
  {
    case 0:
      {
        MCAuto<DataArrayIdType> conn(DataArrayIdType::New());
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

DataArrayIdType *MEDCouplingStructuredMesh::Build1GTNodalConnectivityOfSubLevelMesh(const mcIdType *nodeStBg, const mcIdType *nodeStEnd)
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
 * This method returns the list of ids sorted ascendingly of entities that are in the corner in ghost zone.
 * The ids are returned in a newly created DataArrayIdType having a single component.
 *
 * \param [in] st - The structure \b without ghost cells.
 * \param [in] ghostLev - The size of the ghost zone (>=0)
 * \return DataArrayIdType * - The DataArray containing all the ids the caller is to deallocate.
 */
DataArrayIdType *MEDCouplingStructuredMesh::ComputeCornersGhost(const std::vector<mcIdType>& st, mcIdType ghostLev)
{
  if(ghostLev<0)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ComputeCornersGhost : ghost lev must be >= 0 !");
  std::size_t dim(st.size());
  MCAuto<DataArrayIdType> ret(DataArrayIdType::New());
  switch(dim)
  {
    case 1:
      {
        ret->alloc(2*ghostLev,1);
        mcIdType *ptr(ret->getPointer());
        for(mcIdType i=0;i<ghostLev;i++,ptr++)
          *ptr=i;
        mcIdType offset(st[0]);
        if(offset<0)
          throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ComputeCornersGhost : element in 1D structure must be >= 0 !");
        for(mcIdType i=0;i<ghostLev;i++,ptr++)
          *ptr=offset+ghostLev+i;
        break;
      }
    case 2:
      {
        mcIdType offsetX(st[0]),offsetY(st[1]);
        if(offsetX<0 || offsetY<0)
          throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ComputeCornersGhost : elements in 2D structure must be >= 0 !");
        ret->alloc(4*ghostLev,1);
        mcIdType *ptr(ret->getPointer());
        for(mcIdType i=0;i<ghostLev;i++)
          {
            *ptr++=i*(2*ghostLev+offsetX+1);
            *ptr++=offsetX+2*ghostLev-1+i*(2*ghostLev+offsetX-1);
          }
        for(mcIdType i=0;i<ghostLev;i++)
          {
            *ptr++=(2*ghostLev+offsetX)*(offsetY+ghostLev)+ghostLev-1+i*(2*ghostLev+offsetX-1);
            *ptr++=(2*ghostLev+offsetX)*(offsetY+ghostLev)+offsetX+ghostLev+i*(2*ghostLev+offsetX+1);
          }
        break;
      }
    case 3:
      {
        mcIdType offsetX(st[0]),offsetY(st[1]),offsetZ(st[2]);
        if(offsetX<0 || offsetY<0 || offsetZ<0)
          throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ComputeCornersGhost : elements in 3D structure must be >= 0 !");
        ret->alloc(8*ghostLev,1);
        mcIdType *ptr(ret->getPointer());
        mcIdType zeOffsetZ((offsetX+2*ghostLev)*(offsetY+2*ghostLev));
        for(mcIdType i=0;i<ghostLev;i++)
          {
            *ptr++=i*(2*ghostLev+offsetX+1)+i*zeOffsetZ;
            *ptr++=offsetX+2*ghostLev-1+i*(2*ghostLev+offsetX-1)+i*zeOffsetZ;
            *ptr++=(2*ghostLev+offsetX)*(offsetY+ghostLev)+ghostLev-1+(ghostLev-i-1)*(2*ghostLev+offsetX-1)+i*zeOffsetZ;
            *ptr++=(2*ghostLev+offsetX)*(offsetY+ghostLev)+offsetX+ghostLev+(ghostLev-i-1)*(2*ghostLev+offsetX+1)+i*zeOffsetZ;
          }
        mcIdType j(0),zeOffsetZ2(zeOffsetZ*(offsetZ+ghostLev));
        for(mcIdType i=ghostLev-1;i>=0;i--,j++)
          {
            *ptr++=i*(2*ghostLev+offsetX+1)+j*zeOffsetZ+zeOffsetZ2;
            *ptr++=offsetX+2*ghostLev-1+i*(2*ghostLev+offsetX-1)+j*zeOffsetZ+zeOffsetZ2;
            *ptr++=(2*ghostLev+offsetX)*(offsetY+ghostLev)+ghostLev-1+(ghostLev-i-1)*(2*ghostLev+offsetX-1)+j*zeOffsetZ+zeOffsetZ2;
            *ptr++=(2*ghostLev+offsetX)*(offsetY+ghostLev)+offsetX+ghostLev+(ghostLev-i-1)*(2*ghostLev+offsetX+1)+j*zeOffsetZ+zeOffsetZ2;
          }
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ComputeCornersGhost : Only dimensions 1, 2 and 3 are supported actually !");
  }
  return ret.retn();
}

/*!
 * This method retrieves the number of entities (it can be cells or nodes) given a range in compact standard format
 * used in methods like BuildExplicitIdsFrom,IsPartStructured.
 *
 * \sa BuildExplicitIdsFrom,IsPartStructured
 */
mcIdType MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt(const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat)
{
  mcIdType ret(1);
  std::size_t ii(0);
  for(std::vector< std::pair<mcIdType,mcIdType> >::const_iterator it=partCompactFormat.begin();it!=partCompactFormat.end();it++,ii++)
    {
      mcIdType a((*it).first),b((*it).second);
      if(a<0 || b<0 || b-a<0)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt : invalid input at dimension " << ii << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      ret*=(b-a);
    }
  return ret;
}

mcIdType MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(const std::vector<mcIdType>& st)
{
  mcIdType ret(1);
  bool isFetched(false);
  for(std::size_t i=0;i<st.size();i++)
    {
      if(st[i]<0)
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure : presence of a negative value in structure !");
      ret*=st[i];
      isFetched=true;
    }
  return isFetched?ret:0;
}

void MEDCouplingStructuredMesh::FindTheWidestAxisOfGivenRangeInCompactFrmt(const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat, int& axisId, mcIdType& sizeOfRange)
{
    mcIdType dim(ToIdType(partCompactFormat.size()));
    mcIdType ret(-1);
    for(int i=0;i<dim;i++)
      {
        mcIdType curDelta(partCompactFormat[i].second-partCompactFormat[i].first);
        if(curDelta<0)
          {
            std::ostringstream oss; oss << "MEDCouplingStructuredMesh::FindTheWidestAxisOfGivenRangeInCompactFrmt : at axis #" << i << " the range is invalid (first value < second value) !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
        if(curDelta>ret)
          {
            axisId=i; sizeOfRange=curDelta;
            ret=curDelta;
          }
      }
}

/*!
 * This method is \b NOT wrapped in python because it has no sense in python (for performance reasons).
 * This method starts from a structured mesh with structure \a st on which a boolean field \a crit is set.
 * This method find for such minimalist information of mesh and field which is the part of the mesh, given by the range per axis in output parameter
 * \a partCompactFormat that contains all the True in \a crit. The returned vector of boolean is the field reduced to that part.
 * So the number of True is equal in \a st and in returned vector of boolean.
 *
 * \param [in] minPatchLgth - minimum length that the patch may have for all directions.
 * \param [in] st - The structure per axis of the structured mesh considered.
 * \param [in] crit - The field of boolean (for performance reasons) lying on the mesh defined by \a st.
 * \param [out] partCompactFormat - The minimal part of \a st containing all the true of \a crit.
 * \param [out] reducedCrit - The reduction of \a criterion on \a partCompactFormat.
 * \return - The number of True in \a st (that is equal to those in \a reducedCrit)
 */
mcIdType MEDCouplingStructuredMesh::FindMinimalPartOf(mcIdType minPatchLgth, const std::vector<mcIdType>& st, const std::vector<bool>& crit, std::vector<bool>& reducedCrit, std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat)
{
  if(minPatchLgth<0)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::FindMinimalPartOf : the input minPatchLgth has to be >=0 !");
  if(ToIdType(crit.size())!=DeduceNumberOfGivenStructure(st))
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::FindMinimalPartOf : size of vector of boolean is invalid regarding the declared structure !");
  mcIdType ret(-1);
  switch(st.size())
  {
    case 1:
      {
        ret=FindMinimalPartOf1D(st,crit,partCompactFormat);
        break;
      }
    case 2:
      {
        ret=FindMinimalPartOf2D(st,crit,partCompactFormat);
        break;
      }
    case 3:
      {
        ret=FindMinimalPartOf3D(st,crit,partCompactFormat);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::FindMinimalPartOf : only dimension 1, 2 and 3 are supported actually !");
  }
  std::vector<mcIdType> dims(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(partCompactFormat));
  mcIdType i(0);
  for(std::vector< std::pair<mcIdType,mcIdType> >::iterator it=partCompactFormat.begin();it!=partCompactFormat.end();it++,i++)
    {
      if(st[i]<minPatchLgth)
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::FindMinimalPartOf : the input patch is tinier than the min length constraint !");
      mcIdType start((*it).first),stop((*it).second),middle((start+stop)/2);
      if(stop-start<minPatchLgth)
        {
          (*it).first=middle-minPatchLgth/2;
          (*it).second=middle+minPatchLgth-minPatchLgth/2;
          if((*it).first<0)
            {
              (*it).second+=-(*it).first;
              (*it).first=0;
            }
          if((*it).second>st[i])
            {
              (*it).first-=(*it).second-st[i];
              (*it).second=st[i];
            }
        }
    }
  ExtractFieldOfBoolFrom(st,crit,partCompactFormat,reducedCrit);
  return ret;
}

/*!
 * This method is \b NOT wrapped in python.
 * This method considers \a crit input parameter as a matrix having dimensions specified by \a st. This method returns for each axis
 * the signature, that is to say the number of elems equal to true in \a crit along this axis.
 */
std::vector< std::vector<mcIdType> > MEDCouplingStructuredMesh::ComputeSignaturePerAxisOf(const std::vector<mcIdType>& st, const std::vector<bool>& crit)
{
  std::size_t dim(st.size());
  std::vector< std::vector<mcIdType> > ret(dim);
  switch(dim)
  {
    case 1:
      {
        mcIdType nx(st[0]);
        ret[0].resize(nx);
        std::vector<mcIdType>& retX(ret[0]);
        for(mcIdType i=0;i<nx;i++)
          retX[i]=crit[i]?1:0;
        break;
      }
    case 2:
      {
        mcIdType nx(st[0]),ny(st[1]);
        ret[0].resize(nx); ret[1].resize(ny);
        std::vector<mcIdType>& retX(ret[0]);
        for(mcIdType i=0;i<nx;i++)
          {
            mcIdType cnt(0);
            for(mcIdType j=0;j<ny;j++)
              if(crit[j*nx+i])
                cnt++;
            retX[i]=cnt;
          }
        std::vector<mcIdType>& retY(ret[1]);
        for(mcIdType j=0;j<ny;j++)
          {
            mcIdType cnt(0);
            for(mcIdType i=0;i<nx;i++)
              if(crit[j*nx+i])
                cnt++;
            retY[j]=cnt;
          }
        break;
      }
    case 3:
      {
        mcIdType nx(st[0]),ny(st[1]),nz(st[2]);
        ret[0].resize(nx); ret[1].resize(ny); ret[2].resize(nz);
        std::vector<mcIdType>& retX(ret[0]);
        for(mcIdType i=0;i<nx;i++)
          {
            mcIdType cnt(0);
            for(mcIdType k=0;k<nz;k++)
              {
                mcIdType offz(k*nx*ny+i);
                for(mcIdType j=0;j<ny;j++)
                  if(crit[offz+j*nx])
                    cnt++;
              }
            retX[i]=cnt;
          }
        std::vector<mcIdType>& retY(ret[1]);
        for(mcIdType j=0;j<ny;j++)
          {
            mcIdType cnt(0),offy(j*nx);
            for(mcIdType k=0;k<nz;k++)
              {
                mcIdType offz(k*nx*ny+offy);
                for(mcIdType i=0;i<nx;i++)
                  if(crit[offz+i])
                    cnt++;
              }
            retY[j]=cnt;
          }
        std::vector<mcIdType>& retZ(ret[2]);
        for(mcIdType k=0;k<nz;k++)
          {
            mcIdType cnt(0),offz(k*nx*ny);
            for(mcIdType j=0;j<ny;j++)
              {
                mcIdType offy(offz+j*nx);
                for(mcIdType i=0;i<nx;i++)
                  if(crit[offy+i])
                    cnt++;
              }
            retZ[k]=cnt;
          }
        break;
      }
    default:
       throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ComputeSignatureOf : only dimensions 1, 2 and 3 are supported !");
  }
  return ret;
}

DataArrayIdType *MEDCouplingStructuredMesh::Build1GTNodalConnectivity1D(const mcIdType *nodeStBg)
{
  mcIdType nbOfCells=*nodeStBg-1;
  MCAuto<DataArrayIdType> conn(DataArrayIdType::New());
  conn->alloc(2*nbOfCells,1);
  mcIdType *cp=conn->getPointer();
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      cp[2*i+0]=i;
      cp[2*i+1]=i+1;
    }
  return conn.retn();
}

DataArrayIdType *MEDCouplingStructuredMesh::Build1GTNodalConnectivity2D(const mcIdType *nodeStBg)
{
  mcIdType n1(nodeStBg[0]-1),n2(nodeStBg[1]-1);
  MCAuto<DataArrayIdType> conn(DataArrayIdType::New());
  conn->alloc(4*n1*n2,1);
  mcIdType *cp(conn->getPointer());
  std::size_t pos(0);
  for(mcIdType j=0;j<n2;j++)
    for(mcIdType i=0;i<n1;i++,pos++)
      {
        cp[4*pos+0]=i+1+j*(n1+1);
        cp[4*pos+1]=i+j*(n1+1);
        cp[4*pos+2]=i+(j+1)*(n1+1);
        cp[4*pos+3]=i+1+(j+1)*(n1+1);
      }
  return conn.retn();
}

DataArrayIdType *MEDCouplingStructuredMesh::Build1GTNodalConnectivity3D(const mcIdType *nodeStBg)
{
  mcIdType n1(nodeStBg[0]-1),n2(nodeStBg[1]-1),n3(nodeStBg[2]-1);
  MCAuto<DataArrayIdType> conn(DataArrayIdType::New());
  conn->alloc(8*n1*n2*n3,1);
  mcIdType *cp(conn->getPointer());
  std::size_t pos(0);
  for(mcIdType k=0;k<n3;k++)
    for(mcIdType j=0;j<n2;j++)
      for(mcIdType i=0;i<n1;i++,pos++)
        {
          mcIdType tmp=(n1+1)*(n2+1);
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

DataArrayIdType *MEDCouplingStructuredMesh::Build1GTNodalConnectivityOfSubLevelMesh3D(const mcIdType *nodeStBg)
{
  std::vector<mcIdType> ngs(3);
  mcIdType n0(nodeStBg[0]-1),n1(nodeStBg[1]-1),n2(nodeStBg[2]-1); ngs[0]=n0; ngs[1]=n1; ngs[2]=n2;
  mcIdType off0(nodeStBg[0]),off1(nodeStBg[0]*nodeStBg[1]);
  MCAuto<DataArrayIdType> conn(DataArrayIdType::New());
  conn->alloc(4*GetNumberOfCellsOfSubLevelMesh(ngs,3));
  mcIdType *cp(conn->getPointer());
  //X
  for(mcIdType i=0;i<nodeStBg[0];i++)
    for(mcIdType j=0;j<n1;j++)
      for(mcIdType k=0;k<n2;k++,cp+=4)
        { cp[0]=k*off1+j*off0+i; cp[1]=(k+1)*off1+j*off0+i; cp[2]=(k+1)*off1+(j+1)*off0+i; cp[3]=k*off1+(j+1)*off0+i; }
  //Y
  for(mcIdType j=0;j<nodeStBg[1];j++)
    for(mcIdType i=0;i<n0;i++)
      for(mcIdType k=0;k<n2;k++,cp+=4)
        { cp[0]=k*off1+j*off0+i; cp[1]=(k+1)*off1+j*off0+i; cp[2]=(k+1)*off1+j*off0+(i+1); cp[3]=k*off1+j*off0+(i+1); }
  //Z
  for(mcIdType k=0;k<nodeStBg[2];k++)
    for(mcIdType i=0;i<n0;i++)
      for(mcIdType j=0;j<n1;j++,cp+=4)
        { cp[0]=k*off1+j*off0+i; cp[1]=k*off1+j*off0+(i+1); cp[2]=k*off1+(j+1)*off0+(i+1); cp[3]=k*off1+(j+1)*off0+i; }
  return conn.retn();
}

/*!
 * \sa MEDCouplingStructuredMesh::FindMinimalPartOf
 */
mcIdType MEDCouplingStructuredMesh::FindMinimalPartOf1D(const std::vector<mcIdType>& st, const std::vector<bool>& crit, std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat)
{
  if(st.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::FindMinimalPartOf1D : the input size of st must be equal to 1 !");
  mcIdType nxMin(std::numeric_limits<mcIdType>::max()),nxMax(-std::numeric_limits<mcIdType>::max());
  mcIdType nx(st[0]),ret(0);
  for(mcIdType i=0;i<nx;i++)
    {
      if(crit[i])
        {
          nxMin=std::min(nxMin,i); nxMax=std::max(nxMax,i);
          ret++;
        }
    }
  if(ret==0)
    {
      std::size_t sz(st.size());
      partCompactFormat.resize(sz);
      for(std::size_t i=0;i<sz;i++)
        {
          partCompactFormat[i].first=st[i]/2;
          partCompactFormat[i].second=st[i]/2;
        }
      return ret;
    }
  partCompactFormat.resize(1);
  partCompactFormat[0].first=nxMin; partCompactFormat[0].second=nxMax+1;
  return ret;
}

/*!
 * \sa MEDCouplingStructuredMesh::FindMinimalPartOf
 */
mcIdType MEDCouplingStructuredMesh::FindMinimalPartOf2D(const std::vector<mcIdType>& st, const std::vector<bool>& crit, std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat)
{
  if(st.size()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::FindMinimalPartOf2D : the input size of st must be equal to 2 !");
  mcIdType nxMin(std::numeric_limits<mcIdType>::max()),nxMax(-std::numeric_limits<mcIdType>::max()),nyMin(std::numeric_limits<mcIdType>::max()),nyMax(-std::numeric_limits<mcIdType>::max());
  mcIdType it(0),nx(st[0]),ny(st[1]);
  mcIdType ret(0);
  for(mcIdType i=0;i<ny;i++)
    for(mcIdType j=0;j<nx;j++,it++)
      {
        if(crit[it])
          {
            nxMin=std::min(nxMin,j); nxMax=std::max(nxMax,j);
            nyMin=std::min(nyMin,i); nyMax=std::max(nyMax,i);
            ret++;
          }
      }
  if(ret==0)
    {
      std::size_t sz(st.size());
      partCompactFormat.resize(sz);
      for(std::size_t i=0;i<sz;i++)
        {
          partCompactFormat[i].first=st[i]/2;
          partCompactFormat[i].second=st[i]/2;
        }
      return ret;
    }
  partCompactFormat.resize(2);
  partCompactFormat[0].first=nxMin; partCompactFormat[0].second=nxMax+1;
  partCompactFormat[1].first=nyMin; partCompactFormat[1].second=nyMax+1;
  return ret;
}

/*!
 * \sa MEDCouplingStructuredMesh::FindMinimalPartOf
 */
mcIdType MEDCouplingStructuredMesh::FindMinimalPartOf3D(const std::vector<mcIdType>& st, const std::vector<bool>& crit, std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat)
{
  if(st.size()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::FindMinimalPartOf3D : the input size of st must be equal to 3 !");
  mcIdType nxMin(std::numeric_limits<mcIdType>::max()),nxMax(-std::numeric_limits<mcIdType>::max()),nyMin(std::numeric_limits<mcIdType>::max()),nyMax(-std::numeric_limits<mcIdType>::max()),nzMin(std::numeric_limits<mcIdType>::max()),nzMax(-std::numeric_limits<mcIdType>::max());
  mcIdType it(0),nx(st[0]),ny(st[1]),nz(st[2]);
  mcIdType ret(0);
  for(mcIdType i=0;i<nz;i++)
    for(mcIdType j=0;j<ny;j++)
      for(mcIdType k=0;k<nx;k++,it++)
        {
          if(crit[it])
            {
              nxMin=std::min(nxMin,k); nxMax=std::max(nxMax,k);
              nyMin=std::min(nyMin,j); nyMax=std::max(nyMax,j);
              nzMin=std::min(nzMin,i); nzMax=std::max(nzMax,i);
              ret++;
            }
        }
  if(ret==0)
    {
      std::size_t sz(st.size());
      partCompactFormat.resize(sz);
      for(std::size_t i=0;i<sz;i++)
        {
          partCompactFormat[i].first=st[i]/2;
          partCompactFormat[i].second=st[i]/2;
        }
      return ret;
    }
  partCompactFormat.resize(3);
  partCompactFormat[0].first=nxMin; partCompactFormat[0].second=nxMax+1;
  partCompactFormat[1].first=nyMin; partCompactFormat[1].second=nyMax+1;
  partCompactFormat[2].first=nzMin; partCompactFormat[2].second=nzMax+1;
  return ret;
}

/*!
 * This method computes given the nodal structure defined by [ \a nodeStBg , \a nodeStEnd ) the zipped form.
 * std::distance( \a nodeStBg, \a nodeStEnd ) is equal to the space dimension. The returned value is equal to
 * the meshDimension (or the zipped spaceDimension).
 *
 * \param [out] zipNodeSt - The zipped node structure
 * \return mcIdType - the
 */
int MEDCouplingStructuredMesh::ZipNodeStructure(const mcIdType *nodeStBg, const mcIdType *nodeStEnd, mcIdType zipNodeSt[3])
{
  std::size_t spaceDim(std::distance(nodeStBg,nodeStEnd));
  if(spaceDim>3 || spaceDim<1)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ZipNodeStructure : spaceDim must in [1,2,3] !");
  zipNodeSt[0]=0; zipNodeSt[1]=0; zipNodeSt[2]=0;
  int zippedI(0);
  for(std::size_t i=0;i<spaceDim;i++)
    {
      mcIdType elt(nodeStBg[i]);
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

DataArrayIdType *MEDCouplingStructuredMesh::Build1GTNodalConnectivityOfSubLevelMesh2D(const mcIdType *nodeStBg)
{
  std::vector<mcIdType> ngs(2);
  mcIdType n0(nodeStBg[0]-1),n1(nodeStBg[1]-1); ngs[0]=n0; ngs[1]=n1;
  mcIdType off0(nodeStBg[0]);
  MCAuto<DataArrayIdType> conn(DataArrayIdType::New());
  conn->alloc(2*GetNumberOfCellsOfSubLevelMesh(ngs,2));
  mcIdType *cp(conn->getPointer());
  //X
  for(mcIdType i=0;i<nodeStBg[0];i++)
    for(mcIdType j=0;j<n1;j++,cp+=2)
      { cp[0]=j*off0+i; cp[1]=(j+1)*off0+i; }
  //Y
  for(mcIdType j=0;j<nodeStBg[1];j++)
    for(mcIdType i=0;i<n0;i++,cp+=2)
      { cp[0]=j*off0+i; cp[1]=j*off0+(i+1); }
  return conn.retn();
}

/*!
 * Returns a cell id by its (i,j,k) index. The cell is located between the i-th and
 * ( i + 1 )-th nodes along X axis etc.
 *  \param [in] i - a index of node coordinates array along X axis.
 *  \param [in] j - a index of node coordinates array along Y axis.
 *  \param [in] k - a index of node coordinates array along Z axis.
 *  \return mcIdType - a cell id in \a this mesh.
 */
mcIdType MEDCouplingStructuredMesh::getCellIdFromPos(mcIdType i, mcIdType j, mcIdType k) const
{
  mcIdType tmp[3]={i,j,k};
  mcIdType tmp2[3];
  mcIdType meshDim(getMeshDimension());
  getSplitCellValues(tmp2);
  std::transform(tmp,tmp+meshDim,tmp2,tmp,std::multiplies<mcIdType>());
  return std::accumulate(tmp,tmp+meshDim,0);
}

/*!
 * Returns a node id by its (i,j,k) index.
 *  \param [in] i - a index of node coordinates array along X axis.
 *  \param [in] j - a index of node coordinates array along Y axis.
 *  \param [in] k - a index of node coordinates array along Z axis.
 *  \return mcIdType - a node id in \a this mesh.
 */
mcIdType MEDCouplingStructuredMesh::getNodeIdFromPos(mcIdType i, mcIdType j, mcIdType k) const
{
  mcIdType tmp[3]={i,j,k};
  mcIdType tmp2[3];
  mcIdType spaceDim(getSpaceDimension());
  getSplitNodeValues(tmp2);
  std::transform(tmp,tmp+spaceDim,tmp2,tmp,std::multiplies<mcIdType>());
  return std::accumulate(tmp,tmp+spaceDim,0);
}


mcIdType MEDCouplingStructuredMesh::getNumberOfCells() const
{
  std::vector<mcIdType> ngs(getNodeGridStructure());
  mcIdType ret(1);
  bool isCatched(false);
  std::size_t ii(0);
  for(std::vector<mcIdType>::const_iterator it=ngs.begin();it!=ngs.end();it++,ii++)
    {
      mcIdType elt(*it);
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

mcIdType MEDCouplingStructuredMesh::getNumberOfNodes() const
{
  std::vector<mcIdType> ngs(getNodeGridStructure());
  mcIdType ret(1);
  for(std::vector<mcIdType>::const_iterator it=ngs.begin();it!=ngs.end();it++)
    ret*=*it;
  return ret;
}

/*!
 * This method returns for a cell which id is \a cellId the location (locX,locY,locZ) of this cell in \a this.
 *
 * \param [in] cellId ID of the cell
 * \return - A vector of size this->getMeshDimension()
 * \throw if \a cellId not in [ 0, this->getNumberOfCells() )
 */
std::vector<mcIdType> MEDCouplingStructuredMesh::getLocationFromCellId(mcIdType cellId) const
{
  int meshDim(getMeshDimension());
  std::vector<mcIdType> ret(meshDim);
  std::vector<mcIdType> struc(getCellGridStructure());
  mcIdType nbCells(std::accumulate(struc.begin(),struc.end(),1,std::multiplies<mcIdType>()));
  if(cellId<0 || cellId>=nbCells)
    {
      std::ostringstream oss; oss << "MEDCouplingStructuredMesh::getLocationFromCellId : Input cell id (" << cellId << ") is invalid ! Should be in [0," << nbCells << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::vector<mcIdType> spt(GetSplitVectFromStruct(struc));
  GetPosFromId(cellId,meshDim,&spt[0],&ret[0]);
  return ret;
}

/*!
 * This method returns for a node which id is \a nodeId the location (locX,locY,locZ) of this node in \a this.
 *
 * \param [in] nodeId ID of the node
 * \return - A vector of size this->getSpaceDimension()
 * \throw if \a cellId not in [ 0, this->getNumberOfNodes() )
 */
std::vector<mcIdType> MEDCouplingStructuredMesh::getLocationFromNodeId(mcIdType nodeId) const
{
  int spaceDim(getSpaceDimension());
  std::vector<mcIdType> ret(spaceDim);
  std::vector<mcIdType> struc(getNodeGridStructure());
  mcIdType nbNodes(std::accumulate(struc.begin(),struc.end(),1,std::multiplies<mcIdType>()));
  if(nodeId<0 || nodeId>=nbNodes)
    {
      std::ostringstream oss; oss << "MEDCouplingStructuredMesh::getLocationFromNodeId : Input node id (" << nodeId << ") is invalid ! Should be in [0," << nbNodes << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::vector<mcIdType> spt(GetSplitVectFromStruct(struc));
  GetPosFromId(nodeId,spaceDim,&spt[0],&ret[0]);
  return ret;
}

void MEDCouplingStructuredMesh::GetPosFromId(mcIdType eltId, int meshDim, const mcIdType *split, mcIdType *res)
{
  mcIdType work(eltId);
  for(int i=meshDim-1;i>=0;i--)
    {
      mcIdType pos=work/split[i];
      work=work%split[i];
      res[i]=pos;
    }
}

std::vector<mcIdType> MEDCouplingStructuredMesh::getCellGridStructure() const
{
  std::vector<mcIdType> ret(getNodeGridStructure());
  std::transform(ret.begin(),ret.end(),ret.begin(),std::bind(std::plus<mcIdType>(),std::placeholders::_1,-1));
  return ret;
}

/*!
 * This method returns the squareness of \a this (quadrature). \a this is expected to be with a mesh dimension equal to 2 or 3.
 */
double MEDCouplingStructuredMesh::computeSquareness() const
{
  std::vector<mcIdType> cgs(getCellGridStructure());
  if(cgs.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::computeSquareness : empty mesh !");
  std::size_t dim(cgs.size());
  if(dim==1)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::computeSquareness : A segment cannot be square !");
  if(dim<4)
    {
      mcIdType minAx(cgs[0]),maxAx(cgs[0]);
      for(std::size_t i=1;i<dim;i++)
        {
          minAx=std::min(minAx,cgs[i]);
          maxAx=std::max(maxAx,cgs[i]);
        }
      return (double)minAx/(double)maxAx;
    }
  throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::computeSquareness : only dimension 2 and 3 supported !");
}

/*!
 * Given a struct \a strct it returns a split vector [1,strct[0],strct[0]*strct[1]...]
 * This decomposition allows to quickly find i,j,k given a global id.
 */
std::vector<mcIdType> MEDCouplingStructuredMesh::GetSplitVectFromStruct(const std::vector<mcIdType>& strct)
{
  std::size_t spaceDim(strct.size());
  std::vector<mcIdType> res(spaceDim);
  for(std::size_t l=0;l<spaceDim;l++)
    {
      mcIdType val=1;
      for(std::size_t p=0;p<spaceDim-l-1;p++)
        val*=strct[p];
      res[spaceDim-l-1]=val;
    }
  return res;
}

/*!
 * This method states if given part ids [ \a startIds, \a stopIds) and a structure \a st returns if it can be considered as a structured dataset.
 * If true is returned \a partCompactFormat will contain the information to build the corresponding part.
 *
 * \sa MEDCouplingStructuredMesh::BuildExplicitIdsFrom, MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt
 */
bool MEDCouplingStructuredMesh::IsPartStructured(const mcIdType *startIds, const mcIdType *stopIds, const std::vector<mcIdType>& st, std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat)
{
  int dim((int)st.size());
  partCompactFormat.resize(dim);
  if(dim<1 || dim>3)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::isPartStructured : input structure must be of dimension in [1,2,3] !");
  std::vector<mcIdType> tmp2(dim),tmp(dim),tmp3(dim),tmp4(dim); tmp2[0]=1;
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
  mcIdType szExp(1);
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
  if(szExp!=ToIdType(sz))
    return false;
  const mcIdType *w(startIds);
  switch(dim)
  {
    case 3:
      {
        for(mcIdType i=0;i<tmp4[2];i++)
          {
            mcIdType a=tmp2[2]*(partCompactFormat[2].first+i);
            for(mcIdType j=0;j<tmp4[1];j++)
              {
                mcIdType b=tmp2[1]*(partCompactFormat[1].first+j);
                for(mcIdType k=0;k<tmp4[0];k++,w++)
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
        for(mcIdType j=0;j<tmp4[1];j++)
          {
            mcIdType b=tmp2[1]*(partCompactFormat[1].first+j);
            for(mcIdType k=0;k<tmp4[0];k++,w++)
              {
                if(partCompactFormat[0].first+k+b!=*w)
                  return false;
              }
          }
        return true;
      }
    case 1:
      {
        for(mcIdType k=0;k<tmp4[0];k++,w++)
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
 * This method takes in input a compact format [[Xmax,Xmin),[Ymin,Ymax)] and returns the corresponding dimensions for each axis that is to say
 * [Xmax-Xmin,Ymax-Ymin].
 *
 * \throw if an axis range is so that max<min
 * \sa GetCompactFrmtFromDimensions
 */
std::vector<mcIdType> MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat)
{
  std::vector<mcIdType> ret(partCompactFormat.size());
  for(std::size_t i=0;i<partCompactFormat.size();i++)
    {
      if(partCompactFormat[i].first>partCompactFormat[i].second)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt : For axis #" << i << " end is before start !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      ret[i]=partCompactFormat[i].second-partCompactFormat[i].first;
    }
  return ret;
}

/*!
 * This method takes in input a vector giving the number of entity per axis and returns for each axis a range starting from [0,0...]
 *
 * \throw if there is an axis in \a dims that is < 0.
 * \sa GetDimensionsFromCompactFrmt, ChangeReferenceFromGlobalOfCompactFrmt, ChangeReferenceToGlobalOfCompactFrmt
 */
std::vector< std::pair<mcIdType,mcIdType> > MEDCouplingStructuredMesh::GetCompactFrmtFromDimensions(const std::vector<mcIdType>& dims)
{
  std::size_t sz(dims.size());
  std::vector< std::pair<mcIdType,mcIdType> > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      if(dims[i]<0)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt : For axis #" << i << " dimension < 0 !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      ret[i].first=0;
      ret[i].second=dims[i];
    }
  return ret;
}

/*!
 * This method returns the intersection zone of two ranges (in compact format) \a r1 and \a r2.
 * This method will throw exception if on one axis the intersection is empty.
 *
 * \sa AreRangesIntersect
 */
std::vector< std::pair<mcIdType,mcIdType> > MEDCouplingStructuredMesh::IntersectRanges(const std::vector< std::pair<mcIdType,mcIdType> >& r1, const std::vector< std::pair<mcIdType,mcIdType> >& r2)
{
  std::size_t sz(r1.size());
  if(sz!=r2.size())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::IntersectRanges : the two ranges must have the same dimension !");
  std::vector< std::pair<mcIdType,mcIdType> > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      if(r1[i].first>r1[i].second)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::IntersectRanges : On axis " << i << " of range r1, end is before start !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(r2[i].first>r2[i].second)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::IntersectRanges : On axis " << i << " of range r2, end is before start !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      ret[i].first=std::max(r1[i].first,r2[i].first);
      ret[i].second=std::min(r1[i].second,r2[i].second);
      if(ret[i].first>ret[i].second)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::IntersectRanges : On axis " << i << " the intersection of r1 and r2 is empty !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret;
}

/*!
 * This method states if \a r1 and \a r2 do overlap of not. If yes you can call IntersectRanges to know the intersection area.
 *
 * \sa IntersectRanges
 */
bool MEDCouplingStructuredMesh::AreRangesIntersect(const std::vector< std::pair<mcIdType,mcIdType> >& r1, const std::vector< std::pair<mcIdType,mcIdType> >& r2)
{
  std::size_t sz(r1.size());
  if(sz!=r2.size())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::AreRangesIntersect : the two ranges must have the same dimension !");
  for(std::size_t i=0;i<sz;i++)
    {
      if(r1[i].first>r1[i].second)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::AreRangesIntersect : On axis " << i << " of range r1, end is before start !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(r2[i].first>r2[i].second)
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::AreRangesIntersect : On axis " << i << " of range r2, end is before start !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(r1[i].second<=r2[i].first)
        return false;
      if(r1[i].first>=r2[i].second)
        return false;
    }
  return true;
}

/*!
 * This method is close to BuildExplicitIdsFrom except that instead of returning a DataArrayIdType instance containing explicit ids it
 * enable elems in the vector of booleans (for performance reasons). As it is method for performance, this method is \b not
 * available in python.
 *
 * \param [in] st The entity structure.
 * \param [in] partCompactFormat The compact subpart to be enabled.
 * \param [in,out] vectToSwitchOn Vector which fetched items are enabled.
 *
 * \sa MEDCouplingStructuredMesh::BuildExplicitIdsFrom, ExtractFieldOfBoolFrom
 */
void MEDCouplingStructuredMesh::SwitchOnIdsFrom(const std::vector<mcIdType>& st, const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat, std::vector<bool>& vectToSwitchOn)
{
  if(st.size()!=partCompactFormat.size())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::SwitchOnIdsFrom : input arrays must have the same size !");
  if(ToIdType(vectToSwitchOn.size())!=DeduceNumberOfGivenStructure(st))
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::SwitchOnIdsFrom : invalid size of input vector of boolean regarding the structure !");
  std::vector<mcIdType> dims(GetDimensionsFromCompactFrmt(partCompactFormat));
  switch(st.size())
  {
    case 3:
      {
        for(mcIdType i=0;i<dims[2];i++)
          {
            mcIdType a=(partCompactFormat[2].first+i)*st[0]*st[1];
            for(mcIdType j=0;j<dims[1];j++)
              {
                mcIdType b=(partCompactFormat[1].first+j)*st[0];
                for(mcIdType k=0;k<dims[0];k++)
                  vectToSwitchOn[partCompactFormat[0].first+k+b+a]=true;
              }
          }
        break;
      }
    case 2:
      {
        for(mcIdType j=0;j<dims[1];j++)
          {
            mcIdType b=(partCompactFormat[1].first+j)*st[0];
            for(mcIdType k=0;k<dims[0];k++)
              vectToSwitchOn[partCompactFormat[0].first+k+b]=true;
          }
        break;
      }
    case 1:
      {
        for(mcIdType k=0;k<dims[0];k++)
          vectToSwitchOn[partCompactFormat[0].first+k]=true;
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::SwitchOnIdsFrom : Dimension supported are 1,2 or 3 !");
  }
}

/*!
 * Obviously this method is \b NOT wrapped in python.
 * This method is close to SwitchOnIdsFrom except that here, a sub field \a fieldOut is built starting from the input field \a fieldOfBool having the structure \a st.
 * The extraction is defined by \a partCompactFormat.
 *
 * \param [in] st The entity structure.
 * \param [in] fieldOfBool field of booleans having the size equal to \c MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(st).
 * \param [in] partCompactFormat The compact subpart to be enabled.
 * \param [out] fieldOut the result of the extraction.
 *
 * \sa MEDCouplingStructuredMesh::BuildExplicitIdsFrom, SwitchOnIdsFrom, ExtractFieldOfDoubleFrom
 */
void MEDCouplingStructuredMesh::ExtractFieldOfBoolFrom(const std::vector<mcIdType>& st, const std::vector<bool>& fieldOfBool, const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat, std::vector<bool>& fieldOut)
{
  if(st.size()!=partCompactFormat.size())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ExtractFieldOfBoolFrom : input arrays must have the same size !");
  if(ToIdType(fieldOfBool.size())!=DeduceNumberOfGivenStructure(st))
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ExtractFieldOfBoolFrom : invalid size of input field of boolean regarding the structure !");
  std::vector<mcIdType> dims(GetDimensionsFromCompactFrmt(partCompactFormat));
  mcIdType nbOfTuplesOfOutField(DeduceNumberOfGivenStructure(dims));
  fieldOut.resize(nbOfTuplesOfOutField);
  mcIdType it(0);
  switch(st.size())
  {
    case 3:
      {
        for(mcIdType i=0;i<dims[2];i++)
          {
            mcIdType a=(partCompactFormat[2].first+i)*st[0]*st[1];
            for(mcIdType j=0;j<dims[1];j++)
              {
                mcIdType b=(partCompactFormat[1].first+j)*st[0];
                for(mcIdType k=0;k<dims[0];k++)
                  fieldOut[it++]=fieldOfBool[partCompactFormat[0].first+k+b+a];
              }
          }
        break;
      }
    case 2:
      {
        for(mcIdType j=0;j<dims[1];j++)
          {
            mcIdType b=(partCompactFormat[1].first+j)*st[0];
            for(mcIdType k=0;k<dims[0];k++)
              fieldOut[it++]=fieldOfBool[partCompactFormat[0].first+k+b];
          }
        break;
      }
    case 1:
      {
        for(mcIdType k=0;k<dims[0];k++)
          fieldOut[it++]=fieldOfBool[partCompactFormat[0].first+k];
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ExtractFieldOfBoolFrom : Dimension supported are 1,2 or 3 !");
  }
}

/*!
 * This method is close to SwitchOnIdsFrom except that here, a sub field \a fieldOut is built starting from the input field \a fieldOfDbl having the structure \a st.
 * The extraction is defined by \a partCompactFormat.
 *
 * \param [in] st The entity structure.
 * \param [in] fieldOfDbl field of doubles having a number of tuples equal to \c MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(st).
 * \param [in] partCompactFormat The compact subpart to be enabled.
 * \return DataArrayDouble * -the result of the extraction.
 *
 * \sa MEDCouplingStructuredMesh::BuildExplicitIdsFrom, SwitchOnIdsFrom, ExtractFieldOfBoolFrom
 */
DataArrayDouble *MEDCouplingStructuredMesh::ExtractFieldOfDoubleFrom(const std::vector<mcIdType>& st, const DataArrayDouble *fieldOfDbl, const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat)
{
  if(!fieldOfDbl || !fieldOfDbl->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ExtractFieldOfDoubleFrom : input array of double is NULL or not allocated!");
  if(st.size()!=partCompactFormat.size())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ExtractFieldOfDoubleFrom : input arrays must have the same size !");
  if(fieldOfDbl->getNumberOfTuples()!=DeduceNumberOfGivenStructure(st))
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ExtractFieldOfDoubleFrom : invalid size of input array of double regarding the structure !");
  std::vector<mcIdType> dims(GetDimensionsFromCompactFrmt(partCompactFormat));
  mcIdType nbOfTuplesOfOutField(DeduceNumberOfGivenStructure(dims));
  std::size_t nbComp(fieldOfDbl->getNumberOfComponents());
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New()); ret->alloc(nbOfTuplesOfOutField,nbComp);
  ret->copyStringInfoFrom(*fieldOfDbl);
  double *ptRet(ret->getPointer());
  const double *fieldOfDblPtr(fieldOfDbl->begin());
  switch(st.size())
  {
    case 3:
      {
        for(mcIdType i=0;i<dims[2];i++)
          {
            mcIdType a=(partCompactFormat[2].first+i)*st[0]*st[1];
            for(mcIdType j=0;j<dims[1];j++)
              {
                mcIdType b=(partCompactFormat[1].first+j)*st[0];
                for(mcIdType k=0;k<dims[0];k++)
                  ptRet=std::copy(fieldOfDblPtr+(partCompactFormat[0].first+k+b+a)*nbComp,fieldOfDblPtr+(partCompactFormat[0].first+k+b+a+1)*nbComp,ptRet);
              }
          }
        break;
      }
    case 2:
      {
        for(mcIdType j=0;j<dims[1];j++)
          {
            mcIdType b=(partCompactFormat[1].first+j)*st[0];
            for(mcIdType k=0;k<dims[0];k++)
              ptRet=std::copy(fieldOfDblPtr+(partCompactFormat[0].first+k+b)*nbComp,fieldOfDblPtr+(partCompactFormat[0].first+k+b+1)*nbComp,ptRet);
          }
        break;
      }
    case 1:
      {
        for(mcIdType k=0;k<dims[0];k++)
          ptRet=std::copy(fieldOfDblPtr+(partCompactFormat[0].first+k)*nbComp,fieldOfDblPtr+(partCompactFormat[0].first+k+1)*nbComp,ptRet);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ExtractFieldOfDoubleFrom : Dimension supported are 1,2 or 3 !");
  }
  return ret.retn();
}

/*!
 * This method assign a part of values in \a fieldOfDbl using entirely values of \b other.
 *
 * \param [in] st - the structure of \a fieldOfDbl.
 * \param [in,out] fieldOfDbl - the array that will be partially filled using \a other.
 * \param [in] partCompactFormat - the specification of the part.
 * \param [in] other - the array that will be used to fill \a fieldOfDbl.
 */
void MEDCouplingStructuredMesh::AssignPartOfFieldOfDoubleUsing(const std::vector<mcIdType>& st, DataArrayDouble *fieldOfDbl, const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat, const DataArrayDouble *other)
{//to be optimized
  std::vector<mcIdType> facts(st.size(),1);
  MEDCouplingIMesh::CondenseFineToCoarse(st,other,partCompactFormat,facts,fieldOfDbl);
}

/*!
 * This method changes the reference of a part of structured mesh \a partOfBigInAbs define in absolute reference to a new reference \a bigInAbs.
 * So this method only performs a translation by doing \a partOfBigRelativeToBig = \a partOfBigInAbs - \a bigInAbs
 * This method also checks (if \a check=true) that \a partOfBigInAbs is included in \a bigInAbs.
 * This method is useful to extract a part from a field lying on a big mesh.
 *
 * \sa ChangeReferenceToGlobalOfCompactFrmt, BuildExplicitIdsFrom, SwitchOnIdsFrom, ExtractFieldOfBoolFrom, ExtractFieldOfDoubleFrom
 */
void MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt(const std::vector< std::pair<mcIdType,mcIdType> >& bigInAbs, const std::vector< std::pair<mcIdType,mcIdType> >& partOfBigInAbs, std::vector< std::pair<mcIdType,mcIdType> >& partOfBigRelativeToBig, bool check)
{
  std::size_t dim(bigInAbs.size());
  if(dim!=partOfBigInAbs.size())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt : The size of parts (dimension) must be the same !");
  partOfBigRelativeToBig.resize(dim);
  for(std::size_t i=0;i<dim;i++)
    {
      if(check)
        {
          if(bigInAbs[i].first>bigInAbs[i].second)
            {
              std::ostringstream oss; oss << "MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt : Error at axis #" << i << " the input big part invalid, end before start !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          if(partOfBigInAbs[i].first<bigInAbs[i].first || partOfBigInAbs[i].first>=bigInAbs[i].second)
            {
              std::ostringstream oss; oss << "MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt : Error at axis #" << i << " the part is not included in the big one (start) !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      partOfBigRelativeToBig[i].first=partOfBigInAbs[i].first-bigInAbs[i].first;
      if(check)
        {
          if(partOfBigInAbs[i].second<partOfBigInAbs[i].first || partOfBigInAbs[i].second>bigInAbs[i].second)
            {
              std::ostringstream oss; oss << "MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt : Error at axis #" << i << " the part is not included in the big one (end) !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      partOfBigRelativeToBig[i].second=partOfBigInAbs[i].second-bigInAbs[i].first;
    }
}

/*
 * This method is performs the opposite reference modification than explained in ChangeReferenceFromGlobalOfCompactFrmt.
 *
 * \sa ChangeReferenceFromGlobalOfCompactFrmt
 */
void MEDCouplingStructuredMesh::ChangeReferenceToGlobalOfCompactFrmt(const std::vector< std::pair<mcIdType,mcIdType> >& bigInAbs, const std::vector< std::pair<mcIdType,mcIdType> >& partOfBigRelativeToBig, std::vector< std::pair<mcIdType,mcIdType> >& partOfBigInAbs, bool check)
{
  std::size_t dim(bigInAbs.size());
  if(dim!=partOfBigRelativeToBig.size())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ChangeReferenceToGlobalOfCompactFrmt : The size of parts (dimension) must be the same !");
  partOfBigInAbs.resize(dim);
  for(std::size_t i=0;i<dim;i++)
    {
      if(check)
        {
          if(bigInAbs[i].first>bigInAbs[i].second)
            {
              std::ostringstream oss; oss << "MEDCouplingStructuredMesh::ChangeReferenceToGlobalOfCompactFrmt : Error at axis #" << i << " the input big part invalid, end before start !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          if(partOfBigRelativeToBig[i].first<0 || partOfBigRelativeToBig[i].first>=bigInAbs[i].second-bigInAbs[i].first)
            {
              std::ostringstream oss; oss << "MEDCouplingStructuredMesh::ChangeReferenceToGlobalOfCompactFrmt : Error at axis #" << i << " the start of part is not in the big one !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      partOfBigInAbs[i].first=partOfBigRelativeToBig[i].first+bigInAbs[i].first;
      if(check)
        {
          if(partOfBigRelativeToBig[i].second<partOfBigRelativeToBig[i].first || partOfBigRelativeToBig[i].second>bigInAbs[i].second-bigInAbs[i].first)
            {
              std::ostringstream oss; oss << "MEDCouplingStructuredMesh::ChangeReferenceToGlobalOfCompactFrmt : Error at axis #" << i << " the end of part is not in the big one !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      partOfBigInAbs[i].second=partOfBigRelativeToBig[i].second+bigInAbs[i].first;
    }
}

/*!
 * This method performs a translation (defined by \a translation) of \a part and returns the result of translated part.
 *
 * \sa FindTranslationFrom
 */
std::vector< std::pair<mcIdType,mcIdType> > MEDCouplingStructuredMesh::TranslateCompactFrmt(const std::vector< std::pair<mcIdType,mcIdType> >& part, const std::vector<mcIdType>& translation)
{
  std::size_t sz(part.size());
  if(translation.size()!=sz)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::TranslateCompactFrmt : the sizes are not equal !");
  std::vector< std::pair<mcIdType,mcIdType> > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i].first=part[i].first+translation[i];
      ret[i].second=part[i].second+translation[i];
    }
  return ret;
}

/*!
 * \sa TranslateCompactFrmt
 */
std::vector<mcIdType> MEDCouplingStructuredMesh::FindTranslationFrom(const std::vector< std::pair<mcIdType,mcIdType> >& startingFrom, const std::vector< std::pair<mcIdType,mcIdType> >& goingTo)
{
  std::size_t sz(startingFrom.size());
  if(goingTo.size()!=sz)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::FindTranslationFrom : the sizes are not equal !");
  std::vector< mcIdType > ret(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      ret[i]=goingTo[i].first-startingFrom[i].first;
    }
  return ret;
}

/*!
 * This method builds the explicit entity array from the structure in \a st and the range in \a partCompactFormat.
 * If the range contains invalid values regarding structure an exception will be thrown.
 *
 * \return DataArrayIdType * - a new object.
 * \sa MEDCouplingStructuredMesh::IsPartStructured, MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt, SwitchOnIdsFrom, ExtractFieldOfBoolFrom, ExtractFieldOfDoubleFrom, MultiplyPartOf
 */
DataArrayIdType *MEDCouplingStructuredMesh::BuildExplicitIdsFrom(const std::vector<mcIdType>& st, const std::vector< std::pair<mcIdType,mcIdType> >& partCompactFormat)
{
  if(st.size()!=partCompactFormat.size())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : input arrays must have the same size !");
  mcIdType nbOfItems(1);
  std::vector<mcIdType> dims(st.size());
  for(std::size_t i=0;i<st.size();i++)
    {
      if(partCompactFormat[i].first<0 || partCompactFormat[i].first>st[i])
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : invalid input range 1 !");
      if(partCompactFormat[i].second<0 || partCompactFormat[i].second>st[i])
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : invalid input range 2 !");
      if(partCompactFormat[i].second<partCompactFormat[i].first)
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : invalid input range 3 !");
      dims[i]=partCompactFormat[i].second-partCompactFormat[i].first;
      nbOfItems*=dims[i];
    }
  MCAuto<DataArrayIdType> ret(DataArrayIdType::New());
  ret->alloc(nbOfItems,1);
  mcIdType *pt(ret->getPointer());
  switch(st.size())
  {
    case 3:
      {
        for(mcIdType i=0;i<dims[2];i++)
          {
            mcIdType a=(partCompactFormat[2].first+i)*st[0]*st[1];
            for(mcIdType j=0;j<dims[1];j++)
              {
                mcIdType b=(partCompactFormat[1].first+j)*st[0];
                for(mcIdType k=0;k<dims[0];k++,pt++)
                  *pt=partCompactFormat[0].first+k+b+a;
              }
          }
        break;
      }
    case 2:
      {
        for(mcIdType j=0;j<dims[1];j++)
          {
            mcIdType b=(partCompactFormat[1].first+j)*st[0];
            for(mcIdType k=0;k<dims[0];k++,pt++)
              *pt=partCompactFormat[0].first+k+b;
          }
        break;
      }
    case 1:
      {
        for(mcIdType k=0;k<dims[0];k++,pt++)
          *pt=partCompactFormat[0].first+k;
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : Dimension supported are 1,2 or 3 !");
  }
  return ret.retn();
}

/*!
 * This method multiplies by \a factor values in tuples located by \a part in \a da.
 *
 * \param [in] st - the structure of grid ( \b without considering ghost cells).
 * \param [in] part - the part in the structure ( \b without considering ghost cells) contained in grid whose structure is defined by \a st.
 * \param [in] factor - the factor, the tuples in \a da will be multiply by.
 * \param [in,out] da - The DataArray in which only tuples specified by \a part will be modified.
 *
 * \sa BuildExplicitIdsFrom
 */
void MEDCouplingStructuredMesh::MultiplyPartOf(const std::vector<mcIdType>& st, const std::vector< std::pair<mcIdType,mcIdType> >& part, double factor, DataArrayDouble *da)
{
  if(!da || !da->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::MultiplyPartOf : DataArrayDouble instance must be not NULL and allocated !");
  if(st.size()!=part.size())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::MultiplyPartOf : input arrays must have the same size !");
  std::vector<mcIdType> dims(st.size());
  for(std::size_t i=0;i<st.size();i++)
    {
      if(part[i].first<0 || part[i].first>st[i])
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::MultiplyPartOf : invalid input range 1 !");
      if(part[i].second<0 || part[i].second>st[i])
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::MultiplyPartOf : invalid input range 2 !");
      if(part[i].second<part[i].first)
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::MultiplyPartOf : invalid input range 3 !");
      dims[i]=part[i].second-part[i].first;
    }
  mcIdType nbOfTuplesExp(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(st));
  std::size_t nbCompo(da->getNumberOfComponents());
  if(da->getNumberOfTuples()!=nbOfTuplesExp)
    {
      std::ostringstream oss; oss << "MEDCouplingStructuredMesh::MultiplyPartOf : invalid nb of tuples ! Expected " << nbOfTuplesExp << " having " << da->getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  double *pt(da->getPointer());
  switch(st.size())
  {
    case 3:
      {
        for(mcIdType i=0;i<dims[2];i++)
          {
            mcIdType a=(part[2].first+i)*st[0]*st[1];
            for(mcIdType j=0;j<dims[1];j++)
              {
                mcIdType b=(part[1].first+j)*st[0];
                for(mcIdType k=0;k<dims[0];k++)
                  {
                    mcIdType offset(part[0].first+k+b+a);
                    std::transform(pt+nbCompo*offset,pt+nbCompo*(offset+1),pt+nbCompo*offset,std::bind(std::multiplies<double>(),std::placeholders::_1,factor));
                  }
              }
          }
        break;
      }
    case 2:
      {
        for(mcIdType j=0;j<dims[1];j++)
          {
            mcIdType b=(part[1].first+j)*st[0];
            for(mcIdType k=0;k<dims[0];k++)
              {
                mcIdType offset(part[0].first+k+b);
                std::transform(pt+nbCompo*offset,pt+nbCompo*(offset+1),pt+nbCompo*offset,std::bind(std::multiplies<double>(),std::placeholders::_1,factor));
              }
          }
        break;
      }
    case 1:
      {
        for(mcIdType k=0;k<dims[0];k++)
          {
            mcIdType offset(part[0].first+k);
            std::transform(pt+nbCompo*offset,pt+nbCompo*(offset+1),pt+nbCompo*offset,std::bind(std::multiplies<double>(),std::placeholders::_1,factor));
          }
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::MultiplyPartOf : Dimension supported are 1,2 or 3 !");
  }
}

/*!
 * This method multiplies by \a factor values in tuples located by \a part in \a da.
 *
 * \param [in] st - the structure of grid ( \b without considering ghost cells).
 * \param [in] part - the part in the structure ( \b without considering ghost cells) contained in grid whose structure is defined by \a st.
 * \param [in] ghostSize - \a ghostSize must be >= 0.
 * \param [in] factor - the factor, the tuples in \a da will be multiply by.
 * \param [in,out] da - The DataArray in which only tuples specified by \a part will be modified.
 *
 * \sa MultiplyPartOf, PutInGhostFormat
 */
void MEDCouplingStructuredMesh::MultiplyPartOfByGhost(const std::vector<mcIdType>& st, const std::vector< std::pair<mcIdType,mcIdType> >& part, mcIdType ghostSize, double factor, DataArrayDouble *da)
{
  std::vector<mcIdType> stWG;
  std::vector< std::pair<mcIdType,mcIdType> > partWG;
  PutInGhostFormat(ghostSize,st,part,stWG,partWG);
  MultiplyPartOf(stWG,partWG,factor,da);
}

/*!
 * This method multiplies by \a factor values in tuples located by \a part in \a da.
 *
 * \param [in] st - the structure of grid ( \b without considering ghost cells).
 * \param [in] part - the part in the structure ( \b without considering ghost cells) contained in grid whose structure is defined by \a st.
 * \param [in] ghostSize - \a ghostSize must be >= 0.
 * \param [out] stWithGhost - the structure considering ghost cells.
 * \param [out] partWithGhost - the part considering the ghost cells.
 *
 * \sa MultiplyPartOf, PutInGhostFormat
 */
void MEDCouplingStructuredMesh::PutInGhostFormat(mcIdType ghostSize, const std::vector<mcIdType>& st, const std::vector< std::pair<mcIdType,mcIdType> >& part, std::vector<mcIdType>& stWithGhost, std::vector< std::pair<mcIdType,mcIdType> >&partWithGhost)
{
  if(ghostSize<0)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::PutInGhostFormat : ghost size must be >= 0 !");
  std::size_t dim(part.size());
  if(st.size()!=dim)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::PutInGhostFormat : the dimension of input vectors must be the same !");
  for(std::size_t i=0;i<dim;i++)
    if(part[i].first<0 || part[i].first>part[i].second || part[i].second>st[i])
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::PutInGhostFormat : the specified part is invalid ! The begin must be >= 0 and <= end ! The end must be <= to the size at considered dimension !");
  stWithGhost.resize(st.size());
  std::transform(st.begin(),st.end(),stWithGhost.begin(),std::bind(std::plus<mcIdType>(),std::placeholders::_1,2*ghostSize));
  partWithGhost=part;
  ApplyGhostOnCompactFrmt(partWithGhost,ghostSize);
}

/*!
 * \param [in,out] partBeforeFact - the part of a image mesh in compact format that will be put in ghost reference.
 * \param [in] ghostSize - the ghost size of zone for all axis.
 */
void MEDCouplingStructuredMesh::ApplyGhostOnCompactFrmt(std::vector< std::pair<mcIdType,mcIdType> >& partBeforeFact, mcIdType ghostSize)
{
  if(ghostSize<0)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::ApplyGhostOnCompactFrmt : ghost size must be >= 0 !");
  std::size_t sz(partBeforeFact.size());
  for(std::size_t i=0;i<sz;i++)
    {
      partBeforeFact[i].first+=ghostSize;
      partBeforeFact[i].second+=ghostSize;
    }
}

mcIdType MEDCouplingStructuredMesh::GetNumberOfCellsOfSubLevelMesh(const std::vector<mcIdType>& cgs, int mdim)
{
  mcIdType ret(0);
  for(int i=0;i<mdim;i++)
    {
      mcIdType locRet(1);
      for(int j=0;j<mdim;j++)
        if(j!=i)
          locRet*=cgs[j];
        else
          locRet*=cgs[j]+1;
      ret+=locRet;
    }
  return ret;
}
