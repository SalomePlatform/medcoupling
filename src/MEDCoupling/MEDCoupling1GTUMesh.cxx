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

#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingCMesh.hxx"

#include "SplitterTetra.hxx"
#include "DiameterCalculator.hxx"
#include "OrientationInverter.hxx"
#include "InterpKernelAutoPtr.hxx"

using namespace MEDCoupling;

const int MEDCoupling1SGTUMesh::HEXA8_FACE_PAIRS[6]={0,1,2,4,3,5};

MEDCoupling1GTUMesh::MEDCoupling1GTUMesh():_cm(0)
{
}

MEDCoupling1GTUMesh::MEDCoupling1GTUMesh(const std::string& name, const INTERP_KERNEL::CellModel& cm):_cm(&cm)
{
  setName(name);
}

MEDCoupling1GTUMesh::MEDCoupling1GTUMesh(const MEDCoupling1GTUMesh& other, bool recDeepCpy):MEDCouplingPointSet(other,recDeepCpy),_cm(other._cm)
{
}

MEDCoupling1GTUMesh *MEDCoupling1GTUMesh::New(const std::string& name, INTERP_KERNEL::NormalizedCellType type)
{
  if(type==INTERP_KERNEL::NORM_ERROR)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::New : NORM_ERROR is not a valid type to be used as base geometric type for a mesh !");
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  if(!cm.isDynamic())
    return MEDCoupling1SGTUMesh::New(name,type);
  else
    return MEDCoupling1DGTUMesh::New(name,type);
}

MEDCoupling1GTUMesh *MEDCoupling1GTUMesh::New(const MEDCouplingUMesh *m)
{
  if(!m)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::New : input mesh is null !");
  std::set<INTERP_KERNEL::NormalizedCellType> gts(m->getAllGeoTypes());
  if(gts.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::New : input mesh must have exactly one geometric type !");
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(*gts.begin());
  if(!cm.isDynamic())
    return MEDCoupling1SGTUMesh::New(m);
  else
    return MEDCoupling1DGTUMesh::New(m);
}

const INTERP_KERNEL::CellModel& MEDCoupling1GTUMesh::getCellModel() const
{
  return *_cm;
}

INTERP_KERNEL::NormalizedCellType MEDCoupling1GTUMesh::getCellModelEnum() const
{
  return _cm->getEnum();
}

int MEDCoupling1GTUMesh::getMeshDimension() const
{
  return (int)_cm->getDimension();
}

/*!
 * This method returns a newly allocated array containing cell ids (ascendingly sorted) whose geometric type are equal to type.
 * This method does not throw exception if geometric type \a type is not in \a this.
 * This method throws an INTERP_KERNEL::Exception if meshdimension of \b this is not equal to those of \b type.
 * The coordinates array is not considered here.
 *
 * \param [in] type the geometric type
 * \return cell ids in this having geometric type \a type.
 */
DataArrayInt *MEDCoupling1GTUMesh::giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  if(type==getCellModelEnum())
    ret->alloc(getNumberOfCells(),1);
  else
    ret->alloc(0,1);
  ret->iota();
  return ret.retn();
}

/*!
 * Returns nb of cells having the geometric type \a type. No throw if no cells in \a this has the geometric type \a type.
 */
std::size_t MEDCoupling1GTUMesh::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  return type==getCellModelEnum()?getNumberOfCells():0;
}

/*!
 * Returns a type of a cell by its id.
 *  \param [in] cellId - the id of the cell of interest.
 *  \return INTERP_KERNEL::NormalizedCellType - enumeration item describing the cell type.
 *  \throw If \a cellId is invalid. Valid range is [0, \a this->getNumberOfCells() ).
 */
INTERP_KERNEL::NormalizedCellType MEDCoupling1GTUMesh::getTypeOfCell(std::size_t cellId) const
{
  if(cellId<getNumberOfCells())
    return getCellModelEnum();
  std::ostringstream oss; oss << "MEDCoupling1GTUMesh::getTypeOfCell : Requesting type of cell #" << cellId << " but it should be in [0," << getNumberOfCells() << ") !";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

/*!
 * Returns a set of all cell types available in \a this mesh.
 * \return std::set<INTERP_KERNEL::NormalizedCellType> - the set of cell types.
 * \warning this method does not throw any exception even if \a this is not defined.
 */
std::set<INTERP_KERNEL::NormalizedCellType> MEDCoupling1GTUMesh::getAllGeoTypes() const
{
  std::set<INTERP_KERNEL::NormalizedCellType> ret;
  ret.insert(getCellModelEnum());
  return ret;
}

/*!
 * This method expects that \a this is sorted by types. If not an exception will be thrown.
 * This method returns in the same format as code (see MEDCouplingUMesh::checkTypeConsistencyAndContig or MEDCouplingUMesh::splitProfilePerType) how
 * \a this is composed in cell types.
 * The returned array is of size 3*n where n is the number of different types present in \a this. 
 * For every k in [0,n] ret[3*k+2]==-1 because it has no sense here. 
 * This parameter is kept only for compatibility with other methode listed above.
 */
std::vector<int> MEDCoupling1GTUMesh::getDistributionOfTypes() const
{
  std::vector<int> ret(3);
  ret[0]=(int)getCellModelEnum(); ret[1]=getNumberOfCells(); ret[2]=-1;
  return ret;
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
void MEDCoupling1GTUMesh::splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const
{
  if(!profile)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::splitProfilePerType : input profile is NULL !");
  if(profile->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::splitProfilePerType : input profile should have exactly one component !");
  int nbTuples(profile->getNumberOfTuples()),nbOfCells(getNumberOfCells());
  code.resize(3); idsInPflPerType.resize(1);
  code[0]=(int)getCellModelEnum(); code[1]=nbTuples;
  idsInPflPerType.resize(1);
  if(profile->isIota(nbOfCells))
    {
      code[2]=-1;
      idsInPflPerType[0]=const_cast<DataArrayInt *>(profile); idsInPflPerType[0]->incrRef();
      idsPerType.clear();
      return ;
    }
  code[2]=0;
  profile->checkAllIdsInRange(0,nbOfCells);
  idsPerType.resize(1);
  idsPerType[0]=const_cast<DataArrayInt *>(profile); idsPerType[0]->incrRef();
  idsInPflPerType[0]=DataArrayInt::Range(0,nbTuples,1);
}

/*!
 * This method tries to minimize at most the number of deep copy.
 * So if \a idsPerType is not empty it can be returned directly (without copy, but with ref count incremented) in return.
 * 
 * \sa MEDCouplingUMesh::checkTypeConsistencyAndContig
 */
DataArrayInt *MEDCoupling1GTUMesh::checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const
{
  int nbOfCells=getNumberOfCells();
  if(code.size()!=3)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::checkTypeConsistencyAndContig : invalid input code should be exactly of size 3 !");
  if(code[0]!=(int)getCellModelEnum())
    {
      std::ostringstream oss; oss << "MEDCoupling1GTUMesh::checkTypeConsistencyAndContig : Mismatch of geometric type ! Asking for " << code[0] << " whereas the geometric type is \a this is " << getCellModelEnum() << " (" << _cm->getRepr() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(code[2]==-1)
    {
      if(code[1]==nbOfCells)
        return 0;
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1GTUMesh::checkTypeConsistencyAndContig : mismatch between the number of cells in this (" << nbOfCells << ") and the number of non profile (" << code[1] << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(code[2]!=0)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::checkTypeConsistencyAndContig : single geo type mesh ! 0 or -1 is expected at pos #2 of input code !");
  if(idsPerType.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::checkTypeConsistencyAndContig : input code points to DataArrayInt #0 whereas the size of idsPerType is not equal to 1 !");
  const DataArrayInt *pfl=idsPerType[0];
  if(!pfl)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::checkTypeConsistencyAndContig : the input code points to a NULL DataArrayInt at rank 0 !");
  if(pfl->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::checkTypeConsistencyAndContig : input profile should have exactly one component !");
  pfl->checkAllIdsInRange(0,nbOfCells);
  pfl->incrRef();
  return const_cast<DataArrayInt *>(pfl);
}

void MEDCoupling1GTUMesh::writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const
{
  MCAuto<MEDCouplingUMesh> m=buildUnstructured();
  m->writeVTKLL(ofs,cellData,pointData,byteData);
}

std::string MEDCoupling1GTUMesh::getVTKDataSetType() const
{
  return std::string("UnstructuredGrid");
}

std::string MEDCoupling1GTUMesh::getVTKFileExtension() const
{
  return std::string("vtu");
}

std::size_t MEDCoupling1GTUMesh::getHeapMemorySizeWithoutChildren() const
{
  return MEDCouplingPointSet::getHeapMemorySizeWithoutChildren();
}

bool MEDCoupling1GTUMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const
{
  if(!MEDCouplingPointSet::isEqualIfNotWhy(other,prec,reason))
    return false;
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::isEqualIfNotWhy : input other pointer is null !");
  const MEDCoupling1GTUMesh *otherC=dynamic_cast<const MEDCoupling1GTUMesh *>(other);
  if(!otherC)
    {
      reason="mesh given in input is not castable in MEDCouplingSGTUMesh !";
      return false;
    }
  if(_cm!=otherC->_cm)
    {
      reason="mismatch in geometric type !";
      return false;
    }
  return true;
}

bool MEDCoupling1GTUMesh::isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const
{
  if(!MEDCouplingPointSet::isEqualWithoutConsideringStr(other,prec))
    return false;
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::isEqualWithoutConsideringStr : input other pointer is null !");
  const MEDCoupling1GTUMesh *otherC=dynamic_cast<const MEDCoupling1GTUMesh *>(other);
  if(!otherC)
    return false;
  if(_cm!=otherC->_cm)
    return false;
  return true;
}

void MEDCoupling1GTUMesh::checkConsistencyLight() const
{
  MEDCouplingPointSet::checkConsistencyLight();
}

DataArrayDouble *MEDCoupling1GTUMesh::computeCellCenterOfMass() const
{
  MCAuto<MEDCouplingUMesh> m=buildUnstructured();
  MCAuto<DataArrayDouble> ret=m->computeCellCenterOfMass();
  return ret.retn();
}

MEDCouplingFieldDouble *MEDCoupling1GTUMesh::getMeasureField(bool isAbs) const
{
  MCAuto<MEDCouplingUMesh> m=buildUnstructured();
  MCAuto<MEDCouplingFieldDouble> ret=m->getMeasureField(isAbs);
  ret->setMesh(this);
  return ret.retn();
}

MEDCouplingFieldDouble *MEDCoupling1GTUMesh::getMeasureFieldOnNode(bool isAbs) const
{
  MCAuto<MEDCouplingUMesh> m=buildUnstructured();
  MCAuto<MEDCouplingFieldDouble> ret=m->getMeasureFieldOnNode(isAbs);
  ret->setMesh(this);
  return ret.retn();
}

/*!
 * to improve perf !
 */
int MEDCoupling1GTUMesh::getCellContainingPoint(const double *pos, double eps) const
{
  MCAuto<MEDCouplingUMesh> m(buildUnstructured());
  return m->getCellContainingPoint(pos,eps);
}

/*!
 * to improve perf !
 */
void MEDCoupling1GTUMesh::getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const
{
  MCAuto<MEDCouplingUMesh> m(buildUnstructured());
  return m->getCellsContainingPoint(pos,eps,elts);
}

MEDCouplingFieldDouble *MEDCoupling1GTUMesh::buildOrthogonalField() const
{
  MCAuto<MEDCouplingUMesh> m=buildUnstructured();
  MCAuto<MEDCouplingFieldDouble> ret=m->buildOrthogonalField();
  ret->setMesh(this);
  return ret.retn();
}

DataArrayInt *MEDCoupling1GTUMesh::getCellsInBoundingBox(const double *bbox, double eps) const
{
  MCAuto<MEDCouplingUMesh> m=buildUnstructured();
  return m->getCellsInBoundingBox(bbox,eps);
}

DataArrayInt *MEDCoupling1GTUMesh::getCellsInBoundingBox(const INTERP_KERNEL::DirectedBoundingBox& bbox, double eps)
{
  MCAuto<MEDCouplingUMesh> m=buildUnstructured();
  return m->getCellsInBoundingBox(bbox,eps);
}

MEDCouplingPointSet *MEDCoupling1GTUMesh::buildFacePartOfMySelfNode(const int *start, const int *end, bool fullyIn) const
{
  MCAuto<MEDCouplingUMesh> m=buildUnstructured();
  return m->buildFacePartOfMySelfNode(start,end,fullyIn);
}

DataArrayInt *MEDCoupling1GTUMesh::findBoundaryNodes() const
{
  MCAuto<MEDCouplingUMesh> m=buildUnstructured();
  return m->findBoundaryNodes();
}

MEDCouplingPointSet *MEDCoupling1GTUMesh::buildBoundaryMesh(bool keepCoords) const
{
  MCAuto<MEDCouplingUMesh> m=buildUnstructured();
  return m->buildBoundaryMesh(keepCoords);
}

void MEDCoupling1GTUMesh::findCommonCells(int compType, int startCellId, DataArrayInt *& commonCellsArr, DataArrayInt *& commonCellsIArr) const
{
  MCAuto<MEDCouplingUMesh> m=buildUnstructured();
  m->findCommonCells(compType,startCellId,commonCellsArr,commonCellsIArr);
}

std::size_t MEDCoupling1GTUMesh::getNodalConnectivityLength() const
{
  const DataArrayInt *c1(getNodalConnectivity());
  if(!c1)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::getNodalConnectivityLength : no connectivity set !");
  if(c1->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::getNodalConnectivityLength : Nodal connectivity array set must have exactly one component !");
  if(!c1->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::getNodalConnectivityLength : Nodal connectivity array must be allocated !");
  return c1->getNumberOfTuples();
}

/*!
 * This method aggregates all the meshes in \a parts to put them in a single unstructured mesh (those returned).
 * The order of cells is the returned instance is those in the order of instances in \a parts.
 *
 * \param [in] parts - all not null parts of single geo type meshes to be aggreagated having the same mesh dimension and same coordinates.
 * \return MEDCouplingUMesh * - new object to be dealt by the caller.
 *
 * \throw If one element is null in \a parts.
 * \throw If not all the parts do not have the same mesh dimension.
 * \throw If not all the parts do not share the same coordinates.
 * \throw If not all the parts have their connectivity set properly.
 * \throw If \a parts is empty.
 */
MEDCouplingUMesh *MEDCoupling1GTUMesh::AggregateOnSameCoordsToUMesh(const std::vector< const MEDCoupling1GTUMesh *>& parts)
{
  if(parts.empty())
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::AggregateOnSameCoordsToUMesh : input parts vector is empty !");
  const MEDCoupling1GTUMesh *firstPart(parts[0]);
  if(!firstPart)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::AggregateOnSameCoordsToUMesh : the first instance in input parts is null !");
  const DataArrayDouble *coords(firstPart->getCoords());
  int meshDim(firstPart->getMeshDimension());
  MCAuto<MEDCouplingUMesh> ret(MEDCouplingUMesh::New(firstPart->getName(),meshDim)); ret->setDescription(firstPart->getDescription());
  ret->setCoords(coords);
  int nbOfCells(0),connSize(0);
  for(std::vector< const MEDCoupling1GTUMesh *>::const_iterator it=parts.begin();it!=parts.end();it++)
    {
      if(!(*it))
        throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::AggregateOnSameCoordsToUMesh : presence of null pointer in input vector !");
      if((*it)->getMeshDimension()!=meshDim)
        throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::AggregateOnSameCoordsToUMesh : all the instances in input vector must have same mesh dimension !");
      if((*it)->getCoords()!=coords)
        throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::AggregateOnSameCoordsToUMesh : all the instances must share the same coordinates pointer !");
      nbOfCells+=(*it)->getNumberOfCells();
      connSize+=(*it)->getNodalConnectivityLength();
    }
  MCAuto<DataArrayInt> conn(DataArrayInt::New()),connI(DataArrayInt::New());
  connI->alloc(nbOfCells+1,1); conn->alloc(connSize+nbOfCells,1);
  int *c(conn->getPointer()),*ci(connI->getPointer()); *ci=0;
  for(std::vector< const MEDCoupling1GTUMesh *>::const_iterator it=parts.begin();it!=parts.end();it++)
    {
      int curNbCells((*it)->getNumberOfCells());
      int geoType((int)(*it)->getCellModelEnum());
      const int *cinPtr((*it)->getNodalConnectivity()->begin());
      const MEDCoupling1SGTUMesh *ps(dynamic_cast<const MEDCoupling1SGTUMesh *>(*it));
      const MEDCoupling1DGTUMesh *pd(dynamic_cast<const MEDCoupling1DGTUMesh *>(*it));
      if(ps && !pd)
        {
          int nNodesPerCell(ps->getNumberOfNodesPerCell());
          for(int i=0;i<curNbCells;i++,ci++,cinPtr+=nNodesPerCell)
            {
              *c++=geoType;
              c=std::copy(cinPtr,cinPtr+nNodesPerCell,c);
              ci[1]=ci[0]+nNodesPerCell+1;
            }
        }
      else if(!ps && pd)
        {
          const int *ciinPtr(pd->getNodalConnectivityIndex()->begin());
          for(int i=0;i<curNbCells;i++,ci++,ciinPtr++)
            {
              *c++=geoType;
              c=std::copy(cinPtr+ciinPtr[0],cinPtr+ciinPtr[1],c);
              ci[1]=ci[0]+ciinPtr[1]-ciinPtr[0]+1;
            }
        }
      else
        throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::AggregateOnSameCoordsToUMesh : presence of instance which type is not in [MEDCoupling1SGTUMesh,MEDCoupling1DGTUMesh] !");
    }
  ret->setConnectivity(conn,connI,true);
  return ret.retn();
}

//==

MEDCoupling1SGTUMesh::MEDCoupling1SGTUMesh(const MEDCoupling1SGTUMesh& other, bool recDeepCpy):MEDCoupling1GTUMesh(other,recDeepCpy),_conn(other._conn)
{
  if(recDeepCpy)
    {
      const DataArrayInt *c(other._conn);
      if(c)
        _conn=c->deepCopy();
    }
}

MEDCoupling1SGTUMesh::MEDCoupling1SGTUMesh(const std::string& name, const INTERP_KERNEL::CellModel& cm):MEDCoupling1GTUMesh(name,cm)
{
}

MEDCoupling1SGTUMesh::MEDCoupling1SGTUMesh()
{
}

MEDCoupling1SGTUMesh *MEDCoupling1SGTUMesh::New()
{
  return new MEDCoupling1SGTUMesh;
}

MEDCoupling1SGTUMesh *MEDCoupling1SGTUMesh::New(const std::string& name, INTERP_KERNEL::NormalizedCellType type)
{
  if(type==INTERP_KERNEL::NORM_ERROR)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::New : NORM_ERROR is not a valid type to be used as base geometric type for a mesh !");
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  if(cm.isDynamic())
    {
      std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::New : the input geometric type " << cm.getRepr() << " is dynamic ! Only static types are allowed here !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return new MEDCoupling1SGTUMesh(name,cm);
}

MEDCoupling1SGTUMesh *MEDCoupling1SGTUMesh::New(const MEDCouplingUMesh *m)
{
  if(!m)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::New : input mesh is null !");
  std::set<INTERP_KERNEL::NormalizedCellType> gts(m->getAllGeoTypes());
  if(gts.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::New : input mesh must have exactly one geometric type !");
  int geoType((int)*gts.begin());
  MCAuto<MEDCoupling1SGTUMesh> ret(MEDCoupling1SGTUMesh::New(m->getName(),*gts.begin()));
  ret->setCoords(m->getCoords()); ret->setDescription(m->getDescription());
  int nbCells(m->getNumberOfCells());
  int nbOfNodesPerCell(ret->getNumberOfNodesPerCell());
  MCAuto<DataArrayInt> conn(DataArrayInt::New()); conn->alloc(nbCells*nbOfNodesPerCell,1);
  int *c(conn->getPointer());
  const int *cin(m->getNodalConnectivity()->begin()),*ciin(m->getNodalConnectivityIndex()->begin());
  for(int i=0;i<nbCells;i++,ciin++)
    {
      if(cin[ciin[0]]==geoType)
        {
          if(ciin[1]-ciin[0]==nbOfNodesPerCell+1)
            c=std::copy(cin+ciin[0]+1,cin+ciin[1],c);
          else
            {
              std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::New(const MEDCouplingUMesh *m) : something is wrong in the input mesh at cell #" << i << " ! The size of cell is not those expected (" << nbOfNodesPerCell << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::New(const MEDCouplingUMesh *m) : something is wrong in the input mesh at cell #" << i << " ! The geometric type is not those expected !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  ret->setNodalConnectivity(conn);
  try
  { ret->copyTinyInfoFrom(m); }
  catch(INTERP_KERNEL::Exception&) { }
  return ret.retn();
}

MEDCoupling1SGTUMesh *MEDCoupling1SGTUMesh::clone(bool recDeepCpy) const
{
  return new MEDCoupling1SGTUMesh(*this,recDeepCpy);
}

/*!
 * This method behaves mostly like MEDCoupling1SGTUMesh::deepCopy method, except that only nodal connectivity arrays are deeply copied.
 * The coordinates are shared between \a this and the returned instance.
 * 
 * \return MEDCoupling1SGTUMesh * - A new object instance holding the copy of \a this (deep for connectivity, shallow for coordiantes)
 * \sa MEDCoupling1SGTUMesh::deepCopy
 */
MEDCoupling1SGTUMesh *MEDCoupling1SGTUMesh::deepCopyConnectivityOnly() const
{
  checkConsistencyLight();
  MCAuto<MEDCoupling1SGTUMesh> ret(clone(false));
  MCAuto<DataArrayInt> c(_conn->deepCopy());
  ret->setNodalConnectivity(c);
  return ret.retn();
}

void MEDCoupling1SGTUMesh::shallowCopyConnectivityFrom(const MEDCouplingPointSet *other)
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::shallowCopyConnectivityFrom : input pointer is null !");
  const MEDCoupling1SGTUMesh *otherC=dynamic_cast<const MEDCoupling1SGTUMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::shallowCopyConnectivityFrom : input pointer is not an MEDCoupling1SGTUMesh instance !");
  setNodalConnectivity(otherC->getNodalConnectivity());
}

void MEDCoupling1SGTUMesh::updateTime() const
{
  MEDCoupling1GTUMesh::updateTime();
  const DataArrayInt *c(_conn);
  if(c)
    updateTimeWith(*c);
}

std::size_t MEDCoupling1SGTUMesh::getHeapMemorySizeWithoutChildren() const
{
  return MEDCoupling1GTUMesh::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDCoupling1SGTUMesh::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDCoupling1GTUMesh::getDirectChildrenWithNull());
  ret.push_back((const DataArrayInt *)_conn);
  return ret;
}

MEDCoupling1SGTUMesh *MEDCoupling1SGTUMesh::deepCopy() const
{
  return clone(true);
}

bool MEDCoupling1SGTUMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::isEqualIfNotWhy : input other pointer is null !");
  std::ostringstream oss; oss.precision(15);
  const MEDCoupling1SGTUMesh *otherC=dynamic_cast<const MEDCoupling1SGTUMesh *>(other);
  if(!otherC)
    {
      reason="mesh given in input is not castable in MEDCoupling1SGTUMesh !";
      return false;
    }
  if(!MEDCoupling1GTUMesh::isEqualIfNotWhy(other,prec,reason))
    return false;
  const DataArrayInt *c1(_conn),*c2(otherC->_conn);
  if(c1==c2)
    return true;
  if(!c1 || !c2)
    {
      reason="in connectivity of single static geometric type exactly one among this and other is null !";
      return false;
    }
  if(!c1->isEqualIfNotWhy(*c2,reason))
    {
      reason.insert(0,"Nodal connectivity DataArrayInt differ : ");
      return false;
    }
  return true;
}

bool MEDCoupling1SGTUMesh::isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::isEqualWithoutConsideringStr : input other pointer is null !");
  const MEDCoupling1SGTUMesh *otherC=dynamic_cast<const MEDCoupling1SGTUMesh *>(other);
  if(!otherC)
    return false;
  if(!MEDCoupling1GTUMesh::isEqualWithoutConsideringStr(other,prec))
    return false;
  const DataArrayInt *c1(_conn),*c2(otherC->_conn);
  if(c1==c2)
    return true;
  if(!c1 || !c2)
    return false;
  if(!c1->isEqualWithoutConsideringStr(*c2))
    return false;
  return true;
}

void MEDCoupling1SGTUMesh::checkConsistencyOfConnectivity() const
{
  const DataArrayInt *c1(_conn);
  if(c1)
    {
      if(c1->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("Nodal connectivity array is expected to be with number of components set to one !");
      if(c1->getInfoOnComponent(0)!="")
        throw INTERP_KERNEL::Exception("Nodal connectivity array is expected to have no info on its single component !");
      c1->checkAllocated();
    }
  else
    throw INTERP_KERNEL::Exception("Nodal connectivity array not defined !");
}

void MEDCoupling1SGTUMesh::checkConsistencyLight() const
{
  MEDCouplingPointSet::checkConsistencyLight();
  checkConsistencyOfConnectivity();
}

void MEDCoupling1SGTUMesh::checkConsistency(double eps) const
{
  checkConsistencyLight();
  const DataArrayInt *c1(_conn);
  int nbOfTuples=c1->getNumberOfTuples();
  int nbOfNodesPerCell=(int)_cm->getNumberOfNodes();
  if(nbOfTuples%nbOfNodesPerCell!=0)
    {
      std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::checkConsistency : the nb of tuples in conn is " << nbOfTuples << " and number of nodes per cell is " << nbOfNodesPerCell << ". But " << nbOfTuples << "%" << nbOfNodesPerCell << " !=0 !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int nbOfNodes=getNumberOfNodes();
  int nbOfCells=nbOfTuples/nbOfNodesPerCell;
  const int *w(c1->begin());
  for(int i=0;i<nbOfCells;i++)
    for(int j=0;j<nbOfNodesPerCell;j++,w++)
      {
        if(*w<0 || *w>=nbOfNodes)
          {
            std::ostringstream oss; oss << "At node #" << j << " of  cell #" << i << ", is equal to " << *w << " must be in [0," << nbOfNodes << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
}

std::size_t MEDCoupling1SGTUMesh::getNumberOfCells() const
{
  std::size_t nbOfTuples(getNodalConnectivityLength());
  int nbOfNodesPerCell(getNumberOfNodesPerCell());
  if(nbOfTuples%nbOfNodesPerCell!=0)
    {
      std::ostringstream oss; oss << "MEDCoupling1SGTUMesh:getNumberOfCells: : the nb of tuples in conn is " << nbOfTuples << " and number of nodes per cell is " << nbOfNodesPerCell << ". But " << nbOfTuples << "%" << nbOfNodesPerCell << " !=0 !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return nbOfTuples/nbOfNodesPerCell;
}

int MEDCoupling1SGTUMesh::getNumberOfNodesInCell(int cellId) const
{
  return getNumberOfNodesPerCell();
}

int MEDCoupling1SGTUMesh::getNumberOfNodesPerCell() const
{
  checkNonDynamicGeoType();
  return (int)_cm->getNumberOfNodes();
}

DataArrayInt *MEDCoupling1SGTUMesh::computeNbOfNodesPerCell() const
{
  checkNonDynamicGeoType();
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(getNumberOfCells(),1);
  ret->fillWithValue((int)_cm->getNumberOfNodes());
  return ret.retn();
}

DataArrayInt *MEDCoupling1SGTUMesh::computeNbOfFacesPerCell() const
{
  checkNonDynamicGeoType();
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(getNumberOfCells(),1);
  ret->fillWithValue((int)_cm->getNumberOfSons());
  return ret.retn();
}

DataArrayInt *MEDCoupling1SGTUMesh::computeEffectiveNbOfNodesPerCell() const
{
  checkNonDynamicGeoType();
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  int nbCells(getNumberOfCells());
  ret->alloc(nbCells,1);
  int *retPtr(ret->getPointer());
  int nbNodesPerCell(getNumberOfNodesPerCell());
  const int *conn(_conn->begin());
  for(int i=0;i<nbCells;i++,conn+=nbNodesPerCell,retPtr++)
    {
      std::set<int> s(conn,conn+nbNodesPerCell);
      *retPtr=(int)s.size();
    }
  return ret.retn();
}

void MEDCoupling1SGTUMesh::getNodeIdsOfCell(std::size_t cellId, std::vector<int>& conn) const
{
  int sz=getNumberOfNodesPerCell();
  conn.resize(sz);
  if(cellId<getNumberOfCells())
    std::copy(_conn->begin()+cellId*sz,_conn->begin()+(cellId+1)*sz,conn.begin());
  else
    {
      std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::getNodeIdsOfCell : request for cellId #" << cellId << " must be in [0," << getNumberOfCells() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void MEDCoupling1SGTUMesh::checkNonDynamicGeoType() const
{
  if(_cm->isDynamic())
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::checkNonDynamicGeoType : internal error ! the internal geo type is dynamic ! should be static !");
}

std::string MEDCoupling1SGTUMesh::simpleRepr() const
{
  static const char msg0[]="No coordinates specified !";
  std::ostringstream ret;
  ret << "Single static geometic type (" << _cm->getRepr() << ") unstructured mesh with name : \"" << getName() << "\"\n";
  ret << "Description of mesh : \"" << getDescription() << "\"\n";
  int tmpp1,tmpp2;
  double tt=getTime(tmpp1,tmpp2);
  ret << "Time attached to the mesh [unit] : " << tt << " [" << getTimeUnit() << "]\n";
  ret << "Iteration : " << tmpp1  << " Order : " << tmpp2 << "\n";
  ret << "Mesh dimension : " << getMeshDimension() << "\nSpace dimension : ";
  if(_coords!=0)
    {
      const int spaceDim=getSpaceDimension();
      ret << spaceDim << "\nInfo attached on space dimension : ";
      for(int i=0;i<spaceDim;i++)
        ret << "\"" << _coords->getInfoOnComponent(i) << "\" ";
      ret << "\n";
    }
  else
    ret << msg0 << "\n";
  ret << "Number of nodes : ";
  if(_coords!=0)
    ret << getNumberOfNodes() << "\n";
  else
    ret << msg0 << "\n";
  ret << "Number of cells : ";
  if((const DataArrayInt *)_conn)
    {
      if(_conn->isAllocated())
        {
          if(_conn->getNumberOfComponents()==1)
            ret << getNumberOfCells() << "\n";
          else
            ret << "Nodal connectivity array specified and allocated but with not exactly one component !" << "\n";
        }
      else
        ret << "Nodal connectivity array specified but not allocated !" << "\n";
    }
  else
    ret << "No connectivity specified !" << "\n";
  ret << "Cell type : " << _cm->getRepr() << "\n";
  return ret.str();
}

std::string MEDCoupling1SGTUMesh::advancedRepr() const
{
  std::ostringstream ret;
  ret << simpleRepr();
  ret << "\nCoordinates array : \n___________________\n\n";
  if(_coords)
    _coords->reprWithoutNameStream(ret);
  else
    ret << "No array set !\n";
  ret << "\n\nConnectivity array : \n____________________\n\n";
  //
  if((const DataArrayInt *)_conn)
    {
      if(_conn->isAllocated())
        {
          if(_conn->getNumberOfComponents()==1)
            {
              int nbOfCells=getNumberOfCells();
              int sz=getNumberOfNodesPerCell();
              const int *connPtr=_conn->begin();
              for(int i=0;i<nbOfCells;i++,connPtr+=sz)
                {
                  ret << "Cell #" << i << " : ";
                  std::copy(connPtr,connPtr+sz,std::ostream_iterator<int>(ret," "));
                  ret << "\n";
                }
            }
          else
            ret << "Nodal connectivity array specified and allocated but with not exactly one component !" << "\n";
        }
      else
        ret << "Nodal connectivity array specified but not allocated !" << "\n";
    }
  else
    ret << "No connectivity specified !" << "\n";
  return ret.str();
}

DataArrayDouble *MEDCoupling1SGTUMesh::computeIsoBarycenterOfNodesPerCell() const
{
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  int spaceDim=getSpaceDimension();
  int nbOfCells=getNumberOfCells();//checkConsistencyLight()
  int nbOfNodes=getNumberOfNodes();
  ret->alloc(nbOfCells,spaceDim);
  double *ptToFill=ret->getPointer();
  const double *coor=_coords->begin();
  const int *nodal=_conn->begin();
  int sz=getNumberOfNodesPerCell();
  double coeff=1./(double)sz;
  for(int i=0;i<nbOfCells;i++,ptToFill+=spaceDim)
    {
      std::fill(ptToFill,ptToFill+spaceDim,0.);
      for(int j=0;j<sz;j++,nodal++)
        if(*nodal>=0 && *nodal<nbOfNodes)
          std::transform(coor+spaceDim*nodal[0],coor+spaceDim*(nodal[0]+1),ptToFill,ptToFill,std::plus<double>());
        else
          {
            std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::computeIsoBarycenterOfNodesPerCell : on cell #" << i << " presence of nodeId #" << *nodal << " should be in [0," <<   nbOfNodes << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      std::transform(ptToFill,ptToFill+spaceDim,ptToFill,std::bind2nd(std::multiplies<double>(),coeff));
    }
  return ret.retn();
}

void MEDCoupling1SGTUMesh::renumberCells(const int *old2NewBg, bool check)
{
  int nbCells=getNumberOfCells();
  MCAuto<DataArrayInt> o2n=DataArrayInt::New();
  o2n->useArray(old2NewBg,false,C_DEALLOC,nbCells,1);
  if(check)
    o2n=o2n->checkAndPreparePermutation();
  //
  const int *conn=_conn->begin();
  MCAuto<DataArrayInt> n2o=o2n->invertArrayO2N2N2O(nbCells);
  const int *n2oPtr=n2o->begin();
  MCAuto<DataArrayInt> newConn=DataArrayInt::New();
  newConn->alloc(_conn->getNumberOfTuples(),1);
  newConn->copyStringInfoFrom(*_conn);
  int sz=getNumberOfNodesPerCell();
  //
  int *newC=newConn->getPointer();
  for(int i=0;i<nbCells;i++,newC+=sz)
    {
      int pos=n2oPtr[i];
      std::copy(conn+pos*sz,conn+(pos+1)*sz,newC);
    }
  _conn=newConn;
}

/*!
 * Keeps from \a this only cells which constituing point id are in the ids specified by [\a begin,\a end).
 * The resulting cell ids are stored at the end of the 'cellIdsKept' parameter.
 * Parameter \a fullyIn specifies if a cell that has part of its nodes in ids array is kept or not.
 * If \a fullyIn is true only cells whose ids are \b fully contained in [\a begin,\a end) tab will be kept.
 *
 * \param [in] begin input start of array of node ids.
 * \param [in] end input end of array of node ids.
 * \param [in] fullyIn input that specifies if all node ids must be in [\a begin,\a end) array to consider cell to be in.
 * \param [in,out] cellIdsKeptArr array where all candidate cell ids are put at the end.
 */
void MEDCoupling1SGTUMesh::fillCellIdsToKeepFromNodeIds(const int *begin, const int *end, bool fullyIn, DataArrayInt *&cellIdsKeptArr) const
{
  int nbOfCells=getNumberOfCells();
  MCAuto<DataArrayInt> cellIdsKept=DataArrayInt::New(); cellIdsKept->alloc(0,1);
  int tmp=-1;
  int sz=_conn->getMaxValue(tmp); sz=std::max(sz,0)+1;
  std::vector<bool> fastFinder(sz,false);
  for(const int *work=begin;work!=end;work++)
    if(*work>=0 && *work<sz)
      fastFinder[*work]=true;
  const int *conn=_conn->begin();
  int nbNodesPerCell=getNumberOfNodesPerCell();
  for(int i=0;i<nbOfCells;i++,conn+=nbNodesPerCell)
    {
      int ref=0,nbOfHit=0;
      for(int j=0;j<nbNodesPerCell;j++)
        if(conn[j]>=0)
          {
            ref++;
            if(fastFinder[conn[j]])
              nbOfHit++;
          }
      if((ref==nbOfHit && fullyIn) || (nbOfHit!=0 && !fullyIn))
        cellIdsKept->pushBackSilent(i);
    }
  cellIdsKeptArr=cellIdsKept.retn();
}

MEDCouplingMesh *MEDCoupling1SGTUMesh::mergeMyselfWith(const MEDCouplingMesh *other) const
{
  if(other->getType()!=SINGLE_STATIC_GEO_TYPE_UNSTRUCTURED)
    throw INTERP_KERNEL::Exception("Merge of umesh only available with umesh single static geo type each other !");
  const MEDCoupling1SGTUMesh *otherC=static_cast<const MEDCoupling1SGTUMesh *>(other);
  return Merge1SGTUMeshes(this,otherC);
}

MEDCouplingUMesh *MEDCoupling1SGTUMesh::buildUnstructured() const
{
  MCAuto<MEDCouplingUMesh> ret=MEDCouplingUMesh::New(getName(),getMeshDimension());
  ret->setCoords(getCoords());
  const int *nodalConn=_conn->begin();
  int nbCells=getNumberOfCells();
  int nbNodesPerCell=getNumberOfNodesPerCell();
  int geoType=(int)getCellModelEnum();
  MCAuto<DataArrayInt> c=DataArrayInt::New(); c->alloc(nbCells*(nbNodesPerCell+1),1);
  int *cPtr=c->getPointer();
  for(int i=0;i<nbCells;i++,nodalConn+=nbNodesPerCell)
    {
      *cPtr++=geoType;
      cPtr=std::copy(nodalConn,nodalConn+nbNodesPerCell,cPtr);
    }
  MCAuto<DataArrayInt> cI=DataArrayInt::Range(0,(nbCells+1)*(nbNodesPerCell+1),nbNodesPerCell+1);
  ret->setConnectivity(c,cI,true);
  try
  { ret->copyTinyInfoFrom(this); }
  catch(INTERP_KERNEL::Exception&) { }
  return ret.retn();
}

DataArrayInt *MEDCoupling1SGTUMesh::simplexize(int policy)
{
  switch(policy)
  {
    case 0:
      return simplexizePol0();
    case 1:
      return simplexizePol1();
    case (int) INTERP_KERNEL::PLANAR_FACE_5:
        return simplexizePlanarFace5();
    case (int) INTERP_KERNEL::PLANAR_FACE_6:
        return simplexizePlanarFace6();
    default:
      throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::simplexize : unrecognized policy ! Must be :\n  - 0 or 1 (only available for meshdim=2) \n  - PLANAR_FACE_5, PLANAR_FACE_6  (only for meshdim=3)");
  }
}

/// @cond INTERNAL

struct MEDCouplingAccVisit
{
  MEDCouplingAccVisit():_new_nb_of_nodes(0) { }
  int operator()(int val) { if(val!=-1) return _new_nb_of_nodes++; else return -1; }
  int _new_nb_of_nodes;
};

/// @endcond

/*!
 * This method returns all node ids used in \b this. The data array returned has to be dealt by the caller.
 * The returned node ids are sortes ascendingly. This method is closed to MEDCoupling1SGTUMesh::getNodeIdsInUse except
 * the format of returned DataArrayInt instance.
 *
 * \return a newly allocated DataArrayInt sorted ascendingly of fetched node ids.
 * \sa MEDCoupling1SGTUMesh::getNodeIdsInUse, areAllNodesFetched
 */
DataArrayInt *MEDCoupling1SGTUMesh::computeFetchedNodeIds() const
{
  checkConsistencyOfConnectivity();
  int nbNodes(getNumberOfNodes());
  std::vector<bool> fetchedNodes(nbNodes,false);
  computeNodeIdsAlg(fetchedNodes);
  int sz((int)std::count(fetchedNodes.begin(),fetchedNodes.end(),true));
  MCAuto<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(sz,1);
  int *retPtr(ret->getPointer());
  for(int i=0;i<nbNodes;i++)
    if(fetchedNodes[i])
      *retPtr++=i;
  return ret.retn();
}

/*!
 * Finds nodes not used in any cell and returns an array giving a new id to every node
 * by excluding the unused nodes, for which the array holds -1. The result array is
 * a mapping in "Old to New" mode. 
 *  \param [out] nbrOfNodesInUse - number of node ids present in the nodal connectivity.
 *  \return DataArrayInt * - a new instance of DataArrayInt. Its length is \a
 *          this->getNumberOfNodes(). It holds for each node of \a this mesh either -1
 *          if the node is unused or a new id else. The caller is to delete this
 *          array using decrRef() as it is no more needed.  
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If the nodal connectivity includes an invalid id.
 *  \sa MEDCoupling1SGTUMesh::computeFetchedNodeIds, areAllNodesFetched
 */
DataArrayInt *MEDCoupling1SGTUMesh::getNodeIdsInUse(int& nbrOfNodesInUse) const
{
  nbrOfNodesInUse=-1;
  int nbOfNodes=getNumberOfNodes();
  int nbOfCells=getNumberOfCells();
  MCAuto<DataArrayInt> ret(DataArrayInt::New());
  ret->alloc(nbOfNodes,1);
  int *traducer=ret->getPointer();
  std::fill(traducer,traducer+nbOfNodes,-1);
  const int *conn=_conn->begin();
  int nbNodesPerCell=getNumberOfNodesPerCell();
  for(int i=0;i<nbOfCells;i++)
    for(int j=0;j<nbNodesPerCell;j++,conn++)
      if(*conn>=0 && *conn<nbOfNodes)
        traducer[*conn]=1;
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::getNodeIdsInUse : In cell #" << i  << " presence of node id " <<  conn[j] << " not in [0," << nbOfNodes << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
  nbrOfNodesInUse=(int)std::count(traducer,traducer+nbOfNodes,1);
  std::transform(traducer,traducer+nbOfNodes,traducer,MEDCouplingAccVisit());
  return ret.retn();
}

/*!
 * This method renumbers only nodal connectivity in \a this. The renumbering is only an offset applied. So this method is a specialization of
 * \a renumberNodesInConn. \b WARNING, this method does not check that the resulting node ids in the nodal connectivity is in a valid range !
 *
 * \param [in] offset - specifies the offset to be applied on each element of connectivity.
 *
 * \sa renumberNodesInConn
 */
void MEDCoupling1SGTUMesh::renumberNodesWithOffsetInConn(int offset)
{
  getNumberOfCells();//only to check that all is well defined.
  _conn->applyLin(1,offset);
  updateTime();
}

/*!
 *  Same than renumberNodesInConn(const int *) except that here the format of old-to-new traducer is using map instead
 *  of array. This method is dedicated for renumbering from a big set of nodes the a tiny set of nodes which is the case during extraction
 *  of a big mesh.
 */
void MEDCoupling1SGTUMesh::renumberNodesInConn(const INTERP_KERNEL::HashMap<int,int>& newNodeNumbersO2N)
{
  getNumberOfCells();//only to check that all is well defined.
  int *begPtr(_conn->getPointer());
  int nbElt(_conn->getNumberOfTuples());
  int *endPtr(begPtr+nbElt);
  for(int *it=begPtr;it!=endPtr;it++)
    {
      INTERP_KERNEL::HashMap<int,int>::const_iterator it2(newNodeNumbersO2N.find(*it));
      if(it2!=newNodeNumbersO2N.end())
        {
          *it=(*it2).second;
        }
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::renumberNodesInConn : At pos #" << std::distance(begPtr,it) << " of nodal connectivity value is " << *it << ". Not in map !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  updateTime();
}

/*!
 * Changes ids of nodes within the nodal connectivity arrays according to a permutation
 * array in "Old to New" mode. The node coordinates array is \b not changed by this method.
 * This method is a generalization of shiftNodeNumbersInConn().
 *  \warning This method performs no check of validity of new ids. **Use it with care !**
 *  \param [in] newNodeNumbersO2N - a permutation array, of length \a
 *         this->getNumberOfNodes(), in "Old to New" mode. 
 *         See \ref numbering for more info on renumbering modes.
 *  \throw If the nodal connectivity of cells is not defined.
 */
void MEDCoupling1SGTUMesh::renumberNodesInConn(const int *newNodeNumbersO2N)
{
  getNumberOfCells();//only to check that all is well defined.
  _conn->transformWithIndArr(newNodeNumbersO2N,newNodeNumbersO2N+getNumberOfNodes());
  updateTime();
}

MEDCoupling1SGTUMesh *MEDCoupling1SGTUMesh::Merge1SGTUMeshes(const MEDCoupling1SGTUMesh *mesh1, const MEDCoupling1SGTUMesh *mesh2)
{
  std::vector<const MEDCoupling1SGTUMesh *> tmp(2);
  tmp[0]=const_cast<MEDCoupling1SGTUMesh *>(mesh1); tmp[1]=const_cast<MEDCoupling1SGTUMesh *>(mesh2);
  return Merge1SGTUMeshes(tmp);
}

MEDCoupling1SGTUMesh *MEDCoupling1SGTUMesh::Merge1SGTUMeshes(std::vector<const MEDCoupling1SGTUMesh *>& a)
{
  std::size_t sz=a.size();
  if(sz==0)
    return Merge1SGTUMeshesLL(a);
  for(std::size_t ii=0;ii<sz;ii++)
    if(!a[ii])
      {
        std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::Merge1SGTUMeshes : item #" << ii << " in input array of size "<< sz << " is empty !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  const INTERP_KERNEL::CellModel *cm=&(a[0]->getCellModel());
  for(std::size_t ii=0;ii<sz;ii++)
    if(&(a[ii]->getCellModel())!=cm)
      throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::Merge1SGTUMeshes : all items must have the same geo type !");
  std::vector< MCAuto<MEDCoupling1SGTUMesh> > bb(sz);
  std::vector< const MEDCoupling1SGTUMesh * > aa(sz);
  int spaceDim=-3;
  for(std::size_t i=0;i<sz && spaceDim==-3;i++)
    {
      const MEDCoupling1SGTUMesh *cur=a[i];
      const DataArrayDouble *coo=cur->getCoords();
      if(coo)
        spaceDim=coo->getNumberOfComponents();
    }
  if(spaceDim==-3)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::Merge1SGTUMeshes : no spaceDim specified ! unable to perform merge !");
  for(std::size_t i=0;i<sz;i++)
    {
      bb[i]=a[i]->buildSetInstanceFromThis(spaceDim);
      aa[i]=bb[i];
    }
  return Merge1SGTUMeshesLL(aa);
}

/*!
 * \throw If presence of a null instance in the input vector \a a.
 * \throw If a is empty
 */
MEDCoupling1SGTUMesh *MEDCoupling1SGTUMesh::Merge1SGTUMeshesOnSameCoords(std::vector<const MEDCoupling1SGTUMesh *>& a)
{
  if(a.empty())
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::Merge1SGTUMeshesOnSameCoords : input array must be NON EMPTY !");
  std::vector<const MEDCoupling1SGTUMesh *>::const_iterator it=a.begin();
  if(!(*it))
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::Merge1SGTUMeshesOnSameCoords : null instance in the first element of input vector !");
  std::vector<const DataArrayInt *> ncs(a.size());
  (*it)->getNumberOfCells();//to check that all is OK
  const DataArrayDouble *coords=(*it)->getCoords();
  const INTERP_KERNEL::CellModel *cm=&((*it)->getCellModel());
  ncs[0]=(*it)->getNodalConnectivity();
  it++;
  for(int i=1;it!=a.end();i++,it++)
    {
      if(!(*it))
        throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::Merge1SGTUMeshesOnSameCoords : presence of a null instance in the input vector !");
      if(cm!=&((*it)->getCellModel()))
        throw INTERP_KERNEL::Exception("Geometric types mismatches, Merge1SGTUMeshes impossible !");
      (*it)->getNumberOfCells();//to check that all is OK
      ncs[i]=(*it)->getNodalConnectivity();
      if(coords!=(*it)->getCoords())
        throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::Merge1SGTUMeshesOnSameCoords : not lying on same coords !");
    }
  MCAuto<MEDCoupling1SGTUMesh> ret(new MEDCoupling1SGTUMesh("merge",*cm));
  ret->setCoords(coords);
  ret->_conn=DataArrayInt::Aggregate(ncs);
  return ret.retn();
}

/*!
 * Assume that all instances in \a a are non null. If null it leads to a crash. That's why this method is assigned to be low level (LL)
 */
MEDCoupling1SGTUMesh *MEDCoupling1SGTUMesh::Merge1SGTUMeshesLL(std::vector<const MEDCoupling1SGTUMesh *>& a)
{
  if(a.empty())
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::Merge1SGTUMeshes : input array must be NON EMPTY !");
  std::vector<const MEDCoupling1SGTUMesh *>::const_iterator it=a.begin();
  int nbOfCells=(*it)->getNumberOfCells();
  const INTERP_KERNEL::CellModel *cm=&((*it)->getCellModel());
  int nbNodesPerCell=(*it)->getNumberOfNodesPerCell();
  it++;
  for(;it!=a.end();it++)
    {
      if(cm!=&((*it)->getCellModel()))
        throw INTERP_KERNEL::Exception("Geometric types mismatches, Merge1SGTUMeshes impossible !");
      nbOfCells+=(*it)->getNumberOfCells();
    }
  std::vector<const MEDCouplingPointSet *> aps(a.size());
  std::copy(a.begin(),a.end(),aps.begin());
  MCAuto<DataArrayDouble> pts=MergeNodesArray(aps);
  MCAuto<MEDCoupling1SGTUMesh> ret(new MEDCoupling1SGTUMesh("merge",*cm));
  ret->setCoords(pts);
  MCAuto<DataArrayInt> c=DataArrayInt::New();
  c->alloc(nbOfCells*nbNodesPerCell,1);
  int *cPtr=c->getPointer();
  int offset=0;
  for(it=a.begin();it!=a.end();it++)
    {
      int curConnLgth=(*it)->getNodalConnectivityLength();
      const int *curC=(*it)->_conn->begin();
      cPtr=std::transform(curC,curC+curConnLgth,cPtr,std::bind2nd(std::plus<int>(),offset));
      offset+=(*it)->getNumberOfNodes();
    }
  //
  ret->setNodalConnectivity(c);
  return ret.retn();
}

MEDCouplingPointSet *MEDCoupling1SGTUMesh::buildPartOfMySelfKeepCoords(const int *begin, const int *end) const
{
  int ncell=getNumberOfCells();
  MCAuto<MEDCoupling1SGTUMesh> ret(new MEDCoupling1SGTUMesh(getName(),*_cm));
  ret->setCoords(_coords);
  std::size_t nbOfElemsRet=std::distance(begin,end);
  const int *inConn=_conn->getConstPointer();
  int sz=getNumberOfNodesPerCell();
  MCAuto<DataArrayInt> connRet=DataArrayInt::New(); connRet->alloc((int)nbOfElemsRet*sz,1);
  int *connPtr=connRet->getPointer();
  for(const int *work=begin;work!=end;work++,connPtr+=sz)
    {
      if(*work>=0 && *work<ncell)
        std::copy(inConn+(work[0])*sz,inConn+(work[0]+1)*sz,connPtr);
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::buildPartOfMySelfKeepCoords : On pos #" << std::distance(begin,work) << " input cell id =" << *work << " should be in [0," << ncell << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  ret->_conn=connRet;
  ret->copyTinyInfoFrom(this);
  return ret.retn();
}

MEDCouplingPointSet *MEDCoupling1SGTUMesh::buildPartOfMySelfKeepCoordsSlice(int start, int end, int step) const
{
  int ncell=getNumberOfCells();
  int nbOfElemsRet=DataArray::GetNumberOfItemGivenBESRelative(start,end,step,"MEDCoupling1SGTUMesh::buildPartOfMySelfKeepCoordsSlice : ");
  MCAuto<MEDCoupling1SGTUMesh> ret(new MEDCoupling1SGTUMesh(getName(),*_cm));
  ret->setCoords(_coords);
  const int *inConn=_conn->getConstPointer();
  int sz=getNumberOfNodesPerCell();
  MCAuto<DataArrayInt> connRet=DataArrayInt::New(); connRet->alloc((int)nbOfElemsRet*sz,1);
  int *connPtr=connRet->getPointer();
  int curId=start;
  for(int i=0;i<nbOfElemsRet;i++,connPtr+=sz,curId+=step)
    {
      if(curId>=0 && curId<ncell)
        std::copy(inConn+curId*sz,inConn+(curId+1)*sz,connPtr);
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::buildPartOfMySelfKeepCoordsSlice : On pos #" << i << " input cell id =" << curId  << " should be in [0," << ncell << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  ret->_conn=connRet;
  ret->copyTinyInfoFrom(this);
  return ret.retn();
}

void MEDCoupling1SGTUMesh::computeNodeIdsAlg(std::vector<bool>& nodeIdsInUse) const
{
  int sz((int)nodeIdsInUse.size());
  for(const int *conn=_conn->begin();conn!=_conn->end();conn++)
    {
      if(*conn>=0 && *conn<sz)
       nodeIdsInUse[*conn]=true;
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::computeFetchedNodeIds : At pos #" << std::distance(_conn->begin(),conn) << " value is " << *conn << " must be in [0," << sz << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

MEDCoupling1SGTUMesh *MEDCoupling1SGTUMesh::buildSetInstanceFromThis(int spaceDim) const
{
  MCAuto<MEDCoupling1SGTUMesh> ret(new MEDCoupling1SGTUMesh(getName(),*_cm));
  MCAuto<DataArrayInt> tmp1;
  const DataArrayInt *nodalConn(_conn);
  if(!nodalConn)
    {
      tmp1=DataArrayInt::New(); tmp1->alloc(0,1);
    }
  else
    tmp1=_conn;
  ret->_conn=tmp1;
  if(!_coords)
    {
      MCAuto<DataArrayDouble> coords=DataArrayDouble::New(); coords->alloc(0,spaceDim);
      ret->setCoords(coords);
    }
  else
    ret->setCoords(_coords);
  return ret.retn();
}

DataArrayInt *MEDCoupling1SGTUMesh::simplexizePol0()
{
  int nbOfCells=getNumberOfCells();
  if(getCellModelEnum()!=INTERP_KERNEL::NORM_QUAD4)
    return DataArrayInt::Range(0,nbOfCells,1);
  MCAuto<DataArrayInt> newConn=DataArrayInt::New(); newConn->alloc(2*3*nbOfCells,1);
  MCAuto<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(2*nbOfCells,1);
  const int *c(_conn->begin());
  int *retPtr(ret->getPointer()),*newConnPtr(newConn->getPointer());
  for(int i=0;i<nbOfCells;i++,c+=4,newConnPtr+=6,retPtr+=2)
    {
      newConnPtr[0]=c[0]; newConnPtr[1]=c[1]; newConnPtr[2]=c[2];
      newConnPtr[3]=c[0]; newConnPtr[4]=c[2]; newConnPtr[5]=c[3];
      retPtr[0]=i; retPtr[1]=i;
    }
  _conn=newConn;
  _cm=&INTERP_KERNEL::CellModel::GetCellModel(INTERP_KERNEL::NORM_TRI3);
  updateTime();
  return ret.retn();
}

DataArrayInt *MEDCoupling1SGTUMesh::simplexizePol1()
{
  int nbOfCells=getNumberOfCells();
  if(getCellModelEnum()!=INTERP_KERNEL::NORM_QUAD4)
    return DataArrayInt::Range(0,nbOfCells,1);
  MCAuto<DataArrayInt> newConn=DataArrayInt::New(); newConn->alloc(2*3*nbOfCells,1);
  MCAuto<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(2*nbOfCells,1);
  const int *c(_conn->begin());
  int *retPtr(ret->getPointer()),*newConnPtr(newConn->getPointer());
  for(int i=0;i<nbOfCells;i++,c+=4,newConnPtr+=6,retPtr+=2)
    {
      newConnPtr[0]=c[0]; newConnPtr[1]=c[1]; newConnPtr[2]=c[3];
      newConnPtr[3]=c[1]; newConnPtr[4]=c[2]; newConnPtr[5]=c[3];
      retPtr[0]=i; retPtr[1]=i;
    }
  _conn=newConn;
  _cm=&INTERP_KERNEL::CellModel::GetCellModel(INTERP_KERNEL::NORM_TRI3);
  updateTime();
  return ret.retn();
}

DataArrayInt *MEDCoupling1SGTUMesh::simplexizePlanarFace5()
{
  int nbOfCells=getNumberOfCells();
  if(getCellModelEnum()!=INTERP_KERNEL::NORM_HEXA8)
    return DataArrayInt::Range(0,nbOfCells,1);
  MCAuto<DataArrayInt> newConn=DataArrayInt::New(); newConn->alloc(5*4*nbOfCells,1);
  MCAuto<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(5*nbOfCells,1);
  const int *c(_conn->begin());
  int *retPtr(ret->getPointer()),*newConnPtr(newConn->getPointer());
  for(int i=0;i<nbOfCells;i++,c+=8,newConnPtr+=20,retPtr+=5)
    {
      for(int j=0;j<20;j++)
        newConnPtr[j]=c[INTERP_KERNEL::SPLIT_NODES_5_WO[j]];
      retPtr[0]=i; retPtr[1]=i; retPtr[2]=i; retPtr[3]=i; retPtr[4]=i;
    }
  _conn=newConn;
  _cm=&INTERP_KERNEL::CellModel::GetCellModel(INTERP_KERNEL::NORM_TETRA4);
  updateTime();
  return ret.retn();
}

DataArrayInt *MEDCoupling1SGTUMesh::simplexizePlanarFace6()
{
  int nbOfCells=getNumberOfCells();
  if(getCellModelEnum()!=INTERP_KERNEL::NORM_HEXA8)
    return DataArrayInt::Range(0,nbOfCells,1);
  MCAuto<DataArrayInt> newConn=DataArrayInt::New(); newConn->alloc(6*4*nbOfCells,1);
  MCAuto<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(6*nbOfCells,1);
  const int *c(_conn->begin());
  int *retPtr(ret->getPointer()),*newConnPtr(newConn->getPointer());
  for(int i=0;i<nbOfCells;i++,c+=8,newConnPtr+=24,retPtr+=6)
    {
      for(int j=0;j<24;j++)
        newConnPtr[j]=c[INTERP_KERNEL::SPLIT_NODES_6_WO[j]];
      retPtr[0]=i; retPtr[1]=i; retPtr[2]=i; retPtr[3]=i; retPtr[4]=i; retPtr[5]=i;
    }
  _conn=newConn;
  _cm=&INTERP_KERNEL::CellModel::GetCellModel(INTERP_KERNEL::NORM_TETRA4);
  updateTime();
  return ret.retn();
}

void MEDCoupling1SGTUMesh::reprQuickOverview(std::ostream& stream) const
{
  stream << "MEDCoupling1SGTUMesh C++ instance at " << this << ". Type=" << _cm->getRepr() << ". Name : \"" << getName() << "\".";
  stream << " Mesh dimension : " << getMeshDimension() << ".";
  if(!_coords)
    { stream << " No coordinates set !"; return ; }
  if(!_coords->isAllocated())
    { stream << " Coordinates set but not allocated !"; return ; }
  stream << " Space dimension : " << _coords->getNumberOfComponents() << "." << std::endl;
  stream << "Number of nodes : " << _coords->getNumberOfTuples() << ".";
  if(!(const DataArrayInt *)_conn)
    { stream << std::endl << "Nodal connectivity NOT set !"; return ; }
  if(_conn->isAllocated())
    {
      if(_conn->getNumberOfComponents()==1)
        stream << std::endl << "Number of cells : " << getNumberOfCells() << ".";
    }
}

void MEDCoupling1SGTUMesh::checkFullyDefined() const
{
  if(!((const DataArrayInt *)_conn) || !((const DataArrayDouble *)_coords))
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::checkFullyDefined : part of this is not fully defined.");
}

/*!
 * First step of unserialization process.
 */
bool MEDCoupling1SGTUMesh::isEmptyMesh(const std::vector<int>& tinyInfo) const
{
  throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::isEmptyMesh : not implemented yet !");
}

void MEDCoupling1SGTUMesh::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  int it,order;
  double time=getTime(it,order);
  tinyInfo.clear(); tinyInfoD.clear(); littleStrings.clear();
  //
  littleStrings.push_back(getName());
  littleStrings.push_back(getDescription());
  littleStrings.push_back(getTimeUnit());
  //
  std::vector<std::string> littleStrings2,littleStrings3;
  if((const DataArrayDouble *)_coords)
    _coords->getTinySerializationStrInformation(littleStrings2);
  if((const DataArrayInt *)_conn)
    _conn->getTinySerializationStrInformation(littleStrings3);
  int sz0((int)littleStrings2.size()),sz1((int)littleStrings3.size());
  littleStrings.insert(littleStrings.end(),littleStrings2.begin(),littleStrings2.end());
  littleStrings.insert(littleStrings.end(),littleStrings3.begin(),littleStrings3.end());
  //
  tinyInfo.push_back(getCellModelEnum());
  tinyInfo.push_back(it);
  tinyInfo.push_back(order);
  std::vector<int> tinyInfo2,tinyInfo3;
  if((const DataArrayDouble *)_coords)
    _coords->getTinySerializationIntInformation(tinyInfo2);
  if((const DataArrayInt *)_conn)
    _conn->getTinySerializationIntInformation(tinyInfo3);
  int sz2((int)tinyInfo2.size()),sz3((int)tinyInfo3.size());
  tinyInfo.push_back(sz0); tinyInfo.push_back(sz1); tinyInfo.push_back(sz2); tinyInfo.push_back(sz3);
  tinyInfo.insert(tinyInfo.end(),tinyInfo2.begin(),tinyInfo2.end());
  tinyInfo.insert(tinyInfo.end(),tinyInfo3.begin(),tinyInfo3.end());
  //
  tinyInfoD.push_back(time);
}

void MEDCoupling1SGTUMesh::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const
{
  std::vector<int> tinyInfo2(tinyInfo.begin()+7,tinyInfo.begin()+7+tinyInfo[5]);
  std::vector<int> tinyInfo1(tinyInfo.begin()+7+tinyInfo[5],tinyInfo.begin()+7+tinyInfo[5]+tinyInfo[6]);
  a1->resizeForUnserialization(tinyInfo1);
  a2->resizeForUnserialization(tinyInfo2);
}

void MEDCoupling1SGTUMesh::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const
{
  int sz(0);
  if((const DataArrayInt *)_conn)
    if(_conn->isAllocated())
      sz=_conn->getNbOfElems();
  a1=DataArrayInt::New();
  a1->alloc(sz,1);
  if(sz!=0 && (const DataArrayInt *)_conn)
    std::copy(_conn->begin(),_conn->end(),a1->getPointer());
  sz=0;
  if((const DataArrayDouble *)_coords)
    if(_coords->isAllocated())
      sz=_coords->getNbOfElems();
  a2=DataArrayDouble::New();
  a2->alloc(sz,1);
  if(sz!=0 && (const DataArrayDouble *)_coords)
    std::copy(_coords->begin(),_coords->end(),a2->getPointer());
}

void MEDCoupling1SGTUMesh::unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                                           const std::vector<std::string>& littleStrings)
{
  INTERP_KERNEL::NormalizedCellType gt((INTERP_KERNEL::NormalizedCellType)tinyInfo[0]);
  _cm=&INTERP_KERNEL::CellModel::GetCellModel(gt);
  setName(littleStrings[0]);
  setDescription(littleStrings[1]);
  setTimeUnit(littleStrings[2]);
  setTime(tinyInfoD[0],tinyInfo[1],tinyInfo[2]);
  int sz0(tinyInfo[3]),sz1(tinyInfo[4]),sz2(tinyInfo[5]),sz3(tinyInfo[6]);
  //
  _coords=DataArrayDouble::New();
  std::vector<int> tinyInfo2(tinyInfo.begin()+7,tinyInfo.begin()+7+sz2);
  _coords->resizeForUnserialization(tinyInfo2);
  std::copy(a2->begin(),a2->end(),_coords->getPointer());
  _conn=DataArrayInt::New();
  std::vector<int> tinyInfo3(tinyInfo.begin()+7+sz2,tinyInfo.begin()+7+sz2+sz3);
  _conn->resizeForUnserialization(tinyInfo3);
  std::copy(a1->begin(),a1->end(),_conn->getPointer());
  std::vector<std::string> littleStrings2(littleStrings.begin()+3,littleStrings.begin()+3+sz0);
  _coords->finishUnserialization(tinyInfo2,littleStrings2);
  std::vector<std::string> littleStrings3(littleStrings.begin()+3+sz0,littleStrings.begin()+3+sz0+sz1);
  _conn->finishUnserialization(tinyInfo3,littleStrings3);
}

/*!
 * Checks if \a this and \a other meshes are geometrically equivalent with high
 * probability, else an exception is thrown. The meshes are considered equivalent if
 * (1) meshes contain the same number of nodes and the same number of elements of the
 * same types (2) three cells of the two meshes (first, last and middle) are based
 * on coincident nodes (with a specified precision).
 *  \param [in] other - the mesh to compare with.
 *  \param [in] prec - the precision used to compare nodes of the two meshes.
 *  \throw If the two meshes do not match.
 */
void MEDCoupling1SGTUMesh::checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const
{
  MEDCouplingPointSet::checkFastEquivalWith(other,prec);
  const MEDCoupling1SGTUMesh *otherC=dynamic_cast<const MEDCoupling1SGTUMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::checkFastEquivalWith : Two meshes are not unstructured with single static geometric type !");
  const DataArrayInt *c1(_conn),*c2(otherC->_conn);
  if(c1==c2)
    return;
  if(!c1 || !c2)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::checkFastEquivalWith : presence of nodal connectivity only in one of the 2 meshes !");
  if((c1->isAllocated() && !c2->isAllocated()) || (!c1->isAllocated() && c2->isAllocated()))
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::checkFastEquivalWith : in nodal connectivity, only one is allocated !");
  if(c1->getNumberOfComponents()!=1 || c1->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::checkFastEquivalWith : in nodal connectivity, must have 1 and only 1 component !");
  if(c1->getHashCode()!=c2->getHashCode())
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::checkFastEquivalWith : nodal connectivity differs");
}

MEDCouplingPointSet *MEDCoupling1SGTUMesh::mergeMyselfWithOnSameCoords(const MEDCouplingPointSet *other) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::mergeMyselfWithOnSameCoords : input other is null !");
  const MEDCoupling1SGTUMesh *otherC=dynamic_cast<const MEDCoupling1SGTUMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::mergeMyselfWithOnSameCoords : the input other mesh is not of type single statuc geo type unstructured !");
  std::vector<const MEDCoupling1SGTUMesh *> ms(2);
  ms[0]=this;
  ms[1]=otherC;
  return Merge1SGTUMeshesOnSameCoords(ms);
}

void MEDCoupling1SGTUMesh::getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const
{
  checkFullyDefined();
  int nbOfNodes=getNumberOfNodes();
  int *revNodalIndxPtr=(int *)malloc((nbOfNodes+1)*sizeof(int));
  revNodalIndx->useArray(revNodalIndxPtr,true,C_DEALLOC,nbOfNodes+1,1);
  std::fill(revNodalIndxPtr,revNodalIndxPtr+nbOfNodes+1,0);
  const int *conn=_conn->begin();
  int nbOfCells=getNumberOfCells();
  int nbOfEltsInRevNodal=0;
  int nbOfNodesPerCell=getNumberOfNodesPerCell();
  for(int eltId=0;eltId<nbOfCells;eltId++)
    {
      for(int j=0;j<nbOfNodesPerCell;j++,conn++)
        {
          if(conn[0]>=0 && conn[0]<nbOfNodes)
            {
              nbOfEltsInRevNodal++;
              revNodalIndxPtr[conn[0]+1]++;
            }
          else
            {
              std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::getReverseNodalConnectivity : At cell #" << eltId << " presence of nodeId #" << conn[0] << " should be in [0," << nbOfNodes << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
  std::transform(revNodalIndxPtr+1,revNodalIndxPtr+nbOfNodes+1,revNodalIndxPtr,revNodalIndxPtr+1,std::plus<int>());
  conn=_conn->begin();
  int *revNodalPtr=(int *)malloc((nbOfEltsInRevNodal)*sizeof(int));
  revNodal->useArray(revNodalPtr,true,C_DEALLOC,nbOfEltsInRevNodal,1);
  std::fill(revNodalPtr,revNodalPtr+nbOfEltsInRevNodal,-1);
  for(int eltId=0;eltId<nbOfCells;eltId++)
    {
      for(int j=0;j<nbOfNodesPerCell;j++,conn++)
        {
          *std::find_if(revNodalPtr+revNodalIndxPtr[*conn],revNodalPtr+revNodalIndxPtr[*conn+1],std::bind2nd(std::equal_to<int>(),-1))=eltId;
        }
    }
}

/*!
 * Use \a nodalConn array as nodal connectivity of \a this. The input \a nodalConn pointer can be null.
 */
void MEDCoupling1SGTUMesh::setNodalConnectivity(DataArrayInt *nodalConn)
{
  if(nodalConn)
    nodalConn->incrRef();
  _conn=nodalConn;
  declareAsNew();
}

/*!
 * \return DataArrayInt * - the internal reference to the nodal connectivity. The caller is not reponsible to deallocate it.
 */
DataArrayInt *MEDCoupling1SGTUMesh::getNodalConnectivity() const
{
  const DataArrayInt *ret(_conn);
  return const_cast<DataArrayInt *>(ret);
}

/*!
 * Allocates memory to store an estimation of the given number of cells. Closer is the estimation to the number of cells effectively inserted,
 * less will be the needs to realloc. If the number of cells to be inserted is not known simply put 0 to this parameter.
 * If a nodal connectivity previouly existed before the call of this method, it will be reset.
 *
 *  \param [in] nbOfCells - estimation of the number of cell \a this mesh will contain.
 */
void MEDCoupling1SGTUMesh::allocateCells(int nbOfCells)
{
  if(nbOfCells<0)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::allocateCells : the input number of cells should be >= 0 !");
  _conn=DataArrayInt::New();
  _conn->reserve(getNumberOfNodesPerCell()*nbOfCells);
  declareAsNew();
}

/*!
 * Appends at the end of \a this a cell having nodal connectivity array defined in [ \a nodalConnOfCellBg, \a nodalConnOfCellEnd ).
 *
 * \param [in] nodalConnOfCellBg - the begin (included) of nodal connectivity of the cell to add.
 * \param [in] nodalConnOfCellEnd - the end (excluded) of nodal connectivity of the cell to add.
 * \throw If the length of the input nodal connectivity array of the cell to add is not equal to number of nodes per cell relative to the unique geometric type
 *        attached to \a this.
 * \thow If the nodal connectivity array in \a this is null (call MEDCoupling1SGTUMesh::allocateCells before).
 */
void MEDCoupling1SGTUMesh::insertNextCell(const int *nodalConnOfCellBg, const int *nodalConnOfCellEnd)
{
  int sz=(int)std::distance(nodalConnOfCellBg,nodalConnOfCellEnd);
  int ref=getNumberOfNodesPerCell();
  if(sz==ref)
    {
      DataArrayInt *c(_conn);
      if(c)
        c->pushBackValsSilent(nodalConnOfCellBg,nodalConnOfCellEnd);
      else
        throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::insertNextCell : nodal connectivity array is null ! Call MEDCoupling1SGTUMesh::allocateCells before !");
    }
  else
    {
      std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::insertNextCell : input nodal size (" << sz << ") does not match number of nodes per cell of this (";
      oss << ref << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * This method builds the dual mesh of \a this and returns it.
 * 
 * \return MEDCoupling1SGTUMesh * - newly object created to be managed by the caller.
 * \throw If \a this is not a mesh containing only simplex cells.
 * \throw If \a this is not correctly allocated (coordinates and connectivities have to be correctly set !).
 * \throw If at least one node in \a this is orphan (without any simplex cell lying on it !)
 */
MEDCoupling1GTUMesh *MEDCoupling1SGTUMesh::computeDualMesh() const
{
  const INTERP_KERNEL::CellModel& cm(getCellModel());
  if(!cm.isSimplex())
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::computeDualMesh : this mesh is not a simplex mesh ! Please invoke simplexize of tetrahedrize on this before calling this method !");
  switch(getMeshDimension())
  {
    case 3:
      return computeDualMesh3D();
    case 2:
      return computeDualMesh2D();
    default:
      throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::computeDualMesh : meshdimension must be in [2,3] !");
  }
}

/*!
 * This method explode each NORM_HEXA8 cells in \a this into 6 NORM_QUAD4 cells and put the result into the MEDCoupling1SGTUMesh returned instance.
 * 
 * \return MEDCoupling1SGTUMesh * - a newly allocated instances (to be managed by the caller) storing the result of the explosion.
 * \throw If \a this is not a mesh containing only NORM_HEXA8 cells.
 * \throw If \a this is not properly allocated.
 */
MEDCoupling1SGTUMesh *MEDCoupling1SGTUMesh::explodeEachHexa8To6Quad4() const
{
  const INTERP_KERNEL::CellModel& cm(getCellModel());
  if(cm.getEnum()!=INTERP_KERNEL::NORM_HEXA8)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::explodeEachHexa8To6Quad4 : this method can be applied only on HEXA8 mesh !");
  int nbHexa8(getNumberOfCells());
  const int *inConnPtr(getNodalConnectivity()->begin());
  MCAuto<MEDCoupling1SGTUMesh> ret(MEDCoupling1SGTUMesh::New(getName(),INTERP_KERNEL::NORM_QUAD4));
  MCAuto<DataArrayInt> c(DataArrayInt::New()); c->alloc(nbHexa8*6*4,1);
  int *cPtr(c->getPointer());
  for(int i=0;i<nbHexa8;i++,inConnPtr+=8)
    {
      for(int j=0;j<6;j++,cPtr+=4)
        cm.fillSonCellNodalConnectivity(j,inConnPtr,cPtr);
    }
  ret->setCoords(getCoords());
  ret->setNodalConnectivity(c);
  return ret.retn();
}

/*!
 * This method starts from an unstructured mesh that hides in reality a cartesian mesh.
 * If it is not the case, an exception will be thrown.
 * This method returns three objects : The cartesian mesh geometrically equivalent to \a this (within a precision of \a eps) and a permutation of cells
 * and a permutation of nodes.
 *
 * - this[cellPerm[i]]=ret[i]
 *
 * \param [out] cellPerm the permutation array of size \c this->getNumberOfCells()
 * \param [out] nodePerm the permutation array of size \c this->getNumberOfNodes()
 * \return MEDCouplingCMesh * - a newly allocated mesh that is the result of the structurization of \a this.
 */
MEDCouplingCMesh *MEDCoupling1SGTUMesh::structurizeMe(DataArrayInt *& cellPerm, DataArrayInt *& nodePerm, double eps) const
{
  checkConsistencyLight();
  int spaceDim(getSpaceDimension()),meshDim(getMeshDimension()),nbNodes(getNumberOfNodes());
  if(MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension(meshDim)!=getCellModelEnum())
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::structurizeMe : the unique geo type in this is not compatible with the geometric type regarding mesh dimension !");
  MCAuto<MEDCouplingCMesh> cm(MEDCouplingCMesh::New());
  for(int i=0;i<spaceDim;i++)
    {
      std::vector<int> tmp(1,i);
      MCAuto<DataArrayDouble> elt(static_cast<DataArrayDouble*>(getCoords()->keepSelectedComponents(tmp)));
      elt=elt->getDifferentValues(eps);
      elt->sort(true);
      cm->setCoordsAt(i,elt);
    }
  if(nbNodes!=cm->getNumberOfNodes())
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::structurizeMe : considering the number of nodes after split per components in space this can't be a cartesian mesh ! Maybe your epsilon parameter is invalid ?");
  try
  { cm->copyTinyInfoFrom(this); }
  catch(INTERP_KERNEL::Exception&) { }
  MCAuto<MEDCouplingUMesh> um(cm->buildUnstructured()),self(buildUnstructured());
  self->checkGeoEquivalWith(um,12,eps,cellPerm,nodePerm);
  return cm.retn();
}

/// @cond INTERNAL

bool UpdateHexa8Cell(int validAxis, int neighId, const int *validConnQuad4NeighSide, int *allFacesNodalConn, int *myNeighbours)
{
  static const int TAB[48]={
    0,1,2,3,4,5,6,7,//0
    4,7,6,5,0,3,2,1,//1
    0,3,7,4,1,2,6,5,//2
    4,0,3,7,5,1,2,6,//3
    5,1,0,4,6,2,3,7,//4
    3,7,4,0,2,6,5,1 //5
  };
  static const int TAB2[6]={0,0,3,3,3,3};
  if(myNeighbours[validAxis]==neighId && allFacesNodalConn[4*validAxis+0]==validConnQuad4NeighSide[TAB2[validAxis]])
    return true;
  int oldAxis((int)std::distance(myNeighbours,std::find(myNeighbours,myNeighbours+6,neighId)));
  std::size_t pos(std::distance(MEDCoupling1SGTUMesh::HEXA8_FACE_PAIRS,std::find(MEDCoupling1SGTUMesh::HEXA8_FACE_PAIRS,MEDCoupling1SGTUMesh::HEXA8_FACE_PAIRS+6,oldAxis)));
  std::size_t pos0(pos/2),pos1(pos%2);
  int oldAxisOpp(MEDCoupling1SGTUMesh::HEXA8_FACE_PAIRS[2*pos0+(pos1+1)%2]);
  int oldConn[8],myConn2[8]={-1,-1,-1,-1,-1,-1,-1,-1},myConn[8],edgeConn[2],allFacesTmp[24],neighTmp[6];
  oldConn[0]=allFacesNodalConn[0]; oldConn[1]=allFacesNodalConn[1]; oldConn[2]=allFacesNodalConn[2]; oldConn[3]=allFacesNodalConn[3];
  oldConn[4]=allFacesNodalConn[4]; oldConn[5]=allFacesNodalConn[7]; oldConn[6]=allFacesNodalConn[6]; oldConn[7]=allFacesNodalConn[5];
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(INTERP_KERNEL::NORM_HEXA8));
  for(int i=0;i<4;i++)
    myConn2[i]=validConnQuad4NeighSide[(4-i+TAB2[validAxis])%4];
  for(int i=0;i<4;i++)
    {
      int nodeId(myConn2[i]);//the node id for which the opposite one will be found
      bool found(false);
      INTERP_KERNEL::NormalizedCellType typeOfSon;
      for(int j=0;j<12 && !found;j++)
        {
          cm.fillSonEdgesNodalConnectivity3D(j,oldConn,-1,edgeConn,typeOfSon);
          if(edgeConn[0]==nodeId || edgeConn[1]==nodeId)
            {
              if(std::find(allFacesNodalConn+4*oldAxisOpp,allFacesNodalConn+4*oldAxisOpp+4,edgeConn[0]==nodeId?edgeConn[1]:edgeConn[0])!=allFacesNodalConn+4*oldAxisOpp+4)
                {
                  myConn2[i+4]=edgeConn[0]==nodeId?edgeConn[1]:edgeConn[0];
                  found=true;
                }
            }
        }
      if(!found)
        throw INTERP_KERNEL::Exception("UpdateHexa8Cell : Internal Error !");
    }
  const int *myTab(TAB+8*validAxis);
  for(int i=0;i<8;i++)
    myConn[i]=myConn2[myTab[i]];
  for(int i=0;i<6;i++)
    {
      cm.fillSonCellNodalConnectivity(i,myConn,allFacesTmp+4*i);
      std::set<int> s(allFacesTmp+4*i,allFacesTmp+4*i+4);
      bool found(false);
      for(int j=0;j<6 && !found;j++)
        {
          std::set<int> s1(allFacesNodalConn+4*j,allFacesNodalConn+4*j+4);
          if(s==s1)
            {
              neighTmp[i]=myNeighbours[j];
              found=true;
            }
        }
      if(!found)
        throw INTERP_KERNEL::Exception("UpdateHexa8Cell : Internal Error #2 !");
    }
  std::copy(allFacesTmp,allFacesTmp+24,allFacesNodalConn);
  std::copy(neighTmp,neighTmp+6,myNeighbours);
  return false;
}

/// @endcond

/*!
 * This method expects the \a this contains NORM_HEXA8 cells only. This method will sort each cells in \a this so that their numbering was
 * homogeneous. If it succeeds the result of MEDCouplingUMesh::tetrahedrize will return a conform mesh.
 *
 * \return DataArrayInt * - a newly allocated array (to be managed by the caller) containing renumbered cell ids.
 *
 * \throw If \a this is not a mesh containing only NORM_HEXA8 cells.
 * \throw If \a this is not properly allocated.
 * \sa MEDCouplingUMesh::tetrahedrize, MEDCouplingUMesh::simplexize.
 */
DataArrayInt *MEDCoupling1SGTUMesh::sortHexa8EachOther()
{
  MCAuto<MEDCoupling1SGTUMesh> quads(explodeEachHexa8To6Quad4());//checks that only hexa8
  int nbHexa8(getNumberOfCells()),*cQuads(quads->getNodalConnectivity()->getPointer());
  MCAuto<DataArrayInt> neighOfQuads(DataArrayInt::New()); neighOfQuads->alloc(nbHexa8*6,1); neighOfQuads->fillWithValue(-1);
  int *ptNeigh(neighOfQuads->getPointer());
  {//neighOfQuads tells for each face of each Quad8 which cell (if!=-1) is connected to this face.
    MCAuto<MEDCouplingUMesh> quadsTmp(quads->buildUnstructured());
    MCAuto<DataArrayInt> ccSafe,cciSafe;
    DataArrayInt *cc(0),*cci(0);
    quadsTmp->findCommonCells(3,0,cc,cci);
    ccSafe=cc; cciSafe=cci;
    const int *ccPtr(ccSafe->begin()),nbOfPair(cci->getNumberOfTuples()-1);
    for(int i=0;i<nbOfPair;i++)
      { ptNeigh[ccPtr[2*i+0]]=ccPtr[2*i+1]/6; ptNeigh[ccPtr[2*i+1]]=ccPtr[2*i+0]/6; }
  }
  MCAuto<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
  std::vector<bool> fetched(nbHexa8,false);
  std::vector<bool>::iterator it(std::find(fetched.begin(),fetched.end(),false));
  while(it!=fetched.end())//it will turns as time as number of connected zones
    {
      int cellId((int)std::distance(fetched.begin(),it));//it is the seed of the connected zone.
      std::set<int> s; s.insert(cellId);//s contains already organized.
      while(!s.empty())
        {
          std::set<int> sNext;
          for(std::set<int>::const_iterator it0=s.begin();it0!=s.end();it0++)
            {
              fetched[*it0]=true;
              int *myNeighb(ptNeigh+6*(*it0));
              for(int i=0;i<6;i++)
                {
                  if(myNeighb[i]!=-1 && !fetched[myNeighb[i]])
                    {
                      std::size_t pos(std::distance(HEXA8_FACE_PAIRS,std::find(HEXA8_FACE_PAIRS,HEXA8_FACE_PAIRS+6,i)));
                      std::size_t pos0(pos/2),pos1(pos%2);
                      if(!UpdateHexa8Cell(HEXA8_FACE_PAIRS[2*pos0+(pos1+1)%2],*it0,cQuads+6*4*(*it0)+4*i,cQuads+6*4*myNeighb[i],ptNeigh+6*myNeighb[i]))
                        ret->pushBackSilent(myNeighb[i]);
                      fetched[myNeighb[i]]=true;
                      sNext.insert(myNeighb[i]);
                    }
                }
            }
          s=sNext;
        }
      it=std::find(fetched.begin(),fetched.end(),false);
    }
  if(!ret->empty())
    {
      int *conn(getNodalConnectivity()->getPointer());
      for(const int *pt=ret->begin();pt!=ret->end();pt++)
        {
          int cellId(*pt);
          conn[8*cellId+0]=cQuads[24*cellId+0]; conn[8*cellId+1]=cQuads[24*cellId+1]; conn[8*cellId+2]=cQuads[24*cellId+2]; conn[8*cellId+3]=cQuads[24*cellId+3];
          conn[8*cellId+4]=cQuads[24*cellId+4]; conn[8*cellId+5]=cQuads[24*cellId+7]; conn[8*cellId+6]=cQuads[24*cellId+6]; conn[8*cellId+7]=cQuads[24*cellId+5];
        }
      declareAsNew();
    }
  return ret.retn();
}

MEDCoupling1DGTUMesh *MEDCoupling1SGTUMesh::computeDualMesh3D() const
{
  static const int DUAL_TETRA_0[36]={
    4,1,0, 6,0,3, 7,3,1,
    4,0,1, 5,2,0, 8,1,2,
    6,3,0, 5,0,2, 9,2,3,
    7,1,3, 9,3,2, 8,2,1
  };
  static const int DUAL_TETRA_1[36]={
    8,4,10, 11,5,8, 10,7,11,
    9,4,8, 8,5,12, 12,6,9,
    10,4,9, 9,6,13, 13,7,10,
    12,5,11, 13,6,12, 11,7,13
  };
  static const int FACEID_NOT_SH_NODE[4]={2,3,1,0};
  if(getCellModelEnum()!=INTERP_KERNEL::NORM_TETRA4)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::computeDualMesh3D : only TETRA4 supported !");
  checkFullyDefined();
  MCAuto<MEDCouplingUMesh> thisu(buildUnstructured());
  MCAuto<DataArrayInt> revNodArr(DataArrayInt::New()),revNodIArr(DataArrayInt::New());
  thisu->getReverseNodalConnectivity(revNodArr,revNodIArr);
  const int *revNod(revNodArr->begin()),*revNodI(revNodIArr->begin()),*nodal(_conn->begin());
  MCAuto<DataArrayInt> d1Arr(DataArrayInt::New()),di1Arr(DataArrayInt::New()),rd1Arr(DataArrayInt::New()),rdi1Arr(DataArrayInt::New());
  MCAuto<MEDCouplingUMesh> edges(thisu->explode3DMeshTo1D(d1Arr,di1Arr,rd1Arr,rdi1Arr));
  const int *d1(d1Arr->begin());
  MCAuto<DataArrayInt> d2Arr(DataArrayInt::New()),di2Arr(DataArrayInt::New()),rd2Arr(DataArrayInt::New()),rdi2Arr(DataArrayInt::New());
  MCAuto<MEDCouplingUMesh> faces(thisu->buildDescendingConnectivity(d2Arr,di2Arr,rd2Arr,rdi2Arr));  thisu=0;
  const int *d2(d2Arr->begin()),*rdi2(rdi2Arr->begin());
  MCAuto<DataArrayDouble> edgesBaryArr(edges->computeCellCenterOfMass()),facesBaryArr(faces->computeCellCenterOfMass()),baryArr(computeCellCenterOfMass());
  const int nbOfNodes(getNumberOfNodes()),offset0(nbOfNodes+faces->getNumberOfCells()),offset1(offset0+edges->getNumberOfCells());
  edges=0; faces=0;
  std::vector<const DataArrayDouble *> v(4); v[0]=getCoords(); v[1]=facesBaryArr; v[2]=edgesBaryArr; v[3]=baryArr;
  MCAuto<DataArrayDouble> zeArr(DataArrayDouble::Aggregate(v)); baryArr=0; edgesBaryArr=0; facesBaryArr=0;
  std::string name("DualOf_"); name+=getName();
  MCAuto<MEDCoupling1DGTUMesh> ret(MEDCoupling1DGTUMesh::New(name,INTERP_KERNEL::NORM_POLYHED)); ret->setCoords(zeArr);
  MCAuto<DataArrayInt> cArr(DataArrayInt::New()),ciArr(DataArrayInt::New()); ciArr->alloc(nbOfNodes+1,1); ciArr->setIJ(0,0,0); cArr->alloc(0,1);
  for(int i=0;i<nbOfNodes;i++,revNodI++)
    {
      int nbOfCellsSharingNode(revNodI[1]-revNodI[0]);
      if(nbOfCellsSharingNode==0)
        {
          std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::computeDualMesh3D : Node #" << i << " is orphan !"; 
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      for(int j=0;j<nbOfCellsSharingNode;j++)
        {
          int curCellId(revNod[revNodI[0]+j]);
          const int *connOfCurCell(nodal+4*curCellId);
          std::size_t nodePosInCurCell(std::distance(connOfCurCell,std::find(connOfCurCell,connOfCurCell+4,i)));
          if(j!=0) cArr->pushBackSilent(-1);
          int tmp[14];
          //
          tmp[0]=d1[6*curCellId+DUAL_TETRA_0[nodePosInCurCell*9+0]-4]+offset0; tmp[1]=d2[4*curCellId+DUAL_TETRA_0[nodePosInCurCell*9+1]]+nbOfNodes;
          tmp[2]=curCellId+offset1; tmp[3]=d2[4*curCellId+DUAL_TETRA_0[nodePosInCurCell*9+2]]+nbOfNodes;
          tmp[4]=-1;
          tmp[5]=d1[6*curCellId+DUAL_TETRA_0[nodePosInCurCell*9+3]-4]+offset0; tmp[6]=d2[4*curCellId+DUAL_TETRA_0[nodePosInCurCell*9+4]]+nbOfNodes;
          tmp[7]=curCellId+offset1; tmp[8]=d2[4*curCellId+DUAL_TETRA_0[nodePosInCurCell*9+5]]+nbOfNodes;
          tmp[9]=-1;
          tmp[10]=d1[6*curCellId+DUAL_TETRA_0[nodePosInCurCell*9+6]-4]+offset0; tmp[11]=d2[4*curCellId+DUAL_TETRA_0[nodePosInCurCell*9+7]]+nbOfNodes;
          tmp[12]=curCellId+offset1; tmp[13]=d2[4*curCellId+DUAL_TETRA_0[nodePosInCurCell*9+8]]+nbOfNodes;
          cArr->insertAtTheEnd(tmp,tmp+14);
          int kk(0);
          for(int k=0;k<4;k++)
            {
              if(FACEID_NOT_SH_NODE[nodePosInCurCell]!=k)
                {
                  const int *faceId(d2+4*curCellId+k);
                  if(rdi2[*faceId+1]-rdi2[*faceId]==1)
                    {
                      int tmp2[5]; tmp2[0]=-1; tmp2[1]=i;
                      tmp2[2]=d1[6*curCellId+DUAL_TETRA_1[9*nodePosInCurCell+3*kk+0]-8]+offset0;
                      tmp2[3]=d2[4*curCellId+DUAL_TETRA_1[9*nodePosInCurCell+3*kk+1]-4]+nbOfNodes;
                      tmp2[4]=d1[6*curCellId+DUAL_TETRA_1[9*nodePosInCurCell+3*kk+2]-8]+offset0;
                      cArr->insertAtTheEnd(tmp2,tmp2+5);
                    }
                  kk++;
                }
            }
        }
      ciArr->setIJ(i+1,0,cArr->getNumberOfTuples());
    }
  ret->setNodalConnectivity(cArr,ciArr);
  return ret.retn();
}

MEDCoupling1DGTUMesh *MEDCoupling1SGTUMesh::computeDualMesh2D() const
{
  static const int DUAL_TRI_0[6]={0,2, 1,0, 2,1};
  static const int DUAL_TRI_1[6]={-3,+5, +3,-4, +4,-5};
  static const int FACEID_NOT_SH_NODE[3]={1,2,0};
  if(getCellModelEnum()!=INTERP_KERNEL::NORM_TRI3)
    throw INTERP_KERNEL::Exception("MEDCoupling1SGTUMesh::computeDualMesh2D : only TRI3 supported !");
  checkFullyDefined();
  MCAuto<MEDCouplingUMesh> thisu(buildUnstructured());
  MCAuto<DataArrayInt> revNodArr(DataArrayInt::New()),revNodIArr(DataArrayInt::New());
  thisu->getReverseNodalConnectivity(revNodArr,revNodIArr);
  const int *revNod(revNodArr->begin()),*revNodI(revNodIArr->begin()),*nodal(_conn->begin());
  MCAuto<DataArrayInt> d2Arr(DataArrayInt::New()),di2Arr(DataArrayInt::New()),rd2Arr(DataArrayInt::New()),rdi2Arr(DataArrayInt::New());
  MCAuto<MEDCouplingUMesh> edges(thisu->buildDescendingConnectivity(d2Arr,di2Arr,rd2Arr,rdi2Arr));  thisu=0;
  const int *d2(d2Arr->begin()),*rdi2(rdi2Arr->begin());
  MCAuto<DataArrayDouble> edgesBaryArr(edges->computeCellCenterOfMass()),baryArr(computeCellCenterOfMass());
  const int nbOfNodes(getNumberOfNodes()),offset0(nbOfNodes+edges->getNumberOfCells());
  edges=0;
  std::vector<const DataArrayDouble *> v(3); v[0]=getCoords(); v[1]=edgesBaryArr; v[2]=baryArr;
  MCAuto<DataArrayDouble> zeArr(DataArrayDouble::Aggregate(v)); baryArr=0; edgesBaryArr=0;
  std::string name("DualOf_"); name+=getName();
  MCAuto<MEDCoupling1DGTUMesh> ret(MEDCoupling1DGTUMesh::New(name,INTERP_KERNEL::NORM_POLYGON)); ret->setCoords(zeArr);
  MCAuto<DataArrayInt> cArr(DataArrayInt::New()),ciArr(DataArrayInt::New()); ciArr->alloc(nbOfNodes+1,1); ciArr->setIJ(0,0,0); cArr->alloc(0,1);
  for(int i=0;i<nbOfNodes;i++,revNodI++)
    {
      int nbOfCellsSharingNode(revNodI[1]-revNodI[0]);
      if(nbOfCellsSharingNode==0)
        {
          std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::computeDualMesh2D : Node #" << i << " is orphan !"; 
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      std::vector< std::vector<int> > polyg;
      for(int j=0;j<nbOfCellsSharingNode;j++)
        {
          int curCellId(revNod[revNodI[0]+j]);
          const int *connOfCurCell(nodal+3*curCellId);
          std::size_t nodePosInCurCell(std::distance(connOfCurCell,std::find(connOfCurCell,connOfCurCell+4,i)));
          std::vector<int> locV(3);
          locV[0]=d2[3*curCellId+DUAL_TRI_0[2*nodePosInCurCell+0]]+nbOfNodes; locV[1]=curCellId+offset0; locV[2]=d2[3*curCellId+DUAL_TRI_0[2*nodePosInCurCell+1]]+nbOfNodes;
          polyg.push_back(locV);
          int kk(0);
          for(int k=0;k<3;k++)
            {
              if(FACEID_NOT_SH_NODE[nodePosInCurCell]!=k)
                {
                  const int *edgeId(d2+3*curCellId+k);
                  if(rdi2[*edgeId+1]-rdi2[*edgeId]==1)
                    {
                      std::vector<int> locV2(2);
                      int zeLocEdgeIdRel(DUAL_TRI_1[2*nodePosInCurCell+kk]);
                      if(zeLocEdgeIdRel>0)
                        {  locV2[0]=d2[3*curCellId+zeLocEdgeIdRel-3]+nbOfNodes;  locV2[1]=i; }
                      else
                        {  locV2[0]=i; locV2[1]=d2[3*curCellId-zeLocEdgeIdRel-3]+nbOfNodes; }
                      polyg.push_back(locV2);
                    }
                  kk++;
                }
            }
        }
      std::vector<int> zePolyg(MEDCoupling1DGTUMesh::BuildAPolygonFromParts(polyg));
      cArr->insertAtTheEnd(zePolyg.begin(),zePolyg.end());
      ciArr->setIJ(i+1,0,cArr->getNumberOfTuples());
    }
  ret->setNodalConnectivity(cArr,ciArr);
  return ret.retn();
}

/*!
 * This method aggregate the bbox of each cell and put it into bbox 
 *
 * \param [in] arcDetEps - a parameter specifying in case of 2D quadratic polygon cell the detection limit between linear and arc circle. (By default 1e-12)
 *                         For all other cases this input parameter is ignored.
 * \return DataArrayDouble * - newly created object (to be managed by the caller) \a this number of cells tuples and 2*spacedim components.
 * 
 * \throw If \a this is not fully set (coordinates and connectivity).
 * \throw If a cell in \a this has no valid nodeId.
 */
DataArrayDouble *MEDCoupling1SGTUMesh::getBoundingBoxForBBTree(double arcDetEps) const
{
  int spaceDim(getSpaceDimension()),nbOfCells(getNumberOfCells()),nbOfNodes(getNumberOfNodes()),nbOfNodesPerCell(getNumberOfNodesPerCell());
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New()); ret->alloc(nbOfCells,2*spaceDim);
  double *bbox(ret->getPointer());
  for(int i=0;i<nbOfCells*spaceDim;i++)
    {
      bbox[2*i]=std::numeric_limits<double>::max();
      bbox[2*i+1]=-std::numeric_limits<double>::max();
    }
  const double *coordsPtr(_coords->getConstPointer());
  const int *conn(_conn->getConstPointer());
  for(int i=0;i<nbOfCells;i++)
    {
      int kk(0);
      for(int j=0;j<nbOfNodesPerCell;j++,conn++)
        {
          int nodeId(*conn);
          if(nodeId>=0 && nodeId<nbOfNodes)
            {
              for(int k=0;k<spaceDim;k++)
                {
                  bbox[2*spaceDim*i+2*k]=std::min(bbox[2*spaceDim*i+2*k],coordsPtr[spaceDim*nodeId+k]);
                  bbox[2*spaceDim*i+2*k+1]=std::max(bbox[2*spaceDim*i+2*k+1],coordsPtr[spaceDim*nodeId+k]);
                }
              kk++;
            }
        }
      if(kk==0)
        {
          std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::getBoundingBoxForBBTree : cell #" << i << " contains no valid nodeId !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret.retn();
}

/*!
 * Returns the cell field giving for each cell in \a this its diameter. Diameter means the max length of all possible SEG2 in the cell.
 *
 * \return a new instance of field containing the result. The returned instance has to be deallocated by the caller.
 */
MEDCouplingFieldDouble *MEDCoupling1SGTUMesh::computeDiameterField() const
{
  checkFullyDefined();
  MCAuto<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME));
  int nbCells(getNumberOfCells());
  MCAuto<DataArrayDouble> arr(DataArrayDouble::New());
  arr->alloc(nbCells,1);
  INTERP_KERNEL::AutoCppPtr<INTERP_KERNEL::DiameterCalculator> dc(_cm->buildInstanceOfDiameterCalulator(getSpaceDimension()));
  dc->computeFor1SGTUMeshFrmt(nbCells,_conn->begin(),getCoords()->begin(),arr->getPointer());
  ret->setMesh(this);
  ret->setArray(arr);
  ret->setName("Diameter");
  return ret.retn();
}

/*!
 * This method invert orientation of all cells in \a this. 
 * After calling this method the absolute value of measure of cells in \a this are the same than before calling.
 * This method only operates on the connectivity so coordinates are not touched at all.
 */
void MEDCoupling1SGTUMesh::invertOrientationOfAllCells()
{
  checkConsistencyOfConnectivity();
  INTERP_KERNEL::AutoCppPtr<INTERP_KERNEL::OrientationInverter> oi(INTERP_KERNEL::OrientationInverter::BuildInstanceFrom(getCellModelEnum()));
  int nbOfNodesPerCell((int)_cm->getNumberOfNodes()),nbCells(getNumberOfCells());
  int *conn(_conn->getPointer());
  for(int i=0;i<nbCells;i++)
    oi->operate(conn+i*nbOfNodesPerCell,conn+(i+1)*nbOfNodesPerCell);
  updateTime();
}

//== 

MEDCoupling1DGTUMesh *MEDCoupling1DGTUMesh::New()
{
  return new MEDCoupling1DGTUMesh;
}

MEDCoupling1DGTUMesh *MEDCoupling1DGTUMesh::New(const std::string& name, INTERP_KERNEL::NormalizedCellType type)
{
  if(type==INTERP_KERNEL::NORM_ERROR)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::New : NORM_ERROR is not a valid type to be used as base geometric type for a mesh !");
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  if(!cm.isDynamic())
    {
      std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::New : the input geometric type " << cm.getRepr() << " is static ! Only dynamic types are allowed here !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return new MEDCoupling1DGTUMesh(name,cm);
}

MEDCoupling1DGTUMesh::MEDCoupling1DGTUMesh()
{
}

MEDCoupling1DGTUMesh::MEDCoupling1DGTUMesh(const std::string& name, const INTERP_KERNEL::CellModel& cm):MEDCoupling1GTUMesh(name,cm)
{
}

MEDCoupling1DGTUMesh::MEDCoupling1DGTUMesh(const MEDCoupling1DGTUMesh& other, bool recDeepCpy):MEDCoupling1GTUMesh(other,recDeepCpy),_conn_indx(other._conn_indx),_conn(other._conn)
{
  if(recDeepCpy)
    {
      const DataArrayInt *c(other._conn);
      if(c)
        _conn=c->deepCopy();
      c=other._conn_indx;
      if(c)
        _conn_indx=c->deepCopy();
    }
}

MEDCoupling1DGTUMesh *MEDCoupling1DGTUMesh::clone(bool recDeepCpy) const
{
  return new MEDCoupling1DGTUMesh(*this,recDeepCpy);
}

/*!
 * This method behaves mostly like MEDCoupling1DGTUMesh::deepCopy method, except that only nodal connectivity arrays are deeply copied.
 * The coordinates are shared between \a this and the returned instance.
 * 
 * \return MEDCoupling1DGTUMesh * - A new object instance holding the copy of \a this (deep for connectivity, shallow for coordiantes)
 * \sa MEDCoupling1DGTUMesh::deepCopy
 */
MEDCoupling1DGTUMesh *MEDCoupling1DGTUMesh::deepCopyConnectivityOnly() const
{
  checkConsistencyLight();
  MCAuto<MEDCoupling1DGTUMesh> ret(clone(false));
  MCAuto<DataArrayInt> c(_conn->deepCopy()),ci(_conn_indx->deepCopy());
  ret->setNodalConnectivity(c,ci);
  return ret.retn();
}

void MEDCoupling1DGTUMesh::updateTime() const
{
  MEDCoupling1GTUMesh::updateTime();
  const DataArrayInt *c(_conn);
  if(c)
    updateTimeWith(*c);
  c=_conn_indx;
  if(c)
    updateTimeWith(*c);
}

std::size_t MEDCoupling1DGTUMesh::getHeapMemorySizeWithoutChildren() const
{
  return MEDCoupling1GTUMesh::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDCoupling1DGTUMesh::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDCoupling1GTUMesh::getDirectChildrenWithNull());
  ret.push_back((const DataArrayInt *)_conn);
  ret.push_back((const DataArrayInt *)_conn_indx);
  return ret;
}

MEDCoupling1DGTUMesh *MEDCoupling1DGTUMesh::deepCopy() const
{
  return clone(true);
}

bool MEDCoupling1DGTUMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::isEqualIfNotWhy : input other pointer is null !");
  std::ostringstream oss; oss.precision(15);
  const MEDCoupling1DGTUMesh *otherC=dynamic_cast<const MEDCoupling1DGTUMesh *>(other);
  if(!otherC)
    {
      reason="mesh given in input is not castable in MEDCoupling1DGTUMesh !";
      return false;
    }
  if(!MEDCoupling1GTUMesh::isEqualIfNotWhy(other,prec,reason))
    return false;
  const DataArrayInt *c1(_conn),*c2(otherC->_conn);
  if(c1==c2)
    return true;
  if(!c1 || !c2)
    {
      reason="in connectivity of single dynamic geometric type exactly one among this and other is null !";
      return false;
    }
  if(!c1->isEqualIfNotWhy(*c2,reason))
    {
      reason.insert(0,"Nodal connectivity DataArrayInt differs : ");
      return false;
    }
  c1=_conn_indx; c2=otherC->_conn_indx;
  if(c1==c2)
    return true;
  if(!c1 || !c2)
    {
      reason="in connectivity index of single dynamic geometric type exactly one among this and other is null !";
      return false;
    }
  if(!c1->isEqualIfNotWhy(*c2,reason))
    {
      reason.insert(0,"Nodal connectivity index DataArrayInt differs : ");
      return false;
    }
  return true;
}

bool MEDCoupling1DGTUMesh::isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::isEqualWithoutConsideringStr : input other pointer is null !");
  const MEDCoupling1DGTUMesh *otherC=dynamic_cast<const MEDCoupling1DGTUMesh *>(other);
  if(!otherC)
    return false;
  if(!MEDCoupling1GTUMesh::isEqualWithoutConsideringStr(other,prec))
    return false;
  const DataArrayInt *c1(_conn),*c2(otherC->_conn);
  if(c1==c2)
    return true;
  if(!c1 || !c2)
    return false;
  if(!c1->isEqualWithoutConsideringStr(*c2))
    return false;
  return true;
  c1=_conn_indx; c2=otherC->_conn_indx;
  if(c1==c2)
    return true;
  if(!c1 || !c2)
    return false;
  if(!c1->isEqualWithoutConsideringStr(*c2))
    return false;
  return true;
}

/*!
 * Checks if \a this and \a other meshes are geometrically equivalent with high
 * probability, else an exception is thrown. The meshes are considered equivalent if
 * (1) meshes contain the same number of nodes and the same number of elements of the
 * same types (2) three cells of the two meshes (first, last and middle) are based
 * on coincident nodes (with a specified precision).
 *  \param [in] other - the mesh to compare with.
 *  \param [in] prec - the precision used to compare nodes of the two meshes.
 *  \throw If the two meshes do not match.
 */
void MEDCoupling1DGTUMesh::checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const
{
  MEDCouplingPointSet::checkFastEquivalWith(other,prec);
  const MEDCoupling1DGTUMesh *otherC=dynamic_cast<const MEDCoupling1DGTUMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::checkFastEquivalWith : Two meshes are not unstructured with single dynamic geometric type !");
  const DataArrayInt *c1(_conn),*c2(otherC->_conn);
  if(c1!=c2)
    {
      if(!c1 || !c2)
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::checkFastEquivalWith : presence of nodal connectivity only in one of the 2 meshes !");
      if((c1->isAllocated() && !c2->isAllocated()) || (!c1->isAllocated() && c2->isAllocated()))
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::checkFastEquivalWith : in nodal connectivity, only one is allocated !");
      if(c1->getNumberOfComponents()!=1 || c1->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::checkFastEquivalWith : in nodal connectivity, must have 1 and only 1 component !");
      if(c1->getHashCode()!=c2->getHashCode())
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::checkFastEquivalWith : nodal connectivity differs");
    }
  c1=_conn_indx; c2=otherC->_conn_indx;
  if(c1!=c2)
    {
      if(!c1 || !c2)
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::checkFastEquivalWith : presence of nodal connectivity index only in one of the 2 meshes !");
      if((c1->isAllocated() && !c2->isAllocated()) || (!c1->isAllocated() && c2->isAllocated()))
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::checkFastEquivalWith : in nodal connectivity index, only one is allocated !");
      if(c1->getNumberOfComponents()!=1 || c1->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::checkFastEquivalWith : in nodal connectivity index, must have 1 and only 1 component !");
      if(c1->getHashCode()!=c2->getHashCode())
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::checkFastEquivalWith : nodal connectivity index differs");
    }
}

void MEDCoupling1DGTUMesh::checkConsistencyOfConnectivity() const
{
  const DataArrayInt *c1(_conn);
  if(c1)
    {
      if(c1->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("Nodal connectivity array is expected to be with number of components set to one !");
      if(c1->getInfoOnComponent(0)!="")
        throw INTERP_KERNEL::Exception("Nodal connectivity array is expected to have no info on its single component !");
      c1->checkAllocated();
    }
  else
    throw INTERP_KERNEL::Exception("Nodal connectivity array not defined !");
  //
  int sz2=_conn->getNumberOfTuples();
  c1=_conn_indx;
  if(c1)
    {
      if(c1->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("Nodal connectivity index array is expected to be with number of components set to one !");
      c1->checkAllocated();
      if(c1->getNumberOfTuples()<1)
        throw INTERP_KERNEL::Exception("Nodal connectivity index array is expected to have a a size of 1 at least !");
      if(c1->getInfoOnComponent(0)!="")
        throw INTERP_KERNEL::Exception("Nodal connectivity index array is expected to have no info on its single component !");
      int f=c1->front(),ll=c1->back();
      if(f<0 || (sz2>0 && f>=sz2))
        {
          std::ostringstream oss; oss << "Nodal connectivity index array first value (" << f << ") is expected to be exactly in [0," << sz2 << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(ll<0 || ll>sz2)
        {
          std::ostringstream oss; oss << "Nodal connectivity index array last value (" << ll << ") is expected to be exactly in [0," << sz2 << "] !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(f>ll)
        {
          std::ostringstream oss; oss << "Nodal connectivity index array looks very bad (not increasing monotonic) because front (" << f << ") is greater that back (" << ll << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  else
    throw INTERP_KERNEL::Exception("Nodal connectivity index array not defined !");
  int szOfC1Exp=_conn_indx->back();
  if(sz2<szOfC1Exp)
    {
      std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::checkConsistencyOfConnectivity : The expected length of nodal connectivity array regarding index is " << szOfC1Exp << " but the actual size of it is " << c1->getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * If \a this pass this method, you are sure that connectivity arrays are not null, with exactly one component, no name, no component name, allocated.
 * In addition you are sure that the length of nodal connectivity index array is bigger than or equal to one.
 * In addition you are also sure that length of nodal connectivity is coherent with the content of the last value in the index array.
 */
void MEDCoupling1DGTUMesh::checkConsistencyLight() const
{
  MEDCouplingPointSet::checkConsistencyLight();
  checkConsistencyOfConnectivity();
}

void MEDCoupling1DGTUMesh::checkConsistency(double eps) const
{
  checkConsistencyLight();
  const DataArrayInt *c1(_conn),*c2(_conn_indx);
  if(!c2->isMonotonic(true))
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::checkConsistency : the nodal connectivity index is expected to be increasing monotinic !");
  //
  int nbOfTuples=c1->getNumberOfTuples();
  int nbOfNodes=getNumberOfNodes();
  const int *w(c1->begin());
  for(int i=0;i<nbOfTuples;i++,w++)
    {
      if(*w==-1) continue;
      if(*w<0 || *w>=nbOfNodes)
        {
          std::ostringstream oss; oss << "At pos #" << i << " of nodal connectivity array references to node id #" << *w << " must be in [0," << nbOfNodes << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

std::size_t MEDCoupling1DGTUMesh::getNumberOfCells() const
{
  checkConsistencyOfConnectivity();//do not remove
  return _conn_indx->getNumberOfTuples()-1;
}

/*!
 * This method returns a newly allocated array containing this->getNumberOfCells() tuples and 1 component.
 * For each cell in \b this the number of nodes constituting cell is computed.
 * For each polyhedron cell, the sum of the number of nodes of each face constituting polyhedron cell is returned.
 * So for pohyhedrons some nodes can be counted several times in the returned result.
 * 
 * \return a newly allocated array
 */
DataArrayInt *MEDCoupling1DGTUMesh::computeNbOfNodesPerCell() const
{
  checkConsistencyLight();
  _conn_indx->checkMonotonic(true);
  if(getCellModelEnum()!=INTERP_KERNEL::NORM_POLYHED)
    return _conn_indx->deltaShiftIndex();
  // for polyhedrons
  int nbOfCells=_conn_indx->getNumberOfTuples()-1;
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbOfCells,1);
  int *retPtr=ret->getPointer();
  const int *ci=_conn_indx->begin(),*c=_conn->begin();
  for(int i=0;i<nbOfCells;i++,retPtr++,ci++)
    *retPtr=ci[1]-ci[0]-std::count(c+ci[0],c+ci[1],-1);
  return ret.retn();
}

/*!
 * This method returns a newly allocated array containing this->getNumberOfCells() tuples and 1 component.
 * For each cell in \b this the number of faces constituting (entity of dimension this->getMeshDimension()-1) cell is computed.
 * 
 * \return a newly allocated array
 */
DataArrayInt *MEDCoupling1DGTUMesh::computeNbOfFacesPerCell() const
{
  checkConsistencyLight();
  _conn_indx->checkMonotonic(true);
  if(getCellModelEnum()!=INTERP_KERNEL::NORM_POLYHED && getCellModelEnum()!=INTERP_KERNEL::NORM_QPOLYG)
    return _conn_indx->deltaShiftIndex();
  if(getCellModelEnum()==INTERP_KERNEL::NORM_QPOLYG)
    {
      MCAuto<DataArrayInt> ret=_conn_indx->deltaShiftIndex();
      ret->applyDivideBy(2);
      return ret.retn();
    }
  // for polyhedrons
  int nbOfCells=_conn_indx->getNumberOfTuples()-1;
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbOfCells,1);
  int *retPtr=ret->getPointer();
  const int *ci=_conn_indx->begin(),*c=_conn->begin();
  for(int i=0;i<nbOfCells;i++,retPtr++,ci++)
    *retPtr=std::count(c+ci[0],c+ci[1],-1)+1;
  return ret.retn();
}

/*!
 * This method computes effective number of nodes per cell. That is to say nodes appearing several times in nodal connectivity of a cell,
 * will be counted only once here whereas it will be counted several times in MEDCoupling1DGTUMesh::computeNbOfNodesPerCell method.
 *
 * \return DataArrayInt * - new object to be deallocated by the caller.
 * \sa MEDCoupling1DGTUMesh::computeNbOfNodesPerCell
 */
DataArrayInt *MEDCoupling1DGTUMesh::computeEffectiveNbOfNodesPerCell() const
{
  checkConsistencyLight();
  _conn_indx->checkMonotonic(true);
  int nbOfCells(_conn_indx->getNumberOfTuples()-1);
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbOfCells,1);
  int *retPtr(ret->getPointer());
  const int *ci(_conn_indx->begin()),*c(_conn->begin());
  if(getCellModelEnum()!=INTERP_KERNEL::NORM_POLYHED)
    {
      for(int i=0;i<nbOfCells;i++,retPtr++,ci++)
        {
          std::set<int> s(c+ci[0],c+ci[1]);
          *retPtr=(int)s.size();
        }
    }
  else
    {
      for(int i=0;i<nbOfCells;i++,retPtr++,ci++)
        {
          std::set<int> s(c+ci[0],c+ci[1]); s.erase(-1);
          *retPtr=(int)s.size();
        }
    }
  return ret.retn();
}

void MEDCoupling1DGTUMesh::getNodeIdsOfCell(std::size_t cellId, std::vector<int>& conn) const
{
  std::size_t nbOfCells(getNumberOfCells());//performs checks
  if(cellId<nbOfCells)
    {
      int strt=_conn_indx->getIJ(cellId,0),stp=_conn_indx->getIJ(cellId+1,0);
      int nbOfNodes=stp-strt;
      if(nbOfNodes<0)
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::getNodeIdsOfCell : the index array is invalid ! Should be increasing monotonic !");
      conn.resize(nbOfNodes);
      std::copy(_conn->begin()+strt,_conn->begin()+stp,conn.begin());
    }
  else
    {
      std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::getNodeIdsOfCell : request for cellId #" << cellId << " must be in [0," << nbOfCells << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

int MEDCoupling1DGTUMesh::getNumberOfNodesInCell(int cellId) const
{
  int nbOfCells(getNumberOfCells());//performs checks
  if(cellId>=0 && cellId<nbOfCells)
    {
      const int *conn(_conn->begin());
      int strt=_conn_indx->getIJ(cellId,0),stp=_conn_indx->getIJ(cellId+1,0);
      return stp-strt-std::count(conn+strt,conn+stp,-1);
    }
  else
    {
      std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::getNumberOfNodesInCell : request for cellId #" << cellId << " must be in [0," << nbOfCells << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

std::string MEDCoupling1DGTUMesh::simpleRepr() const
{
  static const char msg0[]="No coordinates specified !";
  std::ostringstream ret;
  ret << "Single dynamic geometic type (" << _cm->getRepr() << ") unstructured mesh with name : \"" << getName() << "\"\n";
  ret << "Description of mesh : \"" << getDescription() << "\"\n";
  int tmpp1,tmpp2;
  double tt=getTime(tmpp1,tmpp2);
  ret << "Time attached to the mesh [unit] : " << tt << " [" << getTimeUnit() << "]\n";
  ret << "Iteration : " << tmpp1  << " Order : " << tmpp2 << "\n";
  ret << "Mesh dimension : " << getMeshDimension() << "\nSpace dimension : ";
  if(_coords!=0)
    {
      const int spaceDim=getSpaceDimension();
      ret << spaceDim << "\nInfo attached on space dimension : ";
      for(int i=0;i<spaceDim;i++)
        ret << "\"" << _coords->getInfoOnComponent(i) << "\" ";
      ret << "\n";
    }
  else
    ret << msg0 << "\n";
  ret << "Number of nodes : ";
  if(_coords!=0)
    ret << getNumberOfNodes() << "\n";
  else
    ret << msg0 << "\n";
  ret << "Number of cells : ";
  bool isOK=true;
  try { checkConsistencyLight(); } catch(INTERP_KERNEL::Exception& /* e */)
  {
      ret << "Nodal connectivity arrays are not set or badly set !\n";
      isOK=false;
  }
  if(isOK)
    ret << getNumberOfCells() << "\n";
  ret << "Cell type : " << _cm->getRepr() << "\n";
  return ret.str();
}

std::string MEDCoupling1DGTUMesh::advancedRepr() const
{
  std::ostringstream ret;
  ret << simpleRepr();
  ret << "\nCoordinates array : \n___________________\n\n";
  if(_coords)
    _coords->reprWithoutNameStream(ret);
  else
    ret << "No array set !\n";
  ret << "\n\nNodal Connectivity : \n____________________\n\n";
  //
  bool isOK=true;
  try { checkConsistency(); } catch(INTERP_KERNEL::Exception& /* e */)
  {
      ret << "Nodal connectivity arrays are not set or badly set !\n";
      isOK=false;
  }
  if(!isOK)
    return ret.str();
  int nbOfCells=getNumberOfCells();
  const int *ci=_conn_indx->begin(),*c=_conn->begin();
  for(int i=0;i<nbOfCells;i++,ci++)
    {
      ret << "Cell #" << i << " : ";
      std::copy(c+ci[0],c+ci[1],std::ostream_iterator<int>(ret," "));
      ret << "\n";
    }
  return ret.str();
}

DataArrayDouble *MEDCoupling1DGTUMesh::computeIsoBarycenterOfNodesPerCell() const
{
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  int spaceDim=getSpaceDimension();
  int nbOfCells=getNumberOfCells();//checkConsistencyLight()
  int nbOfNodes=getNumberOfNodes();
  ret->alloc(nbOfCells,spaceDim);
  double *ptToFill=ret->getPointer();
  const double *coor=_coords->begin();
  const int *nodal=_conn->begin(),*nodali=_conn_indx->begin();
  nodal+=nodali[0];
  if(getCellModelEnum()!=INTERP_KERNEL::NORM_POLYHED)
    {
      for(int i=0;i<nbOfCells;i++,ptToFill+=spaceDim,nodali++)
        {
          std::fill(ptToFill,ptToFill+spaceDim,0.);
          if(nodali[0]<nodali[1])// >= to avoid division by 0.
            {
              for(int j=nodali[0];j<nodali[1];j++,nodal++)
                {
                  if(*nodal>=0 && *nodal<nbOfNodes)
                    std::transform(coor+spaceDim*nodal[0],coor+spaceDim*(nodal[0]+1),ptToFill,ptToFill,std::plus<double>());
                  else
                    {
                      std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::computeIsoBarycenterOfNodesPerCell : on cell #" << i << " presence of nodeId #" << *nodal << " should be in [0," <<   nbOfNodes << ") !";
                      throw INTERP_KERNEL::Exception(oss.str().c_str());
                    }
                  std::transform(ptToFill,ptToFill+spaceDim,ptToFill,std::bind2nd(std::multiplies<double>(),1./(nodali[1]-nodali[0])));
                }
            }
          else
            {
              std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::computeIsoBarycenterOfNodesPerCell : at cell #" << i << " the nodal index array is invalid !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
  else
    {
      for(int i=0;i<nbOfCells;i++,ptToFill+=spaceDim,nodali++)
        {
          std::fill(ptToFill,ptToFill+spaceDim,0.);
          if(nodali[0]<nodali[1])// >= to avoid division by 0.
            {
              int nbOfNod=0;
              for(int j=nodali[0];j<nodali[1];j++,nodal++)
                {
                  if(*nodal==-1) continue;
                  if(*nodal>=0 && *nodal<nbOfNodes)
                    {
                      std::transform(coor+spaceDim*nodal[0],coor+spaceDim*(nodal[0]+1),ptToFill,ptToFill,std::plus<double>());
                      nbOfNod++;
                    }
                  else
                    {
                      std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::computeIsoBarycenterOfNodesPerCell (polyhedron) : on cell #" << i << " presence of nodeId #" << *nodal << " should be in [0," <<   nbOfNodes << ") !";
                      throw INTERP_KERNEL::Exception(oss.str().c_str());
                    }
                }
              if(nbOfNod!=0)
                std::transform(ptToFill,ptToFill+spaceDim,ptToFill,std::bind2nd(std::multiplies<double>(),1./nbOfNod));
              else
                {
                  std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::computeIsoBarycenterOfNodesPerCell (polyhedron) : no nodes in cell #" << i << " !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
          else
            {
              std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::computeIsoBarycenterOfNodesPerCell (polyhedron)  : at cell #" << i << " the nodal index array is invalid !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
  return ret.retn();
}

void MEDCoupling1DGTUMesh::renumberCells(const int *old2NewBg, bool check)
{
  int nbCells=getNumberOfCells();
  MCAuto<DataArrayInt> o2n=DataArrayInt::New();
  o2n->useArray(old2NewBg,false,C_DEALLOC,nbCells,1);
  if(check)
    o2n=o2n->checkAndPreparePermutation();
  //
  const int *o2nPtr=o2n->getPointer();
  const int *conn=_conn->begin(),*conni=_conn_indx->begin();
  MCAuto<DataArrayInt> newConn=DataArrayInt::New();
  MCAuto<DataArrayInt> newConnI=DataArrayInt::New();
  newConn->alloc(_conn->getNumberOfTuples(),1); newConnI->alloc(nbCells,1);
  newConn->copyStringInfoFrom(*_conn); newConnI->copyStringInfoFrom(*_conn_indx);
  //
  int *newC=newConn->getPointer(),*newCI=newConnI->getPointer();
  for(int i=0;i<nbCells;i++)
    {
      int newPos=o2nPtr[i];
      int sz=conni[i+1]-conni[i];
      if(sz>=0)
        newCI[newPos]=sz;
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::renumberCells : the index nodal array is invalid for cell #" << i << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  newConnI->computeOffsetsFull(); newCI=newConnI->getPointer();
  //
  for(int i=0;i<nbCells;i++,conni++)
    {
      int newp=o2nPtr[i];
      std::copy(conn+conni[0],conn+conni[1],newC+newCI[newp]);
    }
  _conn=newConn;
  _conn_indx=newConnI;
}

MEDCouplingMesh *MEDCoupling1DGTUMesh::mergeMyselfWith(const MEDCouplingMesh *other) const
{
  if(other->getType()!=SINGLE_DYNAMIC_GEO_TYPE_UNSTRUCTURED)
    throw INTERP_KERNEL::Exception("Merge of umesh only available with umesh single dynamic geo type each other !");
  const MEDCoupling1DGTUMesh *otherC=static_cast<const MEDCoupling1DGTUMesh *>(other);
  return Merge1DGTUMeshes(this,otherC);
}

MEDCouplingUMesh *MEDCoupling1DGTUMesh::buildUnstructured() const
{
  MCAuto<MEDCouplingUMesh> ret=MEDCouplingUMesh::New(getName(),getMeshDimension());
  ret->setCoords(getCoords());
  const int *nodalConn=_conn->begin(),*nodalConnI=_conn_indx->begin();
  int nbCells=getNumberOfCells();//checkConsistencyLight
  int geoType=(int)getCellModelEnum();
  MCAuto<DataArrayInt> c=DataArrayInt::New(); c->alloc(nbCells+_conn->getNumberOfTuples(),1);
  MCAuto<DataArrayInt> cI=DataArrayInt::New(); cI->alloc(nbCells+1);
  int *cPtr=c->getPointer(),*ciPtr=cI->getPointer();
  ciPtr[0]=0;
  for(int i=0;i<nbCells;i++,ciPtr++)
    {
      int sz=nodalConnI[i+1]-nodalConnI[i];
      if(sz>=0)
        {
          *cPtr++=geoType;
          cPtr=std::copy(nodalConn+nodalConnI[i],nodalConn+nodalConnI[i+1],cPtr);
          ciPtr[1]=ciPtr[0]+sz+1;
        }
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::buildUnstructured : Invalid for nodal index for cell #" << i << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  ret->setConnectivity(c,cI,true);
  try
  { ret->copyTinyInfoFrom(this); }
  catch(INTERP_KERNEL::Exception&) { }
  return ret.retn();
}

/*!
 * Do nothing for the moment, because there is no policy that allows to split polygons, polyhedrons ... into simplexes
 */
DataArrayInt *MEDCoupling1DGTUMesh::simplexize(int policy)
{
  int nbOfCells=getNumberOfCells();
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbOfCells,1);
  ret->iota(0);
  return ret.retn();
}

void MEDCoupling1DGTUMesh::reprQuickOverview(std::ostream& stream) const
{
  stream << "MEDCoupling1DGTUMesh C++ instance at " << this << ". Type=" << _cm->getRepr() << ". Name : \"" << getName() << "\".";
  stream << " Mesh dimension : " << getMeshDimension() << ".";
  if(!_coords)
    { stream << " No coordinates set !"; return ; }
  if(!_coords->isAllocated())
    { stream << " Coordinates set but not allocated !"; return ; }
  stream << " Space dimension : " << _coords->getNumberOfComponents() << "." << std::endl;
  stream << "Number of nodes : " << _coords->getNumberOfTuples() << ".";
  bool isOK=true;
  try { checkConsistencyLight(); } catch(INTERP_KERNEL::Exception&  /* e */)
  {
      stream << std::endl << "Nodal connectivity NOT set properly !\n";
      isOK=false;
  }
  if(isOK)
    stream << std::endl << "Number of cells : " << getNumberOfCells() << ".";
}

void MEDCoupling1DGTUMesh::shallowCopyConnectivityFrom(const MEDCouplingPointSet *other)
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::shallowCopyConnectivityFrom : input pointer is null !");
  const MEDCoupling1DGTUMesh *otherC=dynamic_cast<const MEDCoupling1DGTUMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::shallowCopyConnectivityFrom : input pointer is not an MEDCoupling1DGTUMesh instance !");
  setNodalConnectivity(otherC->getNodalConnectivity(),otherC->getNodalConnectivityIndex());
}

MEDCouplingPointSet *MEDCoupling1DGTUMesh::mergeMyselfWithOnSameCoords(const MEDCouplingPointSet *other) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::mergeMyselfWithOnSameCoords : input other is null !");
  const MEDCoupling1DGTUMesh *otherC=dynamic_cast<const MEDCoupling1DGTUMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::mergeMyselfWithOnSameCoords : the input other mesh is not of type single statuc geo type unstructured !");
  std::vector<const MEDCoupling1DGTUMesh *> ms(2);
  ms[0]=this;
  ms[1]=otherC;
  return Merge1DGTUMeshesOnSameCoords(ms);
}

MEDCouplingPointSet *MEDCoupling1DGTUMesh::buildPartOfMySelfKeepCoords(const int *begin, const int *end) const
{
  checkConsistencyLight();
  MCAuto<MEDCoupling1DGTUMesh> ret(new MEDCoupling1DGTUMesh(getName(),*_cm));
  ret->setCoords(_coords);
  DataArrayInt *c=0,*ci=0;
  MEDCouplingUMesh::ExtractFromIndexedArrays(begin,end,_conn,_conn_indx,c,ci);
  MCAuto<DataArrayInt> cSafe(c),ciSafe(ci);
  ret->setNodalConnectivity(c,ci);
  return ret.retn();
}

MEDCouplingPointSet *MEDCoupling1DGTUMesh::buildPartOfMySelfKeepCoordsSlice(int start, int end, int step) const
{
  checkConsistencyLight();
  MCAuto<MEDCoupling1DGTUMesh> ret(new MEDCoupling1DGTUMesh(getName(),*_cm));
  ret->setCoords(_coords);
  DataArrayInt *c=0,*ci=0;
  MEDCouplingUMesh::ExtractFromIndexedArraysSlice(start,end,step,_conn,_conn_indx,c,ci);
  MCAuto<DataArrayInt> cSafe(c),ciSafe(ci);
  ret->setNodalConnectivity(c,ci);
  return ret.retn();
}

void MEDCoupling1DGTUMesh::computeNodeIdsAlg(std::vector<bool>& nodeIdsInUse) const
{
  checkConsistency();
  int sz((int)nodeIdsInUse.size());
  for(const int *conn=_conn->begin();conn!=_conn->end();conn++)
    {
      if(*conn>=0 && *conn<sz)
        nodeIdsInUse[*conn]=true;
      else
        {
          if(*conn!=-1)
            {
              std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::computeNodeIdsAlg : At pos #" << std::distance(_conn->begin(),conn) << " value is " << *conn << " must be in [0," << sz << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
}

void MEDCoupling1DGTUMesh::getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const
{
  checkFullyDefined();
  int nbOfNodes=getNumberOfNodes();
  int *revNodalIndxPtr=(int *)malloc((nbOfNodes+1)*sizeof(int));
  revNodalIndx->useArray(revNodalIndxPtr,true,C_DEALLOC,nbOfNodes+1,1);
  std::fill(revNodalIndxPtr,revNodalIndxPtr+nbOfNodes+1,0);
  const int *conn=_conn->begin(),*conni=_conn_indx->begin();
  int nbOfCells=getNumberOfCells();
  int nbOfEltsInRevNodal=0;
  for(int eltId=0;eltId<nbOfCells;eltId++)
    {
      int nbOfNodesPerCell=conni[eltId+1]-conni[eltId];
      if(nbOfNodesPerCell>=0)
        {
          for(int j=0;j<nbOfNodesPerCell;j++)
            {
              int nodeId=conn[conni[eltId]+j];
              if(nodeId==-1) continue;            
              if(nodeId>=0 && nodeId<nbOfNodes)
                {
                  nbOfEltsInRevNodal++;
                  revNodalIndxPtr[nodeId+1]++;
                }
              else
                {
                  std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::getReverseNodalConnectivity : At cell #" << eltId << " presence of nodeId #" << conn[0] << " should be in [0," << nbOfNodes << ") !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
        }
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::getReverseNodalConnectivity : At cell #" << eltId << "nodal connectivity is invalid !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  std::transform(revNodalIndxPtr+1,revNodalIndxPtr+nbOfNodes+1,revNodalIndxPtr,revNodalIndxPtr+1,std::plus<int>());
  conn=_conn->begin();
  int *revNodalPtr=(int *)malloc((nbOfEltsInRevNodal)*sizeof(int));
  revNodal->useArray(revNodalPtr,true,C_DEALLOC,nbOfEltsInRevNodal,1);
  std::fill(revNodalPtr,revNodalPtr+nbOfEltsInRevNodal,-1);
  for(int eltId=0;eltId<nbOfCells;eltId++)
    {
      int nbOfNodesPerCell=conni[eltId+1]-conni[eltId];
      for(int j=0;j<nbOfNodesPerCell;j++)
        {
          int nodeId=conn[conni[eltId]+j];
          if(nodeId!=-1)
            *std::find_if(revNodalPtr+revNodalIndxPtr[nodeId],revNodalPtr+revNodalIndxPtr[nodeId+1],std::bind2nd(std::equal_to<int>(),-1))=eltId;
        }
    }
}

void MEDCoupling1DGTUMesh::checkFullyDefined() const
{
  if(!((const DataArrayInt *)_conn) || !((const DataArrayInt *)_conn_indx) || !((const DataArrayDouble *)_coords))
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::checkFullyDefined : part of this is not fully defined.");
}

bool MEDCoupling1DGTUMesh::isEmptyMesh(const std::vector<int>& tinyInfo) const
{
  throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::isEmptyMesh : not implemented yet !");
}

void MEDCoupling1DGTUMesh::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  int it,order;
  double time=getTime(it,order);
  tinyInfo.clear(); tinyInfoD.clear(); littleStrings.clear();
  //
  littleStrings.push_back(getName());
  littleStrings.push_back(getDescription());
  littleStrings.push_back(getTimeUnit());
  //
  std::vector<std::string> littleStrings2,littleStrings3,littleStrings4;
  if((const DataArrayDouble *)_coords)
    _coords->getTinySerializationStrInformation(littleStrings2);
  if((const DataArrayInt *)_conn)
    _conn->getTinySerializationStrInformation(littleStrings3);
  if((const DataArrayInt *)_conn_indx)
    _conn_indx->getTinySerializationStrInformation(littleStrings4);
  int sz0((int)littleStrings2.size()),sz1((int)littleStrings3.size()),sz2((int)littleStrings4.size());
  littleStrings.insert(littleStrings.end(),littleStrings2.begin(),littleStrings2.end());
  littleStrings.insert(littleStrings.end(),littleStrings3.begin(),littleStrings3.end());
  littleStrings.insert(littleStrings.end(),littleStrings4.begin(),littleStrings4.end());
  //
  tinyInfo.push_back(getCellModelEnum());
  tinyInfo.push_back(it);
  tinyInfo.push_back(order);
  std::vector<int> tinyInfo2,tinyInfo3,tinyInfo4;
  if((const DataArrayDouble *)_coords)
    _coords->getTinySerializationIntInformation(tinyInfo2);
  if((const DataArrayInt *)_conn)
    _conn->getTinySerializationIntInformation(tinyInfo3);
  if((const DataArrayInt *)_conn_indx)
    _conn_indx->getTinySerializationIntInformation(tinyInfo4);
  int sz3((int)tinyInfo2.size()),sz4((int)tinyInfo3.size()),sz5((int)tinyInfo4.size());
  tinyInfo.push_back(sz0); tinyInfo.push_back(sz1); tinyInfo.push_back(sz2); tinyInfo.push_back(sz3); tinyInfo.push_back(sz4);  tinyInfo.push_back(sz5);
  tinyInfo.insert(tinyInfo.end(),tinyInfo2.begin(),tinyInfo2.end());
  tinyInfo.insert(tinyInfo.end(),tinyInfo3.begin(),tinyInfo3.end());
  tinyInfo.insert(tinyInfo.end(),tinyInfo4.begin(),tinyInfo4.end());
  //
  tinyInfoD.push_back(time);
}

void MEDCoupling1DGTUMesh::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const
{
  std::vector<int> tinyInfo2(tinyInfo.begin()+9,tinyInfo.begin()+9+tinyInfo[6]);
  std::vector<int> tinyInfo1(tinyInfo.begin()+9+tinyInfo[6],tinyInfo.begin()+9+tinyInfo[6]+tinyInfo[7]);
  std::vector<int> tinyInfo12(tinyInfo.begin()+9+tinyInfo[6]+tinyInfo[7],tinyInfo.begin()+9+tinyInfo[6]+tinyInfo[7]+tinyInfo[8]);
  MCAuto<DataArrayInt> p1(DataArrayInt::New()); p1->resizeForUnserialization(tinyInfo1);
  MCAuto<DataArrayInt> p2(DataArrayInt::New()); p2->resizeForUnserialization(tinyInfo12);
  std::vector<const DataArrayInt *> v(2); v[0]=p1; v[1]=p2;
  p2=DataArrayInt::Aggregate(v);
  a2->resizeForUnserialization(tinyInfo2);
  a1->alloc(p2->getNbOfElems(),1);
}

void MEDCoupling1DGTUMesh::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const
{
  int sz(0);
  if((const DataArrayInt *)_conn)
    if(_conn->isAllocated())
      sz=_conn->getNbOfElems();
  if((const DataArrayInt *)_conn_indx)
    if(_conn_indx->isAllocated())
      sz+=_conn_indx->getNbOfElems();
  a1=DataArrayInt::New();
  a1->alloc(sz,1);
  int *work(a1->getPointer());
  if(sz!=0 && (const DataArrayInt *)_conn)
    work=std::copy(_conn->begin(),_conn->end(),a1->getPointer());
  if(sz!=0 && (const DataArrayInt *)_conn_indx)
    std::copy(_conn_indx->begin(),_conn_indx->end(),work);
  sz=0;
  if((const DataArrayDouble *)_coords)
    if(_coords->isAllocated())
      sz=_coords->getNbOfElems();
  a2=DataArrayDouble::New();
  a2->alloc(sz,1);
  if(sz!=0 && (const DataArrayDouble *)_coords)
    std::copy(_coords->begin(),_coords->end(),a2->getPointer());
}

void MEDCoupling1DGTUMesh::unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                                           const std::vector<std::string>& littleStrings)
{
  INTERP_KERNEL::NormalizedCellType gt((INTERP_KERNEL::NormalizedCellType)tinyInfo[0]);
  _cm=&INTERP_KERNEL::CellModel::GetCellModel(gt);
  setName(littleStrings[0]);
  setDescription(littleStrings[1]);
  setTimeUnit(littleStrings[2]);
  setTime(tinyInfoD[0],tinyInfo[1],tinyInfo[2]);
  int sz0(tinyInfo[3]),sz1(tinyInfo[4]),sz2(tinyInfo[5]),sz3(tinyInfo[6]),sz4(tinyInfo[7]),sz5(tinyInfo[8]);
  //
  _coords=DataArrayDouble::New();
  std::vector<int> tinyInfo2(tinyInfo.begin()+9,tinyInfo.begin()+9+sz3);
  _coords->resizeForUnserialization(tinyInfo2);
  std::copy(a2->begin(),a2->end(),_coords->getPointer());
  _conn=DataArrayInt::New();
  std::vector<int> tinyInfo3(tinyInfo.begin()+9+sz3,tinyInfo.begin()+9+sz3+sz4);
  _conn->resizeForUnserialization(tinyInfo3);
  std::copy(a1->begin(),a1->begin()+_conn->getNbOfElems(),_conn->getPointer());
  _conn_indx=DataArrayInt::New();
  std::vector<int> tinyInfo4(tinyInfo.begin()+9+sz3+sz4,tinyInfo.begin()+9+sz3+sz4+sz5);
  _conn_indx->resizeForUnserialization(tinyInfo4);
  std::copy(a1->begin()+_conn->getNbOfElems(),a1->end(),_conn_indx->getPointer());
  std::vector<std::string> littleStrings2(littleStrings.begin()+3,littleStrings.begin()+3+sz0);
  _coords->finishUnserialization(tinyInfo2,littleStrings2);
  std::vector<std::string> littleStrings3(littleStrings.begin()+3+sz0,littleStrings.begin()+3+sz0+sz1);
  _conn->finishUnserialization(tinyInfo3,littleStrings3);
  std::vector<std::string> littleStrings4(littleStrings.begin()+3+sz0+sz1,littleStrings.begin()+3+sz0+sz1+sz2);
  _conn_indx->finishUnserialization(tinyInfo4,littleStrings4);
}

/*!
 * Finds nodes not used in any cell and returns an array giving a new id to every node
 * by excluding the unused nodes, for which the array holds -1. The result array is
 * a mapping in "Old to New" mode.
 *  \param [out] nbrOfNodesInUse - number of node ids present in the nodal connectivity.
 *  \return DataArrayInt * - a new instance of DataArrayInt. Its length is \a
 *          this->getNumberOfNodes(). It holds for each node of \a this mesh either -1
 *          if the node is unused or a new id else. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If the nodal connectivity includes an invalid id.
 *  \sa MEDCoupling1DGTUMesh::getNodeIdsInUse, areAllNodesFetched
 */
DataArrayInt *MEDCoupling1DGTUMesh::computeFetchedNodeIds() const
{
  checkConsistency();
  int nbNodes(getNumberOfNodes());
  std::vector<bool> fetchedNodes(nbNodes,false);
  computeNodeIdsAlg(fetchedNodes);
  int sz((int)std::count(fetchedNodes.begin(),fetchedNodes.end(),true));
  MCAuto<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(sz,1);
  int *retPtr(ret->getPointer());
  for(int i=0;i<nbNodes;i++)
    if(fetchedNodes[i])
      *retPtr++=i;
  return ret.retn();
}

/*!
 * Finds nodes not used in any cell and returns an array giving a new id to every node
 * by excluding the unused nodes, for which the array holds -1. The result array is
 * a mapping in "Old to New" mode. 
 *  \param [out] nbrOfNodesInUse - number of node ids present in the nodal connectivity.
 *  \return DataArrayInt * - a new instance of DataArrayInt. Its length is \a
 *          this->getNumberOfNodes(). It holds for each node of \a this mesh either -1
 *          if the node is unused or a new id else. The caller is to delete this
 *          array using decrRef() as it is no more needed.  
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If the nodal connectivity includes an invalid id.
 *  \sa MEDCoupling1DGTUMesh::computeFetchedNodeIds, areAllNodesFetched
 */
DataArrayInt *MEDCoupling1DGTUMesh::getNodeIdsInUse(int& nbrOfNodesInUse) const
{
  nbrOfNodesInUse=-1;
  int nbOfNodes=getNumberOfNodes();
  int nbOfCells=getNumberOfCells();//checkConsistencyLight
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbOfNodes,1);
  int *traducer=ret->getPointer();
  std::fill(traducer,traducer+nbOfNodes,-1);
  const int *conn=_conn->begin(),*conni(_conn_indx->begin());
  for(int i=0;i<nbOfCells;i++,conni++)
    {
      int nbNodesPerCell=conni[1]-conni[0];
      for(int j=0;j<nbNodesPerCell;j++)
        {
          int nodeId=conn[conni[0]+j];
          if(nodeId==-1) continue;
          if(nodeId>=0 && nodeId<nbOfNodes)
            traducer[nodeId]=1;
          else
            {
              std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::getNodeIdsInUse : In cell #" << i  << " presence of node id " <<  nodeId << " not in [0," << nbOfNodes << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
  nbrOfNodesInUse=(int)std::count(traducer,traducer+nbOfNodes,1);
  std::transform(traducer,traducer+nbOfNodes,traducer,MEDCouplingAccVisit());
  return ret.retn();
}

/*!
 * This method renumbers only nodal connectivity in \a this. The renumbering is only an offset applied. So this method is a specialization of
 * \a renumberNodesInConn. \b WARNING, this method does not check that the resulting node ids in the nodal connectivity is in a valid range !
 *
 * \param [in] offset - specifies the offset to be applied on each element of connectivity.
 *
 * \sa renumberNodesInConn
 */
void MEDCoupling1DGTUMesh::renumberNodesWithOffsetInConn(int offset)
{
  getNumberOfCells();//only to check that all is well defined.
  //
  int nbOfTuples(_conn->getNumberOfTuples());
  int *pt(_conn->getPointer());
  for(int i=0;i<nbOfTuples;i++,pt++)
    {
      if(*pt==-1) continue;
      *pt+=offset;
    }
  //
  updateTime();
}

/*!
 *  Same than renumberNodesInConn(const int *) except that here the format of old-to-new traducer is using map instead
 *  of array. This method is dedicated for renumbering from a big set of nodes the a tiny set of nodes which is the case during extraction
 *  of a big mesh.
 */
void MEDCoupling1DGTUMesh::renumberNodesInConn(const INTERP_KERNEL::HashMap<int,int>& newNodeNumbersO2N)
{
  getNumberOfCells();//only to check that all is well defined.
  //
  int nbElemsIn(getNumberOfNodes()),nbOfTuples(_conn->getNumberOfTuples());
  int *pt(_conn->getPointer());
  for(int i=0;i<nbOfTuples;i++,pt++)
    {
      if(*pt==-1) continue;
      if(*pt>=0 && *pt<nbElemsIn)
        {
          INTERP_KERNEL::HashMap<int,int>::const_iterator it(newNodeNumbersO2N.find(*pt));
          if(it!=newNodeNumbersO2N.end())
            *pt=(*it).second;
          else
            {
              std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::renumberNodesInConn : At pos #" << i << " of connectivity, node id is " << *pt << ". Not in keys of input map !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::renumberNodesInConn : error on tuple #" << i << " value is " << *pt << " and indirectionnal array as a size equal to " << nbElemsIn;
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  //
  updateTime();
}

/*!
 * Changes ids of nodes within the nodal connectivity arrays according to a permutation
 * array in "Old to New" mode. The node coordinates array is \b not changed by this method.
 * This method is a generalization of shiftNodeNumbersInConn().
 *  \warning This method performs no check of validity of new ids. **Use it with care !**
 *  \param [in] newNodeNumbersO2N - a permutation array, of length \a
 *         this->getNumberOfNodes(), in "Old to New" mode. 
 *         See \ref numbering for more info on renumbering modes.
 *  \throw If the nodal connectivity of cells is not defined.
 */
void MEDCoupling1DGTUMesh::renumberNodesInConn(const int *newNodeNumbersO2N)
{
  getNumberOfCells();//only to check that all is well defined.
  //
  int nbElemsIn(getNumberOfNodes()),nbOfTuples(_conn->getNumberOfTuples());
  int *pt(_conn->getPointer());
  for(int i=0;i<nbOfTuples;i++,pt++)
    {
      if(*pt==-1) continue;
      if(*pt>=0 && *pt<nbElemsIn)
        *pt=newNodeNumbersO2N[*pt];
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::renumberNodesInConn : error on tuple #" << i << " value is " << *pt << " and indirectionnal array as a size equal to " << nbElemsIn;
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  //
  updateTime();
}

/*!
 * Keeps from \a this only cells which constituing point id are in the ids specified by [\a begin,\a end).
 * The resulting cell ids are stored at the end of the 'cellIdsKept' parameter.
 * Parameter \a fullyIn specifies if a cell that has part of its nodes in ids array is kept or not.
 * If \a fullyIn is true only cells whose ids are \b fully contained in [\a begin,\a end) tab will be kept.
 *
 * \param [in] begin input start of array of node ids.
 * \param [in] end input end of array of node ids.
 * \param [in] fullyIn input that specifies if all node ids must be in [\a begin,\a end) array to consider cell to be in.
 * \param [in,out] cellIdsKeptArr array where all candidate cell ids are put at the end.
 */
void MEDCoupling1DGTUMesh::fillCellIdsToKeepFromNodeIds(const int *begin, const int *end, bool fullyIn, DataArrayInt *&cellIdsKeptArr) const
{
  int nbOfCells=getNumberOfCells();
  MCAuto<DataArrayInt> cellIdsKept=DataArrayInt::New(); cellIdsKept->alloc(0,1);
  int tmp=-1;
  int sz=_conn->getMaxValue(tmp); sz=std::max(sz,0)+1;
  std::vector<bool> fastFinder(sz,false);
  for(const int *work=begin;work!=end;work++)
    if(*work>=0 && *work<sz)
      fastFinder[*work]=true;
  const int *conn=_conn->begin(),*conni=_conn_indx->begin();
  for(int i=0;i<nbOfCells;i++,conni++)
    {
      int ref=0,nbOfHit=0;
      int nbNodesPerCell=conni[1]-conni[0];
      if(nbNodesPerCell>=0)
        {
          for(int j=0;j<nbNodesPerCell;j++)
            {
              int nodeId=conn[conni[0]+j];
              if(nodeId>=0)
                {
                  ref++;
                  if(fastFinder[nodeId])
                    nbOfHit++;
                }
            }
        }
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::fillCellIdsToKeepFromNodeIds : invalid index array for cell #" << i << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if((ref==nbOfHit && fullyIn) || (nbOfHit!=0 && !fullyIn))
        cellIdsKept->pushBackSilent(i);
    }
  cellIdsKeptArr=cellIdsKept.retn();
}

void MEDCoupling1DGTUMesh::allocateCells(int nbOfCells)
{
  if(nbOfCells<0)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::allocateCells : the input number of cells should be >= 0 !");
  _conn=DataArrayInt::New();
  _conn->reserve(nbOfCells*3);
  _conn_indx=DataArrayInt::New();
  _conn_indx->reserve(nbOfCells+1); _conn_indx->pushBackSilent(0);
  declareAsNew();
}

/*!
 * Appends at the end of \a this a cell having nodal connectivity array defined in [ \a nodalConnOfCellBg, \a nodalConnOfCellEnd ).
 *
 * \param [in] nodalConnOfCellBg - the begin (included) of nodal connectivity of the cell to add.
 * \param [in] nodalConnOfCellEnd - the end (excluded) of nodal connectivity of the cell to add.
 * \throw If the length of the input nodal connectivity array of the cell to add is not equal to number of nodes per cell relative to the unique geometric type
 *        attached to \a this.
 * \thow If the nodal connectivity array in \a this is null (call MEDCoupling1SGTUMesh::allocateCells before).
 */
void MEDCoupling1DGTUMesh::insertNextCell(const int *nodalConnOfCellBg, const int *nodalConnOfCellEnd)
{
  std::size_t sz(std::distance(nodalConnOfCellBg,nodalConnOfCellEnd));
  DataArrayInt *c(_conn),*c2(_conn_indx);
  if(c && c2)
    {
      int pos=c2->back();
      if(pos==(int)c->getNumberOfTuples())
        {
          c->pushBackValsSilent(nodalConnOfCellBg,nodalConnOfCellEnd);
          c2->pushBackSilent(pos+sz);
        }
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::insertNextCell : The nodal index array (end=" << pos << ") mismatches with nodal array (length=" << c->getNumberOfTuples() << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  else
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::insertNextCell : nodal connectivity array is null ! Call MEDCoupling1DGTUMesh::allocateCells before !");
}

void MEDCoupling1DGTUMesh::setNodalConnectivity(DataArrayInt *nodalConn, DataArrayInt *nodalConnIndex)
{
  if(nodalConn)
    nodalConn->incrRef();
  _conn=nodalConn;
  if(nodalConnIndex)
    nodalConnIndex->incrRef();
  _conn_indx=nodalConnIndex;
  declareAsNew();
}

/*!
 * \return DataArrayInt * - the internal reference to the nodal connectivity. The caller is not reponsible to deallocate it.
 */
DataArrayInt *MEDCoupling1DGTUMesh::getNodalConnectivity() const
{
  const DataArrayInt *ret(_conn);
  return const_cast<DataArrayInt *>(ret);
}

/*!
 * \return DataArrayInt * - the internal reference to the nodal connectivity index. The caller is not reponsible to deallocate it.
 */
DataArrayInt *MEDCoupling1DGTUMesh::getNodalConnectivityIndex() const
{
  const DataArrayInt *ret(_conn_indx);
  return const_cast<DataArrayInt *>(ret);
}

/*!
 * See the definition of the nodal connectivity pack \ref MEDCoupling1DGTUMesh::isPacked "here".
 * This method tries to build a new instance geometrically equivalent to \a this, by limiting at most the number of new object (nodal connectivity).
 * Geometrically the returned mesh is equal to \a this. So if \a this is already packed, the return value is a shallow copy of \a this.
 *
 * Whatever the status of pack of \a this, the coordinates array of the returned newly created instance is the same than those in \a this.
 * 
 * \param [out] isShallowCpyOfNodalConnn - tells if the returned instance share the same pair of nodal connectivity arrays (true) or if nodal
 *              connectivity arrays are different (false)
 * \return a new object to be managed by the caller.
 * 
 * \sa MEDCoupling1DGTUMesh::retrievePackedNodalConnectivity, MEDCoupling1DGTUMesh::isPacked
 */
MEDCoupling1DGTUMesh *MEDCoupling1DGTUMesh::copyWithNodalConnectivityPacked(bool& isShallowCpyOfNodalConnn) const
{
  MCAuto<MEDCoupling1DGTUMesh> ret(new MEDCoupling1DGTUMesh(getName(),*_cm));
  DataArrayInt *nc=0,*nci=0;
  isShallowCpyOfNodalConnn=retrievePackedNodalConnectivity(nc,nci);
  MCAuto<DataArrayInt> ncs(nc),ncis(nci);
  ret->_conn=ncs; ret->_conn_indx=ncis;
  ret->setCoords(getCoords());
  return ret.retn();
}

/*!
 * This method allows to compute, if needed, the packed nodal connectivity pair.
 * Indeed, it is possible to store in \a this a nodal connectivity array bigger than ranges convered by nodal connectivity index array.
 * It is typically the case when nodalConnIndx starts with an id greater than 0, and finishes with id less than number of tuples in \c this->_conn.
 * 
 * If \a this looks packed (the front of nodal connectivity index equal to 0 and back of connectivity index equal to number of tuple of nodal connectivity array)
 * true will be returned and respectively \a this->_conn and \a this->_conn_indx (with ref counter incremented). This is the classical case.
 *
 * If nodal connectivity index points to a subpart of nodal connectivity index the packed pair of arrays will be computed (new objects) and returned and false
 * will be returned.
 * 
 * This method return 3 elements.
 * \param [out] nodalConn - a pointer that can be equal to \a this->_conn if true is returned (general case). Whatever the value of return parameter
 *                          this pointer can be seen as a new object, that is to managed by the caller.
 * \param [out] nodalConnIndx - a pointer that can be equal to \a this->_conn_indx if true is returned (general case). Whatever the value of return parameter
 *                              this pointer can be seen as a new object, that is to managed by the caller.
 * \return bool - an indication of the content of the 2 output parameters. If true, \a this looks packed (general case), if true, \a this is not packed then
 * output parameters are newly created objects.
 *
 * \throw if \a this does not pass MEDCoupling1DGTUMesh::checkConsistencyLight test
 */
bool MEDCoupling1DGTUMesh::retrievePackedNodalConnectivity(DataArrayInt *&nodalConn, DataArrayInt *&nodalConnIndx) const
{
  if(isPacked())//performs the checkConsistencyLight
    {
      const DataArrayInt *c0(_conn),*c1(_conn_indx);
      nodalConn=const_cast<DataArrayInt *>(c0); nodalConnIndx=const_cast<DataArrayInt *>(c1);
      nodalConn->incrRef(); nodalConnIndx->incrRef();
      return true;
    }
  int bg=_conn_indx->front(),end=_conn_indx->back();
  MCAuto<DataArrayInt> nc(_conn->selectByTupleIdSafeSlice(bg,end,1));
  MCAuto<DataArrayInt> nci(_conn_indx->deepCopy());
  nci->applyLin(1,-bg);
  nodalConn=nc.retn(); nodalConnIndx=nci.retn();
  return false;
}

/*
 * If \a this looks packed (the front of nodal connectivity index equal to 0 and back of connectivity index equal to number of tuple of nodal connectivity array)
 * true will be returned and respectively \a this->_conn and \a this->_conn_indx (with ref counter incremented). This is the classical case.
 * If nodal connectivity index points to a subpart of nodal connectivity index false will be returned.
 * \return bool - true if \a this looks packed, false is not.
 *
 * \throw if \a this does not pass MEDCoupling1DGTUMesh::checkConsistencyLight test
 */
bool MEDCoupling1DGTUMesh::isPacked() const
{
  checkConsistencyLight();
  return _conn_indx->front()==0 && _conn_indx->back()==(int)_conn->getNumberOfTuples();
}

MEDCoupling1DGTUMesh *MEDCoupling1DGTUMesh::Merge1DGTUMeshes(const MEDCoupling1DGTUMesh *mesh1, const MEDCoupling1DGTUMesh *mesh2)
{
  std::vector<const MEDCoupling1DGTUMesh *> tmp(2);
  tmp[0]=const_cast<MEDCoupling1DGTUMesh *>(mesh1); tmp[1]=const_cast<MEDCoupling1DGTUMesh *>(mesh2);
  return Merge1DGTUMeshes(tmp);
}

MEDCoupling1DGTUMesh *MEDCoupling1DGTUMesh::Merge1DGTUMeshes(std::vector<const MEDCoupling1DGTUMesh *>& a)
{
  std::size_t sz=a.size();
  if(sz==0)
    return Merge1DGTUMeshesLL(a);
  for(std::size_t ii=0;ii<sz;ii++)
    if(!a[ii])
      {
        std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::Merge1DGTUMeshes : item #" << ii << " in input array of size "<< sz << " is empty !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  const INTERP_KERNEL::CellModel *cm=&(a[0]->getCellModel());
  for(std::size_t ii=0;ii<sz;ii++)
    if(&(a[ii]->getCellModel())!=cm)
      throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::Merge1DGTUMeshes : all items must have the same geo type !");
  std::vector< MCAuto<MEDCoupling1DGTUMesh> > bb(sz);
  std::vector< const MEDCoupling1DGTUMesh * > aa(sz);
  int spaceDim=-3;
  for(std::size_t i=0;i<sz && spaceDim==-3;i++)
    {
      const MEDCoupling1DGTUMesh *cur=a[i];
      const DataArrayDouble *coo=cur->getCoords();
      if(coo)
        spaceDim=coo->getNumberOfComponents();
    }
  if(spaceDim==-3)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::Merge1DGTUMeshes : no spaceDim specified ! unable to perform merge !");
  for(std::size_t i=0;i<sz;i++)
    {
      bb[i]=a[i]->buildSetInstanceFromThis(spaceDim);
      aa[i]=bb[i];
    }
  return Merge1DGTUMeshesLL(aa);
}

/*!
 * \throw If presence of a null instance in the input vector \a a.
 * \throw If a is empty
 */
MEDCoupling1DGTUMesh *MEDCoupling1DGTUMesh::Merge1DGTUMeshesOnSameCoords(std::vector<const MEDCoupling1DGTUMesh *>& a)
{
  if(a.empty())
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::Merge1DGTUMeshesOnSameCoords : input array must be NON EMPTY !");
  std::vector<const MEDCoupling1DGTUMesh *>::const_iterator it=a.begin();
  if(!(*it))
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::Merge1DGTUMeshesOnSameCoords : null instance in the first element of input vector !");
  std::vector< MCAuto<MEDCoupling1DGTUMesh> > objs(a.size());
  std::vector<const DataArrayInt *> ncs(a.size()),ncis(a.size());
  (*it)->getNumberOfCells();//to check that all is OK
  const DataArrayDouble *coords=(*it)->getCoords();
  const INTERP_KERNEL::CellModel *cm=&((*it)->getCellModel());
  bool tmp;
  objs[0]=(*it)->copyWithNodalConnectivityPacked(tmp);
  ncs[0]=objs[0]->getNodalConnectivity(); ncis[0]=objs[0]->getNodalConnectivityIndex();
  it++;
  for(int i=1;it!=a.end();i++,it++)
    {
      if(!(*it))
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::Merge1DGTUMeshesOnSameCoords : presence of null instance !");
      if(cm!=&((*it)->getCellModel()))
        throw INTERP_KERNEL::Exception("Geometric types mismatches, Merge1DGTUMeshes impossible !");
      (*it)->getNumberOfCells();//to check that all is OK
      objs[i]=(*it)->copyWithNodalConnectivityPacked(tmp);
      ncs[i]=objs[i]->getNodalConnectivity(); ncis[i]=objs[i]->getNodalConnectivityIndex();
      if(coords!=(*it)->getCoords())
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::Merge1DGTUMeshesOnSameCoords : not lying on same coords !");
    }
  MCAuto<MEDCoupling1DGTUMesh> ret(new MEDCoupling1DGTUMesh("merge",*cm));
  ret->setCoords(coords);
  ret->_conn=DataArrayInt::Aggregate(ncs);
  ret->_conn_indx=DataArrayInt::AggregateIndexes(ncis);
  return ret.retn();
}

/*!
 * Assume that all instances in \a a are non null. If null it leads to a crash. That's why this method is assigned to be low level (LL)
 */
MEDCoupling1DGTUMesh *MEDCoupling1DGTUMesh::Merge1DGTUMeshesLL(std::vector<const MEDCoupling1DGTUMesh *>& a)
{
  if(a.empty())
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::Merge1DGTUMeshes : input array must be NON EMPTY !");
  std::vector< MCAuto<MEDCoupling1DGTUMesh> > objs(a.size());
  std::vector<const DataArrayInt *> ncs(a.size()),ncis(a.size());
  std::vector<const MEDCoupling1DGTUMesh *>::const_iterator it=a.begin();
  std::vector<int> nbNodesPerElt(a.size());
  int nbOfCells=(*it)->getNumberOfCells();
  bool tmp;
  objs[0]=(*it)->copyWithNodalConnectivityPacked(tmp);
  ncs[0]=objs[0]->getNodalConnectivity(); ncis[0]=objs[0]->getNodalConnectivityIndex();
  nbNodesPerElt[0]=0;
  int prevNbOfNodes=(*it)->getNumberOfNodes();
  const INTERP_KERNEL::CellModel *cm=&((*it)->getCellModel());
  it++;
  for(int i=1;it!=a.end();i++,it++)
    {
      if(cm!=&((*it)->getCellModel()))
        throw INTERP_KERNEL::Exception("Geometric types mismatches, Merge1DGTUMeshes impossible !");
      objs[i]=(*it)->copyWithNodalConnectivityPacked(tmp);
      ncs[i]=objs[i]->getNodalConnectivity(); ncis[i]=objs[i]->getNodalConnectivityIndex();
      nbOfCells+=(*it)->getNumberOfCells();
      nbNodesPerElt[i]=nbNodesPerElt[i-1]+prevNbOfNodes;
      prevNbOfNodes=(*it)->getNumberOfNodes();
    }
  std::vector<const MEDCouplingPointSet *> aps(a.size());
  std::copy(a.begin(),a.end(),aps.begin());
  MCAuto<DataArrayDouble> pts=MergeNodesArray(aps);
  MCAuto<MEDCoupling1DGTUMesh> ret(new MEDCoupling1DGTUMesh("merge",*cm));
  ret->setCoords(pts);
  ret->_conn=AggregateNodalConnAndShiftNodeIds(ncs,nbNodesPerElt);
  ret->_conn_indx=DataArrayInt::AggregateIndexes(ncis);
  return ret.retn();
}

MEDCoupling1DGTUMesh *MEDCoupling1DGTUMesh::buildSetInstanceFromThis(int spaceDim) const
{
  MCAuto<MEDCoupling1DGTUMesh> ret(new MEDCoupling1DGTUMesh(getName(),*_cm));
  MCAuto<DataArrayInt> tmp1,tmp2;
  const DataArrayInt *nodalConn(_conn),*nodalConnI(_conn_indx);
  if(!nodalConn)
    {
      tmp1=DataArrayInt::New(); tmp1->alloc(0,1);
    }
  else
    tmp1=_conn;
  ret->_conn=tmp1;
  //
  if(!nodalConnI)
    {
      tmp2=DataArrayInt::New(); tmp2->alloc(1,1); tmp2->setIJ(0,0,0);
    }
  else
    tmp2=_conn_indx;
  ret->_conn_indx=tmp2;
  //
  if(!_coords)
    {
      MCAuto<DataArrayDouble> coords=DataArrayDouble::New(); coords->alloc(0,spaceDim);
      ret->setCoords(coords);
    }
  else
    ret->setCoords(_coords);
  return ret.retn();
}

/*!
 * This method aggregate the bbox of each cell and put it into bbox parameter.
 * 
 * \param [in] arcDetEps - a parameter specifying in case of 2D quadratic polygon cell the detection limit between linear and arc circle. (By default 1e-12)
 *                         For all other cases this input parameter is ignored.
 * \return DataArrayDouble * - newly created object (to be managed by the caller) \a this number of cells tuples and 2*spacedim components.
 * 
 * \throw If \a this is not fully set (coordinates and connectivity).
 * \throw If a cell in \a this has no valid nodeId.
 */
DataArrayDouble *MEDCoupling1DGTUMesh::getBoundingBoxForBBTree(double arcDetEps) const
{
  checkFullyDefined();
  int spaceDim(getSpaceDimension()),nbOfCells(getNumberOfCells()),nbOfNodes(getNumberOfNodes());
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New()); ret->alloc(nbOfCells,2*spaceDim);
  double *bbox(ret->getPointer());
  for(int i=0;i<nbOfCells*spaceDim;i++)
    {
      bbox[2*i]=std::numeric_limits<double>::max();
      bbox[2*i+1]=-std::numeric_limits<double>::max();
    }
  const double *coordsPtr(_coords->getConstPointer());
  const int *conn(_conn->getConstPointer()),*connI(_conn_indx->getConstPointer());
  for(int i=0;i<nbOfCells;i++)
    {
      int offset=connI[i];
      int nbOfNodesForCell(connI[i+1]-offset),kk(0);
      for(int j=0;j<nbOfNodesForCell;j++)
        {
          int nodeId=conn[offset+j];
          if(nodeId>=0 && nodeId<nbOfNodes)
            {
              for(int k=0;k<spaceDim;k++)
                {
                  bbox[2*spaceDim*i+2*k]=std::min(bbox[2*spaceDim*i+2*k],coordsPtr[spaceDim*nodeId+k]);
                  bbox[2*spaceDim*i+2*k+1]=std::max(bbox[2*spaceDim*i+2*k+1],coordsPtr[spaceDim*nodeId+k]);
                }
              kk++;
            }
        }
      if(kk==0)
        {
          std::ostringstream oss; oss << "MEDCoupling1SGTUMesh::getBoundingBoxForBBTree : cell #" << i << " contains no valid nodeId !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret.retn();
}

/*!
 * Returns the cell field giving for each cell in \a this its diameter. Diameter means the max length of all possible SEG2 in the cell.
 *
 * \return a new instance of field containing the result. The returned instance has to be deallocated by the caller.
 */
MEDCouplingFieldDouble *MEDCoupling1DGTUMesh::computeDiameterField() const
{
  throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::computeDiameterField : not implemented yet for dynamic types !");
}

std::vector<int> MEDCoupling1DGTUMesh::BuildAPolygonFromParts(const std::vector< std::vector<int> >& parts)
{
  std::vector<int> ret;
  if(parts.empty())
    return ret;
  ret.insert(ret.end(),parts[0].begin(),parts[0].end());
  int ref(ret.back());
  std::size_t sz(parts.size()),nbh(1);
  std::vector<bool> b(sz,true); b[0]=false;
  while(nbh<sz)
    {
      std::size_t i(0);
      for(;i<sz;i++) if(b[i] && parts[i].front()==ref) { ret.insert(ret.end(),parts[i].begin()+1,parts[i].end()); nbh++; break; }
      if(i<sz)
        ref=ret.back();
      else
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::BuildAPolygonFromParts : the input vector is not a part of a single polygon !");
    }
  if(ret.back()==ret.front())
    ret.pop_back();
  return ret;
}

/*!
 * This method invert orientation of all cells in \a this. 
 * After calling this method the absolute value of measure of cells in \a this are the same than before calling.
 * This method only operates on the connectivity so coordinates are not touched at all.
 */
void MEDCoupling1DGTUMesh::invertOrientationOfAllCells()
{
  checkConsistencyOfConnectivity();
  INTERP_KERNEL::AutoCppPtr<INTERP_KERNEL::OrientationInverter> oi(INTERP_KERNEL::OrientationInverter::BuildInstanceFrom(getCellModelEnum()));
  int nbCells(getNumberOfCells());
  const int *connI(_conn_indx->begin());
  int *conn(_conn->getPointer());
  for(int i=0;i<nbCells;i++)
    oi->operate(conn+connI[i],conn+connI[i+1]);
  updateTime();
}

/*!
 * This method performs an aggregation of \a nodalConns (as DataArrayInt::Aggregate does) but in addition of that a shift is applied on the 
 * values contained in \a nodalConns using corresponding offset specified in input \a offsetInNodeIdsPerElt.
 * But it also manage the values -1, that have a semantic in MEDCoupling1DGTUMesh class (separator for polyhedron).
 *
 * \param [in] nodalConns - a list of nodal connectivity arrays same size than \a offsetInNodeIdsPerElt.
 * \param [in] offsetInNodeIdsPerElt - a list of offsets to apply.
 * \return DataArrayInt * - A new object (to be managed by the caller) that is the result of the aggregation.
 * \throw If \a nodalConns or \a offsetInNodeIdsPerElt are empty.
 * \throw If \a nodalConns and \a offsetInNodeIdsPerElt have not the same size.
 * \throw If presence of null pointer in \a nodalConns.
 * \throw If presence of not allocated or array with not exactly one component in \a nodalConns.
 */
DataArrayInt *MEDCoupling1DGTUMesh::AggregateNodalConnAndShiftNodeIds(const std::vector<const DataArrayInt *>& nodalConns, const std::vector<int>& offsetInNodeIdsPerElt)
{
  std::size_t sz1(nodalConns.size()),sz2(offsetInNodeIdsPerElt.size());
  if(sz1!=sz2)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::AggregateNodalConnAndShiftNodeIds : input vectors do not have the same size !");
  if(sz1==0)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::AggregateNodalConnAndShiftNodeIds : empty vectors in input !");
  int nbOfTuples=0;
  for(std::vector<const DataArrayInt *>::const_iterator it=nodalConns.begin();it!=nodalConns.end();it++)
    {
      if(!(*it))
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::AggregateNodalConnAndShiftNodeIds : presence of null pointer in input vector !");
      if(!(*it)->isAllocated())
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::AggregateNodalConnAndShiftNodeIds : presence of non allocated array in input vector !");
      if((*it)->getNumberOfComponents()!=1)
        throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::AggregateNodalConnAndShiftNodeIds : presence of array with not exactly one component !");
      nbOfTuples+=(*it)->getNumberOfTuples();
    }
  MCAuto<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(nbOfTuples,1);
  int *pt=ret->getPointer();
  int i=0;
  for(std::vector<const DataArrayInt *>::const_iterator it=nodalConns.begin();it!=nodalConns.end();it++,i++)
    {
      int curNbt=(*it)->getNumberOfTuples();
      const int *inPt=(*it)->begin();
      int offset=offsetInNodeIdsPerElt[i];
      for(int j=0;j<curNbt;j++,pt++)
        {
          if(inPt[j]!=-1)
            *pt=inPt[j]+offset;
          else
            *pt=-1;
        }
    }
  return ret.retn();
}

MEDCoupling1DGTUMesh *MEDCoupling1DGTUMesh::New(const MEDCouplingUMesh *m)
{
  if(!m)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::New : input mesh is null !");
  std::set<INTERP_KERNEL::NormalizedCellType> gts(m->getAllGeoTypes());
  if(gts.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCoupling1DGTUMesh::New : input mesh must have exactly one geometric type !");
  int geoType((int)*gts.begin());
  MCAuto<MEDCoupling1DGTUMesh> ret(MEDCoupling1DGTUMesh::New(m->getName(),*gts.begin()));
  ret->setCoords(m->getCoords()); ret->setDescription(m->getDescription());
  int nbCells(m->getNumberOfCells());
  MCAuto<DataArrayInt> conn(DataArrayInt::New()),connI(DataArrayInt::New());
  conn->alloc(m->getNodalConnectivityArrayLen()-nbCells,1); connI->alloc(nbCells+1,1);
  int *c(conn->getPointer()),*ci(connI->getPointer()); *ci=0;
  const int *cin(m->getNodalConnectivity()->begin()),*ciin(m->getNodalConnectivityIndex()->begin());
  for(int i=0;i<nbCells;i++,ciin++,ci++)
    {
      if(cin[ciin[0]]==geoType)
        {
          if(ciin[1]-ciin[0]>=1)
            {
              c=std::copy(cin+ciin[0]+1,cin+ciin[1],c);
              ci[1]=ci[0]+ciin[1]-ciin[0]-1;
            }
          else
            {
              std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::New(const MEDCouplingUMesh *m) : something is wrong in the input mesh at cell #" << i << " ! The size of cell is not >=0 !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "MEDCoupling1DGTUMesh::New(const MEDCouplingUMesh *m) : something is wrong in the input mesh at cell #" << i << " ! The geometric type is not those expected !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  ret->setNodalConnectivity(conn,connI);
  return ret.retn();
}
