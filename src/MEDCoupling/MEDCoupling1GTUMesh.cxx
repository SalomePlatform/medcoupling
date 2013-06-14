// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#include "MEDCoupling1GTUMesh.hxx"

using namespace ParaMEDMEM;

MEDCoupling1GTUMesh *MEDCoupling1GTUMesh::New(const char *meshName, INTERP_KERNEL::NormalizedCellType type) throw(INTERP_KERNEL::Exception)
{
  if(type==INTERP_KERNEL::NORM_ERROR)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::New : NORM_ERROR is not a valid type to be used as base geometric type for a mesh !");
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  if(!cm.isDynamic())
    return MEDCoupling1SGTUMesh::New(meshName,type);
  throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::New : not implemented yet !");
}

const INTERP_KERNEL::CellModel& MEDCoupling1GTUMesh::getCellModel() const throw(INTERP_KERNEL::Exception)
{
  return _cm;
}

INTERP_KERNEL::NormalizedCellType MEDCoupling1GTUMesh::getCellModelEnum() const throw(INTERP_KERNEL::Exception)
{
  return _cm.getEnum();
}

int MEDCoupling1GTUMesh::getMeshDimension() const
{
  return (int)_cm.getDimension();
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
DataArrayInt *MEDCoupling1GTUMesh::giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
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
int MEDCoupling1GTUMesh::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  return type==getCellModelEnum()?getNumberOfCells():0;
}

/*!
 * Returns a type of a cell by its id.
 *  \param [in] cellId - the id of the cell of interest.
 *  \return INTERP_KERNEL::NormalizedCellType - enumeration item describing the cell type.
 *  \throw If \a cellId is invalid. Valid range is [0, \a this->getNumberOfCells() ).
 */
INTERP_KERNEL::NormalizedCellType MEDCoupling1GTUMesh::getTypeOfCell(int cellId) const
{
  if(cellId>=0 && cellId<getNumberOfCells())
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
 * For every k in [0,n] ret[3*k+2]==0 because it has no sense here. 
 * This parameter is kept only for compatibility with other methode listed above.
 */
std::vector<int> MEDCoupling1GTUMesh::getDistributionOfTypes() const throw(INTERP_KERNEL::Exception)
{
  std::vector<int> ret(3);
  ret[0]=(int)getCellModelEnum(); ret[1]=getNumberOfCells(); ret[2]=0;
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
 */
void MEDCoupling1GTUMesh::splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception)
{
  if(!profile)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::splitProfilePerType : input profile is NULL !");
  if(profile->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCoupling1GTUMesh::splitProfilePerType : input profile should have exactly one component !");
  int nbOfCells=getNumberOfCells();
  
}
