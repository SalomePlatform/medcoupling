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

#include "MEDCouplingField.hxx"
#include "MEDCouplingMesh.hxx"
#include "MEDCouplingFieldDiscretization.hxx"

#include <sstream>

using namespace MEDCoupling;

void MEDCouplingField::checkConsistencyLight() const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Field invalid because no mesh specified !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("MEDCouplingField::checkConsistencyLight : no spatial discretization !");
}

bool MEDCouplingField::isEqualIfNotWhyProtected(const MEDCouplingField *other, double meshPrec, std::string& reason) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingField::isEqualIfNotWhy : other instance is NULL !");
  std::ostringstream oss; oss.precision(15);
  if(_name!=other->_name)
    {
      oss << "Field names differ : this name = \"" << _name << "\" and other name = \"" << other->_name << "\" !";
      reason=oss.str();
      return false;
    }
  if(_desc!=other->_desc)
    {
      oss << "Field descriptions differ : this description = \"" << _desc << "\" and other description = \"" << other->_desc << "\" !";
      reason=oss.str();
      return false;
    }
  if(_nature!=other->_nature)
    {
      oss << "Field nature differ : this nature = \"" << MEDCouplingNatureOfField::GetRepr(_nature) << "\" and other nature = \"" << MEDCouplingNatureOfField::GetRepr(other->_nature) << "\" !";
      reason=oss.str();
      return false;
    }
  if(!_type->isEqualIfNotWhy(other->_type,meshPrec,reason))
    {
      reason.insert(0,"Spatial discretizations differ :");
      return false;
    }
  if(_mesh==0 && other->_mesh==0)
    return true;
  if(_mesh==0 || other->_mesh==0)
    {
      reason="Only one field between the two this and other has its underlying mesh defined !";
      return false;
    }
  if(_mesh==other->_mesh)
    return true;
  bool ret=_mesh->isEqualIfNotWhy(other->_mesh,meshPrec,reason);
  if(!ret)
    reason.insert(0,"Underlying meshes of fields differ for the following reason : ");
  return ret;
}

/*!
 * Checks if \a this and another MEDCouplingField are equal. The textual
 * information like names etc. is not considered.
 *  \param [in] other - the field to compare with \a this one.
 *  \param [in] meshPrec - precision used to compare node coordinates of the underlying mesh.
 *  \return bool - \c true if the two fields are equal, \c false else.
 *  \throw If \a other is NULL.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
bool MEDCouplingField::isEqualWithoutConsideringStrProtected(const MEDCouplingField *other, double meshPrec) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingField::isEqualWithoutConsideringStr : input field is NULL !");
  if(!_type)
    throw INTERP_KERNEL::Exception("MEDCouplingField::isEqualWithoutConsideringStr : spatial discretization of this is NULL !");
  if(!_type->isEqualWithoutConsideringStr(other->_type,meshPrec))
    return false;
  if(_nature!=other->_nature)
    return false;
  if(_mesh==0 && other->_mesh==0)
    return true;
  if(_mesh==0 || other->_mesh==0)
    return false;
  if(_mesh==other->_mesh)
    return true;
  return _mesh->isEqualWithoutConsideringStr(other->_mesh,meshPrec);
}

/*!
 * This method states if 'this' and 'other' are compatibles each other before performing any treatment.
 * This method is good for methods like : mergeFields.
 * This method is not very demanding compared to areStrictlyCompatible that is better for operation on fields.
 */
bool MEDCouplingField::areCompatibleForMerge(const MEDCouplingField *other) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingField::areCompatibleForMerge : input field is NULL !");
  if(!_type->isEqual(other->_type,1.))
    return false;
  if(_nature!=other->_nature)
    return false;
  if(_mesh==other->_mesh)
    return true;
  return _mesh->areCompatibleForMerge(other->_mesh);
}

/*!
 * This method is more strict than MEDCouplingField::areCompatibleForMerge method.
 * This method is used for operation on fields to operate a first check before attempting operation.
 */
bool MEDCouplingField::areStrictlyCompatible(const MEDCouplingField *other) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingField::areStrictlyCompatible : input field is NULL !");
  if(!_type->isEqual(other->_type,1.e-12))
    return false;
  if(_nature!=other->_nature)
    return false;
  return _mesh==other->_mesh;
}

/*!
 * This method is less strict than MEDCouplingField::areStrictlyCompatible method.
 * The difference is that the nature is not checked.
 * This method is used for multiplication and division on fields to operate a first check before attempting operation.
 */
bool MEDCouplingField::areStrictlyCompatibleForMulDiv(const MEDCouplingField *other) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingField::areStrictlyCompatible : input field is NULL !");
  if(!_type->isEqual(other->_type,1.e-12))
    return false;
  return _mesh==other->_mesh;
}


void MEDCouplingField::updateTime() const
{
  if(_mesh)
    updateTimeWith(*_mesh);
  if(_type)
    updateTimeWith(*_type);
}

std::size_t MEDCouplingField::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret=0;
  ret+=_name.capacity();
  ret+=_desc.capacity();
  return ret;
}

std::vector<const BigMemoryObject *> MEDCouplingField::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back(_mesh);
  ret.push_back((const MEDCouplingFieldDiscretization *)_type);
  return ret;
}

/*!
 * Returns a type of \ref MEDCouplingSpatialDisc "spatial discretization" of \a this
 * field in terms of enum MEDCoupling::TypeOfField. 
 * \return MEDCoupling::TypeOfField - the type of \a this field.
 * \throw If the geometric type is empty.
 */
TypeOfField MEDCouplingField::getTypeOfField() const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("MEDCouplingField::getTypeOfField : spatial discretization is null !");
  return _type->getEnum();
}

/*!
 * Returns the nature of \a this field. This information is very important during
 * interpolation process using MEDCoupling::MEDCouplingRemapper or MEDCoupling::InterpKernelDEC.
 * In other context than the two mentioned above, this attribute is unimportant. This
 * attribute is not stored in the MED file.
 * For more information of the semantics and the influence of this attribute to the
 * result of interpolation, see
 * - \ref NatureOfField
 * - \ref TableNatureOfField "How interpolation coefficients depend on Field Nature"
 */
NatureOfField MEDCouplingField::getNature() const
{
  return _nature;
}

/*!
 * Sets the nature of \a this field. This information is very important during
 * interpolation process using MEDCoupling::MEDCouplingRemapper or MEDCoupling::InterpKernelDEC.
 * In other context than the two mentioned above, this attribute is unimportant. This
 * attribute is not stored in the MED file.
 * For more information of the semantics and the influence of this attribute to the
 * result of interpolation, see
 * - \ref NatureOfField
 * - \ref TableNatureOfField "How interpolation coefficients depend on Field Nature"
 *
 *  \param [in] nat - the nature of \a this field.
 *  \throw If \a nat has an invalid value.
 */
void MEDCouplingField::setNature(NatureOfField nat)
{
  MEDCouplingNatureOfField::GetRepr(nat);//generate a throw if nat not recognized
  if(_type)
    _type->checkCompatibilityWithNature(nat);
  _nature=nat;
}

/*!
 * Returns coordinates of field location points that depend on 
 * \ref MEDCouplingSpatialDisc "spatial discretization" of \a this field.
 * - For a field on nodes, returns coordinates of nodes.
 * - For a field on cells, returns barycenters of cells.
 * - For a field on gauss points, returns coordinates of gauss points.
 * 
 *  \return DataArrayDouble * - a new instance of DataArrayDouble. The caller is to
 *          delete this array using decrRef() as it is no more needed. 
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the mesh is not set.
 */
DataArrayDouble *MEDCouplingField::getLocalizationOfDiscr() const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingField::getLocalizationOfDiscr : No mesh set !");
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("MEDCouplingField::getLocalizationOfDiscr : No spatial discretization set !");
  return _type->getLocalizationOfDiscValues(_mesh);
}

/*!
 * Returns a new MEDCouplingFieldDouble containing volumes of cells of a dual mesh whose
 * cells are constructed around field location points (getLocalizationOfDiscr()) of \a this
 * field. (In case of a field on cells, the dual mesh coincides with the underlying mesh).<br>
 * For 1D cells, the returned field contains lengths.<br>
 * For 2D cells, the returned field contains areas.<br>
 * For 3D cells, the returned field contains volumes.
 *  \param [in] isAbs - if \c true, the computed cell volume does not reflect cell
 *         orientation, i.e. the volume is always positive.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble.
 *          The caller is to delete this array using decrRef() as
 *          it is no more needed.
 *  \throw If the mesh is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the spatial discretization of \a this field is not well defined.
 */

MEDCouplingFieldDouble *MEDCouplingField::buildMeasureField(bool isAbs) const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingField::buildMeasureField : no mesh defined !");
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("MEDCouplingField::buildMeasureField : No spatial discretization set !");
  return _type->getMeasureField(_mesh,isAbs);
}

/*!
 * Sets the underlying mesh of \a this field.
 * For examples of field construction, see \ref MEDCouplingFirstSteps3.
 *  \param [in] mesh - the new underlying mesh.
 */
void MEDCouplingField::setMesh(const MEDCouplingMesh *mesh)
{
  if(mesh!=_mesh)
    {
      if(_mesh)
        _mesh->decrRef();
      _mesh=mesh;
      declareAsNew();
      if(_mesh)
        {
          _mesh->incrRef();
          updateTimeWith(*_mesh);
        }
    }
}

/*!
 * Sets localization of Gauss points for a given geometric type of cell.
 *  \param [in] type - the geometric type of cell for which the Gauss localization is set.
 *  \param [in] refCoo - coordinates of points of the reference cell. Size of this vector
 *         must be \c nbOfNodesPerCell * \c dimOfType. 
 *  \param [in] gsCoo - coordinates of Gauss points on the reference cell. Size of this vector
 *         must be  _wg_.size() * \c dimOfType.
 *  \param [in] wg - the weights of Gauss points.
 *  \throw If \a this field is not on Gauss points.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the mesh is not set.
 *  \throw If size of any vector do not match the \a type.
 */
void MEDCouplingField::setGaussLocalizationOnType(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                  const std::vector<double>& gsCoo, const std::vector<double>& wg)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Mesh has to be set before calling setGaussLocalizationOnType method !");
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call setGaussLocalizationOnType method !");
  _type->setGaussLocalizationOnType(_mesh,type,refCoo,gsCoo,wg);
}

/*!
 * Sets localization of Gauss points for given cells specified by their ids.
 *  \param [in] begin - an array of cell ids of interest.
 *  \param [in] end - the end of \a begin, i.e. a pointer to its (last+1)-th element.
 *  \param [in] refCoo - coordinates of points of the reference cell. Size of this vector
 *         must be \c nbOfNodesPerCell * \c dimOfType. 
 *  \param [in] gsCoo - coordinates of Gauss points on the reference cell. Size of this vector
 *         must be  _wg_.size() * \c dimOfType.
 *  \param [in] wg - the weights of Gauss points.
 *  \throw If \a this field is not on Gauss points.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the mesh is not set.
 *  \throw If size of any vector do not match the type of cell # \a begin[0].
 *  \throw If type of any cell in \a begin differs from that of cell # \a begin[0].
 *  \throw If the range [_begin_,_end_) is empty.
 */
void MEDCouplingField::setGaussLocalizationOnCells(const int *begin, const int *end, const std::vector<double>& refCoo,
                                                   const std::vector<double>& gsCoo, const std::vector<double>& wg)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Mesh has to be set before calling setGaussLocalizationOnCells method !");
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call setGaussLocalizationOnCells method !");
  _type->setGaussLocalizationOnCells(_mesh,begin,end,refCoo,gsCoo,wg);
}

/*!
 * Clears data on Gauss points localization.
 *  \throw If \a this field is not on Gauss points.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
void MEDCouplingField::clearGaussLocalizations()
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call clearGaussLocalizations method !");
  _type->clearGaussLocalizations();
}

/*!
 * Returns a reference to the Gauss localization object by its id.
 * \warning This method is not const, so the returned object can be modified without any
 *          problem.
 *  \param [in] locId - the id of the Gauss localization object of interest.
 *         It must be in range <em> 0 <= locId < getNbOfGaussLocalization() </em>.
 *  \return \ref MEDCoupling::MEDCouplingGaussLocalization "MEDCouplingGaussLocalization" & - the
 *  Gauss localization object.
 *  \throw If \a this field is not on Gauss points.
 *  \throw If \a locId is not within the valid range.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
MEDCouplingGaussLocalization& MEDCouplingField::getGaussLocalization(int locId)
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call getGaussLocalization method !");
  return _type->getGaussLocalization(locId);
}

/*!
 * Returns an id of the Gauss localization object corresponding to a given cell type.
 *  \param [in] type - the cell type of interest.
 *  \return int - the id of the Gauss localization object.
 *  \throw If \a this field is not on Gauss points.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If no Gauss localization object found for the given cell \a type.
 *  \throw If more than one Gauss localization object found for the given cell \a type.
 */
int MEDCouplingField::getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call getGaussLocalizationIdOfOneType method !");
  return _type->getGaussLocalizationIdOfOneType(type);
}

/*!
 * Returns ids of Gauss localization objects corresponding to a given cell type.
 *  \param [in] type - the cell type of interest.
 *  \return std::set<int> - ids of the Gauss localization object.
 *  \throw If \a this field is not on Gauss points.
 *  \throw If the spatial discretization of \a this field is NULL
 */
std::set<int> MEDCouplingField::getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call getGaussLocalizationIdsOfOneType method !");
  return _type->getGaussLocalizationIdsOfOneType(type);
}

/*!
 * Returns number of Gauss localization objects available. Implicitly all ids in
 * [0,getNbOfGaussLocalization()) are valid Gauss localization ids. 
 *  \return int - the number of available Gauss localization objects.
 *  \throw If \a this field is not on Gauss points.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
int MEDCouplingField::getNbOfGaussLocalization() const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call getNbOfGaussLocalization method !");
  return _type->getNbOfGaussLocalization();
}

/*!
 * Returns an id of the Gauss localization object corresponding to a type of a given cell.
 *  \param [in] cellId - an id of the cell of interest.
 *  \return int - the id of the Gauss localization object.
 *  \throw If \a this field is not on Gauss points.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If no Gauss localization object found for the given cell.
 */
int MEDCouplingField::getGaussLocalizationIdOfOneCell(int cellId) const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call getGaussLocalizationIdOfOneCell method !");
  return _type->getGaussLocalizationIdOfOneCell(cellId);
}

/*!
 * Returns ids of cells that share the same Gauss localization given by its id.
 *  \param [in] locId - the id of the Gauss localization object of interest. 
 *         It must be in range <em> 0 <= locId < getNbOfGaussLocalization() </em>.
 *  \param [in,out] cellIds - a vector returning ids of found cells. It is cleared before
 *         filling in. It remains empty if no cells found.
 *  \throw If \a this field is not on Gauss points.
 *  \throw If \a locId is not within the valid range.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
void MEDCouplingField::getCellIdsHavingGaussLocalization(int locId, std::vector<int>& cellIds) const
{
  cellIds.clear();
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call getCellIdsHavingGaussLocalization method !");
  _type->getCellIdsHavingGaussLocalization(locId,cellIds);
}

/*!
 * Returns a reference to the Gauss localization object by its id.
 * \warning This method is const, so the returned object is not apt for modification.
 *  \param [in] locId - the id of the Gauss localization object of interest.
 *         It must be in range <em> 0 <= locId < getNbOfGaussLocalization() </em>.
 *  \return const \ref MEDCouplingGaussLocalization & - the Gauss localization object.
 *  \throw If \a this field is not on Gauss points.
 *  \throw If \a locId is not within the valid range.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
const MEDCouplingGaussLocalization& MEDCouplingField::getGaussLocalization(int locId) const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call getGaussLocalization method !");
  return _type->getGaussLocalization(locId);
}

MEDCouplingField::~MEDCouplingField()
{
  if(_mesh)
    _mesh->decrRef();
}

MEDCouplingField::MEDCouplingField(MEDCouplingFieldDiscretization *type, NatureOfField nature):_nature(nature),_mesh(0),_type(type)
{
}

MEDCouplingField::MEDCouplingField(TypeOfField type):_nature(NoNature),_mesh(0),_type(MEDCouplingFieldDiscretization::New(type))
{
}

MEDCouplingField::MEDCouplingField(const MEDCouplingField& other, bool deepCopy):RefCountObject(other),_name(other._name),_desc(other._desc),_nature(other._nature),
    _mesh(0),_type(0)
{
  if(other._mesh)
    {
      _mesh=other._mesh;
      _mesh->incrRef();
    }
  if(deepCopy)
    _type=other._type->clone();
  else
    _type=other._type;
}

/*!
 * Returns a new MEDCouplingMesh constituted by some cells of the underlying mesh of \a
 * this field, and returns ids of entities (nodes, cells, Gauss points) lying on the
 * specified cells. The cells to include to the result mesh are specified by an array of
 * cell ids. The new mesh shares the coordinates array with the underlying mesh. 
 *  \param [in] start - an array of cell ids to include to the result mesh.
 *  \param [in] end - specifies the end of the array \a start, so that
 *              the last value of \a start is \a end[ -1 ].
 *  \param [out] di - a new instance of DataArrayInt holding the ids of entities (nodes,
 *         cells, Gauss points). The caller is to delete this array using decrRef() as it
 *         is no more needed.  
 *  \return MEDCouplingMesh * - a new instance of MEDCouplingMesh. The caller is to
 *         delete this mesh using decrRef() as it is no more needed. 
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the mesh is not set.
 * \sa buildSubMeshDataRange()
 */
MEDCouplingMesh *MEDCouplingField::buildSubMeshData(const int *start, const int *end, DataArrayInt *&di) const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call buildSubMeshData method !");
  return _type->buildSubMeshData(_mesh,start,end,di);
}

/*!
 * This method returns a submesh of 'mesh' instance constituting cell ids defined by a range given by the 3 following inputs \a begin, \a end and \a step.
 * 
 * \param [out] beginOut Valid only if \a di is NULL
 * \param [out] endOut Valid only if \a di is NULL
 * \param [out] stepOut Valid only if \a di is NULL
 * \param [out] di is an array returned that specifies entity ids (nodes, cells, Gauss points... ) in array if no output range is foundable.
 * 
 * \sa MEDCouplingField::buildSubMeshData
 */
MEDCouplingMesh *MEDCouplingField::buildSubMeshDataRange(int begin, int end, int step, int& beginOut, int& endOut, int& stepOut, DataArrayInt *&di) const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call buildSubMeshDataRange method !");
  return _type->buildSubMeshDataRange(_mesh,begin,end,step,beginOut,endOut,stepOut,di);
}

/*!
 * This method returns tuples ids implied by the mesh selection of the  cell ids contained in array defined as an interval [start;end).
 * \return a newly allocated DataArrayInt instance containing tuples ids.
 */
DataArrayInt *MEDCouplingField::computeTupleIdsToSelectFromCellIds(const int *startCellIds, const int *endCellIds) const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call computeTupleIdsToSelectFromCellIds method !");
  return _type->computeTupleIdsToSelectFromCellIds(_mesh,startCellIds,endCellIds);
}

/*!
 * Returns number of tuples expected regarding the spatial discretization of \a this
 * field and number of entities in the underlying mesh. This method behaves exactly as MEDCouplingFieldDouble::getNumberOfTuples.
 *  \return int - the number of expected tuples.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the mesh is not set.
 */
int MEDCouplingField::getNumberOfTuplesExpected() const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call getNumberOfTuplesExpected method !");
  if(_mesh)
    return _type->getNumberOfTuples(_mesh);
  else
    throw INTERP_KERNEL::Exception("MEDCouplingField::getNumberOfTuplesExpected : Empty mesh !");
}

void MEDCouplingField::setDiscretization(MEDCouplingFieldDiscretization *newDisc)
{
  bool needUpdate=(const MEDCouplingFieldDiscretization *)_type!=newDisc;
  _type=newDisc;
  if(newDisc)
    newDisc->incrRef();
  if(needUpdate)
    declareAsNew();
}

/*!
 * Returns number of mesh entities in the underlying mesh of \a this field regarding the
 * spatial discretization.
 *  \return int - the number of mesh entities porting the field values.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the mesh is not set.
 */
int MEDCouplingField::getNumberOfMeshPlacesExpected() const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("Spatial discretization not set ! Impossible to call getNumberOfMeshPlacesExpected method !");
  if(_mesh)
    return _type->getNumberOfMeshPlaces(_mesh);
  else
    throw INTERP_KERNEL::Exception("MEDCouplingField::getNumberOfMeshPlacesExpected : Empty mesh !");
}

/*!
 * Copy tiny info (component names, name, description) but warning the underlying mesh is not renamed (for safety reason).
 */
void MEDCouplingField::copyTinyStringsFrom(const MEDCouplingField *other)
{
  if(other)
    {
      setName(other->_name);
      setDescription(other->_desc);    
    }
}

/*!
 * This method computes the number of tuples a DataArrayDouble instance should have to build a correct MEDCouplingFieldDouble instance starting from a 
 * submesh of a virtual mesh on which a substraction per type had been applied regarding the spatial discretization in \a this.
 * 
 * For spatial discretization \b not equal to ON_GAUSS_NE this method will make the hypothesis that any positive entity id in \a code \a idsPerType is valid.
 * So in those cases attribute \a _mesh of \a this is ignored.
 * 
 * For spatial discretization equal to ON_GAUSS_NE \a _mesh attribute will be taken into account.
 *
 * The input code is those implemented in MEDCouplingUMesh::splitProfilePerType.
 *
 * \param [in] code - a code with format described above.
 * \param [in] idsPerType - a list of subparts
 * \throw If \a this has no spatial discretization set.
 * \throw If input code point to invalid zones in spatial discretization.
 * \throw If spatial discretization in \a this requires a mesh and those mesh is invalid (null,...)
 */
int MEDCouplingField::getNumberOfTuplesExpectedRegardingCode(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const
{
  const MEDCouplingFieldDiscretization *t(_type);
  if(!t)
    throw INTERP_KERNEL::Exception("MEDCouplingField::getNumberOfTuplesExpectedRegardingCode : no spatial discretization set !");
  return t->getNumberOfTuplesExpectedRegardingCode(code,idsPerType);
}
