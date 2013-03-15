// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

#include "MEDCouplingField.hxx"
#include "MEDCouplingMesh.hxx"
#include "MEDCouplingFieldDiscretization.hxx"

#include <sstream>

using namespace ParaMEDMEM;

bool MEDCouplingField::isEqualIfNotWhy(const MEDCouplingField *other, double meshPrec, double valsPrec, std::string& reason) const throw(INTERP_KERNEL::Exception)
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
      oss << "Field nature differ : this nature = \"" << MEDCouplingNatureOfField::getRepr(_nature) << "\" and other nature = \"" << MEDCouplingNatureOfField::getRepr(other->_nature) << "\" !";
      reason=oss.str();
      return false;
    }
  if(!_type->isEqualIfNotWhy(other->_type,valsPrec,reason))
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

bool MEDCouplingField::isEqual(const MEDCouplingField *other, double meshPrec, double valsPrec) const
{
  std::string tmp;
  return isEqualIfNotWhy(other,meshPrec,valsPrec,tmp);
}

bool MEDCouplingField::isEqualWithoutConsideringStr(const MEDCouplingField *other, double meshPrec, double valsPrec) const
{
  if(!_type->isEqualWithoutConsideringStr(other->_type,valsPrec))
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
  if(!_type->isEqual(other->_type,1.e-12))
    return false;
  if(_nature!=other->_nature)
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

std::size_t MEDCouplingField::getHeapMemorySize() const
{
  std::size_t ret=0;
  ret+=_name.capacity();
  ret+=_desc.capacity();
  if(_mesh)
    ret+=_mesh->getHeapMemorySize();
  if((const MEDCouplingFieldDiscretization *)_type)
    ret+=_type->getHeapMemorySize();
  return ret;
}

TypeOfField MEDCouplingField::getTypeOfField() const
{
  return _type->getEnum();
}

/*!
 * This method returns the nature of field. This information is very important during interpolation process using ParaMEDMEM::MEDCouplingRemapper or ParaMEDMEM::InterpKernelDEC.
 * In other context than the two mentioned before this attribute of the field is not sensitive. This attribute is not store in MED file in MEDLoader.
 * More information of the semantic, and the consequence of this attribute in the result of the interpolation, is available \ref NatureOfField "here".
 */
NatureOfField MEDCouplingField::getNature() const
{
  return _nature;
}

/*!
 * This method set the nature of field in \b this.This  information is very important during interpolation process using ParaMEDMEM::MEDCouplingRemapper or ParaMEDMEM::InterpKernelDEC.
 * In other context than the two mentioned before this attribute of the field is not sensitive. This attribute is not store in MED file in MEDLoader.
 * More information of the semantic, and the consequence of this attribute in the result of the interpolation, is available \ref TableNatureOfField "here".
 */
void MEDCouplingField::setNature(NatureOfField nat) throw(INTERP_KERNEL::Exception)
{
  _nature=nat;
}

/*!
 * This method returns is case of success an instance of DataArrayDouble the user is in reponsability to deal with.
 * If 'this->_mesh' is not set an exception will be thrown.
 * For a field on node the array of coords will be returned. For a field on cell a ParaMEDMEM::DataArrayDouble instance
 * containing the barycenter of cells will be returned. And for a field on gauss point the explicit position of gauss points.
 */
DataArrayDouble *MEDCouplingField::getLocalizationOfDiscr() const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingField::getLocalizationOfDiscr : No mesh set !");
  return _type->getLocalizationOfDiscValues(_mesh);
}

/*!
 * This method retrieves the measure field of 'this'. If no '_mesh' is defined an exception will be thrown.
 * Warning the retrieved field life cycle is the responsability of caller.
 */
MEDCouplingFieldDouble *MEDCouplingField::buildMeasureField(bool isAbs) const throw(INTERP_KERNEL::Exception)
{
  if(_mesh==0)
    throw INTERP_KERNEL::Exception("MEDCouplingField::getMeasureField : no mesh defined !!!");
  return _type->getMeasureField(_mesh,isAbs);
}

void MEDCouplingField::setMesh(const MEDCouplingMesh *mesh)
{
  if(mesh!=_mesh)
    {
      if(_mesh)
        _mesh->decrRef();
      _mesh=mesh;
      if(_mesh)
        {
          _mesh->incrRef();
          updateTimeWith(*_mesh);
        }
    }
}

/*!
 * This method sets gauss localization by geometric type.
 * @param type geometric type on which the gauss localization will be set.
 * @param refCoo is the reference coordinates of the specified element. Its size has to be equal to nbOfNodesPerCell*dimOfType
 * @param gsCoo are the coordinates of Gauss points in reference element specified by 'refCoo'. Its size must be equal to wg.size()*dimOfType
 * @param wg are the weights on Gauss points. The size of this array is used to determine the number of Gauss point in the element.
 * @throw when size of 'RefCoo' is not valid regarding 'type' parameter, it throws too when the mesh is not set before or if it is not a field on Gauss points.
 */
void MEDCouplingField::setGaussLocalizationOnType(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                  const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Mesh has to be set before calling setGaussLocalizationOnType method !");
  _type->setGaussLocalizationOnType(_mesh,type,refCoo,gsCoo,wg);
}

/*!
 * This method sets on ids defined by [begin;end) their gauss localization. This method checks the coherency of cells ids in [begin;end) and 'refCoo' size.
 * If an incoherence appears an exception will be thrown and no seting will be performed.
 * An exception is thrown too if [begin,end) has a size lesser than 1.
 * 
 * @param refCoo is the reference coordinates of the specified element. Its size has to be equal to nbOfNodesPerCell*dimOfType
 * @param gsCoo are the coordinates of Gauss points in reference element specified by 'refCoo'. Its size must be equal to wg.size()*dimOfType
 * @param wg are the weights on Gauss points. The size of this array is used to determine the number of Gauss point in the element.
 * @throw when size of 'RefCoo' is not valid regarding cells in [begin,end) parameters, it throws too when the mesh is not set before or if it is not a field on Gauss points.
 */
void MEDCouplingField::setGaussLocalizationOnCells(const int *begin, const int *end, const std::vector<double>& refCoo,
                                                   const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Mesh has to be set before calling setGaussLocalizationOnCells method !");
  _type->setGaussLocalizationOnCells(_mesh,begin,end,refCoo,gsCoo,wg);
}

/*!
 * This method resets all Gauss loalizations if any.
 */
void MEDCouplingField::clearGaussLocalizations()
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Mesh has to be set before calling clearGaussLocalizations method !");
  _type->clearGaussLocalizations();
}

/*!
 * This method returns reference to the Gauss localization object corresponding to 'locId' id.
 * This method throws an exception if there is no mesh, invalid FieldDescription (different from Gauss) and if 'locId' is invalid because out of range given by
 * MEDCouplingField::getNbOfGaussLocalization method.
 * Warning this method is not const, so the returned object could be modified without any problem.
 */
MEDCouplingGaussLocalization& MEDCouplingField::getGaussLocalization(int locId) throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Mesh has to be set before calling getGaussLocalization method !");
  return _type->getGaussLocalization(locId);
}

/*!
 * This method returns reference to the Gauss localization object corresponding to 'locId' id.
 * This method throws an exception if there is no mesh, invalid FieldDescription (different from Gauss) and if several localization ids have been found
 * for a type.
 */
int MEDCouplingField::getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Mesh has to be set before calling getGaussLocalizationIdOfOneType method !");
  return _type->getGaussLocalizationIdOfOneType(type);
}

std::set<int> MEDCouplingField::getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Mesh has to be set before calling getGaussLocalizationIdsOfOneType method !");
  return _type->getGaussLocalizationIdsOfOneType(type);
}

/*!
 * This method returns number of Gauss localization available. Implicitely all ids in [0,getNbOfGaussLocalization()) is a valid Gauss localisation id.
 * This method throws an exception if there is no mesh, invalid FieldDescription (different from Gauss)
 */
int MEDCouplingField::getNbOfGaussLocalization() const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Mesh has to be set before calling getNbOfGaussLocalization method !");
  return _type->getNbOfGaussLocalization();
}

/*!
 * This method returns an id of Gauss localization in [0,getNbOfGaussLocalization()) that corresponds to the localization of the cell specified by its cellId.
 * This methods throws an exception if there is no mesh, invalid FieldDescription (different from Gauss) or if at the cell with id 'cellId' in this->_mesh no
 * Gauss localization has been set.
 */
int MEDCouplingField::getGaussLocalizationIdOfOneCell(int cellId) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Mesh has to be set before calling getGaussLocalizationIdOfOneCell method !");
  return _type->getGaussLocalizationIdOfOneCell(cellId);
}

/*!
 * This method returns all cellIds that share the same Gauss localization given by 'locId' parameter (in range [0,getNbOfGaussLocalization()) ).
 * If no cells fit the Gauss localization given by 'locId' cellIds will be returned empty.
 * @param locId input that specifies the id of Gauss localization.
 * @param cellIds output parameter, that will contain the result if this method succeds. This parameter is systematically cleared when called.
 * @throw  if there is no mesh, invalid FieldDescription (different from Gauss) or if locId not in [0,getNbOfGaussLocalization())
 */
void MEDCouplingField::getCellIdsHavingGaussLocalization(int locId, std::vector<int>& cellIds) const throw(INTERP_KERNEL::Exception)
{
  cellIds.clear();
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Mesh has to be set before calling getGaussLocalizationIdOfOneCell method !");
  _type->getCellIdsHavingGaussLocalization(locId,cellIds);
}

/*!
 * This method returns reference to the Gauss localization object corresponding to 'locId' id.
 * This method throws an exception if there is no mesh, invalid FieldDescription (different from Gauss) and if 'locId' is invalid because out of range given by
 * MEDCouplingField::getNbOfGaussLocalization method.
 * Warning this method is const.
 */
const MEDCouplingGaussLocalization& MEDCouplingField::getGaussLocalization(int locId) const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Mesh has to be set before calling getGaussLocalization method !");
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
 * This method returns a submesh of 'mesh' instance constituting cell ids contained in array defined as an interval [start;end).
 * @param di is an array returned that specifies entity ids (nodes, cells ids...) in mesh 'mesh' of entity in returned submesh.
 */
MEDCouplingMesh *MEDCouplingField::buildSubMeshData(const int *start, const int *end, DataArrayInt *&di) const
{
  return _type->buildSubMeshData(_mesh,start,end,di);
}

/*!
 * This method returns tuples ids implied by the mesh selection of the  cell ids contained in array defined as an interval [start;end).
 * \return a newly allocated DataArrayInt instance containing tuples ids.
 */
DataArrayInt *MEDCouplingField::computeTupleIdsToSelectFromCellIds(const int *startCellIds, const int *endCellIds) const
{
  return _type->computeTupleIdsToSelectFromCellIds(_mesh,startCellIds,endCellIds);
}

/*!
 * This method returns number of tuples expected regarding its discretization and its _mesh attribute.
 * This method expected a not null _mesh instance. If null, an exception will be thrown.
 */
int MEDCouplingField::getNumberOfTuplesExpected() const throw(INTERP_KERNEL::Exception)
{
  if(_mesh)
    return _type->getNumberOfTuples(_mesh);
  else
    throw INTERP_KERNEL::Exception("MEDCouplingField::getNumberOfTuplesExpected : Empty mesh !");
}

/*!
 * This method returns number of mesh placed expected regarding its discretization and its _mesh attribute.
 * This method expected a not null _mesh instance. If null, an exception will be thrown.
 */
int MEDCouplingField::getNumberOfMeshPlacesExpected() const throw(INTERP_KERNEL::Exception)
{
  if(_mesh)
    return _type->getNumberOfMeshPlaces(_mesh);
  else
    throw INTERP_KERNEL::Exception("MEDCouplingField::getNumberOfMeshPlacesExpected : Empty mesh !");
}

/*!
 * Copy tiny info (component names, name, description) but warning the underlying mesh is not renamed (for safety reason).
 */
void MEDCouplingField::copyTinyStringsFrom(const MEDCouplingField *other) throw(INTERP_KERNEL::Exception)
{
  if(other)
    {
      setName(other->_name.c_str());
      setDescription(other->_desc.c_str());    
    }
}
