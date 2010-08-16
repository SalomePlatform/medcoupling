//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#include "MEDCouplingField.hxx"
#include "MEDCouplingMesh.hxx"
#include "MEDCouplingFieldDiscretization.hxx"

using namespace ParaMEDMEM;

bool MEDCouplingField::isEqual(const MEDCouplingField *other, double meshPrec, double valsPrec) const
{
  if(_name!=other->_name)
    return false;
  if(_desc!=other->_desc)
    return false;
  if(!_type->isEqual(other->_type,valsPrec))
    return false;
  if(_mesh==0 && other->_mesh==0)
    return true;
  if(_mesh==0 || other->_mesh==0)
    return false;
  if(_mesh==other->_mesh)
    return true;
  return _mesh->isEqual(other->_mesh,meshPrec);
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
  if(_mesh==other->_mesh)
    return true;
  return _mesh->areCompatibleForMerge(other->_mesh);
}

/*!
 * This method is more strict than MEDCouplingField::areCompatible method.
 * This method is used for operation on fields to operate a first check before attempting operation.
 */
bool MEDCouplingField::areStrictlyCompatible(const MEDCouplingField *other) const
{
  if(!_type->isEqual(other->_type,1.e-12))
    return false;
  return _mesh==other->_mesh;
}

void MEDCouplingField::updateTime()
{
  if(_mesh)
    updateTimeWith(*_mesh);
  if(_type)
    updateTimeWith(*_type);
}

TypeOfField MEDCouplingField::getTypeOfField() const
{
  return _type->getEnum();
}

void MEDCouplingField::setMesh(const MEDCouplingMesh *mesh)
{
  if(mesh!=_mesh)
    {
      if(_mesh)
        ((MEDCouplingMesh *)_mesh)->decrRef();
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
    ((MEDCouplingMesh *)_mesh)->decrRef();
  delete _type;
}

MEDCouplingField::MEDCouplingField(TypeOfField type):_mesh(0),_type(MEDCouplingFieldDiscretization::New(type))
{
}

MEDCouplingField::MEDCouplingField(const MEDCouplingField& other):_name(other._name),_desc(other._desc),
                                                                  _mesh(0),_type(other._type->clone())
{
  if(other._mesh)
    {
      _mesh=other._mesh;
      _mesh->incrRef();
    }
}

/*!
 * This method returns a submesh of 'mesh' instance constituting cell ids contained in array defined as an interval [start;end).
 * @ param di is an array returned that specifies entity ids (nodes, cells ids...) in mesh 'mesh' of entity in returned submesh.
 */
MEDCouplingMesh *MEDCouplingField::buildSubMeshData(const int *start, const int *end, DataArrayInt *&di) const
{
  return _type->buildSubMeshData(start,end,_mesh,di);
}
