//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#include "MEDCouplingRMesh.hxx"
#include "MemArray.hxx"

using namespace ParaMEDMEM;

MEDCouplingRMesh::MEDCouplingRMesh():_x_array(0),_y_array(0),_z_array(0)
{
}

MEDCouplingRMesh::~MEDCouplingRMesh()
{
  if(_x_array)
    _x_array->decrRef();
  if(_y_array)
    _y_array->decrRef();
  if(_z_array)
    _z_array->decrRef();
}

MEDCouplingRMesh *MEDCouplingRMesh::New()
{
  return new MEDCouplingRMesh;
}

void MEDCouplingRMesh::updateTime()
{
  if(_x_array)
    updateTimeWith(*_x_array);
  if(_y_array)
    updateTimeWith(*_y_array);
  if(_z_array)
    updateTimeWith(*_z_array);
}

bool MEDCouplingRMesh::isEqual(const MEDCouplingMesh *other, double prec) const
{
  const MEDCouplingRMesh *otherC=dynamic_cast<const MEDCouplingRMesh *>(other);
  if(!otherC)
    return false;
  return true;
}

void MEDCouplingRMesh::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  const char msg0[]="Invalid ";
  const char msg1[]=" array ! Must contain more than 1 element.";
  if(_x_array)
    if(_x_array->getNbOfElems()<2)
      {
        std::ostringstream os; os << msg0 << 'X' << msg1;
        throw INTERP_KERNEL::Exception(os.str().c_str());
      }
  if(_y_array)
    if(_y_array->getNbOfElems()<2)
      {
        std::ostringstream os; os << msg0 << 'Y' << msg1;
        throw INTERP_KERNEL::Exception(os.str().c_str());
      }
  if(_z_array)
    if(_z_array->getNbOfElems()<2)
      {
        std::ostringstream os; os << msg0 << 'Z' << msg1;
        throw INTERP_KERNEL::Exception(os.str().c_str());
      }
}

bool MEDCouplingRMesh::isStructured() const
{
  return true;
}

int MEDCouplingRMesh::getNumberOfCells() const
{
  int ret=1;
  if(_x_array)
    ret*=_x_array->getNbOfElems()-1;
  if(_y_array)
    ret*=_y_array->getNbOfElems()-1;
  if(_z_array)
    ret*=_z_array->getNbOfElems()-1;
  return ret;
}

int MEDCouplingRMesh::getNumberOfNodes() const
{
  int ret=1;
  if(_x_array)
    ret*=_x_array->getNbOfElems();
  if(_y_array)
    ret*=_y_array->getNbOfElems();
  if(_z_array)
    ret*=_z_array->getNbOfElems();
  return ret;
}

int MEDCouplingRMesh::getSpaceDimension() const
{
  int ret=0;
  if(_x_array)
    ret++;
  if(_y_array)
    ret++;
  if(_z_array)
    ret++;
  return ret;
}

int MEDCouplingRMesh::getMeshDimension() const
{
  int ret=0;
  if(_x_array)
    ret++;
  if(_y_array)
    ret++;
  if(_z_array)
    ret++;
  return ret;
}

DataArrayDouble *MEDCouplingRMesh::getCoordsAt(int i) const throw(INTERP_KERNEL::Exception)
{
  switch(i)
    {
    case 0:
      return _x_array;
    case 1:
      return _y_array;
    case 2:
      return _z_array;
    default:
      throw INTERP_KERNEL::Exception("Invalid rank specified must be 0 or 1 or 2.");
    }
}

void MEDCouplingRMesh::setCoords(DataArrayDouble *coordsX, DataArrayDouble *coordsY, DataArrayDouble *coordsZ)
{
  if(_x_array)
    _x_array->decrRef();
  _x_array=coordsX;
  if(_x_array)
    _x_array->incrRef();
  if(_y_array)
    _y_array->decrRef();
  _y_array=coordsY;
  if(_y_array)
    _y_array->incrRef();
  if(_z_array)
    _z_array->decrRef();
  _z_array=coordsZ;
  if(_z_array)
    _z_array->incrRef();
  declareAsNew();
}

void MEDCouplingRMesh::getBoundingBox(double *bbox) const
{
  //not implemented yet !
}

MEDCouplingFieldDouble *MEDCouplingRMesh::getMeasureField() const
{
  //not implemented yet !
}
