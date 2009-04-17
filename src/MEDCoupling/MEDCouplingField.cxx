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
  if(!_type->isEqual(other->_type))
    return false;
  if(_mesh==0 && other->_mesh==0)
    return true;
  if(_mesh==0 || other->_mesh==0)
    return false;
  if(_mesh==other->_mesh)
    return true;
  return _mesh->isEqual(other->_mesh,meshPrec);
}

void MEDCouplingField::updateTime()
{
  if(_mesh)
    updateTimeWith(*_mesh);
}

TypeOfField MEDCouplingField::getEntity() const
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

MEDCouplingField::~MEDCouplingField()
{
  if(_mesh)
    ((MEDCouplingMesh *)_mesh)->decrRef();
  delete _type;
}

MEDCouplingField::MEDCouplingField(TypeOfField type):_mesh(0),_type(MEDCouplingFieldDiscretization::New(type))
{
}

MEDCouplingField::MEDCouplingField(const MEDCouplingField& other):_name(other._name),_desc(other._name),
                                                                  _mesh(0),_type(other._type->clone())
{
  if(other._mesh)
    {
      _mesh=other._mesh;
      _mesh->incrRef();
    }
}
