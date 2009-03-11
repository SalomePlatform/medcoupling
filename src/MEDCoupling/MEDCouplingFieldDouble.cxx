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
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMesh.hxx"

#include <sstream>

using namespace ParaMEDMEM;

MEDCouplingFieldDouble *MEDCouplingFieldDouble::New(TypeOfField type)
{
  return new MEDCouplingFieldDouble(type);
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::clone(bool recDeepCpy) const
{
  return new MEDCouplingFieldDouble(*this,recDeepCpy);
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(TypeOfField type):MEDCouplingField(type),_array(0)
{
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(const MEDCouplingFieldDouble& other, bool deepCpy):MEDCouplingField(other),_array(0)
{
  if(other._array)
    _array=other._array->performCpy(deepCpy);
}

MEDCouplingFieldDouble::~MEDCouplingFieldDouble()
{
  if(_array)
    _array->decrRef();
}

void MEDCouplingFieldDouble::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Field invalid because no mesh specified !");
  if(!_array)
    throw INTERP_KERNEL::Exception("Field invalid because no values set !");
  if(_type==ON_CELLS)
    {
      if(_mesh->getNumberOfCells()!=_array->getNumberOfTuples())
        {
          std::ostringstream message;
          message << "Field on cells invalid because there are " << _mesh->getNumberOfCells();
          message << " cells in mesh and " << _array->getNumberOfTuples() << " tuples in field !";
          throw INTERP_KERNEL::Exception(message.str().c_str());
        }
    }
  else if(_type==ON_NODES)
    {
      if(_mesh->getNumberOfNodes()!=_array->getNumberOfTuples())
        {
          std::ostringstream message;
          message << "Field on nodes invalid because there are " << _mesh->getNumberOfNodes();
          message << " cells in mesh and " << _array->getNumberOfTuples() << " tuples in field !";
          throw INTERP_KERNEL::Exception(message.str().c_str());
        }
    }
  else
    throw INTERP_KERNEL::Exception("Field of undifined type !!!");
}

void MEDCouplingFieldDouble::applyLin(double a, double b, int compoId)
{
  double *ptr=_array->getPointer();
  ptr+=compoId;
  int nbOfComp=_array->getNumberOfComponents();
  int nbOfTuple=_array->getNumberOfTuples();
  for(int i=0;i<nbOfTuple;i++,ptr+=nbOfComp)
    *ptr=a*(*ptr)+b;
}

int MEDCouplingFieldDouble::getNumberOfComponents() const
{
  return _array->getNumberOfComponents();
}

int MEDCouplingFieldDouble::getNumberOfTuples() const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Impossible to retrieve number of tuples because no mesh specified !");
  if(_type==ON_CELLS)
    return _mesh->getNumberOfCells();
  else if(_type==ON_NODES)
    return _mesh->getNumberOfNodes();
  else
    throw INTERP_KERNEL::Exception("Impossible to retrieve number of tuples because type of entity not implemented yet !");
}

void MEDCouplingFieldDouble::updateTime()
{
  MEDCouplingField::updateTime();
  if(_array)
    updateTimeWith(*_array);
}

void MEDCouplingFieldDouble::setArray(DataArrayDouble *array)
{
  if(array!=_array)
    {
      if(_array)
        _array->decrRef();
      _array=array;
      if(_array)
        _array->incrRef();
      declareAsNew();
    }
}
