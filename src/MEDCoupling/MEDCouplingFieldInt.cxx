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
// Author : Yann Pora (EDF R&D)

#include "MEDCouplingFieldInt.hxx"

using namespace MEDCoupling;

MEDCouplingFieldInt *MEDCouplingFieldInt::New(TypeOfField type, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldInt(type,td);
}

void MEDCouplingFieldInt::checkConsistencyLight() const
{
  MEDCouplingField::checkConsistencyLight();
  if(_array.isNull())
    throw INTERP_KERNEL::Exception("MEDCouplingFieldInt::checkConsistencyLight : array is null !");
  _type->checkCoherencyBetween(_mesh,getArray());
}

std::string MEDCouplingFieldInt::simpleRepr() const
{
  return std::string();
}

void MEDCouplingFieldInt::reprQuickOverview(std::ostream& stream) const
{
}

void MEDCouplingFieldInt::setTimeUnit(const std::string& unit)
{
  _time_discr->setTimeUnit(unit);
}

std::string MEDCouplingFieldInt::getTimeUnit() const
{
  return _time_discr->getTimeUnit();
}

void MEDCouplingFieldInt::setTime(double val, int iteration, int order) 
{ 
  _time_discr->setTime(val,iteration,order); 
}

double MEDCouplingFieldInt::getTime(int& iteration, int& order) const
{
  return _time_discr->getTime(iteration,order);
}

void MEDCouplingFieldInt::setArray(DataArrayInt *array)
{
  MCAuto<DataArrayInt> array2(array);
  if(array2.isNotNull())
    array2->incrRef();
  _array=array2;
  //_time_discr->setArray(array,this);
}

const DataArrayInt *MEDCouplingFieldInt::getArray() const
{
  return _array;
}

DataArrayInt *MEDCouplingFieldInt::getArray()
{
  return _array;
}

MEDCouplingFieldInt::MEDCouplingFieldInt(TypeOfField type, TypeOfTimeDiscretization td):MEDCouplingField(type),_time_discr(MEDCouplingTimeDiscretization::New(td))
{
}

MEDCouplingFieldInt::MEDCouplingFieldInt(const MEDCouplingFieldInt& other, bool deepCopy):MEDCouplingField(other,deepCopy),_time_discr(other._time_discr->performCopyOrIncrRef(deepCopy))
{
  if(other._array.isNull())
    return ;
  if(deepCopy)
    {
      _array=other._array->deepCopy();
    }
  else
    {
      _array=other._array;
    }
}

MEDCouplingFieldInt::MEDCouplingFieldInt(NatureOfField n, MEDCouplingTimeDiscretization *td, MEDCouplingFieldDiscretization *type):MEDCouplingField(type,n),_time_discr(td)
{
}

MEDCouplingFieldInt::~MEDCouplingFieldInt()
{
  delete _time_discr;
}
