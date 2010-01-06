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
#include "MEDCouplingMemArray.txx"

using namespace ParaMEDMEM;

void DataArray::setName(const char *name)
{
  _name=name;
}

void DataArray::copyStringInfoFrom(const DataArray& other) throw(INTERP_KERNEL::Exception)
{
  if(_info_on_compo.size()!=other._info_on_compo.size())
    throw INTERP_KERNEL::Exception("Size of arrays mismatches on copyStringInfoFrom !");
  _name=other._name;
  _info_on_compo=other._info_on_compo;
}

bool DataArray::areInfoEquals(const DataArray& other) const
{
  if(_nb_of_tuples!=other._nb_of_tuples)
    return false;
  if(_name!=other._name)
    return false;
  return _info_on_compo==other._info_on_compo;
}

DataArrayDouble *DataArrayDouble::New()
{
  return new DataArrayDouble;
}

DataArrayDouble *DataArrayDouble::deepCopy() const
{
  return new DataArrayDouble(*this);
}

DataArrayDouble *DataArrayDouble::performCpy(bool deepCpy) const
{
  if(deepCpy)
    return deepCopy();
  else
    {
      incrRef();
      return const_cast<DataArrayDouble *>(this);
    }
}

void DataArrayDouble::alloc(int nbOfTuple, int nbOfCompo)
{
  _nb_of_tuples=nbOfTuple;
  _info_on_compo.resize(nbOfCompo);
  _mem.alloc(nbOfCompo*_nb_of_tuples);
  declareAsNew();
}

bool DataArrayDouble::isEqual(const DataArrayDouble& other, double prec) const
{
  if(!areInfoEquals(other))
    return false;
  return _mem.isEqual(other._mem,prec);
}

void DataArrayDouble::reAlloc(int nbOfTuples)
{
  _mem.reAlloc(_info_on_compo.size()*nbOfTuples);
  _nb_of_tuples=nbOfTuples;
  declareAsNew();
}

void DataArrayDouble::setArrayIn(DataArrayDouble *newArray, DataArrayDouble* &arrayToSet)
{
  if(newArray!=arrayToSet)
    {
      if(arrayToSet)
        arrayToSet->decrRef();
      arrayToSet=newArray;
      if(arrayToSet)
        arrayToSet->incrRef();
    }
}

void DataArrayDouble::useArray(const double *array, bool ownership,  DeallocType type, int nbOfTuple, int nbOfCompo)
{
  _nb_of_tuples=nbOfTuple;
  _info_on_compo.resize(nbOfCompo);
  _mem.useArray(array,ownership,type,nbOfTuple*nbOfCompo);
  declareAsNew();
}

DataArrayInt *DataArrayInt::New()
{
  return new DataArrayInt;
}

DataArrayInt *DataArrayInt::deepCopy() const
{
  return new DataArrayInt(*this);
}

DataArrayInt *DataArrayInt::performCpy(bool deepCpy) const
{
  if(deepCpy)
    return deepCopy();
  else
    {
      incrRef();
      return const_cast<DataArrayInt *>(this);
    }
}

void DataArrayInt::alloc(int nbOfTuple, int nbOfCompo)
{
  _nb_of_tuples=nbOfTuple;
  _info_on_compo.resize(nbOfCompo);
  _mem.alloc(nbOfCompo*_nb_of_tuples);
  declareAsNew();
}

bool DataArrayInt::isEqual(const DataArrayInt& other) const
{
  if(!areInfoEquals(other))
    return false;
  return _mem.isEqual(other._mem,0);
}

void DataArrayInt::useArray(const int *array, bool ownership,  DeallocType type, int nbOfTuple, int nbOfCompo)
{
  _nb_of_tuples=nbOfTuple;
  _info_on_compo.resize(nbOfCompo);
  _mem.useArray(array,ownership,type,nbOfTuple*nbOfCompo);
  declareAsNew();
}

void DataArrayInt::reAlloc(int nbOfTuples)
{
  _mem.reAlloc(_info_on_compo.size()*nbOfTuples);
  _nb_of_tuples=nbOfTuples;
  declareAsNew();
}

void DataArrayInt::setArrayIn(DataArrayInt *newArray, DataArrayInt* &arrayToSet)
{
  if(newArray!=arrayToSet)
    {
      if(arrayToSet)
        arrayToSet->decrRef();
      arrayToSet=newArray;
      if(arrayToSet)
        arrayToSet->incrRef();
    }
}
