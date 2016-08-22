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
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingMesh.hxx"

using namespace MEDCoupling;

MEDCouplingFieldInt *MEDCouplingFieldInt::New(TypeOfField type, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldInt(type,td);
}

MEDCouplingFieldInt *MEDCouplingFieldInt::New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldInt(ft,td);
}

void MEDCouplingFieldInt::checkConsistencyLight() const
{
  MEDCouplingField::checkConsistencyLight();
  _time_discr->checkConsistencyLight();
  _type->checkCoherencyBetween(_mesh,getArray());
}

std::string MEDCouplingFieldInt::simpleRepr() const
{
  std::ostringstream ret;
  ret << "FieldInt with name : \"" << getName() << "\"\n";
  ret << "Description of field is : \"" << getDescription() << "\"\n";
  if(_type)
    { ret << "FieldInt space discretization is : " << _type->getStringRepr() << "\n"; }
  else
    { ret << "FieldInt has no spatial discretization !\n"; }
  if(_time_discr)
    { ret << "FieldInt time discretization is : " << _time_discr->getStringRepr() << "\n"; }
  else
    { ret << "FieldInt has no time discretization !\n"; }
  ret << "FieldInt nature of field is : \"" << MEDCouplingNatureOfField::GetReprNoThrow(_nature) << "\"\n";
  if(getArray())
    {
      if(getArray()->isAllocated())
        {
          int nbOfCompo=getArray()->getNumberOfComponents();
          ret << "FieldInt default array has " << nbOfCompo << " components and " << getArray()->getNumberOfTuples() << " tuples.\n";
          ret << "FieldInt default array has following info on components : ";
          for(int i=0;i<nbOfCompo;i++)
            ret << "\"" << getArray()->getInfoOnComponent(i) << "\" ";
          ret << "\n";
        }
      else
        {
          ret << "Array set but not allocated !\n";
        }
    }
  if(_mesh)
    ret << "Mesh support information :\n__________________________\n" << _mesh->simpleRepr();
  else
    ret << "Mesh support information : No mesh set !\n";
  return ret.str();
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
  _time_discr->setArray(array,this);
}

const DataArrayInt *MEDCouplingFieldInt::getArray() const
{
  return _time_discr->getArray();
}

DataArrayInt *MEDCouplingFieldInt::getArray()
{
  return _time_discr->getArray();
}

MEDCouplingFieldInt::MEDCouplingFieldInt(TypeOfField type, TypeOfTimeDiscretization td):MEDCouplingField(type),_time_discr(MEDCouplingTimeDiscretizationInt::New(td))
{
}

MEDCouplingFieldInt::MEDCouplingFieldInt(const MEDCouplingFieldInt& other, bool deepCopy):MEDCouplingField(other,deepCopy),_time_discr(dynamic_cast<MEDCouplingTimeDiscretizationInt *>(other._time_discr->performCopyOrIncrRef(deepCopy)))
{
}

MEDCouplingFieldInt::MEDCouplingFieldInt(NatureOfField n, MEDCouplingTimeDiscretizationInt *td, MEDCouplingFieldDiscretization *type):MEDCouplingField(type,n),_time_discr(td)
{
}

MEDCouplingFieldInt::~MEDCouplingFieldInt()
{
  delete _time_discr;
}

/*!
 * ** WARINING : This method do not deeply copy neither mesh nor spatial discretization. Only a shallow copy (reference) is done for mesh and spatial discretization ! **
 */
MEDCouplingFieldInt::MEDCouplingFieldInt(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td):MEDCouplingField(ft,false),_time_discr(MEDCouplingTimeDiscretizationInt::New(td))
{
}
