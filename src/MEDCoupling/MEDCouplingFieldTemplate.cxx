// Copyright (C) 2007-2024  CEA, EDF
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
// Author : Anthony Geay (EDF R&D)

#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingMesh.hxx"
#include "MEDCouplingFieldInt32.hxx"
#include "MEDCouplingFieldInt64.hxx"
#include "MEDCouplingFieldFloat.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldDiscretization.hxx"

#include <sstream>

using namespace MEDCoupling;

MEDCouplingFieldTemplate *MEDCouplingFieldTemplate::New(const MEDCouplingFieldDouble& f)
{
  return new MEDCouplingFieldTemplate(f,true);
}

MEDCouplingFieldTemplate *MEDCouplingFieldTemplate::New(const MEDCouplingFieldFloat& f)
{
  return new MEDCouplingFieldTemplate(f,true);
}

MEDCouplingFieldTemplate *MEDCouplingFieldTemplate::New(const MEDCouplingFieldInt32& f)
{
  return new MEDCouplingFieldTemplate(f,true);
}

MEDCouplingFieldTemplate *MEDCouplingFieldTemplate::New(const MEDCouplingFieldInt64& f)
{
  return new MEDCouplingFieldTemplate(f,true);
}

MEDCouplingFieldTemplate *MEDCouplingFieldTemplate::NewWithoutCheck(const MEDCouplingFieldDouble& f)
{
  return new MEDCouplingFieldTemplate(f,false);
}

MEDCouplingFieldTemplate *MEDCouplingFieldTemplate::NewWithoutCheck(const MEDCouplingFieldFloat& f)
{
  return new MEDCouplingFieldTemplate(f,false);
}

MEDCouplingFieldTemplate *MEDCouplingFieldTemplate::NewWithoutCheck(const MEDCouplingFieldInt32& f)
{
  return new MEDCouplingFieldTemplate(f,false);
}

MEDCouplingFieldTemplate *MEDCouplingFieldTemplate::NewWithoutCheck(const MEDCouplingFieldInt64& f)
{
  return new MEDCouplingFieldTemplate(f,false);
}

bool MEDCouplingFieldTemplate::isEqualIfNotWhy(const MEDCouplingFieldTemplate *other, double meshPrec, std::string& reason) const
{
  return isEqualIfNotWhyProtected(other,meshPrec,reason);
}

bool MEDCouplingFieldTemplate::isEqual(const MEDCouplingFieldTemplate *other, double meshPrec) const
{
  std::string tmp;
  return isEqualIfNotWhyProtected(other,meshPrec,tmp);
}

bool MEDCouplingFieldTemplate::isEqualWithoutConsideringStr(const MEDCouplingFieldTemplate *other, double meshPrec) const
{
  return isEqualWithoutConsideringStrProtected(other,meshPrec);
}

/*!
 * The user should \b not use this method. Only useful for CORBA serialization/unserialization.
 */
MEDCouplingFieldTemplate *MEDCouplingFieldTemplate::New(TypeOfField type)
{
  return new MEDCouplingFieldTemplate(type);
}

MEDCouplingFieldTemplate::MEDCouplingFieldTemplate(const MEDCouplingFieldDouble& f, bool isChecked):MEDCouplingField(f,false) 
{
  forceTimeOfThis(f);
  if(isChecked)
    checkConsistencyLight();
}

MEDCouplingFieldTemplate::MEDCouplingFieldTemplate(const MEDCouplingFieldFloat& f, bool isChecked):MEDCouplingField(f,false) 
{
  forceTimeOfThis(f);
  if(isChecked)
    checkConsistencyLight();
}

MEDCouplingFieldTemplate::MEDCouplingFieldTemplate(const MEDCouplingFieldInt32& f, bool isChecked):MEDCouplingField(f,false) 
{
  forceTimeOfThis(f);
  if(isChecked)
    checkConsistencyLight();
}

MEDCouplingFieldTemplate::MEDCouplingFieldTemplate(const MEDCouplingFieldInt64& f, bool isChecked):MEDCouplingField(f,false) 
{
  forceTimeOfThis(f);
  if(isChecked)
    checkConsistencyLight();
}

MEDCouplingFieldTemplate::MEDCouplingFieldTemplate(TypeOfField type):MEDCouplingField(type)
{
}

MEDCouplingFieldTemplate::MEDCouplingFieldTemplate(const MEDCouplingFieldTemplate& other, bool deepCopy):MEDCouplingField(other,deepCopy)
{
}

void MEDCouplingFieldTemplate::checkConsistencyLight() const
{
  if(_mesh==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldTemplate::checkConsistencyLight : Empty mesh !");
}

std::string MEDCouplingFieldTemplate::simpleRepr() const
{
  std::ostringstream ret;
  ret << "FieldTemplate with name : \"" << getName() << "\"\n";
  ret << "Description of field is : \"" << getDescription() << "\"\n";
  if(_type)
    { ret << "FieldTemplate space discretization is : " << _type->getStringRepr() << "\n"; }
  else
    { ret << "FieldTemplate has no spatial discretization !\n"; }
  ret << "FieldTemplate nature of field is : \"" << MEDCouplingNatureOfField::GetReprNoThrow(_nature) << "\"\n";
  if(_mesh)
    ret << "Mesh support information :\n__________________________\n" << _mesh->simpleRepr();
  else
    ret << "Mesh support information : No mesh set !\n";
  return ret.str();
}

std::string MEDCouplingFieldTemplate::advancedRepr() const
{
  return simpleRepr();
}

void MEDCouplingFieldTemplate::getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform getTinySerializationIntInformation !");
  tinyInfo.clear();
  tinyInfo.push_back(ToIdType(_type->getEnum()));
  tinyInfo.push_back(ToIdType(_nature));
  std::vector<mcIdType> tinyInfo2;
  _type->getTinySerializationIntInformation(tinyInfo2);
  tinyInfo.insert(tinyInfo.end(),tinyInfo2.begin(),tinyInfo2.end());
  tinyInfo.push_back(ToIdType(tinyInfo2.size()));
}

void MEDCouplingFieldTemplate::getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform getTinySerializationDbleInformation !");
  tinyInfo.clear();
  _type->getTinySerializationDbleInformation(tinyInfo);
}

void MEDCouplingFieldTemplate::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const
{
  tinyInfo.clear();
  tinyInfo.push_back(_name);
  tinyInfo.push_back(_desc);
}

void MEDCouplingFieldTemplate::resizeForUnserialization(const std::vector<mcIdType>& tinyInfoI, DataArrayIdType *&dataInt)
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform resizeForUnserialization !");
  dataInt=0;
  std::vector<mcIdType> tinyInfoITmp(tinyInfoI.begin()+2,tinyInfoI.end());
  _type->resizeForUnserialization(tinyInfoITmp,dataInt);
}

void MEDCouplingFieldTemplate::finishUnserialization(const std::vector<mcIdType>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS)
{
  if(!((const MEDCouplingFieldDiscretization *)_type))
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform finishUnserialization !");
  _nature=(NatureOfField)tinyInfoI[1];
  _type->finishUnserialization(tinyInfoD);
  _name=tinyInfoS[0];
  _desc=tinyInfoS[1];
}

void MEDCouplingFieldTemplate::serialize(DataArrayIdType *&dataInt) const
{
  _type->getSerializationIntArray(dataInt);
}

void MEDCouplingFieldTemplate::reprQuickOverview(std::ostream& stream) const
{
  stream << "MEDCouplingFieldTemplate C++ instance at " << this << ". Name : \"" << _name << "\"." << std::endl;
  const char *nat=0;
  try
  {
      nat=MEDCouplingNatureOfField::GetRepr(_nature);
      stream << "Nature of field template : " << nat << ".\n";
  }
  catch(INTERP_KERNEL::Exception& /*e*/)
  {  }
  const MEDCouplingFieldDiscretization *fd(_type);
  if(!fd)
    stream << "No spatial discretization set !";
  else
    fd->reprQuickOverview(stream);
  stream << std::endl;
  if(!_mesh)
    stream << "\nNo mesh support defined !";
  else
    {
      std::ostringstream oss;
      _mesh->reprQuickOverview(oss);
      std::string tmp(oss.str());
      stream << "\nMesh info : " << tmp.substr(0,tmp.find('\n'));
    }
}

MCAuto<MEDCouplingFieldTemplate> MEDCouplingFieldTemplate::clone(bool recDeepCpy) const
{
  MCAuto<MEDCouplingFieldTemplate> ret(new MEDCouplingFieldTemplate(*this,recDeepCpy));
  return ret;
}
