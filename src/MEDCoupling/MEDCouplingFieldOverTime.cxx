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

#include "MEDCouplingFieldOverTime.hxx"
#include "MEDCouplingMesh.hxx"

#include <cmath>

using namespace MEDCoupling;

MEDCouplingFieldOverTime *MEDCouplingFieldOverTime::New(const std::vector<MEDCouplingFieldDouble *>& fs)
{
  return new MEDCouplingFieldOverTime(fs);
}

double MEDCouplingFieldOverTime::getTimeTolerance() const
{
  std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();
  if(_fs.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingFieldOverTime::getTimeTolerance : empty set !");
  for(;it!=_fs.end();it++)
    if((const MEDCouplingFieldDouble *)(*it)!=0)
      return (*it)->getTimeTolerance();
  throw INTERP_KERNEL::Exception("MEDCouplingFieldOverTime::getTimeTolerance : only empty fields in this !");
}

void MEDCouplingFieldOverTime::checkConsistencyLight() const
{
  MEDCouplingMultiFields::checkConsistencyLight();
  std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();
  for(;it!=_fs.end();it++)
    if((*it)->getTimeDiscretization()==NO_TIME)
      {
        std::ostringstream oss; oss << "MEDCouplingFieldOverTime::checkConsistencyLight : At rank #" << std::distance(_fs.begin(),it) << " the field has no time !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  if(_fs.empty())
    return ;
  it=_fs.begin();
  const MEDCouplingFieldDouble& ref=*(*(it++));
  int tt1,tt2;
  double reft=ref.getEndTime(tt1,tt2);
  double eps=getTimeTolerance();
  int id=1;
  for(;it!=_fs.end();it++,id++)
    {
      if(!ref.getMesh()->areCompatibleForMerge((*it)->getMesh()))
        {
          std::ostringstream oss; oss << "Field slice at rank #" << id << " is not compatible with the first !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      double curt=(*it)->getStartTime(tt1,tt2);
      if(curt<reft-eps)
        throw INTERP_KERNEL::Exception("MEDCouplingFieldOverTime::checkConsistencyLight : fields are NOT sorted properly in ascending time !");
      reft=(*it)->getEndTime(tt1,tt2);
    }
}

std::string MEDCouplingFieldOverTime::simpleRepr() const
{
  std::ostringstream ret;
  ret << "MEDCouplingFieldOverTime with name : \"" << getName() << "\"\n";
  ret << "Description of MEDCouplingFieldOverTime is : \"" << getDescription() << "\"\n";
  ret << "Number of discretization : " << _fs.size() << "\n";
  ret << "Number of different meshes : ";
  std::vector<MEDCouplingMesh *> ms;
  std::vector<int> refms;
  try
  {
      ms=getDifferentMeshes(refms);
      ret << ms.size() << "\n";
  }
  catch(INTERP_KERNEL::Exception& /*e*/)
  { ret << "Current instance is INVALID !\n"; }
  try
  {
      MEDCouplingDefinitionTime dt=getDefinitionTimeZone();
      dt.appendRepr(ret);
  }
  catch(INTERP_KERNEL::Exception& /*e*/)
  { ret << "Definition zone is INVALID !\n"; }
  return ret.str();
}

bool MEDCouplingFieldOverTime::isEqual(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const
{
  if(!MEDCouplingMultiFields::isEqual(other,meshPrec,valsPrec))
    return false;
  const MEDCouplingFieldOverTime *otherC=dynamic_cast<const MEDCouplingFieldOverTime *>(other);
  if(!otherC)
    return false;
  // to implement
  return true;
}

bool MEDCouplingFieldOverTime::isEqualWithoutConsideringStr(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const
{
  if(!MEDCouplingMultiFields::isEqualWithoutConsideringStr(other,meshPrec,valsPrec))
    return false;
  const MEDCouplingFieldOverTime *otherC=dynamic_cast<const MEDCouplingFieldOverTime *>(other);
  if(!otherC)
    return false;
  // to implement
  return true;
}

std::vector<MEDCouplingMesh *> MEDCouplingFieldOverTime::getMeshes() const
{
  checkConsistencyLight();
  return MEDCouplingMultiFields::getMeshes();
}

std::vector<MEDCouplingMesh *> MEDCouplingFieldOverTime::getDifferentMeshes(std::vector<int>& refs) const
{
  checkConsistencyLight();
  return MEDCouplingMultiFields::getDifferentMeshes(refs);
}

std::vector<DataArrayDouble *> MEDCouplingFieldOverTime::getArrays() const
{
  checkConsistencyLight();
  return MEDCouplingMultiFields::getArrays();
}

std::vector<DataArrayDouble *> MEDCouplingFieldOverTime::getDifferentArrays(std::vector< std::vector<int> >& refs) const
{
  checkConsistencyLight();
  return MEDCouplingMultiFields::getDifferentArrays(refs);
}

MEDCouplingDefinitionTime MEDCouplingFieldOverTime::getDefinitionTimeZone() const
{
  std::vector< std::vector<int> > tmp;
  getDifferentArrays(tmp);
  std::vector<const MEDCouplingFieldDouble *> tmp2(_fs.begin(),_fs.end());
  std::vector<int> tmp3;
  getDifferentMeshes(tmp3);
  return MEDCouplingDefinitionTime(tmp2,tmp3,tmp);
}

MEDCouplingFieldOverTime::MEDCouplingFieldOverTime(const std::vector<MEDCouplingFieldDouble *>& fs):MEDCouplingMultiFields(fs)
{
  checkConsistencyLight();
}

MEDCouplingFieldOverTime::MEDCouplingFieldOverTime()
{
}
