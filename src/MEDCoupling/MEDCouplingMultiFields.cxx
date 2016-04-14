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

#include "MEDCouplingMultiFields.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMesh.hxx"
#include "MCAuto.hxx"

#include <sstream>
#include <algorithm>

using namespace MEDCoupling;

MEDCouplingMultiFields *MEDCouplingMultiFields::New(const std::vector<MEDCouplingFieldDouble *>& fs)
{
  return new MEDCouplingMultiFields(fs);
}

MEDCouplingMultiFields *MEDCouplingMultiFields::New()
{
  return new MEDCouplingMultiFields;
}

MEDCouplingMultiFields *MEDCouplingMultiFields::deepCopy() const
{
  return new MEDCouplingMultiFields(*this);
}

bool MEDCouplingMultiFields::isEqual(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const
{
  std::size_t sz=_fs.size();
  if(sz!=other->_fs.size())
    return false;
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDCouplingFieldDouble *f1=_fs[i];
      const MEDCouplingFieldDouble *f2=other->_fs[i];
      if(f1!=f2)
        {
          if(f1==0 || f2==0)
            return false;
          if(!_fs[i]->isEqual(other->_fs[i],meshPrec,valsPrec))
            return false;
        }
    }
  std::vector<int> refs1,refs2;
  std::vector<MEDCouplingMesh *> ms1=getDifferentMeshes(refs1);
  std::vector<MEDCouplingMesh *> ms2=other->getDifferentMeshes(refs2);
  if(ms1.size()!=ms2.size())
    return false;
  if(refs1!=refs2)
    return false;
  std::vector< std::vector<int> > refs3,refs4;
  std::vector<DataArrayDouble *> das1=getDifferentArrays(refs3);
  std::vector<DataArrayDouble *> das2=getDifferentArrays(refs4);
  if(das1.size()!=das2.size())
    return false;
  if(refs3!=refs4)
    return false;
  return true;
}

std::string MEDCouplingMultiFields::getName() const
{
  std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();
  for(;it!=_fs.end();it++)
    if((const MEDCouplingFieldDouble *)(*it))
      return (*it)->getName();
  return std::string();
}

std::string MEDCouplingMultiFields::getDescription() const
{
  std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();
  for(;it!=_fs.end();it++)
    if((const MEDCouplingFieldDouble *)(*it))
      return (*it)->getDescription();
  return std::string();
}

std::string MEDCouplingMultiFields::getTimeUnit() const
{
  std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();
  for(;it!=_fs.end();it++)
    if((const MEDCouplingFieldDouble *)(*it))
      return (*it)->getTimeUnit();
  return std::string();
}

double MEDCouplingMultiFields::getTimeResolution() const
{
  std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();
  for(;it!=_fs.end();it++)
    if((const MEDCouplingFieldDouble *)(*it))
      return (*it)->getTimeTolerance();
  throw INTERP_KERNEL::Exception("MEDCouplingMultiFields::getTimeResolution : no not null field !");
}

std::string MEDCouplingMultiFields::simpleRepr() const
{
  std::ostringstream ret;
  ret << "MEDCouplingMultiFields with name : \"" << getName() << "\"\n";
  ret << "Description of MEDCouplingMultiFields is : \"" << getDescription() << "\"\n";
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
  return ret.str();
}

std::string MEDCouplingMultiFields::advancedRepr() const
{
  return simpleRepr();
}

bool MEDCouplingMultiFields::isEqualWithoutConsideringStr(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const
{
  std::size_t sz=_fs.size();
  if(sz!=other->_fs.size())
    return false;
  for(std::size_t i=0;i<sz;i++)
    if(!_fs[i]->isEqualWithoutConsideringStr(other->_fs[i],meshPrec,valsPrec))
      return false;
  return true;
}

const MEDCouplingFieldDouble *MEDCouplingMultiFields::getFieldWithId(int id) const
{
  if(id>=(int)_fs.size() || id < 0)
    throw INTERP_KERNEL::Exception("MEDCouplingMultiFields::getFieldWithId : invalid id outside boundaries !");
  return _fs[id];
}

std::vector<const MEDCouplingFieldDouble *> MEDCouplingMultiFields::getFields() const
{
  std::vector<const MEDCouplingFieldDouble *> ret(_fs.size());
  std::copy(_fs.begin(),_fs.end(),ret.begin());
  return ret;
}

int MEDCouplingMultiFields::getNumberOfFields() const
{
  return (int)_fs.size();
}

const MEDCouplingFieldDouble *MEDCouplingMultiFields::getFieldAtPos(int id) const
{
  if(id<0 || id>=(int)_fs.size())
    {
      std::ostringstream oss; oss << "MEDCouplingMultiFields::getFieldAtPos : Invalid given pos : should be >=0 and < " << _fs.size() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return _fs[id];
}

void MEDCouplingMultiFields::updateTime() const
{
  std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();
  for(;it!=_fs.end();it++)
    if((const MEDCouplingFieldDouble *)(*it))
      (*it)->updateTime();
  it=_fs.begin();
  for(;it!=_fs.end();it++)
    if((const MEDCouplingFieldDouble *)(*it))
      updateTimeWith(*(*it));
}

std::size_t MEDCouplingMultiFields::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

std::vector<const BigMemoryObject *> MEDCouplingMultiFields::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  for(std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();it!=_fs.end();it++)
    ret.push_back((const MEDCouplingFieldDouble *)*it);
  return ret;
}

std::vector<MEDCouplingMesh *> MEDCouplingMultiFields::getMeshes() const
{
  std::vector<MEDCouplingMesh *> ms;
  for(std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();it!=_fs.end();it++)
    {
      const MEDCouplingMesh *m=0;
      if((const MEDCouplingFieldDouble *)(*it))
        m=(*it)->getMesh();
      ms.push_back(const_cast<MEDCouplingMesh *>(m));
    }
  return ms;
}

std::vector<MEDCouplingMesh *> MEDCouplingMultiFields::getDifferentMeshes(std::vector<int>& refs) const
{
  refs.resize(_fs.size());
  std::vector<MEDCouplingMesh *> ms;
  int id=0;
  for(std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();it!=_fs.end();it++,id++)
    {
      const MEDCouplingMesh *m=0;
      if((const MEDCouplingFieldDouble *)(*it))
        m=(*it)->getMesh();
      if(m)
        {
          std::vector<MEDCouplingMesh *>::iterator it2=std::find(ms.begin(),ms.end(),m);
          if(it2==ms.end())
            {
              ms.push_back(const_cast<MEDCouplingMesh *>(m));
              refs[id]=(int)ms.size()-1;
            }
          else
            refs[id]=(int)std::distance(ms.begin(),it2);
        }
      else
        refs[id]=-1;
    }
  return ms;
}

std::vector<DataArrayDouble *> MEDCouplingMultiFields::getArrays() const
{
  std::vector<DataArrayDouble *> tmp;
  for(std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();it!=_fs.end();it++)
    {
      std::vector<DataArrayDouble *> tmp2=(*it)->getArrays();
      tmp.insert(tmp.end(),tmp2.begin(),tmp2.end());
    }
  return tmp;
}

std::vector<DataArrayDouble *> MEDCouplingMultiFields::getDifferentArrays(std::vector< std::vector<int> >& refs) const
{
  refs.resize(_fs.size());
  int id=0;
  std::vector<DataArrayDouble *> ret;
  for(std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();it!=_fs.end();it++,id++)
    {
      std::vector<DataArrayDouble *> tmp2;
      if((const MEDCouplingFieldDouble *)(*it))
        {
          tmp2=(*it)->getArrays();
          refs[id].resize(tmp2.size());
        }
      else
        refs[id].clear();
      int id2=0;
      for(std::vector<DataArrayDouble *>::const_iterator it2=tmp2.begin();it2!=tmp2.end();it2++,id2++)
        {
          if(*it2)
            {
              std::vector<DataArrayDouble *>::iterator it3=std::find(ret.begin(),ret.end(),*it2);
              if(it3==ret.end())
                {
                  ret.push_back(*it2);
                  refs[id][id2]=(int)ret.size()-1;
                }
              else
                refs[id][id2]=(int)std::distance(ret.begin(),it3);
            }
          else
            refs[id][id2]=-1;
        }
    }
  return ret;
}

void MEDCouplingMultiFields::checkConsistencyLight() const
{
  std::vector< MCAuto<MEDCouplingFieldDouble> >::const_iterator it=_fs.begin();
  for(;it!=_fs.end();it++)
    {
      if((const MEDCouplingFieldDouble *)(*it)==0)
        throw INTERP_KERNEL::Exception("MEDCouplingMultiFields::checkConsistencyLight : There is an empty Field in array...");
      (*it)->checkConsistencyLight();
    }
}

MEDCouplingMultiFields::MEDCouplingMultiFields(const std::vector<MEDCouplingFieldDouble *>& fs):_fs(fs.size())
{
  int id=0;
  for(std::vector< MEDCouplingFieldDouble * >::const_iterator it=fs.begin();it!=fs.end();it++,id++)
    {
      if(*it)
        (*it)->incrRef();
      else
        throw INTERP_KERNEL::Exception("MEDCouplingMultiFields constructor : empty field found in vector !");
      (*it)->checkConsistencyLight();
      _fs[id]=*it;
    }
}


/*!
 * Performs deepCopy.
 */
MEDCouplingMultiFields::MEDCouplingMultiFields(const MEDCouplingMultiFields& other):RefCountObject(other)
{
  std::size_t sz=other._fs.size();
  _fs.resize(sz);
  std::vector<int> refs;
  std::vector< std::vector<int> > refs2;
  std::vector<MEDCouplingMesh *> ms=other.getDifferentMeshes(refs);
  std::size_t msLgh=ms.size();
  std::vector< MCAuto<MEDCouplingMesh> > ms2(msLgh);
  for(std::size_t i=0;i<msLgh;i++)
    ms2[i]=ms[i]->deepCopy();
  std::vector<DataArrayDouble *> das=other.getDifferentArrays(refs2);
  std::size_t dasLgth=das.size();
  std::vector< MCAuto<DataArrayDouble> > das2(dasLgth);
  for(std::size_t i=0;i<dasLgth;i++)
    das2[i]=das[i]->deepCopy();
  for(std::size_t i=0;i<sz;i++)
    {
      if((const MEDCouplingFieldDouble *)other._fs[i])
        {
          MEDCouplingFieldTemplate *tmp=MEDCouplingFieldTemplate::New(*other._fs[i]);
          _fs[i]=MEDCouplingFieldDouble::New(*tmp,other._fs[i]->getTimeDiscretization());
          tmp->decrRef();
          if(refs[i]!=-1)
            _fs[i]->setMesh(ms2[refs[i]]);
          std::size_t nbOfArr=refs2[i].size();
          std::vector<DataArrayDouble *> tmp2(nbOfArr);
          for(std::size_t j=0;j<nbOfArr;j++)
            {
              if(refs2[i][j]!=-1)
                tmp2[j]=das2[refs2[i][j]];
              else
                tmp2[j]=0;
            }
          _fs[i]->setArrays(tmp2);
          std::vector<int> tinyInfo;
          std::vector<double> tinyInfo2;
          other._fs[i]->getTimeDiscretizationUnderGround()->getTinySerializationIntInformation2(tinyInfo);
          other._fs[i]->getTimeDiscretizationUnderGround()->getTinySerializationDbleInformation2(tinyInfo2);
          _fs[i]->getTimeDiscretizationUnderGround()->finishUnserialization2(tinyInfo,tinyInfo2);
        }
    }
}

MEDCouplingMultiFields::MEDCouplingMultiFields()
{
}

void MEDCouplingMultiFields::getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<double>& tinyInfo2, int& nbOfDiffMeshes, int& nbOfDiffArr) const
{
  std::vector<int> refs;
  std::vector<MEDCouplingMesh *> ms=getDifferentMeshes(refs);
  nbOfDiffMeshes=(int)ms.size();
  std::vector< std::vector<int> > refs2;
  std::vector<DataArrayDouble *> fs=getDifferentArrays(refs2);
  nbOfDiffArr=(int)fs.size();
  //
  std::size_t sz=refs.size();//==_fs.size()
  int sz2=0;
  for(std::size_t i=0;i<sz;i++)
    sz2+=(int)refs2[i].size();
  //
  tinyInfo2.clear();
  std::vector<int> doubleDaInd(sz);
  std::vector<int> timeDiscrInt;
  tinyInfo.resize(sz2+5*sz+3);
  tinyInfo[0]=(int)sz;
  tinyInfo[1]=sz2;
  for(std::size_t i=0;i<sz;i++)
    {
      std::vector<double> tmp;
      std::vector<int> tmp2;
      _fs[i]->getTimeDiscretizationUnderGround()->getTinySerializationDbleInformation2(tmp);
      _fs[i]->getTimeDiscretizationUnderGround()->getTinySerializationIntInformation2(tmp2);
      tinyInfo[3*sz+3+i]=(int)tmp.size();
      tinyInfo[4*sz+3+i]=(int)tmp2.size();
      tinyInfo2.insert(tinyInfo2.end(),tmp.begin(),tmp.end());
      timeDiscrInt.insert(timeDiscrInt.end(),tmp2.begin(),tmp2.end());
    }
  int sz3=(int)timeDiscrInt.size();
  tinyInfo[2]=sz3;
  //
  for(std::size_t i=0;i<sz;i++)
    tinyInfo[i+3]=refs[i];
  for(std::size_t i=0;i<sz;i++)
    tinyInfo[i+sz+3]=(int)refs2[i].size();
  for(std::size_t i=0;i<sz;i++)
    tinyInfo[i+2*sz+3]=(int)_fs[i]->getTimeDiscretization();
  int k=0;
  for(std::size_t i=0;i<sz;i++)
    for(std::vector<int>::const_iterator it=refs2[i].begin();it!=refs2[i].end();it++,k++)
      tinyInfo[5*sz+k+3]=*it;
  tinyInfo.insert(tinyInfo.end(),timeDiscrInt.begin(),timeDiscrInt.end());//tinyInfo has lgth==sz3+sz2+5*sz+3
}

void MEDCouplingMultiFields::finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD,
                                                   const std::vector<MEDCouplingFieldTemplate *>& ft, const std::vector<MEDCouplingMesh *>& ms,
                                                   const std::vector<DataArrayDouble *>& das)
{
  int sz=tinyInfoI[0];
  _fs.resize(sz);
  int sz2=tinyInfoI[1];
  // dealing with ft with no mesh set.
  for(int i=0;i<sz;i++)
    {
      int meshId=tinyInfoI[3+i];
      if(meshId!=-1)
        ft[i]->setMesh(ms[meshId]);
    }
  // dealing with fieldtemplate->fielddouble
  int k=0;
  int offI=0;
  int offD=0;
  for(int i=0;i<sz;i++)
    {
      _fs[i]=MEDCouplingFieldDouble::New(*ft[i],(TypeOfTimeDiscretization)tinyInfoI[2*sz+3+i]);
      int sz3=tinyInfoI[sz+i+3];
      std::vector<DataArrayDouble *> tmp(sz3);
      for(int j=0;j<sz3;j++,k++)
        {
          int daId=tinyInfoI[5*sz+k+3];
          if(daId!=-1)
            tmp[j]=das[daId];
          else
            tmp[j]=0;
        }
      _fs[i]->setArrays(tmp);
      // time discr tiny info
      int lgthI=tinyInfoI[4*sz+3+i];
      int lgthD=tinyInfoI[3*sz+3+i];
      //
      std::vector<int> tdInfoI(tinyInfoI.begin()+sz2+5*sz+3+offI,tinyInfoI.begin()+sz2+5*sz+3+offI+lgthI);
      std::vector<double> tdInfoD(tinyInfoD.begin()+offD,tinyInfoD.begin()+offD+lgthD);
      _fs[i]->getTimeDiscretizationUnderGround()->finishUnserialization2(tdInfoI,tdInfoD);
      //
      offI+=lgthI;
      offD+=lgthD;
    }
}
