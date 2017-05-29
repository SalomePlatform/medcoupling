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

#include "MEDCouplingRefCountObject.hxx"
#include "MEDCoupling_version.h"

#include "InterpKernelException.hxx"

#include <sstream>
#include <algorithm>

using namespace MEDCoupling;

GlobalDict *GlobalDict::UNIQUE_INSTANCE=0;

const char *MEDCoupling::MEDCouplingVersionStr()
{
  return MEDCOUPLING_VERSION_STR;
}

int MEDCoupling::MEDCouplingVersion()
{
  return MEDCOUPLING_VERSION;
}

void MEDCoupling::MEDCouplingVersionMajMinRel(int& maj, int& minor, int& releas)
{
  int ver=MEDCOUPLING_VERSION;
  maj=(ver & 0xFF0000) >> 16;
  minor=(ver & 0xFF00) >> 8;
  releas=(ver & 0xFF);
}

int MEDCoupling::MEDCouplingSizeOfVoidStar()
{
  return 8*sizeof(std::size_t);
}

/*!
 * If true is returned it is a LittleEndian machine.
 * If false it is a BigEndian machine.
 * \return the coding mode of integers of the machine.
 */
bool MEDCoupling::MEDCouplingByteOrder()
{
  unsigned int x(1);
  unsigned char *xc(reinterpret_cast<unsigned char *>(&x));
  return xc[0]==1;
}

const char *MEDCoupling::MEDCouplingByteOrderStr()
{
  static const char LITTLEENDIAN_STR[]="LittleEndian";
  static const char BIGENDIAN_STR[]="BigEndian";
  if(MEDCouplingByteOrder())
    return LITTLEENDIAN_STR;
  else
    return BIGENDIAN_STR;
}

//=

std::size_t BigMemoryObject::getHeapMemorySize() const
{
  std::size_t ret(getHeapMemorySizeWithoutChildren());
  std::vector<const BigMemoryObject *> v(getDirectChildren());
  std::set<const BigMemoryObject *> s1,s2(v.begin(),v.end());
  return ret+GetHeapMemoryOfSet(s1,s2);
}

/*!
 * This method returns all the progeny of \a this (this is \b not included in returned vector).
 * All the progeny means all the subobjects (children), subsubobjects (little children), ... of \a this.
 * The elements in returned array are reported only once even if they appear several times in the progeny of \a this.
 */
std::vector<const BigMemoryObject *> BigMemoryObject::getAllTheProgeny() const
{
  std::vector<const BigMemoryObject *> s1(getDirectChildren());
  std::vector<const BigMemoryObject *> ret;
  while(!s1.empty())
    {
      ret.insert(ret.end(),s1.begin(),s1.end());
      std::vector<const BigMemoryObject *> s3;
      for(std::vector<const BigMemoryObject *>::const_iterator it0=s1.begin();it0!=s1.end();it0++)
        {
          std::vector<const BigMemoryObject *> s2;
          if(*it0)
            s2=(*it0)->getDirectChildren();
          for(std::vector<const BigMemoryObject *>::const_iterator it1=s2.begin();it1!=s2.end();it1++)
            {
              if(*it1)
                if(std::find(ret.begin(),ret.end(),*it1)==ret.end())
                  s3.push_back(*it1);
            }
        }
      s1=s3;
    }
  return ret;
}

/*!
 * This method scan all the progeny of \a this (\a this excluded) to see if \a obj is part of it.
 * If obj is NULL false is returned.
 * \sa BigMemoryObject::getAllTheProgeny
 */
bool BigMemoryObject::isObjectInTheProgeny(const BigMemoryObject *obj) const
{
  if(!obj)
    return false;
  std::vector<const BigMemoryObject *> objs(getAllTheProgeny());
  return std::find(objs.begin(),objs.end(),obj)!=objs.end();
}

std::size_t BigMemoryObject::GetHeapMemorySizeOfObjs(const std::vector<const BigMemoryObject *>& objs)
{
  std::size_t ret(0);
  std::set<const BigMemoryObject *> s1,s2;
  for(std::vector<const BigMemoryObject *>::const_iterator it0=objs.begin();it0!=objs.end();it0++)
    {
      if(*it0)
        if(s1.find(*it0)==s1.end())
          {
            std::vector<const BigMemoryObject *> vTmp((*it0)->getDirectChildren());
            s2.insert(vTmp.begin(),vTmp.end());
            ret+=(*it0)->getHeapMemorySizeWithoutChildren();
            s1.insert(*it0);
          }
    }
  return ret+GetHeapMemoryOfSet(s1,s2);
}

std::size_t BigMemoryObject::GetHeapMemoryOfSet(std::set<const BigMemoryObject *>& s1, std::set<const BigMemoryObject *>& s2)
{
  std::size_t ret(0);
  while(!s2.empty())
    {
      std::set<const BigMemoryObject *> s3;
      for(std::set<const BigMemoryObject *>::const_iterator it=s2.begin();it!=s2.end();it++)
        {
          if(s1.find(*it)==s1.end())
            {
              ret+=(*it)->getHeapMemorySizeWithoutChildren();
              s1.insert(*it);
              std::vector<const BigMemoryObject *> v2((*it)->getDirectChildren());
              for(std::vector<const BigMemoryObject *>::const_iterator it2=v2.begin();it2!=v2.end();it2++)
                if(s1.find(*it2)==s1.end())
                  s3.insert(*it2);
            }
        }
      s2=s3;
    }
  return ret;
}

std::string BigMemoryObject::getHeapMemorySizeStr() const
{
  static const char *UNITS[4]={"B","kB","MB","GB"};
  std::size_t m(getHeapMemorySize());
  std::ostringstream oss; oss.precision(3);
  std::size_t remain(0);
  int i(0);
  for(;i<4;i++)
    {
      if(m<1024)
        {
          oss << m;
          if(remain!=0)
            {
              std::ostringstream oss2; oss2 << std::fixed << ((double)remain)/1024.;
              std::string s(oss2.str());
              s=s.substr(1,4);
              std::size_t pos(s.find_last_not_of('0'));
              if(pos==4)
                oss << s;
              else
                oss << s.substr(0,pos+1);
            }
          oss << " " << UNITS[i];
          break;
        }
      else
        {
          if(i!=3)
            {
              remain=(m%1024);
              m/=1024;
            }
        }
    }
  if(i==4)
    oss << m << " " << UNITS[3];
  return oss.str();
}

std::vector<const BigMemoryObject *> BigMemoryObject::getDirectChildren() const
{
  std::vector<const BigMemoryObject *> ret;
  std::vector<const BigMemoryObject *> retWithNull(getDirectChildrenWithNull());
  for(std::vector<const BigMemoryObject *>::const_iterator it=retWithNull.begin();it!=retWithNull.end();it++)
    if(*it)
      ret.push_back(*it);
  return ret;
}

BigMemoryObject::~BigMemoryObject()
{
}

//=

RefCountObjectOnly::RefCountObjectOnly():_cnt(1)
{
}

RefCountObjectOnly::RefCountObjectOnly(const RefCountObjectOnly& other):_cnt(1)
{
}

bool RefCountObjectOnly::decrRef() const
{
  bool ret=((--_cnt)==0);
  if(ret)
    delete this;
  return ret;
}

void RefCountObjectOnly::incrRef() const
{
  _cnt++;
}

int RefCountObjectOnly::getRCValue() const
{
  return _cnt;
}

RefCountObjectOnly::~RefCountObjectOnly()
{
}

/*!
 * Do nothing here ! It is not a bug ( I hope :) ) because all subclasses that
 * copies using operator= should not copy the ref counter of \a other !
 */
RefCountObjectOnly& RefCountObjectOnly::operator=(const RefCountObjectOnly& other)
{
  return *this;
}

//=

RefCountObject::RefCountObject()
{
}

RefCountObject::RefCountObject(const RefCountObject& other):RefCountObjectOnly(other)
{
}

RefCountObject::~RefCountObject()
{
}

//=

GlobalDict *GlobalDict::GetInstance()
{
  if(!UNIQUE_INSTANCE)
    UNIQUE_INSTANCE=new GlobalDict;
  return UNIQUE_INSTANCE;
}

bool GlobalDict::hasKey(const std::string& key) const
{
  std::map<std::string, std::string>::const_iterator it(_my_map.find(key));
  return it!=_my_map.end();
}

std::string GlobalDict::value(const std::string& key) const
{
  std::map<std::string, std::string>::const_iterator it(_my_map.find(key));
  if(it==_my_map.end())
    {
      std::ostringstream oss;
      oss << "GlobalDict::value : key \"" << key << "\" is not in map !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return (*it).second;
}

std::vector<std::string> GlobalDict::keys() const
{
  std::vector<std::string> ret;
  for(std::map<std::string, std::string>::const_iterator it=_my_map.begin();it!=_my_map.end();it++)
    ret.push_back((*it).first);
  return ret;
}

void GlobalDict::erase(const std::string& key)
{
  std::map<std::string, std::string>::iterator it(_my_map.find(key));
  if(it==_my_map.end())
    {
      std::ostringstream oss;
      oss << "GlobalDict::erase : key \"" << key << "\" is not in map !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _my_map.erase(it);
}

void GlobalDict::clear()
{
  _my_map.clear();
}

void GlobalDict::setKeyValue(const std::string& key, const std::string& val)
{
  std::map<std::string, std::string>::const_iterator it(_my_map.find(key));
  if(it!=_my_map.end())
    {
      std::ostringstream oss;
      oss << "GlobalDict::setKeyValue : key \"" << key << "\" already exists !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _my_map[key]=val;
}

void GlobalDict::setKeyValueForce(const std::string& key, const std::string& val)
{
  _my_map[key]=val;
}

std::string GlobalDict::printSelf() const
{
  std::ostringstream oss;
  for(std::map<std::string, std::string>::const_iterator it=_my_map.begin();it!=_my_map.end();it++)
    {
      oss << "(" << (*it).first << "," << (*it).second << ")" << std::endl;
    }
  return oss.str();
}
