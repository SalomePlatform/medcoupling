// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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
#include "MED_version.h"

#include <sstream>

using namespace ParaMEDMEM;

const char *ParaMEDMEM::MEDCouplingVersionStr()
{
  return SALOMEMED_VERSION_STR;
}

int ParaMEDMEM::MEDCouplingVersion()
{
  return SALOMEMED_VERSION;
}

void ParaMEDMEM::MEDCouplingVersionMajMinRel(int& maj, int& minor, int& releas)
{
  int ver=SALOMEMED_VERSION;
  maj=(ver & 0xFF0000) >> 16;
  minor=(ver & 0xFF00) >> 8;
  releas=(ver & 0xFF);
}

int ParaMEDMEM::MEDCouplingSizeOfVoidStar()
{
  return 8*sizeof(std::size_t);
}

/*!
 * If true is returned it is a LittleEndian machine.
 * If false it is a BigEndian machine.
 * \return the coding mode of integers of the machine.
 */
bool ParaMEDMEM::MEDCouplingByteOrder()
{
  unsigned int x(1);
  unsigned char *xc(reinterpret_cast<unsigned char *>(&x));
  return xc[0]==1;
}

const char *ParaMEDMEM::MEDCouplingByteOrderStr()
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
