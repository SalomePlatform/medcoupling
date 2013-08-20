// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

RefCountObject::RefCountObject():_cnt(1)
{
}

RefCountObject::RefCountObject(const RefCountObject& other):_cnt(1)
{
}

/*!
 * Do nothing here ! It is not a bug ( I hope :) ) because all subclasses that
 * copies using operator= should not copy the ref counter of \a other !
 */
RefCountObject& RefCountObject::operator=(const RefCountObject& other)
{
  return *this;
}

bool RefCountObject::decrRef() const
{
  bool ret=((--_cnt)==0);
  if(ret)
    delete this;
  return ret;
}

void RefCountObject::incrRef() const
{
  _cnt++;
}

int RefCountObject::getRCValue() const
{
  return _cnt;
}

RefCountObject::~RefCountObject()
{
}
