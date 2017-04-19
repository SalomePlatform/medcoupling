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

#include "MEDFileMeshReadSelector.hxx"

using namespace MEDCoupling;

MEDFileMeshReadSelector::MEDFileMeshReadSelector():_code(0xFFFFFFFF)
{
}

MEDFileMeshReadSelector::MEDFileMeshReadSelector(unsigned int code):_code(code)
{
}

unsigned int MEDFileMeshReadSelector::getCode() const
{
  return _code;
}

void MEDFileMeshReadSelector::setCode(unsigned int newCode)
{
  _code=newCode;
}

bool MEDFileMeshReadSelector::isCellFamilyFieldReading() const
{
  return _code & 0x00000001;
}

bool MEDFileMeshReadSelector::isNodeFamilyFieldReading() const
{
  return _code & 0x00000002;
}

bool MEDFileMeshReadSelector::isCellNameFieldReading() const
{
  return _code & 0x00000004;
}

bool MEDFileMeshReadSelector::isNodeNameFieldReading() const
{
  return _code & 0x00000008;
}

bool MEDFileMeshReadSelector::isCellNumFieldReading() const
{
  return _code & 0x00000010;
}

bool MEDFileMeshReadSelector::isNodeNumFieldReading() const
{
  return _code & 0x00000020;
}

bool MEDFileMeshReadSelector::isGlobalNodeNumFieldReading() const
{
  return _code & 0x00000040;
}

void MEDFileMeshReadSelector::setCellFamilyFieldReading(bool b)
{
  unsigned int code(_code & 0xFFFFFFFE);
  unsigned int b2=b?1:0;
  //b2<<=0;
  code+=b2;
  _code=code;
}

void MEDFileMeshReadSelector::setNodeFamilyFieldReading(bool b)
{
  unsigned int code(_code & 0xFFFFFFFD);
  unsigned int b2=b?1:0;
  b2<<=1;
  code+=b2;
  _code=code;
}

void MEDFileMeshReadSelector::setCellNameFieldReading(bool b)
{
  unsigned int code(_code & 0xFFFFFFFB);
  unsigned int b2=b?1:0;
  b2<<=2;
  code+=b2;
  _code=code;
}

void MEDFileMeshReadSelector::setNodeNameFieldReading(bool b)
{
  unsigned int code(_code & 0xFFFFFFF7);
  unsigned int b2=b?1:0;
  b2<<=3;
  code+=b2;
  _code=code;
}

void MEDFileMeshReadSelector::setCellNumFieldReading(bool b)
{
  unsigned int code(_code & 0xFFFFFFEF);
  unsigned int b2=b?1:0;
  b2<<=4;
  code+=b2;
  _code=code;
}

void MEDFileMeshReadSelector::setNodeNumFieldReading(bool b)
{
  unsigned int code(_code & 0xFFFFFFDF);
  unsigned int b2=b?1:0;
  b2<<=5;
  code+=b2;
  _code=code;
}

void MEDFileMeshReadSelector::setGlobalNodeNumFieldReading(bool b)
{
  unsigned int code(_code & 0xFFFFFFBF);
  unsigned int b2=b?1:0;
  b2<<=6;
  code+=b2;
  _code=code;
}

void MEDFileMeshReadSelector::reprAll(std::ostream& str) const
{
  str << "MEDFileMeshReadSelector (code=" << _code << ") : \n";
  str << "Read family field on cells : " << ReprStatus(isCellFamilyFieldReading()) << std::endl;
  str << "Read family field on nodes : " << ReprStatus(isNodeFamilyFieldReading()) << std::endl;
  str << "Read name field on cells : " << ReprStatus(isCellNameFieldReading()) << std::endl;
  str << "Read name field on nodes : " << ReprStatus(isNodeNameFieldReading()) << std::endl;
  str << "Read number field on cells : " << ReprStatus(isCellNumFieldReading()) << std::endl;
  str << "Read number field name on nodes : " << ReprStatus(isNodeNumFieldReading()) << std::endl;
  str << "Read global number field name on nodes : " << ReprStatus(isGlobalNodeNumFieldReading());
}

std::string MEDFileMeshReadSelector::ReprStatus(bool v)
{
  if(v)
    return std::string("ON");
  else
    return std::string("OFF");
}

