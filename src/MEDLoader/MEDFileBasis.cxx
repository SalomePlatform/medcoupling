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

#include "MEDFileBasis.hxx"

#include <cstring>

using namespace MEDCoupling;

MEDFileString::MEDFileString(int maxLgth):_max_lgth(maxLgth),_content(new char[maxLgth+1])
{
  std::fill(_content,_content+maxLgth+1,'\0');
}

MEDFileString::~MEDFileString()
{
  delete [] _content;
}

void MEDFileString::clear()
{
  std::fill(_content,_content+_max_lgth+1,'\0');
}

void MEDFileString::set(const char *s)
{
  if((int)strlen(s)>_max_lgth)
    throw INTERP_KERNEL::Exception("Name is too long to be stored in MEDfile !");
  clear();
  strcpy(_content,s);
}

std::string MEDFileString::getRepr() const
{
  return std::string(_content);
}

