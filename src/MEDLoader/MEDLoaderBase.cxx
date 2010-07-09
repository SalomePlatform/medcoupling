//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#include "MEDLoaderBase.hxx"
#include "InterpKernelException.hxx"

#include <fstream>

int MEDLoaderBase::getStatusOfFile(const char *fileName)
{
  std::ifstream ifs;
  ifs.open(fileName);
  unsigned int res=0;
  if((ifs.rdstate() & std::ifstream::failbit)!=0)
    {
      res+=1;
      ifs.close();
    }
  std::ofstream ofs(fileName,std::ios_base::app);
  if((ofs.rdstate() & std::ofstream::failbit)!=0)
    {
      ofs.close();
      res+=2;
    }
  switch(res)
    {
    case 0:
      return EXIST_RW;
    case 1:
      {
        std::ifstream ifs2;
        ifs2.open(fileName);
        if((ifs2.rdstate() & std::ifstream::failbit)!=0)
          return EXIST_WRONLY;
        else
          return NOT_EXIST;
      }
    case 2:
      return EXIST_RDONLY;
    case 3:
      return DIR_LOCKED;
    default:
      throw INTERP_KERNEL::Exception("Internal error !");
    }
}

char *MEDLoaderBase::buildEmptyString(int lgth)
{
  char *ret=new char[lgth+1];
  std::fill(ret,ret+lgth,' ');
  ret[lgth]='\0';
  return ret;
}
