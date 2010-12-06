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

#include "MEDFileUtilities.hxx"
#include "MEDLoaderBase.hxx"

#include <sstream>

med_mode_acces MEDFileUtilities::TraduceWriteMode(int medloaderwritemode) throw(INTERP_KERNEL::Exception)
{
  switch(medloaderwritemode)
    {
    case 2:
      return MED_CREATION;
    case 1:
      return MED_LECTURE_AJOUT;
    case 0:
      return MED_LECTURE_ECRITURE;
    default:
      throw INTERP_KERNEL::Exception("Invalid write mode specified ! must be 0(write with no question), 1(append) or 2(creation)");
    }
}

void MEDFileUtilities::CheckMEDCode(int code, med_idt fid, const char *msg) throw(INTERP_KERNEL::Exception)
{
  if(code!=0)
    {
      std::ostringstream oss;
      oss << "MEDFile has returned an error code (" << code <<") : " << msg;
      MEDfermer(fid);
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void MEDFileUtilities::CheckFileForRead(const char *fileName) throw(INTERP_KERNEL::Exception)
{
  int status=MEDLoaderBase::getStatusOfFile(fileName);
  std::ostringstream oss;
  oss << " File : \"" << fileName << "\"";
  switch(status)
    {
    case MEDLoaderBase::DIR_LOCKED:
      {
        oss << " has been detected as unreadable : impossible to read anything !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    case MEDLoaderBase::NOT_EXIST:
      {
        oss << " has been detected as NOT EXISTING : impossible to read anything !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    case MEDLoaderBase::EXIST_WRONLY:
      {
        oss << " has been detected as WRITE ONLY : impossible to read anything !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    }
  int fid=MEDouvrir((char *)fileName,MED_LECTURE);
  if(fid<0)
    {
      oss << " has been detected as unreadable by MED file : impossible to read anything !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  oss << " has been detected readable but ";
  int major,minor,release;
  MEDversionLire(fid,&major,&minor,&release);
  if(major<2 || (major==2 && minor<2))
    {
      oss << "version of MED file is < 2.2 : impossible to read anything !";
      MEDfermer(fid);
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MEDfermer(fid);
}
