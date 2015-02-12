// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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

#include "MEDFileUtilities.hxx"
#include "MEDLoaderBase.hxx"

#include <sstream>

med_access_mode MEDFileUtilities::TraduceWriteMode(int medloaderwritemode)
{
  switch(medloaderwritemode)
  {
    case 2:
      return MED_ACC_CREAT;
    case 1:
      return MED_ACC_RDEXT;
    case 0:
      return MED_ACC_RDWR;
    default:
      throw INTERP_KERNEL::Exception("Invalid write mode specified ! must be 0(write with no question), 1(append) or 2(creation)");
  }
}

const char *MEDFileUtilities::GetReadableMEDFieldType(med_field_type ft)
{
  static const char medFloat64[]="MED_FLOAT64";
  static const char medInt32[]="MED_INT32";
  static const char medInt64[]="MED_INT64";
  switch(ft)
  {
    case MED_FLOAT64:
      return medFloat64;
    case MED_INT32:
      return medInt32;
    case MED_INT64:
      return medInt64;
    default:
      throw INTERP_KERNEL::Exception("Non supported field type ! Should be FLOAT64, INT32 or INT64 !");
  }
}

void MEDFileUtilities::CheckMEDCode(int code, med_idt fid, const std::string& msg)
{
  if(code<0)
    {
      std::ostringstream oss;
      oss << "MEDFile has returned an error code (" << code <<") : " << msg;
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void MEDFileUtilities::CheckFileForRead(const std::string& fileName)
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
  AutoFid fid=MEDfileOpen(fileName.c_str(),MED_ACC_RDONLY);
  if(fid<0)
    {
      oss << " has been detected as unreadable by MED file : impossible to read anything !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  oss << " has been detected readable but ";
  int major,minor,release;
  MEDfileNumVersionRd(fid,&major,&minor,&release);
  if(major<2 || (major==2 && minor<2))
    {
      oss << "version of MED file is < 2.2 : impossible to read anything !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

MEDFileUtilities::AutoFid::AutoFid(med_idt fid):_fid(fid)
{
}

MEDFileUtilities::AutoFid::~AutoFid()
{
  MEDfileClose(_fid);
}

ParaMEDMEM::MEDFileWritable::MEDFileWritable():_too_long_str(0),_zipconn_pol(2)
{
}

void ParaMEDMEM::MEDFileWritable::copyOptionsFrom(const MEDFileWritable& other) const
{
  _too_long_str=other._too_long_str;
  _zipconn_pol=other._zipconn_pol;
}

int ParaMEDMEM::MEDFileWritable::getTooLongStrPolicy() const
{
  return _too_long_str;
}

void ParaMEDMEM::MEDFileWritable::setTooLongStrPolicy(int newVal)
{
  if(newVal!=2 && newVal!=1 && newVal!=0)
    throw INTERP_KERNEL::Exception("MEDFileWritable::setTooLongStrPolicy : invalid policy should be in 0,1 or 2 !");
  _too_long_str=newVal;
}

int ParaMEDMEM::MEDFileWritable::getZipConnPolicy()
{
  return _zipconn_pol;
}

void ParaMEDMEM::MEDFileWritable::setZipConnPolicy(int newVal)
{
  _zipconn_pol=newVal;
}
