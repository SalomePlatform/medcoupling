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

#ifndef __MEDFILEUTILITIES_HXX__
#define __MEDFILEUTILITIES_HXX__

#include "InterpKernelException.hxx"
#include "MEDLoaderDefines.hxx"

#include "med.h"

namespace MEDFileUtilities
{
  med_access_mode TraduceWriteMode(int medloaderwritemode);
  const char *GetReadableMEDFieldType(med_field_type ft);
  void CheckMEDCode(int code, med_idt fid, const std::string& msg);
  void CheckFileForRead(const std::string& fileName);

  class AutoFid
  {
  public:
    AutoFid(med_idt fid);
    operator med_idt() const;
    ~AutoFid();
  private:
    med_idt _fid;
  };
}
  
namespace ParaMEDMEM
{
  class MEDLOADER_EXPORT MEDFileWritable
  {
  public:
    MEDFileWritable();
    void copyOptionsFrom(const MEDFileWritable& other) const;
    int getTooLongStrPolicy() const;
    void setTooLongStrPolicy(int newVal);
    int getZipConnPolicy();
    void setZipConnPolicy(int newVal);
  protected://policies on write
    mutable int _too_long_str;
    mutable int _zipconn_pol;
  };
}

#endif
