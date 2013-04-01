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

#ifndef __MEDFILEUTILITIES_HXX__
#define __MEDFILEUTILITIES_HXX__

#include "InterpKernelException.hxx"
#include "MEDLoaderDefines.hxx"

#include "med.h"

namespace MEDFileUtilities
{
  med_access_mode TraduceWriteMode(int medloaderwritemode) throw(INTERP_KERNEL::Exception);
  int TraduceFieldType(med_field_type ft) throw(INTERP_KERNEL::Exception);
  void CheckMEDCode(int code, med_idt fid, const char *msg) throw(INTERP_KERNEL::Exception);
  void CheckFileForRead(const char *fileName) throw(INTERP_KERNEL::Exception);

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
    int getTooLongStrPolicy() const throw(INTERP_KERNEL::Exception);
    void setTooLongStrPolicy(int newVal) throw(INTERP_KERNEL::Exception);
    int getZipConnPolicy() throw(INTERP_KERNEL::Exception);
    void setZipConnPolicy(int newVal) throw(INTERP_KERNEL::Exception);
  protected://policies on write
    mutable int _too_long_str;
    mutable int _zipconn_pol;
  };
}

#endif
