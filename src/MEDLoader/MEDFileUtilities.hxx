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

#ifndef __MEDFILEUTILITIES_HXX__
#define __MEDFILEUTILITIES_HXX__

#include "InterpKernelException.hxx"
#include "MEDLoaderDefines.hxx"

#include "MCAuto.hxx"
#include "MEDCouplingMemArray.hxx"

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
    operator med_idt() const { return _fid; }
    ~AutoFid();
  private:
    med_idt _fid;
  };
}

namespace MEDCoupling
{
  class MEDLOADER_EXPORT MEDFileWritable
  {
  public:
    MEDFileWritable();
    virtual ~MEDFileWritable() {}
    void copyOptionsFrom(const MEDFileWritable& other) const;
    int getTooLongStrPolicy() const;
    void setTooLongStrPolicy(int newVal);
    int getZipConnPolicy();
    void setZipConnPolicy(int newVal);
    static std::string FileNameFromFID(med_idt fid);
  protected://policies on write
    mutable int _too_long_str;
    mutable int _zipconn_pol;
  };

  class MEDFileWritableStandAlone : public MEDFileWritable
  {
  public:
    MEDLOADER_EXPORT virtual void writeLL(med_idt fid) const = 0;
    MEDLOADER_EXPORT virtual void write(const std::string& fileName, int mode) const;
    MEDLOADER_EXPORT virtual void write30(const std::string& fileName, int mode) const;
    MEDLOADER_EXPORT MCAuto<DataArrayByte> serialize() const;
    MEDLOADER_EXPORT static std::string GenerateUniqueDftFileNameInMem();
  public:
    MEDLOADER_EXPORT static const char DFT_FILENAME_IN_MEM[];
    template<class T>
    static T *BuildFromMemoryChunk(DataArrayByte *db);
  };
  
  MEDFileUtilities::AutoFid OpenMEDFileForRead(const std::string& fileName);
}

#endif
