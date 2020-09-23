// Copyright (C) 2016-2020  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

#pragma once

#include "MEDLoaderDefines.hxx"
#include "MCType.hxx"

namespace MEDCoupling
{
  template<class T>
  struct MEDLOADER_EXPORT MLFieldTraits
  {
    typedef T EltType;
  };

  class MEDFileFieldMultiTS;
  class MEDFileField1TS;
  class MEDFileInt32FieldMultiTS;
  class MEDFileInt64FieldMultiTS;
  class MEDFileInt32Field1TS;
  class MEDFileInt64Field1TS;
  class MEDFileFloatFieldMultiTS;
  class MEDFileFloatField1TS;
  class MEDFileField1TSWithoutSDA;
  class MEDFileInt32Field1TSWithoutSDA;
  class MEDFileInt64Field1TSWithoutSDA;
  class MEDFileFloatField1TSWithoutSDA;
  class MEDFileFieldMultiTSWithoutSDA;
  class MEDFileInt32FieldMultiTSWithoutSDA;
  class MEDFileInt64FieldMultiTSWithoutSDA;
  class MEDFileFloatFieldMultiTSWithoutSDA;
  
  template<>
  struct MEDLOADER_EXPORT MLFieldTraits<double>
  {
    typedef MEDFileFieldMultiTSWithoutSDA FMTSWSDAType;
    typedef MEDFileFieldMultiTS FMTSType;
    typedef MEDFileField1TS F1TSType;
    typedef MEDFileField1TSWithoutSDA F1TSWSDAType;
  };

  template<>
  struct MEDLOADER_EXPORT MLFieldTraits<float>
  {
    typedef MEDFileFloatFieldMultiTSWithoutSDA FMTSWSDAType;
    typedef MEDFileFloatFieldMultiTS FMTSType;
    typedef MEDFileFloatField1TS F1TSType;
    typedef MEDFileFloatField1TSWithoutSDA F1TSWSDAType;
  };
  
  template<>
  struct MEDLOADER_EXPORT MLFieldTraits<Int32>
  {
    typedef MEDFileInt32FieldMultiTSWithoutSDA FMTSWSDAType;
    typedef MEDFileInt32FieldMultiTS FMTSType;
    typedef MEDFileInt32Field1TS F1TSType;
    typedef MEDFileInt32Field1TSWithoutSDA F1TSWSDAType;
  };

  template<>
  struct MEDLOADER_EXPORT MLFieldTraits<Int64>
  {
    typedef MEDFileInt64FieldMultiTSWithoutSDA FMTSWSDAType;
    typedef MEDFileInt64FieldMultiTS FMTSType;
    typedef MEDFileInt64Field1TS F1TSType;
    typedef MEDFileInt64Field1TSWithoutSDA F1TSWSDAType;
  };

  template<class T>
  struct MEDLOADER_EXPORT MLMeshTraits
  {
  };
  
  class MEDFileUMesh;
  class MEDFileCMesh;
  class MEDFileCurveLinearMesh;
  
  template<>
  struct MEDLOADER_EXPORT MLMeshTraits<MEDFileUMesh>
  {
    static const char ClassName[];
  };
  
  template<>
  struct MEDLOADER_EXPORT MLMeshTraits<MEDFileCMesh>
  {
    static const char ClassName[];
  };

  template<>
  struct MEDLOADER_EXPORT MLMeshTraits<MEDFileCurveLinearMesh>
  {
    static const char ClassName[];
  };
}
