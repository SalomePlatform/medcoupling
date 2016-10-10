// Copyright (C) 2016  CEA/DEN, EDF R&D
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

#ifndef __MEDLOADERTRAITS_HXX__
#define __MEDLOADERTRAITS_HXX__

#include "MEDLoaderDefines.hxx"

namespace MEDCoupling
{
  template<class T>
  struct MEDLOADER_EXPORT MLFieldTraits
  {
    typedef T EltType;
  };

  class MEDFileFieldMultiTS;
  class MEDFileField1TS;
  class MEDFileIntFieldMultiTS;
  class MEDFileIntField1TS;
  class MEDFileField1TSWithoutSDA;
  class MEDFileIntField1TSWithoutSDA;
  
  template<>
  struct MEDLOADER_EXPORT MLFieldTraits<double>
  {
    typedef MEDFileFieldMultiTS FMTSType;
    typedef MEDFileField1TS F1TSType;
    typedef MEDFileField1TSWithoutSDA F1TSWSDAType;
  };
  
  template<>
  struct MEDLOADER_EXPORT MLFieldTraits<int>
  {
    typedef MEDFileIntFieldMultiTS FMTSType;
    typedef MEDFileIntField1TS F1TSType;
    typedef MEDFileIntField1TSWithoutSDA F1TSWSDAType;
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

#endif
