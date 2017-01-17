// Copyright (C) 2007-2017  CEA/DEN, EDF R&D
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

#ifndef __MEDFILEMESHSUPPORT_HXX__
#define __MEDFILEMESHSUPPORT_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDFileUtilities.txx"

#include "MEDCouplingRefCountObject.hxx"

namespace MEDCoupling
{
  class MEDFileMeshSupport : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT static MEDFileMeshSupport *New(const std::string& fileName);
    MEDLOADER_EXPORT static MEDFileMeshSupport *New(const std::string& fileName, int smid);
    MEDLOADER_EXPORT static MEDFileMeshSupport *New(med_idt fid, int smid);
    MEDLOADER_EXPORT static MEDFileMeshSupport *New();
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT int getSpaceDim() const;
    MEDLOADER_EXPORT void setSpaceDim(int dim);
    MEDLOADER_EXPORT int getMeshDim() const;
    MEDLOADER_EXPORT void setMeshDim(int dim);
  private:
    void writeLL(med_idt fid) const;
    MEDFileMeshSupport(med_idt fid, int smid);
    MEDFileMeshSupport();
  private:
    int _space_dim;
    int _mesh_dim;
    std::string _description;
  };
}

#endif
