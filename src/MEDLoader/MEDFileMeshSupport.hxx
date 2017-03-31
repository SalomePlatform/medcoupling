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
#include "MEDFileMesh.hxx"

#include "MEDCouplingRefCountObject.hxx"

namespace MEDCoupling
{
  class MEDFileMeshSupports : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT static MEDFileMeshSupports *New(const std::string& fileName);
    MEDLOADER_EXPORT static MEDFileMeshSupports *New(med_idt fid);
    MEDLOADER_EXPORT static MEDFileMeshSupports *New();
  public:
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    MEDLOADER_EXPORT std::vector<std::string> getSupMeshNames() const;
    MEDLOADER_EXPORT const MEDFileUMesh *getSupMeshWithName(const std::string& name) const;
    MEDLOADER_EXPORT int getNumberOfNodesInConnOf(TypeOfField entity, const std::string& name) const;
  private:
    MEDFileMeshSupports(med_idt fid);
    MEDFileMeshSupports();
    ~MEDFileMeshSupports();
  private:
    std::vector< MCAuto<MEDFileUMesh> > _supports;
  };
}

#endif
