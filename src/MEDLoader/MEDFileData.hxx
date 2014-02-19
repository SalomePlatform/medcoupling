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

#ifndef __MEDFILEDATA_HXX__
#define __MEDFILEDATA_HXX__

#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDFileParameter.hxx"
#include "MEDFileField.hxx"
#include "MEDFileMesh.hxx"

namespace ParaMEDMEM
{
  /*!
   * User class.
   */
  class MEDFileData : public RefCountObject, public MEDFileWritable
  {
  public:
    MEDLOADER_EXPORT static MEDFileData *New(const std::string& fileName);
    MEDLOADER_EXPORT static MEDFileData *New();
    MEDLOADER_EXPORT MEDFileData *deepCpy() const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDLOADER_EXPORT MEDFileFields *getFields() const;
    MEDLOADER_EXPORT MEDFileMeshes *getMeshes() const;
    MEDLOADER_EXPORT MEDFileParameters *getParams() const;
    MEDLOADER_EXPORT void setFields(MEDFileFields *fields);
    MEDLOADER_EXPORT void setMeshes(MEDFileMeshes *meshes);
    MEDLOADER_EXPORT void setParams(MEDFileParameters *params);
    MEDLOADER_EXPORT int getNumberOfFields() const;
    MEDLOADER_EXPORT int getNumberOfMeshes() const;
    MEDLOADER_EXPORT int getNumberOfParams() const;
    MEDLOADER_EXPORT std::string simpleRepr() const;
    //
    MEDLOADER_EXPORT bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT bool changeMeshName(const std::string& oldMeshName, const std::string& newMeshName);
    MEDLOADER_EXPORT bool unPolyzeMeshes();
    //
    MEDLOADER_EXPORT void write(const std::string& fileName, int mode) const;
  private:
    MEDFileData();
    MEDFileData(const std::string& fileName);
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDFileFields> _fields;
    MEDCouplingAutoRefCountObjectPtr<MEDFileMeshes> _meshes;
    MEDCouplingAutoRefCountObjectPtr<MEDFileParameters> _params;
  };
}

#endif
