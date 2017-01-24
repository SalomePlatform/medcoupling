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

#ifndef __MEDFILEDATA_HXX__
#define __MEDFILEDATA_HXX__

#include "MCAuto.hxx"
#include "MEDFileParameter.hxx"
#include "MEDFileField.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileMeshSupport.hxx"
#include "MEDFileStructureElement.hxx"

namespace MEDCoupling
{
  /*!
   * User class.
   */
  class MEDFileData : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT static MEDFileData *New(const std::string& fileName);
    MEDLOADER_EXPORT static MEDFileData *New(med_idt fid);
    MEDLOADER_EXPORT static MEDFileData *New();
    MEDLOADER_EXPORT static MEDFileData *New(DataArrayByte *db) { return BuildFromMemoryChunk<MEDFileData>(db); }
    MEDLOADER_EXPORT MEDFileData *deepCopy() const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
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
    MEDLOADER_EXPORT std::string getHeader() const;
    MEDLOADER_EXPORT void setHeader(const std::string& header);
    //
    MEDLOADER_EXPORT bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT bool changeMeshName(const std::string& oldMeshName, const std::string& newMeshName);
    MEDLOADER_EXPORT bool unPolyzeMeshes();
    MEDLOADER_EXPORT void dealWithStructureElements();
    MEDLOADER_EXPORT static MCAuto<MEDFileData> Aggregate(const std::vector<const MEDFileData *>& mfds);
    //
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
  private:
    MEDFileData();
    MEDFileData(med_idt fid);
    void readHeader(med_idt fid);
    void writeHeader(med_idt fid) const;
    void readMeshSupports(med_idt fid);
  private:
    MCAuto<MEDFileFields> _fields;
    MCAuto<MEDFileMeshes> _meshes;
    MCAuto<MEDFileParameters> _params;
    MCAuto<MEDFileMeshSupports> _mesh_supports;
    MCAuto<MEDFileStructureElements> _struct_elems;
    std::string _header;
  };
}

#endif
