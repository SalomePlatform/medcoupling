// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

#ifndef __MEDFILEDATA_HXX__
#define __MEDFILEDATA_HXX__

#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDFileField.hxx"
#include "MEDFileMesh.hxx"

namespace ParaMEDMEM
{
  /*!
   * User class.
   */
  class MEDLOADER_EXPORT MEDFileData : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileData *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileData *New();
    MEDFileFields *getFields() const;
    MEDFileMeshes *getMeshes() const;
    void setFields(MEDFileFields *fields) throw(INTERP_KERNEL::Exception);
    void setMeshes(MEDFileMeshes *meshes) throw(INTERP_KERNEL::Exception);
    int getNumberOfFields() const throw(INTERP_KERNEL::Exception);
    int getNumberOfMeshes() const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    //
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    bool changeMeshName(const char *oldMeshName, const char *newMeshName) throw(INTERP_KERNEL::Exception);
    //
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
  private:
    MEDFileData();
    MEDFileData(const char *fileName) throw(INTERP_KERNEL::Exception);
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDFileFields> _fields;
    MEDCouplingAutoRefCountObjectPtr<MEDFileMeshes> _meshes;
  };
}

#endif
