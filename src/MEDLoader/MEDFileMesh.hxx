//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#ifndef __MEDFILEMESH_HXX__
#define __MEDFILEMESH_HXX__

#include "MEDFileMeshLL.hxx"

#include <map>

namespace ParaMEDMEM
{
  class MEDFileMesh : public RefCountObject
  {
  public:
    static MEDFileMesh *New(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception);
    void setName(const char *name) { _name=name; }
    const char *getName() const { return _name.c_str(); }
    void setUnivName(const char *name) { _univ_name=name; }
    const char *getUnivName() const { return _univ_name.c_str(); }
    void setDescription(const char *name) { _desc_name=name; }
    const char *getDescription() const { return _desc_name.c_str(); }
  protected:
    std::string _name;
    std::string _univ_name;
    std::string _desc_name;
  };

  class MEDFileUMesh : public MEDFileMesh
  {
  public:
    static MEDFileUMesh *New(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception);
    static MEDFileUMesh *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileUMesh *New();
    ~MEDFileUMesh();
    //
    void removeGroup(const char *name) throw(INTERP_KERNEL::Exception);
    void removeFamily(const char *name) throw(INTERP_KERNEL::Exception);
    void changeGroupName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception);
    void changeFamilyName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception);
    //
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    int getNumberOfLevels() const { return _ms.size(); }
    int getNumberOfNonEmptyLevels() const;
    int getFamilyId(const char *name) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamiliesIds(const std::vector<std::string>& famNames) const throw(INTERP_KERNEL::Exception);
    std::string getFamilyNameGivenId(int id) const throw(INTERP_KERNEL::Exception);
    int getMeshDimension() const;
    std::vector<std::string> getFamiliesOnGroup(const char *name) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getGroupsNames() const;
    std::vector<std::string> getFamiliesNames() const;
    std::vector<int> getNonEmptyLevels() const;
    DataArrayDouble *getCoords() const;
    MEDCouplingUMesh *getGroup(int meshDimRelToMax, const char *grp, bool renum=true) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getGroupArr(int meshDimRelToMax, const char *grp, bool renum=true) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getGroups(int meshDimRelToMax, const std::vector<std::string>& grps, bool renum=true) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getGroupsArr(int meshDimRelToMax, const std::vector<std::string>& grps, bool renum=true) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getFamily(int meshDimRelToMax, const char *fam, bool renum=true) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getFamilyArr(int meshDimRelToMax, const char *fam, bool renum=true) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getFamilies(int meshDimRelToMax, const std::vector<std::string>& fams, bool renum=true) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getFamiliesArr(int meshDimRelToMax, const std::vector<std::string>& fams, bool renum=true) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getNodeGroupArr(const char *grp, bool renum=true) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getNodeGroupsArr(const std::vector<std::string>& grps, bool renum=true) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getNodeFamilyArr(const char *fam, bool renum=true) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum=true) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getMeshAtRank(int meshDimRelToMax, bool renum=true) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getRank0Mesh(bool renum=true) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getRankM1Mesh(bool renum=true) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getRankM2Mesh(bool renum=true) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getRankM3Mesh(bool renum=true) const throw(INTERP_KERNEL::Exception);
    //
    void setMeshAtRank(int meshDimRelToMax, MEDCouplingUMesh *m) throw(INTERP_KERNEL::Exception);
    void setGroupsFromScratch(int meshDimRelToMax, const std::vector<MEDCouplingUMesh *>& ms) throw(INTERP_KERNEL::Exception);
    void setGroupsOnSetMesh(int meshDimRelToMax, const std::vector<MEDCouplingUMesh *>& ms) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileUMesh();
    MEDFileUMesh(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception);
    const MEDFileUMeshSplitL1 *getMeshAtLevSafe(int meshDimRelToMax) const throw(INTERP_KERNEL::Exception);
    void checkMeshDimCoherency(int meshDim, int meshDimRelToMax) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *checkMultiMesh(const std::vector<MEDCouplingUMesh *>& ms) const throw(INTERP_KERNEL::Exception);
  private:
    std::map<std::string, std::vector<std::string> > _groups;
    std::map<std::string,int> _families;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> > _ms;
    MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> _coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam_coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num_coords;
    int _too_long_str;
  };


  class MEDFileCMesh : public MEDFileMesh
  {
  };
}

#endif
