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

#ifndef __MEDFILEMESHLL_HXX__
#define __MEDFILEMESHLL_HXX__

#include "MEDFileBasis.hxx"
#include "MEDFileMeshElt.hxx"

#include "MEDCouplingAutoRefCountObjectPtr.hxx"

extern "C"
{
#include "med.h"
}

#include <map>

namespace ParaMEDMEM
{
  class MEDCouplingUMesh;

  class MEDFileMeshL2 : public RefCountObject
  {
  public:
    MEDFileMeshL2();
    const char *getName() const { return _name.getReprForWrite(); }
    const char *getDescription() const { return _description.getReprForWrite(); } 
  protected:
    MEDFileString _name;
    MEDFileString _description;
  };

  class MEDFileUMeshL2 : public MEDFileMeshL2
  {
  public:
    MEDFileUMeshL2();
    void loadAll(med_idt fid, int mId, const char *mName);
    void loadConnectivity(med_idt fid, int mdim, const char *mName);
    void loadCoords(med_idt fid, int mId, int mdim, const char *mName) throw(INTERP_KERNEL::Exception);
    int getNumberOfLevels() const { return _per_type_mesh.size(); }
    bool emptyLev(int levId) const { return _per_type_mesh[levId].empty(); }
    const std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> >& getLev(int levId) const { return _per_type_mesh[levId]; }
    bool isFamDefinedOnLev(int levId) const;
    bool isNumDefinedOnLev(int levId) const;
    MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> getCoords() const { return _coords; }
    static int getMeshIdFromName(med_idt fid, const char *mname) throw(INTERP_KERNEL::Exception);
    static void readFamiliesAndGrps(med_idt fid, const char *mname, std::map<std::string,int>& fams, std::map<std::string, std::vector<std::string> >& grps);
    static void writeFamiliesAndGrps(med_idt fid, const char *mname, const std::map<std::string,int>& fams, const std::map<std::string, std::vector<std::string> >& grps);
    static void writeCoords(med_idt fid, const char *mname, const DataArrayDouble *coords);
  private:
    void sortTypes();
  private:
    std::vector< std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> > > _per_type_mesh;
    MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> _coords;
  };

  class MEDFileUMeshL2CMesh : public MEDFileMeshL2
  {
  };

  class MEDFileUMeshSplitL1 : public RefCountObject
  {
  public:
    MEDFileUMeshSplitL1(const MEDFileUMeshL2& l2, const char *mName, int id);
    MEDFileUMeshSplitL1(MEDCouplingUMesh *m);
    bool empty() const;
    int getMeshDimension() const;
    MEDCouplingUMesh *getFamilyPart(const std::vector<int>& ids) const;
    MEDCouplingUMesh *getWholeMesh() const;
    void setGroupsFromScratch(const std::vector<MEDCouplingUMesh *>& ms, std::map<std::string,int>& familyIds,
                              std::map<std::string, std::vector<std::string> >& groups) throw(INTERP_KERNEL::Exception);
    void write(med_idt fid, const char *mName, int mdim) const;
    static std::vector<int> getNewFamiliesNumber(int nb, const std::map<std::string,int>& families);
    static void traduceFamilyNumber(const std::vector< std::vector<int> >& fidsGrps, std::map<std::string,int>& familyIds,
                                    std::map<int,int>& famIdTrad, std::map<int,std::string>& newfams);
  private:
    MEDCouplingUMesh *renumIfNeeded(MEDCouplingUMesh *m, const int *cellIds) const;
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> _m_by_types;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num;
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> _m;
  };
}

#endif
