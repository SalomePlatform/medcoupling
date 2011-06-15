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
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> getCoordsFamily() const { return _fam_coords; }
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> getCoordsNum() const { return _num_coords; }
    static int getMeshIdFromName(med_idt fid, const char *mname) throw(INTERP_KERNEL::Exception);
    static void readFamiliesAndGrps(med_idt fid, const char *mname, std::map<std::string,int>& fams, std::map<std::string, std::vector<std::string> >& grps);
    static void writeFamiliesAndGrps(med_idt fid, const char *mname, const std::map<std::string,int>& fams, const std::map<std::string, std::vector<std::string> >& grps, int tooLongStrPol);
    static void writeCoords(med_idt fid, const char *mname, const DataArrayDouble *coords, const DataArrayInt *famCoords, const DataArrayInt *numCoords);
  private:
    void sortTypes();
  private:
    std::vector< std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> > > _per_type_mesh;
    MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> _coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam_coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num_coords;
  };

  class MEDFileUMeshL2CMesh : public MEDFileMeshL2
  {
  };

  class MEDFileUMeshSplitL1;

  class MEDFileUMeshPermCompute
  {
  public:
    MEDFileUMeshPermCompute(const MEDFileUMeshSplitL1* st);
    operator MEDCouplingUMesh *() const;
    void operator=(MEDCouplingUMesh *m);
    void updateTime() const;
  private:
    const MEDFileUMeshSplitL1 *_st;
    mutable unsigned int _mpt_time;
    mutable unsigned int _num_time;
    mutable MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> _m;
  };

  class MEDFileUMeshSplitL1 : public RefCountObject
  {
    friend class MEDFileUMeshPermCompute;
  public:
    MEDFileUMeshSplitL1(const MEDFileUMeshL2& l2, const char *mName, int id);
    MEDFileUMeshSplitL1(MEDCouplingUMesh *m);
    MEDFileUMeshSplitL1(MEDCouplingUMesh *m, bool newOrOld);
    void assignMesh(MEDCouplingUMesh *m, bool newOrOld) throw(INTERP_KERNEL::Exception);
    bool empty() const;
    int getMeshDimension() const;
    int getSize() const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getFamilyPart(const std::vector<int>& ids, bool renum) const;
    DataArrayInt *getFamilyPartArr(const std::vector<int>& ids, bool renum) const;
    MEDCouplingUMesh *getWholeMesh(bool renum) const;
    const DataArrayInt *getFamilyField() const;
    const DataArrayInt *getNumberField() const;
    const DataArrayInt *getRevNumberField() const;
    void eraseFamilyField();
    void setGroupsFromScratch(const std::vector<const MEDCouplingUMesh *>& ms, std::map<std::string,int>& familyIds,
                              std::map<std::string, std::vector<std::string> >& groups) throw(INTERP_KERNEL::Exception);
    void write(med_idt fid, const char *mName, int mdim) const;
    //
    void setFamilyArr(DataArrayInt *famArr);
    void setRenumArr(DataArrayInt *renumArr);
    //
    static std::vector<int> getNewFamiliesNumber(int nb, const std::map<std::string,int>& families);
    static void traduceFamilyNumber(const std::vector< std::vector<int> >& fidsGrps, std::map<std::string,int>& familyIds,
                                    std::map<int,int>& famIdTrad, std::map<int,std::string>& newfams);
    static DataArrayInt *Renumber(const DataArrayInt *renum, const DataArrayInt *da);
    static MEDCouplingUMesh *Renumber2(const DataArrayInt *renum, MEDCouplingUMesh *m, const int *cellIds);
  private:
    MEDCouplingUMesh *renumIfNeeded(MEDCouplingUMesh *m, const int *cellIds) const;
    DataArrayInt *renumIfNeededArr(const DataArrayInt *da) const;
    void computeRevNum() const;
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> _m_by_types;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num;
    mutable MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _rev_num;
    MEDFileUMeshPermCompute _m;
  };
}

#endif
