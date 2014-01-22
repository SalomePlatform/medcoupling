// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#ifndef __MEDFILEMESHLL_HXX__
#define __MEDFILEMESHLL_HXX__

#include "MEDFileBasis.hxx"
#include "MEDFileMeshElt.hxx"

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingCurveLinearMesh.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include "med.h"

#include <map>

namespace ParaMEDMEM
{
  class MEDFileMeshReadSelector;
  
  class MEDFileMeshL2 : public RefCountObject
  {
  public:
    MEDFileMeshL2();
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
    const char *getName() const { return _name.getReprForWrite(); }
    const char *getDescription() const { return _description.getReprForWrite(); }
    const char *getUnivName() const { return _univ_name.getReprForWrite(); }
    const char *getTimeUnit() const { return _dt_unit.getReprForWrite(); }
    int getIteration() const { return _iteration; }
    int getOrder() const { return _order; }
    double getTime() { return _time; }
    std::vector<std::string> getAxisInfoOnMesh(med_idt fid, int mId, const std::string& mName, ParaMEDMEM::MEDCouplingMeshType& meshType, int& nstep, int& Mdim);
    static int GetMeshIdFromName(med_idt fid, const std::string& mName, ParaMEDMEM::MEDCouplingMeshType& meshType, int& dt, int& it, std::string& dtunit1);
    static double CheckMeshTimeStep(med_idt fid, const std::string& mname, int nstep, int dt, int it);
    static void ReadFamiliesAndGrps(med_idt fid, const std::string& mname, std::map<std::string,int>& fams, std::map<std::string, std::vector<std::string> >& grps, MEDFileMeshReadSelector *mrs);
    static void WriteFamiliesAndGrps(med_idt fid, const std::string& mname, const std::map<std::string,int>& fams, const std::map<std::string, std::vector<std::string> >& grps, int tooLongStrPol);
  protected:
    MEDFileString _name;
    MEDFileString _description;
    MEDFileString _univ_name;
    MEDFileString _dt_unit;
    int _iteration;
    int _order;
    double _time;
  };

  class MEDFileUMeshL2 : public MEDFileMeshL2
  {
  public:
    MEDFileUMeshL2();
    void loadAll(med_idt fid, int mId, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadConnectivity(med_idt fid, int mdim, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadCoords(med_idt fid, int mId, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it);
    int getNumberOfLevels() const { return _per_type_mesh.size(); }
    bool emptyLev(int levId) const { return _per_type_mesh[levId].empty(); }
    const std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> >& getLev(int levId) const { return _per_type_mesh[levId]; }
    bool isFamDefinedOnLev(int levId) const;
    bool isNumDefinedOnLev(int levId) const;
    bool isNamesDefinedOnLev(int levId) const;
    MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> getCoords() const { return _coords; }
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> getCoordsFamily() const { return _fam_coords; }
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> getCoordsNum() const { return _num_coords; }
    MEDCouplingAutoRefCountObjectPtr<DataArrayAsciiChar> getCoordsName() const { return _name_coords; }
    static void WriteCoords(med_idt fid, const std::string& mname, int dt, int it, double time, const DataArrayDouble *coords, const DataArrayInt *famCoords, const DataArrayInt *numCoords, const DataArrayAsciiChar *nameCoords);
  private:
    void sortTypes();
  private:
    std::vector< std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshPerType> > > _per_type_mesh;
    MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> _coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam_coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num_coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayAsciiChar> _name_coords;
  };

  class MEDFileStrMeshL2 : public MEDFileMeshL2
  {
  };

  class MEDFileCMeshL2 : public MEDFileStrMeshL2
  {
  public:
    MEDFileCMeshL2();
    void loadAll(med_idt fid, int mId, const std::string& mName, int dt, int it);
    MEDCouplingCMesh *getMesh() { return _cmesh; }
  private:
    static med_data_type GetDataTypeCorrespondingToSpaceId(int id);
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingCMesh> _cmesh;
  };
  
  class MEDFileCLMeshL2 : public MEDFileStrMeshL2
  {
  public:
    MEDFileCLMeshL2();
    void loadAll(med_idt fid, int mId, const std::string& mName, int dt, int it);
    MEDCouplingCurveLinearMesh *getMesh() { return _clmesh; }
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingCurveLinearMesh> _clmesh;
  };

  class MEDFileMesh;
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
    mutable std::size_t _mpt_time;
    mutable std::size_t _num_time;
    mutable MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> _m;
  };

  class MEDFileUMeshAggregateCompute : public BigMemoryObject
  {
  public:
    MEDFileUMeshAggregateCompute();
    void assignParts(const std::vector< const MEDCoupling1GTUMesh * >& mParts);
    void assignUMesh(MEDCouplingUMesh *m);
    MEDCouplingUMesh *getUmesh() const;
    std::vector<MEDCoupling1GTUMesh *> getParts() const;
    std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypes() const;
    std::vector<MEDCoupling1GTUMesh *> getPartsWithoutComputation() const;
    MEDCoupling1GTUMesh *getPartWithoutComputation(INTERP_KERNEL::NormalizedCellType gt) const;
    void getStartStopOfGeoTypeWithoutComputation(INTERP_KERNEL::NormalizedCellType gt, int& start, int& stop) const;
    std::size_t getTimeOfThis() const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDFileUMeshAggregateCompute deepCpy(DataArrayDouble *coords) const;
    bool isEqual(const MEDFileUMeshAggregateCompute& other, double eps, std::string& what) const;
    void clearNonDiscrAttributes() const;
    void synchronizeTinyInfo(const MEDFileMesh& master) const;
    bool empty() const;
    int getMeshDimension() const;
    std::vector<int> getDistributionOfTypes() const;
    int getSize() const;
    void setCoords(DataArrayDouble *coords);
    void forceComputationOfPartsFromUMesh() const;
  private:
    std::size_t getTimeOfParts() const;
    std::size_t getTimeOfUMesh() const;
  private:
    mutable std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCoupling1GTUMesh> > _m_parts;
    mutable std::size_t _mp_time;
    mutable std::size_t _m_time;
    mutable MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> _m;
  };
  
  class MEDFileUMeshSplitL1 : public RefCountObject
  {
    friend class MEDFileUMeshPermCompute;
  public:
    MEDFileUMeshSplitL1(const MEDFileUMeshSplitL1& other);
    MEDFileUMeshSplitL1(const MEDFileUMeshL2& l2, const std::string& mName, int id);
    MEDFileUMeshSplitL1(MEDCoupling1GTUMesh *m);
    MEDFileUMeshSplitL1(MEDCouplingUMesh *m);
    MEDFileUMeshSplitL1(MEDCouplingUMesh *m, bool newOrOld);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDFileUMeshSplitL1 *deepCpy(DataArrayDouble *coords) const;
    void setCoords(DataArrayDouble *coords);
    bool isEqual(const MEDFileUMeshSplitL1 *other, double eps, std::string& what) const;
    void clearNonDiscrAttributes() const;
    void synchronizeTinyInfo(const MEDFileMesh& master) const;
    void assignMesh(MEDCouplingUMesh *m, bool newOrOld);
    void assignParts(const std::vector< const MEDCoupling1GTUMesh * >& mParts);
    void forceComputationOfParts() const;
    bool empty() const;
    bool presenceOfOneFams(const std::vector<int>& ids) const;
    int getMeshDimension() const;
    void simpleRepr(std::ostream& oss) const;
    int getSize() const;
    MEDCouplingUMesh *getFamilyPart(const int *idsBg, const int *idsEnd, bool renum) const;
    DataArrayInt *getFamilyPartArr(const int *idsBg, const int *idsEnd, bool renum) const;
    MEDCouplingUMesh *getWholeMesh(bool renum) const;
    std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypes() const;
    std::vector<MEDCoupling1GTUMesh *> getDirectUndergroundSingleGeoTypeMeshes() const { return _m_by_types.getPartsWithoutComputation(); }
    MEDCoupling1GTUMesh *getDirectUndergroundSingleGeoTypeMesh(INTERP_KERNEL::NormalizedCellType gt) const { return _m_by_types.getPartWithoutComputation(gt); }
    DataArrayInt *extractFamilyFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    DataArrayInt *extractNumberFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    std::vector<int> getDistributionOfTypes() const { return _m_by_types.getDistributionOfTypes(); }
    DataArrayInt *getOrCreateAndGetFamilyField();
    const DataArrayInt *getFamilyField() const;
    const DataArrayInt *getNumberField() const;
    const DataArrayAsciiChar *getNameField() const;
    const DataArrayInt *getRevNumberField() const;
    void eraseFamilyField();
    void setGroupsFromScratch(const std::vector<const MEDCouplingUMesh *>& ms, std::map<std::string,int>& familyIds,
                              std::map<std::string, std::vector<std::string> >& groups) throw(INTERP_KERNEL::Exception);
    void write(med_idt fid, const std::string& mName, int mdim) const;
    //
    void setFamilyArr(DataArrayInt *famArr);
    void setRenumArr(DataArrayInt *renumArr);
    void setNameArr(DataArrayAsciiChar *nameArr);
    void changeFamilyIdArr(int oldId, int newId);
    //
    void renumberNodesInConn(const int *newNodeNumbersO2N);
    //
    static void ClearNonDiscrAttributes(const MEDCouplingMesh *tmp);
    static std::vector<int> GetNewFamiliesNumber(int nb, const std::map<std::string,int>& families);
    static void TraduceFamilyNumber(const std::vector< std::vector<int> >& fidsGrps, std::map<std::string,int>& familyIds,
                                    std::map<int,int>& famIdTrad, std::map<int,std::string>& newfams);
    static DataArrayInt *Renumber(const DataArrayInt *renum, const DataArrayInt *da);
    static MEDCouplingUMesh *Renumber2(const DataArrayInt *renum, MEDCouplingUMesh *m, const int *cellIds);
  private:
    void assignCommonPart();
    MEDCouplingUMesh *renumIfNeeded(MEDCouplingUMesh *m, const int *cellIds) const;
    DataArrayInt *renumIfNeededArr(const DataArrayInt *da) const;
    void computeRevNum() const;
  private:
    MEDFileUMeshAggregateCompute _m_by_types;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num;
    MEDCouplingAutoRefCountObjectPtr<DataArrayAsciiChar> _names;
    mutable MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _rev_num;
    MEDFileUMeshPermCompute _m;
  };
}

#endif
