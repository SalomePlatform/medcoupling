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

#ifndef __MEDFILEMESHLL_HXX__
#define __MEDFILEMESHLL_HXX__

#include "MEDFileBasis.hxx"
#include "MEDFileMeshElt.hxx"

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingPartDefinition.hxx"
#include "MEDCouplingCurveLinearMesh.hxx"
#include "MCAuto.hxx"

#include "InterpKernelAutoPtr.hxx"

#include "med.h"

#include <map>

namespace MEDCoupling
{
  class MEDFileMeshReadSelector;

  class MeshOrStructMeshCls
  {
  protected:
    MeshOrStructMeshCls(int mid):_mid(mid) { }
  public:
    virtual ~MeshOrStructMeshCls() {}
    int getID() const { return _mid; }
    virtual std::vector<std::string> getAxisInfoOnMesh(med_idt fid, const std::string& mName, MEDCoupling::MEDCouplingMeshType& meshType, MEDCoupling::MEDCouplingAxisType& axType, int& nstep, int& Mdim, MEDFileString& description, MEDFileString& dtunit, MEDFileString& univName) const = 0;
    virtual double checkMeshTimeStep(med_idt fid, const std::string& mName, int nstep, int dt, int it) const = 0;
  private:
    int _mid;
  };

  class MeshCls : public MeshOrStructMeshCls
  {
  public:
    MeshCls(int mid):MeshOrStructMeshCls(mid) { }
    std::vector<std::string> getAxisInfoOnMesh(med_idt fid, const std::string& mName, MEDCoupling::MEDCouplingMeshType& meshType, MEDCoupling::MEDCouplingAxisType& axType, int& nstep, int& Mdim, MEDFileString& description, MEDFileString& dtunit, MEDFileString& univName) const;
    double checkMeshTimeStep(med_idt fid, const std::string& mName, int nstep, int dt, int it) const;
  };

  class StructMeshCls : public MeshOrStructMeshCls
  {
  public:
    StructMeshCls(int mid):MeshOrStructMeshCls(mid) { }
    std::vector<std::string> getAxisInfoOnMesh(med_idt fid, const std::string& mName, MEDCoupling::MEDCouplingMeshType& meshType, MEDCoupling::MEDCouplingAxisType& axType, int& nstep, int& Mdim, MEDFileString& description, MEDFileString& dtunit, MEDFileString& univName) const;
    double checkMeshTimeStep(med_idt fid, const std::string& mName, int nstep, int dt, int it) const;
  };
  
  class MEDFileMeshL2 : public RefCountObject
  {
  public:
    MEDFileMeshL2();
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    const char *getName() const { return _name.getReprForWrite(); }
    const char *getDescription() const { return _description.getReprForWrite(); }
    const char *getUnivName() const { return _univ_name.getReprForWrite(); }
    const char *getTimeUnit() const { return _dt_unit.getReprForWrite(); }
    int getIteration() const { return _iteration; }
    int getOrder() const { return _order; }
    double getTime() const { return _time; }
    MCAuto<PartDefinition> getPartDefOfCoo() const { return _part_coords; }
    std::vector<std::string> getAxisInfoOnMesh(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, MEDCoupling::MEDCouplingMeshType& meshType, MEDCoupling::MEDCouplingAxisType& axType, int& nstep, int& Mdim);
    static INTERP_KERNEL::AutoCppPtr<MeshOrStructMeshCls> GetMeshIdFromName(med_idt fid, const std::string& mName, MEDCoupling::MEDCouplingMeshType& meshType, MEDCoupling::MEDCouplingAxisType& axType, int& dt, int& it, std::string& dtunit1);
    static void ReadFamiliesAndGrps(med_idt fid, const std::string& mname, std::map<std::string,int>& fams, std::map<std::string, std::vector<std::string> >& grps, MEDFileMeshReadSelector *mrs);
    static void WriteFamiliesAndGrps(med_idt fid, const std::string& mname, const std::map<std::string,int>& fams, const std::map<std::string, std::vector<std::string> >& grps, int tooLongStrPol);
    static bool RenameFamiliesFromFileToMem(std::vector< std::string >& famNames);
    static bool RenameFamiliesFromMemToFile(std::vector< std::string >& famNames);
    static MEDCoupling::MEDCouplingAxisType TraduceAxisType(med_axis_type at);
    static MEDCoupling::MEDCouplingAxisType TraduceAxisTypeStruct(med_grid_type gt);
    static med_axis_type TraduceAxisTypeRev(MEDCoupling::MEDCouplingAxisType at);
    static med_grid_type TraduceAxisTypeRevStruct(MEDCoupling::MEDCouplingAxisType at);
  private:
    typedef bool (*RenameFamiliesPatternFunc)(std::vector< std::string >&);
    static void RenameFamiliesPatternInternal(std::vector< std::pair<std::string,std::pair<int,std::vector<std::string> > > >& crudeFams, RenameFamiliesPatternFunc func);
    static void RenameFamiliesFromFileToMemInternal(std::vector< std::pair<std::string,std::pair<int,std::vector<std::string> > > >& crudeFams);
    static void RenameFamiliesFromMemToFileInternal(std::vector< std::pair<std::string,std::pair<int,std::vector<std::string> > > >& crudeFams);
  public:
    static const char ZE_SEP_FOR_FAMILY_KILLERS[];
    static int ZE_SEP2_FOR_FAMILY_KILLERS;
  protected:
    MEDFileString _name;
    MEDFileString _description;
    MEDFileString _univ_name;
    MEDFileString _dt_unit;
    int _iteration;
    int _order;
    double _time;
    MCAuto<PartDefinition> _part_coords;
  };

  class MEDFileUMeshL2 : public MEDFileMeshL2
  {
  public:
    MEDFileUMeshL2();
    std::vector<std::string> loadCommonPart(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, int dt, int it, int& Mdim);
    void loadAll(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadPart(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<int>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadConnectivity(med_idt fid, int mdim, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadPartOfConnectivity(med_idt fid, int mdim, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<int>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it);
    void loadPartCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, int nMin, int nMax);
    int getNumberOfLevels() const { return _per_type_mesh.size(); }
    bool emptyLev(int levId) const { return _per_type_mesh[levId].empty(); }
    const std::vector< MCAuto<MEDFileUMeshPerType> >& getLev(int levId) const { return _per_type_mesh[levId]; }
    bool isFamDefinedOnLev(int levId) const;
    bool isNumDefinedOnLev(int levId) const;
    bool isNamesDefinedOnLev(int levId) const;
    MCAuto<DataArrayDouble> getCoords() const { return _coords; }
    MCAuto<DataArrayInt> getCoordsFamily() const { return _fam_coords; }
    MCAuto<DataArrayInt> getCoordsNum() const { return _num_coords; }
    MCAuto<DataArrayInt> getCoordsGlobalNum() const { return _global_num_coords; }
    MCAuto<DataArrayAsciiChar> getCoordsName() const { return _name_coords; }
    static void WriteCoords(med_idt fid, const std::string& mname, int dt, int it, double time, const DataArrayDouble *coords, const DataArrayInt *famCoords, const DataArrayInt *numCoords, const DataArrayAsciiChar *nameCoords, const DataArrayInt *globalNumCoords);
  private:
    void sortTypes();
  private:
    std::vector< std::vector< MCAuto<MEDFileUMeshPerType> > > _per_type_mesh;
    MCAuto<DataArrayDouble> _coords;
    MCAuto<DataArrayInt> _fam_coords;
    MCAuto<DataArrayInt> _num_coords;
    MCAuto<DataArrayInt> _global_num_coords;
    MCAuto<DataArrayAsciiChar> _name_coords;
  };

  class MEDFileStrMeshL2 : public MEDFileMeshL2
  {
  };

  class MEDFileCMeshL2 : public MEDFileStrMeshL2
  {
  public:
    MEDFileCMeshL2();
    void loadAll(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, int dt, int it);
    MEDCouplingCMesh *getMesh() { return _cmesh; }
    MEDCoupling::MEDCouplingAxisType getAxisType() const { return _ax_type; }
  private:
    static med_data_type GetDataTypeCorrespondingToSpaceId(int id);
  private:
    MCAuto<MEDCouplingCMesh> _cmesh;
    MEDCoupling::MEDCouplingAxisType _ax_type;
  };

  class MEDFileCLMeshL2 : public MEDFileStrMeshL2
  {
  public:
    MEDFileCLMeshL2();
    void loadAll(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, int dt, int it);
    MEDCouplingCurveLinearMesh *getMesh() { return _clmesh; }
  private:
    MCAuto<MEDCouplingCurveLinearMesh> _clmesh;
  };

  class MEDFileMesh;
  class MEDFileUMeshSplitL1;

  class MEDFileUMeshPermCompute : public BigMemoryObject
  {
  public:
    MEDFileUMeshPermCompute(const MEDFileUMeshSplitL1* st);
    operator MEDCouplingUMesh *() const;
    void operator=(MEDCouplingUMesh *m);
    void updateTime() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
  private:
    const MEDFileUMeshSplitL1 *_st;
    mutable std::size_t _mpt_time;
    mutable std::size_t _num_time;
    mutable MCAuto<MEDCouplingUMesh> _m;
  };

  class MEDFileUMeshAggregateCompute : public BigMemoryObject
  {
  public:
    MEDFileUMeshAggregateCompute();
    void setName(const std::string& name);
    void assignParts(const std::vector< const MEDCoupling1GTUMesh * >& mParts);
    void assignDefParts(const std::vector<const PartDefinition *>& partDefs);
    void assignUMesh(MEDCouplingUMesh *m);
    MEDCouplingUMesh *getUmesh() const;
    int getNumberOfCells() const;
    std::vector<MEDCoupling1GTUMesh *> getParts() const;
    std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypes() const;
    int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType ct) const;
    std::vector<MEDCoupling1GTUMesh *> retrievePartsWithoutComputation() const;
    MEDCoupling1GTUMesh *retrievePartWithoutComputation(INTERP_KERNEL::NormalizedCellType gt) const;
    void getStartStopOfGeoTypeWithoutComputation(INTERP_KERNEL::NormalizedCellType gt, int& start, int& stop) const;
    void renumberNodesInConnWithoutComputation(const int *newNodeNumbersO2N);
    bool isStoredSplitByType() const;
    std::size_t getTimeOfThis() const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileUMeshAggregateCompute deepCopy(DataArrayDouble *coords) const;
    void shallowCpyMeshes();
    bool isEqual(const MEDFileUMeshAggregateCompute& other, double eps, std::string& what) const;
    void checkConsistency() const;
    void clearNonDiscrAttributes() const;
    void synchronizeTinyInfo(const MEDFileMesh& master) const;
    bool empty() const;
    int getMeshDimension() const;
    std::vector<int> getDistributionOfTypes() const;
    int getSize() const;
    void setCoords(DataArrayDouble *coords);
    void forceComputationOfPartsFromUMesh() const;
    const PartDefinition *getPartDefOfWithoutComputation(INTERP_KERNEL::NormalizedCellType gt) const;
    void serialize(std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI) const;
    void unserialize(const std::string& name, DataArrayDouble *coo, std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI);
  private:
    std::size_t getTimeOfParts() const;
    std::size_t getTimeOfUMesh() const;
  private:
    mutable std::vector< MCAuto<MEDCoupling1GTUMesh> > _m_parts;
    mutable std::size_t _mp_time;
    mutable std::size_t _m_time;
    mutable MCAuto<MEDCouplingUMesh> _m;
    mutable std::vector< MCAuto<PartDefinition> > _part_def;
  };

  class MEDFileUMeshSplitL1 : public RefCountObject
  {
    friend class MEDFileUMeshPermCompute;
    friend class MEDFileUMesh;
  public:
    MEDFileUMeshSplitL1(const MEDFileUMeshSplitL1& other);
    MEDFileUMeshSplitL1(const MEDFileUMeshL2& l2, const std::string& mName, int id);
    MEDFileUMeshSplitL1(MEDCoupling1GTUMesh *m);
    MEDFileUMeshSplitL1(MEDCouplingUMesh *m);
    MEDFileUMeshSplitL1(MEDCouplingUMesh *m, bool newOrOld);
    void setName(const std::string& name);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileUMeshSplitL1 *shallowCpyUsingCoords(DataArrayDouble *coords) const;
    MEDFileUMeshSplitL1 *deepCopy(DataArrayDouble *coords) const;
    void checkConsistency() const;
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
    int getNumberOfCells() const;
    bool isMeshStoredSplitByType() const { return _m_by_types.isStoredSplitByType(); }
    std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypes() const;
    int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType ct) const;
    std::vector<MEDCoupling1GTUMesh *> getDirectUndergroundSingleGeoTypeMeshes() const { return _m_by_types.retrievePartsWithoutComputation(); }
    MEDCoupling1GTUMesh *getDirectUndergroundSingleGeoTypeMesh(INTERP_KERNEL::NormalizedCellType gt) const { return _m_by_types.retrievePartWithoutComputation(gt); }
    DataArrayInt *extractFamilyFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    DataArrayInt *extractNumberFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    std::vector<int> getDistributionOfTypes() const { return _m_by_types.getDistributionOfTypes(); }
    DataArrayInt *getOrCreateAndGetFamilyField();
    const DataArrayInt *getFamilyField() const;
    const DataArrayInt *getNumberField() const;
    const DataArrayAsciiChar *getNameField() const;
    const DataArrayInt *getRevNumberField() const;
    const PartDefinition *getPartDef(INTERP_KERNEL::NormalizedCellType gt) const;
    void eraseFamilyField();
    void setGroupsFromScratch(const std::vector<const MEDCouplingUMesh *>& ms, std::map<std::string,int>& familyIds,
                              std::map<std::string, std::vector<std::string> >& groups);
    void write(med_idt fid, const std::string& mName, int mdim) const;
    //
    void setFamilyArr(DataArrayInt *famArr);
    DataArrayInt *getFamilyField();
    void setRenumArr(DataArrayInt *renumArr);
    void setNameArr(DataArrayAsciiChar *nameArr);
    void changeFamilyIdArr(int oldId, int newId);
    //
    void renumberNodesInConn(const int *newNodeNumbersO2N);
    //
    void serialize(std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI) const;
    void unserialize(const std::string& name, DataArrayDouble *coo, std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI);
    //
    static void ClearNonDiscrAttributes(const MEDCouplingMesh *tmp);
    static std::vector<int> GetNewFamiliesNumber(int nb, const std::map<std::string,int>& families);
    static void TraduceFamilyNumber(const std::vector< std::vector<int> >& fidsGrps, std::map<std::string,int>& familyIds,
                                    std::map<int,int>& famIdTrad, std::map<int,std::string>& newfams);
    static DataArrayInt *Renumber(const DataArrayInt *renum, const DataArrayInt *da);
    static MEDCouplingUMesh *Renumber2(const DataArrayInt *renum, MEDCouplingUMesh *m, const int *cellIds);
    static MEDFileUMeshSplitL1 *Unserialize(const std::string& name, DataArrayDouble *coo, std::vector<int>& tinyInt, std::vector< MCAuto<DataArrayInt> >& bigArraysI);
  private:
    MEDFileUMeshSplitL1();
    void assignCommonPart();
    MEDCouplingUMesh *renumIfNeeded(MEDCouplingUMesh *m, const int *cellIds) const;
    DataArrayInt *renumIfNeededArr(const DataArrayInt *da) const;
    void computeRevNum() const;
  private:
    MEDFileUMeshAggregateCompute _m_by_types;
    MCAuto<DataArrayInt> _fam;
    MCAuto<DataArrayInt> _num;
    MCAuto<DataArrayInt> _global_num;
    MCAuto<DataArrayAsciiChar> _names;
    mutable MCAuto<DataArrayInt> _rev_num;
    MEDFileUMeshPermCompute _m;
  };

  class MEDFileEltStruct4Mesh : public RefCountObject
  {
  public:
    static MEDFileEltStruct4Mesh *New(med_idt fid, const std::string& mName, int dt, int it, int iterOnStEltOfMesh, MEDFileMeshReadSelector *mrs);
    std::string getGeoTypeName() const { return _geo_type_name; }
    MCAuto<DataArrayInt> getConn() const { return _conn; }
    MCAuto<MEDFileUMeshPerTypeCommon> getMeshDef() const { return _common; }
    const std::vector< MCAuto<DataArray> >& getVars() const { return _vars; }
  private:
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const MEDCoupling::BigMemoryObject*> getDirectChildrenWithNull() const;
  private:
    ~MEDFileEltStruct4Mesh() { }
  private:
    MEDFileEltStruct4Mesh(med_idt fid, const std::string& mName, int dt, int it, int iterOnStEltOfMesh, MEDFileMeshReadSelector *mrs);
  private:
    std::string _geo_type_name;
    int _geo_type;
    MCAuto<DataArrayInt> _conn;
    MCAuto<MEDFileUMeshPerTypeCommon> _common;
    std::vector< MCAuto<DataArray> > _vars;
  };
}

#endif
