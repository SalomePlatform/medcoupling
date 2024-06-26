// Copyright (C) 2007-2024  CEA, EDF
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
    static void ReadFamiliesAndGrps(med_idt fid, const std::string& mname, std::map<std::string,mcIdType>& fams, std::map<std::string, std::vector<std::string> >& grps, MEDFileMeshReadSelector *mrs);
    static void WriteFamiliesAndGrps(med_idt fid, const std::string& mname, const std::map<std::string,mcIdType>& fams, const std::map<std::string, std::vector<std::string> >& grps, int tooLongStrPol);
    static bool RenameFamiliesFromFileToMem(std::vector< std::string >& famNames);
    static bool RenameFamiliesFromMemToFile(std::vector< std::string >& famNames);
    static MEDCoupling::MEDCouplingAxisType TraduceAxisType(med_axis_type at);
    static MEDCoupling::MEDCouplingAxisType TraduceAxisTypeStruct(med_grid_type gt);
    static med_axis_type TraduceAxisTypeRev(MEDCoupling::MEDCouplingAxisType at);
    static med_grid_type TraduceAxisTypeRevStruct(MEDCoupling::MEDCouplingAxisType at);
  private:
    typedef bool (*RenameFamiliesPatternFunc)(std::vector< std::string >&);
    static void RenameFamiliesPatternInternal(std::vector< std::pair<std::string,std::pair<mcIdType,std::vector<std::string> > > >& crudeFams, RenameFamiliesPatternFunc func);
    static void RenameFamiliesFromFileToMemInternal(std::vector< std::pair<std::string,std::pair<mcIdType,std::vector<std::string> > > >& crudeFams);
    static void RenameFamiliesFromMemToFileInternal(std::vector< std::pair<std::string,std::pair<mcIdType,std::vector<std::string> > > >& crudeFams);
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
    std::string getClassName() const override { return std::string("MEDFileUMeshL2"); }
    std::vector<std::string> loadCommonPart(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, int dt, int it, int& Mdim);
    void loadAll(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadPart(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<mcIdType>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadPartFromUserDistrib(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, const std::map<INTERP_KERNEL::NormalizedCellType,std::vector<mcIdType>>& distrib, int dt, int it, MEDFileMeshReadSelector *mrs);
    void dealWithCoordsInLoadPart(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, const std::vector<std::string>& infosOnComp, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<mcIdType>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs);
    std::vector<std::string> loadPartConnectivityOnly(med_idt fid, const MeshOrStructMeshCls *mId, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<mcIdType>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs, int& Mdim);
    void loadConnectivity(med_idt fid, int mdim, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadPartOfConnectivity(med_idt fid, int mdim, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<mcIdType>& slicPerTyp, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadPartOfConnectivityFromUserDistrib(med_idt fid, int mdim, const std::string& mName, const std::map<INTERP_KERNEL::NormalizedCellType,std::vector<mcIdType>>& distrib, int dt, int it, MEDFileMeshReadSelector *mrs);

    void loadCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it);
    void loadPartCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, mcIdType nMin, mcIdType nMax);
    void loadPartCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, const std::vector<mcIdType>& distribNodes);
    void loadPartCoordsSlice(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, const DataArrayIdType *nodeIds, mcIdType nbOfCoordLS);
    int getNumberOfLevels() const { return (int)_per_type_mesh.size(); }
    bool emptyLev(int levId) const { return _per_type_mesh[levId].empty(); }
    const std::vector< MCAuto<MEDFileUMeshPerType> >& getLev(int levId) const { return _per_type_mesh[levId]; }
    bool isFamDefinedOnLev(int levId) const;
    bool isNumDefinedOnLev(int levId) const;
    bool isNamesDefinedOnLev(int levId) const;
    MCAuto<DataArrayDouble> getCoords() const { return _coords; }
    MCAuto<DataArrayIdType> getCoordsFamily() const { return _fam_coords; }
    MCAuto<DataArrayIdType> getCoordsNum() const { return _num_coords; }
    MCAuto<DataArrayIdType> getCoordsGlobalNum() const { return _global_num_coords; }
    MCAuto<DataArrayAsciiChar> getCoordsName() const { return _name_coords; }
    static void WriteCoords(med_idt fid, const std::string& mname, int dt, int it, double time, const DataArrayDouble *coords, const DataArrayIdType *famCoords, const DataArrayIdType *numCoords, const DataArrayAsciiChar *nameCoords, const DataArrayIdType *globalNumCoords);
    static void LoadPartCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, const std::vector<mcIdType>& distribNodes,
    MCAuto<DataArrayDouble>& _coords, MCAuto<PartDefinition>& _part_coords, MCAuto<DataArrayIdType>& _fam_coords, MCAuto<DataArrayIdType>& _num_coords, MCAuto<DataArrayAsciiChar>& _name_coords);
    static void LoadPartCoords(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, mcIdType nMin, mcIdType nMax,
    MCAuto<DataArrayDouble>& _coords, MCAuto<PartDefinition>& _part_coords, MCAuto<DataArrayIdType>& _fam_coords, MCAuto<DataArrayIdType>& _num_coords, MCAuto<DataArrayAsciiChar>& _name_coords);
    static void LoadPartCoordsArray(med_idt fid, const std::vector<std::string>& infosOnComp, const std::string& mName, int dt, int it, const DataArrayIdType *nodeIds,
MCAuto<DataArrayDouble>& _coords, MCAuto<DataArrayIdType>& _fam_coords, MCAuto<DataArrayIdType>& _num_coords, MCAuto<DataArrayAsciiChar>& _name_coords);
    static void allocCoordsPartCoords(mcIdType spaceDim, mcIdType nMin, mcIdType nMax, MCAuto<DataArrayDouble>& _coords, MCAuto<PartDefinition>& _part_coords);
    static void allocCoordsPartCoords(mcIdType spaceDim, const std::vector<mcIdType>& nodeIds, MCAuto<DataArrayDouble>& _coords, MCAuto<PartDefinition>& _part_coords);
    static void fillPartCoords(med_idt fid, mcIdType spaceDim, const std::string& mName, int dt, int it, const PartDefinition *partCoords, MCAuto<DataArrayDouble>& _coords, MCAuto<DataArrayIdType>& _fam_coords, MCAuto<DataArrayIdType>& _num_coords, MCAuto<DataArrayAsciiChar>& _name_coords);
  private:
    void sortTypes();
  private:
    std::vector< std::vector< MCAuto<MEDFileUMeshPerType> > > _per_type_mesh;
    MCAuto<DataArrayDouble> _coords;
    MCAuto<DataArrayIdType> _fam_coords;
    MCAuto<DataArrayIdType> _num_coords;
    MCAuto<DataArrayIdType> _global_num_coords;
    MCAuto<DataArrayAsciiChar> _name_coords;
  };

  class MEDFileStrMeshL2 : public MEDFileMeshL2
  {
  };

  class MEDFileCMeshL2 : public MEDFileStrMeshL2
  {
  public:
    MEDFileCMeshL2();
    std::string getClassName() const override { return std::string("MEDFileCMeshL2"); }
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
    std::string getClassName() const override { return std::string("MEDFileCLMeshL2"); }
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
    std::string getClassName() const override { return std::string("MEDFileUMeshPermCompute"); }
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
    std::string getClassName() const override { return std::string("MEDFileUMeshAggregateCompute"); }
    void setName(const std::string& name);
    void assignParts(const std::vector< const MEDCoupling1GTUMesh * >& mParts);
    void assignDefParts(const std::vector<const PartDefinition *>& partDefs);
    void assignUMesh(MEDCouplingUMesh *m);
    MEDCouplingUMesh *getUmesh() const;
    mcIdType getNumberOfCells() const;
    void highlightUsedNodes(std::vector<bool>& nodesToBeHighlighted) const;
    std::vector<MEDCoupling1GTUMesh *> getParts() const;
    std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypes() const;
    mcIdType getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType ct) const;
    std::vector<MEDCoupling1GTUMesh *> retrievePartsWithoutComputation() const;
    MEDCoupling1GTUMesh *retrievePartWithoutComputation(INTERP_KERNEL::NormalizedCellType gt) const;
    void getStartStopOfGeoTypeWithoutComputation(INTERP_KERNEL::NormalizedCellType gt, mcIdType& start, mcIdType& stop) const;
    void renumberNodesInConnWithoutComputation(const mcIdType *newNodeNumbersO2N);
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
    std::vector<mcIdType> getDistributionOfTypes() const;
    mcIdType getSize() const;
    void setCoords(DataArrayDouble *coords);
    void forceComputationOfPartsFromUMesh() const;
    void declarePartsUpdated() const;
    const PartDefinition *getPartDefOfWithoutComputation(INTERP_KERNEL::NormalizedCellType gt) const;
    void serialize(std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI) const;
    void unserialize(const std::string& name, DataArrayDouble *coo, std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI);
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
    std::string getClassName() const override { return std::string("MEDFileUMeshSplitL1"); }
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
    void declarePartsUpdated() const;
    bool empty() const;
    bool presenceOfOneFams(const std::vector<mcIdType>& ids) const;
    int getMeshDimension() const;
    void simpleRepr(std::ostream& oss) const;
    mcIdType getSize() const;
    MEDCouplingUMesh *getFamilyPart(const mcIdType *idsBg, const mcIdType *idsEnd, bool renum) const;
    DataArrayIdType *getFamilyPartArr(const mcIdType *idsBg, const mcIdType *idsEnd, bool renum) const;
    MEDCouplingUMesh *getWholeMesh(bool renum) const;
    mcIdType getNumberOfCells() const;
    bool isMeshStoredSplitByType() const { return _m_by_types.isStoredSplitByType(); }
    std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypes() const;
    mcIdType getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType ct) const;
    std::vector<MEDCoupling1GTUMesh *> getDirectUndergroundSingleGeoTypeMeshes() const { return _m_by_types.retrievePartsWithoutComputation(); }
    MEDCoupling1GTUMesh *getDirectUndergroundSingleGeoTypeMesh(INTERP_KERNEL::NormalizedCellType gt) const { return _m_by_types.retrievePartWithoutComputation(gt); }
    DataArrayIdType *extractFamilyFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    DataArrayIdType *extractNumberFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    std::vector<mcIdType> getDistributionOfTypes() const { return _m_by_types.getDistributionOfTypes(); }
    DataArrayIdType *getOrCreateAndGetFamilyField();
    const DataArrayIdType *getFamilyField() const;
    const DataArrayIdType *getNumberField() const;
    const DataArrayAsciiChar *getNameField() const;
    const DataArrayIdType *getRevNumberField() const;
    void highlightUsedNodes(std::vector<bool>& nodesToBeHighlighted) const { _m_by_types.highlightUsedNodes(nodesToBeHighlighted); }
    const PartDefinition *getPartDef(INTERP_KERNEL::NormalizedCellType gt) const;
    void eraseFamilyField();
    void setGroupsFromScratch(const std::vector<const MEDCouplingUMesh *>& ms, std::map<std::string,mcIdType>& familyIds,
                              std::map<std::string, std::vector<std::string> >& groups);
    void checkCoordsConsistency(const DataArrayDouble *coords) const;
    void write(med_idt fid, const std::string& mName, int mdim) const;
    //
    void setFamilyArr(DataArrayIdType *famArr);
    DataArrayIdType *getFamilyField();
    void setRenumArr(DataArrayIdType *renumArr);
    void setNameArr(DataArrayAsciiChar *nameArr);
    void changeFamilyIdArr(mcIdType oldId, mcIdType newId);
    //
    void renumberNodesInConn(const mcIdType *newNodeNumbersO2N);
    //
    void serialize(std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI) const;
    void unserialize(const std::string& name, DataArrayDouble *coo, std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI);
    //
    static void ClearNonDiscrAttributes(const MEDCouplingMesh *tmp);
    static std::vector<mcIdType> GetNewFamiliesNumber(mcIdType nb, const std::map<std::string,mcIdType>& families);
    static void TraduceFamilyNumber(const std::vector< std::vector<mcIdType> >& fidsGrps, std::map<std::string,mcIdType>& familyIds,
                                    std::map<mcIdType,mcIdType>& famIdTrad, std::map<mcIdType,std::string>& newfams);
    static DataArrayIdType *Renumber(const DataArrayIdType *renum, const DataArrayIdType *da);
    static MEDCouplingUMesh *Renumber2(const DataArrayIdType *renum, MEDCouplingUMesh *m, const mcIdType *cellIds);
    static MEDFileUMeshSplitL1 *Unserialize(const std::string& name, DataArrayDouble *coo, std::vector<mcIdType>& tinyInt, std::vector< MCAuto<DataArrayIdType> >& bigArraysI);
  private:
    MEDFileUMeshSplitL1();
    void assignCommonPart();
    MEDCouplingUMesh *renumIfNeeded(MEDCouplingUMesh *m, const mcIdType *cellIds) const;
    DataArrayIdType *renumIfNeededArr(const DataArrayIdType *da) const;
    void computeRevNum() const;
  private:
    MEDFileUMeshAggregateCompute _m_by_types;
    MCAuto<DataArrayIdType> _fam;
    MCAuto<DataArrayIdType> _num;
    MCAuto<DataArrayIdType> _global_num;
    MCAuto<DataArrayAsciiChar> _names;
    mutable MCAuto<DataArrayIdType> _rev_num;
    MEDFileUMeshPermCompute _m;
  };

  class MEDFileEltStruct4Mesh : public RefCountObject
  {
  public:
    static MEDFileEltStruct4Mesh *New(med_idt fid, const std::string& mName, int dt, int it, int iterOnStEltOfMesh, MEDFileMeshReadSelector *mrs);
    std::string getClassName() const override { return std::string("MEDFileEltStruct4Mesh"); }
    std::string getGeoTypeName() const { return _geo_type_name; }
    MCAuto<DataArrayIdType> getConn() const { return _conn; }
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
    MCAuto<DataArrayIdType> _conn;
    MCAuto<MEDFileUMeshPerTypeCommon> _common;
    std::vector< MCAuto<DataArray> > _vars;
  };
}

#endif
