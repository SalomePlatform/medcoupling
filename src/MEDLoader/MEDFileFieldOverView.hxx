// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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

#ifndef __MEDFILEFIELDOVERVIEW_HXX__
#define __MEDFILEFIELDOVERVIEW_HXX__

#include "MEDLoaderDefines.hxx"

#include "MCAuto.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDCoupling1GTUMesh.hxx"

#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"

#include <vector>

namespace MEDCoupling
{
  class DataArrayInt;
  class MEDCouplingMesh;
  class MEDFileMesh;
  class MEDFileUMesh;
  class MEDFileCMesh;
  class MEDFileStructuredMesh;
  class MEDFileCurveLinearMesh;
  class MEDFileFieldGlobs;
  class MEDFileFieldGlobsReal;
  class MEDFileAnyTypeField1TS;
  class MEDFileAnyTypeFieldMultiTS;

  class MEDFileMeshStruct : public RefCountObject
  {
  public:
    MEDLOADER_EXPORT static MEDFileMeshStruct *New(const MEDFileMesh *mesh);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    const MEDFileMesh *getTheMesh() const { return _mesh; }
    mcIdType getNumberOfNodes() const { return _nb_nodes; }
    bool doesManageGeoType(INTERP_KERNEL::NormalizedCellType t) const;
    mcIdType getNumberOfElemsOfGeoType(INTERP_KERNEL::NormalizedCellType t) const;
    int getLevelOfGeoType(INTERP_KERNEL::NormalizedCellType t) const;
    int getNumberOfLevs() const;
    int getNumberOfGeoTypesInLev(int relativeLev) const;
    // non const methods
    void appendIfImplicitType(INTERP_KERNEL::NormalizedCellType t);
  private:
    MEDFileMeshStruct(const MEDFileMesh *mesh);
  private:
    const MEDFileMesh *_mesh;
    std::string _name;
    mcIdType _nb_nodes;
    std::vector< std::vector<mcIdType> > _geo_types_distrib;
  }; 

  class MEDFileField1TSStructItem;

  class MEDMeshMultiLev : public RefCountObject
  {
  public:
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
  public:
    static MEDMeshMultiLev *New(const MEDFileMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayIdType *>& pfls, const std::vector<mcIdType>& nbEntities);
    static MEDMeshMultiLev *New(const MEDFileMesh *m, const std::vector<int>& levs);
    static MEDMeshMultiLev *NewOnlyOnNode(const MEDFileMesh *m, const DataArrayIdType *pflOnNode);
    void setNodeReduction(const DataArrayIdType *nr);
    void setCellReduction(const DataArrayIdType *cr);
    bool isFastlyTheSameStruct(const MEDFileField1TSStructItem& fst, const MEDFileFieldGlobsReal *globs) const;
    MEDLOADER_EXPORT DataArray *buildDataArray(const MEDFileField1TSStructItem& fst, const MEDFileFieldGlobsReal *globs, const DataArray *vals) const;
    MEDLOADER_EXPORT void retrieveFamilyIdsOnCells(DataArrayIdType *& famIds, bool& isWithoutCopy) const;
    MEDLOADER_EXPORT void retrieveNumberIdsOnCells(DataArrayIdType *& numIds, bool& isWithoutCopy) const;
    MEDLOADER_EXPORT void retrieveFamilyIdsOnNodes(DataArrayIdType *& famIds, bool& isWithoutCopy) const;
    MEDLOADER_EXPORT void retrieveNumberIdsOnNodes(DataArrayIdType *& numIds, bool& isWithoutCopy) const;
    MEDLOADER_EXPORT DataArrayIdType *retrieveGlobalNodeIdsIfAny() const;
    MEDLOADER_EXPORT std::vector< INTERP_KERNEL::NormalizedCellType > getGeoTypes() const;
    void setFamilyIdsOnCells(DataArrayIdType *famIds);
    void setNumberIdsOnCells(DataArrayIdType *numIds);
    void setFamilyIdsOnNodes(DataArrayIdType *famIds);
    void setNumberIdsOnNodes(DataArrayIdType *numIds);
    virtual void selectPartOfNodes(const DataArrayIdType *pflNodes) = 0;
    virtual MEDMeshMultiLev *prepare() const = 0;
    mcIdType getNumberOfCells(INTERP_KERNEL::NormalizedCellType t) const;
    mcIdType getNumberOfNodes() const;
  protected:
    std::string getPflNameOfId(int id) const;
    DataArray *constructDataArray(const MEDFileField1TSStructItem& fst, const MEDFileFieldGlobsReal *globs, const DataArray *vals) const;
    virtual void appendVertices(const DataArrayIdType *verticesToAdd, DataArrayIdType *nr);
  protected:
    MEDMeshMultiLev(const MEDFileMesh *mesh);
    MEDMeshMultiLev(const MEDMeshMultiLev& other);
    MEDMeshMultiLev(const MEDFileMesh *mesh, mcIdType nbNodes, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayIdType *>& pfls, const std::vector<mcIdType>& nbEntities);
  protected:
    const MEDFileMesh *_mesh;
    std::vector< MCAuto<DataArrayIdType> > _pfls;
    std::vector< INTERP_KERNEL::NormalizedCellType > _geo_types;
    std::vector<mcIdType> _nb_entities;
    MCAuto<DataArrayIdType> _node_reduction;
    mcIdType _nb_nodes;
    //
    MCAuto<DataArrayIdType> _cell_fam_ids;
    MCAuto<DataArrayIdType> _cell_num_ids;
    MCAuto<DataArrayIdType> _node_fam_ids;
    MCAuto<DataArrayIdType> _node_num_ids;
  public:
    MEDLOADER_EXPORT static const int PARAMEDMEM_2_VTKTYPE_LGTH=MEDCOUPLING2VTKTYPETRADUCER_LGTH;
    MEDLOADER_EXPORT static const unsigned char *PARAMEDMEM_2_VTKTYPE;
    MEDLOADER_EXPORT static const unsigned char HEXA27_PERM_ARRAY[27];
  };

  class MEDStructuredMeshMultiLev;

  class MEDUMeshMultiLev : public MEDMeshMultiLev
  {
  public:
    static MEDUMeshMultiLev *New(const MEDFileUMesh *m, const std::vector<int>& levs);
    static MEDUMeshMultiLev *New(const MEDFileUMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayIdType *>& pfls, const std::vector<mcIdType>& nbEntities);
    void selectPartOfNodes(const DataArrayIdType *pflNodes);
    MEDMeshMultiLev *prepare() const;
    MEDUMeshMultiLev(const MEDStructuredMeshMultiLev& other, const MCAuto<MEDCoupling1GTUMesh>& part);
    MEDLOADER_EXPORT bool buildVTUArrays(DataArrayDouble *& coords, DataArrayByte *&types, DataArrayIdType *&cellLocations, DataArrayIdType *& cells, DataArrayIdType *&faceLocations, DataArrayIdType *&faces) const;
  protected:
    void appendVertices(const DataArrayIdType *verticesToAdd, DataArrayIdType *nr);
  private:
    void reorderNodesIfNecessary(MCAuto<DataArrayDouble>& coords, DataArrayIdType *nodalConnVTK, DataArrayIdType *polyhedNodalConnVTK) const;
  private:
    MEDUMeshMultiLev(const MEDUMeshMultiLev& other);
    MEDUMeshMultiLev(const MEDFileUMesh *m, const std::vector<int>& levs);
    MEDUMeshMultiLev(const MEDFileUMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayIdType *>& pfls, const std::vector<mcIdType>& nbEntities);
  private:
    std::vector< MCAuto<MEDCoupling1GTUMesh> > _parts;
    //! this attribute is used only for mesh with no cells but having coordinates. For classical umeshes those pointer is equal to pointer of coordinates of instances in this->_parts.
    MCAuto<DataArrayDouble> _coords;
  };

  class MEDStructuredMeshMultiLev : public MEDMeshMultiLev
  {
  public:
    void selectPartOfNodes(const DataArrayIdType *pflNodes);
    virtual std::vector<mcIdType> getNodeGridStructure() const = 0;
  protected:
    MEDStructuredMeshMultiLev(const MEDStructuredMeshMultiLev& other);
    MEDStructuredMeshMultiLev(const MEDFileStructuredMesh *m, const std::vector<int>& lev);
    MEDStructuredMeshMultiLev(const MEDFileStructuredMesh *m, mcIdType nbOfNodes, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayIdType *>& pfls, const std::vector<mcIdType>& nbEntities);
    void dealWithImplicitUnstructuredMesh(const MEDFileMesh *m);
  protected:
    void moveFaceToCell() const;
    bool prepareForImplicitUnstructuredMeshCase(MEDMeshMultiLev *&ret) const;
  private:
    void initStdFieldOfIntegers(const MEDFileStructuredMesh *m);
  protected:
    bool _is_internal;
    MCAuto<DataArrayIdType> _face_fam_ids;
    MCAuto<DataArrayIdType> _face_num_ids;
  };

  class MEDCMeshMultiLev : public MEDStructuredMeshMultiLev
  {
  public:
    static MEDCMeshMultiLev *New(const MEDFileCMesh *m, const std::vector<int>& levs);
    static MEDCMeshMultiLev *New(const MEDFileCMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayIdType *>& pfls, const std::vector<mcIdType>& nbEntities);
    std::vector<mcIdType> getNodeGridStructure() const;
    MEDMeshMultiLev *prepare() const;
    MEDLOADER_EXPORT std::vector< DataArrayDouble * > buildVTUArrays(bool& isInternal) const;
  private:
    MEDCMeshMultiLev(const MEDCMeshMultiLev& other);
    MEDCMeshMultiLev(const MEDFileCMesh *m, const std::vector<int>& levs);
    MEDCMeshMultiLev(const MEDFileCMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayIdType *>& pfls, const std::vector<mcIdType>& nbEntities);
  private:
    std::vector< MCAuto<DataArrayDouble> > _coords;
  };

  class MEDCurveLinearMeshMultiLev : public MEDStructuredMeshMultiLev
  {
  public:
    static MEDCurveLinearMeshMultiLev *New(const MEDFileCurveLinearMesh *m, const std::vector<int>& levs);
    static MEDCurveLinearMeshMultiLev *New(const MEDFileCurveLinearMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayIdType *>& pfls , const std::vector<mcIdType>& nbEntities);
    std::vector<mcIdType> getNodeGridStructure() const;
    MEDMeshMultiLev *prepare() const;
    MEDLOADER_EXPORT void buildVTUArrays(DataArrayDouble *&coords, std::vector<mcIdType>& nodeStrct, bool& isInternal) const;
  private:
    MEDCurveLinearMeshMultiLev(const MEDCurveLinearMeshMultiLev& other);
    MEDCurveLinearMeshMultiLev(const MEDFileCurveLinearMesh *m, const std::vector<int>& levs);
    MEDCurveLinearMeshMultiLev(const MEDFileCurveLinearMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayIdType *>& pfls, const std::vector<mcIdType>& nbEntities);
  private:
    MCAuto<DataArrayDouble> _coords;
    std::vector<mcIdType> _structure;
  };

  class MEDFileField1TSStructItem2 : public BigMemoryObject
  {
  public:
    MEDFileField1TSStructItem2();
    MEDFileField1TSStructItem2(INTERP_KERNEL::NormalizedCellType a, const std::pair<mcIdType,mcIdType>& b, const std::string& pfl, const std::string& loc);
    void checkWithMeshStructForCells(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs);
    void checkWithMeshStructForGaussNE(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs);
    void checkWithMeshStructForGaussPT(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs);
    //
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    //
    const DataArrayIdType *getPfl(const MEDFileFieldGlobsReal *globs) const;
    INTERP_KERNEL::NormalizedCellType getGeo() const { return _geo_type; }
    mcIdType getNbEntity() const { return _nb_of_entity; }
    const std::pair<mcIdType,mcIdType>& getStartStop() const { return _start_end; }
    std::string getPflName() const;
    int getNbOfIntegrationPts(const MEDFileFieldGlobsReal *globs) const;
    //! warning this method also set _nb_of_entity attribute !
    void checkInRange(mcIdType nbOfEntity, int nip, const MEDFileFieldGlobsReal *globs);
    bool isFastlyEqual(mcIdType& startExp, INTERP_KERNEL::NormalizedCellType gt, const std::string& pflName) const;
    bool operator==(const MEDFileField1TSStructItem2& other) const;
    bool isCellSupportEqual(const MEDFileField1TSStructItem2& other, const MEDFileFieldGlobsReal *globs) const;
    bool isNodeSupportEqual(const MEDFileField1TSStructItem2& other, const MEDFileFieldGlobsReal *globs) const;
    static MEDFileField1TSStructItem2 BuildAggregationOf(const std::vector<const MEDFileField1TSStructItem2 *>& objs, const MEDFileFieldGlobsReal *globs);
  public:
    static const char NEWLY_CREATED_PFL_NAME[];
  private:
    INTERP_KERNEL::NormalizedCellType _geo_type;
    std::pair<mcIdType,mcIdType> _start_end;
    MCAuto<DataArrayIdType> _pfl;
    std::string _loc;
    mcIdType _nb_of_entity;
  };

  class MEDFileField1TSStructItem : public BigMemoryObject
  {
  public:
    MEDFileField1TSStructItem():_computed(false),_type(ON_CELLS) { }
    MEDFileField1TSStructItem(TypeOfField a, const std::vector< MEDFileField1TSStructItem2 >& b);
    void checkWithMeshStruct(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs);
    bool operator==(const MEDFileField1TSStructItem& other) const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    bool isEntityCell() const;
    bool isComputed() const { return _computed; }
    TypeOfField getType() const { return _type; }
    std::size_t getNumberOfItems() const { return _items.size(); }
    const MEDFileField1TSStructItem2& operator[](std::size_t i) const;
    //
    bool isCellSupportEqual(const MEDFileField1TSStructItem& other, const MEDFileFieldGlobsReal *globs) const;
    bool isNodeSupportEqual(const MEDFileField1TSStructItem& other, const MEDFileFieldGlobsReal *globs) const;
    MEDFileField1TSStructItem simplifyMeOnCellEntity(const MEDFileFieldGlobsReal *globs) const;
    bool isCompatibleWithNodesDiscr(const MEDFileField1TSStructItem& other, const MEDFileMeshStruct *meshSt, const MEDFileFieldGlobsReal *globs) const;
    bool isFullyOnOneLev(const MEDFileMeshStruct *meshSt, int& theFirstLevFull) const;
    std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypes(const MEDFileMesh *m) const;
    MEDLOADER_EXPORT MEDMeshMultiLev *buildFromScratchDataSetSupportOnCells(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs) const;
    MEDLOADER_EXPORT static MEDFileField1TSStructItem BuildItemFrom(const MEDFileAnyTypeField1TS *ref, const MEDFileMeshStruct *meshSt);
  private:
    bool _computed;
    TypeOfField _type;
    std::vector< MEDFileField1TSStructItem2 > _items;
  };

  class MEDFileField1TSStruct : public RefCountObject
  {
  public:
    static MEDFileField1TSStruct *New(const MEDFileAnyTypeField1TS *ref, MEDFileMeshStruct *mst);
    void checkWithMeshStruct(MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    bool isEqualConsideringThePast(const MEDFileAnyTypeField1TS *other, const MEDFileMeshStruct *mst) const;
    bool isSupportSameAs(const MEDFileAnyTypeField1TS *other, const MEDFileMeshStruct *meshSt);
    bool isCompatibleWithNodesDiscr(const MEDFileAnyTypeField1TS *other, const MEDFileMeshStruct *meshSt);
    MEDLOADER_EXPORT MEDMeshMultiLev *buildFromScratchDataSetSupport(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs) const;
    bool isDataSetSupportFastlyEqualTo(const MEDFileField1TSStruct& other, const MEDFileFieldGlobsReal *globs) const;
    std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypes(const MEDFileMesh *m) const;
  private:
    MEDFileField1TSStruct(const MEDFileAnyTypeField1TS *ref, MEDFileMeshStruct *mst);
    bool presenceOfCellDiscr(int& pos) const;
    bool presenceOfPartialNodeDiscr(int& pos) const;
  private:
    std::vector<MEDFileField1TSStructItem> _already_checked;
  };

  class MEDFileFastCellSupportComparator : public RefCountObject
  {
  public:
    MEDLOADER_EXPORT static MEDFileFastCellSupportComparator *New(const MEDFileMeshStruct *m, const MEDFileAnyTypeFieldMultiTS *ref);
    MEDLOADER_EXPORT MEDMeshMultiLev *buildFromScratchDataSetSupport(int timeStepId, const MEDFileFieldGlobsReal *globs) const;
    MEDLOADER_EXPORT bool isDataSetSupportEqualToThePreviousOne(int timeStepId, const MEDFileFieldGlobsReal *globs) const;
    MEDLOADER_EXPORT int getNumberOfTS() const;
    MEDLOADER_EXPORT std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypesAt(int timeStepId, const MEDFileMesh *m) const;
    bool isEqual(const MEDFileAnyTypeFieldMultiTS *other);
    bool isCompatibleWithNodesDiscr(const MEDFileAnyTypeFieldMultiTS *other);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
  private:
    MEDFileFastCellSupportComparator(const MEDFileMeshStruct *m, const MEDFileAnyTypeFieldMultiTS *ref);
  private:
    MCAuto<MEDFileMeshStruct> _mesh_comp;
    std::vector< MCAuto<MEDFileField1TSStruct> > _f1ts_cmps;
  };
}

#endif
