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

#ifndef __MEDFILEFIELDOVERVIEW_HXX__
#define __MEDFILEFIELDOVERVIEW_HXX__

#include "MEDLoaderDefines.hxx"

#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDCoupling1GTUMesh.hxx"

#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"

#include <vector>

namespace ParaMEDMEM
{
  class DataArrayInt;
  class MEDCouplingMesh;
  class MEDFileMesh;
  class MEDFileUMesh;
  class MEDFileCMesh;
  class MEDFileCurveLinearMesh;
  class MEDFileFieldGlobs;
  class MEDFileFieldGlobsReal;
  class MEDFileAnyTypeField1TS;
  class MEDFileAnyTypeFieldMultiTS;

  class MEDFileMeshStruct : public RefCountObject
  {
  public:
    MEDLOADER_EXPORT static MEDFileMeshStruct *New(const MEDFileMesh *mesh);
    std::size_t getHeapMemorySize() const;
    const MEDFileMesh *getTheMesh() const { return _mesh; }
    int getNumberOfNodes() const { return _nb_nodes; }
    int getNumberOfElemsOfGeoType(INTERP_KERNEL::NormalizedCellType t) const throw(INTERP_KERNEL::Exception);
    int getLevelOfGeoType(INTERP_KERNEL::NormalizedCellType t) const throw(INTERP_KERNEL::Exception);
    int getNumberOfLevs() const;
    int getNumberOfGeoTypesInLev(int relativeLev) const throw(INTERP_KERNEL::Exception);
  private:
    MEDFileMeshStruct(const MEDFileMesh *mesh);
  private:
    const MEDFileMesh *_mesh;
    std::string _name;
    int _nb_nodes;
    std::vector< std::vector<int> > _geo_types_distrib;
  }; 
  
  class MEDFileField1TSStructItem;
  
  class MEDMeshMultiLev : public RefCountObject
  {
  public:
    std::size_t getHeapMemorySize() const;
  public:
    static MEDMeshMultiLev *New(const MEDFileMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities) throw(INTERP_KERNEL::Exception);
    static MEDMeshMultiLev *New(const MEDFileMesh *m, const std::vector<int>& levs) throw(INTERP_KERNEL::Exception);
    static MEDMeshMultiLev *NewOnlyOnNode(const MEDFileMesh *m, const DataArrayInt *pflOnNode) throw(INTERP_KERNEL::Exception);
    void setNodeReduction(const DataArrayInt *nr);
    bool isFastlyTheSameStruct(const MEDFileField1TSStructItem& fst, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    MEDLOADER_EXPORT DataArray *buildDataArray(const MEDFileField1TSStructItem& fst, const MEDFileFieldGlobsReal *globs, const DataArray *vals) const throw(INTERP_KERNEL::Exception);
    virtual void selectPartOfNodes(const DataArrayInt *pflNodes) throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDMeshMultiLev *prepare() const throw(INTERP_KERNEL::Exception) = 0;
    int getNumberOfCells(INTERP_KERNEL::NormalizedCellType t) const throw(INTERP_KERNEL::Exception);
    int getNumberOfNodes() const throw(INTERP_KERNEL::Exception);
  protected:
    std::string getPflNameOfId(int id) const;
    DataArray *constructDataArray(const MEDFileField1TSStructItem& fst, const MEDFileFieldGlobsReal *globs, const DataArray *vals) const throw(INTERP_KERNEL::Exception);
  protected:
    MEDMeshMultiLev();
    MEDMeshMultiLev(const MEDMeshMultiLev& other);
    MEDMeshMultiLev(int nbNodes, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities);
  protected:
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > _pfls;
    std::vector< INTERP_KERNEL::NormalizedCellType > _geo_types;
    std::vector<int> _nb_entities;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _node_reduction;
    int _nb_nodes;
  public:
    static const int PARAMEDMEM_2_VTKTYPE_LGTH=34;
    static const unsigned char PARAMEDMEM_2_VTKTYPE[PARAMEDMEM_2_VTKTYPE_LGTH];
  };
  
  class MEDStructuredMeshMultiLev;
  
  class MEDUMeshMultiLev : public MEDMeshMultiLev
  {
  public:
    static MEDUMeshMultiLev *New(const MEDFileUMesh *m, const std::vector<int>& levs) throw(INTERP_KERNEL::Exception);
    static MEDUMeshMultiLev *New(const MEDFileUMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities) throw(INTERP_KERNEL::Exception);
    void selectPartOfNodes(const DataArrayInt *pflNodes) throw(INTERP_KERNEL::Exception);
    MEDMeshMultiLev *prepare() const throw(INTERP_KERNEL::Exception);
    MEDUMeshMultiLev(const MEDStructuredMeshMultiLev& other, const MEDCouplingAutoRefCountObjectPtr<MEDCoupling1GTUMesh>& part);
    MEDLOADER_EXPORT void buildVTUArrays(DataArrayDouble *& coords, DataArrayByte *&types, DataArrayInt *&cellLocations, DataArrayInt *& cells, DataArrayInt *&faceLocations, DataArrayInt *&faces) const throw(INTERP_KERNEL::Exception);
  private:
    void reorderNodesIfNecessary(MEDCouplingAutoRefCountObjectPtr<DataArrayDouble>& coords, DataArrayInt *nodalConnVTK, DataArrayInt *polyhedNodalConnVTK) const throw(INTERP_KERNEL::Exception);
  private:
    MEDUMeshMultiLev(const MEDUMeshMultiLev& other);
    MEDUMeshMultiLev(const MEDFileUMesh *m, const std::vector<int>& levs);
    MEDUMeshMultiLev(const MEDFileUMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCoupling1GTUMesh> > _parts;
  };

  class MEDStructuredMeshMultiLev : public MEDMeshMultiLev
  {
  public:
    void selectPartOfNodes(const DataArrayInt *pflNodes) throw(INTERP_KERNEL::Exception);
    virtual std::vector<int> getNodeGridStructure() const throw(INTERP_KERNEL::Exception) = 0;
  protected:
    MEDStructuredMeshMultiLev();
    MEDStructuredMeshMultiLev(const MEDStructuredMeshMultiLev& other);
    MEDStructuredMeshMultiLev(int nbOfNodes, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities);
  };
  
  class MEDCMeshMultiLev : public MEDStructuredMeshMultiLev
  {
  public:
    static MEDCMeshMultiLev *New(const MEDFileCMesh *m, const std::vector<int>& levs) throw(INTERP_KERNEL::Exception);
    static MEDCMeshMultiLev *New(const MEDFileCMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities) throw(INTERP_KERNEL::Exception);
    std::vector<int> getNodeGridStructure() const throw(INTERP_KERNEL::Exception);
    MEDMeshMultiLev *prepare() const throw(INTERP_KERNEL::Exception);
    MEDLOADER_EXPORT std::vector< DataArrayDouble * > buildVTUArrays() const throw(INTERP_KERNEL::Exception);
  private:
    MEDCMeshMultiLev(const MEDCMeshMultiLev& other);
    MEDCMeshMultiLev(const MEDFileCMesh *m, const std::vector<int>& levs);
    MEDCMeshMultiLev(const MEDFileCMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > _coords;
  };

  class MEDCurveLinearMeshMultiLev : public MEDStructuredMeshMultiLev
  {
  public:
    static MEDCurveLinearMeshMultiLev *New(const MEDFileCurveLinearMesh *m, const std::vector<int>& levs) throw(INTERP_KERNEL::Exception);
    static MEDCurveLinearMeshMultiLev *New(const MEDFileCurveLinearMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls , const std::vector<int>& nbEntities) throw(INTERP_KERNEL::Exception);
    std::vector<int> getNodeGridStructure() const throw(INTERP_KERNEL::Exception);
    MEDMeshMultiLev *prepare() const throw(INTERP_KERNEL::Exception);
    MEDLOADER_EXPORT void buildVTUArrays(DataArrayDouble *&coords, std::vector<int>& nodeStrct) const throw(INTERP_KERNEL::Exception);
  private:
    MEDCurveLinearMeshMultiLev(const MEDCurveLinearMeshMultiLev& other);
    MEDCurveLinearMeshMultiLev(const MEDFileCurveLinearMesh *m, const std::vector<int>& levs);
    MEDCurveLinearMeshMultiLev(const MEDFileCurveLinearMesh *m, const std::vector<INTERP_KERNEL::NormalizedCellType>& gts, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& nbEntities);
  private:
    MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> _coords;
    std::vector<int> _structure;
  };

  class MEDFileField1TSStructItem2 : public RefCountObject
  {
  public:
    MEDFileField1TSStructItem2();
    MEDFileField1TSStructItem2(INTERP_KERNEL::NormalizedCellType a, const std::pair<int,int>& b, const std::string& pfl, const std::string& loc);
    void checkWithMeshStructForCells(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs) throw(INTERP_KERNEL::Exception);
    void checkWithMeshStructForGaussNE(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs) throw(INTERP_KERNEL::Exception);
    void checkWithMeshStructForGaussPT(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs) throw(INTERP_KERNEL::Exception);
    //
    MEDLOADER_EXPORT std::size_t getHeapMemorySize() const;
    //
    const DataArrayInt *getPfl(const MEDFileFieldGlobsReal *globs) const;
    INTERP_KERNEL::NormalizedCellType getGeo() const { return _geo_type; }
    int getNbEntity() const { return _nb_of_entity; }
    const std::pair<int,int>& getStartStop() const { return _start_end; }
    std::string getPflName() const;
    int getNbOfIntegrationPts(const MEDFileFieldGlobsReal *globs) const;
    //! warning this method also set _nb_of_entity attribute !
    void checkInRange(int nbOfEntity, int nip, const MEDFileFieldGlobsReal *globs) throw(INTERP_KERNEL::Exception);
    bool isFastlyEqual(int& startExp, INTERP_KERNEL::NormalizedCellType gt, const char *pflName) const;
    bool operator==(const MEDFileField1TSStructItem2& other) const throw(INTERP_KERNEL::Exception);
    bool isCellSupportEqual(const MEDFileField1TSStructItem2& other, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    bool isNodeSupportEqual(const MEDFileField1TSStructItem2& other, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    static MEDFileField1TSStructItem2 BuildAggregationOf(const std::vector<const MEDFileField1TSStructItem2 *>& objs, const MEDFileFieldGlobsReal *globs) throw(INTERP_KERNEL::Exception);
  public:
    static const char NEWLY_CREATED_PFL_NAME[];
  private:
    INTERP_KERNEL::NormalizedCellType _geo_type;
    std::pair<int,int> _start_end;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _pfl;
    std::string _loc;
    int _nb_of_entity;
  };

  class MEDFileField1TSStructItem : public RefCountObject
  {
  public:
    MEDFileField1TSStructItem() { }
    MEDFileField1TSStructItem(TypeOfField a, const std::vector< MEDFileField1TSStructItem2 >& b);
    void checkWithMeshStruct(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs) throw(INTERP_KERNEL::Exception);
    bool operator==(const MEDFileField1TSStructItem& other) const throw(INTERP_KERNEL::Exception);
    MEDLOADER_EXPORT std::size_t getHeapMemorySize() const;
    bool isEntityCell() const;
    bool isComputed() const { return _computed; }
    TypeOfField getType() const { return _type; }
    std::size_t getNumberOfItems() const { return _items.size(); }
    const MEDFileField1TSStructItem2& operator[](std::size_t i) const throw(INTERP_KERNEL::Exception);
    //
    bool isCellSupportEqual(const MEDFileField1TSStructItem& other, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    bool isNodeSupportEqual(const MEDFileField1TSStructItem& other, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    MEDFileField1TSStructItem simplifyMeOnCellEntity(const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    bool isCompatibleWithNodesDiscr(const MEDFileField1TSStructItem& other, const MEDFileMeshStruct *meshSt, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    bool isFullyOnOneLev(const MEDFileMeshStruct *meshSt, int& theFirstLevFull) const throw(INTERP_KERNEL::Exception);
    MEDLOADER_EXPORT MEDMeshMultiLev *buildFromScratchDataSetSupportOnCells(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    MEDLOADER_EXPORT static MEDFileField1TSStructItem BuildItemFrom(const MEDFileAnyTypeField1TS *ref, const MEDFileMeshStruct *meshSt);
  private:
    bool _computed;
    TypeOfField _type;
    std::vector< MEDFileField1TSStructItem2 > _items;
  };

  class MEDFileField1TSStruct : public RefCountObject
  {
  public:
    static MEDFileField1TSStruct *New(const MEDFileAnyTypeField1TS *ref, MEDFileMeshStruct *mst) throw(INTERP_KERNEL::Exception);
    void checkWithMeshStruct(MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    bool isEqualConsideringThePast(const MEDFileAnyTypeField1TS *other, const MEDFileMeshStruct *mst) const throw(INTERP_KERNEL::Exception);
    bool isSupportSameAs(const MEDFileAnyTypeField1TS *other, const MEDFileMeshStruct *meshSt) throw(INTERP_KERNEL::Exception);
    bool isCompatibleWithNodesDiscr(const MEDFileAnyTypeField1TS *other, const MEDFileMeshStruct *meshSt) throw(INTERP_KERNEL::Exception);
    MEDLOADER_EXPORT MEDMeshMultiLev *buildFromScratchDataSetSupport(const MEDFileMeshStruct *mst, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    bool isDataSetSupportFastlyEqualTo(const MEDFileField1TSStruct& other, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
  private:
    MEDFileField1TSStruct(const MEDFileAnyTypeField1TS *ref, MEDFileMeshStruct *mst);
    bool presenceOfCellDiscr(int& pos) const throw(INTERP_KERNEL::Exception);
    bool presenceOfPartialNodeDiscr(int& pos) const throw(INTERP_KERNEL::Exception);
  private:
    std::vector<MEDFileField1TSStructItem> _already_checked;
  };

  class MEDFileFastCellSupportComparator : public RefCountObject
  {
  public:
    MEDLOADER_EXPORT static MEDFileFastCellSupportComparator *New(const MEDFileMeshStruct *m, const MEDFileAnyTypeFieldMultiTS *ref) throw(INTERP_KERNEL::Exception);
    MEDLOADER_EXPORT MEDMeshMultiLev *buildFromScratchDataSetSupport(int timeStepId, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    MEDLOADER_EXPORT bool isDataSetSupportEqualToThePreviousOne(int timeStepId, const MEDFileFieldGlobsReal *globs) const throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDFileAnyTypeFieldMultiTS *other) throw(INTERP_KERNEL::Exception);
    bool isCompatibleWithNodesDiscr(const MEDFileAnyTypeFieldMultiTS *other) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
  private:
    MEDFileFastCellSupportComparator(const MEDFileMeshStruct *m, const MEDFileAnyTypeFieldMultiTS *ref);
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDFileMeshStruct> _mesh_comp;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSStruct> > _f1ts_cmps;
  };
}

#endif
