// Copyright (C) 2007-2025  CEA, EDF
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
// Author : Anthony Geay (EDF R&D)

#ifndef __MEDFILEFIELDINTERNAL_HXX__
#define __MEDFILEFIELDINTERNAL_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDFileUtilities.hxx"
#include "NormalizedGeometricTypes"
#include "InterpKernelAutoPtr.hxx"
#include "MCAuto.hxx"

#include "med.h"

#include <string>
#include <list>

namespace MEDCoupling
{
  class MEDFileGTKeeper
  {
  public:
    virtual MEDFileGTKeeper *deepCopy() const = 0;
    virtual INTERP_KERNEL::NormalizedCellType getGeoType() const = 0;
    virtual std::string getRepr() const = 0;
    virtual bool isEqual(const MEDFileGTKeeper *other) const = 0;
    virtual ~MEDFileGTKeeper();
  };

  class MEDFileGTKeeperSta : public MEDFileGTKeeper
  {
  public:
    MEDFileGTKeeperSta(INTERP_KERNEL::NormalizedCellType gt):_geo_type(gt) { }
    MEDFileGTKeeper *deepCopy() const;
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    std::string getRepr() const;
    bool isEqual(const MEDFileGTKeeper *other) const;
  private:
    INTERP_KERNEL::NormalizedCellType _geo_type;
  };


  class MEDFileStructureElement;
  class MEDFileUMesh;

  class MEDFileGTKeeperDyn : public MEDFileGTKeeper
  {
  public:
    MEDFileGTKeeperDyn(const MEDFileUMesh *mesh, const MEDFileUMesh *section, const MEDFileStructureElement *se);
    MEDFileGTKeeper *deepCopy() const;
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    std::string getRepr() const;
    bool isEqual(const MEDFileGTKeeper *other) const;
    const MEDFileUMesh *getMesh() const { return _mesh; }
    const MEDFileUMesh *getSection() const { return _section; }
    const MEDFileStructureElement *getSE() const { return _se; }
  private:
    MCConstAuto<MEDFileUMesh> _mesh;
    MCConstAuto<MEDFileUMesh> _section;
    MCConstAuto<MEDFileStructureElement> _se;
  };

  class MEDFileEntities;

  class MEDFileFieldLoc : public RefCountObject
  {
  public:
    MEDLOADER_EXPORT void simpleRepr(std::ostream& oss) const;
    MEDLOADER_EXPORT std::string getName() const { return _name; }
    MEDLOADER_EXPORT void setName(const std::string& name);
    std::string getClassName() const override { return std::string("MEDFileFieldLoc"); }
    static MEDFileFieldLoc *New(med_idt fid, const std::string& locName);
    static MEDFileFieldLoc *New(med_idt fid, int i, const MEDFileEntities *entities);
    static MEDFileFieldLoc *New(const std::string& locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileFieldLoc *deepCopy() const;
    bool isOnStructureElement() const;
    const MEDFileGTKeeper *getUndergroundGTKeeper() const { return _gt; }
    MEDLOADER_EXPORT int getNbOfGaussPtPerCell() const { return _nb_gauss_pt; }
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    MEDLOADER_EXPORT std::string repr() const;
    MEDLOADER_EXPORT bool isName(const std::string& name) const { return _name==name; }
    MEDLOADER_EXPORT int getDimension() const { return _dim; }
    MEDLOADER_EXPORT int getNumberOfGaussPoints() const { return _nb_gauss_pt; }
    MEDLOADER_EXPORT int getNumberOfPointsInCells() const { return _nb_node_per_cell; }
    MEDLOADER_EXPORT const std::vector<double>& getRefCoords() const { return _ref_coo; }
    MEDLOADER_EXPORT const std::vector<double>& getGaussCoords() const { return _gs_coo; }
    MEDLOADER_EXPORT const std::vector<double>& getGaussWeights() const { return _w; }
    MEDLOADER_EXPORT INTERP_KERNEL::NormalizedCellType getGeoType() const { return _gt->getGeoType(); }
    MEDLOADER_EXPORT bool isEqual(const MEDFileFieldLoc& other, double eps) const;
  private:
    MEDFileFieldLoc(const MEDFileFieldLoc& other);
    MEDFileFieldLoc(med_idt fid, const std::string& locName);
    MEDFileFieldLoc(med_idt fid, int id, const MEDFileEntities *entities);
    MEDFileFieldLoc(const std::string& locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w);
  private:
    int _dim;
    int _nb_gauss_pt;
    INTERP_KERNEL::AutoCppPtr<MEDFileGTKeeper> _gt;
    int _nb_node_per_cell;
    std::string _name;
    std::vector<double> _ref_coo;
    std::vector<double> _gs_coo;
    std::vector<double> _w;
  };

  /// @cond INTERNAL
  class MEDFileFieldPerMeshPerTypeCommon;
  class MEDFileFieldPerMeshPerType;
  class MEDCouplingFieldTemplate;
  class MEDFileFieldNameScope;
  class MEDFileFieldGlobsReal;
  class MEDCouplingMesh;

  class MEDFileFieldPerMeshPerTypePerDisc : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMeshPerTypePerDisc *NewOnRead(MEDFileFieldPerMeshPerTypeCommon *fath, TypeOfField type, int profileIt, const PartDefinition *pd);
    static MEDFileFieldPerMeshPerTypePerDisc *New(MEDFileFieldPerMeshPerTypeCommon *fath, TypeOfField type, mcIdType locId);
    static MEDFileFieldPerMeshPerTypePerDisc *New(const MEDFileFieldPerMeshPerTypePerDisc& other);
    std::string getClassName() const override { return std::string("MEDFileFieldPerMeshPerTypePerDisc"); }
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileFieldPerMeshPerTypePerDisc *deepCopy(MEDFileFieldPerMeshPerTypeCommon *father) const;
    void assignFieldNoProfile(mcIdType& start, mcIdType offset, mcIdType nbOfCells, const MEDCouplingFieldTemplate *field, const DataArray *arrr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void assignFieldProfile(bool isPflAlone, mcIdType& start, const DataArrayIdType *multiTypePfl, const DataArrayIdType *idsInPfl, DataArrayIdType *locIds, mcIdType nbOfEltsInWholeMesh, const MEDCouplingFieldTemplate *field, const DataArray *arrr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void assignNodeFieldNoProfile(mcIdType& start, const MEDCouplingFieldTemplate *field, const DataArray *arrr, MEDFileFieldGlobsReal& glob);
    void getCoarseData(TypeOfField& type, std::pair<mcIdType,mcIdType>& dad, std::string& pfl, std::string& loc) const;
    void writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const;
    const MEDFileFieldPerMeshPerTypeCommon *getFather() const;
    void loadOnlyStructureOfDataRecursively(med_idt fid, mcIdType& start, const MEDFileFieldNameScope& nasc);
    void loadBigArray(med_idt fid, const MEDFileFieldNameScope& nasc);
    void setNewStart(mcIdType newValueOfStart);
    int getIteration() const;
    int getOrder() const;
    double getTime() const;
    std::string getMeshName() const;
    TypeOfField getType() const;
    void simpleRepr(int bkOffset, std::ostream& oss, int id) const;
    void fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const;
    void setType(TypeOfField newType);
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    INTERP_KERNEL::NormalizedCellType getGeoTypeStatic() const;
    std::size_t getNumberOfComponents() const;
    mcIdType getNumberOfTuples() const;
    mcIdType getStart() const { return _start; }
    mcIdType getEnd() const { return _end; }
    void setEnd(mcIdType endd) { _end=endd; }
    mcIdType getNumberOfVals() const { return _nval; }
    void incrementNbOfVals(mcIdType deltaNbVal);
    DataArray *getOrCreateAndGetArray();
    const DataArray *getOrCreateAndGetArray() const;
    const std::vector<std::string>& getInfo() const;
    std::string getProfile() const;
    void setProfile(const std::string& newPflName);
    std::string getLocalization() const;
    void setLocalization(const std::string& newLocName);
    mcIdType getLocId() const { return _loc_id; }
    void setLocId(mcIdType newId) const { _loc_id=newId; }
    void setFather(MEDFileFieldPerMeshPerTypeCommon *newFather) { _father=newFather; }
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    void getFieldAtLevel(TypeOfField type, const MEDFileFieldGlobsReal *glob, std::vector< std::pair<mcIdType,mcIdType> >& dads, std::vector<const DataArrayIdType *>& pfls, std::vector<int>& locs,
                         std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const;
    void fillValues(int discId, mcIdType& startEntryId, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<mcIdType,mcIdType> > >& entries) const;
    mcIdType fillEltIdsFromCode(mcIdType offset, const std::vector<mcIdType>& codeOfMesh, const MEDFileFieldGlobsReal& glob, mcIdType *ptToFill) const;
    mcIdType fillTupleIds(mcIdType *ptToFill) const;
    static int ConvertType(TypeOfField type, mcIdType locId);
    static std::vector< std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> > SplitPerDiscretization(const std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>& entries);
    static bool RenumberChunks(mcIdType offset, const std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>& entriesOnSameDisc,
                               const DataArrayIdType *explicitIdsInMesh, const std::vector<mcIdType>& newCode,
                               MEDFileFieldGlobsReal& glob, DataArrayDouble *arr, std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> >& result);
    static MEDFileFieldPerMeshPerTypePerDisc *NewObjectOnSameDiscThanPool(TypeOfField typeF, INTERP_KERNEL::NormalizedCellType geoType, DataArrayIdType *idsOfMeshElt,
                                                                          bool isPfl, mcIdType nbi, mcIdType offset, std::list< const MEDFileFieldPerMeshPerTypePerDisc *>& entriesOnSameDisc,
                                                                          MEDFileFieldGlobsReal& glob, bool &notInExisting);
    static MCAuto<MEDFileFieldPerMeshPerTypePerDisc> Aggregate(mcIdType &start, const std::vector<std::pair<int,const MEDFileFieldPerMeshPerTypePerDisc *> >& pms, const std::vector< std::vector< std::pair<int,mcIdType> > >& dts, TypeOfField tof, MEDFileFieldPerMeshPerType *father, std::vector<std::pair< int, std::pair<mcIdType,mcIdType> > >& extractInfo);
    MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerTypeCommon *fath, TypeOfField type):_type(type),_father(fath),_start(-1),_end(-1),_nval(-1),_loc_id(-5),_profile_it(-1) { }
  private:
    MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerTypeCommon *fath, TypeOfField type, mcIdType profileIt, const PartDefinition *pd);
    MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerTypeCommon *fath, TypeOfField type, mcIdType profileIt, const std::string& dummy);
    MEDFileFieldPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc& other);
    MEDFileFieldPerMeshPerTypePerDisc();
  private:
    void goReadZeValuesInFile(med_idt fid, const std::string& fieldName, int nbOfCompo, int iteration, int order, med_entity_type menti, med_geometry_type mgeoti, unsigned char *startFeedingPtr);
  private:
    TypeOfField _type;
    MEDFileFieldPerMeshPerTypeCommon *_father;
    mcIdType _start;
    mcIdType _end;
    //! _nval is different than end-start in case of ON_GAUSS_PT and ON_GAUSS_NE ! (_nval=(_end-_start)/nbi)
    mcIdType _nval;
    std::string _profile;
    std::string _localization;
    //! only on assignment -3 : ON_NODES, -2 : ON_CELLS, -1 : ON_GAUSS_NE, 0..* : ON_GAUSS_PT
    mutable mcIdType _loc_id;
    mutable mcIdType _profile_it;
    MCAuto<PartDefinition> _pd;
  public:
    mutable mcIdType _tmp_work1;
  };

  class MEDFileFieldVisitor;
  class MEDFileFieldPerMesh;

  class MEDFileFieldPerMeshPerTypeCommon : public RefCountObject, public MEDFileWritable
  {
  public:
    std::string getClassName() const override { return std::string("MEDFileFieldPerMeshPerTypeCommon"); }
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    void assignFieldNoProfile(mcIdType& start, mcIdType offset, mcIdType nbOfCells, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void assignFieldProfile(bool isPflAlone, mcIdType& start, const DataArrayIdType *multiTypePfl, const DataArrayIdType *idsInPfl, DataArrayIdType *locIds, mcIdType nbOfEltsInWholeMesh, const MEDCouplingFieldTemplate *field, const DataArray *arr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void assignNodeFieldNoProfile(mcIdType& start, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob);
    void assignNodeFieldProfile(mcIdType& start, const DataArrayIdType *pfl, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    const MEDFileFieldPerMesh *getFather() const;
    void loadOnlyStructureOfDataRecursively(med_idt fid, mcIdType &start, const MEDFileFieldNameScope& nasc);
    void loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc);
    void writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const;
    bool isUniqueLevel(int& dim) const;
    void fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const;
    void fillFieldSplitedByType(std::vector< std::pair<mcIdType,mcIdType> >& dads, std::vector<TypeOfField>& types, std::vector<std::string>& pfls, std::vector<std::string>& locs) const;
    int getIteration() const;
    int getOrder() const;
    double getTime() const;
    std::string getMeshName() const;
    void getSizes(mcIdType& globalSz, mcIdType& nbOfEntries) const;
    std::size_t getNumberOfComponents() const;
    bool presenceOfMultiDiscPerGeoType() const;
    void pushDiscretization(MEDFileFieldPerMeshPerTypePerDisc *disc);
    DataArray *getOrCreateAndGetArray();
    const DataArray *getOrCreateAndGetArray() const;
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenLocId(mcIdType locId);
    const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenLocId(mcIdType locId) const;
    int getNumberOfLoc() const { return (int)_field_pm_pt_pd.size(); }
    int locIdOfLeaf(const MEDFileFieldPerMeshPerTypePerDisc *leaf) const;
    void fillValues(mcIdType& startEntryId, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<mcIdType,mcIdType> > >& entries) const;
    void setLeaves(const std::vector< MCAuto< MEDFileFieldPerMeshPerTypePerDisc > >& leaves);
    bool keepOnlySpatialDiscretization(TypeOfField tof, mcIdType &globalNum, std::vector< std::pair<mcIdType,mcIdType> >& its);
    bool keepOnlyGaussDiscretization(std::size_t idOfDisc, mcIdType &globalNum, std::vector< std::pair<mcIdType,mcIdType> >& its);
    static med_entity_type ConvertIntoMEDFileType(TypeOfField ikType, INTERP_KERNEL::NormalizedCellType ikGeoType, med_geometry_type& medfGeoType);
    MEDFileFieldPerMeshPerTypeCommon(MEDFileFieldPerMesh *father):_father(father) { }
    void setFather(MEDFileFieldPerMesh *father);
    void accept(MEDFileFieldVisitor& visitor) const;
  public:
    virtual ~MEDFileFieldPerMeshPerTypeCommon();
    virtual void getDimension(int& dim) const = 0;
    virtual INTERP_KERNEL::NormalizedCellType getGeoType() const = 0;
    virtual INTERP_KERNEL::NormalizedCellType getGeoTypeStatic() const = 0;
    virtual void entriesForMEDfile(TypeOfField mct, med_geometry_type& gt, med_entity_type& ent) const = 0;
    virtual void simpleRepr(int bkOffset, std::ostream& oss, int id) const = 0;
    virtual std::string getGeoTypeRepr() const = 0;
    virtual MEDFileFieldPerMeshPerTypeCommon *deepCopy(MEDFileFieldPerMesh *father) const = 0;
    virtual void getFieldAtLevel(int meshDim, TypeOfField type, const MEDFileFieldGlobsReal *glob, std::vector< std::pair<mcIdType,mcIdType> >& dads, std::vector<const DataArrayIdType *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const = 0;
  protected:
    void deepCopyElements();
    std::vector<mcIdType> addNewEntryIfNecessary(const MEDCouplingFieldTemplate *field, mcIdType offset, mcIdType nbOfCells);
    std::vector<mcIdType> addNewEntryIfNecessaryGauss(const MEDCouplingFieldTemplate *field, mcIdType offset, mcIdType nbOfCells);
    std::vector<mcIdType> addNewEntryIfNecessary(const MEDCouplingFieldTemplate *field, const DataArrayIdType *subCells);
    std::vector<mcIdType> addNewEntryIfNecessaryGauss(const MEDCouplingFieldTemplate *field, const DataArrayIdType *subCells);
  private:
    MEDFileFieldPerMesh *_father;
  protected:
    std::vector< MCAuto<MEDFileFieldPerMeshPerTypePerDisc> > _field_pm_pt_pd;
  };

  class MEDFileFieldPerMeshPerType : public MEDFileFieldPerMeshPerTypeCommon
  {
  public:
    std::string getClassName() const override { return std::string("MEDFileFieldPerMeshPerType"); }
    static MEDFileFieldPerMeshPerType *New(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType);
    static MEDFileFieldPerMeshPerType *NewOnRead(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, const MEDFileFieldNameScope& nasc, const PartDefinition *pd);
    static MCAuto<MEDFileFieldPerMeshPerType> Aggregate(mcIdType &start, const std::vector< std::pair<int,const MEDFileFieldPerMeshPerType *> >& pms, const std::vector< std::vector< std::pair<int,mcIdType> > >& dts, INTERP_KERNEL::NormalizedCellType gt, MEDFileFieldPerMesh *father, std::vector<std::pair< int, std::pair<mcIdType,mcIdType> > >& extractInfo);
  public:// overload of abstract methods
    void getDimension(int& dim) const;
    INTERP_KERNEL::NormalizedCellType getGeoType() const override;
    INTERP_KERNEL::NormalizedCellType getGeoTypeStatic() const override;
    void entriesForMEDfile(TypeOfField mct, med_geometry_type& gt, med_entity_type& ent) const;
    void simpleRepr(int bkOffset, std::ostream& oss, int id) const;
    std::string getGeoTypeRepr() const;
    MEDFileFieldPerMeshPerType *deepCopy(MEDFileFieldPerMesh *father) const;
    void getFieldAtLevel(int meshDim, TypeOfField type, const MEDFileFieldGlobsReal *glob, std::vector< std::pair<mcIdType,mcIdType> >& dads, std::vector<const DataArrayIdType *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const;
  private:
    MEDFileFieldPerMeshPerType(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, const MEDFileFieldNameScope& nasc, const PartDefinition *pd);
    MEDFileFieldPerMeshPerType(MEDFileFieldPerMesh *father, INTERP_KERNEL::NormalizedCellType gt);
  private:
    INTERP_KERNEL::NormalizedCellType _geo_type;
  };

  class MEDFileFieldPerMeshPerTypeDyn : public MEDFileFieldPerMeshPerTypeCommon
  {
  public:
    static MEDFileFieldPerMeshPerTypeDyn *NewOnRead(med_idt fid, MEDFileFieldPerMesh *fath, const MEDFileEntities *entities, int idGT, const MEDFileFieldNameScope& nasc);
    std::string getClassName() const override { return std::string("MEDFileFieldPerMeshPerTypeDyn"); }
    int getDynGT() const;
    std::string getModelName() const;
  public:
    void getDimension(int& dim) const;
    INTERP_KERNEL::NormalizedCellType getGeoType() const override;
    INTERP_KERNEL::NormalizedCellType getGeoTypeStatic() const override;
    void entriesForMEDfile(TypeOfField mct, med_geometry_type& gt, med_entity_type& ent) const;
    void simpleRepr(int bkOffset, std::ostream& oss, int id) const;
    std::string getGeoTypeRepr() const;
    MEDFileFieldPerMeshPerTypeDyn *deepCopy(MEDFileFieldPerMesh *father) const;
    void getFieldAtLevel(int meshDim, TypeOfField type, const MEDFileFieldGlobsReal *glob, std::vector< std::pair<mcIdType,mcIdType> >& dads, std::vector<const DataArrayIdType *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const;
  private:
    MEDFileFieldPerMeshPerTypeDyn(med_idt fid, MEDFileFieldPerMesh *fath, const MEDFileStructureElement *se, const MEDFileFieldNameScope& nasc);
  private:
    MCConstAuto<MEDFileStructureElement> _se;
  };

  class MEDFileMesh;
  class MEDFileAnyTypeField1TSWithoutSDA;
  class MEDFileField1TSWithoutSDA;

  class MEDFileFieldPerMesh : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMesh *New(MEDFileAnyTypeField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh);
    static MEDFileFieldPerMesh *NewOnRead(med_idt fid, MEDFileAnyTypeField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder, const MEDFileFieldNameScope& nasc, const MEDFileMesh *mm, const MEDFileEntities *entities);
    static MEDFileFieldPerMesh *NewOnRead(med_idt fid, MEDFileAnyTypeField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder, const MEDFileFieldNameScope& nasc, const PartDefinition *pd, const MEDFileEntities *entities);
    std::string getClassName() const override { return std::string("MEDFileFieldPerMesh"); }
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileFieldPerMesh *deepCopy(MEDFileAnyTypeField1TSWithoutSDA *father) const;
    void simpleRepr(int bkOffset,std::ostream& oss, int id) const;
    void copyTinyInfoFrom(const MEDCouplingMesh *mesh);
    void assignFieldProfile(mcIdType& start, const DataArrayIdType *multiTypePfl, const std::vector<mcIdType>& code, const std::vector<mcIdType>& code2, const std::vector<DataArrayIdType *>& idsInPflPerType, const std::vector<DataArrayIdType *>& idsPerType, const MEDCouplingFieldTemplate *field, const DataArray *arr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void assignFieldNoProfileNoRenum(mcIdType& start, const std::vector<mcIdType>& code, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void assignNodeFieldNoProfile(mcIdType& start, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob);
    void assignNodeFieldProfile(mcIdType& start, const DataArrayIdType *pfl, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void loadOnlyStructureOfDataRecursively(med_idt fid, mcIdType &start, const MEDFileFieldNameScope& nasc);
    void loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc);
    void writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const;
    void fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const;
    std::vector< std::vector< std::pair<mcIdType,mcIdType> > > getFieldSplitedByType(std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    void accept(MEDFileFieldVisitor& visitor) const;
    void getDimension(int& dim) const;
    bool isUniqueLevel(int& dim) const;
    double getTime() const;
    int getIteration() const;
    int getOrder() const;
    int getMeshIteration() const { return _mesh_iteration; }
    int getMeshOrder() const { return _mesh_order; }
    std::string getMeshName() const;
    void setMeshName(const std::string& meshName);
    std::size_t getNumberOfComponents() const;
    bool presenceOfMultiDiscPerGeoType() const;
    bool presenceOfStructureElements() const;
    bool onlyStructureElements() const;
    void killStructureElements();
    void keepOnlyStructureElements();
    void keepOnlyOnSE(const std::string& seName);
    void getMeshSENames(std::vector< std::pair<std::string,std::string> >& ps) const;
    DataArray *getOrCreateAndGetArray();
    const DataArray *getOrCreateAndGetArray() const;
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    void convertMedBallIntoClassic();
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    bool renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<mcIdType>& oldCode, const std::vector<mcIdType>& newCode, const DataArrayIdType *renumO2N, MEDFileFieldGlobsReal& glob);
    void keepOnlySpatialDiscretization(TypeOfField tof, mcIdType &globalNum, std::vector< std::pair<mcIdType,mcIdType> >& its);
    void keepOnlyGaussDiscretization(std::size_t idOfDisc, mcIdType &globalNum, std::vector< std::pair<mcIdType,mcIdType> >& its);
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, bool& isPfl, MCAuto<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    DataArray *getFieldOnMeshAtLevelWithPfl(TypeOfField type, const MEDCouplingMesh *mesh, DataArrayIdType *&pfl, const MEDFileFieldGlobsReal *glob, const MEDFileFieldNameScope& nasc) const;
    void getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<mcIdType,mcIdType> > >& entries) const;
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, mcIdType locId);
    const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, mcIdType locId) const;
    static MCAuto<MEDFileFieldPerMesh> Aggregate(mcIdType &start, const std::vector<const MEDFileFieldPerMesh *>& pms, const std::vector< std::vector< std::pair<int,mcIdType> > >& dts, MEDFileAnyTypeField1TSWithoutSDA *father, std::vector<std::pair< int, std::pair<mcIdType,mcIdType> > >& extractInfo);
  private:
    int addNewEntryIfNecessary(INTERP_KERNEL::NormalizedCellType type);
    MEDCouplingFieldDouble *finishField(TypeOfField type, const MEDFileFieldGlobsReal *glob,
                                        const std::vector< std::pair<mcIdType,mcIdType> >& dads, const std::vector<int>& locs, const MEDCouplingMesh *mesh, bool& isPfl, MCAuto<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    MEDCouplingFieldDouble *finishField2(TypeOfField type, const MEDFileFieldGlobsReal *glob,
                                         const std::vector< std::pair<mcIdType,mcIdType> >& dads, const std::vector<int>& locs,
                                         const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes,
                                         const MEDCouplingMesh *mesh, const DataArrayIdType *da, bool& isPfl, MCAuto<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    MEDCouplingFieldDouble *finishFieldNode2(const MEDFileFieldGlobsReal *glob,
                                             const std::vector< std::pair<mcIdType,mcIdType> >& dads, const std::vector<int>& locs,
                                             const MEDCouplingMesh *mesh, const DataArrayIdType *da, bool& isPfl, MCAuto<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    DataArray *finishField4(const std::vector< std::pair<mcIdType,mcIdType> >& dads, const DataArrayIdType *pflIn, mcIdType nbOfElems, DataArrayIdType *&pflOut) const;
    void assignNewLeaves(const std::vector< MCAuto< MEDFileFieldPerMeshPerTypePerDisc > >& leaves);
    static void SortArraysPerType(const MEDFileFieldGlobsReal *glob, TypeOfField type,
                                  const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<mcIdType,mcIdType> >& dads, const std::vector<const DataArrayIdType *>& pfls, const std::vector<int>& locs,
                                  std::vector<mcIdType>& code, std::vector<DataArrayIdType *>& notNullPfls);
    static mcIdType ComputeNbOfElems(const MEDFileFieldGlobsReal *glob, TypeOfField type, const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<mcIdType,mcIdType> >& dads, const std::vector<int>& locs);
    MEDFileFieldPerMesh(med_idt fid, MEDFileAnyTypeField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder, const MEDFileFieldNameScope& nasc, const MEDFileMesh *mm, const MEDFileEntities *entities);
    MEDFileFieldPerMesh(med_idt fid, MEDFileAnyTypeField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder, const MEDFileFieldNameScope& nasc, const PartDefinition *pd, const MEDFileEntities *entities);
    MEDFileFieldPerMesh(MEDFileAnyTypeField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh);
    MEDFileFieldPerMesh(MEDFileAnyTypeField1TSWithoutSDA *fath, int meshIt, int meshOrd):_mesh_iteration(meshIt),_mesh_order(meshOrd),_father(fath) { }
  private:
    int _mesh_iteration;
    int _mesh_order;
    MEDFileAnyTypeField1TSWithoutSDA *_father;
    std::vector< MCAuto< MEDFileFieldPerMeshPerTypeCommon > > _field_pm_pt;
  };
}

#endif
