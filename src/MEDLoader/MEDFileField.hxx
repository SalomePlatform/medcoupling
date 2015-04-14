// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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

#ifndef __MEDFILEFIELD_HXX__
#define __MEDFILEFIELD_HXX__

#include "MEDLoaderDefines.hxx"

#include "MEDFileFieldOverView.hxx"
#include "MEDFileUtilities.hxx"

#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingMemArray.hxx"

#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"

#include <vector>
#include <string>
#include <list>
#include <set>

#include "med.h"

namespace ParaMEDMEM
{
  class MEDFileFieldGlobs;
  class MEDCouplingMesh;
  class MEDCouplingFieldDouble;
  class MEDFileMesh;

  class MEDFileFieldLoc : public RefCountObject
  {
  public:
    MEDLOADER_EXPORT void simpleRepr(std::ostream& oss) const;
    MEDLOADER_EXPORT std::string getName() const { return _name; }
    MEDLOADER_EXPORT void setName(const std::string& name);
    static MEDFileFieldLoc *New(med_idt fid, const std::string& locName);
    static MEDFileFieldLoc *New(med_idt fid, int id);
    static MEDFileFieldLoc *New(const std::string& locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileFieldLoc *deepCpy() const;
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
    MEDLOADER_EXPORT INTERP_KERNEL::NormalizedCellType getGeoType() const { return _geo_type; }
    MEDLOADER_EXPORT bool isEqual(const MEDFileFieldLoc& other, double eps) const;
  private:
    MEDFileFieldLoc(med_idt fid, const std::string& locName);
    MEDFileFieldLoc(med_idt fid, int id);
    MEDFileFieldLoc(const std::string& locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w);
  private:
    int _dim;
    int _nb_gauss_pt;
    int _nb_node_per_cell;
    std::string _name;
    INTERP_KERNEL::NormalizedCellType _geo_type;
    std::vector<double> _ref_coo;
    std::vector<double> _gs_coo;
    std::vector<double> _w;
  };

  /// @cond INTERNAL
  class MEDFileAnyTypeField1TSWithoutSDA;
  class MEDFileFieldPerMeshPerType;
  class MEDFileField1TSWithoutSDA;
  class MEDFileFieldNameScope;
  class MEDFileFieldGlobsReal;
  class MEDFileFieldPerMesh;
  class PartDefinition;

  class MEDFileFieldPerMeshPerTypePerDisc : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMeshPerTypePerDisc *NewOnRead(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int profileIt, const PartDefinition *pd);
    static MEDFileFieldPerMeshPerTypePerDisc *New(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int locId);
    static MEDFileFieldPerMeshPerTypePerDisc *New(const MEDFileFieldPerMeshPerTypePerDisc& other);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileFieldPerMeshPerTypePerDisc *deepCpy(MEDFileFieldPerMeshPerType *father) const;
    void assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, const DataArray *arrr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void assignFieldProfile(bool isPflAlone, int& start, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, DataArrayInt *locIds, int nbOfEltsInWholeMesh, const MEDCouplingFieldDouble *field, const DataArray *arrr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, const DataArray *arrr, MEDFileFieldGlobsReal& glob);
    void getCoarseData(TypeOfField& type, std::pair<int,int>& dad, std::string& pfl, std::string& loc) const;
    void writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const;
    const MEDFileFieldPerMeshPerType *getFather() const;
    void loadOnlyStructureOfDataRecursively(med_idt fid, int& start, const MEDFileFieldNameScope& nasc);
    void loadBigArray(med_idt fid, const MEDFileFieldNameScope& nasc);
    void setNewStart(int newValueOfStart);
    int getIteration() const;
    int getOrder() const;
    double getTime() const;
    std::string getMeshName() const;
    TypeOfField getType() const;
    void simpleRepr(int bkOffset, std::ostream& oss, int id) const;
    void fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const;
    void setType(TypeOfField newType);
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    int getNumberOfComponents() const;
    int getNumberOfTuples() const;
    int getStart() const { return _start; }
    int getEnd() const { return _end; }
    DataArray *getOrCreateAndGetArray();
    const DataArray *getOrCreateAndGetArray() const;
    const std::vector<std::string>& getInfo() const;
    std::string getProfile() const;
    void setProfile(const std::string& newPflName);
    std::string getLocalization() const;
    void setLocalization(const std::string& newLocName);
    int getLocId() const { return _loc_id; }
    void setLocId(int newId) const { _loc_id=newId; }
    void setFather(MEDFileFieldPerMeshPerType *newFather) { _father=newFather; }
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    void getFieldAtLevel(TypeOfField type, const MEDFileFieldGlobsReal *glob, std::vector< std::pair<int,int> >& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs,
                         std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const;
    void fillValues(int discId, int& startEntryId, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    int fillEltIdsFromCode(int offset, const std::vector<int>& codeOfMesh, const MEDFileFieldGlobsReal& glob, int *ptToFill) const;
    int fillTupleIds(int *ptToFill) const;
    static int ConvertType(TypeOfField type, int locId);
    static std::vector< std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> > SplitPerDiscretization(const std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>& entries);
    static bool RenumberChunks(int offset, const std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>& entriesOnSameDisc,
                               const DataArrayInt *explicitIdsInMesh, const std::vector<int>& newCode,
                               MEDFileFieldGlobsReal& glob, DataArrayDouble *arr, std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >& result);
    static MEDFileFieldPerMeshPerTypePerDisc *NewObjectOnSameDiscThanPool(TypeOfField typeF, INTERP_KERNEL::NormalizedCellType geoType, DataArrayInt *idsOfMeshElt,
                                                                          bool isPfl, int nbi, int offset, std::list< const MEDFileFieldPerMeshPerTypePerDisc *>& entriesOnSameDisc,
                                                                          MEDFileFieldGlobsReal& glob, bool &notInExisting);
  private:
    MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int profileIt, const PartDefinition *pd);
    MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int profileIt, const std::string& dummy);
    MEDFileFieldPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc& other);
    MEDFileFieldPerMeshPerTypePerDisc();
  private:
    void goReadZeValuesInFile(med_idt fid, const std::string& fieldName, int nbOfCompo, int iteration, int order, med_entity_type menti, med_geometry_type mgeoti, unsigned char *startFeedingPtr);
  private:
    TypeOfField _type;
    MEDFileFieldPerMeshPerType *_father;
    int _start;
    int _end;
    //! _nval is different than end-start in case of ON_GAUSS_PT and ON_GAUSS_NE ! (_nval=(_end-_start)/nbi)
    int _nval;
    std::string _profile;
    std::string _localization;
    //! only on assignement -3 : ON_NODES, -2 : ON_CELLS, -1 : ON_GAUSS_NE, 0..* : ON_GAUSS_PT
    mutable int _loc_id;
    mutable int _profile_it;
    MEDCouplingAutoRefCountObjectPtr<PartDefinition> _pd;
  public:
    mutable int _tmp_work1;
  };

  class MEDFileFieldPerMeshPerType : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMeshPerType *New(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType);
    static MEDFileFieldPerMeshPerType *NewOnRead(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, const MEDFileFieldNameScope& nasc, const PartDefinition *pd);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileFieldPerMeshPerType *deepCpy(MEDFileFieldPerMesh *father) const;
    void assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void assignFieldProfile(bool isPflAlone, int& start, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, DataArrayInt *locIds, int nbOfEltsInWholeMesh, const MEDCouplingFieldDouble *field, const DataArray *arr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob);
    void assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    const MEDFileFieldPerMesh *getFather() const;
    void loadOnlyStructureOfDataRecursively(med_idt fid, int &start, const MEDFileFieldNameScope& nasc);
    void loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc);
    void writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const;
    void getDimension(int& dim) const;
    void fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const;
    void fillFieldSplitedByType(std::vector< std::pair<int,int> >& dads, std::vector<TypeOfField>& types, std::vector<std::string>& pfls, std::vector<std::string>& locs) const;
    int getIteration() const;
    int getOrder() const;
    double getTime() const;
    std::string getMeshName() const;
    void simpleRepr(int bkOffset, std::ostream& oss, int id) const;
    void getSizes(int& globalSz, int& nbOfEntries) const;
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    int getNumberOfComponents() const;
    bool presenceOfMultiDiscPerGeoType() const;
    DataArray *getOrCreateAndGetArray();
    const DataArray *getOrCreateAndGetArray() const;
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenLocId(int locId);
    const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenLocId(int locId) const;
    void getFieldAtLevel(int meshDim, TypeOfField type, const MEDFileFieldGlobsReal *glob, std::vector< std::pair<int,int> >& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const;
    void fillValues(int& startEntryId, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    void setLeaves(const std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerTypePerDisc > >& leaves);
    bool keepOnlySpatialDiscretization(TypeOfField tof, int &globalNum, std::vector< std::pair<int,int> >& its);
    bool keepOnlyGaussDiscretization(std::size_t idOfDisc, int &globalNum, std::vector< std::pair<int,int> >& its);
    static med_entity_type ConvertIntoMEDFileType(TypeOfField ikType, INTERP_KERNEL::NormalizedCellType ikGeoType, med_geometry_type& medfGeoType);
  private:
    std::vector<int> addNewEntryIfNecessary(const MEDCouplingFieldDouble *field, int offset, int nbOfCells);
    std::vector<int> addNewEntryIfNecessaryGauss(const MEDCouplingFieldDouble *field, int offset, int nbOfCells);
    std::vector<int> addNewEntryIfNecessary(const MEDCouplingFieldDouble *field, const DataArrayInt *subCells);
    std::vector<int> addNewEntryIfNecessaryGauss(const MEDCouplingFieldDouble *field, const DataArrayInt *subCells);
    MEDFileFieldPerMeshPerType(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, const MEDFileFieldNameScope& nasc, const PartDefinition *pd);
    MEDFileFieldPerMeshPerType(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType);
  private:
    MEDFileFieldPerMesh *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> > _field_pm_pt_pd;
    INTERP_KERNEL::NormalizedCellType _geo_type;
  };

  class MEDFileMesh;

  class MEDFileFieldPerMesh : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMesh *New(MEDFileAnyTypeField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh);
    static MEDFileFieldPerMesh *NewOnRead(med_idt fid, MEDFileAnyTypeField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder, const MEDFileFieldNameScope& nasc, const MEDFileMesh *mm, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileFieldPerMesh *deepCpy(MEDFileAnyTypeField1TSWithoutSDA *father) const;
    void simpleRepr(int bkOffset,std::ostream& oss, int id) const;
    void copyTinyInfoFrom(const MEDCouplingMesh *mesh);
    void assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const std::vector<int>& code, const std::vector<int>& code2, const std::vector<DataArrayInt *>& idsInPflPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, const DataArray *arr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void assignFieldNoProfileNoRenum(int& start, const std::vector<int>& code, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob);
    void assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    void loadOnlyStructureOfDataRecursively(med_idt fid, int &start, const MEDFileFieldNameScope& nasc);
    void loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc);
    void writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const;
    void fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const;
    std::vector< std::vector< std::pair<int,int> > > getFieldSplitedByType(std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    void getDimension(int& dim) const;
    double getTime() const;
    int getIteration() const;
    int getOrder() const;
    int getMeshIteration() const { return _mesh_iteration; }
    int getMeshOrder() const { return _mesh_order; }
    std::string getMeshName() const { return _mesh_name; }
    int getNumberOfComponents() const;
    bool presenceOfMultiDiscPerGeoType() const;
    DataArray *getOrCreateAndGetArray();
    const DataArray *getOrCreateAndGetArray() const;
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    bool renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N, MEDFileFieldGlobsReal& glob);
    void keepOnlySpatialDiscretization(TypeOfField tof, int &globalNum, std::vector< std::pair<int,int> >& its);
    void keepOnlyGaussDiscretization(std::size_t idOfDisc, int &globalNum, std::vector< std::pair<int,int> >& its);
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, bool& isPfl, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    DataArray *getFieldOnMeshAtLevelWithPfl(TypeOfField type, const MEDCouplingMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob, const MEDFileFieldNameScope& nasc) const;
    void getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, int locId);
    const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, int locId) const;
  private:
    int addNewEntryIfNecessary(INTERP_KERNEL::NormalizedCellType type);
    MEDCouplingFieldDouble *finishField(TypeOfField type, const MEDFileFieldGlobsReal *glob,
                                        const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs, const MEDCouplingMesh *mesh, bool& isPfl, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    MEDCouplingFieldDouble *finishField2(TypeOfField type, const MEDFileFieldGlobsReal *glob,
                                         const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs,
                                         const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes,
                                         const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    MEDCouplingFieldDouble *finishFieldNode2(const MEDFileFieldGlobsReal *glob,
                                             const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs,
                                             const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    DataArray *finishField4(const std::vector< std::pair<int,int> >& dads, const DataArrayInt *pflIn, int nbOfElems, DataArrayInt *&pflOut) const;
    void assignNewLeaves(const std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerTypePerDisc > >& leaves);
    static void SortArraysPerType(const MEDFileFieldGlobsReal *glob, TypeOfField type, 
                                  const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& locs,
                                  std::vector<int>& code, std::vector<DataArrayInt *>& notNullPfls);
    static int ComputeNbOfElems(const MEDFileFieldGlobsReal *glob, TypeOfField type, const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs);
    MEDFileFieldPerMesh(med_idt fid, MEDFileAnyTypeField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder, const MEDFileFieldNameScope& nasc, const MEDFileMesh *mm, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    MEDFileFieldPerMesh(MEDFileAnyTypeField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh);
  private:
    std::string _mesh_name;
    int _mesh_iteration;
    int _mesh_order;
    MEDFileAnyTypeField1TSWithoutSDA *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > > _field_pm_pt;
  };

  class MEDFileFieldGlobsReal;

  class MEDFileFieldGlobs : public RefCountObject
  {
  public:
    static MEDFileFieldGlobs *New(const std::string& fname);
    static MEDFileFieldGlobs *New();
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileFieldGlobs *deepCpy() const;
    MEDFileFieldGlobs *shallowCpyPart(const std::vector<std::string>& pfls, const std::vector<std::string>& locs) const;
    MEDFileFieldGlobs *deepCpyPart(const std::vector<std::string>& pfls, const std::vector<std::string>& locs) const;
    void simpleRepr(std::ostream& oss) const;
    void appendGlobs(const MEDFileFieldGlobs& other, double eps);
    void checkGlobsPflsPartCoherency(const std::vector<std::string>& pflsUsed) const;
    void checkGlobsLocsPartCoherency(const std::vector<std::string>& locsUsed) const;
    void loadProfileInFile(med_idt fid, int id, const std::string& pflName);
    void loadProfileInFile(med_idt fid, int id);
    void loadGlobals(med_idt fid, const MEDFileFieldGlobsReal& real);
    void loadAllGlobals(med_idt fid);
    void writeGlobals(med_idt fid, const MEDFileWritable& opt) const;
    std::vector<std::string> getPfls() const;
    std::vector<std::string> getLocs() const;
    bool existsPfl(const std::string& pflName) const;
    bool existsLoc(const std::string& locName) const;
    std::string createNewNameOfPfl() const;
    std::string createNewNameOfLoc() const;
    std::vector< std::vector<int> > whichAreEqualProfiles() const;
    std::vector< std::vector<int> > whichAreEqualLocs(double eps) const;
    void setFileName(const std::string& fileName);
    void changePflsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    void changeLocsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    int getNbOfGaussPtPerCell(int locId) const;
    int getLocalizationId(const std::string& loc) const;
    std::string getFileName() const { return _file_name; }
    const MEDFileFieldLoc& getLocalizationFromId(int locId) const;
    const MEDFileFieldLoc& getLocalization(const std::string& locName) const;
    const DataArrayInt *getProfileFromId(int pflId) const;
    const DataArrayInt *getProfile(const std::string& pflName) const;
    MEDFileFieldLoc& getLocalizationFromId(int locId);
    MEDFileFieldLoc& getLocalization(const std::string& locName);
    DataArrayInt *getProfile(const std::string& pflName);
    DataArrayInt *getProfileFromId(int pflId);
    void killProfileIds(const std::vector<int>& pflIds);
    void killLocalizationIds(const std::vector<int>& locIds);
    //
    void appendProfile(DataArrayInt *pfl);
    void appendLoc(const std::string& locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w);
    //
    static std::string CreateNewNameNotIn(const std::string& prefix, const std::vector<std::string>& namesToAvoid);
  protected:
    MEDFileFieldGlobs(const std::string& fname);
    MEDFileFieldGlobs();
    ~MEDFileFieldGlobs();
  protected:
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > _pfls;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> > _locs;
    std::string _file_name;
  };

  /// @endcond INTERNAL

  class MEDFileFieldGlobsReal
  {
  public:
    MEDLOADER_EXPORT MEDFileFieldGlobsReal(const std::string& fname);
    MEDLOADER_EXPORT MEDFileFieldGlobsReal();
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT void simpleReprGlobs(std::ostream& oss) const;
    MEDLOADER_EXPORT void resetContent();
    MEDLOADER_EXPORT void shallowCpyGlobs(const MEDFileFieldGlobsReal& other);
    MEDLOADER_EXPORT void deepCpyGlobs(const MEDFileFieldGlobsReal& other);
    MEDLOADER_EXPORT void shallowCpyOnlyUsedGlobs(const MEDFileFieldGlobsReal& other);
    MEDLOADER_EXPORT void deepCpyOnlyUsedGlobs(const MEDFileFieldGlobsReal& other);
    MEDLOADER_EXPORT void appendGlobs(const MEDFileFieldGlobsReal& other, double eps);
    MEDLOADER_EXPORT void checkGlobsCoherency() const;
    MEDLOADER_EXPORT void checkGlobsPflsPartCoherency() const;
    MEDLOADER_EXPORT void checkGlobsLocsPartCoherency() const;
    MEDLOADER_EXPORT virtual std::vector<std::string> getPflsReallyUsed() const = 0;
    MEDLOADER_EXPORT virtual std::vector<std::string> getLocsReallyUsed() const = 0;
    MEDLOADER_EXPORT virtual std::vector<std::string> getPflsReallyUsedMulti() const = 0;
    MEDLOADER_EXPORT virtual std::vector<std::string> getLocsReallyUsedMulti() const = 0;
    MEDLOADER_EXPORT virtual void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) = 0;
    MEDLOADER_EXPORT virtual void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) = 0;
    MEDLOADER_EXPORT virtual ~MEDFileFieldGlobsReal();
    //
    MEDLOADER_EXPORT void loadProfileInFile(med_idt fid, int id, const std::string& pflName);
    MEDLOADER_EXPORT void loadProfileInFile(med_idt fid, int id);
    MEDLOADER_EXPORT void loadGlobals(med_idt fid);
    MEDLOADER_EXPORT void loadAllGlobals(med_idt fid);
    MEDLOADER_EXPORT void writeGlobals(med_idt fid, const MEDFileWritable& opt) const;
    MEDLOADER_EXPORT std::vector<std::string> getPfls() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocs() const;
    MEDLOADER_EXPORT bool existsPfl(const std::string& pflName) const;
    MEDLOADER_EXPORT bool existsLoc(const std::string& locName) const;
    MEDLOADER_EXPORT std::string createNewNameOfPfl() const;
    MEDLOADER_EXPORT std::string createNewNameOfLoc() const;
    MEDLOADER_EXPORT std::vector< std::vector<int> > whichAreEqualProfiles() const;
    MEDLOADER_EXPORT std::vector< std::vector<int> > whichAreEqualLocs(double eps) const;
    MEDLOADER_EXPORT void setFileName(const std::string& fileName);
    MEDLOADER_EXPORT void changePflsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changeLocsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changePflsNames(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changeLocsNames(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changePflName(const std::string& oldName, const std::string& newName);
    MEDLOADER_EXPORT void changeLocName(const std::string& oldName, const std::string& newName);
    MEDLOADER_EXPORT std::vector< std::pair<std::vector<std::string>, std::string > > zipPflsNames();
    MEDLOADER_EXPORT std::vector< std::pair<std::vector<std::string>, std::string > > zipLocsNames(double eps);
    MEDLOADER_EXPORT int getNbOfGaussPtPerCell(int locId) const;
    MEDLOADER_EXPORT int getLocalizationId(const std::string& loc) const;
    MEDLOADER_EXPORT std::string getFileName() const;
    MEDLOADER_EXPORT const MEDFileFieldLoc& getLocalizationFromId(int locId) const;
    MEDLOADER_EXPORT const MEDFileFieldLoc& getLocalization(const std::string& locName) const;
    MEDLOADER_EXPORT MEDFileFieldLoc& getLocalizationFromId(int locId);
    MEDLOADER_EXPORT MEDFileFieldLoc& getLocalization(const std::string& locName);
    MEDLOADER_EXPORT const DataArrayInt *getProfile(const std::string& pflName) const;
    MEDLOADER_EXPORT const DataArrayInt *getProfileFromId(int pflId) const;
    MEDLOADER_EXPORT DataArrayInt *getProfile(const std::string& pflName);
    MEDLOADER_EXPORT DataArrayInt *getProfileFromId(int pflId);
    MEDLOADER_EXPORT void killProfileIds(const std::vector<int>& pflIds);
    MEDLOADER_EXPORT void killLocalizationIds(const std::vector<int>& locIds);
    //
    MEDLOADER_EXPORT void appendProfile(DataArrayInt *pfl);
    MEDLOADER_EXPORT void appendLoc(const std::string& locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w);
  protected:
    MEDFileFieldGlobs *contentNotNull();
    const MEDFileFieldGlobs *contentNotNull() const;
  protected:
    MEDCouplingAutoRefCountObjectPtr< MEDFileFieldGlobs > _globals;
  };

  class MEDFileFieldNameScope
  {
  public:
    MEDLOADER_EXPORT MEDFileFieldNameScope();
    MEDLOADER_EXPORT MEDFileFieldNameScope(const std::string& fieldName);
    MEDLOADER_EXPORT std::string getName() const;
    MEDLOADER_EXPORT void setName(const std::string& fieldName);
    MEDLOADER_EXPORT std::string getDtUnit() const;
    MEDLOADER_EXPORT void setDtUnit(const std::string& dtUnit);
    MEDLOADER_EXPORT void copyNameScope(const MEDFileFieldNameScope& other);
  protected:
    std::string _name;
    std::string _dt_unit;
  };

  class MEDFileMeshes;

  /*!
   * SDA is for Shared Data Arrays such as profiles.
   */
  class MEDFileAnyTypeField1TSWithoutSDA : public RefCountObject, public MEDFileFieldNameScope
  {
  public:
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA();
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA(const std::string& fieldName, int csit, int iteration, int order);
    MEDLOADER_EXPORT int getIteration() const { return _iteration; }
    MEDLOADER_EXPORT int getOrder() const { return _order; }
    MEDLOADER_EXPORT double getTime(int& iteration, int& order) const { iteration=_iteration; order=_order; return _dt; }
    MEDLOADER_EXPORT void setTime(int iteration, int order, double val) { _dt=val; _iteration=iteration; _order=order; }
    MEDLOADER_EXPORT int getDimension() const;
    MEDLOADER_EXPORT std::string getMeshName() const;
    MEDLOADER_EXPORT void setMeshName(const std::string& newMeshName);
    MEDLOADER_EXPORT bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT int getMeshIteration() const;
    MEDLOADER_EXPORT int getMeshOrder() const;
    MEDLOADER_EXPORT bool isDealingTS(int iteration, int order) const;
    MEDLOADER_EXPORT std::pair<int,int> getDtIt() const;
    MEDLOADER_EXPORT void fillIteration(std::pair<int,int>& p) const;
    MEDLOADER_EXPORT void fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const;
    MEDLOADER_EXPORT std::vector<TypeOfField> getTypesOfFieldAvailable() const;
    //
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsed2() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsed2() const;
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsedMulti2() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsedMulti2() const;
    MEDLOADER_EXPORT void changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    //
    MEDLOADER_EXPORT int getNonEmptyLevels(const std::string& mname, std::vector<int>& levs) const;
    MEDLOADER_EXPORT std::vector< std::vector<std::pair<int,int> > > getFieldSplitedByType(const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    //
    MEDLOADER_EXPORT MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId);
    MEDLOADER_EXPORT const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const;
    MEDLOADER_EXPORT void deepCpyLeavesFrom(const MEDFileAnyTypeField1TSWithoutSDA& other);
  public:
    MEDLOADER_EXPORT int getNumberOfComponents() const;
    MEDLOADER_EXPORT const std::vector<std::string>& getInfo() const;
    MEDLOADER_EXPORT std::vector<std::string>& getInfo();
    MEDLOADER_EXPORT bool presenceOfMultiDiscPerGeoType() const;
    MEDLOADER_EXPORT void setInfo(const std::vector<std::string>& infos);
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT int copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr);
    MEDLOADER_EXPORT void setFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    MEDLOADER_EXPORT void setFieldProfile(const MEDCouplingFieldDouble *field, const DataArray *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    MEDLOADER_EXPORT virtual void simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeField1TSWithoutSDA *deepCpy() const = 0;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeField1TSWithoutSDA *shallowCpy() const = 0;
    MEDLOADER_EXPORT virtual std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > splitComponents() const;
    MEDLOADER_EXPORT virtual const char *getTypeStr() const = 0;
    MEDLOADER_EXPORT virtual DataArray *getUndergroundDataArray() const = 0;
    MEDLOADER_EXPORT virtual DataArray *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const = 0;
    MEDLOADER_EXPORT virtual void setArray(DataArray *arr) = 0;
    MEDLOADER_EXPORT virtual DataArray *createNewEmptyDataArrayInstance() const = 0;
    MEDLOADER_EXPORT virtual DataArray *getOrCreateAndGetArray() = 0;
    MEDLOADER_EXPORT virtual const DataArray *getOrCreateAndGetArray() const = 0;
  public:
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, const std::string& mName, int renumPol, const MEDFileFieldGlobsReal *glob, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDFileMesh *mesh, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, const std::string& mName, int renumPol, const MEDFileFieldGlobsReal *glob, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, const DataArrayInt *cellRenum, const DataArrayInt *nodeRenum, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    DataArray *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob, const MEDFileFieldNameScope& nasc) const;
  public:
    MEDLOADER_EXPORT bool renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N, MEDFileFieldGlobsReal& glob);
    MEDLOADER_EXPORT std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > splitDiscretizations() const;
    MEDLOADER_EXPORT std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > splitMultiDiscrPerGeoTypes() const;
    MEDLOADER_EXPORT int keepOnlySpatialDiscretization(TypeOfField tof, std::vector< std::pair<int,int> >& its);
    MEDLOADER_EXPORT int keepOnlyGaussDiscretization(std::size_t idOfDisc, std::vector< std::pair<int,int> >& its);
  public:
    MEDLOADER_EXPORT void allocNotFromFile(int newNbOfTuples);
    MEDLOADER_EXPORT bool allocIfNecessaryTheArrayToReceiveDataFromFile();
    MEDLOADER_EXPORT void loadOnlyStructureOfDataRecursively(med_idt fid, const MEDFileFieldNameScope& nasc, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    MEDLOADER_EXPORT void loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc);
    MEDLOADER_EXPORT void loadBigArraysRecursivelyIfNecessary(med_idt fid, const MEDFileFieldNameScope& nasc);
    MEDLOADER_EXPORT void loadStructureAndBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    MEDLOADER_EXPORT void unloadArrays();
    MEDLOADER_EXPORT void writeLL(med_idt fid, const MEDFileWritable& opts, const MEDFileFieldNameScope& nasc) const;
  protected:
    int getMeshIdFromMeshName(const std::string& mName) const;
    int addNewEntryIfNecessary(const MEDCouplingMesh *mesh);
    void updateData(int newLgth, const std::vector< std::pair<int,int> >& oldStartStops);
  protected:
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > > _field_per_mesh;
    int _iteration;
    int _order;
    double _dt;
  public:
    //! only useable on reading
    mutable int _csit;
    // -3 means allocated and build from scratch
    // -2 means allocated and read from a file
    // -1 means not allocated and build from scratch
    // >=0 means not allocated and read from a file
    mutable int _nb_of_tuples_to_be_allocated;
  };

  class MEDFileIntField1TSWithoutSDA;

  /*!
   * SDA is for Shared Data Arrays such as profiles.
   */
  class MEDFileField1TSWithoutSDA : public MEDFileAnyTypeField1TSWithoutSDA
  {
  public:
    MEDLOADER_EXPORT const char *getTypeStr() const;
    MEDLOADER_EXPORT DataArray *getUndergroundDataArray() const;
    MEDLOADER_EXPORT DataArray *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDLOADER_EXPORT DataArrayDouble *getUndergroundDataArrayDouble() const;
    MEDLOADER_EXPORT DataArrayDouble *getUndergroundDataArrayDoubleExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDLOADER_EXPORT std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    MEDLOADER_EXPORT static void CheckMeshDimRel(int meshDimRelToMax);
    MEDLOADER_EXPORT static std::vector<int> CheckSBTMesh(const MEDCouplingMesh *mesh);
    MEDLOADER_EXPORT static MEDFileField1TSWithoutSDA *New(const std::string& fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos);
  public:
    MEDLOADER_EXPORT MEDFileField1TSWithoutSDA();
    MEDLOADER_EXPORT MEDFileField1TSWithoutSDA(const std::string& fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos);
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA *shallowCpy() const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA *deepCpy() const;
    MEDLOADER_EXPORT void setArray(DataArray *arr);
    MEDLOADER_EXPORT DataArray *createNewEmptyDataArrayInstance() const;
    MEDLOADER_EXPORT DataArray *getOrCreateAndGetArray();
    MEDLOADER_EXPORT const DataArray *getOrCreateAndGetArray() const;
    MEDLOADER_EXPORT DataArrayDouble *getOrCreateAndGetArrayDouble();
    MEDLOADER_EXPORT const DataArrayDouble *getOrCreateAndGetArrayDouble() const;
    MEDLOADER_EXPORT MEDFileIntField1TSWithoutSDA *convertToInt() const;
  protected:
    MEDCouplingAutoRefCountObjectPtr< DataArrayDouble > _arr;
  public:
    static const char TYPE_STR[];
  };

  /*!
   * SDA is for Shared Data Arrays such as profiles.
   */
  class MEDFileIntField1TSWithoutSDA : public MEDFileAnyTypeField1TSWithoutSDA
  {
  public:
    MEDLOADER_EXPORT MEDFileIntField1TSWithoutSDA();
    MEDLOADER_EXPORT static MEDFileIntField1TSWithoutSDA *New(const std::string& fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos);
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA *deepCpy() const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA *shallowCpy() const;
    MEDLOADER_EXPORT const char *getTypeStr() const;
    MEDLOADER_EXPORT DataArray *getUndergroundDataArray() const;
    MEDLOADER_EXPORT DataArray *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDLOADER_EXPORT void setArray(DataArray *arr);
    MEDLOADER_EXPORT DataArray *createNewEmptyDataArrayInstance() const;
    MEDLOADER_EXPORT DataArray *getOrCreateAndGetArray();
    MEDLOADER_EXPORT const DataArray *getOrCreateAndGetArray() const;
    MEDLOADER_EXPORT DataArrayInt *getOrCreateAndGetArrayInt();
    MEDLOADER_EXPORT const DataArrayInt *getOrCreateAndGetArrayInt() const;
    MEDLOADER_EXPORT DataArrayInt *getUndergroundDataArrayInt() const;
    MEDLOADER_EXPORT DataArrayInt *getUndergroundDataArrayIntExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDLOADER_EXPORT MEDFileField1TSWithoutSDA *convertToDouble() const;
  protected:
    MEDFileIntField1TSWithoutSDA(const std::string& fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos);
  protected:
    MEDCouplingAutoRefCountObjectPtr< DataArrayInt > _arr;
  public:
    MEDLOADER_EXPORT static const char TYPE_STR[];
  };

  /*!
   * User class.
   */
  class MEDFileAnyTypeField1TS : public RefCountObject, public MEDFileWritable, public MEDFileFieldGlobsReal
  {
  protected:
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS();
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS(const std::string& fileName, bool loadAll, const MEDFileMeshes *ms);
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS(const std::string& fileName, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms);
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms);
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS(const MEDFileAnyTypeField1TSWithoutSDA& other, bool shallowCopyOfContent);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *BuildNewInstanceFromContent(MEDFileAnyTypeField1TSWithoutSDA *c, const std::string& fileName);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TSWithoutSDA *BuildContentFrom(med_idt fid, const std::string& fileName, bool loadAll, const MEDFileMeshes *ms);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TSWithoutSDA *BuildContentFrom(med_idt fid, const std::string& fileName, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TSWithoutSDA *BuildContentFrom(med_idt fid, const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms);
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    // direct forwarding to MEDFileAnyTypeField1TSWithoutSDA instance _content
  public:
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *New(const std::string& fileName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true);
    MEDLOADER_EXPORT int getDimension() const;
    MEDLOADER_EXPORT int getIteration() const;
    MEDLOADER_EXPORT int getOrder() const;
    MEDLOADER_EXPORT double getTime(int& iteration, int& order) const;
    MEDLOADER_EXPORT void setTime(int iteration, int order, double val);
    MEDLOADER_EXPORT std::string getName() const;
    MEDLOADER_EXPORT void setName(const std::string& name);
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT void simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const;
    MEDLOADER_EXPORT std::string getDtUnit() const;
    MEDLOADER_EXPORT void setDtUnit(const std::string& dtUnit);
    MEDLOADER_EXPORT std::string getMeshName() const;
    MEDLOADER_EXPORT void setMeshName(const std::string& newMeshName);
    MEDLOADER_EXPORT bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT int getMeshIteration() const;
    MEDLOADER_EXPORT int getMeshOrder() const;
    MEDLOADER_EXPORT int getNumberOfComponents() const;
    MEDLOADER_EXPORT bool isDealingTS(int iteration, int order) const;
    MEDLOADER_EXPORT std::pair<int,int> getDtIt() const;
    MEDLOADER_EXPORT void fillIteration(std::pair<int,int>& p) const;
    MEDLOADER_EXPORT void fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const;
    MEDLOADER_EXPORT void setInfo(const std::vector<std::string>& infos);
    MEDLOADER_EXPORT const std::vector<std::string>& getInfo() const;
    MEDLOADER_EXPORT std::vector<std::string>& getInfo();
    MEDLOADER_EXPORT bool presenceOfMultiDiscPerGeoType() const;
    MEDLOADER_EXPORT std::vector<TypeOfField> getTypesOfFieldAvailable() const;
    MEDLOADER_EXPORT std::vector< std::vector<std::pair<int,int> > > getFieldSplitedByType(const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
        std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    MEDLOADER_EXPORT MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId);
    MEDLOADER_EXPORT const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const;
    MEDLOADER_EXPORT int getNonEmptyLevels(const std::string& mname, std::vector<int>& levs) const;
  public:
    MEDLOADER_EXPORT void write(const std::string& fileName, int mode) const;
    MEDLOADER_EXPORT void loadArrays();
    MEDLOADER_EXPORT void loadArraysIfNecessary();
    MEDLOADER_EXPORT void unloadArrays();
    MEDLOADER_EXPORT void unloadArraysWithoutDataLoss();
    MEDLOADER_EXPORT std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TS > > splitComponents() const;
    MEDLOADER_EXPORT std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TS > > splitDiscretizations() const;
    MEDLOADER_EXPORT std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TS > > splitMultiDiscrPerGeoTypes() const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *deepCpy() const;
    MEDLOADER_EXPORT int copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr);
    MEDLOADER_EXPORT virtual MEDFileAnyTypeField1TS *shallowCpy() const = 0;
  public:
    //! underground method see MEDFileField1TSWithoutSDA::setProfileNameOnLeaf
    MEDLOADER_EXPORT void setProfileNameOnLeaf(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newPflName, bool forceRenameOnGlob=false);
    //! underground method see MEDFileField1TSWithoutSDA::setLocNameOnLeaf
    MEDLOADER_EXPORT void setLocNameOnLeaf(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newLocName, bool forceRenameOnGlob=false);
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsed() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsed() const;
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsedMulti() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsedMulti() const;
    MEDLOADER_EXPORT void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
  public:
    MEDLOADER_EXPORT static int LocateField2(med_idt fid, const std::string& fileName, int fieldIdCFormat, bool checkFieldId, std::string& fieldName, med_field_type& typcha, std::vector<std::string>& infos, std::string& dtunitOut);
    MEDLOADER_EXPORT static int LocateField(med_idt fid, const std::string& fileName, const std::string& fieldName, int& posCFormat, med_field_type& typcha, std::vector<std::string>& infos, std::string& dtunitOut);
  public:
    MEDLOADER_EXPORT virtual med_field_type getMEDFileFieldType() const = 0;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA *contentNotNullBase();
    MEDLOADER_EXPORT const MEDFileAnyTypeField1TSWithoutSDA *contentNotNullBase() const;
  protected:
    MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> _content;
  };

  class MEDFileIntField1TS;

  /*!
   * User class.
   */
  class MEDFileField1TS : public MEDFileAnyTypeField1TS
  {
  public:
    MEDLOADER_EXPORT static MEDFileField1TS *New(const std::string& fileName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileField1TS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileField1TS *New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileField1TS *New(const MEDFileField1TSWithoutSDA& other, bool shallowCopyOfContent);
    MEDLOADER_EXPORT static MEDFileField1TS *New();
    MEDLOADER_EXPORT MEDFileIntField1TS *convertToInt(bool isDeepCpyGlobs=true) const;
    //
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, const std::string& mname, int meshDimRelToMax, int renumPol=0) const;
    MEDLOADER_EXPORT DataArrayDouble *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const;
    //
    MEDLOADER_EXPORT void setFieldNoProfileSBT(const MEDCouplingFieldDouble *field);
    MEDLOADER_EXPORT void setFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    // direct forwarding to MEDFileField1TSWithoutSDA instance _content
  public:
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *shallowCpy() const;
    MEDLOADER_EXPORT DataArrayDouble *getUndergroundDataArray() const;
    MEDLOADER_EXPORT DataArrayDouble *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;

    MEDLOADER_EXPORT std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
        std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
  public:
    MEDLOADER_EXPORT static void SetDataArrayDoubleInField(MEDCouplingFieldDouble *f, MEDCouplingAutoRefCountObjectPtr<DataArray>& arr);
    MEDLOADER_EXPORT static DataArrayDouble *ReturnSafelyDataArrayDouble(MEDCouplingAutoRefCountObjectPtr<DataArray>& arr);
  private:
    med_field_type getMEDFileFieldType() const { return MED_FLOAT64; }
    const MEDFileField1TSWithoutSDA *contentNotNull() const;
    MEDFileField1TSWithoutSDA *contentNotNull();
  private:
    ~MEDFileField1TS() { }
    MEDFileField1TS(const std::string& fileName, bool loadAll, const MEDFileMeshes *ms);
    MEDFileField1TS(const std::string& fileName, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms);
    MEDFileField1TS(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms);
    MEDFileField1TS(const MEDFileField1TSWithoutSDA& other, bool shallowCopyOfContent);
    MEDFileField1TS();
  };

  class MEDFileIntField1TS : public MEDFileAnyTypeField1TS
  {
  public:
    MEDLOADER_EXPORT static MEDFileIntField1TS *New();
    MEDLOADER_EXPORT static MEDFileIntField1TS *New(const std::string& fileName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileIntField1TS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileIntField1TS *New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileIntField1TS *New(const MEDFileIntField1TSWithoutSDA& other, bool shallowCopyOfContent);
    MEDLOADER_EXPORT MEDFileField1TS *convertToDouble(bool isDeepCpyGlobs=true) const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *shallowCpy() const;
    //
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, DataArrayInt* &arrOut, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, DataArrayInt* &arrOut, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt* &arrOut, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, DataArrayInt* &arrOut, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, const std::string& mname, int meshDimRelToMax, DataArrayInt* &arrOut, int renumPol=0) const;
    MEDLOADER_EXPORT DataArrayInt *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const;
    //
    MEDLOADER_EXPORT void setFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals);
    MEDLOADER_EXPORT void setFieldProfile(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    MEDLOADER_EXPORT DataArrayInt *getUndergroundDataArray() const;
  public:
    MEDLOADER_EXPORT static DataArrayInt *ReturnSafelyDataArrayInt(MEDCouplingAutoRefCountObjectPtr<DataArray>& arr);
  private:
    med_field_type getMEDFileFieldType() const { return MED_INT32; }
    const MEDFileIntField1TSWithoutSDA *contentNotNull() const;
    MEDFileIntField1TSWithoutSDA *contentNotNull();
  private:
    ~MEDFileIntField1TS() { }
    MEDFileIntField1TS();
    MEDFileIntField1TS(const std::string& fileName, bool loadAll, const MEDFileMeshes *ms);
    MEDFileIntField1TS(const std::string& fileName, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms);
    MEDFileIntField1TS(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms);
    MEDFileIntField1TS(const MEDFileIntField1TSWithoutSDA& other, bool shallowCopyOfContent);
  };

  class MEDFileAnyTypeFieldMultiTSWithoutSDA : public RefCountObject, public MEDFileFieldNameScope
  {
  protected:
    MEDFileAnyTypeFieldMultiTSWithoutSDA();
    MEDFileAnyTypeFieldMultiTSWithoutSDA(const std::string& fieldName);
    MEDFileAnyTypeFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    MEDFileAnyTypeFieldMultiTSWithoutSDA(med_idt fid, const std::string& fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
  public:
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeFieldMultiTSWithoutSDA *deepCpy() const;
    MEDLOADER_EXPORT virtual std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > splitComponents() const;
    MEDLOADER_EXPORT virtual std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > splitDiscretizations() const;
    MEDLOADER_EXPORT virtual std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > splitMultiDiscrPerGeoTypes() const;
    MEDLOADER_EXPORT virtual const char *getTypeStr() const = 0;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeFieldMultiTSWithoutSDA *shallowCpy() const = 0;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeFieldMultiTSWithoutSDA *createNew() const = 0;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeField1TSWithoutSDA *createNew1TSWithoutSDAEmptyInstance() const = 0;
    MEDLOADER_EXPORT virtual void checkCoherencyOfType(const MEDFileAnyTypeField1TSWithoutSDA *f1ts) const = 0;
    MEDLOADER_EXPORT const std::vector<std::string>& getInfo() const;
    MEDLOADER_EXPORT bool presenceOfMultiDiscPerGeoType() const;
    MEDLOADER_EXPORT void setInfo(const std::vector<std::string>& info);
    MEDLOADER_EXPORT int getTimeStepPos(int iteration, int order) const;
    MEDLOADER_EXPORT const MEDFileAnyTypeField1TSWithoutSDA& getTimeStepEntry(int iteration, int order) const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA& getTimeStepEntry(int iteration, int order);
    MEDLOADER_EXPORT std::string getMeshName() const;
    MEDLOADER_EXPORT void setMeshName(const std::string& newMeshName);
    MEDLOADER_EXPORT bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT int getNumberOfTS() const;
    MEDLOADER_EXPORT void eraseEmptyTS();
    MEDLOADER_EXPORT void eraseTimeStepIds(const int *startIds, const int *endIds);
    MEDLOADER_EXPORT void eraseTimeStepIds2(int bg, int end, int step);
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *buildFromTimeStepIds(const int *startIds, const int *endIds) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *buildFromTimeStepIds2(int bg, int end, int step) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *partOfThisLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *partOfThisNotLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const;
    MEDLOADER_EXPORT int getPosOfTimeStep(int iteration, int order) const;
    MEDLOADER_EXPORT int getPosGivenTime(double time, double eps=1e-8) const;
    MEDLOADER_EXPORT std::vector< std::pair<int,int> > getIterations() const;
    MEDLOADER_EXPORT std::vector< std::pair<int,int> > getTimeSteps(std::vector<double>& ret1) const;
    MEDLOADER_EXPORT void pushBackTimeStep(MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA>& tse);
    MEDLOADER_EXPORT void synchronizeNameScope();
    MEDLOADER_EXPORT void simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const;
    MEDLOADER_EXPORT int getNonEmptyLevels(int iteration, int order, const std::string& mname, std::vector<int>& levs) const;
    MEDLOADER_EXPORT void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob);
    MEDLOADER_EXPORT void appendFieldProfile(const MEDCouplingFieldDouble *field, const DataArray *arr, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob);
    MEDLOADER_EXPORT std::vector< std::vector< std::pair<int,int> > > getFieldSplitedByType(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    MEDLOADER_EXPORT std::vector< std::vector<TypeOfField> > getTypesOfFieldAvailable() const;
    MEDLOADER_EXPORT DataArray *getUndergroundDataArray(int iteration, int order) const;
    MEDLOADER_EXPORT DataArray *getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDLOADER_EXPORT bool renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N, MEDFileFieldGlobsReal& glob);
    MEDLOADER_EXPORT void loadStructureOrStructureAndBigArraysRecursively(med_idt fid, int nbPdt, med_field_type fieldTyp, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    MEDLOADER_EXPORT void writeLL(med_idt fid, const MEDFileWritable& opts) const;
    MEDLOADER_EXPORT void loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc);
    MEDLOADER_EXPORT void loadBigArraysRecursivelyIfNecessary(med_idt fid, const MEDFileFieldNameScope& nasc);
    MEDLOADER_EXPORT void unloadArrays();
  public:
    MEDLOADER_EXPORT const MEDFileAnyTypeField1TSWithoutSDA *getTimeStepAtPos2(int pos) const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA *getTimeStepAtPos2(int pos);
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsed2() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsed2() const;
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsedMulti2() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsedMulti2() const;
    MEDLOADER_EXPORT void changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void setIteration(int i, MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> ts);
  protected:
    virtual med_field_type getMEDFileFieldType() const = 0;
    void copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr);
    void checkCoherencyOfTinyInfo(const MEDCouplingFieldDouble *field, const DataArray *arr) const;
    void checkThatComponentsMatch(const std::vector<std::string>& compos) const;
    void checkThatNbOfCompoOfTSMatchThis() const;
  protected:
    std::vector<std::string> _infos;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > _time_steps;
  };

  class MEDFileIntFieldMultiTSWithoutSDA;

  class MEDFileFieldMultiTSWithoutSDA : public MEDFileAnyTypeFieldMultiTSWithoutSDA
  {
  public:
    MEDLOADER_EXPORT static MEDFileFieldMultiTSWithoutSDA *New(med_idt fid, const std::string& fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    MEDLOADER_EXPORT MEDFileFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    MEDLOADER_EXPORT const char *getTypeStr() const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *shallowCpy() const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *createNew() const;
    MEDLOADER_EXPORT std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    MEDLOADER_EXPORT MEDFileIntFieldMultiTSWithoutSDA *convertToInt() const;
  protected:
    MEDFileFieldMultiTSWithoutSDA(const std::string& fieldName);
    MEDFileFieldMultiTSWithoutSDA(med_idt fid, const std::string& fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    med_field_type getMEDFileFieldType() const { return MED_FLOAT64; }
    MEDFileAnyTypeField1TSWithoutSDA *createNew1TSWithoutSDAEmptyInstance() const;
    void checkCoherencyOfType(const MEDFileAnyTypeField1TSWithoutSDA *f1ts) const;
  public:
    MEDLOADER_EXPORT MEDFileFieldMultiTSWithoutSDA();
  };

  class MEDFileIntFieldMultiTSWithoutSDA : public MEDFileAnyTypeFieldMultiTSWithoutSDA
  {
  public:
    MEDLOADER_EXPORT static MEDFileIntFieldMultiTSWithoutSDA *New(med_idt fid, const std::string& fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    MEDLOADER_EXPORT MEDFileIntFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    MEDLOADER_EXPORT const char *getTypeStr() const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *shallowCpy() const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *createNew() const;
    MEDLOADER_EXPORT MEDFileFieldMultiTSWithoutSDA *convertToDouble() const;
  protected:
    MEDFileIntFieldMultiTSWithoutSDA(const std::string& fieldName);
    MEDFileIntFieldMultiTSWithoutSDA(med_idt fid, const std::string& fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    med_field_type getMEDFileFieldType() const { return MED_INT32; }
    MEDFileAnyTypeField1TSWithoutSDA *createNew1TSWithoutSDAEmptyInstance() const;
    void checkCoherencyOfType(const MEDFileAnyTypeField1TSWithoutSDA *f1ts) const;
  public:
    MEDLOADER_EXPORT MEDFileIntFieldMultiTSWithoutSDA();
  };

  class MEDFileAnyTypeFieldMultiTSIterator;
  class MEDFileFastCellSupportComparator;
  /*!
   * User class.
   */
  class MEDFileAnyTypeFieldMultiTS : public RefCountObject, public MEDFileWritable, public MEDFileFieldGlobsReal
  {
  protected:
    MEDFileAnyTypeFieldMultiTS();
    MEDFileAnyTypeFieldMultiTS(const std::string& fileName, bool loadAll, const MEDFileMeshes *ms);
    MEDFileAnyTypeFieldMultiTS(const std::string& fileName, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities=0);
    MEDFileAnyTypeFieldMultiTS(const MEDFileAnyTypeFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent);
    static MEDFileAnyTypeFieldMultiTS *BuildNewInstanceFromContent(MEDFileAnyTypeFieldMultiTSWithoutSDA *c, const std::string& fileName);
    static MEDFileAnyTypeFieldMultiTSWithoutSDA *BuildContentFrom(med_idt fid, const std::string& fileName, bool loadAll, const MEDFileMeshes *ms);
    static MEDFileAnyTypeFieldMultiTSWithoutSDA *BuildContentFrom(med_idt fid, const std::string& fileName, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
  public:
    MEDLOADER_EXPORT static MEDFileAnyTypeFieldMultiTS *New(const std::string& fileName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeFieldMultiTS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT void loadArrays();
    MEDLOADER_EXPORT void loadArraysIfNecessary();
    MEDLOADER_EXPORT void unloadArrays();
    MEDLOADER_EXPORT void unloadArraysWithoutDataLoss();
    MEDLOADER_EXPORT void write(const std::string& fileName, int mode) const;
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeFieldMultiTS *deepCpy() const;
    MEDLOADER_EXPORT std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTS > > splitComponents() const;
    MEDLOADER_EXPORT std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTS > > splitDiscretizations() const;
    MEDLOADER_EXPORT std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTS > > splitMultiDiscrPerGeoTypes() const;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeFieldMultiTS *shallowCpy() const = 0;
    MEDLOADER_EXPORT virtual void checkCoherencyOfType(const MEDFileAnyTypeField1TS *f1ts) const = 0;
    //
    MEDLOADER_EXPORT virtual MEDFileAnyTypeField1TS *getTimeStepAtPos(int pos) const = 0;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *getTimeStep(int iteration, int order) const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *getTimeStepGivenTime(double time, double eps=1e-8) const;
    MEDLOADER_EXPORT static std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > SplitIntoCommonTimeSeries(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS);
    MEDLOADER_EXPORT static std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > SplitPerCommonSupport(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS, const MEDFileMesh *mesh, std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFastCellSupportComparator> >& fsc);
    MEDLOADER_EXPORT static int CheckSupportAcrossTime(MEDFileAnyTypeFieldMultiTS *f0, MEDFileAnyTypeFieldMultiTS *f1, const MEDFileMesh *mesh, TypeOfField& tof0, TypeOfField& tof1);
  public:// direct forwarding to MEDFileField1TSWithoutSDA instance _content
    MEDLOADER_EXPORT std::string getName() const;
    MEDLOADER_EXPORT void setName(const std::string& name);
    MEDLOADER_EXPORT std::string getDtUnit() const;
    MEDLOADER_EXPORT void setDtUnit(const std::string& dtUnit);
    MEDLOADER_EXPORT std::string getMeshName() const;
    MEDLOADER_EXPORT void setMeshName(const std::string& newMeshName);
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT void simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const;
    MEDLOADER_EXPORT int getNumberOfTS() const;
    MEDLOADER_EXPORT void eraseEmptyTS();
    MEDLOADER_EXPORT void eraseTimeStepIds(const int *startIds, const int *endIds);
    MEDLOADER_EXPORT void eraseTimeStepIds2(int bg, int end, int step);
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *buildSubPart(const int *startIds, const int *endIds) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *buildSubPartSlice(int bg, int end, int step) const;
    MEDLOADER_EXPORT std::vector< std::pair<int,int> > getTimeSteps(std::vector<double>& ret1) const;
    MEDLOADER_EXPORT std::vector< std::pair<int,int> > getIterations() const;
    MEDLOADER_EXPORT void pushBackTimeSteps(const std::vector<MEDFileAnyTypeField1TS *>& f1ts);
    MEDLOADER_EXPORT void pushBackTimeSteps(MEDFileAnyTypeFieldMultiTS *fmts);
    MEDLOADER_EXPORT void pushBackTimeStep(MEDFileAnyTypeField1TS *f1ts);
    MEDLOADER_EXPORT void synchronizeNameScope();
    MEDLOADER_EXPORT int getPosOfTimeStep(int iteration, int order) const;
    MEDLOADER_EXPORT int getPosGivenTime(double time, double eps=1e-8) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSIterator *iterator();
    MEDLOADER_EXPORT bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT const std::vector<std::string>& getInfo() const;
    MEDLOADER_EXPORT bool presenceOfMultiDiscPerGeoType() const;
    MEDLOADER_EXPORT void setInfo(const std::vector<std::string>& info);
    MEDLOADER_EXPORT int getNumberOfComponents() const;
    MEDLOADER_EXPORT int getNonEmptyLevels(int iteration, int order, const std::string& mname, std::vector<int>& levs) const;
    MEDLOADER_EXPORT std::vector< std::vector<TypeOfField> > getTypesOfFieldAvailable() const;
    MEDLOADER_EXPORT std::vector< std::vector< std::pair<int,int> > > getFieldSplitedByType(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    MEDLOADER_EXPORT MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> getContent();
  public:
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsed() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsed() const;
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsedMulti() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsedMulti() const;
    MEDLOADER_EXPORT void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
  protected:
    MEDFileAnyTypeFieldMultiTSWithoutSDA *contentNotNullBase();
    const MEDFileAnyTypeFieldMultiTSWithoutSDA *contentNotNullBase() const;
  private:
    static std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > SplitPerCommonSupportNotNodesAlg(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS, const MEDFileMesh *mesh, std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFastCellSupportComparator> >& cmps);
  protected:
    MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> _content;
  };

  class MEDFileIntFieldMultiTS;

  /*!
   * User class.
   */
  class MEDFileFieldMultiTS : public MEDFileAnyTypeFieldMultiTS
  {
  public:
    MEDLOADER_EXPORT static MEDFileFieldMultiTS *New();
    MEDLOADER_EXPORT static MEDFileFieldMultiTS *New(const std::string& fileName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileFieldMultiTS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileFieldMultiTS *New(const MEDFileFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent);
    MEDLOADER_EXPORT static MEDFileFieldMultiTS *LoadSpecificEntities(const std::string& fileName, const std::string& fieldName, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >& entities, bool loadAll=true);
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *shallowCpy() const;
    MEDLOADER_EXPORT void checkCoherencyOfType(const MEDFileAnyTypeField1TS *f1ts) const;
    MEDLOADER_EXPORT MEDFileIntFieldMultiTS *convertToInt(bool isDeepCpyGlobs=true) const;
    //
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *getTimeStepAtPos(int pos) const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *getTimeStep(int iteration, int order) const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *getTimeStepGivenTime(double time, double eps=1e-8) const;
    //
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, const std::string& mname, int iteration, int order, int meshDimRelToMax, int renumPol=0) const;
    MEDLOADER_EXPORT DataArrayDouble *getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const;
    //
    MEDLOADER_EXPORT void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field);
    MEDLOADER_EXPORT void appendFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    MEDLOADER_EXPORT std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    MEDLOADER_EXPORT DataArrayDouble *getUndergroundDataArray(int iteration, int order) const;
    MEDLOADER_EXPORT DataArrayDouble *getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
  private:
    const MEDFileFieldMultiTSWithoutSDA *contentNotNull() const;
    MEDFileFieldMultiTSWithoutSDA *contentNotNull();
  private:
    ~MEDFileFieldMultiTS() { }
    MEDFileFieldMultiTS();
    MEDFileFieldMultiTS(const MEDFileFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent);
    MEDFileFieldMultiTS(const std::string& fileName, bool loadAll, const MEDFileMeshes *ms);
    MEDFileFieldMultiTS(const std::string& fileName, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities=0);
  };

  /*!
   * User class.
   */
  class MEDFileIntFieldMultiTS : public MEDFileAnyTypeFieldMultiTS
  {
  public:
    MEDLOADER_EXPORT static MEDFileIntFieldMultiTS *New();
    MEDLOADER_EXPORT static MEDFileIntFieldMultiTS *New(const std::string& fileName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileIntFieldMultiTS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileIntFieldMultiTS *New(const MEDFileIntFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent);
    MEDLOADER_EXPORT static MEDFileIntFieldMultiTS *LoadSpecificEntities(const std::string& fileName, const std::string& fieldName, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >& entities, bool loadAll=true);
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *shallowCpy() const;
    MEDLOADER_EXPORT void checkCoherencyOfType(const MEDFileAnyTypeField1TS *f1ts) const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *getTimeStepAtPos(int pos) const;
    MEDLOADER_EXPORT MEDFileFieldMultiTS *convertToDouble(bool isDeepCpyGlobs=true) const;
    //
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, DataArrayInt* &arrOut, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int iteration, int order, DataArrayInt* &arrOut, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt* &arrOut, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, DataArrayInt* &arrOut, int renumPol=0) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, int iteration, int order, const std::string& mname, int meshDimRelToMax, DataArrayInt* &arrOut, int renumPol=0) const;
    MEDLOADER_EXPORT DataArrayInt *getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const;
    //
    MEDLOADER_EXPORT void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals);
    MEDLOADER_EXPORT void appendFieldProfile(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    //
    MEDLOADER_EXPORT DataArrayInt *getUndergroundDataArray(int iteration, int order) const;
  private:
    const MEDFileIntFieldMultiTSWithoutSDA *contentNotNull() const;
    MEDFileIntFieldMultiTSWithoutSDA *contentNotNull();
  private:
    ~MEDFileIntFieldMultiTS() { }
    MEDFileIntFieldMultiTS();
    MEDFileIntFieldMultiTS(const MEDFileIntFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent);
    MEDFileIntFieldMultiTS(const std::string& fileName, bool loadAll, const MEDFileMeshes *ms);
    MEDFileIntFieldMultiTS(const std::string& fileName, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities=0);
  };

  class MEDFileAnyTypeFieldMultiTSIterator
  {
  public:
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSIterator(MEDFileAnyTypeFieldMultiTS *fmts);
    MEDLOADER_EXPORT ~MEDFileAnyTypeFieldMultiTSIterator();
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *nextt();
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTS> _fmts;
    int _iter_id;
    int _nb_iter;
  };

  class MEDFileFieldsIterator;

  /*!
   * Use class.
   */
  class MEDFileFields : public RefCountObject, public MEDFileFieldGlobsReal, public MEDFileWritable
  {
  public:
    MEDLOADER_EXPORT static MEDFileFields *New();
    MEDLOADER_EXPORT static MEDFileFields *New(const std::string& fileName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileFields *LoadPartOf(const std::string& fileName, bool loadAll=true, const MEDFileMeshes *ms=0);
    MEDLOADER_EXPORT static MEDFileFields *LoadSpecificEntities(const std::string& fileName, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >& entities, bool loadAll=true);
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT MEDFileFields *deepCpy() const;
    MEDLOADER_EXPORT MEDFileFields *shallowCpy() const;
    MEDLOADER_EXPORT void write(const std::string& fileName, int mode) const;
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    MEDLOADER_EXPORT void loadArrays();
    MEDLOADER_EXPORT void loadArraysIfNecessary();
    MEDLOADER_EXPORT void unloadArrays();
    MEDLOADER_EXPORT void unloadArraysWithoutDataLoss();
    MEDLOADER_EXPORT int getNumberOfFields() const;
    MEDLOADER_EXPORT std::vector< std::pair<int,int> > getCommonIterations(bool& areThereSomeForgottenTS) const;
    MEDLOADER_EXPORT std::vector<std::string> getFieldsNames() const;
    MEDLOADER_EXPORT std::vector<std::string> getMeshesNames() const;
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT void simpleRepr(int bkOffset, std::ostream& oss) const;
    //
    MEDLOADER_EXPORT void resize(int newSize);
    MEDLOADER_EXPORT void pushField(MEDFileAnyTypeFieldMultiTS *field);
    MEDLOADER_EXPORT void pushFields(const std::vector<MEDFileAnyTypeFieldMultiTS *>& fields);
    MEDLOADER_EXPORT void setFieldAtPos(int i, MEDFileAnyTypeFieldMultiTS *field);
    MEDLOADER_EXPORT int getPosFromFieldName(const std::string& fieldName) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *getFieldAtPos(int i) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *getFieldWithName(const std::string& fieldName) const;
    MEDLOADER_EXPORT MEDFileFields *buildSubPart(const int *startIds, const int *endIds) const;
    MEDLOADER_EXPORT bool removeFieldsWithoutAnyTimeStep();
    MEDLOADER_EXPORT MEDFileFields *partOfThisLyingOnSpecifiedMeshName(const std::string& meshName) const;
    MEDLOADER_EXPORT MEDFileFields *partOfThisLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const;
    MEDLOADER_EXPORT MEDFileFields *partOfThisNotLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const;
    MEDLOADER_EXPORT MEDFileFieldsIterator *iterator();
    MEDLOADER_EXPORT void destroyFieldAtPos(int i);
    MEDLOADER_EXPORT void destroyFieldsAtPos(const int *startIds, const int *endIds);
    MEDLOADER_EXPORT void destroyFieldsAtPos2(int bg, int end, int step);
    MEDLOADER_EXPORT bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT bool renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N);
  public:
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsed() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsed() const;
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsedMulti() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsedMulti() const;
    MEDLOADER_EXPORT void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
  private:
    ~MEDFileFields() { }
    MEDFileFields();
    MEDFileFields(const std::string& fileName, bool loadAll, const MEDFileMeshes *ms, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > _fields;
  };

  class MEDFileFieldsIterator
  {
  public:
    MEDLOADER_EXPORT MEDFileFieldsIterator(MEDFileFields *fs);
    MEDLOADER_EXPORT ~MEDFileFieldsIterator();
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *nextt();
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDFileFields> _fs;
    int _iter_id;
    int _nb_iter;
  };
}

#endif
