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

#ifndef __MEDFILEFIELD_HXX__
#define __MEDFILEFIELD_HXX__

#include "MEDLoaderDefines.hxx"

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
    void MEDLOADER_EXPORT simpleRepr(std::ostream& oss) const;
    std::string MEDLOADER_EXPORT getName() const { return _name; }
    void MEDLOADER_EXPORT setName(const char *name);
    static MEDFileFieldLoc *New(med_idt fid, const char *locName);
    static MEDFileFieldLoc *New(med_idt fid, int id);
    static MEDFileFieldLoc *New(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w);
    std::size_t getHeapMemorySize() const;
    MEDFileFieldLoc *deepCpy() const;
    int MEDLOADER_EXPORT getNbOfGaussPtPerCell() const { return _nb_gauss_pt; }
    void MEDLOADER_EXPORT writeLL(med_idt fid) const;
    std::string MEDLOADER_EXPORT repr() const;
    bool MEDLOADER_EXPORT isName(const char *name) const { return _name==name; }
    int MEDLOADER_EXPORT getDimension() const { return _dim; }
    int MEDLOADER_EXPORT getNumberOfGaussPoints() const { return _nb_gauss_pt; }
    int MEDLOADER_EXPORT getNumberOfPointsInCells() const { return _nb_node_per_cell; }
    const MEDLOADER_EXPORT std::vector<double>& getRefCoords() const { return _ref_coo; }
    const MEDLOADER_EXPORT std::vector<double>& getGaussCoords() const { return _gs_coo; }
    const MEDLOADER_EXPORT std::vector<double>& getGaussWeights() const { return _w; }
    bool MEDLOADER_EXPORT isEqual(const MEDFileFieldLoc& other, double eps) const;
  private:
    MEDFileFieldLoc(med_idt fid, const char *locName);
    MEDFileFieldLoc(med_idt fid, int id);
    MEDFileFieldLoc(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w);
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

  class MEDFileFieldPerMeshPerTypePerDisc : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMeshPerTypePerDisc *NewOnRead(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int profileIt) throw(INTERP_KERNEL::Exception);
    static MEDFileFieldPerMeshPerTypePerDisc *New(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int locId);
    static MEDFileFieldPerMeshPerTypePerDisc *New(const MEDFileFieldPerMeshPerTypePerDisc& other);
    std::size_t getHeapMemorySize() const;
    MEDFileFieldPerMeshPerTypePerDisc *deepCpy(MEDFileFieldPerMeshPerType *father) const throw(INTERP_KERNEL::Exception);
    void assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, const DataArray *arrr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, DataArrayInt *locIds, int nbOfEltsInWholeMesh, const MEDCouplingFieldDouble *field, const DataArray *arrr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, const DataArray *arrr, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void getCoarseData(TypeOfField& type, std::pair<int,int>& dad, std::string& pfl, std::string& loc) const throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerType *getFather() const;
    void loadOnlyStructureOfDataRecursively(med_idt fid, int& start, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void loadBigArray(med_idt fid, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void setNewStart(int newValueOfStart) throw(INTERP_KERNEL::Exception);
    int getIteration() const;
    int getOrder() const;
    double getTime() const;
    std::string getMeshName() const;
    TypeOfField getType() const;
    void simpleRepr(int bkOffset, std::ostream& oss, int id) const;
    void fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const throw(INTERP_KERNEL::Exception);
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
    void setProfile(const char *newPflName);
    std::string getLocalization() const;
    void setLocalization(const char *newLocName);
    int getLocId() const { return _loc_id; }
    void setLocId(int newId) const { _loc_id=newId; }
    void setFather(MEDFileFieldPerMeshPerType *newFather) { _father=newFather; }
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void getFieldAtLevel(TypeOfField type, const MEDFileFieldGlobsReal *glob, std::vector< std::pair<int,int> >& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs,
                         std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const;
    void fillValues(int discId, int& startEntryId, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    int fillEltIdsFromCode(int offset, const std::vector<int>& codeOfMesh, const MEDFileFieldGlobsReal& glob, int *ptToFill) const throw(INTERP_KERNEL::Exception);
    int fillTupleIds(int *ptToFill) const throw(INTERP_KERNEL::Exception);
    static int ConvertType(TypeOfField type, int locId) throw(INTERP_KERNEL::Exception);
    static std::vector< std::vector< const MEDFileFieldPerMeshPerTypePerDisc *> > SplitPerDiscretization(const std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>& entries);
    static bool RenumberChunks(int offset, const std::vector< const MEDFileFieldPerMeshPerTypePerDisc *>& entriesOnSameDisc,
                               const DataArrayInt *explicitIdsInMesh, const std::vector<int>& newCode,
                               MEDFileFieldGlobsReal& glob, DataArrayDouble *arr, std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> >& result);
    static MEDFileFieldPerMeshPerTypePerDisc *NewObjectOnSameDiscThanPool(TypeOfField typeF, INTERP_KERNEL::NormalizedCellType geoType, DataArrayInt *idsOfMeshElt,
                                                                          bool isPfl, int nbi, int offset, std::list< const MEDFileFieldPerMeshPerTypePerDisc *>& entriesOnSameDisc,
                                                                          MEDFileFieldGlobsReal& glob, bool &notInExisting) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int profileIt) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int profileIt, const std::string& dummy);
    MEDFileFieldPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc& other);
    MEDFileFieldPerMeshPerTypePerDisc();
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
  public:
    mutable int _tmp_work1;
  };

  class MEDFileFieldPerMeshPerType : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMeshPerType *New(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception);
    static MEDFileFieldPerMeshPerType *NewOnRead(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    MEDFileFieldPerMeshPerType *deepCpy(MEDFileFieldPerMesh *father) const throw(INTERP_KERNEL::Exception);
    void assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, DataArrayInt *locIds, int nbOfEltsInWholeMesh, const MEDCouplingFieldDouble *field, const DataArray *arr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMesh *getFather() const;
    void loadOnlyStructureOfDataRecursively(med_idt fid, int &start, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
    void getDimension(int& dim) const;
    void fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const throw(INTERP_KERNEL::Exception);
    void fillFieldSplitedByType(std::vector< std::pair<int,int> >& dads, std::vector<TypeOfField>& types, std::vector<std::string>& pfls, std::vector<std::string>& locs) const throw(INTERP_KERNEL::Exception);
    int getIteration() const;
    int getOrder() const;
    double getTime() const;
    std::string getMeshName() const;
    void simpleRepr(int bkOffset, std::ostream& oss, int id) const;
    void getSizes(int& globalSz, int& nbOfEntries) const;
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    int getNumberOfComponents() const;
    DataArray *getOrCreateAndGetArray();
    const DataArray *getOrCreateAndGetArray() const;
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenLocId(int locId) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenLocId(int locId) const throw(INTERP_KERNEL::Exception);
    void getFieldAtLevel(int meshDim, TypeOfField type, const MEDFileFieldGlobsReal *glob, std::vector< std::pair<int,int> >& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const;
    void fillValues(int& startEntryId, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    void setLeaves(const std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerTypePerDisc > >& leaves) throw(INTERP_KERNEL::Exception);
    bool keepOnlySpatialDiscretization(TypeOfField tof, int &globalNum, std::vector< std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception);
    static med_entity_type ConvertIntoMEDFileType(TypeOfField ikType, INTERP_KERNEL::NormalizedCellType ikGeoType, med_geometry_type& medfGeoType);
  private:
    std::vector<int> addNewEntryIfNecessary(const MEDCouplingFieldDouble *field, int offset, int nbOfCells) throw(INTERP_KERNEL::Exception);
    std::vector<int> addNewEntryIfNecessaryGauss(const MEDCouplingFieldDouble *field, int offset, int nbOfCells) throw(INTERP_KERNEL::Exception);
    std::vector<int> addNewEntryIfNecessary(const MEDCouplingFieldDouble *field, const DataArrayInt *subCells) throw(INTERP_KERNEL::Exception);
    std::vector<int> addNewEntryIfNecessaryGauss(const MEDCouplingFieldDouble *field, const DataArrayInt *subCells) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerType(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerType(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFieldPerMesh *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> > _field_pm_pt_pd;
    INTERP_KERNEL::NormalizedCellType _geo_type;
  };

  class MEDFileFieldPerMesh : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMesh *New(MEDFileAnyTypeField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh);
    static MEDFileFieldPerMesh *NewOnRead(med_idt fid, MEDFileAnyTypeField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    MEDFileFieldPerMesh *deepCpy(MEDFileAnyTypeField1TSWithoutSDA *father) const throw(INTERP_KERNEL::Exception);
    void simpleRepr(int bkOffset,std::ostream& oss, int id) const;
    void copyTinyInfoFrom(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    void assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const std::vector<int>& code, const std::vector<int>& code2, const std::vector<DataArrayInt *>& idsInPflPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, const DataArray *arr, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void assignFieldNoProfileNoRenum(int& start, const std::vector<int>& code, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void loadOnlyStructureOfDataRecursively(med_idt fid, int &start, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
    void fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector< std::pair<int,int> > > getFieldSplitedByType(std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    void getDimension(int& dim) const;
    double getTime() const;
    int getIteration() const;
    int getOrder() const;
    int getMeshIteration() const { return _mesh_iteration; }
    int getMeshOrder() const { return _mesh_order; }
    std::string getMeshName() const { return _mesh_name; }
    int getNumberOfComponents() const;
    DataArray *getOrCreateAndGetArray();
    const DataArray *getOrCreateAndGetArray() const;
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    bool renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void keepOnlySpatialDiscretization(TypeOfField tof, int &globalNum, std::vector< std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception);
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, bool& isPfl, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
    DataArray *getFieldOnMeshAtLevelWithPfl(TypeOfField type, const MEDCouplingMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
    void getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception);
  private:
    int addNewEntryIfNecessary(INTERP_KERNEL::NormalizedCellType type);
    MEDCouplingFieldDouble *finishField(TypeOfField type, const MEDFileFieldGlobsReal *glob,
                                        const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs, const MEDCouplingMesh *mesh, bool& isPfl, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *finishField2(TypeOfField type, const MEDFileFieldGlobsReal *glob,
                                         const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs,
                                         const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes,
                                         const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *finishFieldNode2(const MEDFileFieldGlobsReal *glob,
                                             const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs,
                                             const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
    DataArray *finishField4(const std::vector< std::pair<int,int> >& dads, const DataArrayInt *pflIn, int nbOfElems, DataArrayInt *&pflOut) const throw(INTERP_KERNEL::Exception);
    void assignNewLeaves(const std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerTypePerDisc > >& leaves) throw(INTERP_KERNEL::Exception);
    static void SortArraysPerType(const MEDFileFieldGlobsReal *glob, TypeOfField type, 
                                  const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& locs,
                                  std::vector<int>& code, std::vector<DataArrayInt *>& notNullPfls);
    static int ComputeNbOfElems(const MEDFileFieldGlobsReal *glob, TypeOfField type, const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMesh(med_idt fid, MEDFileAnyTypeField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMesh(MEDFileAnyTypeField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh);
  private:
    std::string _mesh_name;
    int _mesh_iteration;
    int _mesh_order;
    int _mesh_csit;
    MEDFileAnyTypeField1TSWithoutSDA *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > > _field_pm_pt;
  };

  class MEDFileFieldGlobsReal;

  class MEDFileFieldGlobs : public RefCountObject
  {
  public:
    static MEDFileFieldGlobs *New(const char *fname);
    static MEDFileFieldGlobs *New();
    std::size_t getHeapMemorySize() const;
    MEDFileFieldGlobs *deepCpy() const throw(INTERP_KERNEL::Exception);
    MEDFileFieldGlobs *shallowCpyPart(const std::vector<std::string>& pfls, const std::vector<std::string>& locs) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldGlobs *deepCpyPart(const std::vector<std::string>& pfls, const std::vector<std::string>& locs) const throw(INTERP_KERNEL::Exception);
    void simpleRepr(std::ostream& oss) const;
    void appendGlobs(const MEDFileFieldGlobs& other, double eps) throw(INTERP_KERNEL::Exception);
    void checkGlobsPflsPartCoherency(const std::vector<std::string>& pflsUsed) const throw(INTERP_KERNEL::Exception);
    void checkGlobsLocsPartCoherency(const std::vector<std::string>& locsUsed) const throw(INTERP_KERNEL::Exception);
    void loadProfileInFile(med_idt fid, int id, const char *pflName) throw(INTERP_KERNEL::Exception);
    void loadProfileInFile(med_idt fid, int id);
    void loadGlobals(med_idt fid, const MEDFileFieldGlobsReal& real) throw(INTERP_KERNEL::Exception);
    void loadAllGlobals(med_idt fid) throw(INTERP_KERNEL::Exception);
    void writeGlobals(med_idt fid, const MEDFileWritable& opt) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPfls() const;
    std::vector<std::string> getLocs() const;
    bool existsPfl(const char *pflName) const;
    bool existsLoc(const char *locName) const;
    std::string createNewNameOfPfl() const throw(INTERP_KERNEL::Exception);
    std::string createNewNameOfLoc() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<int> > whichAreEqualProfiles() const;
    std::vector< std::vector<int> > whichAreEqualLocs(double eps) const;
    void setFileName(const char *fileName);
    void changePflsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    int getNbOfGaussPtPerCell(int locId) const throw(INTERP_KERNEL::Exception);
    int getLocalizationId(const char *loc) const throw(INTERP_KERNEL::Exception);
    const char *getFileName() const { return _file_name.c_str(); }
    std::string getFileName2() const { return _file_name; }
    const MEDFileFieldLoc& getLocalizationFromId(int locId) const throw(INTERP_KERNEL::Exception);
    const MEDFileFieldLoc& getLocalization(const char *locName) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getProfileFromId(int pflId) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getProfile(const char *pflName) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldLoc& getLocalizationFromId(int locId) throw(INTERP_KERNEL::Exception);
    MEDFileFieldLoc& getLocalization(const char *locName) throw(INTERP_KERNEL::Exception);
    DataArrayInt *getProfile(const char *pflName) throw(INTERP_KERNEL::Exception);
    DataArrayInt *getProfileFromId(int pflId) throw(INTERP_KERNEL::Exception);
    void killProfileIds(const std::vector<int>& pflIds) throw(INTERP_KERNEL::Exception);
    void killLocalizationIds(const std::vector<int>& locIds) throw(INTERP_KERNEL::Exception);
    //
    void appendProfile(DataArrayInt *pfl) throw(INTERP_KERNEL::Exception);
    void appendLoc(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w) throw(INTERP_KERNEL::Exception);
    //
    static std::string CreateNewNameNotIn(const char *prefix, const std::vector<std::string>& namesToAvoid) throw(INTERP_KERNEL::Exception);
  protected:
    MEDFileFieldGlobs(const char *fname);
    MEDFileFieldGlobs();
    ~MEDFileFieldGlobs();
  protected:
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > _pfls;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> > _locs;
    std::string _file_name;
  };

/// @endcond INTERNAL

  class MEDLOADER_EXPORT MEDFileFieldGlobsReal
  {
  public:
    MEDFileFieldGlobsReal(const char *fname);
    MEDFileFieldGlobsReal();
    std::size_t getHeapMemorySize() const;
    void simpleReprGlobs(std::ostream& oss) const;
    void resetContent();
    void shallowCpyGlobs(const MEDFileFieldGlobsReal& other);
    void deepCpyGlobs(const MEDFileFieldGlobsReal& other);
    void shallowCpyOnlyUsedGlobs(const MEDFileFieldGlobsReal& other) throw(INTERP_KERNEL::Exception);
    void deepCpyOnlyUsedGlobs(const MEDFileFieldGlobsReal& other) throw(INTERP_KERNEL::Exception);
    void appendGlobs(const MEDFileFieldGlobsReal& other, double eps) throw(INTERP_KERNEL::Exception);
    void checkGlobsCoherency() const throw(INTERP_KERNEL::Exception);
    void checkGlobsPflsPartCoherency() const throw(INTERP_KERNEL::Exception);
    void checkGlobsLocsPartCoherency() const throw(INTERP_KERNEL::Exception);
    virtual std::vector<std::string> getPflsReallyUsed() const = 0;
    virtual std::vector<std::string> getLocsReallyUsed() const = 0;
    virtual std::vector<std::string> getPflsReallyUsedMulti() const = 0;
    virtual std::vector<std::string> getLocsReallyUsedMulti() const = 0;
    virtual void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception) = 0;
    virtual void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception) = 0;
    virtual ~MEDFileFieldGlobsReal();
    //
    void loadProfileInFile(med_idt fid, int id, const char *pflName) throw(INTERP_KERNEL::Exception);
    void loadProfileInFile(med_idt fid, int id);
    void loadGlobals(med_idt fid) throw(INTERP_KERNEL::Exception);
    void loadAllGlobals(med_idt fid) throw(INTERP_KERNEL::Exception);
    void writeGlobals(med_idt fid, const MEDFileWritable& opt) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPfls() const;
    std::vector<std::string> getLocs() const;
    bool existsPfl(const char *pflName) const;
    bool existsLoc(const char *locName) const;
    std::string createNewNameOfPfl() const throw(INTERP_KERNEL::Exception);
    std::string createNewNameOfLoc() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<int> > whichAreEqualProfiles() const;
    std::vector< std::vector<int> > whichAreEqualLocs(double eps) const;
    void setFileName(const char *fileName);
    void changePflsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changePflsNames(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsNames(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changePflName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception);
    void changeLocName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception);
    std::vector< std::pair<std::vector<std::string>, std::string > > zipPflsNames() throw(INTERP_KERNEL::Exception);
    std::vector< std::pair<std::vector<std::string>, std::string > > zipLocsNames(double eps) throw(INTERP_KERNEL::Exception);
    int getNbOfGaussPtPerCell(int locId) const throw(INTERP_KERNEL::Exception);
    int getLocalizationId(const char *loc) const throw(INTERP_KERNEL::Exception);
    const char *getFileName() const;
    std::string getFileName2() const;
    const MEDFileFieldLoc& getLocalizationFromId(int locId) const throw(INTERP_KERNEL::Exception);
    const MEDFileFieldLoc& getLocalization(const char *locName) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldLoc& getLocalizationFromId(int locId) throw(INTERP_KERNEL::Exception);
    MEDFileFieldLoc& getLocalization(const char *locName) throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getProfile(const char *pflName) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getProfileFromId(int pflId) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getProfile(const char *pflName) throw(INTERP_KERNEL::Exception);
    DataArrayInt *getProfileFromId(int pflId) throw(INTERP_KERNEL::Exception);
    void killProfileIds(const std::vector<int>& pflIds) throw(INTERP_KERNEL::Exception);
    void killLocalizationIds(const std::vector<int>& locIds) throw(INTERP_KERNEL::Exception);
    //
    void appendProfile(DataArrayInt *pfl) throw(INTERP_KERNEL::Exception);
    void appendLoc(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w) throw(INTERP_KERNEL::Exception);
  protected:
    MEDFileFieldGlobs *contentNotNull() throw(INTERP_KERNEL::Exception);
    const MEDFileFieldGlobs *contentNotNull() const throw(INTERP_KERNEL::Exception);
  protected:
    MEDCouplingAutoRefCountObjectPtr< MEDFileFieldGlobs > _globals;
  };

  class MEDLOADER_EXPORT MEDFileFieldNameScope
  {
  public:
    MEDFileFieldNameScope();
    MEDFileFieldNameScope(const char *fieldName);
    std::string getName() const throw(INTERP_KERNEL::Exception);
    void setName(const char *fieldName) throw(INTERP_KERNEL::Exception);
    std::string getDtUnit() const throw(INTERP_KERNEL::Exception);
    void setDtUnit(const char *dtUnit) throw(INTERP_KERNEL::Exception);
    void copyNameScope(const MEDFileFieldNameScope& other);
  protected:
    std::string _name;
    std::string _dt_unit;
  };

  /*!
   * SDA is for Shared Data Arrays such as profiles.
   */
  class MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA : public RefCountObject, public MEDFileFieldNameScope
  {
  public:
    MEDFileAnyTypeField1TSWithoutSDA();
    MEDFileAnyTypeField1TSWithoutSDA(const char *fieldName, int csit, int iteration, int order);
    int getIteration() const { return _iteration; }
    int getOrder() const { return _order; }
    double getTime(int& iteration, int& order) const { iteration=_iteration; order=_order; return _dt; }
    void setTime(int iteration, int order, double val) { _dt=val; _iteration=iteration; _order=order; }
    int getDimension() const;
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    void setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception);
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    int getMeshIteration() const throw(INTERP_KERNEL::Exception);
    int getMeshOrder() const throw(INTERP_KERNEL::Exception);
    bool isDealingTS(int iteration, int order) const;
    std::pair<int,int> getDtIt() const;
    void fillIteration(std::pair<int,int>& p) const;
    void fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const throw(INTERP_KERNEL::Exception);
    std::vector<TypeOfField> getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception);
    //
    std::vector<std::string> getPflsReallyUsed2() const;
    std::vector<std::string> getLocsReallyUsed2() const;
    std::vector<std::string> getPflsReallyUsedMulti2() const;
    std::vector<std::string> getLocsReallyUsedMulti2() const;
    void changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    //
    int getNonEmptyLevels(const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<std::pair<int,int> > > getFieldSplitedByType(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    //
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception);
     const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception);
     void deepCpyLeavesFrom(const MEDFileAnyTypeField1TSWithoutSDA& other) throw(INTERP_KERNEL::Exception);
  public:
    int getNumberOfComponents() const;
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string>& getInfo();
    void setInfo(const std::vector<std::string>& infos) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    int copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr) throw(INTERP_KERNEL::Exception);
    void setFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void setFieldProfile(const MEDCouplingFieldDouble *field, const DataArray *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    virtual void simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const;
    virtual MEDFileAnyTypeField1TSWithoutSDA *deepCpy() const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDFileAnyTypeField1TSWithoutSDA *shallowCpy() const throw(INTERP_KERNEL::Exception) = 0;
    virtual std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > splitComponents() const throw(INTERP_KERNEL::Exception);
    virtual const char *getTypeStr() const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArray *getUndergroundDataArray() const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArray *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void setArray(DataArray *arr) throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArray *createNewEmptyDataArrayInstance() const = 0;
    virtual DataArray *getOrCreateAndGetArray() = 0;
    virtual const DataArray *getOrCreateAndGetArray() const = 0;
  public:
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, const char *mName, int renumPol, const MEDFileFieldGlobsReal *glob, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDFileMesh *mesh, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, const char *mName, int renumPol, const MEDFileFieldGlobsReal *glob, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, const DataArrayInt *cellRenum, const DataArrayInt *nodeRenum, MEDCouplingAutoRefCountObjectPtr<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
    DataArray *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
  public:
    bool renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > splitDiscretizations() const throw(INTERP_KERNEL::Exception);
    int keepOnlySpatialDiscretization(TypeOfField tof, std::vector< std::pair<int,int> >& its) throw(INTERP_KERNEL::Exception);
  public:
    void allocNotFromFile(int newNbOfTuples) throw(INTERP_KERNEL::Exception);
    bool allocIfNecessaryTheArrayToReceiveDataFromFile() throw(INTERP_KERNEL::Exception);
    void loadOnlyStructureOfDataRecursively(med_idt fid, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void loadBigArraysRecursivelyIfNecessary(med_idt fid, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void loadStructureAndBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void unloadArrays() throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid, const MEDFileWritable& opts, const MEDFileFieldNameScope& nasc) const throw(INTERP_KERNEL::Exception);
  protected:
    int getMeshIdFromMeshName(const char *mName) const throw(INTERP_KERNEL::Exception);
    int addNewEntryIfNecessary(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    void updateData(int newLgth, const std::vector< std::pair<int,int> >& oldStartStops) throw(INTERP_KERNEL::Exception);
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
  class MEDLOADER_EXPORT MEDFileField1TSWithoutSDA : public MEDFileAnyTypeField1TSWithoutSDA
  {
  public:
    const char *getTypeStr() const throw(INTERP_KERNEL::Exception);
    DataArray *getUndergroundDataArray() const throw(INTERP_KERNEL::Exception);
    DataArray *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArrayDouble() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArrayDoubleExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    static void CheckMeshDimRel(int meshDimRelToMax) throw(INTERP_KERNEL::Exception);
    static std::vector<int> CheckSBTMesh(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    static MEDFileField1TSWithoutSDA *New(const char *fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos);
  public:
    MEDFileField1TSWithoutSDA();
    MEDFileField1TSWithoutSDA(const char *fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos);
    MEDFileAnyTypeField1TSWithoutSDA *shallowCpy() const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TSWithoutSDA *deepCpy() const throw(INTERP_KERNEL::Exception);
    void setArray(DataArray *arr) throw(INTERP_KERNEL::Exception);
    DataArray *createNewEmptyDataArrayInstance() const;
    DataArray *getOrCreateAndGetArray();
    const DataArray *getOrCreateAndGetArray() const;
    DataArrayDouble *getOrCreateAndGetArrayDouble();
    const DataArrayDouble *getOrCreateAndGetArrayDouble() const;
    MEDFileIntField1TSWithoutSDA *convertToInt() const throw(INTERP_KERNEL::Exception);
  protected:
    MEDCouplingAutoRefCountObjectPtr< DataArrayDouble > _arr;
  public:
    static const char TYPE_STR[];
  };

  /*!
   * SDA is for Shared Data Arrays such as profiles.
   */
  class MEDLOADER_EXPORT MEDFileIntField1TSWithoutSDA : public MEDFileAnyTypeField1TSWithoutSDA
  {
  public:
    MEDFileIntField1TSWithoutSDA();
    static MEDFileIntField1TSWithoutSDA *New(const char *fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos);
    MEDFileAnyTypeField1TSWithoutSDA *deepCpy() const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TSWithoutSDA *shallowCpy() const throw(INTERP_KERNEL::Exception);
    const char *getTypeStr() const throw(INTERP_KERNEL::Exception);
    DataArray *getUndergroundDataArray() const throw(INTERP_KERNEL::Exception);
    DataArray *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
    void setArray(DataArray *arr) throw(INTERP_KERNEL::Exception);
    DataArray *createNewEmptyDataArrayInstance() const;
    DataArray *getOrCreateAndGetArray();
    const DataArray *getOrCreateAndGetArray() const;
    DataArrayInt *getOrCreateAndGetArrayInt();
    const DataArrayInt *getOrCreateAndGetArrayInt() const;
    DataArrayInt *getUndergroundDataArrayInt() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getUndergroundDataArrayIntExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
    MEDFileField1TSWithoutSDA *convertToDouble() const throw(INTERP_KERNEL::Exception);
  protected:
    MEDFileIntField1TSWithoutSDA(const char *fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos);
  protected:
    MEDCouplingAutoRefCountObjectPtr< DataArrayInt > _arr;
  public:
    static const char TYPE_STR[];
  };

  /*!
   * User class.
   */
  class MEDLOADER_EXPORT MEDFileAnyTypeField1TS : public RefCountObject, public MEDFileWritable, public MEDFileFieldGlobsReal
  {
  protected:
    MEDFileAnyTypeField1TS();
    MEDFileAnyTypeField1TS(const char *fileName, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TS(const char *fileName, const char *fieldName, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TS(const char *fileName, const char *fieldName, int iteration, int order, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TS(const MEDFileAnyTypeField1TSWithoutSDA& other, bool shallowCopyOfContent);
    static MEDFileAnyTypeField1TS *BuildNewInstanceFromContent(MEDFileAnyTypeField1TSWithoutSDA *c, const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileAnyTypeField1TSWithoutSDA *BuildContentFrom(med_idt fid, const char *fileName, bool loadAll) throw(INTERP_KERNEL::Exception);
    static MEDFileAnyTypeField1TSWithoutSDA *BuildContentFrom(med_idt fid, const char *fileName, const char *fieldName, bool loadAll) throw(INTERP_KERNEL::Exception);
    static MEDFileAnyTypeField1TSWithoutSDA *BuildContentFrom(med_idt fid, const char *fileName, const char *fieldName, int iteration, int order, bool loadAll) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    // direct forwarding to MEDFileAnyTypeField1TSWithoutSDA instance _content
  public:
    static MEDFileAnyTypeField1TS *New(const char *fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileAnyTypeField1TS *New(const char *fileName, const char *fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileAnyTypeField1TS *New(const char *fileName, const char *fieldName, int iteration, int order, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    int getDimension() const;
    int getIteration() const;
    int getOrder() const;
    double getTime(int& iteration, int& order) const;
    void setTime(int iteration, int order, double val);
    std::string getName() const;
    void setName(const char *name);
    std::string simpleRepr() const;
    void simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const;
    std::string getDtUnit() const throw(INTERP_KERNEL::Exception);
    void setDtUnit(const char *dtUnit) throw(INTERP_KERNEL::Exception);
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    void setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception);
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    int getMeshIteration() const throw(INTERP_KERNEL::Exception);
    int getMeshOrder() const throw(INTERP_KERNEL::Exception);
    int getNumberOfComponents() const;
    bool isDealingTS(int iteration, int order) const;
    std::pair<int,int> getDtIt() const;
    void fillIteration(std::pair<int,int>& p) const;
    void fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const throw(INTERP_KERNEL::Exception);
    void setInfo(const std::vector<std::string>& infos) throw(INTERP_KERNEL::Exception);
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string>& getInfo();
    std::vector<TypeOfField> getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<std::pair<int,int> > > getFieldSplitedByType(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
                                                                          std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception);
    int getNonEmptyLevels(const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception);
  public:
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void loadArrays() throw(INTERP_KERNEL::Exception);
    void loadArraysIfNecessary() throw(INTERP_KERNEL::Exception);
    void unloadArrays() throw(INTERP_KERNEL::Exception);
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TS > > splitComponents() const throw(INTERP_KERNEL::Exception);
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeField1TS > > splitDiscretizations() const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TS *deepCpy() const throw(INTERP_KERNEL::Exception);
    int copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr) throw(INTERP_KERNEL::Exception);
    virtual MEDFileAnyTypeField1TS *shallowCpy() const throw(INTERP_KERNEL::Exception) = 0;
  public:
    //! underground method see MEDFileField1TSWithoutSDA::setProfileNameOnLeaf
    void setProfileNameOnLeaf(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newPflName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception);
    //! underground method see MEDFileField1TSWithoutSDA::setLocNameOnLeaf
    void setLocNameOnLeaf(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newLocName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
  public:
    static int LocateField2(med_idt fid, const char *fileName, int fieldIdCFormat, bool checkFieldId, std::string& fieldName, med_field_type& typcha, std::vector<std::string>& infos, std::string& dtunitOut) throw(INTERP_KERNEL::Exception);
    static int LocateField(med_idt fid, const char *fileName, const char *fieldName, int& posCFormat, med_field_type& typcha, std::vector<std::string>& infos, std::string& dtunitOut) throw(INTERP_KERNEL::Exception);
  public:
    virtual med_field_type getMEDFileFieldType() const = 0;
    MEDFileAnyTypeField1TSWithoutSDA *contentNotNullBase() throw(INTERP_KERNEL::Exception);
    const MEDFileAnyTypeField1TSWithoutSDA *contentNotNullBase() const throw(INTERP_KERNEL::Exception);
  protected:
    MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> _content;
  };
  
  class MEDFileIntField1TS;
  
  /*!
   * User class.
   */
  class MEDLOADER_EXPORT MEDFileField1TS : public MEDFileAnyTypeField1TS
  {
  public:
    static MEDFileField1TS *New(const char *fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileField1TS *New(const char *fileName, const char *fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileField1TS *New(const char *fileName, const char *fieldName, int iteration, int order, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileField1TS *New(const MEDFileField1TSWithoutSDA& other, bool shallowCopyOfContent);
    static MEDFileField1TS *New();
    MEDFileIntField1TS *convertToInt(bool deepCpyGlobs=true) const throw(INTERP_KERNEL::Exception);
    //
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, const char *mname, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception);
    //
    void setFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    void setFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception);
    // direct forwarding to MEDFileField1TSWithoutSDA instance _content
  public:
    MEDFileAnyTypeField1TS *shallowCpy() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArray() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
    
    std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
                                                                         std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
  public:
    static void SetDataArrayDoubleInField(MEDCouplingFieldDouble *f, MEDCouplingAutoRefCountObjectPtr<DataArray>& arr) throw(INTERP_KERNEL::Exception);
    static DataArrayDouble *ReturnSafelyDataArrayDouble(MEDCouplingAutoRefCountObjectPtr<DataArray>& arr) throw(INTERP_KERNEL::Exception);
  private:
    med_field_type getMEDFileFieldType() const { return MED_FLOAT64; }
    const MEDFileField1TSWithoutSDA *contentNotNull() const throw(INTERP_KERNEL::Exception);
    MEDFileField1TSWithoutSDA *contentNotNull() throw(INTERP_KERNEL::Exception);
  private:
    MEDFileField1TS(const char *fileName, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileField1TS(const char *fileName, const char *fieldName, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileField1TS(const MEDFileField1TSWithoutSDA& other, bool shallowCopyOfContent);
    MEDFileField1TS();
  };

  class MEDLOADER_EXPORT MEDFileIntField1TS : public MEDFileAnyTypeField1TS
  {
  public:
    static MEDFileIntField1TS *New();
    static MEDFileIntField1TS *New(const char *fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileIntField1TS *New(const char *fileName, const char *fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileIntField1TS *New(const char *fileName, const char *fieldName, int iteration, int order, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileIntField1TS *New(const MEDFileIntField1TSWithoutSDA& other, bool shallowCopyOfContent);
    MEDFileField1TS *convertToDouble(bool deepCpyGlobs=true) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TS *shallowCpy() const throw(INTERP_KERNEL::Exception);
    //
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, DataArrayInt* &arrOut, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, DataArrayInt* &arrOut, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt* &arrOut, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, DataArrayInt* &arrOut, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, const char *mname, int meshDimRelToMax, DataArrayInt* &arrOut, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception);
    //
    void setFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals) throw(INTERP_KERNEL::Exception);
    void setFieldProfile(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception);
    DataArrayInt *getUndergroundDataArray() const throw(INTERP_KERNEL::Exception);
  public:
    static DataArrayInt *ReturnSafelyDataArrayInt(MEDCouplingAutoRefCountObjectPtr<DataArray>& arr) throw(INTERP_KERNEL::Exception);
  private:
    med_field_type getMEDFileFieldType() const { return MED_INT32; }
    const MEDFileIntField1TSWithoutSDA *contentNotNull() const throw(INTERP_KERNEL::Exception);
    MEDFileIntField1TSWithoutSDA *contentNotNull() throw(INTERP_KERNEL::Exception);
  private:
    MEDFileIntField1TS();
    MEDFileIntField1TS(const char *fileName, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileIntField1TS(const char *fileName, const char *fieldName, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileIntField1TS(const char *fileName, const char *fieldName, int iteration, int order, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileIntField1TS(const MEDFileIntField1TSWithoutSDA& other, bool shallowCopyOfContent);
  };
  
  class MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA : public RefCountObject, public MEDFileFieldNameScope
  {
  protected:
    MEDFileAnyTypeFieldMultiTSWithoutSDA();
    MEDFileAnyTypeFieldMultiTSWithoutSDA(const char *fieldName);
    MEDFileAnyTypeFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTSWithoutSDA(med_idt fid, const char *fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll) throw(INTERP_KERNEL::Exception);
  public:
    std::size_t getHeapMemorySize() const;
    virtual MEDFileAnyTypeFieldMultiTSWithoutSDA *deepCpy() const throw(INTERP_KERNEL::Exception);
    virtual std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > splitComponents() const throw(INTERP_KERNEL::Exception);
    virtual std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > splitDiscretizations() const throw(INTERP_KERNEL::Exception);
    virtual const char *getTypeStr() const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDFileAnyTypeFieldMultiTSWithoutSDA *shallowCpy() const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDFileAnyTypeFieldMultiTSWithoutSDA *createNew() const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDFileAnyTypeField1TSWithoutSDA *createNew1TSWithoutSDAEmptyInstance() const throw(INTERP_KERNEL::Exception) = 0;
    virtual void checkCoherencyOfType(const MEDFileAnyTypeField1TSWithoutSDA *f1ts) const throw(INTERP_KERNEL::Exception) = 0;
    const std::vector<std::string>& getInfo() const throw(INTERP_KERNEL::Exception);
    void setInfo(const std::vector<std::string>& info) throw(INTERP_KERNEL::Exception);
    int getTimeStepPos(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    const MEDFileAnyTypeField1TSWithoutSDA& getTimeStepEntry(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TSWithoutSDA& getTimeStepEntry(int iteration, int order) throw(INTERP_KERNEL::Exception);
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    void setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception);
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    int getNumberOfTS() const;
    void eraseEmptyTS() throw(INTERP_KERNEL::Exception);
    void eraseTimeStepIds(const int *startIds, const int *endIds) throw(INTERP_KERNEL::Exception);
    void eraseTimeStepIds2(int bg, int end, int step) throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTSWithoutSDA *buildFromTimeStepIds(const int *startIds, const int *endIds) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTSWithoutSDA *buildFromTimeStepIds2(int bg, int end, int step) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTSWithoutSDA *partOfThisLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTSWithoutSDA *partOfThisNotLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const throw(INTERP_KERNEL::Exception);
    int getPosOfTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    int getPosGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    std::vector< std::pair<int,int> > getIterations() const;
    std::vector< std::pair<int,int> > getTimeSteps(std::vector<double>& ret1) const throw(INTERP_KERNEL::Exception);
    void pushBackTimeStep(MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA>& tse) throw(INTERP_KERNEL::Exception);
    void synchronizeNameScope() throw(INTERP_KERNEL::Exception);
    void simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const;
    int getNonEmptyLevels(int iteration, int order, const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception);
    void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void appendFieldProfile(const MEDCouplingFieldDouble *field, const DataArray *arr, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    std::vector< std::vector< std::pair<int,int> > > getFieldSplitedByType(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<TypeOfField> > getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception);
    DataArray *getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    DataArray *getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
    bool renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void loadStructureOrStructureAndBigArraysRecursively(med_idt fid, int nbPdt, med_field_type fieldTyp, bool loadAll) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid, const MEDFileWritable& opts) const throw(INTERP_KERNEL::Exception);
    void loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void loadBigArraysRecursivelyIfNecessary(med_idt fid, const MEDFileFieldNameScope& nasc) throw(INTERP_KERNEL::Exception);
    void unloadArrays() throw(INTERP_KERNEL::Exception);
  public:
    const MEDFileAnyTypeField1TSWithoutSDA *getTimeStepAtPos2(int pos) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TSWithoutSDA *getTimeStepAtPos2(int pos) throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPflsReallyUsed2() const;
    std::vector<std::string> getLocsReallyUsed2() const;
    std::vector<std::string> getPflsReallyUsedMulti2() const;
    std::vector<std::string> getLocsReallyUsedMulti2() const;
    void changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void setIteration(int i, MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> ts) throw(INTERP_KERNEL::Exception);
  protected:
    virtual med_field_type getMEDFileFieldType() const = 0;
    void copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr) throw(INTERP_KERNEL::Exception);
    void checkCoherencyOfTinyInfo(const MEDCouplingFieldDouble *field, const DataArray *arr) const throw(INTERP_KERNEL::Exception);
    void checkThatComponentsMatch(const std::vector<std::string>& compos) const throw(INTERP_KERNEL::Exception);
    void checkThatNbOfCompoOfTSMatchThis() const throw(INTERP_KERNEL::Exception);
  protected:
    std::vector<std::string> _infos;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeField1TSWithoutSDA> > _time_steps;
  };

  class MEDFileIntFieldMultiTSWithoutSDA;

  class MEDLOADER_EXPORT MEDFileFieldMultiTSWithoutSDA : public MEDFileAnyTypeFieldMultiTSWithoutSDA
  {
  public:
    static MEDFileFieldMultiTSWithoutSDA *New(med_idt fid, const char *fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll) throw(INTERP_KERNEL::Exception);
    const char *getTypeStr() const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTSWithoutSDA *shallowCpy() const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTSWithoutSDA *createNew() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    MEDFileIntFieldMultiTSWithoutSDA *convertToInt() const throw(INTERP_KERNEL::Exception);
  protected:
    MEDFileFieldMultiTSWithoutSDA(const char *fieldName);
    MEDFileFieldMultiTSWithoutSDA(med_idt fid, const char *fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll) throw(INTERP_KERNEL::Exception);
    med_field_type getMEDFileFieldType() const { return MED_FLOAT64; }
    MEDFileAnyTypeField1TSWithoutSDA *createNew1TSWithoutSDAEmptyInstance() const throw(INTERP_KERNEL::Exception);
    void checkCoherencyOfType(const MEDFileAnyTypeField1TSWithoutSDA *f1ts) const throw(INTERP_KERNEL::Exception);
  public:
    MEDFileFieldMultiTSWithoutSDA();
  };

  class MEDLOADER_EXPORT MEDFileIntFieldMultiTSWithoutSDA : public MEDFileAnyTypeFieldMultiTSWithoutSDA
  {
  public:
    static MEDFileIntFieldMultiTSWithoutSDA *New(med_idt fid, const char *fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileIntFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll) throw(INTERP_KERNEL::Exception);
    const char *getTypeStr() const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTSWithoutSDA *shallowCpy() const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTSWithoutSDA *createNew() const throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTSWithoutSDA *convertToDouble() const throw(INTERP_KERNEL::Exception);
  protected:
    MEDFileIntFieldMultiTSWithoutSDA(const char *fieldName);
    MEDFileIntFieldMultiTSWithoutSDA(med_idt fid, const char *fieldName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll) throw(INTERP_KERNEL::Exception);
    med_field_type getMEDFileFieldType() const { return MED_INT32; }
    MEDFileAnyTypeField1TSWithoutSDA *createNew1TSWithoutSDAEmptyInstance() const throw(INTERP_KERNEL::Exception);
    void checkCoherencyOfType(const MEDFileAnyTypeField1TSWithoutSDA *f1ts) const throw(INTERP_KERNEL::Exception);
  public:
    MEDFileIntFieldMultiTSWithoutSDA();
  };

  class MEDFileAnyTypeFieldMultiTSIterator;
  class MEDFileFastCellSupportComparator;
  /*!
   * User class.
   */
  class MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS : public RefCountObject, public MEDFileWritable, public MEDFileFieldGlobsReal
  {
  protected:
    MEDFileAnyTypeFieldMultiTS();
    MEDFileAnyTypeFieldMultiTS(const char *fileName, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTS(const char *fileName, const char *fieldName, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTS(const MEDFileAnyTypeFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent);
    static MEDFileAnyTypeFieldMultiTS *BuildNewInstanceFromContent(MEDFileAnyTypeFieldMultiTSWithoutSDA *c, const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileAnyTypeFieldMultiTSWithoutSDA *BuildContentFrom(med_idt fid, const char *fileName, bool loadAll) throw(INTERP_KERNEL::Exception);
    static MEDFileAnyTypeFieldMultiTSWithoutSDA *BuildContentFrom(med_idt fid, const char *fileName, const char *fieldName, bool loadAll) throw(INTERP_KERNEL::Exception);
  public:
    static MEDFileAnyTypeFieldMultiTS *New(const char *fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileAnyTypeFieldMultiTS *New(const char *fileName, const char *fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    void loadArrays() throw(INTERP_KERNEL::Exception);
    void loadArraysIfNecessary() throw(INTERP_KERNEL::Exception);
    void unloadArrays() throw(INTERP_KERNEL::Exception);
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    virtual MEDFileAnyTypeFieldMultiTS *deepCpy() const throw(INTERP_KERNEL::Exception);
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTS > > splitComponents() const throw(INTERP_KERNEL::Exception);
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileAnyTypeFieldMultiTS > > splitDiscretizations() const throw(INTERP_KERNEL::Exception);
    virtual MEDFileAnyTypeFieldMultiTS *shallowCpy() const throw(INTERP_KERNEL::Exception) = 0;
    virtual void checkCoherencyOfType(const MEDFileAnyTypeField1TS *f1ts) const throw(INTERP_KERNEL::Exception) = 0;
    //
    virtual MEDFileAnyTypeField1TS *getTimeStepAtPos(int pos) const throw(INTERP_KERNEL::Exception) = 0;
    MEDFileAnyTypeField1TS *getTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TS *getTimeStepGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    static std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > SplitIntoCommonTimeSeries(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS) throw(INTERP_KERNEL::Exception);
    static std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > SplitPerCommonSupport(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS, const MEDFileMesh *mesh, std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFastCellSupportComparator> >& fsc) throw(INTERP_KERNEL::Exception);
    static int CheckSupportAcrossTime(MEDFileAnyTypeFieldMultiTS *f0, MEDFileAnyTypeFieldMultiTS *f1, const MEDFileMesh *mesh, TypeOfField& tof0, TypeOfField& tof1) throw(INTERP_KERNEL::Exception);
  public:// direct forwarding to MEDFileField1TSWithoutSDA instance _content
    std::string getName() const;
    void setName(const char *name);
    std::string getDtUnit() const throw(INTERP_KERNEL::Exception);
    void setDtUnit(const char *dtUnit) throw(INTERP_KERNEL::Exception);
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    void setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    void simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const;
    int getNumberOfTS() const;
    void eraseEmptyTS() throw(INTERP_KERNEL::Exception);
    void eraseTimeStepIds(const int *startIds, const int *endIds) throw(INTERP_KERNEL::Exception);
    void eraseTimeStepIds2(int bg, int end, int step) throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTS *buildSubPart(const int *startIds, const int *endIds) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTS *buildSubPartSlice(int bg, int end, int step) const throw(INTERP_KERNEL::Exception);
    std::vector< std::pair<int,int> > getTimeSteps(std::vector<double>& ret1) const throw(INTERP_KERNEL::Exception);
    std::vector< std::pair<int,int> > getIterations() const;
    void pushBackTimeSteps(const std::vector<MEDFileAnyTypeField1TS *>& f1ts) throw(INTERP_KERNEL::Exception);
    void pushBackTimeStep(MEDFileAnyTypeField1TS *f1ts) throw(INTERP_KERNEL::Exception);
    void synchronizeNameScope() throw(INTERP_KERNEL::Exception);
    int getPosOfTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    int getPosGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTSIterator *iterator() throw(INTERP_KERNEL::Exception);
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    const std::vector<std::string>& getInfo() const throw(INTERP_KERNEL::Exception);
    void setInfo(const std::vector<std::string>& info) throw(INTERP_KERNEL::Exception);
    int getNumberOfComponents() const throw(INTERP_KERNEL::Exception);
    int getNonEmptyLevels(int iteration, int order, const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<TypeOfField> > getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector< std::pair<int,int> > > getFieldSplitedByType(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> getContent();
  public:
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
  protected:
    MEDFileAnyTypeFieldMultiTSWithoutSDA *contentNotNullBase() throw(INTERP_KERNEL::Exception);
    const MEDFileAnyTypeFieldMultiTSWithoutSDA *contentNotNullBase() const throw(INTERP_KERNEL::Exception);
  private:
    static std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > SplitPerCommonSupportNotNodesAlg(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS, const MEDFileMesh *mesh, std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFastCellSupportComparator> >& cmps) throw(INTERP_KERNEL::Exception);
  protected:
    MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> _content;
  };

  class MEDFileIntFieldMultiTS;

  /*!
   * User class.
   */
  class MEDLOADER_EXPORT MEDFileFieldMultiTS : public MEDFileAnyTypeFieldMultiTS
  {
  public:
    static MEDFileFieldMultiTS *New();
    static MEDFileFieldMultiTS *New(const char *fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileFieldMultiTS *New(const char *fileName, const char *fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileFieldMultiTS *New(const MEDFileFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent);
    MEDFileAnyTypeFieldMultiTS *shallowCpy() const throw(INTERP_KERNEL::Exception);
    void checkCoherencyOfType(const MEDFileAnyTypeField1TS *f1ts) const throw(INTERP_KERNEL::Exception);
    MEDFileIntFieldMultiTS *convertToInt(bool deepCpyGlobs=true) const throw(INTERP_KERNEL::Exception);
    //
    MEDFileAnyTypeField1TS *getTimeStepAtPos(int pos) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TS *getTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TS *getTimeStepGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    //
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, const char *mname, int iteration, int order, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception);
    //
    void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    void appendFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception);
std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
  private:
    const MEDFileFieldMultiTSWithoutSDA *contentNotNull() const throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTSWithoutSDA *contentNotNull() throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFieldMultiTS();
    MEDFileFieldMultiTS(const MEDFileFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent);
    MEDFileFieldMultiTS(const char *fileName, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTS(const char *fileName, const char *fieldName, bool loadAll) throw(INTERP_KERNEL::Exception);
  };

  /*!
   * User class.
   */
  class MEDLOADER_EXPORT MEDFileIntFieldMultiTS : public MEDFileAnyTypeFieldMultiTS
  {
  public:
    static MEDFileIntFieldMultiTS *New();
    static MEDFileIntFieldMultiTS *New(const char *fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileIntFieldMultiTS *New(const char *fileName, const char *fieldName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    static MEDFileIntFieldMultiTS *New(const MEDFileIntFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent);
    MEDFileAnyTypeFieldMultiTS *shallowCpy() const throw(INTERP_KERNEL::Exception);
    void checkCoherencyOfType(const MEDFileAnyTypeField1TS *f1ts) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeField1TS *getTimeStepAtPos(int pos) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTS *convertToDouble(bool deepCpyGlobs=true) const throw(INTERP_KERNEL::Exception);
    //
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, DataArrayInt* &arrOut, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int iteration, int order, DataArrayInt* &arrOut, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt* &arrOut, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, DataArrayInt* &arrOut, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, int iteration, int order, const char *mname, int meshDimRelToMax, DataArrayInt* &arrOut, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception);
    //
    void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals) throw(INTERP_KERNEL::Exception);
    void appendFieldProfile(const MEDCouplingFieldDouble *field, const DataArrayInt *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception);
    //
    DataArrayInt *getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception);
  private:
    const MEDFileIntFieldMultiTSWithoutSDA *contentNotNull() const throw(INTERP_KERNEL::Exception);
    MEDFileIntFieldMultiTSWithoutSDA *contentNotNull() throw(INTERP_KERNEL::Exception);
  private:
    MEDFileIntFieldMultiTS();
    MEDFileIntFieldMultiTS(const MEDFileIntFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent);
    MEDFileIntFieldMultiTS(const char *fileName, bool loadAll) throw(INTERP_KERNEL::Exception);
    MEDFileIntFieldMultiTS(const char *fileName, const char *fieldName, bool loadAll) throw(INTERP_KERNEL::Exception);
  };

  class MEDCOUPLING_EXPORT MEDFileAnyTypeFieldMultiTSIterator
  {
  public:
    MEDFileAnyTypeFieldMultiTSIterator(MEDFileAnyTypeFieldMultiTS *fmts);
    ~MEDFileAnyTypeFieldMultiTSIterator();
    MEDFileAnyTypeField1TS *nextt() throw(INTERP_KERNEL::Exception);
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTS> _fmts;
     int _iter_id;
     int _nb_iter;
  };

  class MEDFileFieldsIterator;

  /*!
   * Use class.
   */
  class MEDLOADER_EXPORT MEDFileFields : public RefCountObject, public MEDFileFieldGlobsReal, public MEDFileWritable
  {
  public:
    static MEDFileFields *New();
    static MEDFileFields *New(const char *fileName, bool loadAll=true) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    MEDFileFields *deepCpy() const throw(INTERP_KERNEL::Exception);
    MEDFileFields *shallowCpy() const throw(INTERP_KERNEL::Exception);
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    void loadArrays() throw(INTERP_KERNEL::Exception);
    void loadArraysIfNecessary() throw(INTERP_KERNEL::Exception);
    void unloadArrays() throw(INTERP_KERNEL::Exception);
    int getNumberOfFields() const;
    std::vector< std::pair<int,int> > getCommonIterations(bool& areThereSomeForgottenTS) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getFieldsNames() const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getMeshesNames() const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    void simpleRepr(int bkOffset, std::ostream& oss) const;
    //
    void resize(int newSize) throw(INTERP_KERNEL::Exception);
    void pushField(MEDFileAnyTypeFieldMultiTS *field) throw(INTERP_KERNEL::Exception);
    void pushFields(const std::vector<MEDFileAnyTypeFieldMultiTS *>& fields) throw(INTERP_KERNEL::Exception);
    void setFieldAtPos(int i, MEDFileAnyTypeFieldMultiTS *field) throw(INTERP_KERNEL::Exception);
    int getPosFromFieldName(const char *fieldName) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTS *getFieldAtPos(int i) const throw(INTERP_KERNEL::Exception);
    MEDFileAnyTypeFieldMultiTS *getFieldWithName(const char *fieldName) const throw(INTERP_KERNEL::Exception);
    MEDFileFields *buildSubPart(const int *startIds, const int *endIds) const throw(INTERP_KERNEL::Exception);
    MEDFileFields *partOfThisLyingOnSpecifiedMeshName(const char *meshName) const throw(INTERP_KERNEL::Exception);
    MEDFileFields *partOfThisLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const throw(INTERP_KERNEL::Exception);
    MEDFileFields *partOfThisNotLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldsIterator *iterator() throw(INTERP_KERNEL::Exception);
    void destroyFieldAtPos(int i) throw(INTERP_KERNEL::Exception);
    void destroyFieldsAtPos(const int *startIds, const int *endIds) throw(INTERP_KERNEL::Exception);
    void destroyFieldsAtPos2(int bg, int end, int step) throw(INTERP_KERNEL::Exception);
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    bool renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N) throw(INTERP_KERNEL::Exception);
  public:
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFields();
    MEDFileFields(const char *fileName, bool loadAll) throw(INTERP_KERNEL::Exception);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileAnyTypeFieldMultiTSWithoutSDA> > _fields;
  };

  class MEDCOUPLING_EXPORT MEDFileFieldsIterator
  {
  public:
    MEDFileFieldsIterator(MEDFileFields *fs);
    ~MEDFileFieldsIterator();
    MEDFileAnyTypeFieldMultiTS *nextt();
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDFileFields> _fs;
     int _iter_id;
     int _nb_iter;
  };
}

#endif
