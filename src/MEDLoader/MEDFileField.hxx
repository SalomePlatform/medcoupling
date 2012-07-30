// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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
    const MEDLOADER_EXPORT std::string& getName() const { return _name; }
    void MEDLOADER_EXPORT setName(const char *name);
    static MEDFileFieldLoc *New(med_idt fid, const char *locName);
    static MEDFileFieldLoc *New(med_idt fid, int id);
    static MEDFileFieldLoc *New(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w);
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
  class MEDFileFieldPerMeshPerType;
  class MEDFileField1TSWithoutSDA;
  class MEDFileFieldGlobsReal;
  class MEDFileFieldPerMesh;

  class MEDFileFieldPerMeshPerTypePerDisc : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMeshPerTypePerDisc *NewOnRead(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int profileIt) throw(INTERP_KERNEL::Exception);
    static MEDFileFieldPerMeshPerTypePerDisc *New(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int locId);
    static MEDFileFieldPerMeshPerTypePerDisc *New(const MEDFileFieldPerMeshPerTypePerDisc& other);
    void assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignFieldProfile(int& start, const char *pflName, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void getCoarseData(TypeOfField& type, std::pair<int,int>& dad, std::string& pfl, std::string& loc) const throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerType *getFather() const;
    void prepareLoading(med_idt fid, int profileIt, int& start) throw(INTERP_KERNEL::Exception);
    void finishLoading(med_idt fid, int profileIt, int ft) throw(INTERP_KERNEL::Exception);
    void setNewStart(int newValueOfStart) throw(INTERP_KERNEL::Exception);
    int getIteration() const;
    int getOrder() const;
    double getTime() const;
    std::string getName() const;
    std::string getMeshName() const;
    TypeOfField getType() const;
    void simpleRepr(int bkOffset, std::ostream& oss, int id) const;
    void fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const throw(INTERP_KERNEL::Exception);
    void setType(TypeOfField newType);
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    int getNumberOfComponents() const;
    int getNumberOfTuples() const;
    int getStart() const { return _start; }
    DataArrayDouble *getArray();
    const DataArrayDouble *getArray() const;
    const std::vector<std::string>& getInfo() const;
    std::string getProfile() const;
    void setProfile(const char *newPflName);
    std::string getLocalization() const;
    void setLocalization(const char *newLocName);
    int getLocId() const { return _loc_id; }
    void setLocId(int newId) { _loc_id=newId; }
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
                                                                          bool isPfl, int offset, std::list< const MEDFileFieldPerMeshPerTypePerDisc *>& entriesOnSameDisc,
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
  };

  class MEDFileFieldPerMeshPerType : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMeshPerType *New(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception);
    static MEDFileFieldPerMeshPerType *NewOnRead(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception);
    void assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, DataArrayInt *locIds, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMesh *getFather() const;
    void prepareLoading(med_idt fid, int &start) throw(INTERP_KERNEL::Exception);
    void finishLoading(med_idt fid, int ft) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    void getDimension(int& dim) const;
    void fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const throw(INTERP_KERNEL::Exception);
    void fillFieldSplitedByType(std::vector< std::pair<int,int> >& dads, std::vector<TypeOfField>& types, std::vector<std::string>& pfls, std::vector<std::string>& locs) const throw(INTERP_KERNEL::Exception);
    int getIteration() const;
    int getOrder() const;
    double getTime() const;
    std::string getName() const;
    std::string getMeshName() const;
    void simpleRepr(int bkOffset, std::ostream& oss, int id) const;
    void getSizes(int& globalSz, int& nbOfEntries) const;
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    int getNumberOfComponents() const;
    DataArrayDouble *getArray();
    const DataArrayDouble *getArray() const;
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
    static med_entity_type ConvertIntoMEDFileType(TypeOfField ikType, INTERP_KERNEL::NormalizedCellType ikGeoType, med_geometry_type& medfGeoType);
  private:
    std::vector<int> addNewEntryIfNecessary(const MEDCouplingFieldDouble *field, int offset, int nbOfCells) throw(INTERP_KERNEL::Exception);
    std::vector<int> addNewEntryIfNecessaryGauss(const MEDCouplingFieldDouble *field, int offset, int nbOfCells) throw(INTERP_KERNEL::Exception);
    std::vector<int> addNewEntryIfNecessary(const MEDCouplingFieldDouble *field, const DataArrayInt *subCells) throw(INTERP_KERNEL::Exception);
    std::vector<int> addNewEntryIfNecessaryGauss(const MEDCouplingFieldDouble *field, const DataArrayInt *subCells) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerType(med_idt fid, MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerType(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFieldPerMesh *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> > _field_pm_pt_pd;
    INTERP_KERNEL::NormalizedCellType _geo_type;
  };

  class MEDFileFieldPerMesh : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMesh *New(MEDFileField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh);
    static MEDFileFieldPerMesh *NewOnRead(med_idt fid, MEDFileField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder) throw(INTERP_KERNEL::Exception);
    void simpleRepr(int bkOffset,std::ostream& oss, int id) const;
    void copyTinyInfoFrom(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    void assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const std::vector<int>& code, const std::vector<DataArrayInt *>& idsInPflPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignFieldProfileGeneral(int& start, const DataArrayInt *multiTypePfl, const std::vector<int>& code, const std::vector<DataArrayInt *>& idsInPflPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignFieldNoProfileNoRenum(int& start, const std::vector<int>& code, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void prepareLoading(med_idt fid, int &start) throw(INTERP_KERNEL::Exception);
    void finishLoading(med_idt fid, int ft) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    void fillTypesOfFieldAvailable(std::set<TypeOfField>& types) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector< std::pair<int,int> > > getFieldSplitedByType(std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    void getDimension(int& dim) const;
    double getTime() const;
    int getIteration() const;
    int getOrder() const;
    int getMeshIteration() const { return _mesh_iteration; }
    int getMeshOrder() const { return _mesh_order; }
    const std::string& getDtUnit() const;
    std::string getName() const;
    std::string getMeshName() const { return _mesh_name; }
    int getNumberOfComponents() const;
    DataArrayDouble *getArray();
    const DataArrayDouble *getArray() const;
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    bool renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, bool& isPfl) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getFieldOnMeshAtLevelWithPfl(TypeOfField type, const MEDCouplingMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception);
  private:
    int addNewEntryIfNecessary(INTERP_KERNEL::NormalizedCellType type);
    MEDCouplingFieldDouble *finishField(TypeOfField type, const MEDFileFieldGlobsReal *glob,
                                        const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs, const MEDCouplingMesh *mesh, bool& isPfl) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *finishField2(TypeOfField type, const MEDFileFieldGlobsReal *glob,
                                         const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs,
                                         const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes,
                                         const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *finishField3(const MEDFileFieldGlobsReal *glob,
                                         const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs,
                                         const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *finishField4(const std::vector< std::pair<int,int> >& dads, const DataArrayInt *pflIn, int nbOfElems, DataArrayInt *&pflOut) const throw(INTERP_KERNEL::Exception);
    void assignNewLeaves(const std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerTypePerDisc > >& leaves) throw(INTERP_KERNEL::Exception);
    static void SortArraysPerType(const MEDFileFieldGlobsReal *glob, TypeOfField type, 
                                  const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& locs,
                                  std::vector<int>& code, std::vector<DataArrayInt *>& notNullPfls);
    static int ComputeNbOfElems(const MEDFileFieldGlobsReal *glob, TypeOfField type, const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMesh(med_idt fid, MEDFileField1TSWithoutSDA *fath, int meshCsit, int meshIteration, int meshOrder) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMesh(MEDFileField1TSWithoutSDA *fath, const MEDCouplingMesh *mesh);
  private:
    std::string _mesh_name;
    int _mesh_iteration;
    int _mesh_order;
    int _mesh_csit;
    MEDFileField1TSWithoutSDA *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > > _field_pm_pt;
  };

  class MEDFileFieldGlobsReal;

  class MEDFileFieldGlobs : public RefCountObject
  {
  public:
    static MEDFileFieldGlobs *New(const char *fname);
    static MEDFileFieldGlobs *New();
    void simpleRepr(std::ostream& oss) const;
    void appendGlobs(const MEDFileFieldGlobs& other, double eps) throw(INTERP_KERNEL::Exception);
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
    void simpleRepr(std::ostream& oss) const;
    void shallowCpyGlobs(const MEDFileFieldGlobsReal& other);
    void appendGlobs(const MEDFileFieldGlobsReal& other, double eps) throw(INTERP_KERNEL::Exception);
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
    MEDCouplingAutoRefCountObjectPtr< MEDFileFieldGlobs > _globals;
  };

  /*!
   * SDA is for Shared Data Arrays such as profiles.
   */
  class MEDLOADER_EXPORT MEDFileField1TSWithoutSDA : public RefCountObject
  {
  public:
    int copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    int getDimension() const;
    int getIteration() const { return _iteration; }
    int getOrder() const { return _order; }
    double getTime(int& iteration, int& order) const { iteration=_iteration; order=_order; return _dt; }
    void setTime(int iteration, int order, double val) { _dt=val; _iteration=iteration; _order=order; }
    std::string getName() const;
    void setName(const char *name);
    void simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const;
    const std::string& getDtUnit() const { return _dt_unit; }
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
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string>& getInfo();
    //
    std::vector<std::string> getPflsReallyUsed2() const;
    std::vector<std::string> getLocsReallyUsed2() const;
    std::vector<std::string> getPflsReallyUsedMulti2() const;
    std::vector<std::string> getLocsReallyUsedMulti2() const;
    void changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArray() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
    int getNonEmptyLevels(const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception);
    std::vector<TypeOfField> getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<std::pair<int,int> > > getFieldSplitedByType(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
     std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
     MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception);
     const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception);
    static void CheckMeshDimRel(int meshDimRelToMax) throw(INTERP_KERNEL::Exception);
    static std::vector<int> CheckSBTMesh(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    static MEDFileField1TSWithoutSDA *New(const char *fieldName, int csit, int fieldtype, int iteration, int order, const std::vector<std::string>& infos);
  public:
    void finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception);
    virtual void writeLL(med_idt fid, const MEDFileWritable& opts) const throw(INTERP_KERNEL::Exception);
    void setFieldNoProfileSBT(const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void setFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    bool renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, const char *mName, int renumPol, const MEDFileFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDFileMesh *mesh) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, const char *mName, int renumPol, const MEDFileFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, const DataArrayInt *cellRenum, const DataArrayInt *nodeRenum) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception);
  protected:
    int addNewEntryIfNecessary(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    int getMeshIdFromMeshName(const char *mName) const throw(INTERP_KERNEL::Exception);
  public:
    MEDFileField1TSWithoutSDA();
    MEDFileField1TSWithoutSDA(const char *fieldName, int csit, int fieldtype, int iteration, int order, const std::vector<std::string>& infos);
    DataArrayDouble *getOrCreateAndGetArray();
    const DataArrayDouble *getOrCreateAndGetArray() const;
  protected:
    std::string _dt_unit;
    MEDCouplingAutoRefCountObjectPtr< DataArrayDouble > _arr;
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > > _field_per_mesh;
    int _iteration;
    int _order;
    double _dt;
  public:
    //! only useable on reading
    mutable int _csit;
    //! only useable on reading. 0 is for float, 1 for int32, 2 for int64
    mutable int _field_type;
  };

  /*!
   * User class.
   */
  class MEDLOADER_EXPORT MEDFileField1TS : public RefCountObject, public MEDFileWritable, public MEDFileFieldGlobsReal
  {
  public:
    static MEDFileField1TS *New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
    static MEDFileField1TS *New(const MEDFileField1TSWithoutSDA& other, bool deepCpy);
    static MEDFileField1TS *New();
    std::string simpleRepr() const;
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
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
    int copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    int getDimension() const;
    int getIteration() const;
    int getOrder() const;
    double getTime(int& iteration, int& order) const;
    void setTime(int iteration, int order, double val);
    std::string getName() const;
    void setName(const char *name);
    void simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const;
    const std::string& getDtUnit() const;
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
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string>& getInfo();
    DataArrayDouble *getUndergroundDataArray() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception);
    int getNonEmptyLevels(const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception);
    std::vector<TypeOfField> getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<std::pair<int,int> > > getFieldSplitedByType(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
                                                                          std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
                                                                         std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
  public:
    //! underground method see MEDFileField1TSWithoutSDA::setProfileNameOnLeaf
    void setProfileNameOnLeaf(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newPflName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception);
    //! underground method see MEDFileField1TSWithoutSDA::setLocNameOnLeaf
    void setLocNameOnLeaf(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newLocName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception);
  private:
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
    MEDFileField1TS(const MEDFileField1TSWithoutSDA& other, bool deepCpy);
    MEDFileField1TS();
  protected:
    MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA> _content;
  };
  
  class MEDLOADER_EXPORT MEDFileFieldMultiTSWithoutSDA : public RefCountObject
  {
  public:
    static MEDFileFieldMultiTSWithoutSDA *New(med_idt fid, const char *fieldName, int id, int ft, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTSWithoutSDA(med_idt fid, int fieldId) throw(INTERP_KERNEL::Exception);
    int getNumberOfTS() const;
    void eraseEmptyTS() throw(INTERP_KERNEL::Exception);
    void eraseTimeStepIds(const int *startIds, const int *endIds) throw(INTERP_KERNEL::Exception);
    std::vector< std::pair<int,int> > getIterations() const;
    int getPosOfTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    int getPosGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    int getNonEmptyLevels(int iteration, int order, const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<TypeOfField> > getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector< std::pair<int,int> > > getFieldSplitedByType(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    virtual void writeLL(med_idt fid, const MEDFileWritable& opts) const throw(INTERP_KERNEL::Exception);
    std::string getName() const;
    void setName(const char *name);
    void simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const;
    std::vector< std::pair<int,int> > getTimeSteps(std::vector<double>& ret1) const throw(INTERP_KERNEL::Exception);
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    void setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception);
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    const std::vector<std::string>& getInfo() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
    bool renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void appendFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
  public:
    const MEDFileField1TSWithoutSDA *getTimeStepAtPos2(int pos) const throw(INTERP_KERNEL::Exception);
    MEDFileField1TSWithoutSDA *getTimeStepAtPos2(int pos) throw(INTERP_KERNEL::Exception);
    const MEDFileField1TSWithoutSDA& getTimeStepEntry(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPflsReallyUsed2() const;
    std::vector<std::string> getLocsReallyUsed2() const;
    std::vector<std::string> getPflsReallyUsedMulti2() const;
    std::vector<std::string> getLocsReallyUsedMulti2() const;
    void changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
  protected:
    MEDFileField1TSWithoutSDA& getTimeStepEntry(int iteration, int order) throw(INTERP_KERNEL::Exception);
    std::string getDtUnit() const throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTSWithoutSDA(const char *fieldName);
    MEDFileFieldMultiTSWithoutSDA(med_idt fid, const char *fieldName, int id, int ft, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception);
    void finishLoading(med_idt fid, int nbPdt) throw(INTERP_KERNEL::Exception);
    void copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    void checkCoherencyOfTinyInfo(const MEDCouplingFieldDouble *field) const throw(INTERP_KERNEL::Exception);
  public:
    MEDFileFieldMultiTSWithoutSDA();
  protected:
    std::string _name;
    std::vector<std::string> _infos;
    //! only useable on reading. 0 is for float, 1 for int32, 2 for int64
    mutable int _field_type;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutSDA> > _time_steps;
  };

  class MEDFileFieldMultiTSIterator;

  /*!
   * User class.
   */
  class MEDLOADER_EXPORT MEDFileFieldMultiTS : public RefCountObject, public MEDFileWritable, public MEDFileFieldGlobsReal
  {
  public:
    static MEDFileFieldMultiTS *New();
    static MEDFileFieldMultiTS *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileFieldMultiTS *New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception);
    static MEDFileFieldMultiTS *New(const MEDFileFieldMultiTSWithoutSDA& other, bool deepCpy);
    //
    MEDFileField1TS *getTimeStepAtPos(int pos) const throw(INTERP_KERNEL::Exception);
    MEDFileField1TS *getTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    MEDFileField1TS *getTimeStepGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTSIterator *iterator() throw(INTERP_KERNEL::Exception);
    //
    std::string simpleRepr() const;
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevelOld(TypeOfField type, const char *mname, int iteration, int order, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const throw(INTERP_KERNEL::Exception);
    //
    void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    void appendFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile) throw(INTERP_KERNEL::Exception);
  public:// direct forwarding to MEDFileField1TSWithoutSDA instance _content
    int getNumberOfTS() const;
    void eraseEmptyTS() throw(INTERP_KERNEL::Exception);
    void eraseTimeStepIds(const int *startIds, const int *endIds) throw(INTERP_KERNEL::Exception);
    std::vector< std::pair<int,int> > getIterations() const;
    int getPosOfTimeStep(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    int getPosGivenTime(double time, double eps=1e-8) const throw(INTERP_KERNEL::Exception);
    int getNonEmptyLevels(int iteration, int order, const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<TypeOfField> > getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector< std::pair<int,int> > > getFieldSplitedByType(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    std::string getName() const;
    void setName(const char *name);
    void simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const;
    std::vector< std::pair<int,int> > getTimeSteps(std::vector<double>& ret1) const throw(INTERP_KERNEL::Exception);
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    void setMeshName(const char *newMeshName) throw(INTERP_KERNEL::Exception);
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    const std::vector<std::string>& getInfo() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
  public:
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
  public:
    MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutSDA> getContent();
  private:
    MEDFileFieldMultiTS();
    MEDFileFieldMultiTS(const MEDFileFieldMultiTSWithoutSDA& other, bool deepCpy);
    MEDFileFieldMultiTS(const char *fileName) throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  protected:
    MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutSDA> _content;
  };

  class MEDCOUPLING_EXPORT MEDFileFieldMultiTSIterator
  {
  public:
    MEDFileFieldMultiTSIterator(MEDFileFieldMultiTS *fmts);
    ~MEDFileFieldMultiTSIterator();
    MEDFileField1TS *nextt();
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> _fmts;
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
    static MEDFileFields *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    int getNumberOfFields() const;
    std::vector<std::string> getFieldsNames() const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getMeshesNames() const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    void simpleRepr(int bkOffset, std::ostream& oss) const;
    //
    void resize(int newSize) throw(INTERP_KERNEL::Exception);
    void pushField(MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception);
    void setFieldAtPos(int i, MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTS *getFieldAtPos(int i) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTS *getFieldWithName(const char *fieldName) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldsIterator *iterator() throw(INTERP_KERNEL::Exception);
    void destroyFieldAtPos(int i) throw(INTERP_KERNEL::Exception);
    bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    bool renumberEntitiesLyingOnMesh(const char *meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N) throw(INTERP_KERNEL::Exception);
  private:
    int getPosFromFieldName(const char *fieldName) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
    void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFields();
    MEDFileFields(const char *fileName) throw(INTERP_KERNEL::Exception);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutSDA> > _fields;
  };

  class MEDCOUPLING_EXPORT MEDFileFieldsIterator
  {
  public:
    MEDFileFieldsIterator(MEDFileFields *fs);
    ~MEDFileFieldsIterator();
    MEDFileFieldMultiTS *nextt();
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDFileFields> _fs;
     int _iter_id;
     int _nb_iter;
  };
}

#endif
