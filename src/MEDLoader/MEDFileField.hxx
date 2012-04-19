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

extern "C"
{
#include "med.h"
}

namespace ParaMEDMEM
{
  class MEDFieldFieldGlobs;
  class MEDCouplingMesh;
  class MEDCouplingFieldDouble;
  class MEDFileMesh;
/// @cond INTERNAL
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

  class MEDFileFieldPerMeshPerType;
  class MEDFileField1TSWithoutDAS;
  class MEDFieldFieldGlobsReal;
  class MEDFileFieldPerMesh;

  class MEDFileFieldPerMeshPerTypePerDisc : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMeshPerTypePerDisc *NewOnRead(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int profileIt) throw(INTERP_KERNEL::Exception);
    static MEDFileFieldPerMeshPerTypePerDisc *New(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int locId);
    void assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignFieldProfile(int& start, const char *pflName, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void getCoarseData(TypeOfField& type, std::pair<int,int>& dad, std::string& pfl, std::string& loc) const throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerType *getFather() const;
    void prepareLoading(med_idt fid, int profileIt, int& start) throw(INTERP_KERNEL::Exception);
    void finishLoading(med_idt fid, int profileIt, int ft) throw(INTERP_KERNEL::Exception);
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
    DataArrayDouble *getArray();
    const DataArrayDouble *getArray() const;
    const std::vector<std::string>& getInfo() const;
    std::string getProfile() const;
    void setProfile(const char *newPflName);
    std::string getLocalization() const;
    void setLocalization(const char *newLocName);
    int getLocId() const { return _loc_id; }
    void getFieldAtLevel(TypeOfField type, const MEDFieldFieldGlobsReal *glob, std::vector< std::pair<int,int> >& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs,
                         std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const;
    void fillValues(int discId, int& startEntryId, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    static int ConvertType(TypeOfField type, int locId) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int profileIt) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, TypeOfField type, int profileIt, const std::string& dummy);
  private:
    TypeOfField _type;
    MEDFileFieldPerMeshPerType *_father;
    int _start;
    int _end;
    //! _nval is different than end-start in case of ON_GAUSS_PT and ON_GAUSS_NE ! (_nval=(_end-_start)/nbi)
    int _nval;
    int _profile_it;
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
    void assignFieldNoProfile(int& start, int offset, int nbOfCells, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const DataArrayInt *idsInPfl, DataArrayInt *locIds, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
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
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenLocId(int locId) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenLocId(int locId) const throw(INTERP_KERNEL::Exception);
    void getFieldAtLevel(int meshDim, TypeOfField type, const MEDFieldFieldGlobsReal *glob, std::vector< std::pair<int,int> >& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const;
    void fillValues(int& startEntryId, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
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
    static MEDFileFieldPerMesh *New(MEDFileField1TSWithoutDAS *fath, const MEDCouplingMesh *mesh);
    static MEDFileFieldPerMesh *NewOnRead(med_idt fid, MEDFileField1TSWithoutDAS *fath, int meshCsit, int meshIteration, int meshOrder) throw(INTERP_KERNEL::Exception);
    void simpleRepr(int bkOffset,std::ostream& oss, int id) const;
    void copyTinyInfoFrom(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    void assignFieldProfile(int& start, const DataArrayInt *multiTypePfl, const std::vector<int>& code, const std::vector<DataArrayInt *>& idsInPflPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignFieldProfileGeneral(int& start, const DataArrayInt *multiTypePfl, const std::vector<int>& code, const std::vector<DataArrayInt *>& idsInPflPerType, const std::vector<DataArrayInt *>& idsPerType, const MEDCouplingFieldDouble *field, const MEDCouplingMesh *mesh, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignFieldNoProfileNoRenum(int& start, const std::vector<int>& code, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldNoProfile(int& start, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void assignNodeFieldProfile(int& start, const DataArrayInt *pfl, const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
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
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDFieldFieldGlobsReal *glob, const MEDCouplingMesh *mesh, bool& isPfl) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getFieldOnMeshAtLevelWithPfl(TypeOfField type, const MEDCouplingMesh *mesh, DataArrayInt *&pfl, const MEDFieldFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenTypeAndLocId(INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception);
  private:
    int addNewEntryIfNecessary(INTERP_KERNEL::NormalizedCellType type);
    MEDCouplingFieldDouble *finishField(TypeOfField type, const MEDFieldFieldGlobsReal *glob,
                                        const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs, const MEDCouplingMesh *mesh, bool& isPfl) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *finishField2(TypeOfField type, const MEDFieldFieldGlobsReal *glob,
                                         const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs,
                                         const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes,
                                         const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *finishField3(const MEDFieldFieldGlobsReal *glob,
                                         const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs,
                                         const MEDCouplingMesh *mesh, const DataArrayInt *da, bool& isPfl) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *finishField4(const std::vector< std::pair<int,int> >& dads, const DataArrayInt *pflIn, int nbOfElems, DataArrayInt *&pflOut) const throw(INTERP_KERNEL::Exception);
    static void SortArraysPerType(const MEDFieldFieldGlobsReal *glob, TypeOfField type, 
                                  const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& locs,
                                  std::vector<int>& code, std::vector<DataArrayInt *>& notNullPfls);
    static int ComputeNbOfElems(const MEDFieldFieldGlobsReal *glob, TypeOfField type, const std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, const std::vector< std::pair<int,int> >& dads, const std::vector<int>& locs) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMesh(med_idt fid, MEDFileField1TSWithoutDAS *fath, int meshCsit, int meshIteration, int meshOrder) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMesh(MEDFileField1TSWithoutDAS *fath, const MEDCouplingMesh *mesh);
  private:
    std::string _mesh_name;
    int _mesh_iteration;
    int _mesh_order;
    int _mesh_csit;
    MEDFileField1TSWithoutDAS *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > > _field_pm_pt;
  };

  class MEDFieldFieldGlobsReal;

  class MEDFieldFieldGlobs : public RefCountObject
  {
  public:
    static MEDFieldFieldGlobs *New(const char *fname);
    static MEDFieldFieldGlobs *New();
    void simpleRepr(std::ostream& oss) const;
    void appendGlobs(const MEDFieldFieldGlobs& other, double eps) throw(INTERP_KERNEL::Exception);
    void loadProfileInFile(med_idt fid, int id, const char *pflName) throw(INTERP_KERNEL::Exception);
    void loadProfileInFile(med_idt fid, int id);
    void loadGlobals(med_idt fid, const MEDFieldFieldGlobsReal& real) throw(INTERP_KERNEL::Exception);
    void loadAllGlobals(med_idt fid) throw(INTERP_KERNEL::Exception);
    void writeGlobals(med_idt fid, const MEDFileWritable& opt) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPfls() const;
    std::vector<std::string> getLocs() const;
    bool existsPfl(const char *pflName) const;
    bool existsLoc(const char *locName) const;
    void setFileName(const char *fileName);
    int getNbOfGaussPtPerCell(int locId) const throw(INTERP_KERNEL::Exception);
    int getLocalizationId(const char *loc) const throw(INTERP_KERNEL::Exception);
    const char *getFileName() const { return _file_name.c_str(); }
    std::string getFileName2() const { return _file_name; }
    const MEDFileFieldLoc& getLocalizationFromId(int locId) const throw(INTERP_KERNEL::Exception);
    const MEDFileFieldLoc& getLocalization(const char *locName) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getProfile(const char *pflName) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldLoc& getLocalizationFromId(int locId) throw(INTERP_KERNEL::Exception);
    MEDFileFieldLoc& getLocalization(const char *locName) throw(INTERP_KERNEL::Exception);
    DataArrayInt *getProfile(const char *pflName) throw(INTERP_KERNEL::Exception);
    //
    void appendProfile(DataArrayInt *pfl) throw(INTERP_KERNEL::Exception);
    void appendLoc(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w) throw(INTERP_KERNEL::Exception);
  protected:
    MEDFieldFieldGlobs(const char *fname);
    MEDFieldFieldGlobs();
    ~MEDFieldFieldGlobs();
  protected:
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > _pfls;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> > _locs;
    std::string _file_name;
  };

/// @endcond INTERNAL

  class MEDLOADER_EXPORT MEDFieldFieldGlobsReal
  {
  public:
    MEDFieldFieldGlobsReal(const char *fname);
    MEDFieldFieldGlobsReal();
    void simpleRepr(std::ostream& oss) const;
    void shallowCpyGlobs(const MEDFieldFieldGlobsReal& other);
    void appendGlobs(const MEDFieldFieldGlobsReal& other, double eps) throw(INTERP_KERNEL::Exception);
    virtual std::vector<std::string> getPflsReallyUsed() const = 0;
    virtual std::vector<std::string> getLocsReallyUsed() const = 0;
    virtual std::vector<std::string> getPflsReallyUsedMulti() const = 0;
    virtual std::vector<std::string> getLocsReallyUsedMulti() const = 0;
    virtual ~MEDFieldFieldGlobsReal();
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
    void setFileName(const char *fileName);
    int getNbOfGaussPtPerCell(int locId) const throw(INTERP_KERNEL::Exception);
    int getLocalizationId(const char *loc) const throw(INTERP_KERNEL::Exception);
    const char *getFileName() const;
    std::string getFileName2() const;
    const MEDFileFieldLoc& getLocalizationFromId(int locId) const throw(INTERP_KERNEL::Exception);
    const MEDFileFieldLoc& getLocalization(const char *locName) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getProfile(const char *pflName) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldLoc& getLocalizationFromId(int locId) throw(INTERP_KERNEL::Exception);
    MEDFileFieldLoc& getLocalization(const char *locName) throw(INTERP_KERNEL::Exception);
    DataArrayInt *getProfile(const char *pflName) throw(INTERP_KERNEL::Exception);
    //
    void appendProfile(DataArrayInt *pfl) throw(INTERP_KERNEL::Exception);
    void appendLoc(const char *locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w) throw(INTERP_KERNEL::Exception);
  protected:
    MEDCouplingAutoRefCountObjectPtr< MEDFieldFieldGlobs > _globals;
  };

  /*!
   * DAS is for Shared Data Arrays such as profiles.
   */
  class MEDLOADER_EXPORT MEDFileField1TSWithoutDAS : public RefCountObject, public MEDFileWritable
  {
  public:
    int copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    int getDimension() const;
    int getIteration() const { return _iteration; }
    int getOrder() const { return _order; }
    double getTime(int& iteration, int& order) const { iteration=_iteration; order=_order; return _dt; }
    void setTime(double val, int iteration, int order) { _dt=val; _iteration=iteration; _order=order; }
    std::string getName() const;
    void simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const;
    const std::string& getDtUnit() const { return _dt_unit; }
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    int getMeshIteration() const throw(INTERP_KERNEL::Exception);
    int getMeshOrder() const throw(INTERP_KERNEL::Exception);
    int getNumberOfComponents() const;
    void simpleReprWithoutHeader(std::ostream& oss) const;
    bool isDealingTS(int iteration, int order) const;
    std::pair<int,int> getDtIt() const;
    void fillIteration(std::pair<int,int>& p) const;
    void fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const throw(INTERP_KERNEL::Exception);
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string>& getInfo();
    //
    static MEDFileField1TSWithoutDAS *New(const char *fieldName, int csit, int fieldtype, int iteration, int order, const std::vector<std::string>& infos);
    void finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception);
    virtual void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPflsReallyUsed2() const;
    std::vector<std::string> getLocsReallyUsed2() const;
    std::vector<std::string> getPflsReallyUsedMulti2() const;
    std::vector<std::string> getLocsReallyUsedMulti2() const;
    static void CheckMeshDimRel(int meshDimRelToMax) throw(INTERP_KERNEL::Exception);
    static std::vector<int> CheckSBTMesh(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    //
    void setFieldNoProfileSBT(const MEDCouplingFieldDouble *field, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
    void setFieldProfile(const MEDCouplingFieldDouble *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFieldFieldGlobsReal& glob) throw(INTERP_KERNEL::Exception);
  public:
    int getNonEmptyLevels(const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception);
    std::vector<TypeOfField> getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<std::pair<int,int> > > getFieldSplitedByType(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, const char *mName, int renumPol, const MEDFieldFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol, const MEDFieldFieldGlobsReal *glob, const MEDFileMesh *mesh) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, const char *mName, int renumPol, const MEDFieldFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int renumPol, const MEDFieldFieldGlobsReal *glob, const MEDCouplingMesh *mesh, const DataArrayInt *cellRenum, const DataArrayInt *nodeRenum) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl, const MEDFieldFieldGlobsReal *glob) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArray() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
  protected:
    int addNewEntryIfNecessary(const MEDCouplingMesh *mesh) throw(INTERP_KERNEL::Exception);
    int getMeshIdFromMeshName(const char *mName) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const throw(INTERP_KERNEL::Exception);
    MEDFileField1TSWithoutDAS(const char *fieldName, int csit, int fieldtype, int iteration, int order, const std::vector<std::string>& infos);
  public:
    MEDFileField1TSWithoutDAS();
    DataArrayDouble *getOrCreateAndGetArray();
    const DataArrayDouble *getOrCreateAndGetArray() const;
  protected:
    std::string _dt_unit;
    MEDCouplingAutoRefCountObjectPtr< DataArrayDouble > _arr;
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > > _field_per_mesh;
    //! only useable on reading
    mutable int _csit;
    //! only useable on reading. 0 is for float, 1 for int32, 2 for int64
    mutable int _field_type;
    int _iteration;
    int _order;
    double _dt;
  };

  /*!
   * User class.
   */
  class MEDLOADER_EXPORT MEDFileField1TS : public MEDFileField1TSWithoutDAS, public MEDFieldFieldGlobsReal
  {
  public:
    static MEDFileField1TS *New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
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
  public:
    //! underground method see MEDFileField1TSWithoutDAS::setProfileNameOnLeaf
    void setProfileNameOnLeaf(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newPflName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception);
    //! underground method see MEDFileField1TSWithoutDAS::setLocNameOnLeaf
    void setLocNameOnLeaf(const char *mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const char *newLocName, bool forceRenameOnGlob=false) throw(INTERP_KERNEL::Exception);
  private:
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
    MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
    MEDFileField1TS();
  };
  
  class MEDLOADER_EXPORT MEDFileFieldMultiTSWithoutDAS : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldMultiTSWithoutDAS *New(med_idt fid, const char *fieldName, int id, int ft, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception);
    int getNumberOfTS() const;
    std::vector< std::pair<int,int> > getIterations() const;
    int getNonEmptyLevels(int iteration, int order, const char *mname, std::vector<int>& levs) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<TypeOfField> > getTypesOfFieldAvailable() const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector< std::pair<int,int> > > getFieldSplitedByType(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(int iteration, int order, const char *mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const throw(INTERP_KERNEL::Exception);
    virtual void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    std::string getName() const;
    void simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const;
    std::vector< std::pair<int,int> > getTimeSteps(std::vector<double>& ret1) const throw(INTERP_KERNEL::Exception);
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    const std::vector<std::string>& getInfo() const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArray(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const throw(INTERP_KERNEL::Exception);
  public:
    std::vector<std::string> getPflsReallyUsed2() const;
    std::vector<std::string> getLocsReallyUsed2() const;
    std::vector<std::string> getPflsReallyUsedMulti2() const;
    std::vector<std::string> getLocsReallyUsedMulti2() const;
  protected:
    const MEDFileField1TSWithoutDAS& getTimeStepEntry(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    MEDFileField1TSWithoutDAS& getTimeStepEntry(int iteration, int order) throw(INTERP_KERNEL::Exception);
    std::string getDtUnit() const throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTSWithoutDAS();
    MEDFileFieldMultiTSWithoutDAS(const char *fieldName);
    MEDFileFieldMultiTSWithoutDAS(med_idt fid, const char *fieldName, int id, int ft, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception);
    void finishLoading(med_idt fid, int nbPdt) throw(INTERP_KERNEL::Exception);
    void copyTinyInfoFrom(const MEDCouplingFieldDouble *field) throw(INTERP_KERNEL::Exception);
    void checkCoherencyOfTinyInfo(const MEDCouplingFieldDouble *field) const throw(INTERP_KERNEL::Exception);
  protected:
    std::string _name;
    std::vector<std::string> _infos;
    //! only useable on reading. 0 is for float, 1 for int32, 2 for int64
    mutable int _field_type;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutDAS>  > _time_steps;
  };

  /*!
   * User class.
   */
  class MEDLOADER_EXPORT MEDFileFieldMultiTS : public MEDFileFieldMultiTSWithoutDAS, public MEDFieldFieldGlobsReal
  {
  public:
    static MEDFileFieldMultiTS *New();
    static MEDFileFieldMultiTS *New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception);
    static MEDFileFieldMultiTS *New(const MEDFileFieldMultiTSWithoutDAS& other);
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
  private:
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
  private:
    MEDFileFieldMultiTS();
    MEDFileFieldMultiTS(const MEDFileFieldMultiTSWithoutDAS& other);
    MEDFileFieldMultiTS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  };

  /*!
   * Use class.
   */
  class MEDLOADER_EXPORT MEDFileFields : public RefCountObject, public MEDFieldFieldGlobsReal, public MEDFileWritable
  {
  public:
    static MEDFileFields *New();
    static MEDFileFields *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    int getNumberOfFields() const;
    std::vector<std::string> getFieldsNames() const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    void simpleRepr(int bkOffset, std::ostream& oss) const;
    //
    void resize(int newSize) throw(INTERP_KERNEL::Exception);
    void pushField(MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception);
    void setFieldAtPos(int i, MEDFileFieldMultiTS *field) throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTS *getFieldAtPos(int i) const throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTS *getFieldWithName(const char *fieldName) const throw(INTERP_KERNEL::Exception);
    void destroyFieldAtPos(int i) throw(INTERP_KERNEL::Exception);
  private:
    int getPosFromFieldName(const char *fieldName) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<std::string> getPflsReallyUsedMulti() const;
    std::vector<std::string> getLocsReallyUsedMulti() const;
  private:
    MEDFileFields();
    MEDFileFields(const char *fileName) throw(INTERP_KERNEL::Exception);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutDAS> > _fields;
  };
}

#endif
