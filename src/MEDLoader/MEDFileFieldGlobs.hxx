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
// Author : Anthony Geay (EDF R&D)

#ifndef __MEDFILEFIELDGLOBS_HXX__
#define __MEDFILEFIELDGLOBS_HXX__

#include "MEDLoaderDefines.hxx"

#include "NormalizedGeometricTypes"
#include "MEDCouplingMemArray.hxx"
#include "MCAuto.hxx"

#include "med.h"

namespace MEDCoupling
{
  class MEDFileFieldGlobsReal;
  class MEDFileEntities;
  class MEDFileWritable;
  class MEDFileFieldLoc;
  
  class MEDFileFieldGlobs : public RefCountObject
  {
  public:
    static MEDFileFieldGlobs *New(med_idt fid);
    static MEDFileFieldGlobs *New();
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDFileFieldGlobs *deepCopy() const;
    MEDFileFieldGlobs *shallowCpyPart(const std::vector<std::string>& pfls, const std::vector<std::string>& locs) const;
    MEDFileFieldGlobs *deepCpyPart(const std::vector<std::string>& pfls, const std::vector<std::string>& locs) const;
    void simpleRepr(std::ostream& oss) const;
    void appendGlobs(const MEDFileFieldGlobs& other, double eps);
    void checkGlobsPflsPartCoherency(const std::vector<std::string>& pflsUsed) const;
    void checkGlobsLocsPartCoherency(const std::vector<std::string>& locsUsed) const;
    void loadProfileInFile(med_idt fid, int id, const std::string& pflName);
    void loadProfileInFile(med_idt fid, int id);
    void loadGlobals(med_idt fid, const MEDFileFieldGlobsReal& real);
    void loadAllGlobals(med_idt fid, const MEDFileEntities *entities);
    void writeGlobals(med_idt fid, const MEDFileWritable& opt) const;
    std::vector<std::string> getPfls() const;
    std::vector<std::string> getLocs() const;
    bool existsPfl(const std::string& pflName) const;
    bool existsLoc(const std::string& locName) const;
    std::string createNewNameOfPfl() const;
    std::string createNewNameOfLoc() const;
    std::vector< std::vector<int> > whichAreEqualProfiles() const;
    std::vector< std::vector<int> > whichAreEqualLocs(double eps) const;
    void setFileName(const std::string& fileName) { _file_name=fileName; }
    void changePflsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    void changeLocsNamesInStruct(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    int getNbOfGaussPtPerCell(int locId) const;
    int getLocalizationId(const std::string& loc) const;
    int getProfileId(const std::string& pfl) const;
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
    void killStructureElementsInGlobs();
    //
    void appendProfile(DataArrayInt *pfl);
    void appendLoc(const std::string& locName, INTERP_KERNEL::NormalizedCellType geoType, const std::vector<double>& refCoo, const std::vector<double>& gsCoo, const std::vector<double>& w);
    //
    static std::string CreateNewNameNotIn(const std::string& prefix, const std::vector<std::string>& namesToAvoid);
  protected:
    MEDFileFieldGlobs(med_idt fid);
    MEDFileFieldGlobs();
    ~MEDFileFieldGlobs();
  protected:
    std::vector< MCAuto<DataArrayInt> > _pfls;
    std::vector< MCAuto<MEDFileFieldLoc> > _locs;
    std::string _file_name;
  };

  /// @endcond INTERNAL

  class MEDFileFieldGlobsReal
  {
  public:
    MEDLOADER_EXPORT MEDFileFieldGlobsReal(med_idt fid);
    MEDLOADER_EXPORT MEDFileFieldGlobsReal();
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT void simpleReprGlobs(std::ostream& oss) const;
    MEDLOADER_EXPORT void resetContent();
    MEDLOADER_EXPORT void killStructureElementsInGlobs();
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
    MEDLOADER_EXPORT void loadAllGlobals(med_idt fid, const MEDFileEntities *entities=0);
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
    MEDLOADER_EXPORT int getProfileId(const std::string& pfl) const;
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
    MCAuto< MEDFileFieldGlobs > _globals;
  };

  class MEDFileFieldNameScope
  {
  public:
    MEDLOADER_EXPORT MEDFileFieldNameScope();
    MEDLOADER_EXPORT MEDFileFieldNameScope(const std::string& fieldName, const std::string& meshName);
    MEDLOADER_EXPORT std::string getName() const;
    MEDLOADER_EXPORT void setName(const std::string& fieldName);
    MEDLOADER_EXPORT std::string getDtUnit() const;
    MEDLOADER_EXPORT void setDtUnit(const std::string& dtUnit);
    MEDLOADER_EXPORT void copyNameScope(const MEDFileFieldNameScope& other);
    MEDLOADER_EXPORT std::string getMeshName() const;
    MEDLOADER_EXPORT void setMeshName(const std::string& meshName);
  protected:
    std::string _name;
    std::string _dt_unit;
    std::string _mesh_name;
  };
}

#endif
