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
// Author : Anthony Geay

#ifndef __MEDCOUPLINGAMRATTRIBUTE_HXX__
#define __MEDCOUPLINGAMRATTRIBUTE_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingNatureOfFieldEnum"
#include "MEDCouplingCartesianAMRMesh.hxx"

namespace MEDCoupling
{
  /// @cond INTERNAL
  class DataArrayDoubleCollection : public RefCountObject, public TimeLabel
  {
  public:
    static DataArrayDoubleCollection *New(const std::vector< std::pair<std::string,int> >& fieldNames);
    std::string getClassName() const override { return std::string("DataArrayDoubleCollection"); }
    DataArrayDoubleCollection *deepCopy() const;
    void allocTuples(mcIdType nbOfTuples);
    void dellocTuples();
    void copyFrom(const DataArrayDoubleCollection& other);
    void spillInfoOnComponents(const std::vector< std::vector<std::string> >& compNames);
    void spillNatures(const std::vector<NatureOfField>& nfs);
    std::vector< std::pair<std::string, std::vector<std::string> > > getInfoOnComponents() const;
    std::vector<NatureOfField> getNatures() const;
    std::vector<DataArrayDouble *> retrieveFields() const;
    const DataArrayDouble *getFieldWithName(const std::string& name) const;
    DataArrayDouble *getFieldWithName(const std::string& name);
    DataArrayDouble *at(mcIdType pos);
    const DataArrayDouble *at(mcIdType pos) const;
    mcIdType size() const;
    static void SynchronizeFineToCoarse(mcIdType ghostLev, const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh, mcIdType patchId, const DataArrayDoubleCollection *fine, DataArrayDoubleCollection *coarse);
    static void SynchronizeCoarseToFine(mcIdType ghostLev, const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh, mcIdType patchId, const DataArrayDoubleCollection *coarse, DataArrayDoubleCollection *fine);
    static void SynchronizeFineEachOther(mcIdType patchId, mcIdType ghostLev, const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh, const std::vector<const MEDCouplingCartesianAMRMeshGen *>& children, const std::vector<DataArrayDoubleCollection *>& fieldsOnFine);
    static void SynchronizeCoarseToFineOnlyInGhostZone(mcIdType ghostLev, const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh, mcIdType patchId, const DataArrayDoubleCollection *coarse, DataArrayDoubleCollection *fine);
    static void SynchronizeGhostZoneOfOneUsingTwo(mcIdType ghostLev, const MEDCouplingCartesianAMRPatch *p1, const DataArrayDoubleCollection *p1dac, const MEDCouplingCartesianAMRPatch *p2, const DataArrayDoubleCollection *p2dac);
    void synchronizeMyGhostZoneUsing(mcIdType ghostLev, const DataArrayDoubleCollection& other, const MEDCouplingCartesianAMRPatch *thisp, const MEDCouplingCartesianAMRPatch *otherp, const MEDCouplingCartesianAMRMeshGen *father) const;
    void synchronizeMyGhostZoneUsingExt(mcIdType ghostLev, const DataArrayDoubleCollection& other, const MEDCouplingCartesianAMRPatch *thisp, const MEDCouplingCartesianAMRPatch *otherp) const;
  private:
    DataArrayDoubleCollection(const std::vector< std::pair<std::string,int> >& fieldNames);
    DataArrayDoubleCollection(const DataArrayDoubleCollection& other);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    void updateTime() const;
    static void CheckDiscriminantNames(const std::vector<std::string>& names);
    static bool IsConservativeNature(NatureOfField n);
    static void CheckSameNatures(NatureOfField n1, NatureOfField n2);
    static void CheckValidNature(NatureOfField n);
  private:
    std::vector< std::pair< MCAuto<DataArrayDouble>, NatureOfField > > _arrs;
  };

  class MEDCouplingGridCollection : public RefCountObject, public TimeLabel
  {
  public:
    static MEDCouplingGridCollection *New(const std::vector<const MEDCouplingCartesianAMRMeshGen *>& ms, const std::vector< std::pair<std::string,int> >& fieldNames);
    std::string getClassName() const override { return std::string("MEDCouplingGridCollection"); }
    MEDCouplingGridCollection *deepCopy(const MEDCouplingCartesianAMRMeshGen *newGf, const MEDCouplingCartesianAMRMeshGen *oldGf) const;
    void alloc(mcIdType ghostLev);
    void dealloc();
    void spillInfoOnComponents(const std::vector< std::vector<std::string> >& compNames);
    void spillNatures(const std::vector<NatureOfField>& nfs);
    std::vector< std::pair<std::string, std::vector<std::string> > > getInfoOnComponents() const;
    std::vector<NatureOfField> getNatures() const;
    bool presenceOf(const MEDCouplingCartesianAMRMeshGen *m, mcIdType& pos) const;
    const DataArrayDoubleCollection& getFieldsAt(mcIdType pos) const;
    DataArrayDoubleCollection& getFieldsAt(mcIdType pos);
    void copyOverlappedZoneFrom(mcIdType ghostLev, const MEDCouplingGridCollection& other);
    static void SynchronizeFineToCoarse(mcIdType ghostLev, const MEDCouplingGridCollection *fine, const MEDCouplingGridCollection *coarse);
    static void SynchronizeCoarseToFine(mcIdType ghostLev, const MEDCouplingGridCollection *coarse, const MEDCouplingGridCollection *fine);
    void synchronizeFineEachOther(mcIdType ghostLev, const std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> >& ps) const;
    void synchronizeFineEachOtherExt(mcIdType ghostLev, const std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> >& ps) const;
    std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > findNeighbors(mcIdType ghostLev) const;
    static void SynchronizeCoarseToFineOnlyInGhostZone(mcIdType ghostLev, const MEDCouplingGridCollection *coarse, const MEDCouplingGridCollection *fine);
    void fillIfInTheProgenyOf(const std::string& fieldName, const MEDCouplingCartesianAMRMeshGen *head, std::vector<const DataArrayDouble *>& recurseArrs) const;
  private:
    MEDCouplingGridCollection(const std::vector<const MEDCouplingCartesianAMRMeshGen *>& ms, const std::vector< std::pair<std::string,int> >& fieldNames);
    MEDCouplingGridCollection(const MEDCouplingGridCollection& other, const MEDCouplingCartesianAMRMeshGen *newGf, const MEDCouplingCartesianAMRMeshGen *oldGf);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    void updateTime() const;
  private:
    std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MCAuto<DataArrayDoubleCollection> > > _map_of_dadc;
  };

  /// @endcond

  class MEDCouplingDataForGodFather : public RefCountObject
  {
    friend class MEDCouplingCartesianAMRMesh;
  public:
    MEDCOUPLING_EXPORT MEDCouplingCartesianAMRMesh *getMyGodFather();
    std::string getClassName() const override { return std::string("MEDCouplingDataForGodFather"); }
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMesh *getMyGodFather() const;
    MEDCOUPLING_EXPORT virtual void synchronizeFineToCoarse() = 0;
    MEDCOUPLING_EXPORT virtual void synchronizeFineToCoarseBetween(mcIdType fromLev, mcIdType toLev) = 0;
    MEDCOUPLING_EXPORT virtual void synchronizeCoarseToFine() = 0;
    MEDCOUPLING_EXPORT virtual void synchronizeCoarseToFineBetween(mcIdType fromLev, mcIdType toLev) = 0;
    MEDCOUPLING_EXPORT virtual void synchronizeAllGhostZones() = 0;
    MEDCOUPLING_EXPORT virtual void synchronizeAllGhostZonesOfDirectChidrenOf(const MEDCouplingCartesianAMRMeshGen *mesh) = 0;
    MEDCOUPLING_EXPORT virtual void synchronizeAllGhostZonesAtASpecifiedLevel(mcIdType level) = 0;
    MEDCOUPLING_EXPORT virtual void synchronizeAllGhostZonesAtASpecifiedLevelUsingOnlyFather(mcIdType level) = 0;
    MEDCOUPLING_EXPORT virtual void alloc() = 0;
    MEDCOUPLING_EXPORT virtual void dealloc() = 0;
  protected:
    MEDCouplingDataForGodFather(MEDCouplingCartesianAMRMesh *gf);
    void checkGodFatherFrozen() const;
  protected:
    virtual bool changeGodFather(MEDCouplingCartesianAMRMesh *gf);
    MEDCouplingDataForGodFather(const MEDCouplingDataForGodFather& other, bool deepCpyGF);
  protected:
    MCAuto<MEDCouplingCartesianAMRMesh> _gf;
    TimeLabelConstOverseer _tlc;
  };

  class MEDCouplingAMRAttribute : public MEDCouplingDataForGodFather, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingAMRAttribute *New(MEDCouplingCartesianAMRMesh *gf, const std::vector< std::pair<std::string,int> >& fieldNames, mcIdType ghostLev);
    MEDCOUPLING_EXPORT static MEDCouplingAMRAttribute *New(MEDCouplingCartesianAMRMesh *gf, const std::vector< std::pair<std::string, std::vector<std::string> > >& fieldNames, mcIdType ghostLev);
    std::string getClassName() const override { return std::string("MEDCouplingAMRAttribute"); }
    MEDCOUPLING_EXPORT void spillInfoOnComponents(const std::vector< std::vector<std::string> >& compNames);
    MEDCOUPLING_EXPORT void spillNatures(const std::vector<NatureOfField>& nfs);
    MEDCOUPLING_EXPORT MEDCouplingAMRAttribute *deepCopy() const;
    MEDCOUPLING_EXPORT MEDCouplingAMRAttribute *deepCpyWithoutGodFather() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfLevels() const;
    MEDCOUPLING_EXPORT std::vector<DataArrayDouble *> retrieveFieldsOn(MEDCouplingCartesianAMRMeshGen *mesh) const;
    MEDCOUPLING_EXPORT const DataArrayDouble *getFieldOn(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const;
    MEDCOUPLING_EXPORT DataArrayDouble *getFieldOn(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName);
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildCellFieldOnRecurseWithoutOverlapWithoutGhost(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildCellFieldOnWithGhost(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildCellFieldOnWithoutGhost(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const;
    MEDCOUPLING_EXPORT std::string writeVTHB(const std::string& fileName) const;
    //
    MEDCOUPLING_EXPORT MEDCouplingAMRAttribute *projectTo(MEDCouplingCartesianAMRMesh *targetGF) const;
    //
    MEDCOUPLING_EXPORT void synchronizeFineToCoarse();
    MEDCOUPLING_EXPORT void synchronizeFineToCoarseBetween(mcIdType fromLev, mcIdType toLev);
    MEDCOUPLING_EXPORT void synchronizeCoarseToFine();
    MEDCOUPLING_EXPORT void synchronizeCoarseToFineBetween(mcIdType fromLev, mcIdType toLev);
    MEDCOUPLING_EXPORT void synchronizeAllGhostZones();
    MEDCOUPLING_EXPORT void synchronizeAllGhostZonesOfDirectChidrenOf(const MEDCouplingCartesianAMRMeshGen *mesh);
    MEDCOUPLING_EXPORT void synchronizeAllGhostZonesAtASpecifiedLevel(mcIdType level);
    MEDCOUPLING_EXPORT void synchronizeAllGhostZonesAtASpecifiedLevelUsingOnlyFather(mcIdType level);
    //
    MEDCOUPLING_EXPORT void alloc();
    MEDCOUPLING_EXPORT void dealloc();
    MEDCOUPLING_EXPORT bool changeGodFather(MEDCouplingCartesianAMRMesh *gf);
    //
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT void updateTime() const;
  private:
    MEDCouplingAMRAttribute(MEDCouplingCartesianAMRMesh *gf, const std::vector< std::pair<std::string,int> >& fieldNames, mcIdType ghostLev);
    MEDCouplingAMRAttribute(const MEDCouplingAMRAttribute& other, bool deepCpyGF);
    const DataArrayDoubleCollection& findCollectionAttachedTo(const MEDCouplingCartesianAMRMeshGen *m) const;
    void synchronizeFineToCoarseByOneLevel(mcIdType level);
    void synchronizeCoarseToFineByOneLevel(mcIdType level);
  private:
    mcIdType _ghost_lev;
    std::vector< MCAuto<MEDCouplingGridCollection> > _levs;
    std::vector< std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > > _neighbors;
    std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > _mixed_lev_neighbors;
    std::vector< std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > > _cross_lev_neighbors;
  };
}

#endif
