// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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
#include "MEDCouplingCartesianAMRMesh.hxx"

namespace ParaMEDMEM
{
  /// @cond INTERNAL
  class DataArrayDoubleCollection : public RefCountObject, public TimeLabel
  {
  public:
    static DataArrayDoubleCollection *New(const std::vector< std::pair<std::string,int> >& fieldNames);
    void allocTuples(int nbOfTuples);
    void dellocTuples();
    void spillInfoOnComponents(const std::vector< std::vector<std::string> >& compNames);
    std::vector<DataArrayDouble *> retrieveFields() const;
    DataArrayDouble *retrieveFieldWithName(const std::string& name) const;
    static void SynchronizeFineToCoarse(int ghostLev, const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh, int patchId, const DataArrayDoubleCollection *fine, DataArrayDoubleCollection *coarse);
    static void SynchronizeCoarseToFine(int ghostLev, const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh, int patchId, const DataArrayDoubleCollection *coarse, DataArrayDoubleCollection *fine);
    static void SynchronizeFineEachOther(int patchId, int ghostLev, const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh, const std::vector<const MEDCouplingCartesianAMRMeshGen *>& children, const std::vector<DataArrayDoubleCollection *>& fieldsOnFine);
    static void SynchronizeCoarseToFineOnlyInGhostZone(int ghostLev, const MEDCouplingCartesianAMRMeshGen *fatherOfFineMesh, int patchId, const DataArrayDoubleCollection *coarse, DataArrayDoubleCollection *fine);
  private:
    DataArrayDoubleCollection(const std::vector< std::pair<std::string,int> >& fieldNames);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
    void updateTime() const;
    static void CheckDiscriminantNames(const std::vector<std::string>& names);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > _arrs;
  };

  class MEDCouplingGridCollection : public RefCountObject, public TimeLabel
  {
  public:
    static MEDCouplingGridCollection *New(const std::vector<const MEDCouplingCartesianAMRMeshGen *>& ms, const std::vector< std::pair<std::string,int> >& fieldNames);
    void alloc(int ghostLev);
    void dealloc();
    void spillInfoOnComponents(const std::vector< std::vector<std::string> >& compNames);
    bool presenceOf(const MEDCouplingCartesianAMRMeshGen *m, int& pos) const;
    const DataArrayDoubleCollection& retrieveFieldsAt(int pos) const;
    static void SynchronizeFineToCoarse(int ghostLev, const MEDCouplingGridCollection *fine, const MEDCouplingGridCollection *coarse);
    static void SynchronizeCoarseToFine(int ghostLev, const MEDCouplingGridCollection *coarse, const MEDCouplingGridCollection *fine);
    void synchronizeFineEachOther(int ghostLev) const;
    static void SynchronizeCoarseToFineOnlyInGhostZone(int ghostLev, const MEDCouplingGridCollection *coarse, const MEDCouplingGridCollection *fine);
  private:
    MEDCouplingGridCollection(const std::vector<const MEDCouplingCartesianAMRMeshGen *>& ms, const std::vector< std::pair<std::string,int> >& fieldNames);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
    void updateTime() const;
  private:
    std::vector< std::pair<const MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > > _map_of_dadc;
  };
  /// @endcond

  class MEDCouplingAMRAttribute : public MEDCouplingDataForGodFather, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingAMRAttribute *New(MEDCouplingCartesianAMRMesh *gf, const std::vector< std::pair<std::string,int> >& fieldNames);
    MEDCOUPLING_EXPORT static MEDCouplingAMRAttribute *New(MEDCouplingCartesianAMRMesh *gf, const std::vector< std::pair<std::string, std::vector<std::string> > >& fieldNames);
    MEDCOUPLING_EXPORT void spillInfoOnComponents(const std::vector< std::vector<std::string> >& compNames);
    MEDCOUPLING_EXPORT std::vector<DataArrayDouble *> retrieveFieldsOn(MEDCouplingCartesianAMRMeshGen *mesh) const;
    MEDCOUPLING_EXPORT DataArrayDouble *retrieveFieldOn(MEDCouplingCartesianAMRMeshGen *mesh, const std::string& fieldName) const;
    //
    MEDCOUPLING_EXPORT void synchronizeFineToCoarse(int ghostLev);
    MEDCOUPLING_EXPORT void synchronizeCoarseToFine(int ghostLev);
    MEDCOUPLING_EXPORT void synchronizeCoarseToFineOnlyInGhostZone(int ghostLev);
    MEDCOUPLING_EXPORT void synchronizeFineEachOtherInGhostZone(int ghostLev);
    MEDCOUPLING_EXPORT void alloc(int ghostLev);
    MEDCOUPLING_EXPORT void dealloc();
    MEDCOUPLING_EXPORT bool changeGodFather(MEDCouplingCartesianAMRMesh *gf);
    //
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDCOUPLING_EXPORT void updateTime() const;
  private:
    MEDCouplingAMRAttribute(MEDCouplingCartesianAMRMesh *gf, const std::vector< std::pair<std::string,int> >& fieldNames);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingGridCollection> > _levs;
  };
}

#endif
