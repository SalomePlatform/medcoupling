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
  class DataArrayDoubleCollection : public RefCountObject
  {
  public:
    static DataArrayDoubleCollection *New(const std::vector< std::pair<std::string,int> >& fieldNames);
    void allocTuples(int nbOfTuples);
    void dellocTuples();
  private:
    DataArrayDoubleCollection(const std::vector< std::pair<std::string,int> >& fieldNames);
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
    static void CheckDiscriminantNames(const std::vector<std::string>& names);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> > _arrs;
  };

  class MEDCouplingGridCollection : public RefCountObject
  {
  private:
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
  private:
    std::vector< std::pair<MEDCouplingCartesianAMRMeshGen *,MEDCouplingAutoRefCountObjectPtr<DataArrayDoubleCollection> > > _map_of_dadc;
  };
  /// @endcond

  class MEDCouplingAMRAttribute : public MEDCouplingDataForGodFather, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingAMRAttribute *New(MEDCouplingCartesianAMRMesh *gf, const std::vector< std::pair<std::string,int> >& fieldNames);
    //
    MEDCOUPLING_EXPORT void alloc();
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
