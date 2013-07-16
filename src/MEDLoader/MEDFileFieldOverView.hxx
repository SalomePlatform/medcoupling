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

#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDCouplingRefCountObject.hxx"

#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"

#include <vector>

namespace ParaMEDMEM
{
  class DataArrayInt;
  class MEDFileMesh;
  class MEDFileFieldGlobs;
  class MEDFileAnyTypeField1TS;
  class MEDFileAnyTypeFieldMultiTS;

  class MEDFileMeshStruct : public RefCountObject
  {
  public:
    static MEDFileMeshStruct *New(const MEDFileMesh *mesh);
    std::size_t getHeapMemorySize() const;
    int getNumberOfNodes() const { return _nb_nodes; }
    int getNumberOfElemsOfGeoType(INTERP_KERNEL::NormalizedCellType t) const throw(INTERP_KERNEL::Exception);
  private:
    MEDFileMeshStruct(const MEDFileMesh *mesh);
  private:
    const MEDFileMesh *_mesh;
    std::string _name;
    int _nb_nodes;
    std::vector< std::vector<int> > _geo_types_distrib;
  };

  class MEDFileField1TSStructItem2
  {
  public:
    MEDFileField1TSStructItem2(INTERP_KERNEL::NormalizedCellType a, const std::pair<int,int>& b, const std::string& pfl, const std::string& loc);
    void checkWithMeshStructForCells(MEDFileMeshStruct *mst, const MEDFileFieldGlobs *globs) throw(INTERP_KERNEL::Exception);
    void checkWithMeshStructForGaussNE(MEDFileMeshStruct *mst, const MEDFileFieldGlobs *globs) throw(INTERP_KERNEL::Exception);
    void checkWithMeshStructForGaussPT(MEDFileMeshStruct *mst, const MEDFileFieldGlobs *globs) throw(INTERP_KERNEL::Exception);
    //
    INTERP_KERNEL::NormalizedCellType getGeo() const { return _geo_type; }
    //! warning this method also set _nb_of_entity attribute !
    void checkInRange(int nbOfEntity, int nip, const MEDFileFieldGlobs *globs) throw(INTERP_KERNEL::Exception);
    bool operator==(const MEDFileField1TSStructItem2& other) const throw(INTERP_KERNEL::Exception);
  private:
    INTERP_KERNEL::NormalizedCellType _geo_type;
    std::pair<int,int> _start_end;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _pfl;
    std::string _loc;
    int _nb_of_entity;
  };

  class MEDFileField1TSStructItem
  {
  public:
    MEDFileField1TSStructItem(TypeOfField a, const std::vector< MEDFileField1TSStructItem2 >& b);
    void checkWithMeshStruct(MEDFileMeshStruct *mst, const MEDFileFieldGlobs *globs) throw(INTERP_KERNEL::Exception);
    bool operator==(const MEDFileField1TSStructItem& other) const throw(INTERP_KERNEL::Exception);
    bool isEntityCell() const;
    //
    MEDFileField1TSStructItem simplifyMeOnCellEntity(const MEDFileFieldGlobs *globs) const throw(INTERP_KERNEL::Exception);
  private:
    TypeOfField _type;
    std::vector< MEDFileField1TSStructItem2 > _items;
  };

  class MEDFileField1TSStruct : public RefCountObject
  {
  public:
    static MEDFileField1TSStruct *New(const MEDFileAnyTypeField1TS *ref) throw(INTERP_KERNEL::Exception);
    void checkWithMeshStruct(MEDFileMeshStruct *mst, const MEDFileFieldGlobs *globs) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    bool isEqualConsideringThePast(const MEDFileAnyTypeField1TS *other) throw(INTERP_KERNEL::Exception);
    bool isSupportSameAs(const MEDFileAnyTypeField1TS *other) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileField1TSStruct(const MEDFileAnyTypeField1TS *ref);
    static MEDFileField1TSStructItem BuildItemFrom(const MEDFileAnyTypeField1TS *ref);
  private:
    std::vector<MEDFileField1TSStructItem> _already_checked;
  };

  class MEDFileFastCellSupportComparator : public RefCountObject
  {
  public:
    static MEDFileFastCellSupportComparator *New(const MEDFileMesh *m, const MEDFileAnyTypeFieldMultiTS *ref) throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDFileAnyTypeFieldMultiTS *other) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
  private:
    MEDFileFastCellSupportComparator(const MEDFileMesh *m, const MEDFileAnyTypeFieldMultiTS *ref);
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDFileMeshStruct> _mesh_comp;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSStruct> > _f1ts_cmps;
  };
}

#endif
