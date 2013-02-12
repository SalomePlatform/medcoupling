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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGSTRUCTUREDMESH_HXX__
#define __PARAMEDMEM_MEDCOUPLINGSTRUCTUREDMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingMesh.hxx"

namespace ParaMEDMEM
{
  class MEDCOUPLING_EXPORT MEDCouplingStructuredMesh : public MEDCouplingMesh
  {
  public:
    INTERP_KERNEL::NormalizedCellType getTypeOfCell(int cellId) const;
    std::set<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const;
    int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const;
    static void GetPosFromId(int nodeId, int meshDim, const int *split, int *res);
    void getNodeIdsOfCell(int cellId, std::vector<int>& conn) const;
    std::size_t getHeapMemorySize() const;
    void copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception);
    bool isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception);
    //tools
    std::vector<int> getDistributionOfTypes() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception);
    void splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *buildUnstructured() const throw(INTERP_KERNEL::Exception);
    MEDCouplingMesh *buildPart(const int *start, const int *end) const;
    MEDCouplingMesh *buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const;
    DataArrayInt *simplexize(int policy) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildOrthogonalField() const;
    void fill1DUnstructuredMesh(MEDCouplingUMesh *m) const;
    void fill2DUnstructuredMesh(MEDCouplingUMesh *m) const;
    void fill3DUnstructuredMesh(MEDCouplingUMesh *m) const;
    //some useful methods
    int getCellIdFromPos(int i, int j, int k) const;
    int getNodeIdFromPos(int i, int j, int k) const;
    virtual void getNodeGridStructure(int *res) const = 0;
    virtual void getSplitCellValues(int *res) const = 0;
    virtual void getSplitNodeValues(int *res) const = 0;
  protected:
    MEDCouplingStructuredMesh();
    MEDCouplingStructuredMesh(const MEDCouplingStructuredMesh& other, bool deepCpy);
    ~MEDCouplingStructuredMesh();
  };
}

#endif
