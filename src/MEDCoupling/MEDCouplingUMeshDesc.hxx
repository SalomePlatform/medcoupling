//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef __PARAMEDMEM_MEDCOUPLINGUMESHDESC_HXX__
#define __PARAMEDMEM_MEDCOUPLINGUMESHDESC_HXX__

#include "MEDCouplingPointSet.hxx"
#include "MEDCoupling.hxx"
#include "NormalizedUnstructuredMesh.hxx"

#include <set>

namespace ParaMEDMEM
{
  class MEDCOUPLING_EXPORT MEDCouplingUMeshDesc : public MEDCouplingPointSet
  {
  public:
    static MEDCouplingUMeshDesc *New();
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    void setMeshDimension(unsigned meshDim);
    int getNumberOfCells() const;
    int getNumberOfFaces() const;
    int getCellMeshLength() const;
    int getFaceMeshLength() const;
    int getMeshDimension() const { return _mesh_dim; }
    void setConnectivity(DataArrayInt *descConn, DataArrayInt *descConnIndex, DataArrayInt *nodalFaceConn, DataArrayInt *nodalFaceConnIndx);
    //tools to overload
    void getTinySerializationInformation(std::vector<int>& tinyInfo) const;
    void resizeForSerialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2);
    void serialize(DataArrayInt *&a1, DataArrayDouble *&a2);
    MEDCouplingPointSet *buildObjectFromUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2);
    void giveElemsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems);
    MEDCouplingPointSet *buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const;
    MEDCouplingFieldDouble *getMeasureField() const;
  private:
    MEDCouplingUMeshDesc();
    ~MEDCouplingUMeshDesc();
    void computeTypes();
  private:
    unsigned _mesh_dim;
    DataArrayInt *_desc_connec;
    DataArrayInt *_desc_connec_index;
    DataArrayInt *_nodal_connec_face;
    DataArrayInt *_nodal_connec_face_index;
    std::set<INTERP_KERNEL::NormalizedCellType> _types;
  };
}

#endif
