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
#ifndef __PARAMEDMEM_MEDCOUPLINGPOINTSET_HXX__
#define __PARAMEDMEM_MEDCOUPLINGPOINTSET_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingMesh.hxx"

#include <vector>

namespace ParaMEDMEM
{
  class DataArrayInt;
  class DataArrayDouble;

  class MEDCOUPLING_EXPORT MEDCouplingPointSet : public MEDCouplingMesh
  {
  protected:
    MEDCouplingPointSet();
    MEDCouplingPointSet(const MEDCouplingPointSet& other, bool deepCpy);
    ~MEDCouplingPointSet();
  public:
    void updateTime();
    bool isStructured() const;
    int getNumberOfNodes() const;
    int getSpaceDimension() const;
    void setCoords(DataArrayDouble *coords);
    DataArrayDouble *getCoords() const { return _coords; }
    bool areCoordsEqual(const MEDCouplingPointSet& other, double prec) const;
    void getBoundingBox(double *bbox) const;
    void zipCoords();
    void rotate(const double *center, const double *vector, double angle);
    void translate(const double *vector);
    static MEDCouplingPointSet *buildInstanceFromMeshType(MEDCouplingMeshType type);
    virtual MEDCouplingPointSet *buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const = 0;
    virtual MEDCouplingPointSet *buildPartOfMySelfNode(const int *start, const int *end, bool fullyIn) const = 0;
    virtual void findBoundaryNodes(std::vector<int>& nodes) const = 0;
    virtual MEDCouplingPointSet *buildBoundaryMesh(bool keepCoords) const = 0;
    virtual void renumberConnectivity(const int *newNodeNumbers) = 0;
    //! size of returned tinyInfo must be always the same.
    virtual void getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    virtual void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings);
    virtual void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    virtual void unserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2,
                                 const std::vector<std::string>& littleStrings);
    virtual void giveElemsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems) = 0;
    virtual DataArrayInt *zipCoordsTraducer() = 0;
  protected:
    virtual void checkFullyDefined() const throw(INTERP_KERNEL::Exception) = 0;
    static bool intersectsBoundingBox(const double* bb1, const double* bb2, int dim, double eps);
    void rotate2D(const double *center, double angle);
    void rotate3D(const double *center, const double *vect, double angle);
  protected:
    DataArrayDouble *_coords;
  };
}

#endif
