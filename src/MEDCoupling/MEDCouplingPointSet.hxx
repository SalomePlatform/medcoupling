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
    static MEDCouplingPointSet *buildInstanceFromMeshType(MEDCouplingMeshType type);
    virtual MEDCouplingPointSet *buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const = 0;
    //! size of returned tinyInfo must be always the same.
    virtual void getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const;
    virtual void resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings);
    virtual void serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const;
    virtual void unserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2,
                                 const std::vector<std::string>& littleStrings);
    virtual void giveElemsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems) = 0;
  protected:
    static bool intersectsBoundingBox(const double* bb1, const double* bb2, int dim, double eps);
  protected:
    DataArrayDouble *_coords;
  };
}

#endif
