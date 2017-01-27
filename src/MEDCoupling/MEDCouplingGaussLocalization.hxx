// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGGAUSSLOCALIZATION_HXX__
#define __PARAMEDMEM_MEDCOUPLINGGAUSSLOCALIZATION_HXX__

#include "MEDCoupling.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "MEDCouplingMemArray.hxx"
#include "InterpKernelException.hxx"

#include <vector>

namespace MEDCoupling
{
  class MEDCouplingMesh;
  class MEDCouplingUMesh;

  class MEDCouplingGaussLocalization
  {
  public:
    MEDCOUPLING_EXPORT MEDCouplingGaussLocalization(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                    const std::vector<double>& gsCoo, const std::vector<double>& w);
    MEDCOUPLING_EXPORT MEDCouplingGaussLocalization(INTERP_KERNEL::NormalizedCellType typ);
    MEDCOUPLING_EXPORT INTERP_KERNEL::NormalizedCellType getType() const { return _type; }
    MEDCOUPLING_EXPORT void setType(INTERP_KERNEL::NormalizedCellType typ);
    MEDCOUPLING_EXPORT int getNumberOfGaussPt() const { return (int)_weight.size(); }
    MEDCOUPLING_EXPORT int getDimension() const;
    MEDCOUPLING_EXPORT int getNumberOfPtsInRefCell() const;
    MEDCOUPLING_EXPORT std::string getStringRepr() const;
    MEDCOUPLING_EXPORT std::size_t getMemorySize() const;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT bool isEqual(const MEDCouplingGaussLocalization& other, double eps) const;
    MEDCOUPLING_EXPORT void pushTinySerializationIntInfo(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void pushTinySerializationDblInfo(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT const double *fillWithValues(const double *vals);
    //
    MEDCOUPLING_EXPORT MCAuto<DataArrayDouble> localizePtsInRefCooForEachCell(const DataArrayDouble *ptsInRefCoo, const MEDCouplingUMesh *mesh) const;
    MEDCOUPLING_EXPORT MCAuto<MEDCouplingUMesh> buildRefCell() const;
    //
    MEDCOUPLING_EXPORT const std::vector<double>& getRefCoords() const { return _ref_coord; }
    MEDCOUPLING_EXPORT double getRefCoord(int ptIdInCell, int comp) const;
    MEDCOUPLING_EXPORT const std::vector<double>& getGaussCoords() const { return _gauss_coord; }
    MEDCOUPLING_EXPORT double getGaussCoord(int gaussPtIdInCell, int comp) const;
    MEDCOUPLING_EXPORT const std::vector<double>& getWeights() const { return _weight; }
    MEDCOUPLING_EXPORT double getWeight(int gaussPtIdInCell, double newVal) const;
    MEDCOUPLING_EXPORT void setRefCoord(int ptIdInCell, int comp, double newVal);
    MEDCOUPLING_EXPORT void setGaussCoord(int gaussPtIdInCell, int comp, double newVal);
    MEDCOUPLING_EXPORT void setWeight(int gaussPtIdInCell, double newVal);
    MEDCOUPLING_EXPORT void setRefCoords(const std::vector<double>& refCoo);
    MEDCOUPLING_EXPORT void setGaussCoords(const std::vector<double>& gsCoo);
    MEDCOUPLING_EXPORT void setWeights(const std::vector<double>& w);
    //
    MEDCOUPLING_EXPORT static MEDCouplingGaussLocalization BuildNewInstanceFromTinyInfo(int dim, const std::vector<int>& tinyData);
    MEDCOUPLING_EXPORT static bool AreAlmostEqual(const std::vector<double>& v1, const std::vector<double>& v2, double eps);
  private:
    int checkCoherencyOfRequest(int gaussPtIdInCell, int comp) const;
  private:
    INTERP_KERNEL::NormalizedCellType _type;
    std::vector<double> _ref_coord;
    std::vector<double> _gauss_coord;
    std::vector<double> _weight;
  };
}

#endif
