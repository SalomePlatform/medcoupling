// Copyright (C) 2022-2025  CEA, EDF
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

#pragma once

#include "MEDCouplingFieldDiscretization.hxx"

#include <functional>

namespace MEDCoupling
{
  /*!
   * Class in charge to implement FE functions with shape functions
   */
  class MEDCouplingFieldDiscretizationOnNodesFE : public MEDCouplingFieldDiscretizationOnNodes
  {
    public:
      TypeOfField getEnum() const override { return TYPE; }
      std::string getClassName() const override { return std::string("MEDCouplingFieldDiscretizationOnNodesFE"); }
      MEDCOUPLING_EXPORT  const char *getRepr() const override;
      MEDCOUPLING_EXPORT std::string getStringRepr() const override;
      MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const override;
      MEDCOUPLING_EXPORT MCAuto<MEDCouplingFieldDiscretization> aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const override;
      MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const override;
      MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization *clone() const override;
      MEDCOUPLING_EXPORT void checkCompatibilityWithNature(NatureOfField nat) const override;
      MEDCOUPLING_EXPORT MEDCouplingFieldDouble *getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const override;
      MEDCOUPLING_EXPORT void getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const override;
      MEDCOUPLING_EXPORT DataArrayDouble *getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const override;
    public:
      MEDCOUPLING_EXPORT MCAuto<DataArrayDouble> getCooInRefElement(const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const;
    public:
      MEDCOUPLING_EXPORT static void GetRefCoordOfListOf3DPtsIn3D(const MEDCouplingUMesh *umesh, const double *ptsCoo, mcIdType nbOfPts,
      std::function<void(const MEDCouplingGaussLocalization&, const std::vector<mcIdType>&)> customFunc);
    private:
      const MEDCouplingUMesh *checkConfig3D(const MEDCouplingMesh *mesh) const;
    public:
      static const char REPR[];
      static constexpr TypeOfField TYPE = ON_NODES_FE;
  };
}
