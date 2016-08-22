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

#ifndef __PARAMEDMEM_MEDCOUPLINGFIELDTEMPLATE_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELDTEMPLATE_HXX__

#include "MEDCouplingField.hxx"

namespace MEDCoupling
{
  class MEDCouplingFieldInt;
  class MEDCouplingFieldDouble;
  /*!
   * \brief A field template can be seen as a field without the array of values.
   *
   * A field template aggregates a MEDCouplingMesh and a spatial discretization object (instance of
   * MEDCouplingFieldDiscretization).
   * 
   * MEDCouplingFieldTemplate is the most appropriate type for the preparation of matrix using
   * MEDCouplingRemapper::prepareEx, since it contains the minimal information requireds to prepare
   * the interpolation matrix.
   */
  class MEDCouplingFieldTemplate : public MEDCouplingField
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingFieldTemplate *New(const MEDCouplingFieldDouble& f);
    MEDCOUPLING_EXPORT static MEDCouplingFieldTemplate *New(const MEDCouplingFieldInt& f);
    MEDCOUPLING_EXPORT static MEDCouplingFieldTemplate *New(TypeOfField type);
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT std::string advancedRepr() const;
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    //
    MEDCOUPLING_EXPORT void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    MEDCOUPLING_EXPORT void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    MEDCOUPLING_EXPORT void resizeForUnserialization(const std::vector<int>& tinyInfoI, DataArrayInt *&dataInt);
    MEDCOUPLING_EXPORT void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    MEDCOUPLING_EXPORT void serialize(DataArrayInt *&dataInt) const;
    //
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
  private:
    MEDCouplingFieldTemplate(const MEDCouplingFieldDouble& f);
    MEDCouplingFieldTemplate(const MEDCouplingFieldInt& f);
    MEDCouplingFieldTemplate(TypeOfField type);
  };
}

#endif
