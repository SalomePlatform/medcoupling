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

#ifndef __PARAMEDMEM_MEDCOUPLINGFIELDOVERTIME_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELDOVERTIME_HXX__

#include "MEDCouplingMultiFields.hxx"
#include "MEDCouplingDefinitionTime.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include <vector>

namespace MEDCoupling
{
  class MEDCouplingFieldOverTime : public MEDCouplingMultiFields
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingFieldOverTime *New(const std::vector<MEDCouplingFieldDouble *>& fs);
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT double getTimeTolerance() const;
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT bool isEqual(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStr(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    //void getIdsToFetch(double time, int& fieldId, int& arrId, int& meshId) const;
    //void setFieldOnId(int fieldId, MEDCouplingFieldDouble *f);
    //void dispatchPointers();
    MEDCOUPLING_EXPORT std::vector<MEDCouplingMesh *> getMeshes() const;
    MEDCOUPLING_EXPORT std::vector<MEDCouplingMesh *> getDifferentMeshes(std::vector<int>& refs) const;
    MEDCOUPLING_EXPORT std::vector<DataArrayDouble *> getArrays() const;
    MEDCOUPLING_EXPORT std::vector<DataArrayDouble *> getDifferentArrays(std::vector< std::vector<int> >& refs) const;
    MEDCOUPLING_EXPORT MEDCouplingDefinitionTime getDefinitionTimeZone() const;
  protected:
    MEDCOUPLING_EXPORT MEDCouplingFieldOverTime();
  private:
    MEDCouplingFieldOverTime(const std::vector<MEDCouplingFieldDouble *>& fs);
  };
}

#endif
