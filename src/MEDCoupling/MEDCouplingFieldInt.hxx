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
// Author : Yann Pora (EDF R&D)

#ifndef __PARAMEDMEM_MEDCOUPLINGFIELDINT_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELDINT_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingField.hxx"
#include "MEDCouplingTimeDiscretization.hxx"
#include "MEDCouplingMemArray.hxx"

#include <string>

namespace MEDCoupling
{
  class MEDCouplingFieldInt : public MEDCouplingField
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingFieldInt *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    MEDCOUPLING_EXPORT void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT std::string simpleRepr() const;
    MEDCOUPLING_EXPORT void reprQuickOverview(std::ostream& stream) const;
    MEDCOUPLING_EXPORT void setTimeUnit(const std::string& unit);
    MEDCOUPLING_EXPORT std::string getTimeUnit() const;
    MEDCOUPLING_EXPORT void setTime(double val, int iteration, int order);
    MEDCOUPLING_EXPORT double getTime(int& iteration, int& order) const;
    MEDCOUPLING_EXPORT void setArray(DataArrayInt *array);
    MEDCOUPLING_EXPORT const DataArrayInt *getArray() const;
    MEDCOUPLING_EXPORT DataArrayInt *getArray();
  protected:
    MEDCouplingFieldInt(TypeOfField type, TypeOfTimeDiscretization td);
    MEDCouplingFieldInt(const MEDCouplingFieldInt& other, bool deepCopy);
    MEDCouplingFieldInt(NatureOfField n, MEDCouplingTimeDiscretization *td, MEDCouplingFieldDiscretization *type);
    ~MEDCouplingFieldInt();
  private:
    MEDCouplingTimeDiscretization *_time_discr;
    MCAuto<DataArrayInt> _array;// agy : don't panic ! this is temporary ! templatization of time discr is planned !
  };
}

#endif
