// Copyright (C) 2020-2025  CEA, EDF
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

#include "MEDCoupling.hxx"
#include "MEDCouplingFieldT.hxx"
#include "MEDCouplingMemArray.hxx"

#include <string>

namespace MEDCoupling
{
  class MEDCouplingFieldDouble;
  class MEDCouplingFieldTemplate;
  
  class MEDCouplingFieldInt64 : public MEDCouplingFieldT<Int64>
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingFieldInt64 *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    MEDCOUPLING_EXPORT static MEDCouplingFieldInt64 *New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
    MEDCOUPLING_EXPORT MEDCouplingFieldInt64 *deepCopy() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldInt64 *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *convertToDblField() const;
    std::string getClassName() const override { return std::string("MEDCouplingFieldInt64"); }
  protected:
    MEDCouplingFieldInt64(TypeOfField type, TypeOfTimeDiscretization td);
    MEDCouplingFieldInt64(const MEDCouplingFieldInt64& other, bool deepCopy);
    MEDCouplingFieldInt64(NatureOfField n, MEDCouplingTimeDiscretizationInt64 *td, MEDCouplingFieldDiscretization *type);
    MEDCouplingFieldInt64(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td);
    ~MEDCouplingFieldInt64() { }
  };
}
