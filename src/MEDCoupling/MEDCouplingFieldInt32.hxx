// Copyright (C) 2007-2026  CEA, EDF
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

#pragma once

#include "MEDCoupling.hxx"
#include "MEDCouplingFieldT.hxx"
#include "MEDCouplingMemArray.hxx"

#include <string>

namespace MEDCoupling
{
class MEDCouplingFieldDouble;
class MEDCouplingFieldTemplate;

class MEDCouplingFieldInt32 : public MEDCouplingFieldT<Int32>
{
   public:
    MEDCOUPLING_EXPORT static MEDCouplingFieldInt32 *New(TypeOfField type, TypeOfTimeDiscretization td = ONE_TIME);
    MEDCOUPLING_EXPORT static MEDCouplingFieldInt32 *New(
        const MEDCouplingFieldTemplate &ft, TypeOfTimeDiscretization td = ONE_TIME
    );
    MEDCOUPLING_EXPORT MEDCouplingFieldInt32 *deepCopy() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldInt32 *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *convertToDblField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldFloat *convertToFloatField() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldInt64 *convertToInt64Field() const;
    std::string getClassName() const override { return std::string("MEDCouplingFieldInt32"); }

   protected:
    MEDCouplingFieldInt32(TypeOfField type, TypeOfTimeDiscretization td);
    MEDCouplingFieldInt32(const MEDCouplingFieldInt32 &other, bool deepCopy);
    MEDCouplingFieldInt32(
        NatureOfField n, MEDCouplingTimeDiscretizationInt32 *td, MEDCouplingFieldDiscretization *type
    );
    MEDCouplingFieldInt32(const MEDCouplingFieldTemplate &ft, TypeOfTimeDiscretization td);
    ~MEDCouplingFieldInt32() {}
};
}  // namespace MEDCoupling
