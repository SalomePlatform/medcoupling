// Copyright (C) 2007-2017  CEA/DEN, EDF R&D
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
// Author : Anthony GEAY (EDF R&D)

#ifndef __MEDCOUPLINGFIELDFLOAT_HXX__
#define __MEDCOUPLINGFIELDFLOAT_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingFieldT.hxx"
#include "MEDCouplingMemArray.hxx"

#include <string>

namespace MEDCoupling
{
  class MEDCouplingFieldDouble;
  class MEDCouplingFieldTemplate;
  
  class MEDCouplingFieldFloat : public MEDCouplingFieldT<float>
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingFieldFloat *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    MEDCOUPLING_EXPORT static MEDCouplingFieldFloat *New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
    MEDCOUPLING_EXPORT MEDCouplingFieldFloat *deepCopy() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldFloat *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *convertToDblField() const;
  protected:
    MEDCouplingFieldFloat(TypeOfField type, TypeOfTimeDiscretization td);
    MEDCouplingFieldFloat(const MEDCouplingFieldFloat& other, bool deepCpy);
    MEDCouplingFieldFloat(NatureOfField n, MEDCouplingTimeDiscretizationFloat *td, MEDCouplingFieldDiscretization *type);
    MEDCouplingFieldFloat(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td);
    ~MEDCouplingFieldFloat() { }
  };
}

#endif
