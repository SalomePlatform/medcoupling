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

#ifndef __MEDCOUPLINGFIELDINT_HXX__
#define __MEDCOUPLINGFIELDINT_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingFieldT.hxx"
#include "MEDCouplingMemArray.hxx"

#include <string>

namespace MEDCoupling
{
  class MEDCouplingFieldDouble;
  class MEDCouplingFieldTemplate;
  
  class MEDCouplingFieldInt : public MEDCouplingFieldT<int>
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingFieldInt *New(TypeOfField type, TypeOfTimeDiscretization td=ONE_TIME);
    MEDCOUPLING_EXPORT static MEDCouplingFieldInt *New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td=ONE_TIME);
    MEDCOUPLING_EXPORT MEDCouplingFieldInt *deepCopy() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldInt *clone(bool recDeepCpy) const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *convertToDblField() const;
  protected:
    MEDCouplingFieldInt(TypeOfField type, TypeOfTimeDiscretization td);
    MEDCouplingFieldInt(const MEDCouplingFieldInt& other, bool deepCopy);
    MEDCouplingFieldInt(NatureOfField n, MEDCouplingTimeDiscretizationInt *td, MEDCouplingFieldDiscretization *type);
    MEDCouplingFieldInt(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td);
    ~MEDCouplingFieldInt() { }
  };
}

#endif
