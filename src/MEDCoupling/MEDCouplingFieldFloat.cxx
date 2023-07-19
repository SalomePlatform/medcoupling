// Copyright (C) 2007-2023  CEA, EDF
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

#include "MEDCouplingFieldFloat.hxx"
#include "MEDCouplingFieldT.txx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldInt.hxx"
#include "MEDCouplingFieldInt64.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingMesh.hxx"
#include "MEDCouplingMemArray.txx"

using namespace MEDCoupling;

template class MEDCoupling::MEDCouplingFieldT<float>;

MEDCouplingFieldFloat *MEDCouplingFieldFloat::New(TypeOfField type, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldFloat(type,td);
}

MEDCouplingFieldFloat *MEDCouplingFieldFloat::New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldFloat(ft,td);
}

MEDCouplingFieldFloat::MEDCouplingFieldFloat(TypeOfField type, TypeOfTimeDiscretization td):MEDCouplingFieldT<float>(type,MEDCouplingTimeDiscretizationFloat::New(td))
{
}

MEDCouplingFieldFloat::MEDCouplingFieldFloat(const MEDCouplingFieldFloat& other, bool deepCpy):MEDCouplingFieldT<float>(other,deepCpy)
{
}

MEDCouplingFieldFloat::MEDCouplingFieldFloat(NatureOfField n, MEDCouplingTimeDiscretizationFloat *td, MEDCouplingFieldDiscretization *type):MEDCouplingFieldT<float>(type,n,td)
{
}

/*!
 * ** WARINING : This method do not deeply copy neither mesh nor spatial discretization. Only a shallow copy (reference) is done for mesh and spatial discretization ! **
 */
MEDCouplingFieldFloat::MEDCouplingFieldFloat(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td):MEDCouplingFieldT<float>(ft,MEDCouplingTimeDiscretizationFloat::New(td),false)
{
}

MEDCouplingFieldFloat *MEDCouplingFieldFloat::deepCopy() const
{
  return cloneWithMesh(true);
}

MEDCouplingFieldFloat *MEDCouplingFieldFloat::clone(bool recDeepCpy) const
{
  return new MEDCouplingFieldFloat(*this,recDeepCpy);
}

template<class U>
typename Traits<U>::FieldType *ConvertToUField(const MEDCouplingFieldFloat *self)
{
  MCAuto<MEDCouplingFieldTemplate> tmp(MEDCouplingFieldTemplate::New(*self));
  int t1,t2;
  double t0(self->getTime(t1,t2));
  MCAuto<typename Traits<U>::FieldType> ret(Traits<U>::FieldType::New(*tmp,self->getTimeDiscretization()));
  ret->setTime(t0,t1,t2);
  if(self->getArray())
    {
      MCAuto<typename Traits<U>::ArrayType> arr(self->getArray()->convertToOtherTypeOfArr<U>());
      ret->setArray(arr);
    }
  return ret.retn();
}

MEDCouplingFieldDouble * MEDCouplingFieldFloat::convertToDblField() const
{
  return ConvertToUField<double>(this);
}

MEDCouplingFieldInt32 * MEDCouplingFieldFloat::convertToIntField() const
{
  return ConvertToUField<int>(this);
}

MEDCouplingFieldInt64 * MEDCouplingFieldFloat::convertToInt64Field() const
{
  return ConvertToUField<Int64>(this);
}
