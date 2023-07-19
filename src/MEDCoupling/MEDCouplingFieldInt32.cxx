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
// Author : Yann Pora (EDF R&D)

#include "MEDCouplingFieldInt32.hxx"
#include "MEDCouplingFieldInt64.hxx"
#include "MEDCouplingFieldT.txx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldFloat.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingMesh.hxx"
#include "MEDCouplingMemArray.txx"

using namespace MEDCoupling;

template class MEDCoupling::MEDCouplingFieldT<Int32>;

MEDCouplingFieldInt32 *MEDCouplingFieldInt32::New(TypeOfField type, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldInt32(type,td);
}

MEDCouplingFieldInt32 *MEDCouplingFieldInt32::New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldInt32(ft,td);
}

MEDCouplingFieldInt32::MEDCouplingFieldInt32(TypeOfField type, TypeOfTimeDiscretization td):MEDCouplingFieldT<Int32>(type,MEDCouplingTimeDiscretizationInt32::New(td))
{
}

MEDCouplingFieldInt32::MEDCouplingFieldInt32(const MEDCouplingFieldInt32& other, bool deepCpy):MEDCouplingFieldT<Int32>(other,deepCpy)
{
}

MEDCouplingFieldInt32::MEDCouplingFieldInt32(NatureOfField n, MEDCouplingTimeDiscretizationInt32 *td, MEDCouplingFieldDiscretization *type):MEDCouplingFieldT<Int32>(type,n,td)
{
}

/*!
 * ** WARINING : This method do not deeply copy neither mesh nor spatial discretization. Only a shallow copy (reference) is done for mesh and spatial discretization ! **
 */
MEDCouplingFieldInt32::MEDCouplingFieldInt32(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td):MEDCouplingFieldT<Int32>(ft,MEDCouplingTimeDiscretizationInt32::New(td),false)
{
}

MEDCouplingFieldInt32 *MEDCouplingFieldInt32::deepCopy() const
{
  return cloneWithMesh(true);
}

MEDCouplingFieldInt32 *MEDCouplingFieldInt32::clone(bool recDeepCpy) const
{
  return new MEDCouplingFieldInt32(*this,recDeepCpy);
}

template<class U>
typename Traits<U>::FieldType *ConvertToUField(const MEDCouplingFieldInt32 *self)
{
  MCAuto<MEDCouplingFieldTemplate> tmp(MEDCouplingFieldTemplate::New(*self));
  int t1,t2;
  double t0(self->getTime(t1,t2));
  MCAuto<typename Traits<U>::FieldType > ret(Traits<U>::FieldType::New(*tmp,self->getTimeDiscretization()));
  ret->setTime(t0,t1,t2);
  if(self->getArray())
    {
      MCAuto<typename Traits<U>::ArrayType> arr(self->getArray()->convertToOtherTypeOfArr<U>());
      ret->setArray(arr);
    }
  return ret.retn();
}

MEDCouplingFieldDouble *MEDCouplingFieldInt32::convertToDblField() const
{
  return ConvertToUField<double>(this);
}

MEDCouplingFieldInt64 *MEDCouplingFieldInt32::convertToInt64Field() const
{
  return ConvertToUField<Int64>(this);
}

MEDCouplingFieldFloat *MEDCouplingFieldInt32::convertToFloatField() const
{
  return ConvertToUField<float>(this);
}