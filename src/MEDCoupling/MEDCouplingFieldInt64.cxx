// Copyright (C) 2020-2024  CEA, EDF
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

#include "MEDCouplingFieldInt64.hxx"
#include "MEDCouplingFieldT.txx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingMesh.hxx"
#include "MEDCouplingMemArray.txx"

using namespace MEDCoupling;

template class MEDCoupling::MEDCouplingFieldT<Int64>;

MEDCouplingFieldInt64 *MEDCouplingFieldInt64::New(TypeOfField type, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldInt64(type,td);
}

MEDCouplingFieldInt64 *MEDCouplingFieldInt64::New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldInt64(ft,td);
}

MEDCouplingFieldInt64::MEDCouplingFieldInt64(TypeOfField type, TypeOfTimeDiscretization td):MEDCouplingFieldT<Int64>(type,MEDCouplingTimeDiscretizationInt64::New(td))
{
}

MEDCouplingFieldInt64::MEDCouplingFieldInt64(const MEDCouplingFieldInt64& other, bool deepCpy):MEDCouplingFieldT<Int64>(other,deepCpy)
{
}

MEDCouplingFieldInt64::MEDCouplingFieldInt64(NatureOfField n, MEDCouplingTimeDiscretizationInt64 *td, MEDCouplingFieldDiscretization *type):MEDCouplingFieldT<Int64>(type,n,td)
{
}

/*!
 * ** WARINING : This method do not deeply copy neither mesh nor spatial discretization. Only a shallow copy (reference) is done for mesh and spatial discretization ! **
 */
MEDCouplingFieldInt64::MEDCouplingFieldInt64(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td):MEDCouplingFieldT<Int64>(ft,MEDCouplingTimeDiscretizationInt64::New(td),false)
{
}

MEDCouplingFieldInt64 *MEDCouplingFieldInt64::deepCopy() const
{
  return cloneWithMesh(true);
}

MEDCouplingFieldInt64 *MEDCouplingFieldInt64::clone(bool recDeepCpy) const
{
  return new MEDCouplingFieldInt64(*this,recDeepCpy);
}

MEDCouplingFieldDouble *MEDCouplingFieldInt64::convertToDblField() const
{
  MCAuto<MEDCouplingFieldTemplate> tmp(MEDCouplingFieldTemplate::New(*this));
  int t1,t2;
  double t0(getTime(t1,t2));
  MCAuto<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(*tmp,getTimeDiscretization()));
  ret->setTime(t0,t1,t2);
  if(getArray())
    {
      MCAuto<DataArrayDouble> arr(getArray()->convertToDblArr());
      ret->setArray(arr);
    }
  return ret.retn();
}
