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

#include "MEDCouplingFieldInt.hxx"
#include "MEDCouplingFieldT.txx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingMesh.hxx"

using namespace MEDCoupling;

template class MEDCoupling::MEDCouplingFieldT<int>;

MEDCouplingFieldInt *MEDCouplingFieldInt::New(TypeOfField type, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldInt(type,td);
}

MEDCouplingFieldInt *MEDCouplingFieldInt::New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldInt(ft,td);
}

MEDCouplingFieldInt::MEDCouplingFieldInt(TypeOfField type, TypeOfTimeDiscretization td):MEDCouplingFieldT<int>(type,MEDCouplingTimeDiscretizationInt::New(td))
{
}

MEDCouplingFieldInt::MEDCouplingFieldInt(const MEDCouplingFieldInt& other, bool deepCpy):MEDCouplingFieldT<int>(other,deepCpy)
{
}

MEDCouplingFieldInt::MEDCouplingFieldInt(NatureOfField n, MEDCouplingTimeDiscretizationInt *td, MEDCouplingFieldDiscretization *type):MEDCouplingFieldT<int>(type,n,td)
{
}

/*!
 * ** WARINING : This method do not deeply copy neither mesh nor spatial discretization. Only a shallow copy (reference) is done for mesh and spatial discretization ! **
 */
MEDCouplingFieldInt::MEDCouplingFieldInt(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td):MEDCouplingFieldT<int>(ft,MEDCouplingTimeDiscretizationInt::New(td),false)
{
}

MEDCouplingFieldInt *MEDCouplingFieldInt::deepCopy() const
{
  return cloneWithMesh(true);
}

MEDCouplingFieldInt *MEDCouplingFieldInt::clone(bool recDeepCpy) const
{
  return new MEDCouplingFieldInt(*this,recDeepCpy);
}

MEDCouplingFieldDouble *MEDCouplingFieldInt::convertToDblField() const
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
