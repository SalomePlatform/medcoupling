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
// Author : Anthony Geay (CEA/DEN)

%module MEDLoader

%include "MEDLoaderCommon.i"

%pythoncode %{
def MEDCouplingDataArrayDoubleIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDouble____iadd___(self, self, *args)
def MEDCouplingDataArrayDoubleIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDouble____isub___(self, self, *args)
def MEDCouplingDataArrayDoubleImul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDouble____imul___(self, self, *args)
def MEDCouplingDataArrayDoubleIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDouble____idiv___(self, self, *args)
def MEDCouplingDataArrayDoubleIpow(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDouble____ipow___(self, self, *args)
def MEDCouplingFieldDoubleIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldDouble____iadd___(self, self, *args)
def MEDCouplingFieldDoubleIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldDouble____isub___(self, self, *args)
def MEDCouplingFieldDoubleImul(self,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldDouble____imul___(self, self, *args)
def MEDCouplingFieldDoubleIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldDouble____idiv___(self, self, *args)
def MEDCouplingFieldDoubleIpow(self,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldDouble____ipow___(self, self, *args)
def MEDCouplingDataArrayFloatIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayFloat____iadd___(self, self, *args)
def MEDCouplingDataArrayFloatIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayFloat____isub___(self, self, *args)
def MEDCouplingDataArrayFloatImul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayFloat____imul___(self, self, *args)
def MEDCouplingDataArrayFloatIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayFloat____idiv___(self, self, *args)
def MEDCouplingDataArrayInt32Iadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt32____iadd___(self, self, *args)
def MEDCouplingDataArrayInt32Isub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt32____isub___(self, self, *args)
def MEDCouplingDataArrayInt32Imul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt32____imul___(self, self, *args)
def MEDCouplingDataArrayInt32Idiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt32____idiv___(self, self, *args)
def MEDCouplingDataArrayInt32Imod(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt32____imod___(self, self, *args)
def MEDCouplingDataArrayInt32Ipow(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt32____ipow___(self, self, *args)
def MEDCouplingDataArrayInt64Iadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt64____iadd___(self, self, *args)
def MEDCouplingDataArrayInt64Isub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt64____isub___(self, self, *args)
def MEDCouplingDataArrayInt64Imul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt64____imul___(self, self, *args)
def MEDCouplingDataArrayInt64Idiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt64____idiv___(self, self, *args)
def MEDCouplingDataArrayInt64Imod(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt64____imod___(self, self, *args)
def MEDCouplingDataArrayInt64Ipow(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt64____ipow___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDoubleTuple____iadd___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDoubleTuple____isub___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleImul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDoubleTuple____imul___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDoubleTuple____idiv___(self, self, *args)
def MEDCouplingDataArrayInt32TupleIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt32Tuple____iadd___(self, self, *args)
def MEDCouplingDataArrayInt32TupleIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt32Tuple____isub___(self, self, *args)
def MEDCouplingDataArrayInt32TupleImul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt32Tuple____imul___(self, self, *args)
def MEDCouplingDataArrayInt32TupleIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt32Tuple____idiv___(self, self, *args)
def MEDCouplingDataArrayInt32TupleImod(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt32Tuple____imod___(self, self, *args)
def MEDCouplingDataArrayInt64TupleIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt64Tuple____iadd___(self, self, *args)
def MEDCouplingDataArrayInt64TupleIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt64Tuple____isub___(self, self, *args)
def MEDCouplingDataArrayInt64TupleImul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt64Tuple____imul___(self, self, *args)
def MEDCouplingDataArrayInt64TupleIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt64Tuple____idiv___(self, self, *args)
def MEDCouplingDataArrayInt64TupleImod(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt64Tuple____imod___(self, self, *args)
def MEDCouplingDenseMatrixIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DenseMatrix____iadd___(self, self, *args)
def MEDCouplingDenseMatrixIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DenseMatrix____isub___(self, self, *args)
%}

%include "MEDCouplingFinalize.i"

%include "MEDLoaderFinalize.i"
