// Copyright (C) 2007-2025  CEA, EDF
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

%module MEDRenumber

%include "MEDRenumberCommon.i"

%pythoncode %{
def MEDCouplingDataArrayDoubleIadd(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayDouble____iadd___(self, self, *args)
def MEDCouplingDataArrayDoubleIsub(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayDouble____isub___(self, self, *args)
def MEDCouplingDataArrayDoubleImul(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayDouble____imul___(self, self, *args)
def MEDCouplingDataArrayDoubleIdiv(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayDouble____idiv___(self, self, *args)
def MEDCouplingDataArrayDoubleIpow(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayDouble____ipow___(self, self, *args)
def MEDCouplingDataArrayInt32Iadd(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt32____iadd___(self, self, *args)
def MEDCouplingDataArrayInt32Isub(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt32____isub___(self, self, *args)
def MEDCouplingDataArrayInt32Imul(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt32____imul___(self, self, *args)
def MEDCouplingDataArrayInt32Idiv(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt32____idiv___(self, self, *args)
def MEDCouplingDataArrayInt32Imod(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt32____imod___(self, self, *args)
def MEDCouplingDataArrayInt32Ipow(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt32____ipow___(self, self, *args)
def MEDCouplingDataArrayInt64Iadd(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt64____iadd___(self, self, *args)
def MEDCouplingDataArrayInt64Isub(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt64____isub___(self, self, *args)
def MEDCouplingDataArrayInt64Imul(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt64____imul___(self, self, *args)
def MEDCouplingDataArrayInt64Idiv(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt64____idiv___(self, self, *args)
def MEDCouplingDataArrayInt64Imod(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt64____imod___(self, self, *args)
def MEDCouplingDataArrayInt64Ipow(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt64____ipow___(self, self, *args)
def MEDCouplingDataArrayFloatIadd(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayFloat____iadd___(self, self, *args)
def MEDCouplingDataArrayFloatIsub(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayFloat____isub___(self, self, *args)
def MEDCouplingDataArrayFloatImul(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayFloat____imul___(self, self, *args)
def MEDCouplingDataArrayFloatIdiv(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayFloat____idiv___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIadd(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayDoubleTuple____iadd___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIsub(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayDoubleTuple____isub___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleImul(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayDoubleTuple____imul___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIdiv(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayDoubleTuple____idiv___(self, self, *args)
def MEDCouplingDataArrayInt32TupleIadd(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt32Tuple____iadd___(self, self, *args)
def MEDCouplingDataArrayInt32TupleIsub(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt32Tuple____isub___(self, self, *args)
def MEDCouplingDataArrayInt32TupleImul(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt32Tuple____imul___(self, self, *args)
def MEDCouplingDataArrayInt32TupleIdiv(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt32Tuple____idiv___(self, self, *args)
def MEDCouplingDataArrayInt32TupleImod(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt32Tuple____imod___(self, self, *args)
def MEDCouplingDataArrayInt64TupleIadd(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt64Tuple____iadd___(self, self, *args)
def MEDCouplingDataArrayInt64TupleIsub(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt64Tuple____isub___(self, self, *args)
def MEDCouplingDataArrayInt64TupleImul(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt64Tuple____imul___(self, self, *args)
def MEDCouplingDataArrayInt64TupleIdiv(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt64Tuple____idiv___(self, self, *args)
def MEDCouplingDataArrayInt64TupleImod(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt64Tuple____imod___(self, self, *args)
%}


%pythoncode %{
InterpKernelException.__reduce__=INTERPKERNELExceptionReduce
DataArrayDouble.__reduce__=MEDCouplingDataArrayDoubleReduce
DataArrayDouble.__iadd__=MEDCouplingDataArrayDoubleIadd
DataArrayDouble.__isub__=MEDCouplingDataArrayDoubleIsub
DataArrayDouble.__imul__=MEDCouplingDataArrayDoubleImul
DataArrayDouble.__idiv__=MEDCouplingDataArrayDoubleIdiv
DataArrayDouble.__ipow__=MEDCouplingDataArrayDoubleIpow

DataArrayInt32.__reduce__=MEDCouplingDataArrayInt32Reduce
DataArrayInt32.__iadd__=MEDCouplingDataArrayInt32Iadd
DataArrayInt32.__isub__=MEDCouplingDataArrayInt32Isub
DataArrayInt32.__imul__=MEDCouplingDataArrayInt32Imul
DataArrayInt32.__idiv__=MEDCouplingDataArrayInt32Idiv
DataArrayInt32.__imod__=MEDCouplingDataArrayInt32Imod
DataArrayInt32.__ipow__=MEDCouplingDataArrayInt32Ipow
DataArrayInt64.__reduce__=MEDCouplingDataArrayInt64Reduce
DataArrayInt64.__iadd__=MEDCouplingDataArrayInt64Iadd
DataArrayInt64.__isub__=MEDCouplingDataArrayInt64Isub
DataArrayInt64.__imul__=MEDCouplingDataArrayInt64Imul
DataArrayInt64.__idiv__=MEDCouplingDataArrayInt64Idiv
DataArrayInt64.__imod__=MEDCouplingDataArrayInt64Imod
DataArrayInt64.__ipow__=MEDCouplingDataArrayInt64Ipow

DataArrayDoubleTuple.__iadd__=MEDCouplingDataArrayDoubleTupleIadd
DataArrayDoubleTuple.__isub__=MEDCouplingDataArrayDoubleTupleIsub
DataArrayDoubleTuple.__imul__=MEDCouplingDataArrayDoubleTupleImul
DataArrayDoubleTuple.__idiv__=MEDCouplingDataArrayDoubleTupleIdiv

DataArrayInt32Tuple.__iadd__=MEDCouplingDataArrayInt32TupleIadd
DataArrayInt32Tuple.__isub__=MEDCouplingDataArrayInt32TupleIsub
DataArrayInt32Tuple.__imul__=MEDCouplingDataArrayInt32TupleImul
DataArrayInt32Tuple.__idiv__=MEDCouplingDataArrayInt32TupleIdiv
DataArrayInt32Tuple.__itruediv__=MEDCouplingDataArrayInt32TupleIdiv
DataArrayInt32Tuple.__ifloordiv__=MEDCouplingDataArrayInt32TupleIdiv
DataArrayInt32Tuple.__imod__=MEDCouplingDataArrayInt32TupleImod

DataArrayInt64Tuple.__iadd__=MEDCouplingDataArrayInt64TupleIadd
DataArrayInt64Tuple.__isub__=MEDCouplingDataArrayInt64TupleIsub
DataArrayInt64Tuple.__imul__=MEDCouplingDataArrayInt64TupleImul
DataArrayInt64Tuple.__idiv__=MEDCouplingDataArrayInt64TupleIdiv
DataArrayInt64Tuple.__itruediv__=MEDCouplingDataArrayInt64TupleIdiv
DataArrayInt64Tuple.__ifloordiv__=MEDCouplingDataArrayInt64TupleIdiv
DataArrayInt64Tuple.__imod__=MEDCouplingDataArrayInt64TupleImod




del INTERPKERNELExceptionReduce
del MEDCouplingDataArrayDoubleIadd
del MEDCouplingDataArrayDoubleIdiv
del MEDCouplingDataArrayDoubleImul
del MEDCouplingDataArrayDoubleIpow
del MEDCouplingDataArrayDoubleIsub
del MEDCouplingDataArrayDoubleReduce
del MEDCouplingDataArrayDoubleTupleIadd
del MEDCouplingDataArrayDoubleTupleIdiv
del MEDCouplingDataArrayDoubleTupleImul
del MEDCouplingDataArrayDoubleTupleIsub
del MEDCouplingDataArrayInt32Iadd
del MEDCouplingDataArrayInt32Idiv
del MEDCouplingDataArrayInt32Imod
del MEDCouplingDataArrayInt32Imul
del MEDCouplingDataArrayInt32Ipow
del MEDCouplingDataArrayInt32Isub
del MEDCouplingDataArrayInt32Reduce
del MEDCouplingDataArrayInt32TupleIadd
del MEDCouplingDataArrayInt32TupleIdiv
del MEDCouplingDataArrayInt32TupleImod
del MEDCouplingDataArrayInt32TupleImul
del MEDCouplingDataArrayInt32TupleIsub
del MEDCouplingDataArrayInt64Iadd
del MEDCouplingDataArrayInt64Idiv
del MEDCouplingDataArrayInt64Imod
del MEDCouplingDataArrayInt64Imul
del MEDCouplingDataArrayInt64Ipow
del MEDCouplingDataArrayInt64Isub
del MEDCouplingDataArrayInt64Reduce
del MEDCouplingDataArrayInt64TupleIadd
del MEDCouplingDataArrayInt64TupleIdiv
del MEDCouplingDataArrayInt64TupleImod
del MEDCouplingDataArrayInt64TupleImul
del MEDCouplingDataArrayInt64TupleIsub

%}
