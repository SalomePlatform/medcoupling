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

%include "MEDRenumberCommon.i"

%pythoncode %{
def MEDCouplingDataArrayDoublenew(cls,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayDouble____new___(cls,args)
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
def MEDCouplingDataArrayIntnew(cls,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt____new___(cls,args)
def MEDCouplingDataArrayIntIadd(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt____iadd___(self, self, *args)
def MEDCouplingDataArrayIntIsub(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt____isub___(self, self, *args)
def MEDCouplingDataArrayIntImul(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt____imul___(self, self, *args)
def MEDCouplingDataArrayIntIdiv(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt____idiv___(self, self, *args)
def MEDCouplingDataArrayIntImod(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt____imod___(self, self, *args)
def MEDCouplingDataArrayIntIpow(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayInt____ipow___(self, self, *args)
def MEDCouplingDataArrayBytenew(cls,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayByte____new___(cls,args)
def MEDCouplingDataArrayFloatnew(cls,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayFloat____new___(cls,args)
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
def MEDCouplingDataArrayIntTupleIadd(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayIntTuple____iadd___(self, self, *args)
def MEDCouplingDataArrayIntTupleIsub(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayIntTuple____isub___(self, self, *args)
def MEDCouplingDataArrayIntTupleImul(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayIntTuple____imul___(self, self, *args)
def MEDCouplingDataArrayIntTupleIdiv(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayIntTuple____idiv___(self, self, *args)
def MEDCouplingDataArrayIntTupleImod(self,*args):
    import _MEDRenumber
    return _MEDRenumber.DataArrayIntTuple____imod___(self, self, *args)
%}


%pythoncode %{
InterpKernelException.__reduce__=INTERPKERNELExceptionReduce
DataArrayDouble.__new__=classmethod(MEDCouplingDataArrayDoublenew)
DataArrayDouble.__iadd__=MEDCouplingDataArrayDoubleIadd
DataArrayDouble.__isub__=MEDCouplingDataArrayDoubleIsub
DataArrayDouble.__imul__=MEDCouplingDataArrayDoubleImul
DataArrayDouble.__idiv__=MEDCouplingDataArrayDoubleIdiv
DataArrayDouble.__ipow__=MEDCouplingDataArrayDoubleIpow

DataArrayInt.__new__=classmethod(MEDCouplingDataArrayIntnew)
DataArrayInt.__iadd__=MEDCouplingDataArrayIntIadd
DataArrayInt.__isub__=MEDCouplingDataArrayIntIsub
DataArrayInt.__imul__=MEDCouplingDataArrayIntImul
DataArrayInt.__idiv__=MEDCouplingDataArrayIntIdiv
DataArrayInt.__imod__=MEDCouplingDataArrayIntImod
DataArrayInt.__ipow__=MEDCouplingDataArrayIntIpow

DataArrayDoubleTuple.__iadd__=MEDCouplingDataArrayDoubleTupleIadd
DataArrayDoubleTuple.__isub__=MEDCouplingDataArrayDoubleTupleIsub
DataArrayDoubleTuple.__imul__=MEDCouplingDataArrayDoubleTupleImul
DataArrayDoubleTuple.__idiv__=MEDCouplingDataArrayDoubleTupleIdiv

DataArrayIntTuple.__iadd__=MEDCouplingDataArrayIntTupleIadd
DataArrayIntTuple.__isub__=MEDCouplingDataArrayIntTupleIsub
DataArrayIntTuple.__imul__=MEDCouplingDataArrayIntTupleImul
DataArrayIntTuple.__idiv__=MEDCouplingDataArrayIntTupleIdiv
DataArrayIntTuple.__imod__=MEDCouplingDataArrayIntTupleImod

del INTERPKERNELExceptionReduce
del MEDCouplingDataArrayDoublenew
del MEDCouplingDataArrayDoubleIadd
del MEDCouplingDataArrayDoubleIsub
del MEDCouplingDataArrayDoubleImul
del MEDCouplingDataArrayDoubleIdiv
del MEDCouplingDataArrayIntnew
del MEDCouplingDataArrayIntIadd
del MEDCouplingDataArrayIntIsub
del MEDCouplingDataArrayIntImul
del MEDCouplingDataArrayIntIdiv
del MEDCouplingDataArrayIntImod
del MEDCouplingDataArrayDoubleTupleIadd
del MEDCouplingDataArrayDoubleTupleIsub
del MEDCouplingDataArrayDoubleTupleImul
del MEDCouplingDataArrayDoubleTupleIdiv
del MEDCouplingDataArrayIntTupleIadd
del MEDCouplingDataArrayIntTupleIsub
del MEDCouplingDataArrayIntTupleImul
del MEDCouplingDataArrayIntTupleIdiv
del MEDCouplingDataArrayIntTupleImod
%}
