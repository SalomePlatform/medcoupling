// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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
def ParaMEDMEMDataArrayDoubleIadd(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDouble____iadd___(self, self, *args)
def ParaMEDMEMDataArrayDoubleIsub(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDouble____isub___(self, self, *args)
def ParaMEDMEMDataArrayDoubleImul(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDouble____imul___(self, self, *args)
def ParaMEDMEMDataArrayDoubleIdiv(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDouble____idiv___(self, self, *args)
def ParaMEDMEMDataArrayDoubleIpow(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDouble____ipow___(self, self, *args)
def ParaMEDMEMDataArrayIntIadd(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____iadd___(self, self, *args)
def ParaMEDMEMDataArrayIntIsub(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____isub___(self, self, *args)
def ParaMEDMEMDataArrayIntImul(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____imul___(self, self, *args)
def ParaMEDMEMDataArrayIntIdiv(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____idiv___(self, self, *args)
def ParaMEDMEMDataArrayIntImod(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____imod___(self, self, *args)
def ParaMEDMEMDataArrayIntIpow(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____ipow___(self, self, *args)
def ParaMEDMEMDataArrayDoubleTupleIadd(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDoubleTuple____iadd___(self, self, *args)
def ParaMEDMEMDataArrayDoubleTupleIsub(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDoubleTuple____isub___(self, self, *args)
def ParaMEDMEMDataArrayDoubleTupleImul(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDoubleTuple____imul___(self, self, *args)
def ParaMEDMEMDataArrayDoubleTupleIdiv(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDoubleTuple____idiv___(self, self, *args)
def ParaMEDMEMDataArrayIntTupleIadd(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayIntTuple____iadd___(self, self, *args)
def ParaMEDMEMDataArrayIntTupleIsub(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayIntTuple____isub___(self, self, *args)
def ParaMEDMEMDataArrayIntTupleImul(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayIntTuple____imul___(self, self, *args)
def ParaMEDMEMDataArrayIntTupleIdiv(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayIntTuple____idiv___(self, self, *args)
def ParaMEDMEMDataArrayIntTupleImod(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayIntTuple____imod___(self, self, *args)
%}

%pythoncode %{
DataArrayDouble.__iadd__=ParaMEDMEMDataArrayDoubleIadd
DataArrayDouble.__isub__=ParaMEDMEMDataArrayDoubleIsub
DataArrayDouble.__imul__=ParaMEDMEMDataArrayDoubleImul
DataArrayDouble.__idiv__=ParaMEDMEMDataArrayDoubleIdiv
DataArrayDouble.__ipow__=ParaMEDMEMDataArrayDoubleIpow

DataArrayInt.__iadd__=ParaMEDMEMDataArrayIntIadd
DataArrayInt.__isub__=ParaMEDMEMDataArrayIntIsub
DataArrayInt.__imul__=ParaMEDMEMDataArrayIntImul
DataArrayInt.__idiv__=ParaMEDMEMDataArrayIntIdiv
DataArrayInt.__imod__=ParaMEDMEMDataArrayIntImod
DataArrayInt.__ipow__=ParaMEDMEMDataArrayIntIpow

DataArrayDoubleTuple.__iadd__=ParaMEDMEMDataArrayDoubleTupleIadd
DataArrayDoubleTuple.__isub__=ParaMEDMEMDataArrayDoubleTupleIsub
DataArrayDoubleTuple.__imul__=ParaMEDMEMDataArrayDoubleTupleImul
DataArrayDoubleTuple.__idiv__=ParaMEDMEMDataArrayDoubleTupleIdiv

DataArrayIntTuple.__iadd__=ParaMEDMEMDataArrayIntTupleIadd
DataArrayIntTuple.__isub__=ParaMEDMEMDataArrayIntTupleIsub
DataArrayIntTuple.__imul__=ParaMEDMEMDataArrayIntTupleImul
DataArrayIntTuple.__idiv__=ParaMEDMEMDataArrayIntTupleIdiv
DataArrayIntTuple.__imod__=ParaMEDMEMDataArrayIntTupleImod

del ParaMEDMEMDataArrayDoubleIadd
del ParaMEDMEMDataArrayDoubleIsub
del ParaMEDMEMDataArrayDoubleImul
del ParaMEDMEMDataArrayDoubleIdiv
del ParaMEDMEMDataArrayIntIadd
del ParaMEDMEMDataArrayIntIsub
del ParaMEDMEMDataArrayIntImul
del ParaMEDMEMDataArrayIntIdiv
del ParaMEDMEMDataArrayIntImod
del ParaMEDMEMDataArrayDoubleTupleIadd
del ParaMEDMEMDataArrayDoubleTupleIsub
del ParaMEDMEMDataArrayDoubleTupleImul
del ParaMEDMEMDataArrayDoubleTupleIdiv
del ParaMEDMEMDataArrayIntTupleIadd
del ParaMEDMEMDataArrayIntTupleIsub
del ParaMEDMEMDataArrayIntTupleImul
del ParaMEDMEMDataArrayIntTupleIdiv
del ParaMEDMEMDataArrayIntTupleImod
%}
