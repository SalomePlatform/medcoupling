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

%include "MEDCouplingCommon.i"

%pythoncode %{
def MEDCouplingDataArrayDoublenew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDouble____new___(cls,args)
def MEDCouplingDataArrayDoubleIadd(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDouble____iadd___(self, self, *args)
def MEDCouplingDataArrayDoubleIsub(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDouble____isub___(self, self, *args)
def MEDCouplingDataArrayDoubleImul(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDouble____imul___(self, self, *args)
def MEDCouplingDataArrayDoubleIdiv(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDouble____idiv___(self, self, *args)
def MEDCouplingDataArrayDoubleIpow(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDouble____ipow___(self, self, *args)
def MEDCouplingFieldDoublenew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldDouble____new___(cls,args)
def MEDCouplingFieldDoubleIadd(self,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldDouble____iadd___(self, self, *args)
def MEDCouplingFieldDoubleIsub(self,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldDouble____isub___(self, self, *args)
def MEDCouplingFieldDoubleImul(self,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldDouble____imul___(self, self, *args)
def MEDCouplingFieldDoubleIdiv(self,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldDouble____idiv___(self, self, *args)
def MEDCouplingFieldDoubleIpow(self,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldDouble____ipow___(self, self, *args)
def MEDCouplingFieldIntnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldInt____new___(cls,args)
def MEDCouplingFieldFloatnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldFloat____new___(cls,args)
def MEDCouplingDataArrayIntnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____new___(cls,args)
def MEDCouplingDataArrayIntIadd(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____iadd___(self, self, *args)
def MEDCouplingDataArrayIntIsub(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____isub___(self, self, *args)
def MEDCouplingDataArrayIntImul(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____imul___(self, self, *args)
def MEDCouplingDataArrayIntIdiv(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____idiv___(self, self, *args)
def MEDCouplingDataArrayIntImod(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____imod___(self, self, *args)
def MEDCouplingDataArrayIntIpow(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayInt____ipow___(self, self, *args)
def MEDCouplingDataArrayBytenew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayByte____new___(cls,args)
def MEDCouplingDataArrayFloatnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayFloat____new___(cls,args)
def MEDCouplingDataArrayFloatIadd(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayFloat____iadd___(self, self, *args)
def MEDCouplingDataArrayFloatIsub(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayFloat____isub___(self, self, *args)
def MEDCouplingDataArrayFloatImul(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayFloat____imul___(self, self, *args)
def MEDCouplingDataArrayFloatIdiv(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayFloat____idiv___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIadd(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDoubleTuple____iadd___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIsub(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDoubleTuple____isub___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleImul(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDoubleTuple____imul___(self, self, *args)
def MEDCouplingDataArrayDoubleTupleIdiv(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayDoubleTuple____idiv___(self, self, *args)
def MEDCouplingDataArrayIntTupleIadd(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayIntTuple____iadd___(self, self, *args)
def MEDCouplingDataArrayIntTupleIsub(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayIntTuple____isub___(self, self, *args)
def MEDCouplingDataArrayIntTupleImul(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayIntTuple____imul___(self, self, *args)
def MEDCouplingDataArrayIntTupleIdiv(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayIntTuple____idiv___(self, self, *args)
def MEDCouplingDataArrayIntTupleImod(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DataArrayIntTuple____imod___(self, self, *args)
def ParaMEDMEMDenseMatrixIadd(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DenseMatrix____iadd___(self, self, *args)
def ParaMEDMEMDenseMatrixIsub(self,*args):
    import _MEDCoupling
    return _MEDCoupling.DenseMatrix____isub___(self, self, *args)
def MEDCouplingUMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingUMesh____new___(cls,args)
def MEDCoupling1DGTUMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCoupling1DGTUMesh____new___(cls,args)
def MEDCoupling1SGTUMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCoupling1SGTUMesh____new___(cls,args)
def MEDCouplingCurveLinearMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingCurveLinearMesh____new___(cls,args)
def MEDCouplingCMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingCMesh____new___(cls,args)
def MEDCouplingIMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingIMesh____new___(cls,args)
def MEDCouplingExtrudedMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingMappedExtrudedMesh____new___(cls,args)
%}

%include "MEDCouplingFinalize.i"
