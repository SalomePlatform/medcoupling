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
def ParaMEDMEMMEDCouplingFieldDoublenew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldDouble____new___(cls,args)
def ParaMEDMEMMEDCouplingFieldDoubleIadd(self,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldDouble____iadd___(self, self, *args)
def ParaMEDMEMMEDCouplingFieldDoubleIsub(self,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldDouble____isub___(self, self, *args)
def ParaMEDMEMMEDCouplingFieldDoubleImul(self,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldDouble____imul___(self, self, *args)
def ParaMEDMEMMEDCouplingFieldDoubleIdiv(self,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldDouble____idiv___(self, self, *args)
def ParaMEDMEMMEDCouplingFieldDoubleIpow(self,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingFieldDouble____ipow___(self, self, *args)
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
def ParaMEDMEMMEDCouplingUMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingUMesh____new___(cls,args)
def ParaMEDMEMMEDCoupling1DGTUMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCoupling1DGTUMesh____new___(cls,args)
def ParaMEDMEMMEDCoupling1SGTUMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCoupling1SGTUMesh____new___(cls,args)
def ParaMEDMEMMEDCouplingCurveLinearMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingCurveLinearMesh____new___(cls,args)
def ParaMEDMEMMEDCouplingCMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingCMesh____new___(cls,args)
def ParaMEDMEMMEDCouplingIMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingIMesh____new___(cls,args)
def ParaMEDMEMMEDCouplingExtrudedMeshnew(cls,*args):
    import _MEDCoupling
    return _MEDCoupling.MEDCouplingExtrudedMesh____new___(cls,args)
%}

%include "MEDCouplingFinalize.i"
