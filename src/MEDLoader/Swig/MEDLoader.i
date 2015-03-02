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
// Author : Anthony Geay (CEA/DEN)

%include "MEDLoaderCommon.i"

%pythoncode %{
def ParaMEDMEMDataArrayDoublenew(cls,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDouble____new___(cls,args)
def ParaMEDMEMDataArrayDoubleIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDouble____iadd___(self, self, *args)
def ParaMEDMEMDataArrayDoubleIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDouble____isub___(self, self, *args)
def ParaMEDMEMDataArrayDoubleImul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDouble____imul___(self, self, *args)
def ParaMEDMEMDataArrayDoubleIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDouble____idiv___(self, self, *args)
def ParaMEDMEMDataArrayDoubleIpow(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDouble____ipow___(self, self, *args)
def ParaMEDMEMMEDCouplingFieldDoublenew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldDouble____new___(cls,args)
def ParaMEDMEMMEDCouplingFieldDoubleIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldDouble____iadd___(self, self, *args)
def ParaMEDMEMMEDCouplingFieldDoubleIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldDouble____isub___(self, self, *args)
def ParaMEDMEMMEDCouplingFieldDoubleImul(self,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldDouble____imul___(self, self, *args)
def ParaMEDMEMMEDCouplingFieldDoubleIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldDouble____idiv___(self, self, *args)
def ParaMEDMEMMEDCouplingFieldDoubleIpow(self,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldDouble____ipow___(self, self, *args)
def ParaMEDMEMDataArrayIntnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____new___(cls,args)
def ParaMEDMEMDataArrayIntIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____iadd___(self, self, *args)
def ParaMEDMEMDataArrayIntIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____isub___(self, self, *args)
def ParaMEDMEMDataArrayIntImul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____imul___(self, self, *args)
def ParaMEDMEMDataArrayIntIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____idiv___(self, self, *args)
def ParaMEDMEMDataArrayIntImod(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____imod___(self, self, *args)
def ParaMEDMEMDataArrayIntIpow(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____ipow___(self, self, *args)
def ParaMEDMEMDataArrayDoubleTupleIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDoubleTuple____iadd___(self, self, *args)
def ParaMEDMEMDataArrayDoubleTupleIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDoubleTuple____isub___(self, self, *args)
def ParaMEDMEMDataArrayDoubleTupleImul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDoubleTuple____imul___(self, self, *args)
def ParaMEDMEMDataArrayDoubleTupleIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDoubleTuple____idiv___(self, self, *args)
def ParaMEDMEMDataArrayIntTupleIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayIntTuple____iadd___(self, self, *args)
def ParaMEDMEMDataArrayIntTupleIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayIntTuple____isub___(self, self, *args)
def ParaMEDMEMDataArrayIntTupleImul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayIntTuple____imul___(self, self, *args)
def ParaMEDMEMDataArrayIntTupleIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayIntTuple____idiv___(self, self, *args)
def ParaMEDMEMDataArrayIntTupleImod(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayIntTuple____imod___(self, self, *args)
def ParaMEDMEMDenseMatrixIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DenseMatrix____iadd___(self, self, *args)
def ParaMEDMEMDenseMatrixIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DenseMatrix____isub___(self, self, *args)
def ParaMEDMEMMEDCouplingUMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingUMesh____new___(cls,args)
def ParaMEDMEMMEDCoupling1DGTUMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCoupling1DGTUMesh____new___(cls,args)
def ParaMEDMEMMEDCoupling1SGTUMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCoupling1SGTUMesh____new___(cls,args)
def ParaMEDMEMMEDCouplingCurveLinearMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingCurveLinearMesh____new___(cls,args)
def ParaMEDMEMMEDCouplingCMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingCMesh____new___(cls,args)
def ParaMEDMEMMEDCouplingIMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingIMesh____new___(cls,args)
def ParaMEDMEMMEDCouplingExtrudedMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingExtrudedMesh____new___(cls,args)
%}

%include "MEDCouplingFinalize.i"
