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
// Author : Anthony Geay (CEA/DEN)

%include "MEDLoaderCommon.i"

%pythoncode %{
def MEDCouplingDataArrayDoublenew(cls,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayDouble____new___(cls,args)
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
def MEDCouplingFieldDoublenew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldDouble____new___(cls,args)
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
def MEDCouplingFieldIntnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldInt____new___(cls,args)
def MEDCouplingFieldFloatnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingFieldFloat____new___(cls,args)
def MEDCouplingDataArrayBytenew(cls,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayByte____new___(cls,args)
def MEDCouplingDataArrayFloatnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayFloat____new___(cls,args)
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
def MEDCouplingDataArrayIntnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____new___(cls,args)
def MEDCouplingDataArrayIntIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____iadd___(self, self, *args)
def MEDCouplingDataArrayIntIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____isub___(self, self, *args)
def MEDCouplingDataArrayIntImul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____imul___(self, self, *args)
def MEDCouplingDataArrayIntIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____idiv___(self, self, *args)
def MEDCouplingDataArrayIntImod(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____imod___(self, self, *args)
def MEDCouplingDataArrayIntIpow(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayInt____ipow___(self, self, *args)
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
def MEDCouplingDataArrayIntTupleIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayIntTuple____iadd___(self, self, *args)
def MEDCouplingDataArrayIntTupleIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayIntTuple____isub___(self, self, *args)
def MEDCouplingDataArrayIntTupleImul(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayIntTuple____imul___(self, self, *args)
def MEDCouplingDataArrayIntTupleIdiv(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayIntTuple____idiv___(self, self, *args)
def MEDCouplingDataArrayIntTupleImod(self,*args):
    import _MEDLoader
    return _MEDLoader.DataArrayIntTuple____imod___(self, self, *args)
def ParaMEDMEMDenseMatrixIadd(self,*args):
    import _MEDLoader
    return _MEDLoader.DenseMatrix____iadd___(self, self, *args)
def ParaMEDMEMDenseMatrixIsub(self,*args):
    import _MEDLoader
    return _MEDLoader.DenseMatrix____isub___(self, self, *args)
def MEDCouplingUMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingUMesh____new___(cls,args)
def MEDCoupling1DGTUMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCoupling1DGTUMesh____new___(cls,args)
def MEDCoupling1SGTUMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCoupling1SGTUMesh____new___(cls,args)
def MEDCouplingCurveLinearMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingCurveLinearMesh____new___(cls,args)
def MEDCouplingCMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingCMesh____new___(cls,args)
def MEDCouplingIMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingIMesh____new___(cls,args)
def MEDCouplingExtrudedMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDCouplingMappedExtrudedMesh____new___(cls,args)
%}

%pythoncode %{
def MEDCouplingMEDFileUMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileUMesh____new___(cls,args)
def MEDCouplingMEDFileCMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileCMesh____new___(cls,args)
def MEDCouplingMEDFileCurveLinearMeshnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileCurveLinearMesh____new___(cls,args)
def MEDCouplingMEDFileMeshesnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileMeshes____new___(cls,args)
def MEDCouplingMEDFileDatanew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileData____new___(cls,args)
def MEDCouplingMEDFileFieldsnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileFields____new___(cls,args)
def MEDCouplingMEDFileField1TSnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileField1TS____new___(cls,args)
def MEDCouplingMEDFileFieldMultiTSnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileFieldMultiTS____new___(cls,args)
def MEDCouplingMEDFileIntField1TSnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileIntField1TS____new___(cls,args)
def MEDCouplingMEDFileIntFieldMultiTSnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileIntFieldMultiTS____new___(cls,args)
def MEDCouplingMEDFileFloatField1TSnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileFloatField1TS____new___(cls,args)
def MEDCouplingMEDFileFloatFieldMultiTSnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileFloatFieldMultiTS____new___(cls,args)
def MEDCouplingMEDFileParametersnew(cls,*args):
    import _MEDLoader
    return _MEDLoader.MEDFileParameters____new___(cls,args)
%}

%include "MEDCouplingFinalize.i"

%pythoncode %{
MEDFileMeshesIterator.__next__ = MEDFileMeshesIterator.next
MEDFileAnyTypeFieldMultiTSIterator.__next__ = MEDFileAnyTypeFieldMultiTSIterator.next
MEDFileFieldsIterator.__next__ = MEDFileFieldsIterator.next
%}

%pythoncode %{
MEDFileUMesh.__new__=classmethod(MEDCouplingMEDFileUMeshnew)
del MEDCouplingMEDFileUMeshnew
MEDFileCMesh.__new__=classmethod(MEDCouplingMEDFileCMeshnew)
del MEDCouplingMEDFileCMeshnew
MEDFileCurveLinearMesh.__new__=classmethod(MEDCouplingMEDFileCurveLinearMeshnew)
del MEDCouplingMEDFileCurveLinearMeshnew
MEDFileData.__new__=classmethod(MEDCouplingMEDFileDatanew)
del MEDCouplingMEDFileDatanew
MEDFileMeshes.__new__=classmethod(MEDCouplingMEDFileMeshesnew)
del MEDCouplingMEDFileMeshesnew
MEDFileFields.__new__=classmethod(MEDCouplingMEDFileFieldsnew)
del MEDCouplingMEDFileFieldsnew
MEDFileField1TS.__new__=classmethod(MEDCouplingMEDFileField1TSnew)
del MEDCouplingMEDFileField1TSnew
MEDFileFieldMultiTS.__new__=classmethod(MEDCouplingMEDFileFieldMultiTSnew)
del MEDCouplingMEDFileFieldMultiTSnew
MEDFileIntField1TS.__new__=classmethod(MEDCouplingMEDFileIntField1TSnew)
del MEDCouplingMEDFileIntField1TSnew
MEDFileIntFieldMultiTS.__new__=classmethod(MEDCouplingMEDFileIntFieldMultiTSnew)
del MEDCouplingMEDFileIntFieldMultiTSnew
MEDFileFloatField1TS.__new__=classmethod(MEDCouplingMEDFileFloatField1TSnew)
del MEDCouplingMEDFileFloatField1TSnew
MEDFileFloatFieldMultiTS.__new__=classmethod(MEDCouplingMEDFileFloatFieldMultiTSnew)
del MEDCouplingMEDFileFloatFieldMultiTSnew
MEDFileParameters.__new__=classmethod(MEDCouplingMEDFileParametersnew)
del MEDCouplingMEDFileParametersnew
%}
