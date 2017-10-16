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

%pythoncode %{
InterpKernelException.__reduce__=INTERPKERNELExceptionReduce
DataArrayDouble.__new__=classmethod(MEDCouplingDataArrayDoublenew)
DataArrayDouble.__iadd__=MEDCouplingDataArrayDoubleIadd
DataArrayDouble.__isub__=MEDCouplingDataArrayDoubleIsub
DataArrayDouble.__imul__=MEDCouplingDataArrayDoubleImul
DataArrayDouble.__idiv__=MEDCouplingDataArrayDoubleIdiv
DataArrayDouble.__ipow__=MEDCouplingDataArrayDoubleIpow
DataArrayDouble.__truediv__=DataArrayDouble.__div__
DataArrayDouble.__itruediv__=MEDCouplingDataArrayDoubleIdiv
DataArrayDouble.__rtruediv__=DataArrayDouble.__rdiv__
DataArrayDouble.__floordiv__=DataArrayDouble.__div__
DataArrayDouble.__ifloordiv__=MEDCouplingDataArrayDoubleIdiv
DataArrayDouble.__rfloordiv__=DataArrayDouble.__rdiv__
  
DataArrayInt.__reduce__=MEDCouplingDataArrayIntReduce
DataArrayInt.__iadd__=MEDCouplingDataArrayIntIadd
DataArrayInt.__isub__=MEDCouplingDataArrayIntIsub
DataArrayInt.__imul__=MEDCouplingDataArrayIntImul
DataArrayInt.__idiv__=MEDCouplingDataArrayIntIdiv
DataArrayInt.__imod__=MEDCouplingDataArrayIntImod
DataArrayInt.__ipow__=MEDCouplingDataArrayIntIpow
DataArrayInt.__truediv__=DataArrayInt.__div__
DataArrayInt.__itruediv__=MEDCouplingDataArrayIntIdiv
DataArrayInt.__rtruediv__=DataArrayInt.__rdiv__
DataArrayInt.__floordiv__=DataArrayInt.__div__
DataArrayInt.__ifloordiv__=MEDCouplingDataArrayIntIdiv
DataArrayInt.__rfloordiv__=DataArrayInt.__rdiv__

DataArrayByte.__new__=classmethod(MEDCouplingDataArrayBytenew)

DataArrayFloat.__new__=classmethod(MEDCouplingDataArrayFloatnew)
DataArrayFloat.__iadd__=MEDCouplingDataArrayFloatIadd
DataArrayFloat.__isub__=MEDCouplingDataArrayFloatIsub
DataArrayFloat.__imul__=MEDCouplingDataArrayFloatImul
DataArrayFloat.__idiv__=MEDCouplingDataArrayFloatIdiv
DataArrayFloat.__itruediv__=MEDCouplingDataArrayFloatIdiv
DataArrayFloat.__ifloordiv__=MEDCouplingDataArrayFloatIdiv

MEDCouplingFieldDouble.__new__=classmethod(MEDCouplingFieldDoublenew)
MEDCouplingFieldDouble.__iadd__=MEDCouplingFieldDoubleIadd
MEDCouplingFieldDouble.__isub__=MEDCouplingFieldDoubleIsub
MEDCouplingFieldDouble.__imul__=MEDCouplingFieldDoubleImul
MEDCouplingFieldDouble.__idiv__=MEDCouplingFieldDoubleIdiv
MEDCouplingFieldDouble.__ipow__=MEDCouplingFieldDoubleIpow
MEDCouplingFieldDouble.__truediv__=MEDCouplingFieldDouble.__div__
MEDCouplingFieldDouble.__rtruediv__=MEDCouplingFieldDouble.__rdiv__
MEDCouplingFieldDouble.__itruediv__=MEDCouplingFieldDoubleIdiv
MEDCouplingFieldDouble.__floordiv__=MEDCouplingFieldDouble.__div__
MEDCouplingFieldDouble.__rfloordiv__=MEDCouplingFieldDouble.__rdiv__
MEDCouplingFieldDouble.__ifloordiv__=MEDCouplingFieldDoubleIdiv

MEDCouplingFieldInt.__new__=classmethod(MEDCouplingFieldIntnew)
  
MEDCouplingFieldFloat.__new__=classmethod(MEDCouplingFieldFloatnew)

DataArrayDoubleTuple.__iadd__=MEDCouplingDataArrayDoubleTupleIadd
DataArrayDoubleTuple.__isub__=MEDCouplingDataArrayDoubleTupleIsub
DataArrayDoubleTuple.__imul__=MEDCouplingDataArrayDoubleTupleImul
DataArrayDoubleTuple.__idiv__=MEDCouplingDataArrayDoubleTupleIdiv
DataArrayDoubleTuple.__itruediv__=MEDCouplingDataArrayDoubleTupleIdiv
DataArrayDoubleTuple.__ifloordiv__=MEDCouplingDataArrayDoubleTupleIdiv

DataArrayIntTuple.__iadd__=MEDCouplingDataArrayIntTupleIadd
DataArrayIntTuple.__isub__=MEDCouplingDataArrayIntTupleIsub
DataArrayIntTuple.__imul__=MEDCouplingDataArrayIntTupleImul
DataArrayIntTuple.__idiv__=MEDCouplingDataArrayIntTupleIdiv
DataArrayIntTuple.__itruediv__=MEDCouplingDataArrayIntTupleIdiv
DataArrayIntTuple.__ifloordiv__=MEDCouplingDataArrayIntTupleIdiv
DataArrayIntTuple.__imod__=MEDCouplingDataArrayIntTupleImod

DenseMatrix.__iadd__=ParaMEDMEMDenseMatrixIadd
DenseMatrix.__isub__=ParaMEDMEMDenseMatrixIsub

MEDCouplingUMesh.__new__=classmethod(MEDCouplingUMeshnew)
MEDCoupling1DGTUMesh.__new__=classmethod(MEDCoupling1DGTUMeshnew)
MEDCoupling1SGTUMesh.__new__=classmethod(MEDCoupling1SGTUMeshnew)
MEDCouplingCurveLinearMesh.__new__=classmethod(MEDCouplingCurveLinearMeshnew)
MEDCouplingCMesh.__new__=classmethod(MEDCouplingCMeshnew)
MEDCouplingIMesh.__new__=classmethod(MEDCouplingIMeshnew)
MEDCouplingMappedExtrudedMesh.__new__=classmethod(MEDCouplingExtrudedMeshnew)

DataArrayAsciiCharIterator.__next__ = DataArrayAsciiCharIterator.next
DataArrayIntIterator.__next__ = DataArrayIntIterator.next
DataArrayDoubleIterator.__next__ = DataArrayDoubleIterator.next
MEDCouplingUMeshCellIterator.__next__ = MEDCouplingUMeshCellIterator.next
MEDCouplingUMeshCellByTypeIterator.__next__ = MEDCouplingUMeshCellByTypeIterator.next

del INTERPKERNELExceptionReduce
del MEDCouplingDataArrayDoublenew
del MEDCouplingDataArrayDoubleIadd
del MEDCouplingDataArrayDoubleIsub
del MEDCouplingDataArrayDoubleImul
del MEDCouplingDataArrayDoubleIdiv
del MEDCouplingFieldDoubleIadd
del MEDCouplingFieldDoubleIsub
del MEDCouplingFieldDoubleImul
del MEDCouplingFieldDoubleIdiv
del MEDCouplingFieldDoubleIpow
del MEDCouplingDataArrayIntnew
del MEDCouplingDataArrayIntIadd
del MEDCouplingDataArrayIntIsub
del MEDCouplingDataArrayIntImul
del MEDCouplingDataArrayIntIdiv
del MEDCouplingDataArrayIntImod
del MEDCouplingDataArrayBytenew
del MEDCouplingDataArrayFloatnew
del MEDCouplingDataArrayFloatIadd
del MEDCouplingDataArrayFloatIsub
del MEDCouplingDataArrayFloatImul
del MEDCouplingDataArrayFloatIdiv
del MEDCouplingDataArrayDoubleTupleIadd
del MEDCouplingDataArrayDoubleTupleIsub
del MEDCouplingDataArrayDoubleTupleImul
del MEDCouplingDataArrayDoubleTupleIdiv
del MEDCouplingDataArrayIntTupleIadd
del MEDCouplingDataArrayIntTupleIsub
del MEDCouplingDataArrayIntTupleImul
del MEDCouplingDataArrayIntTupleIdiv
del MEDCouplingDataArrayIntTupleImod
del ParaMEDMEMDenseMatrixIadd
del ParaMEDMEMDenseMatrixIsub
del MEDCouplingUMeshnew
del MEDCoupling1DGTUMeshnew
del MEDCoupling1SGTUMeshnew
del MEDCouplingCurveLinearMeshnew
del MEDCouplingCMeshnew
del MEDCouplingIMeshnew
del MEDCouplingExtrudedMeshnew
del MEDCouplingFieldDoublenew
%}
