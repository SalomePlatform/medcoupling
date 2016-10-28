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

DataArrayInt.__new__=classmethod(MEDCouplingDataArrayIntnew)
DataArrayInt.__iadd__=MEDCouplingDataArrayIntIadd
DataArrayInt.__isub__=MEDCouplingDataArrayIntIsub
DataArrayInt.__imul__=MEDCouplingDataArrayIntImul
DataArrayInt.__idiv__=MEDCouplingDataArrayIntIdiv
DataArrayInt.__imod__=MEDCouplingDataArrayIntImod
DataArrayInt.__ipow__=MEDCouplingDataArrayIntIpow

DataArrayByte.__new__=classmethod(MEDCouplingDataArrayBytenew)

MEDCouplingFieldDouble.__iadd__=MEDCouplingFieldDoubleIadd
MEDCouplingFieldDouble.__isub__=MEDCouplingFieldDoubleIsub
MEDCouplingFieldDouble.__imul__=MEDCouplingFieldDoubleImul
MEDCouplingFieldDouble.__idiv__=MEDCouplingFieldDoubleIdiv
MEDCouplingFieldDouble.__ipow__=MEDCouplingFieldDoubleIpow

DataArrayDoubleTuple.__iadd__=MEDCouplingDataArrayDoubleTupleIadd
DataArrayDoubleTuple.__isub__=MEDCouplingDataArrayDoubleTupleIsub
DataArrayDoubleTuple.__imul__=MEDCouplingDataArrayDoubleTupleImul
DataArrayDoubleTuple.__idiv__=MEDCouplingDataArrayDoubleTupleIdiv

DataArrayIntTuple.__iadd__=MEDCouplingDataArrayIntTupleIadd
DataArrayIntTuple.__isub__=MEDCouplingDataArrayIntTupleIsub
DataArrayIntTuple.__imul__=MEDCouplingDataArrayIntTupleImul
DataArrayIntTuple.__idiv__=MEDCouplingDataArrayIntTupleIdiv
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
MEDCouplingFieldDouble.__new__=classmethod(MEDCouplingFieldDoublenew)

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
