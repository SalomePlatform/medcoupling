// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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

MEDCouplingFieldDouble.__iadd__=ParaMEDMEMMEDCouplingFieldDoubleIadd
MEDCouplingFieldDouble.__isub__=ParaMEDMEMMEDCouplingFieldDoubleIsub
MEDCouplingFieldDouble.__imul__=ParaMEDMEMMEDCouplingFieldDoubleImul
MEDCouplingFieldDouble.__idiv__=ParaMEDMEMMEDCouplingFieldDoubleIdiv
MEDCouplingFieldDouble.__ipow__=ParaMEDMEMMEDCouplingFieldDoubleIpow

DataArrayDoubleTuple.__iadd__=ParaMEDMEMDataArrayDoubleTupleIadd
DataArrayDoubleTuple.__isub__=ParaMEDMEMDataArrayDoubleTupleIsub
DataArrayDoubleTuple.__imul__=ParaMEDMEMDataArrayDoubleTupleImul
DataArrayDoubleTuple.__idiv__=ParaMEDMEMDataArrayDoubleTupleIdiv

DataArrayIntTuple.__iadd__=ParaMEDMEMDataArrayIntTupleIadd
DataArrayIntTuple.__isub__=ParaMEDMEMDataArrayIntTupleIsub
DataArrayIntTuple.__imul__=ParaMEDMEMDataArrayIntTupleImul
DataArrayIntTuple.__idiv__=ParaMEDMEMDataArrayIntTupleIdiv
DataArrayIntTuple.__imod__=ParaMEDMEMDataArrayIntTupleImod

DenseMatrix.__iadd__=ParaMEDMEMDenseMatrixIadd
DenseMatrix.__isub__=ParaMEDMEMDenseMatrixIsub

del ParaMEDMEMDataArrayDoubleIadd
del ParaMEDMEMDataArrayDoubleIsub
del ParaMEDMEMDataArrayDoubleImul
del ParaMEDMEMDataArrayDoubleIdiv
del ParaMEDMEMMEDCouplingFieldDoubleIadd
del ParaMEDMEMMEDCouplingFieldDoubleIsub
del ParaMEDMEMMEDCouplingFieldDoubleImul
del ParaMEDMEMMEDCouplingFieldDoubleIdiv
del ParaMEDMEMMEDCouplingFieldDoubleIpow
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
del ParaMEDMEMDenseMatrixIadd
del ParaMEDMEMDenseMatrixIsub
%}
