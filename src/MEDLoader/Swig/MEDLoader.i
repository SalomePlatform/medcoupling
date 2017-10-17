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
%}

%include "MEDCouplingFinalize.i"

%pythoncode %{
MEDFileMeshesIterator.__next__ = MEDFileMeshesIterator.next
MEDFileAnyTypeFieldMultiTSIterator.__next__ = MEDFileAnyTypeFieldMultiTSIterator.next
MEDFileFieldsIterator.__next__ = MEDFileFieldsIterator.next
%}

%pythoncode %{
def MEDCouplingMEDFileUMeshReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileUMesh,((),(self.__getstate__()),))
MEDFileUMesh.__reduce__=MEDCouplingMEDFileUMeshReduce
del MEDCouplingMEDFileUMeshReduce
def MEDCouplingMEDFileCMeshReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileCMesh,((self.serialize(),),(self.__getstate__()),))
MEDFileCMesh.__reduce__=MEDCouplingMEDFileCMeshReduce
del MEDCouplingMEDFileCMeshReduce
def MEDCouplingMEDFileCurveLinearMeshReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileCurveLinearMesh,((self.serialize(),),(self.__getstate__()),))
MEDFileCurveLinearMesh.__reduce__=MEDCouplingMEDFileCurveLinearMeshReduce
del MEDCouplingMEDFileCurveLinearMeshReduce
def MEDCouplingMEDFileDataReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileData,((self.serialize(),),(self.__getstate__()),))
MEDFileData.__reduce__=MEDCouplingMEDFileDataReduce
del MEDCouplingMEDFileDataReduce
def MEDCouplingMEDFileMeshesReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileMeshes,((self.serialize(),),(self.__getstate__()),))
MEDFileMeshes.__reduce__=MEDCouplingMEDFileMeshesReduce
del MEDCouplingMEDFileMeshesReduce
def MEDCouplingMEDFileFieldsReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileFields,((self.serialize(),),(self.__getstate__()),))
MEDFileFields.__reduce__=MEDCouplingMEDFileFieldsReduce
del MEDCouplingMEDFileFieldsReduce
def MEDCouplingMEDFileField1TSReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileField1TS,((self.serialize(),),(self.__getstate__()),))
MEDFileField1TS.__reduce__=MEDCouplingMEDFileField1TSReduce
del MEDCouplingMEDFileField1TSReduce
def MEDCouplingMEDFileFieldMultiTSReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileFieldMultiTS,((self.serialize(),),(self.__getstate__()),))
MEDFileFieldMultiTS.__reduce__=MEDCouplingMEDFileFieldMultiTSReduce
del MEDCouplingMEDFileFieldMultiTSReduce
def MEDCouplingMEDFileIntField1TSReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileIntField1TS,((self.serialize(),),(self.__getstate__()),))
MEDFileIntField1TS.__reduce__=MEDCouplingMEDFileIntField1TSReduce
def MEDCouplingMEDFileIntFieldMultiTSReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileIntFieldMultiTS,((self.serialize(),),(self.__getstate__()),))
MEDFileIntFieldMultiTS.__reduce__=MEDCouplingMEDFileIntFieldMultiTSReduce
del MEDCouplingMEDFileIntFieldMultiTSReduce
def MEDCouplingMEDFileFloatField1TSReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileFloatField1TS,((self.serialize(),),(self.__getstate__()),))
MEDFileFloatField1TS.__reduce__=MEDCouplingMEDFileFloatField1TSReduce
def MEDCouplingMEDFileFloatFieldMultiTSReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileFloatFieldMultiTS,((self.serialize(),),(self.__getstate__()),))
MEDFileFloatFieldMultiTS.__reduce__=MEDCouplingMEDFileFloatFieldMultiTSReduce
del MEDCouplingMEDFileFloatFieldMultiTSReduce
def MEDCouplingMEDFileParametersReduce(self):
  return MEDCouplingStdReduceFunct,(MEDFileParameters,((self.serialize(),),(self.__getstate__()),))
MEDFileParameters.__reduce__=MEDCouplingMEDFileParametersReduce
del MEDCouplingMEDFileParametersReduce
%}
