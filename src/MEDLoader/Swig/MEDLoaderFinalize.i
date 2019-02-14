// Copyright (C) 2017-2019  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

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
