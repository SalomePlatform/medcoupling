#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2012  CEA/DEN, EDF R&D, OPEN CASCADE
#
# Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
# CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

###################################################################################
# This Python script uses the wrapped C++ class MESHING to buid a mesh from only
# primitive data like coordinates (Pythoin double array) and connectivity (Python
# integer arrays). It is the Python equivalent of the C++ program
# test_MEDMEM_Meshing.cxx in the ../MEDMEM directory of the SALOME distribution
###################################################################################
#
from libMEDMEM_Swig import *

# files name to save the generated MESH(ING) in different format
# Med V2.1 Med V2.2 and vtk

med21FileName = "toto21.med"

med22FileName = "toto22.med"

vtkFileName = "toto.vtk"

myMeshing = MESHING()

myMeshing.setName("meshing")

# definition of the coordinates

spaceDimension = 3

numberOfNodes = 19

coordinates = [
    0.0, 0.0, 0.0  ,
    0.0, 0.0, 1.0  ,
    2.0, 0.0, 1.0  ,
    0.0, 2.0, 1.0  ,
    -2.0, 0.0, 1.0 ,
    0.0, -2.0, 1.0 ,
    1.0, 1.0, 2.0  ,
    -1.0, 1.0, 2.0 ,
    -1.0, -1.0, 2.0,
    1.0, -1.0, 2.0 ,
    1.0, 1.0, 3.0  ,
    -1.0, 1.0, 3.0 ,
    -1.0, -1.0, 3.0,
    1.0, -1.0, 3.0 ,
    1.0, 1.0, 4.0  ,
    -1.0, 1.0, 4.0 ,
    -1.0, -1.0, 4.0,
    1.0, -1.0, 4.0 ,
    0.0, 0.0, 5.0]

myMeshing.setCoordinates(spaceDimension,numberOfNodes,coordinates,"CARTESIAN",MED_FULL_INTERLACE)

for i in range(spaceDimension):
    unit = "cm      "
    if (i == 0):
        name = "X       "
    elif (i == 1):
        name = "Y       "
    elif (i == 2):
        name = "Z       "

    myMeshing.setCoordinateName(name,i)
    myMeshing.setCoordinateUnit(unit,i)

# definition of connectivities
# cell part

numberOfTypes = 3
entity = MED_CELL

types = []
numberOfElements = []

types.append(MED_TETRA4)
numberOfElements.append(12)

types.append(MED_PYRA5)
numberOfElements.append(2)

types.append(MED_HEXA8)
numberOfElements.append(2)

myMeshing.setNumberOfTypes(numberOfTypes,entity)
myMeshing.setTypes(types,entity)
myMeshing.setNumberOfElements(numberOfElements,entity)

connectivityTetra = [
    1,2,3,6 ,
    1,2,4,3 ,
    1,2,5,4 ,
    1,2,6,5 ,
    2,7,4,3 ,
    2,8,5,4 ,
    2,9,6,5 ,
    2,10,3,6,
    2,7,3,10,
    2,8,4,7 ,
    2,9,5,8 ,
    2,10,6,9]

myMeshing.setConnectivity(entity,types[0],connectivityTetra)

connectivityPyra = [
    7,8,9,10,2,
    15,18,17,16,19]

myMeshing.setConnectivity(entity,types[1],connectivityPyra)

connectivityHexa = [
    11,12,13,14,7,8,9,10,
    15,16,17,18,11,12,13,14]

myMeshing.setConnectivity(entity,types[2],connectivityPyra)

# face part

numberOfTypes = 2
entity = MED_FACE

types = []
numberOfElements = []

types.append(MED_TRIA3)
numberOfElements.append(4)

types.append(MED_QUAD4)
numberOfElements.append(4)

myMeshing.setNumberOfTypes(numberOfTypes,entity)
myMeshing.setTypes(types,entity)
myMeshing.setNumberOfElements(numberOfElements,entity)

connectivityTria = [
    1,4,3,
    1,5,4,
    1,6,5,
    1,3,6]

myMeshing.setConnectivity(entity,types[0],connectivityPyra)

connectivityQuad = [
    7,8,9,10   ,
    11,12,13,14,
    11,7,8,12  ,
    12,8,9,13]

myMeshing.setConnectivity(entity,types[1],connectivityQuad)

# edge part

# adding GROUPs
# on Node

myGroup = GROUP()
myGroup.setName("SomeNodes")
myGroup.setMesh(myMeshing)
myGroup.setEntity(MED_NODE)
myGroup.setNumberOfGeometricType(1)

myTypes = [MED_NONE]
myGroup.setGeometricType(myTypes)

myNumberOfElements = [4]
myGroup.setNumberOfElements(myNumberOfElements)

index = [1,5]
values = [1,4,5,7]
myGroup.setNumber(index,values)

myMeshing.addGroup(myGroup)

myGroup = GROUP()
myGroup.setName("OtherNodes")
myGroup.setMesh(myMeshing)
myGroup.setEntity(MED_NODE)
myGroup.setNumberOfGeometricType(1)

myTypes = [MED_NONE]
myGroup.setGeometricType(myTypes)

myNumberOfElements = [3]
myGroup.setNumberOfElements(myNumberOfElements)

index = [1,4]
values = [2,3,6]
myGroup.setNumber(index,values)

myMeshing.addGroup(myGroup)

# on Cell

myGroup = GROUP()
myGroup.setName("SomeCells")
myGroup.setMesh(myMeshing)
myGroup.setEntity(MED_CELL)
myGroup.setNumberOfGeometricType(3)

myTypes = [MED_TETRA4,MED_PYRA5,MED_HEXA8]
myGroup.setGeometricType(myTypes)

myNumberOfElements = [4,1,2]
myGroup.setNumberOfElements(myNumberOfElements)

index = [1,5,6,8]
values = [
    2,7,8,12,
    13,
    15,16
    ]
myGroup.setNumber(index,values)

myMeshing.addGroup(myGroup)

myGroup = GROUP()
myGroup.setName("OtherCells")
myGroup.setMesh(myMeshing)
myGroup.setEntity(MED_CELL)
myGroup.setNumberOfGeometricType(2)

myTypes = [MED_TETRA4,MED_PYRA5]
myGroup.setGeometricType(myTypes)

myNumberOfElements = [4,1]
myGroup.setNumberOfElements(myNumberOfElements)

index = [1,5,6]
values = [
    3,4,5,9,
    14
    ]
myGroup.setNumber(index,values)

myMeshing.addGroup(myGroup)

# on Face

myGroup = GROUP()
myGroup.setName("SomeFaces")
myGroup.setMesh(myMeshing)
myGroup.setEntity(MED_FACE)
myGroup.setNumberOfGeometricType(2)

myTypes = [MED_TRIA3,MED_QUAD4]
myGroup.setGeometricType(myTypes)

myNumberOfElements = [2,3]
myGroup.setNumberOfElements(myNumberOfElements)

index = [1,3,6]
values = [
    2,4,
    5,6,8
    ]
myGroup.setNumber(index,values)

myMeshing.addGroup(myGroup)

myGroup = GROUP()
myGroup.setName("OtherFaces")
myGroup.setMesh(myMeshing)
myGroup.setEntity(MED_FACE)
myGroup.setNumberOfGeometricType(1)

myTypes = [MED_TRIA3]
myGroup.setGeometricType(myTypes)

myNumberOfElements = [2]
myGroup.setNumberOfElements(myNumberOfElements)

index = [1,3]
values = [
    1,3
    ]
myGroup.setNumber(index,values)

myMeshing.addGroup(myGroup)

# saving of the generated mesh in MED and VTK format

myMeshing.write(MED_DRIVER,med22FileName)

myMeshing.write(VTK_DRIVER,vtkFileName)

# we build now 8 fields : 4 fields double (integer) :
#                         2 fields on nodes (cells) :
#                         1 scalar (vector)

supportOnNodes = myMeshing.getSupportOnAll(MED_NODE)
numberOfNodes = supportOnNodes.getNumberOfElements(MED_ALL_ELEMENTS)

supportOnCells = myMeshing.getSupportOnAll(MED_CELL)
numberOfCells = supportOnCells.getNumberOfElements(MED_ALL_ELEMENTS)

fieldDoubleScalarOnNodes = FIELDDOUBLE(supportOnNodes,1)
fieldDoubleScalarOnNodes.setName("fieldScalarDoubleNode")
fieldDoubleScalarOnNodes.setIterationNumber(-1)
fieldDoubleScalarOnNodes.setOrderNumber(-1)
fieldDoubleScalarOnNodes.setTime(0.0)

fieldDoubleScalarOnNodes.setComponentName(1,"Vx")
fieldDoubleScalarOnNodes.setComponentDescription(1,"comp1")
fieldDoubleScalarOnNodes.setMEDComponentUnit(1,"unit1")

fieldDoubleVectorOnNodes = FIELDDOUBLE(supportOnNodes,spaceDimension)
fieldDoubleVectorOnNodes.setName("fieldVectorDoubleNode")
fieldDoubleVectorOnNodes.setIterationNumber(-1)
fieldDoubleVectorOnNodes.setOrderNumber(-1)
fieldDoubleVectorOnNodes.setTime(0.0)

fieldDoubleVectorOnNodes.setComponentName(1,"Vx")
fieldDoubleVectorOnNodes.setComponentDescription(1,"comp1")
fieldDoubleVectorOnNodes.setMEDComponentUnit(1,"unit1")
fieldDoubleVectorOnNodes.setComponentName(2,"Vy")
fieldDoubleVectorOnNodes.setComponentDescription(2,"comp2")
fieldDoubleVectorOnNodes.setMEDComponentUnit(2,"unit2")
fieldDoubleVectorOnNodes.setComponentName(3,"Vz")
fieldDoubleVectorOnNodes.setComponentDescription(3,"comp3")
fieldDoubleVectorOnNodes.setMEDComponentUnit(3,"unit3")

fieldDoubleScalarOnCells = FIELDDOUBLE(supportOnCells,1)
fieldDoubleScalarOnCells.setName("fieldScalarDoubleCell")
fieldDoubleScalarOnCells.setIterationNumber(-1)
fieldDoubleScalarOnCells.setOrderNumber(-1)
fieldDoubleScalarOnCells.setTime(0.0)

fieldDoubleScalarOnCells.setComponentName(1,"Vx")
fieldDoubleScalarOnCells.setComponentDescription(1,"comp1")
fieldDoubleScalarOnCells.setMEDComponentUnit(1,"unit1")

fieldDoubleVectorOnCells = FIELDDOUBLE(supportOnCells,spaceDimension)
fieldDoubleVectorOnCells.setName("fieldVectorrDoubleCell")
fieldDoubleVectorOnCells.setIterationNumber(-1)
fieldDoubleVectorOnCells.setOrderNumber(-1)
fieldDoubleVectorOnCells.setTime(0.0)

fieldDoubleVectorOnCells.setComponentName(1,"Vx")
fieldDoubleVectorOnCells.setComponentDescription(1,"comp1")
fieldDoubleVectorOnCells.setMEDComponentUnit(1,"unit1")
fieldDoubleVectorOnCells.setComponentName(2,"Vy")
fieldDoubleVectorOnCells.setComponentDescription(2,"comp2")
fieldDoubleVectorOnCells.setMEDComponentUnit(2,"unit2")
fieldDoubleVectorOnCells.setComponentName(3,"Vz")
fieldDoubleVectorOnCells.setComponentDescription(3,"comp3")
fieldDoubleVectorOnCells.setMEDComponentUnit(3,"unit3")

fieldIntScalarOnNodes = FIELDINT(supportOnNodes,1)
fieldIntScalarOnNodes.setName("fieldScalarIntNode")
fieldIntScalarOnNodes.setIterationNumber(-1)
fieldIntScalarOnNodes.setOrderNumber(-1)
fieldIntScalarOnNodes.setTime(0.0)

fieldIntScalarOnNodes.setComponentName(1,"Vx")
fieldIntScalarOnNodes.setComponentDescription(1,"comp1")
fieldIntScalarOnNodes.setMEDComponentUnit(1,"unit1")

fieldIntVectorOnNodes = FIELDINT(supportOnNodes,spaceDimension)
fieldIntVectorOnNodes.setName("fieldVectorIntNode")
fieldIntVectorOnNodes.setIterationNumber(-1)
fieldIntVectorOnNodes.setOrderNumber(-1)
fieldIntVectorOnNodes.setTime(0.0)

fieldIntVectorOnNodes.setComponentName(1,"Vx")
fieldIntVectorOnNodes.setComponentDescription(1,"comp1")
fieldIntVectorOnNodes.setMEDComponentUnit(1,"unit1")
fieldIntVectorOnNodes.setComponentName(2,"Vy")
fieldIntVectorOnNodes.setComponentDescription(2,"comp2")
fieldIntVectorOnNodes.setMEDComponentUnit(2,"unit2")
fieldIntVectorOnNodes.setComponentName(3,"Vz")
fieldIntVectorOnNodes.setComponentDescription(3,"comp3")
fieldIntVectorOnNodes.setMEDComponentUnit(3,"unit3")

fieldIntScalarOnCells = FIELDINT(supportOnCells,1)
fieldIntScalarOnCells.setName("fieldScalarIntCell")
fieldIntScalarOnCells.setIterationNumber(-1)
fieldIntScalarOnCells.setOrderNumber(-1)
fieldIntScalarOnCells.setTime(0.0)

fieldIntScalarOnCells.setComponentName(1,"Vx")
fieldIntScalarOnCells.setComponentDescription(1,"comp1")
fieldIntScalarOnCells.setMEDComponentUnit(1,"unit1")

fieldIntVectorOnCells = FIELDINT(supportOnCells,spaceDimension)
fieldIntVectorOnCells.setName("fieldVectorrIntCell")
fieldIntVectorOnCells.setIterationNumber(-1)
fieldIntVectorOnCells.setOrderNumber(-1)
fieldIntVectorOnCells.setTime(0.0)

fieldIntVectorOnCells.setComponentName(1,"Vx")
fieldIntVectorOnCells.setComponentDescription(1,"comp1")
fieldIntVectorOnCells.setMEDComponentUnit(1,"unit1")
fieldIntVectorOnCells.setComponentName(2,"Vy")
fieldIntVectorOnCells.setComponentDescription(2,"comp2")
fieldIntVectorOnCells.setMEDComponentUnit(2,"unit2")
fieldIntVectorOnCells.setComponentName(3,"Vz")
fieldIntVectorOnCells.setComponentDescription(3,"comp3")
fieldIntVectorOnCells.setMEDComponentUnit(3,"unit3")

for i in range(numberOfNodes):
    valueInt1 = i+1
    valueInt2 = i+2
    valueInt3 = i+3
    valueDbl1 = valueInt1*0.1
    valueDbl2 = valueInt2*0.1
    valueDbl3 = valueInt3*0.1
    fieldDoubleScalarOnNodes.setValueIJ(i+1,1,valueDbl1)

    fieldIntScalarOnNodes.setValueIJ(i+1,1,valueInt1)

    fieldDoubleVectorOnNodes.setValueIJ(i+1,1,valueDbl1)
    fieldDoubleVectorOnNodes.setValueIJ(i+1,2,valueDbl2)
    fieldDoubleVectorOnNodes.setValueIJ(i+1,3,valueDbl3)

    fieldIntVectorOnNodes.setValueIJ(i+1,1,valueInt1)
    fieldIntVectorOnNodes.setValueIJ(i+1,2,valueInt2)
    fieldIntVectorOnNodes.setValueIJ(i+1,3,valueInt3)

for i in range(numberOfCells):
    valueInt1 = i+1
    valueInt2 = i+2
    valueInt3 = i+3
    valueDbl1 = valueInt1*0.1
    valueDbl2 = valueInt2*0.1
    valueDbl3 = valueInt3*0.1
    fieldDoubleScalarOnCells.setValueIJ(i+1,1,valueDbl1)

    fieldIntScalarOnCells.setValueIJ(i+1,1,valueInt1)

    fieldDoubleVectorOnCells.setValueIJ(i+1,1,valueDbl1)
    fieldDoubleVectorOnCells.setValueIJ(i+1,2,valueDbl2)
    fieldDoubleVectorOnCells.setValueIJ(i+1,3,valueDbl3)

    fieldIntVectorOnCells.setValueIJ(i+1,1,valueInt1)
    fieldIntVectorOnCells.setValueIJ(i+1,2,valueInt2)
    fieldIntVectorOnCells.setValueIJ(i+1,3,valueInt3)

fieldIntScalarOnNodes.write(MED_DRIVER,med21FileName)
fieldDoubleVectorOnNodes.write(MED_DRIVER,med21FileName)
fieldIntVectorOnNodes.write(MED_DRIVER,med21FileName)
fieldDoubleScalarOnCells.write(MED_DRIVER,med21FileName)
fieldIntScalarOnCells.write(MED_DRIVER,med21FileName)
fieldDoubleVectorOnCells.write(MED_DRIVER,med21FileName)
fieldIntVectorOnCells.write(MED_DRIVER,med21FileName)


idVtk = fieldDoubleScalarOnNodes.addDriver(VTK_DRIVER,vtkFileName,fieldDoubleScalarOnNodes.getName())
fieldDoubleScalarOnNodes.writeAppend(idVtk)

idVtk = fieldIntScalarOnNodes.addDriver(VTK_DRIVER,vtkFileName,fieldIntScalarOnNodes.getName())
fieldIntScalarOnNodes.writeAppend(idVtk)

idVtk = fieldDoubleVectorOnNodes.addDriver(VTK_DRIVER,vtkFileName,fieldDoubleVectorOnNodes.getName())
fieldDoubleVectorOnNodes.writeAppend(idVtk)

idVtk = fieldIntVectorOnNodes.addDriver(VTK_DRIVER,vtkFileName,fieldIntVectorOnNodes.getName())
fieldIntVectorOnNodes.writeAppend(idVtk)

idVtk = fieldDoubleScalarOnCells.addDriver(VTK_DRIVER,vtkFileName,fieldDoubleScalarOnCells.getName())
fieldDoubleScalarOnCells.writeAppend(idVtk)

idVtk = fieldIntScalarOnCells.addDriver(VTK_DRIVER,vtkFileName,fieldIntScalarOnCells.getName())
fieldIntScalarOnCells.writeAppend(idVtk)

idVtk = fieldDoubleVectorOnCells.addDriver(VTK_DRIVER,vtkFileName,fieldDoubleVectorOnCells.getName())
fieldDoubleVectorOnCells.writeAppend(idVtk)

idVtk = fieldIntVectorOnCells.addDriver(VTK_DRIVER,vtkFileName,fieldIntVectorOnCells.getName())
fieldIntVectorOnCells.writeAppend(idVtk)
