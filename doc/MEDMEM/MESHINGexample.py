#  -*- coding: iso-8859-1 -*-
#  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
#
#  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
#  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
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

coordinates = []

coordinate = [0.0, 0.0, 0.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [0.0, 0.0, 1.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [2.0, 0.0, 1.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [0.0, 2.0, 1.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [-2.0, 0.0, 1.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [0.0, -2.0, 1.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [1.0, 1.0, 2.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [-1.0, 1.0, 2.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [-1.0, -1.0, 2.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [1.0, -1.0, 2.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [1.0, 1.0, 3.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [-1.0, 1.0, 3.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [-1.0, -1.0, 3.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [1.0, -1.0, 3.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [1.0, 1.0, 4.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [-1.0, 1.0, 4.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [-1.0, -1.0, 4.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [1.0, -1.0, 4.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])
coordinate = [0.0, 0.0, 5.0]
coordinates.append(coordinate[0])
coordinates.append(coordinate[1])
coordinates.append(coordinate[2])

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

connectivityTetra = []

connectivity =  [1,2,3,6]
connectivityTetra.append(connectivity[0])
connectivityTetra.append(connectivity[1])
connectivityTetra.append(connectivity[2])
connectivityTetra.append(connectivity[3])
connectivity =  [1,2,4,3]
connectivityTetra.append(connectivity[0])
connectivityTetra.append(connectivity[1])
connectivityTetra.append(connectivity[2])
connectivityTetra.append(connectivity[3])
connectivity =  [1,2,5,4]
connectivityTetra.append(connectivity[0])
connectivityTetra.append(connectivity[1])
connectivityTetra.append(connectivity[2])
connectivityTetra.append(connectivity[3])
connectivity =  [1,2,6,5]
connectivityTetra.append(connectivity[0])
connectivityTetra.append(connectivity[1])
connectivityTetra.append(connectivity[2])
connectivityTetra.append(connectivity[3])
connectivity =  [2,7,4,3]
connectivityTetra.append(connectivity[0])
connectivityTetra.append(connectivity[1])
connectivityTetra.append(connectivity[2])
connectivityTetra.append(connectivity[3])
connectivity =  [2,8,5,4]
connectivityTetra.append(connectivity[0])
connectivityTetra.append(connectivity[1])
connectivityTetra.append(connectivity[2])
connectivityTetra.append(connectivity[3])
connectivity =  [2,9,6,5]
connectivityTetra.append(connectivity[0])
connectivityTetra.append(connectivity[1])
connectivityTetra.append(connectivity[2])
connectivityTetra.append(connectivity[3])
connectivity =  [2,10,3,6]
connectivityTetra.append(connectivity[0])
connectivityTetra.append(connectivity[1])
connectivityTetra.append(connectivity[2])
connectivityTetra.append(connectivity[3])
connectivity =  [2,7,3,10]
connectivityTetra.append(connectivity[0])
connectivityTetra.append(connectivity[1])
connectivityTetra.append(connectivity[2])
connectivityTetra.append(connectivity[3])
connectivity =  [2,8,4,7]
connectivityTetra.append(connectivity[0])
connectivityTetra.append(connectivity[1])
connectivityTetra.append(connectivity[2])
connectivityTetra.append(connectivity[3])
connectivity =  [2,9,5,8]
connectivityTetra.append(connectivity[0])
connectivityTetra.append(connectivity[1])
connectivityTetra.append(connectivity[2])
connectivityTetra.append(connectivity[3])
connectivity =  [2,10,6,9]
connectivityTetra.append(connectivity[0])
connectivityTetra.append(connectivity[1])
connectivityTetra.append(connectivity[2])
connectivityTetra.append(connectivity[3])

myMeshing.setConnectivity(connectivityTetra,entity,types[0])

connectivityPyra = []
connectivity =  [7,8,9,10,2]
connectivityPyra.append(connectivity[0])
connectivityPyra.append(connectivity[1])
connectivityPyra.append(connectivity[2])
connectivityPyra.append(connectivity[3])
connectivityPyra.append(connectivity[4])
connectivity =  [15,18,17,16,19]
connectivityPyra.append(connectivity[0])
connectivityPyra.append(connectivity[1])
connectivityPyra.append(connectivity[2])
connectivityPyra.append(connectivity[3])
connectivityPyra.append(connectivity[4])

myMeshing.setConnectivity(connectivityPyra,entity,types[1])

connectivityHexa = []
connectivity =  [11,12,13,14,7,8,9,10]
connectivityHexa.append(connectivity[0])
connectivityHexa.append(connectivity[1])
connectivityHexa.append(connectivity[2])
connectivityHexa.append(connectivity[3])
connectivityHexa.append(connectivity[4])
connectivityHexa.append(connectivity[5])
connectivityHexa.append(connectivity[6])
connectivityHexa.append(connectivity[7])
connectivity =  [15,16,17,18,11,12,13,14]
connectivityHexa.append(connectivity[0])
connectivityHexa.append(connectivity[1])
connectivityHexa.append(connectivity[2])
connectivityHexa.append(connectivity[3])
connectivityHexa.append(connectivity[4])
connectivityHexa.append(connectivity[5])
connectivityHexa.append(connectivity[6])
connectivityHexa.append(connectivity[7])

myMeshing.setConnectivity(connectivityHexa,entity,types[2])

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

connectivityTria = []
connectivity =  [1,4,3]
connectivityTria.append(connectivity[0])
connectivityTria.append(connectivity[1])
connectivityTria.append(connectivity[2])
connectivity =  [1,5,4]
connectivityTria.append(connectivity[0])
connectivityTria.append(connectivity[1])
connectivityTria.append(connectivity[2])
connectivity =  [1,6,5]
connectivityTria.append(connectivity[0])
connectivityTria.append(connectivity[1])
connectivityTria.append(connectivity[2])
connectivity =  [1,3,6]
connectivityTria.append(connectivity[0])
connectivityTria.append(connectivity[1])
connectivityTria.append(connectivity[2])

myMeshing.setConnectivity(connectivityTria,entity,types[0])

connectivityQuad = []
connectivity =  [7,8,9,10]
connectivityQuad.append(connectivity[0])
connectivityQuad.append(connectivity[1])
connectivityQuad.append(connectivity[2])
connectivityQuad.append(connectivity[3])
connectivity =  [11,12,13,14]
connectivityQuad.append(connectivity[0])
connectivityQuad.append(connectivity[1])
connectivityQuad.append(connectivity[2])
connectivityQuad.append(connectivity[3])
connectivity =  [11,7,8,12]
connectivityQuad.append(connectivity[0])
connectivityQuad.append(connectivity[1])
connectivityQuad.append(connectivity[2])
connectivityQuad.append(connectivity[3])
connectivity =  [12,8,9,13]
connectivityQuad.append(connectivity[0])
connectivityQuad.append(connectivity[1])
connectivityQuad.append(connectivity[2])
connectivityQuad.append(connectivity[3])

myMeshing.setConnectivity(connectivityQuad,entity,types[1])

meshDimension = spaceDimension # because there 3D cells in the mesh
myMeshing.setMeshDimension(meshDimension)

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

# saving of the generated mesh in MED 2.1, 2.2 and VTK format

medFileVersion = getMedFileVersionForWriting()
print "Med File Version For Writing ",medFileVersion

if (medFileVersion == V22):
    setMedFileVersionForWriting(V21)

idMedV21 = myMeshing.addDriver(MED_DRIVER,med21FileName,myMeshing.getName())
myMeshing.write(idMedV21)

medFileVersion = getMedFileVersionForWriting()
if (medFileVersion == V21):
    setMedFileVersionForWriting(V22)

idMedV22 = myMeshing.addDriver(MED_DRIVER,med22FileName,myMeshing.getName())
myMeshing.write(idMedV22)

idVtk = myMeshing.addDriver(VTK_DRIVER,vtkFileName,myMeshing.getName())
myMeshing.write(idVtk)

# we build now 8 fields : 4 fields double (integer) :
#                         2 fields on nodes (cells) :
#                         1 scalar (vector)

supportOnNodes = SUPPORT(myMeshing,"On_All_Nodes",MED_NODE)
numberOfNodes = supportOnNodes.getNumberOfElements(MED_ALL_ELEMENTS)

supportOnCells = SUPPORT(myMeshing,"On_All_Cells",MED_CELL)
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

medFileVersion = getMedFileVersionForWriting()
print "Med File Version For Writing ",medFileVersion

if (medFileVersion == V22):
    setMedFileVersionForWriting(V21)

idMedV21 = fieldDoubleScalarOnNodes.addDriver(MED_DRIVER,med21FileName,fieldDoubleScalarOnNodes.getName())
fieldDoubleScalarOnNodes.write(idMedV21)

idMedV21 = fieldIntScalarOnNodes.addDriver(MED_DRIVER,med21FileName,fieldIntScalarOnNodes.getName())
fieldIntScalarOnNodes.write(idMedV21)

idMedV21 = fieldDoubleVectorOnNodes.addDriver(MED_DRIVER,med21FileName,fieldDoubleVectorOnNodes.getName())
fieldDoubleVectorOnNodes.write(idMedV21)

idMedV21 = fieldIntVectorOnNodes.addDriver(MED_DRIVER,med21FileName,fieldIntVectorOnNodes.getName())
fieldIntVectorOnNodes.write(idMedV21)

idMedV21 = fieldDoubleScalarOnCells.addDriver(MED_DRIVER,med21FileName,fieldDoubleScalarOnCells.getName())
fieldDoubleScalarOnCells.write(idMedV21)

idMedV21 = fieldIntScalarOnCells.addDriver(MED_DRIVER,med21FileName,fieldIntScalarOnCells.getName())
fieldIntScalarOnCells.write(idMedV21)

idMedV21 = fieldDoubleVectorOnCells.addDriver(MED_DRIVER,med21FileName,fieldDoubleVectorOnCells.getName())
fieldDoubleVectorOnCells.write(idMedV21)

idMedV21 = fieldIntVectorOnCells.addDriver(MED_DRIVER,med21FileName,fieldIntVectorOnCells.getName())
fieldIntVectorOnCells.write(idMedV21)

medFileVersion = getMedFileVersionForWriting()
if (medFileVersion == V21):
    setMedFileVersionForWriting(V22)

idMedV22 = fieldDoubleScalarOnNodes.addDriver(MED_DRIVER,med22FileName,fieldDoubleScalarOnNodes.getName())
fieldDoubleScalarOnNodes.write(idMedV22)

idMedV22 = fieldIntScalarOnNodes.addDriver(MED_DRIVER,med22FileName,fieldIntScalarOnNodes.getName())
fieldIntScalarOnNodes.write(idMedV22)

idMedV22 = fieldDoubleVectorOnNodes.addDriver(MED_DRIVER,med22FileName,fieldDoubleVectorOnNodes.getName())
fieldDoubleVectorOnNodes.write(idMedV22)

idMedV22 = fieldIntVectorOnNodes.addDriver(MED_DRIVER,med22FileName,fieldIntVectorOnNodes.getName())
fieldIntVectorOnNodes.write(idMedV22)

idMedV22 = fieldDoubleScalarOnCells.addDriver(MED_DRIVER,med22FileName,fieldDoubleScalarOnCells.getName())
fieldDoubleScalarOnCells.write(idMedV22)

idMedV22 = fieldIntScalarOnCells.addDriver(MED_DRIVER,med22FileName,fieldIntScalarOnCells.getName())
fieldIntScalarOnCells.write(idMedV22)

idMedV22 = fieldDoubleVectorOnCells.addDriver(MED_DRIVER,med22FileName,fieldDoubleVectorOnCells.getName())
fieldDoubleVectorOnCells.write(idMedV22)

idMedV22 = fieldIntVectorOnCells.addDriver(MED_DRIVER,med22FileName,fieldIntVectorOnCells.getName())
fieldIntVectorOnCells.write(idMedV22)

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
