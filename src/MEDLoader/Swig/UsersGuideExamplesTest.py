#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2019  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
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

from MEDLoader import *
import os

from MEDLoaderDataForTest import MEDLoaderDataForTest
m = MEDLoaderDataForTest.build1DMesh_1()
m.setName("mesh2")
#! [UG_ReadMeshFromFile_3]
m.checkConsecutiveCellTypesForMEDFileFrmt()
#! [UG_ReadMeshFromFile_3]
#! [UG_ReadMeshFromFile_0]
from MEDLoader import WriteMesh
WriteMesh("file2.med",m,True)
#! [UG_ReadMeshFromFile_0]
#! [UG_ReadMeshFromFile_1]
from MEDLoader import ReadMeshFromFile
m=ReadMeshFromFile("file2.med")
#! [UG_ReadMeshFromFile_1]
#! [UG_ReadMeshFromFile_2]
m=ReadMeshFromFile("file2.med","mesh2")
assert(m.getName()=="mesh2")
#! [UG_ReadMeshFromFile_2]

mesh3D=MEDLoaderDataForTest.build3DMesh_1()
mesh2D=mesh3D.computeSkin()
mesh1D=mesh2D.computeSkin()
#! [UG_ReadMeshFromFile_4]
from MEDLoader import MEDFileUMesh
mm=MEDFileUMesh.New()
mm.setMeshAtLevel(0,mesh3D)
mm.setMeshAtLevel(-1,mesh2D)
#! [UG_ReadMeshFromFile_4]
otherCoordArray=mesh3D.getCoords()
#! [UG_ReadMeshFromFile_5]
mm.setCoords(otherCoordArray)
#! [UG_ReadMeshFromFile_5]
#! [UG_ReadMeshFromFile_6]
groupNodes=DataArrayInt([1,3,4,5]);  groupNodes.setName("myNodes")
groupFaces=DataArrayInt([12,13,15]); groupFaces.setName("myFaces")
mm.addGroup(1,groupNodes)
mm.addGroup(-1,groupFaces)
#! [UG_ReadMeshFromFile_6]
oldName,newName="myNodes","myNodes2"
oldFamName,newFamName="Family_2","Family_3"
#! [UG_ReadMeshFromFile_7]
mm.changeGroupName(oldName,newName)
mm.changeFamilyName(oldFamName,newFamName)
#! [UG_ReadMeshFromFile_7]
#! [UG_ReadMeshFromFile_8]
mm.write("file.med",2)
#! [UG_ReadMeshFromFile_8]
mm=MEDFileUMesh.New()
#! [UG_ReadMeshFromFile_9]
mm.setMeshAtLevel(0,mesh3D,True)
#! [UG_ReadMeshFromFile_9]
pass

from MEDLoaderDataForTest import MEDLoaderDataForTest
f=MEDLoaderDataForTest.buildVecFieldOnGauss_1();
f.setName("Field1")
#! [UG_ReadField_0]
from MEDLoader import WriteField
WriteField("file.med",f,True)
#! [UG_ReadField_0]
#! [UG_ReadField_1]
from MEDLoader import ReadField
f=ReadField("file.med")
#! [UG_ReadField_1]
#! [UG_ReadField_2]
from MEDLoader import GetAllFieldNames
print(GetAllFieldNames("file.med"))
#! [UG_ReadField_2]
#! [UG_ReadField_3]
f=ReadField("file.med","Field1")
#! [UG_ReadField_3]
#! [UG_ReadField_4]
from MEDLoader import GetAllFieldIterations
print(GetAllFieldIterations("file.med","Field1"))
#! [UG_ReadField_4]
#! [UG_ReadField_5]
ts0,ts1=1,5
f=ReadField("file.med","Field1",ts0,ts1)
#! [UG_ReadField_5]
fs = [ f ]
#! [UG_ReadField_6]
m=fs[0].getMesh()
WriteMesh("file5.med",m,True)
for f in fs:
    assert(f.getMesh().getHiddenCppPointer()==m.getHiddenCppPointer())
    # extra line to insist on the fact that
    WriteFieldUsingAlreadyWrittenMesh("file5.med",f)
#! [UG_ReadField_6]

from MEDLoaderDataForTest import MEDLoaderDataForTest
fname="PyExamples1.med"
meshName="mesh"
fieldName="FieldOnAll"
iteration=3
order=4
m=MEDLoaderDataForTest.build2DMesh_3()
m.setName(meshName)
f=m.getMeasureField(False)
f=f.buildNewTimeReprFromThis(ONE_TIME,False)
f.setTime(5.5,iteration,order)
f.setName(fieldName)
mesh=m
field=f
level=0
#! [UG_RWFieldAdv_0]
from MEDLoader import MEDFileUMesh, MEDFileField1TS
mm=MEDFileUMesh.New()
mm.setMeshAtLevel(0,mesh)
ff=MEDFileField1TS.New()
ff.setFieldNoProfileSBT(field)
mm.write(fname,2)
ff.write(fname,0)
#! [UG_RWFieldAdv_0]
#! [UG_RWFieldAdv_1]
profile=DataArrayInt([1,3,7]); profile.setName("pfl137")
fieldPartial=field[profile]
fieldPartial.setName("fieldPartial")
ff.setFieldProfile(fieldPartial,mm,level,profile)
ff.write(fname,0)
#! [UG_RWFieldAdv_1]
#! [UG_RWFieldAdv_2]
ff=MEDFileField1TS.New(fname,fieldName,iteration,order)
mm=MEDFileMesh.New(fname)
# you can choose an appropriate method
field=ff.field(mm)
field=ff.getFieldAtLevel(ON_CELLS,level)
field=ff.getFieldOnMeshAtLevel(ON_CELLS,level,mm)
#! [UG_RWFieldAdv_2]
#! [UG_RWFieldAdv_3]
maxDim,maxRelDims=ff.getNonEmptyLevels()
#! [UG_RWFieldAdv_3]

fieldTS1 = f
fieldTS2 = f.deepCopy()
fieldTS2.setTime(4.5,iteration+1,order)
fieldPartialTS1 = fieldPartial
fieldPartialTS2 = fieldPartial.deepCopy()
fieldPartialTS1.setTime(4.0,iteration+2,order)
fieldPartialTS2.setTime(3.5,iteration+3,order)
fieldPartialTS1.setName( fieldTS1.getName() )
fieldPartialTS2.setName( fieldTS1.getName() )
fname="PyExamples2.med"
mm.write(fname,2)
#mm=MEDFileMesh.New(fname)
#! [UG_RWFieldAdv_4]
ff=MEDFileFieldMultiTS.New()
ff.appendFieldNoProfileSBT(fieldTS1)
ff.appendFieldNoProfileSBT(fieldTS2)
ff.appendFieldProfile(fieldPartialTS1,mm,level,profile)
ff.appendFieldProfile(fieldPartialTS2,mm,level,profile)
ff.write(fname,0)
#! [UG_RWFieldAdv_4]

#! [UG_RWFieldAdv_5]
mm=MEDFileMesh.New(fname)
ff=MEDFileFieldMultiTS.New(fname,fieldName)
for ff1TS in ff:
    iteration,order,time=ff1TS.getTime()
    # you can choose an appropriate method
    field=ff1TS.field(mm)
    field=ff1TS.getFieldAtLevel(ON_CELLS,level)
    field=ff1TS.getFieldOnMeshAtLevel(ON_CELLS,level,mm)
#! [UG_RWFieldAdv_5]
