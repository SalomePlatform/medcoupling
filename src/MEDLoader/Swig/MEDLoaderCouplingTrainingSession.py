#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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
# Author : Anthony GEAY (CEA/DEN/DM2S/STMF/LGLS)

from MEDLoader import *
from MEDCouplingRemapper import *
import math, os

d=DataArrayDouble.New(6,2)
d[:,0]=3.
d[:, 1] = list(range(6))
d[:,1]*=math.pi/3.
d=d.fromPolarToCart()
d.setInfoOnComponents(["X [m]","Y [m]"])
print(d.getValues())
print(d)
print(d.magnitude().isUniform(3.,1e-12))
#
radius=3.
translationToPerform=[[0.,0.],[3./2.*radius,-radius*math.sqrt(3.)/2],[3./2.*radius,radius*math.sqrt(3.)/2],[0.,radius*math.sqrt(3.)],[-3./2.*radius,radius*math.sqrt(3.)/2],[-3./2.*radius,-radius*math.sqrt(3.)/2],[0.,-radius*math.sqrt(3.)]]
ds=len(translationToPerform)*[None]
for pos,t in enumerate(translationToPerform):
  ds[pos]=d[:]
  ds[pos]+=t
  pass
#
d2=DataArrayDouble.Aggregate(ds)
oldNbOfTuples=d2.getNumberOfTuples()
c,cI=d2.findCommonTuples(1e-12)
tmp=c[cI[0]:cI[0+1]]
print(tmp)
a=cI.deltaShiftIndex()
b=a-1
myNewNbOfTuples=oldNbOfTuples-sum(b.getValues())
o2n,newNbOfTuples=DataArrayInt.ConvertIndexArrayToO2N(oldNbOfTuples,c,cI)
print("Ai je trouve le bon resultat ? %s"%(str(myNewNbOfTuples==newNbOfTuples))) ; assert myNewNbOfTuples==newNbOfTuples
#
d3=d2.renumberAndReduce(o2n,newNbOfTuples)
n2o=o2n.invertArrayO2N2N2O(newNbOfTuples)
d3_bis=d2[n2o]
print("Ai je trouve le bon resultat (2) ? %s"%(str(d3.isEqual(d3_bis,1e-12)))) ; assert d3.isEqual(d3_bis,1e-12)
#
d3+=[3.3,4.4]
# d3 contains coordinates
m=MEDCouplingUMesh.New("My7hexagons",2)
m.setCoords(d3)
m.allocateCells(7)
for i in range(7):
  m.insertNextCell(NORM_POLYGON,o2n[6*i:6*(i+1)].getValues())
  pass
m.finishInsertingCells()
m.checkConsistencyLight()
#
m.writeVTK("My7hexagons.vtu")

########

coords=[0.,0.,0., 1.,1.,0., 1.,1.25,0., 1.,0.,0., 1.,1.5,0., 2.,0.,0., 2.,1.,0., 1.,2.,0., 0.,2.,0., 3.,1.,0.,
        3.,2.,0., 0.,1.,0., 1.,3.,0., 2.,2.,0., 2.,3.,0.,
        0.,0.,1., 1.,1.,1., 1.,1.25,1., 1.,0.,1., 1.,1.5,1., 2.,0.,1., 2.,1.,1., 1.,2.,1., 0.,2.,1., 3.,1.,1.,
        3.,2.,1., 0.,1.,1., 1.,3.,1., 2.,2.,1., 2.,3.,1.,
        0.,0.,2., 1.,1.,2., 1.,1.25,2., 1.,0.,2., 1.,1.5,2., 2.,0.,2., 2.,1.,2., 1.,2.,2., 0.,2.,2., 3.,1.,2.,
        3.,2.,2., 0.,1.,2., 1.,3.,2., 2.,2.,2., 2.,3.,2.,
        0.,0.,3., 1.,1.,3., 1.,1.25,3., 1.,0.,3., 1.,1.5,3., 2.,0.,3., 2.,1.,3., 1.,2.,3., 0.,2.,3., 3.,1.,3.,
        3.,2.,3., 0.,1.,3., 1.,3.,3., 2.,2.,3., 2.,3.,3.]
conn=[0,11,1,3,15,26,16,18,   1,2,4,7,13,6,-1,1,16,21,6,-1,6,21,28,13,-1,13,7,22,28,-1,7,4,19,22,-1,4,2,17,19,-1,2,1,16,17,-1,16,21,28,22,19,17,
      1,6,5,3,16,21,20,18,   13,10,9,6,28,25,24,21, 11,8,7,4,2,1,-1,11,26,16,1,-1,1,16,17,2,-1,2,17,19,4,-1,4,19,22,7,-1,7,8,23,22,-1,8,11,26,23,-1,26,16,17,19,22,23,
      7,12,14,13,22,27,29,28,  15,26,16,18,30,41,31,33, 16,17,19,22,28,21,-1,16,31,36,21,-1,21,36,43,28,-1,28,22,37,43,-1,22,19,34,37,-1,19,17,32,34,-1,17,16,31,32,-1,31,36,43,37,34,32,
      16,21,20,18,31,36,35,33,   28,25,24,21,43,40,39,36, 26,23,22,19,17,16,-1,26,41,31,16,-1,16,31,32,17,-1,17,32,34,19,-1,19,34,37,22,-1,22,23,38,37,-1,23,26,41,38,-1,41,31,32,34,37,38,
      22,27,29,28,37,42,44,43, 30,41,31,33,45,56,46,48,  31,32,34,37,43,36,-1,31,46,51,36,-1,36,51,58,43,-1,43,37,52,58,-1,37,34,49,52,-1,34,32,47,49,-1,32,31,46,47,-1,46,51,58,52,49,47,
      31,36,35,33,46,51,50,48,  43,40,39,36,58,55,54,51, 41,38,37,34,32,31,-1,41,56,46,31,-1,31,46,47,32,-1,32,47,49,34,-1,34,49,52,37,-1,37,38,53,52,-1,38,41,56,53,-1,56,46,47,49,52,53,
      37,42,44,43,52,57,59,58]
mesh3D=MEDCouplingUMesh.New("mesh3D",3);
mesh3D.allocateCells(18);
mesh3D.insertNextCell(NORM_HEXA8,conn[0:8]); mesh3D.insertNextCell(NORM_POLYHED,conn[8:51]); mesh3D.insertNextCell(NORM_HEXA8,conn[51:59]); mesh3D.insertNextCell(NORM_HEXA8,conn[59:67]); mesh3D.insertNextCell(NORM_POLYHED,conn[67:110]); mesh3D.insertNextCell(NORM_HEXA8,conn[110:118]);
mesh3D.insertNextCell(NORM_HEXA8,conn[118:126]); mesh3D.insertNextCell(NORM_POLYHED,conn[126:169]); mesh3D.insertNextCell(NORM_HEXA8,conn[169:177]); mesh3D.insertNextCell(NORM_HEXA8,conn[177:185]); mesh3D.insertNextCell(NORM_POLYHED,conn[185:228]); mesh3D.insertNextCell(NORM_HEXA8,conn[228:236]);
mesh3D.insertNextCell(NORM_HEXA8,conn[236:244]); mesh3D.insertNextCell(NORM_POLYHED,conn[244:287]); mesh3D.insertNextCell(NORM_HEXA8,conn[287:295]); mesh3D.insertNextCell(NORM_HEXA8,conn[295:303]); mesh3D.insertNextCell(NORM_POLYHED,conn[303:346]); mesh3D.insertNextCell(NORM_HEXA8,conn[346:354]);
mesh3D.finishInsertingCells();
myCoords=DataArrayDouble.New(coords,60,3);
myCoords.setInfoOnComponents(["X [m]","Y [m]","Z [m]"])
mesh3D.setCoords(myCoords);
mesh3D.orientCorrectlyPolyhedrons()
mesh3D.sortCellsInMEDFileFrmt()
mesh3D.checkConsistencyLight()
renum = DataArrayInt.New(60) ; renum[:15] = list(range(15, 30)) ; renum[15:30] = list(range(15)) ; renum[30:45] = list(range(45, 60)) ; renum[45:] = list(range(30, 45))
mesh3D.renumberNodes(renum,60)
#
mesh3D.getCoords()[:]*=100.
mesh3D.getCoords().setInfoOnComponents(["X [cm]","Y [cm]","Z [cm]"])
#
zLev=mesh3D.getCoords()[:,2]
zLev = zLev.getDifferentValues(1e-12)
zLev.sort()
#
tmp,cellIdsSol1=mesh3D.buildSlice3D([0.,0.,(zLev[1]+zLev[2])/2],[0.,0.,1.],1e-12)
bary=mesh3D.computeCellCenterOfMass()
baryZ=bary[:,2]
cellIdsSol2=baryZ.findIdsInRange(zLev[1],zLev[2])
nodeIds=mesh3D.findNodesOnPlane([0.,0.,zLev[0]],[0.,0.,1.],1e-10)
mesh2D=mesh3D.buildFacePartOfMySelfNode(nodeIds,True)
extMesh=MEDCouplingMappedExtrudedMesh.New(mesh3D,mesh2D,0)
cellIdsSol3=extMesh.getMesh3DIds()[mesh2D.getNumberOfCells():2*mesh2D.getNumberOfCells()]
for i in range(3):
  exec("print( cellIdsSol%s.getValues())"%(i+1))
#
mesh3DPart=mesh3D[cellIdsSol2] # equivalent to mesh3DPart=mesh3D.buildPartOfMySelf(cellIdsSol2,True)
mesh3DPart.zipCoords()
print(mesh3DPart.checkConsecutiveCellTypesAndOrder([NORM_HEXA8,NORM_POLYHED])) ; assert mesh3DPart.checkConsecutiveCellTypesAndOrder([NORM_HEXA8,NORM_POLYHED])
print(mesh3DPart.checkConsecutiveCellTypes()) ; assert mesh3DPart.checkConsecutiveCellTypes()
#print mesh3DPart.advancedRepr()
#
baryXY=bary[:,[0,1]]
baryXY-=[250.,150.]
magn=baryXY.magnitude()
cellIds2Sol1=magn.findIdsInRange(0.,1e-12)
#
bary2=mesh2D.computeCellCenterOfMass()[:,[0,1]]
bary2-=[250.,150.]
magn=bary2.magnitude()
ids=magn.findIdsInRange(0.,1e-12)
idStart=int(ids) # ids is assumed to contain only one value, if not an exception is thrown
cellIds2Sol2 = extMesh.getMesh3DIds()[list(range(idStart, mesh3D.getNumberOfCells(), mesh2D.getNumberOfCells()))]
#
mesh3DSlice2=mesh3D[cellIds2Sol1]
mesh3DSlice2.zipCoords()
#
mesh3DSlice2bis=mesh3DSlice2.deepCopy()
mesh3DSlice2bis.translate([0.,1000.,0.])
mesh3DSlice2All=MEDCouplingUMesh.MergeUMeshes([mesh3DSlice2,mesh3DSlice2bis])
mesh3DSlice2All.writeVTK("mesh3DSlice2All.vtu")
#
mesh3DSurf,desc,descIndx,revDesc,revDescIndx=mesh3D.buildDescendingConnectivity()
numberOf3DCellSharing=revDescIndx.deltaShiftIndex()
cellIds=numberOf3DCellSharing.findIdsNotEqual(1)
mesh3DSurfInside=mesh3DSurf[cellIds]
mesh3DSurfInside.writeVTK("mesh3DSurfInside.vtu")

######

xarr=DataArrayDouble.New(11,1)
xarr.iota(0.)
cmesh=MEDCouplingCMesh.New()
cmesh.setCoords(xarr,xarr,xarr)
mesh=cmesh.buildUnstructured()
mesh.convertToPolyTypes(DataArrayInt.Range(0,mesh.getNumberOfCells(),2))
#

f=mesh.fillFromAnalytic(ON_CELLS,1,"(x-5.)*(x-5.)+(y-5.)*(y-5.)+(z-5.)*(z-5.)")
f.setName("MyField")
#
f2=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
f2.setMesh(mesh)
f2.setName("MyField2")
f2.fillFromAnalytic(1,"(x-5.)*(x-5.)+(y-5.)*(y-5.)+(z-5.)*(z-5.)")
print("f and f2 are equal : %s"%(f.isEqualWithoutConsideringStr(f2,1e-13,1e-12))) ; assert f.isEqualWithoutConsideringStr(f2,1e-13,1e-12)
#
ids1=f.getArray().findIdsInRange(0.,5.)
fPart1=f.buildSubPart(ids1)
ids2=f.getArray().findIdsInRange(50.,1.e300)
fPart2=f.buildSubPart(ids2)
#Renumbering cells to follow MED file
fPart1Cpy=fPart1.deepCopy()
o2n=fPart1Cpy.getMesh().sortCellsInMEDFileFrmt()
fPart1Cpy.getArray().renumberInPlace(o2n)
#Check that fPart1Cpy and fPart1 are the same
fPart1Cpy.substractInPlaceDM(fPart1,12,1e-12)
fPart1Cpy.getArray().abs()
print("Fields are the same ? %s"%(fPart1Cpy.getArray().accumulate()[0]<1e-12)) ; assert fPart1Cpy.getArray().accumulate()[0]<1e-12
#
fPart12=MEDCouplingFieldDouble.MergeFields([fPart1,fPart2])
# evaluation on points
bary=fPart12.getMesh().computeCellCenterOfMass()
arr1=fPart12.getValueOnMulti(bary)
arr2=f.getValueOnMulti(bary)
delta=arr1-arr2
delta.abs()
print("Check OK : %s"%(delta.accumulate()[0]<1e-12)) ; assert delta.accumulate()[0]<1e-12
#
print(abs(fPart12.integral(0,True)-fPart12.getArray().accumulate()[0])<1e-10) ; assert abs(fPart12.integral(0,True)-fPart12.getArray().accumulate()[0])<1e-10
fPart12.getMesh().scale([0.,0.,0.],1.2)
print(abs(fPart12.integral(0,True)-fPart12.getArray().accumulate()[0]*1.2*1.2*1.2)<1e-8) ; assert abs(fPart12.integral(0,True)-fPart12.getArray().accumulate()[0]*1.2*1.2*1.2)<1e-8
# Explosion of field
fVec=mesh.fillFromAnalytic(ON_CELLS,3,"(x-5.)*IVec+(y-5.)*JVec+(z-5.)*KVec")
fVecPart1=fVec.buildSubPart(ids1)
fVecPart1.setName("fVecPart1")
cells=fPart1.getMesh().getNumberOfCells()*[None]
for icell,vec in enumerate(fVecPart1.getArray()):
  m=fPart1.getMesh()[[icell]]
  m.zipCoords()
  m.translate(vec)
  cells[icell]=m
  pass
meshFVecPart1Exploded=MEDCouplingUMesh.MergeUMeshes(cells)
fPart1.setMesh(meshFVecPart1Exploded)

####

targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ];
targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4];
targetMesh=MEDCouplingUMesh.New("MyMesh",2);
targetMesh.allocateCells(5);
targetMesh.insertNextCell(NORM_TRI3,3,targetConn[4:7]);
targetMesh.insertNextCell(NORM_TRI3,3,targetConn[7:10]);
targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[0:4]);
targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[10:14]);
targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[14:18]);
targetMesh.finishInsertingCells();
myCoords=DataArrayDouble.New(targetCoords,9,2);
myCoords.setInfoOnComponents(["X [km]","YY [mm]"])
targetMesh.setCoords(myCoords);
#
WriteUMesh("TargetMesh.med",targetMesh,True)
#
meshRead=ReadUMeshFromFile("TargetMesh.med",targetMesh.getName(),0)
print("Is the mesh read in file equals targetMesh ? %s"%(meshRead.isEqual(targetMesh,1e-12))) ; assert meshRead.isEqual(targetMesh,1e-12)
#
f=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
f.setTime(5.6,7,8)
f.setArray(targetMesh.computeCellCenterOfMass())
f.setMesh(targetMesh)
f.setName("AFieldName")
WriteField("MyFirstField.med",f,True)
#
f2=ReadFieldCell("MyFirstField.med",f.getMesh().getName(),0,f.getName(),7,8)
print("Is the field read in file equals f ? %s"%(f2.isEqual(f,1e-12,1e-12))) ; assert f2.isEqual(f,1e-12,1e-12)
#
WriteUMesh("MySecondField.med",f.getMesh(),True)
WriteFieldUsingAlreadyWrittenMesh("MySecondField.med",f)
#
f2=f.clone(True)
f2.getArray()[:]*=2.0
f2.setTime(7.8,9,10)
WriteFieldUsingAlreadyWrittenMesh("MySecondField.med",f2)
#
f3=ReadFieldCell("MySecondField.med",f.getMesh().getName(),0,f.getName(),7,8)
print("Is the field read in file equals f ? %s"%(f.isEqual(f3,1e-12,1e-12))) ; assert f.isEqual(f3,1e-12,1e-12)
f4=ReadFieldCell("MySecondField.med",f.getMesh().getName(),0,f.getName(),9,10)
print("Is the field read in file equals f ? %s"%(f2.isEqual(f4,1e-12,1e-12))) ; assert f2.isEqual(f4,1e-12,1e-12)

#####

targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ];
targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4];
targetMesh=MEDCouplingUMesh.New("MyMesh",2);
targetMesh.allocateCells(5);
targetMesh.insertNextCell(NORM_TRI3,3,targetConn[4:7]);
targetMesh.insertNextCell(NORM_TRI3,3,targetConn[7:10]);
targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[0:4]);
targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[10:14]);
targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[14:18]);
targetMesh.finishInsertingCells();
myCoords=DataArrayDouble.New(targetCoords,9,2);
myCoords.setInfoOnComponents(["X [km]","YY [mm]"])
targetMesh.setCoords(myCoords);
#
targetMeshConsti=targetMesh.buildDescendingConnectivity()[0]
targetMesh1=targetMeshConsti[[3,4,7,8]]
targetMesh1.setName(targetMesh.getName())
#
meshMEDFile=MEDFileUMesh.New()
meshMEDFile.setMeshAtLevel(0,targetMesh)
meshMEDFile.setMeshAtLevel(-1,targetMesh1)
# Some groups on cells Level 0
grp0_0=DataArrayInt.New([0,1,3]) ; grp0_0.setName("grp0_Lev0")
grp1_0=DataArrayInt.New([1,2,3,4]) ; grp1_0.setName("grp1_Lev0")
meshMEDFile.setGroupsAtLevel(0,[grp0_0,grp1_0])
# Some groups on cells Level -1
grp0_M1=DataArrayInt.New([0,1]) ; grp0_M1.setName("grp0_LevM1")
grp1_M1=DataArrayInt.New([0,1,2]) ; grp1_M1.setName("grp1_LevM1")
grp2_M1=DataArrayInt.New([1,2,3]) ; grp2_M1.setName("grp2_LevM1")
meshMEDFile.setGroupsAtLevel(-1,[grp0_M1,grp1_M1,grp2_M1])
#
meshMEDFile.write("TargetMesh2.med",2) # 2 stands for write from scratch
#
meshMEDFileRead=MEDFileMesh.New("TargetMesh2.med")
meshRead0=meshMEDFileRead.getMeshAtLevel(0)
meshRead1=meshMEDFileRead.getMeshAtLevel(-1)
print("Is the mesh at level 0 read in file equals targetMesh ? %s"%(meshRead0.isEqual(targetMesh,1e-12))) ; assert meshRead0.isEqual(targetMesh,1e-12)
print("Is the mesh at level -1 read in file equals targetMesh ? %s"%(meshRead1.isEqual(targetMesh1,1e-12))) ; assert meshRead1.isEqual(targetMesh1,1e-12)
#
print(meshMEDFileRead.getGrpNonEmptyLevels("grp0_Lev0"))
grp0_0_read=meshMEDFileRead.getGroupArr(0,"grp0_Lev0")
print("Is group \"grp0_Lev0\" are the same ? %s"%(grp0_0_read.isEqual(grp0_0))) ; assert grp0_0_read.isEqual(grp0_0)
#
# Fields
#
f=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
f.setTime(5.6,7,8)
f.setArray(targetMesh.computeCellCenterOfMass())
f.setMesh(targetMesh)
f.setName("AFieldName")
#
fMEDFile=MEDFileField1TS.New()
fMEDFile.setFieldNoProfileSBT(f)
#
fMEDFile.write("TargetMesh2.med",0) # 0 is very important here because we want to append to TargetMesh2.med and not to scratch it
#
fMEDFileRead=MEDFileField1TS.New("TargetMesh2.med",f.getName(),7,8)
fRead1=fMEDFileRead.getFieldOnMeshAtLevel(ON_CELLS,0,meshMEDFileRead) # fastest method. No read in file.
fRead2=fMEDFileRead.getFieldAtLevel(ON_CELLS,0) # basic method like, mesh is reread in file...
print("Does the field f remains the same using fast method ? %s"%(fRead1.isEqual(f,1e-12,1e-12))) ; assert fRead1.isEqual(f,1e-12,1e-12)
print("Does the field f remains the same using slow method ? %s"%(fRead2.isEqual(f,1e-12,1e-12))) ; assert fRead2.isEqual(f,1e-12,1e-12)
#
# Writing and Reading fields on profile using MEDLoader advanced API
#
pfl=DataArrayInt.New([1,2,3]) ; pfl.setName("My1stPfl")
fPart=f.buildSubPart(pfl)
fPart.setName("fPart")
#
fMEDFile2=MEDFileField1TS.New()
fMEDFile2.setFieldProfile(fPart,meshMEDFileRead,0,pfl)
fMEDFile2.write("TargetMesh2.med",0) # 0 is very important here because we want to append to TargetMesh2.med and not to scratch it
#
fMEDFileRead2=MEDFileField1TS.New("TargetMesh2.med",fPart.getName(),7,8)
fPartRead,pflRead=fMEDFileRead2.getFieldWithProfile(ON_CELLS,0,meshMEDFileRead)
print(fPartRead.isEqualWithoutConsideringStr(fPart.getArray(),1e-12)) ; assert fPartRead.isEqualWithoutConsideringStr(fPart.getArray(),1e-12)
print(pflRead.isEqualWithoutConsideringStr(pfl)) ; assert pflRead.isEqualWithoutConsideringStr(pfl)

#####

m0=MEDCouplingCMesh()
arr=DataArrayDouble(31,1) ; arr.iota(0.)
m0.setCoords(arr,arr)
m0=m0.buildUnstructured()
m00=m0[::2] ; m00.simplexize(0) ; m01=m0[1::2]
m0=MEDCouplingUMesh.MergeUMeshes([m00,m01])
m0.getCoords()[:]*=1/15.
m0.setName("mesh")
#
CellField=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; CellField.setTime(5.6,5,6) ; CellField.setMesh(m0)
CellField.setName("CellField")
CellField.fillFromAnalytic(1,"exp(-((x-1)*(x-1)+(y-1)*(y-1)))") ; CellField.getArray().setInfoOnComponent(0,"powercell [W]")
NodeField=MEDCouplingFieldDouble(ON_NODES,ONE_TIME) ; NodeField.setTime(5.6,5,6) ; NodeField.setMesh(m0)
NodeField.setName("NodeField")
NodeField.fillFromAnalytic(1,"exp(-((x-1)*(x-1)+(y-1)*(y-1)))") ; NodeField.getArray().setInfoOnComponent(0,"powernode [W]")
#
proc0=m0.getCellsInBoundingBox([(0.,0.4),(0.,0.4)],1e-10)
proc1=proc0.buildComplement(m0.getNumberOfCells())
#
NodeField0=NodeField[proc0] ; NodeField0.getMesh().setName(m0.getName()) ; CellField0=CellField[proc0] ; CellField0.setMesh(NodeField0.getMesh())
NodeField1=NodeField[proc1] ; NodeField1.getMesh().setName(m0.getName()) ; CellField1=CellField[proc1] ; CellField1.setMesh(NodeField1.getMesh())
#
proc0_fname="proc0.med"
WriteField(proc0_fname,NodeField0,True)
WriteFieldUsingAlreadyWrittenMesh(proc0_fname,CellField0)
proc1_fname="proc1.med"
WriteField(proc1_fname,NodeField1,True)
WriteFieldUsingAlreadyWrittenMesh(proc1_fname,CellField1)
#
CellField0_read=ReadFieldCell("proc0.med","mesh",0,"CellField",5,6)
CellField1_read=ReadFieldCell("proc1.med","mesh",0,"CellField",5,6)
CellField_read=MEDCouplingFieldDouble.MergeFields([CellField0_read,CellField1_read])
CellFieldCpy=CellField.deepCopy()
CellFieldCpy.substractInPlaceDM(CellField_read,10,1e-12)
CellFieldCpy.getArray().abs()
print(CellFieldCpy.getArray().isUniform(0.,1e-12))
#
NodeField0_read=ReadFieldNode("proc0.med","mesh",0,"NodeField",5,6)
NodeField1_read=ReadFieldNode("proc1.med","mesh",0,"NodeField",5,6)
NodeField_read=MEDCouplingFieldDouble.MergeFields([NodeField0_read,NodeField1_read])
NodeField_read.mergeNodes(1e-10)
NodeFieldCpy=NodeField.deepCopy()
NodeFieldCpy.mergeNodes(1e-10)
NodeFieldCpy.substractInPlaceDM(NodeField_read,10,1e-12)
print(NodeFieldCpy.getArray().isUniform(0.,1e-12)) ; assert NodeFieldCpy.getArray().isUniform(0.,1e-12)
#
fileNames=["proc0.med","proc1.med"]
msML=[MEDFileMesh.New(fname) for fname in fileNames]
fsML=[MEDFileFields.New(fname) for fname in fileNames]
mergeMLMesh=MEDFileUMesh()
mergeMLFields=MEDFileFields()
for lev in msML[0].getNonEmptyLevels():
    o2nML=len(msML[0].getNonEmptyLevels())*[None]
    cs=[mML.getCoords() for mML in msML]
    mergeMLMesh.setCoords(DataArrayDouble.Aggregate(cs))
    ms=[mML.getMeshAtLevel(lev) for mML in msML]
    m=MEDCouplingUMesh.MergeUMeshes(ms) ; m.setCoords(mergeMLMesh.getCoords())
    o2nML[lev]=m.sortCellsInMEDFileFrmt()
    mergeMLMesh.setMeshAtLevel(lev,m)
    pass
#
for fieldName in fsML[0].getFieldsNames():
    fmts=[fML[fieldName] for fML in fsML]
    mergeField=MEDFileFieldMultiTS()
    for dt,it,tim in fmts[0].getTimeSteps():
        fts=[fmt[dt,it] for fmt in fmts]
        arrs=len(fts)*[None]
        for typp in fts[0].getTypesOfFieldAvailable():
            arr1s=[]
            if typp==ON_CELLS:
               for ft in fts:
                   for geoTyp,smth in ft.getFieldSplitedByType():
                       if geoTyp!=NORM_ERROR:
                           smth1=[x for x in smth if x[0]==ON_CELLS]
                           arr2s=[ft.getUndergroundDataArray()[elt[1][0]:elt[1][1]] for elt in smth1]
                           arr1s.append(DataArrayDouble.Aggregate(arr2s))
                           pass
                       pass
                   pass
               pass
            else:
                for ft in fts:
                    smth=[x for x in ft.getFieldSplitedByType() if x[0]==NORM_ERROR]
                    arr2=DataArrayDouble.Aggregate([ft.getUndergroundDataArray()[elt[1][0][1][0]:elt[1][0][1][1]] for elt in smth])
                    arr1s.append(arr2)
                    pass
                pass
            arr=DataArrayDouble.Aggregate(arr1s)
            if typp==ON_CELLS:
               arr.renumberInPlace(o2nML[lev])
            mcf=MEDCouplingFieldDouble(typp,ONE_TIME) ; mcf.setName(fieldName) ; mcf.setTime(tim,dt,it) ; mcf.setArray(arr)
            mcf.setMesh(mergeMLMesh.getMeshAtLevel(lev)) ; mcf.checkConsistencyLight()
            mergeField.appendFieldNoProfileSBT(mcf)
            pass
        pass
    mergeMLFields.pushField(mergeField)
    pass
mergeMLMesh.write("merge.med",2)
mergeMLFields.write("merge.med",0)

#####

arr=DataArrayDouble(11) ; arr.iota(0)
trgMesh=MEDCouplingCMesh() ; trgMesh.setCoords(arr,arr) ; trgMesh=trgMesh.buildUnstructured()
#
arr=DataArrayDouble(21) ; arr.iota(0) ; arr*=0.5
srcMesh=MEDCouplingCMesh() ; srcMesh.setCoords(arr,arr) ; srcMesh=srcMesh.buildUnstructured()
#
tmp=srcMesh[:20] ; tmp.simplexize(0)
srcMesh=MEDCouplingUMesh.MergeUMeshes([tmp,srcMesh[20:]])
#
remap=MEDCouplingRemapper()
remap.prepare(srcMesh,trgMesh,"P0P0")
#
myMatrix=remap.getCrudeMatrix()
print(myMatrix) # pour voir a quoi elle ressemble
sumByRows=DataArrayDouble(len(myMatrix))
for i,wIt in enumerate(sumByRows):
  su=0.
  for it in myMatrix[i]:
    su+=myMatrix[i][it]
  wIt[0]=su
print("Does interpolation look OK ? %s"%(str(sumByRows.isUniform(1.,1e-12)))) ; assert sumByRows.isUniform(1.,1e-12)
#
srcField=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; srcField.setMesh(srcMesh)
srcField.fillFromAnalytic(1,"7-sqrt((x-5.)*(x-5.)+(y-5.)*(y-5.))") ; CellField.getArray().setInfoOnComponent(0,"powercell [W]")
#
#remap.transferField(srcField,1e300)
srcField.setNature(IntensiveMaximum)
trgFieldCV=remap.transferField(srcField,1e300)
#
print("IntensiveMaximum %lf == %lf"%(srcField.integral(True)[0],trgFieldCV.integral(True)[0])) ; assert abs(srcField.integral(True)[0]-trgFieldCV.integral(True)[0])<1e-6
print("IntensiveMaximum %lf != %lf"%(srcField.getArray().accumulate()[0],trgFieldCV.getArray().accumulate()[0])) ; assert abs(srcField.getArray().accumulate()[0]-trgFieldCV.getArray().accumulate()[0])>1e-6
#
srcField.setNature(ExtensiveMaximum)
trgFieldI=remap.transferField(srcField,1e300)
#
print("ExtensiveConservation %lf != %lf"%(srcField.integral(True)[0],trgFieldI.integral(True)[0])) ; assert abs(srcField.integral(True)[0]-trgFieldI.integral(True)[0])>1e-6
print("ExtensiveConservation %lf == %lf"%(srcField.getArray().accumulate()[0],trgFieldI.getArray().accumulate()[0])) ; assert abs(srcField.getArray().accumulate()[0]-trgFieldI.getArray().accumulate()[0])<1e-6

######

from numpy import *
from math import acos

med_root_dir=os.getenv("MEDCOUPLING_ROOT_DIR")
agitateur_file = ""
if med_root_dir:
  agitateur_file = os.path.join(os.getenv("MEDCOUPLING_ROOT_DIR"),"share","resources","med","agitateur.med")
if not os.path.exists(agitateur_file):
  current_dir = os.path.dirname(os.path.realpath(__file__))
  agitateur_file=os.path.join(current_dir, "..", "..", "..", "resources","agitateur.med")
pass
data=MEDFileData(agitateur_file)
ts=data.getFields()[0].getTimeSteps()
print(ts)
#
fMts=data.getFields()["DISTANCE_INTERFACE_ELEM_BODY_ELEM_DOM"]
f1ts=fMts[(2,-1)]
fMc=f1ts.getFieldAtLevel(ON_CELLS,0)
arr=fMc.getArray()
arr.getMinMaxPerComponent() # juste pour voir la plage de variation du champ par compo
ids=arr.findIdsInRange(0.,1.)
f2Mc=fMc[ids]
#
pressMts=data.getFields()["PRESSION_ELEM_DOM"]
press1ts=pressMts[(2,-1)]
pressMc=press1ts.getFieldAtLevel(ON_CELLS,0)
pressOnAgitateurMc=pressMc[ids]
#
pressOnAgitateurMc.getMesh().zipCoords()
#
agitateurMesh3DMc=pressOnAgitateurMc.getMesh()
m3DSurf,desc,descI,revDesc,revDescI=agitateurMesh3DMc.buildDescendingConnectivity()
nbOf3DCellSharing=revDescI.deltaShiftIndex()
ids2=nbOf3DCellSharing.findIdsEqual(1)
agitateurSkinMc=m3DSurf[ids2]
OffsetsOfTupleIdsInField=revDescI[ids2]
tupleIdsInField=revDesc[OffsetsOfTupleIdsInField]
pressOnSkinAgitateurMc=pressOnAgitateurMc[tupleIdsInField]
pressOnSkinAgitateurMc.setMesh(agitateurSkinMc)
#
pressSkin=pressOnSkinAgitateurMc.getArray()
pressSkin*=1e5
areaSkin=agitateurSkinMc.getMeasureField(True).getArray()
forceSkin=pressSkin*areaSkin
normalSkin=agitateurSkinMc.buildOrthogonalField().getArray()
forceVectSkin=forceSkin*normalSkin
#
singlePolyhedron=agitateurMesh3DMc.buildSpreadZonesWithPoly()
singlePolyhedron.orientCorrectlyPolyhedrons()
centerOfMass=singlePolyhedron.computeCellCenterOfMass()

barySkin=agitateurSkinMc.computeCellCenterOfMass()
posSkin=barySkin-centerOfMass

torquePerCellOnSkin=DataArrayDouble.CrossProduct(posSkin,forceVectSkin)

zeTorque=torquePerCellOnSkin.accumulate()
print("couple = %r N.m"%(zeTorque[2])) ; assert abs(zeTorque[2]-0.37)<1e-2

speedMts=data.getFields()["VITESSE_ELEM_DOM"]
speed1ts=speedMts[(2,-1)]
speedMc=speed1ts.getFieldAtLevel(ON_CELLS,0)
speedOnSkin=speedMc.getArray()[tupleIdsInField]
powerSkin=DataArrayDouble.Dot(forceVectSkin,speedOnSkin)
power=powerSkin.accumulate()[0]
print("power = %r W"%(power)) ; assert abs(power-4.22)<1e-2

x2=posSkin[:,0]*posSkin[:,0] ; x2=x2.accumulate()[0]
y2=posSkin[:,1]*posSkin[:,1] ; y2=y2.accumulate()[0]
xy=posSkin[:,0]*posSkin[:,1] ; xy=xy.accumulate()[0]
inertiaSkin=matrix([[x2,xy],[xy,y2]])
inertiaSkinValues,inertiaSkinVects=linalg.eig(inertiaSkin)
pos=max(enumerate(inertiaSkinValues),key=lambda x: x[1])[0]
vect0=inertiaSkinVects[pos].tolist()[0]
print(vect0)

def computeAngle(locAgitateur1ts):
    fMc=locAgitateur1ts.getFieldAtLevel(ON_CELLS,0)
    arr=fMc.getArray()
    ids=arr.findIdsInRange(0.,1.)
    f2Mc=fMc[ids]
    m3DSurf,desc,descI,revDesc,revDescI=f2Mc.getMesh().buildDescendingConnectivity()
    nbOf3DCellSharing=revDescI.deltaShiftIndex()
    ids2=nbOf3DCellSharing.findIdsEqual(1)
    agitateurSkinMc=m3DSurf[ids2]
    #
    singlePolyhedron=agitateurMesh3DMc.buildSpreadZonesWithPoly()
    singlePolyhedron.orientCorrectlyPolyhedrons()
    centerOfMass=singlePolyhedron.computeCellCenterOfMass()
    bary=agitateurSkinMc.computeCellCenterOfMass()
    posSkin=bary-centerOfMass
    x2=posSkin[:,0]*posSkin[:,0] ; x2=x2.accumulate()[0]
    y2=posSkin[:,1]*posSkin[:,1] ; y2=y2.accumulate()[0]
    xy=posSkin[:,0]*posSkin[:,1] ; xy=xy.accumulate()[0]
    inertiaSkin=matrix([[x2,xy],[xy,y2]])
    inertiaSkinValues,inertiaSkinVects=linalg.eig(inertiaSkin)
    pos=max(enumerate(inertiaSkinValues),key=lambda x: x[1])[0]
    vect0=inertiaSkinVects[pos].tolist()[0]
    return vect0

vects=len(ts)*[None]
for itts,locAgitateur1ts in zip(ts,data.getFields()["DISTANCE_INTERFACE_ELEM_BODY_ELEM_DOM"]):
    angle=computeAngle(locAgitateur1ts)
    vects[itts[0]]=angle
    pass

angle2=len(ts)*[0.]
for pos in range(2, len(vects)):
    norm1=sqrt(vects[pos-1][0]*vects[pos-1][0]+vects[pos-1][1]*vects[pos-1][1])
    norm2=sqrt(vects[pos][0]*vects[pos][0]+vects[pos][1]*vects[pos][1])
    crs=vects[pos-1][0]*vects[pos][0]+vects[pos-1][1]*vects[pos][1]
    crs/=norm1 ; crs/=norm2 ; crs=min(crs,1.)
    angle2[pos]=acos(crs)#/(ts[pos][2]-ts[pos-1][2])
    pass

omega=sum(angle2)/(ts[-1][2]-ts[0][2])
print(sum(angle2)) ; assert abs(sum(angle2)-1.12)<1e-2
print("Au pdt (%d,%d) a %r s le couple est de : %r N.m, power/omega=%r N.m"%(ts[2][0],ts[2][1],ts[2][2],zeTorque[2],power/omega))
assert abs(power/omega-0.37)<1e-2
