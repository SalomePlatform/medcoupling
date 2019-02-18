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


import sys
if sys.platform == "win32":
    from MEDCouplingCompat import *
else:
    from medcoupling import *
from math import pi, sqrt

# ! [PySnippetUMeshStdBuild1_1]
coords=[-0.3,-0.3,0.,   0.2,-0.3,0.,   0.7,-0.3,0.,   -0.3,0.2,0.,   0.2,0.2,0.,
        0.7,0.2,0.,    -0.3,0.7,0.,    0.2,0.7,0.,     0.7,0.7,0. ]
nodalConnPerCell=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
# ! [PySnippetUMeshStdBuild1_1]
# ! [PySnippetUMeshStdBuild1_2]
mesh=MEDCouplingUMesh("My2DMesh",2)
# ! [PySnippetUMeshStdBuild1_2]
# ! [PySnippetUMeshStdBuild1_3]
mesh.allocateCells(5)#You can put more than 5 if you want but not less.
# adding cells
mesh.insertNextCell(NORM_QUAD4,nodalConnPerCell[:4])
mesh.insertNextCell(NORM_TRI3,nodalConnPerCell[4:7])
mesh.insertNextCell(NORM_TRI3,nodalConnPerCell[7:10])
mesh.insertNextCell(NORM_QUAD4,nodalConnPerCell[10:14])
mesh.insertNextCell(NORM_QUAD4,nodalConnPerCell[14:])
# compacting
mesh.finishInsertingCells()
# ! [PySnippetUMeshStdBuild1_3]
# ! [PySnippetUMeshStdBuild1_4]
coordsArr=DataArrayDouble(coords,9,3)#here coordsArr are declared to have 3 components, mesh will deduce that its spaceDim==3.
mesh.setCoords(coordsArr)#coordsArr contains 9 tuples, that is to say mesh contains 9 nodes.
# ! [PySnippetUMeshStdBuild1_4]
# ! [PySnippetUMeshStdBuild1_5]
mesh.checkConsistencyLight()
# ! [PySnippetUMeshStdBuild1_5]

# ! [PySnippetCMeshStdBuild1_1]
XCoords=[-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22] # 9 values along X
YCoords=[0.,0.1,0.37,0.45,0.47,0.49,1.007] # 7 values along Y
arrX=DataArrayDouble(XCoords)
arrX.setInfoOnComponent(0,"X [m]")
arrY=DataArrayDouble(YCoords)
arrY.setInfoOnComponent(0,"Y [m]")
# ! [PySnippetCMeshStdBuild1_1]
# ! [PySnippetCMeshStdBuild1_2]
mesh=MEDCouplingCMesh("My2D_CMesh")
mesh.setCoords(arrX,arrY)
# ! [PySnippetCMeshStdBuild1_2]

nodalConnPerCell=list(range(4*4))
# ! [GU_MEDCoupling1SGTUMesh_0]
mesh=MEDCoupling1SGTUMesh("myQuadMesh",NORM_QUAD4)
mesh.allocateCells(3)
mesh.insertNextCell(nodalConnPerCell[:4])
mesh.insertNextCell(nodalConnPerCell[4:8])
mesh.insertNextCell(nodalConnPerCell[8:12])
# ! [GU_MEDCoupling1SGTUMesh_0]

# ! [GU_MEDCoupling1SGTUMesh_1]
polymesh=MEDCoupling1DGTUMesh("myPolyhedra",NORM_POLYHED)
polymesh.allocateCells(1)
polymesh.insertNextCell([0,1,2,3,-1,7,6,5,4,-1,0,4,5,1,-1,1,5,6,2,-1,3,2,6,7,-1,0,3,7,4])
# ! [GU_MEDCoupling1SGTUMesh_1]

#! [UG_DataArrayDouble_0]
d=DataArrayDouble([1,2,3,4,5,6],3,2)
#! [UG_DataArrayDouble_0]
#! [UG_DataArrayDouble_1]
d=DataArrayDouble([(1,2),(3,4),(5,6)])
#! [UG_DataArrayDouble_1]
#! [UG_DataArrayDouble_2]
d=DataArrayDouble([(1,2,3),(4,5,6)])
#! [UG_DataArrayDouble_2]
#! [UG_DataArrayDouble_3]
d.rearrange(2)
#! [UG_DataArrayDouble_3]
#! [UG_DataArrayDouble_4]
i=DataArrayInt([(1,2,3),(4,5,6)])
f=DataArrayFloat([(1,2,3),(4,5,6)])
#! [UG_DataArrayDouble_4]

#! [UG_MEDCouplingCurveLinearMesh_0]
m=MEDCouplingCurveLinearMesh("myCurveLinearMesh")
m.setNodeGridStructure([2,3])
#! [UG_MEDCouplingCurveLinearMesh_0]
#! [UG_MEDCouplingCurveLinearMesh_1]
coords=DataArrayDouble([0.,0., 2.,0., 0.,1., 1.9,1.1, 0.3,1.9, 2.2,2.1],6,2)
coords.setInfoOnComponents(["X [m]","Y [m]"])
m.setCoords(coords)
#! [UG_MEDCouplingCurveLinearMesh_1]
#! [UG_MEDCouplingCurveLinearMesh_2]
m.checkConsistencyLight()
#! [UG_MEDCouplingCurveLinearMesh_2]

XCoords=[-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22];  arrX=DataArrayDouble(XCoords)
YCoords=[0.,0.1,0.37,0.45,0.47,0.49,1.007];  arrY=DataArrayDouble(YCoords)
mesh=MEDCouplingCMesh("My2D_CMesh")
mesh.setCoords(arrX,arrY)
#! [UG_MEDCouplingFieldDouble_0]
fieldOnCells=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
fieldOnCells.setName("MyTensorFieldOnCellOneTime")
fieldOnCells.setMesh(mesh)
#! [UG_MEDCouplingFieldDouble_0]
#! [UG_MEDCouplingFieldDouble_1]
fieldOnCells.setTimeUnit("ms") # Time unit is ms.
fieldOnCells.setTime(4.22,2,-1) # Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
#! [UG_MEDCouplingFieldDouble_1]
#! [UG_MEDCouplingFieldDouble_2]
array=DataArrayDouble()
array.alloc(fieldOnCells.getMesh().getNumberOfCells(),2) # Implicitly fieldOnCells will be a 2 components field.
array.fillWithValue(7.)
fieldOnCells.setArray(array)
fieldOnCells.checkConsistencyLight()
# fieldOnCells is now usable
# ...
#! [UG_MEDCouplingFieldDouble_2]
#
#! [UG_MEDCouplingFieldDouble_3]
fieldOnNodes=MEDCouplingFieldDouble(ON_NODES,ONE_TIME)
fieldOnNodes.setName("MyScalarFieldOnNodeOneTime")
fieldOnNodes.setMesh(mesh)
fieldOnNodes.setTimeUnit("ms") # Time unit is ms.
fieldOnNodes.setTime(4.22,2,-1) # Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
array=DataArrayDouble()
array.alloc(fieldOnNodes.getMesh().getNumberOfNodes(),1) # Implicitly fieldOnNodes will be a 1 component field.
array.fillWithValue(7.)
fieldOnNodes.setArray(array)
fieldOnNodes.checkConsistencyLight()
# fieldOnNodes is now usable
# ...
#! [UG_MEDCouplingFieldDouble_3]
field=fieldOnCells
#! [UG_MEDCouplingFieldDouble_4]
print(mesh.getHeapMemorySizeStr())
print(field.getHeapMemorySizeStr())
print(array.getHeapMemorySizeStr())
#! [UG_MEDCouplingFieldDouble_4]

f=fieldOnCells
nbComp=3
val=5.
#! [UG_MEDCouplingFieldDouble_5]
f.applyFunc(nbComp,val)
#! [UG_MEDCouplingFieldDouble_5]
#! [UG_MEDCouplingFieldDouble_6]
f.applyFunc(1,"sqrt(X*X+Y*Y+Z*Z)")
#! [UG_MEDCouplingFieldDouble_6]
f.applyFunc(nbComp,val)
#! [UG_MEDCouplingFieldDouble_7]
f.applyFunc(4,"IVec*y+JVec*x+KVec*z+LVec*sqrt(x*x+y*y+z*z)")
#! [UG_MEDCouplingFieldDouble_7]

field=fieldOnCells
#! [UG_MEDCouplingFieldDouble_8]
points=DataArrayDouble([(0.,0.),(1,1)])
values=field.getValueOnMulti(points)
#! [UG_MEDCouplingFieldDouble_8]
#! [UG_MEDCouplingFieldDouble_9]
field.integral(True)
field.integral(0,True)
#! [UG_MEDCouplingFieldDouble_9]
field.applyFunc(6,1.)
#! [UG_MEDCouplingFieldDouble_10]
assert(field.getNumberOfComponents()==6)
diviatorfield=field.deviator()
#! [UG_MEDCouplingFieldDouble_10]


from MEDCouplingDataForTest import MEDCouplingDataForTest
mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
#! [UG_MEDCouplingGaussPointField_0]
fieldGauss=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME);
fieldGauss.setMesh(mesh);
fieldGauss.setName("MyFirstFieldOnGaussPoint");
fieldGauss.setTimeUnit("ms") # Time unit is ms.
fieldGauss.setTime(4.22,2,-1) # Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
#! [UG_MEDCouplingGaussPointField_0]
#! [UG_MEDCouplingGaussPointField_1]
tria3CooRef=[ 0.0, 0.0, 1.0 , 0.0, 0.0, 1.0 ]
tria3CooGauss=[ 0.1, 0.8, 0.2, 0.7 ]
wg3=[0.3,0.3];
fieldGauss.setGaussLocalizationOnType(NORM_TRI3,tria3CooRef,tria3CooGauss,wg3);
#
quad4CooRef=[-1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0]
quad4CooGauss=[ 0.3, 0.2, 0.2, 0.1, 0.2, 0.4, 0.15, 0.27 ]
wg4=[0.3,0.3,0.3,0.3];
fieldGauss.setGaussLocalizationOnType(NORM_QUAD4,quad4CooRef,quad4CooGauss,wg4);
#! [UG_MEDCouplingGaussPointField_1]
#! [UG_MEDCouplingGaussPointField_2]
nbTuples=mesh.getNumberOfCellsWithType(NORM_TRI3)*len(wg3)+mesh.getNumberOfCellsWithType(NORM_QUAD4)*len(wg4)
array=DataArrayDouble.New();
values=[float(i+1) for i in range(nbTuples)]
array.setValues(values,nbTuples,1);
fieldGauss.setArray(array);
fieldGauss.checkConsistencyLight();
#! [UG_MEDCouplingGaussPointField_2]

from MEDCouplingDataForTest import MEDCouplingDataForTest
m=MEDCouplingDataForTest.build2DTargetMesh_1();
Ids=list(range(1,3))
#! [UG_ExtractForMeshes_0]
part=m[Ids]
#! [UG_ExtractForMeshes_0]
#! [UG_ExtractForMeshes_1]
subNodeIds=part.computeFetchedNodeIds()
#! [UG_ExtractForMeshes_1]
#! [UG_ExtractForMeshes_2]
m.getCoords()[subNodeIds]
#! [UG_ExtractForMeshes_2]
#! [UG_ExtractForMeshes_3]
part.zipCoords()
#! [UG_ExtractForMeshes_3]
#! [UG_ExtractForMeshes_4]
o2n=part.zipCoordsTraducer()
#! [UG_ExtractForMeshes_4]

m2=MEDCouplingDataForTest.build3DExtrudedUMesh_1()[0]
#! [UG_ExtractForMeshes_5]
bn = m2.findBoundaryNodes()
#! [UG_ExtractForMeshes_5]
#! [UG_ExtractForMeshes_6]
bc = m2.getCellIdsLyingOnNodes(bn,False)
#! [UG_ExtractForMeshes_6]

#! [UG_ExtractForMeshes_7]
m2.translate([1.,2.,3.])
#! [UG_ExtractForMeshes_7]
#! [UG_ExtractForMeshes_8]
m2.getCoords()[:]+=DataArrayDouble([1.,2.,3.],1,3)
#! [UG_ExtractForMeshes_8]
import math
#! [UG_ExtractForMeshes_9]
m2.rotate([1,2,1],[0,1,0],math.pi/3)
#! [UG_ExtractForMeshes_9]
#! [UG_ExtractForMeshes_10]
MEDCouplingPointSet.Rotate3DAlg([1,2,1],[0,1,0],math.pi/3,m2.getCoords())
#! [UG_ExtractForMeshes_10]

#! [UG_ExtractForMeshes_11]
volPerCell=m2.getMeasureField(True)
#! [UG_ExtractForMeshes_11]
#! [UG_ExtractForMeshes_12]
volPerCell.getArray().accumulate()
#! [UG_ExtractForMeshes_12]
t1=-1
#! [UG_ExtractForMeshes_13]
part=volPerCell.getArray().findIdsGreaterOrEqualTo(t1)
#! [UG_ExtractForMeshes_13]
#! [UG_ExtractForMeshes_14]
m2[part]
#! [UG_ExtractForMeshes_14]
#! [UG_ExtractForMeshes_15]
centers=m2.computeCellCenterOfMass()
#! [UG_ExtractForMeshes_15]
#! [UG_ExtractForMeshes_16]
(centers*volPerCell.getArray()).accumulate()/DataArrayDouble(volPerCell.accumulate())
#! [UG_ExtractForMeshes_16]

#! [UG_ExtractForMeshes_17]
m2.scale( [1,2,4], 6. )
#! [UG_ExtractForMeshes_17]

#! [UG_ExtractForMeshes_18]
ortho_field=m.buildOrthogonalField()
#! [UG_ExtractForMeshes_18]

#! [UG_ExtractForMeshes_19]
ibc=m2.computeIsoBarycenterOfNodesPerCell()
#! [UG_ExtractForMeshes_19]

#! [UG_ExtractForMeshes_20]
# make a structured mesh 1x5
coords=DataArrayDouble(list(range(6)))
cmesh=MEDCouplingCMesh("cmesh")
cmesh.setCoords(coords,coords[:2])

# make a mesh with two zones
zmesh=cmesh.buildUnstructured()[0,1,3,4]

# get cells ids of zones
zoneArrays=zmesh.partitionBySpreadZone()
print([ ids.getValues() for ids in zoneArrays])
#! [UG_ExtractForMeshes_20]

coordsArr=DataArrayDouble(list(range(6)))
mesh2d=MEDCouplingCMesh("mesh2d")
mesh2d.setCoords(coordsArr,coordsArr[:2])
mesh2d=mesh2d.buildUnstructured()
mesh1d=MEDCouplingCMesh("mesh1d")
mesh1d.setCoords(coordsArr,coordsArr[:1])
mesh1d=mesh1d.buildUnstructured()
mesh1d.rotate( [2.3,0], math.radians( 25 ))
mesh1d.translate( [0.2,0.4] )
#! [UG_ExtractForMeshes_21]
m2d,m1d,a2d,a1d=MEDCouplingUMesh.Intersect2DMeshWith1DLine( mesh2d, mesh1d, 1e-12 )
#! [UG_ExtractForMeshes_21]
#print ("a2d",a2d.getValues())
#print (a1d.getValues())

mesh1=mesh2d.deepCopy()
mesh1.rotate( [2.3,0], math.radians( 25 ))
mesh2=mesh2d
#! [UG_ExtractForMeshes_22]
m,a1,a2=MEDCouplingUMesh.Intersect2DMeshes( mesh1, mesh2, 1e-12 )
#! [UG_ExtractForMeshes_22]


points=DataArrayDouble([(0,0,0)])
mesh=m2.computeSkin()
#! [UG_ExtractForMeshes_23]
dist,cells=mesh.distanceToPoints(points)
#! [UG_ExtractForMeshes_23]

arr=DataArrayDouble([1,2,3,4,5,6],3,2)
a,b=2,5
#! [UG_ExtractForArrays_0]
tupleIds = arr[:,0].findIdsInRange(a,b)
#! [UG_ExtractForArrays_0]
c,d=0,7
#! [UG_ExtractForArrays_1]
tupleIds1 = arr.magnitude().findIdsInRange(c,d)
#! [UG_ExtractForArrays_1]
#! [UG_ExtractForArrays_2]
tupleIds2 = DataArrayInt.buildSubstraction(tupleIds,tupleIds1)
#! [UG_ExtractForArrays_2]

valsArr1=DataArrayDouble(list(range(9*2)),9,2)
field4 = MEDCouplingFieldDouble(ON_NODES)
field4.setArray(valsArr1)
mesh=MEDCouplingCMesh("My2D_CMesh")
coo=DataArrayDouble([0,1,2])
mesh.setCoords(coo,coo)
field4.setMesh(mesh)
ids4=[1,2]
#! [UG_ExtractForFields_0]
subField = field4[ids4]
#! [UG_ExtractForFields_0]

m4=MEDCouplingCMesh("box")
coo=DataArrayDouble(list(range(7)))
m4.setCoords(coo[:5],coo[:5],coo)
m4=m4.buildUnstructured()
valsArr1=m4.computeCellCenterOfMass()
valsArr1.applyFunc(1,"sqrt(X*X+Y*Y+Z*Z)")
field5 = MEDCouplingFieldDouble(ON_CELLS)
field5.setArray(valsArr1)
field5.setMesh(m4)
#! [UG_ExtractForFields_1]
origin=[0,0,2]
normvec=[-1,-1,6]
slice5=field5.extractSlice3D(origin,normvec,1e-10)
#! [UG_ExtractForFields_1]

from MEDCouplingDataForTest import MEDCouplingDataForTest
m1=MEDCouplingDataForTest.build2DTargetMesh_1();
m2=m1
eps=1e-12
#! [UG_MeshComparison_0]
m1.isEqual(m2,eps)
#! [UG_MeshComparison_0]
#! [UG_MeshComparison_1]
m1.isEqualWithoutConsideringStr(m2,eps)
#! [UG_MeshComparison_1]

arr=DataArrayDouble(5) ; arr.iota()
m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
m=m.buildUnstructured()
m1=m.deepCopy()
arr=DataArrayInt([0,3,1,4,2])
m2=m1[arr]
#! [UG_MeshComparison_2]
a,_=m1.checkGeoEquivalWith(m2,20,1e-12)
assert(m1[a].isEqualWithoutConsideringStr(m2,1e-12))
#! [UG_MeshComparison_2]
m3=m1
m3.translate([100.,0])
m1=MEDCouplingUMesh.MergeUMeshes([m3,m])
m2=MEDCouplingUMesh.MergeUMeshes([m,m3])
#! [UG_MeshComparison_3]
a,b=m1.checkGeoEquivalWith(m2,12,1e-12)
m2.renumberNodes(b,len(b))
assert(m1[a].isEqualWithoutConsideringStr(m2,1e-12))
#! [UG_MeshComparison_3]


coords = [0.,2.,4.]
coordsArr=DataArrayDouble(coords,3,1)
m4=MEDCouplingCMesh()
m4.setCoords(coordsArr,coordsArr,coordsArr)
m4=m4.buildUnstructured()
field1 = m4.fillFromAnalytic(ON_NODES,1,"x+y")
pts = [0.5,0.5,0.5,1.,1,1.]
#! [UG_CommonHandlingMesh_0]
assert(field1.getTypeOfField()==ON_NODES)
field1.getMesh().simplexize(PLANAR_FACE_5)
field1.getValueOnMulti(pts)
#! [UG_CommonHandlingMesh_0]

from MEDCouplingDataForTest import MEDCouplingDataForTest
m2=MEDCouplingDataForTest.build2DTargetMesh_1();
m2.changeSpaceDimension(3);
m1=MEDCouplingDataForTest.buildCU1DMesh_U();
m1.changeSpaceDimension(3);
center=[0.,0.,0.]
vector=[0.,1.,0.]
m1.rotate(center,vector,-pi/2.);
#! [UG_CommonHandlingMesh_1]
m3=m2.buildExtrudedMesh(m1,0);
#! [UG_CommonHandlingMesh_1]

#! [UG_CommonHandlingMesh_2]
m5=MEDCouplingUMesh.MergeUMeshes([m3,m4])
#! [UG_CommonHandlingMesh_2]

#! [UG_CommonHandlingMesh_3]
m1.convertLinearCellsToQuadratic(0)
#! [UG_CommonHandlingMesh_3]

#! [UG_CommonHandlingMesh_4]
skin = m1.computeSkin()
#! [UG_CommonHandlingMesh_4]
#! [UG_CommonHandlingMesh_5]
#bc = m1.getCellIdsLyingOnNodes(bn,False)
#! [UG_CommonHandlingMesh_5]

mesh3d=m3
#! [UG_CommonHandlingMesh_6]
mesh1d,d,di,r,ri=mesh3d.explodeIntoEdges()
#! [UG_CommonHandlingMesh_6]
#! [UG_CommonHandlingMesh_7]
mesh2d,d,di,r,ri=mesh3d.buildDescendingConnectivity()
#! [UG_CommonHandlingMesh_7]

mesh2d=MEDCouplingCMesh()
mesh2d.setCoords(coordsArr,coordsArr)
mesh2d=mesh2d.buildUnstructured()
eps=1e-12
#! [UG_CommonHandlingMesh_8]
changedCells=mesh2d.conformize2D(eps)
#! [UG_CommonHandlingMesh_8]

#! [UG_CommonHandlingMesh_9]
mesh2d.duplicateNodes([3,4])
#! [UG_CommonHandlingMesh_9]

m1d=MEDCouplingCMesh()
m1d.setCoords(coordsArr)
m1d=m1d.buildUnstructured()
#! [UG_CommonHandlingMesh_10]
m1d.renumberCells(m1d.orderConsecutiveCells1D().invertArrayN2O2O2N(m1d.getNumberOfCells()))
#! [UG_CommonHandlingMesh_10]

skin=m3.computeSkin()
#! [UG_CommonHandlingMesh_11]
vec=[0,0,-1]
skin.orientCorrectly2DCells(vec,False)
#! [UG_CommonHandlingMesh_11]

#! [UG_CommonHandlingMesh_12]
m3.orientCorrectlyPolyhedrons()
#! [UG_CommonHandlingMesh_12]

#! [UG_CommonHandlingMesh_13]
m3.renumberCells(m3.rearrange2ConsecutiveCellTypes())
m3.sortCellsInMEDFileFrmt()
#! [UG_CommonHandlingMesh_13]

m=MEDCouplingCMesh()
m.setCoords(coordsArr[:2],coordsArr[:2])
m=m.buildUnstructured()
#! [UG_CommonHandlingMesh_14]
m.renumberNodes([2,1,0,-1],3)
#! [UG_CommonHandlingMesh_14]

#! [UG_CommonHandlingMesh_15]
mtet,n2ocells,np=m3.tetrahedrize(PLANAR_FACE_5)
#! [UG_CommonHandlingMesh_15]
#! [UG_CommonHandlingMesh_16]
m.mergeNodes(1e-12)
#! [UG_CommonHandlingMesh_16]


#! [UG_Projection_0]
srcCoo=DataArrayDouble([(0,0),(1,0),(3,0),(0,1),(1,1),(3,1)])
src=MEDCouplingUMesh("src",2)
src.setCoords(srcCoo)
src.allocateCells()
src.insertNextCell(NORM_QUAD4,[0,3,4,1])
src.insertNextCell(NORM_QUAD4,[1,4,5,2])
#
trgCoo=DataArrayDouble([(0.5,0.5),(1.5,0.5),(1.5,1.5)])
trg=MEDCouplingUMesh("trg",2)
trg.setCoords(trgCoo)
trg.allocateCells()
trg.insertNextCell(NORM_TRI3,[0,2,1])
#! [UG_Projection_0]
from MEDCouplingRemapper import MEDCouplingRemapper
#! [UG_Projection_1]
rem=MEDCouplingRemapper()
rem.prepare(src,trg,"P0P0")
print(rem.getCrudeMatrix())
#! [UG_Projection_1]
#! [UG_Projection_2]
srcF=MEDCouplingFieldDouble(ON_CELLS)
srcF.setMesh(src)
srcF.setArray(DataArrayDouble([3,4]))
srcF.setNature(IntensiveMaximum)
#
trgF=rem.transferField(srcF,-1)
#! [UG_Projection_2]
#! [UG_Projection_3]
rem=MEDCouplingRemapper()
rem.prepare(src,trg,"P0P1")
print(rem.getCrudeMatrix())
#! [UG_Projection_3]
#! [UG_Projection_4]
rem=MEDCouplingRemapper()
rem.prepare(src,trg,"P1P0")
print(rem.getCrudeMatrix())
#! [UG_Projection_4]


#! [UG_Projection_5]
coo=DataArrayDouble([(0.,0.,0.), (1,0,0), (0,1,0)])
src=MEDCouplingUMesh("src",2)
src.setCoords(coo)
src.allocateCells(1)
src.insertNextCell(NORM_TRI3,[0,1,2])
tgt = src.deepCopy()
rem=MEDCouplingRemapper()
rem.prepare(src,tgt,"P0P0")
print(rem.getCrudeMatrix())
#! [UG_Projection_5]
#! [UG_Projection_6]
src.translate([0,0,1e-3])
rem.prepare(src,tgt,"P0P0")
print(rem.getCrudeMatrix())
#! [UG_Projection_6]
#! [UG_Projection_7]
rem.setBoundingBoxAdjustmentAbs( 1e-3 )
rem.prepare(src,tgt,"P0P0")
print(rem.getCrudeMatrix())
#! [UG_Projection_7]

import math
#! [UG_Projection_8]
src.rotate([0,0,0],[0,1,0],math.pi/4)
rem.prepare(src,tgt,"P0P0")
print(rem.getCrudeMatrix())
#! [UG_Projection_8]
#! [UG_Projection_9]
rem.setMaxDistance3DSurfIntersect( 0.1 )
rem.prepare(src,tgt,"P0P0")
print(rem.getCrudeMatrix())
#! [UG_Projection_9]

#! [UG_Projection_10]
rem.setMaxDistance3DSurfIntersect( -1 ) # switch it off
rem.setMinDotBtwPlane3DSurfIntersect( 0.8 )
rem.prepare(src,tgt,"P0P0")
print(rem.getCrudeMatrix())
#! [UG_Projection_10]
