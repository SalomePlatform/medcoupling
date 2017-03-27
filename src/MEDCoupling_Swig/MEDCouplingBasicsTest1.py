#  -*- coding: utf-8 -*-
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

from MEDCoupling import *
import unittest
from math import pi,e,sqrt,cos,sin
from datetime import datetime
from MEDCouplingDataForTest import MEDCouplingDataForTest
import rlcompleter,readline # this line has to be here, to ensure a usability of MEDCoupling/MEDLoader. B4 removing it please notify to anthony.geay@edf.fr

class MEDCouplingBasicsTest1(unittest.TestCase):
    def testArray2(self):
        arr=DataArrayDouble.New()
        arr.setValues([12.,11.,10.,9.,8.,7.,6.,5.,4.,3.,2.,1.],3,4)
        arr.setInfoOnComponent(0,"ggg");
        arr.setInfoOnComponent(1,"hhhh");
        arr.setInfoOnComponent(2,"jj");
        arr.setInfoOnComponent(3,"kkkkkk");
        arr2=arr.convertToIntArr();
        arr3=arr2.convertToDblArr();
        self.assertTrue(arr.isEqual(arr3,1e-14))
        pass

    def testArray3(self):
        arr1=DataArrayInt.New();
        arr1Ref=[0,10,1,11,2,12,3,13,4,14,5,15,6,16]
        arr1.setValues(arr1Ref,7,2);
        self.assertEqual(7,arr1.getNumberOfTuples());
        self.assertEqual(2,arr1.getNumberOfComponents());
        self.assertEqual(arr1Ref,list(arr1.getValues()));
        arr2=arr1.subArray(3);
        self.assertEqual(4,arr2.getNumberOfTuples());
        self.assertEqual(2,arr2.getNumberOfComponents());
        self.assertEqual(arr1Ref[6:],list(arr2.getValues()));
        arr3=arr1.subArray(2,5);
        self.assertEqual(3,arr3.getNumberOfTuples());
        self.assertEqual(2,arr3.getNumberOfComponents());
        self.assertEqual(arr1Ref[4:10],list(arr3.getValues()));
        #
        arr4=DataArrayDouble.New();
        arr4Ref=[0.8,10.8,1.9,11.9,2.1,12.1,3.2,13.2,4.3,14.3,5.4,15.4,6.5,16.5]
        arr4.setValues(arr4Ref,7,2);
        self.assertEqual(7,arr4.getNumberOfTuples());
        self.assertEqual(2,arr4.getNumberOfComponents());
        tmp=arr4.getValues()
        for i in range(14):
            self.assertTrue(abs(arr4Ref[i]-tmp[i])<1e-14);
            pass
        arr5=arr4.subArray(3);
        self.assertEqual(4,arr5.getNumberOfTuples());
        self.assertEqual(2,arr5.getNumberOfComponents());
        tmp=arr5.getValues()
        for i in range(8):
            self.assertTrue(abs(arr4Ref[6+i]-tmp[i])<1e-14);
            pass
        arr6=arr4.subArray(2,5);
        self.assertEqual(3,arr6.getNumberOfTuples());
        self.assertEqual(2,arr6.getNumberOfComponents());
        tmp=arr6.getValues()
        for i in range(6):
            self.assertTrue(abs(arr4Ref[4+i]-tmp[i])<1e-14);
            pass
        pass

    def testMesh(self):
        tab4=[1, 2, 8, 7, 2, 3, 9, 8, 3,
              4, 10, 9, 4, 5, 11, 10, 5,
              0, 6, 11, 0, 1, 7, 6 ]
        nbOfNodes=12
        nbOfCells=6
        coords=[ 0.024155, 0.04183768725682622, -0.305, 0.04831000000000001, -1.015761910347357e-17,
                 -0.305, 0.09662000000000001, -1.832979297858306e-18, -0.305, 0.120775, 0.04183768725682623,
                 -0.305, 0.09662000000000001, 0.08367537451365245, -0.305, 0.04831000000000001, 0.08367537451365246,
                 -0.305, 0.024155, 0.04183768725682622, -0.2863, 0.04831000000000001, -1.015761910347357e-17, -0.2863, 
                 0.09662000000000001, -1.832979297858306e-18, -0.2863, 0.120775, 0.04183768725682623, -0.2863, 0.09662000000000001,
                 0.08367537451365245, -0.2863, 0.04831000000000001, 0.08367537451365246, -0.2863 ]
        self.assertEqual(MEDCouplingMesh.GetNumberOfNodesOfGeometricType(NORM_TRI3),3)
        self.assertTrue(MEDCouplingMesh.IsStaticGeometricType(NORM_TRI3))
        self.assertTrue(MEDCouplingMesh.IsLinearGeometricType(NORM_TRI3))
        self.assertEqual(MEDCouplingMesh.GetDimensionOfGeometricType(NORM_TRI3),2)
        self.assertEqual(MEDCouplingMesh.GetReprOfGeometricType(NORM_TRI3),"NORM_TRI3")
        self.assertRaises(InterpKernelException,MEDCouplingMesh.GetNumberOfNodesOfGeometricType,NORM_POLYGON)
        self.assertTrue(not MEDCouplingMesh.IsStaticGeometricType(NORM_POLYGON))
        self.assertTrue(MEDCouplingMesh.IsLinearGeometricType(NORM_POLYGON))
        self.assertEqual(MEDCouplingMesh.GetDimensionOfGeometricType(NORM_POLYGON),2)
        self.assertEqual(MEDCouplingMesh.GetReprOfGeometricType(NORM_POLYGON),"NORM_POLYGON")
        self.assertEqual(MEDCouplingMesh.GetNumberOfNodesOfGeometricType(NORM_TRI6),6)
        self.assertTrue(MEDCouplingMesh.IsStaticGeometricType(NORM_TRI6))
        self.assertTrue(not MEDCouplingMesh.IsLinearGeometricType(NORM_TRI6))
        self.assertEqual(MEDCouplingMesh.GetDimensionOfGeometricType(NORM_TRI6),2)
        self.assertEqual(MEDCouplingMesh.GetReprOfGeometricType(NORM_TRI6),"NORM_TRI6")
        mesh=MEDCouplingUMesh.New()
        mesh.setMeshDimension(2)
        mesh.allocateCells(8);
        mesh.setName("mesh1")
        self.assertTrue(mesh.getName()=="mesh1")
        for i in range(nbOfCells):
            mesh.insertNextCell(NORM_QUAD4,4,tab4[4*i:4*(i+1)]);
            pass
        mesh.finishInsertingCells()
        self.assertTrue(mesh.getNumberOfCells()==nbOfCells)
        self.assertTrue(mesh.getNodalConnectivity().getNbOfElems()==30)
        self.assertTrue(mesh.getNodalConnectivityIndex().getNbOfElems()==nbOfCells+1)
        myCoords=DataArrayDouble.New()
        myCoords.setValues(coords,nbOfNodes,3);
        self.assertTrue(myCoords.getIJ(3,2)==-0.305)
        mesh.setCoords(myCoords);
        mesh.checkConsistencyLight();
        self.assertTrue(mesh.getAllGeoTypes()==[4])
        myFalseConn=DataArrayInt.New()
        myFalseConn.setValues(tab4,6,4)
        self.assertTrue(myFalseConn.getIJ(1,1)==3)
        #
        field=MEDCouplingFieldDouble.New(ON_CELLS)
        field.setMesh(mesh)
        field.setNature(ExtensiveMaximum)
        myCoords=DataArrayDouble.New()
        sampleTab=[]
        for i in range(nbOfCells * 9):
            sampleTab.append(float(i))
        myCoords.setValues(sampleTab,nbOfCells,9);
        field.setArray(myCoords)
        self.assertTrue(3==mesh.getSpaceDimension())
        field.checkConsistencyLight()
        mesh2=mesh.clone(False)
        mesh3=mesh.clone(True)
        mesh3=0
        mesh2=0
        ## deep full recursively copy of field -> both field and mesh underneath copied
        field2=field.clone(True)
        field2.setMesh(field.getMesh().clone(True))
        mesh3=mesh.clone(True)
        field3=mesh3.fillFromAnalytic(ON_CELLS,2,"x*IVec+(y+z)*JVec")
        field3.applyFunc("u*u*u+cos(u)")
        pass
        
    def testMeshPointsCloud(self):
        targetCoords=[-0.3,-0.3,0.5, 0.2,-0.3,1., 0.7,-0.3,1.5,
                      -0.3,0.2,0.5, 0.2,0.2,1., 0.7,0.2,1.5, -0.3,0.7,0.5, 0.2,0.7,1., 0.7,0.7,1.5]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(0);
        targetMesh.allocateCells(8);
        targetMesh.insertNextCell(NORM_POINT1,1,[0]);
        targetMesh.insertNextCell(NORM_POINT1,1,[1]);
        targetMesh.insertNextCell(NORM_POINT1,1,[2]);
        targetMesh.insertNextCell(NORM_POINT1,1,[3]);
        targetMesh.insertNextCell(NORM_POINT1,1,[4]);
        targetMesh.insertNextCell(NORM_POINT1,1,[5]);
        targetMesh.insertNextCell(NORM_POINT1,1,[7]);
        targetMesh.insertNextCell(NORM_POINT1,1,[6]);
        targetMesh.finishInsertingCells();
        self.assertRaises(InterpKernelException,targetMesh.checkConsistencyLight);
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,9,3);
        targetMesh.setCoords(myCoords);
        self.assertEqual(targetMesh.getSpaceDimension(),3)
        self.assertEqual(targetMesh.getNumberOfCells(),8)
        self.assertEqual(targetMesh.getNumberOfNodes(),9)
        self.assertEqual(targetMesh.getMeshDimension(),0)
        pass

    def testMeshM1D(self):
        meshM1D=MEDCouplingUMesh.New();
        self.assertRaises(InterpKernelException,meshM1D.getMeshDimension);
        self.assertRaises(InterpKernelException,meshM1D.getNumberOfNodes);
        self.assertRaises(InterpKernelException,meshM1D.getNumberOfCells);
        self.assertRaises(InterpKernelException,meshM1D.setMeshDimension,-2)
        self.assertRaises(InterpKernelException,meshM1D.setMeshDimension,-10)
        meshM1D.setMeshDimension(-1);
        meshM1D.checkConsistencyLight();
        self.assertEqual(meshM1D.getMeshDimension(),-1);
        self.assertEqual(meshM1D.getNumberOfCells(),1);
        self.assertRaises(InterpKernelException,meshM1D.getNumberOfNodes);
        self.assertRaises(InterpKernelException,meshM1D.getSpaceDimension);
        cpy=meshM1D.clone(True);
        self.assertTrue(cpy.isEqual(meshM1D,1e-12));
        fieldOnCells=MEDCouplingFieldDouble.New(ON_CELLS);
        fieldOnCells.setMesh(meshM1D);
        array=DataArrayDouble.New();
        array.setValues(6*[7.],1,6);
        fieldOnCells.setArray(array);
        fieldOnCells.checkConsistencyLight();
        pass
    
    def testDeepCopy(self):
        array=DataArrayDouble.New();
        array.setValues(5*3*[7.],5,3);
        self.assertEqual(array.getIJ(3,2),7.);
        array2=array.deepCopy();
        self.assertEqual(array2.getIJ(3,2),7.)
        #
        array3=DataArrayInt.New();
        array3.setValues(5*3*[17],5,3);
        self.assertEqual(array3.getIJ(3,2),17);
        array4=array3.deepCopy();
        self.assertEqual(array4.getIJ(3,2),17);
        pass
    
    def testRevNodal(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1()
        revNodal,revNodalIndx=mesh.getReverseNodalConnectivity();
        revNodalExpected=[0,0,1,1,2,0,3,0,1,2,3,4,2,4,3,3,4,4];
        revNodalIndexExpected=[0,1,3,5,7,12,14,15,17,18];
        self.assertEqual(revNodal.getNbOfElems(),18)
        self.assertEqual(revNodalIndx.getNbOfElems(),10)
        self.assertEqual(list(revNodal.getValues()),revNodalExpected)
        self.assertEqual(list(revNodalIndx.getValues()),revNodalIndexExpected)
        pass
    
    def testConvertToPolyTypes(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        elts=[1,3];
        mesh.convertToPolyTypes(elts);
        mesh.checkConsistencyLight();
        self.assertEqual(5,mesh.getNumberOfCells());
        self.assertEqual(23,mesh.getNodalConnectivity().getNumberOfTuples());
        expected1=[4, 0, 3, 4, 1, 5, 1, 4, 2, 3, 4, 5, 2, 5, 6, 7, 4, 3, 4, 7, 8, 5, 4]
        self.assertEqual(expected1,list(mesh.getNodalConnectivity().getValues()));
        #
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        mesh.convertToPolyTypes(elts);
        mesh.checkConsistencyLight();
        self.assertEqual(8,mesh.getNumberOfCells());
        self.assertEqual(114,mesh.getNodalConnectivity().getNumberOfTuples());
        mesh.convertToPolyTypes(elts);
        mesh.checkConsistencyLight();
        self.assertEqual(8,mesh.getNumberOfCells());
        self.assertEqual(114,mesh.getNodalConnectivity().getNumberOfTuples());
        pass

    def testDescConn2D(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        desc=DataArrayInt.New();
        descIndx=DataArrayInt.New();
        revDesc=DataArrayInt.New();
        revDescIndx=DataArrayInt.New();
        mesh2=mesh.buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
        mesh2.checkConsistencyLight();
        self.assertEqual(1,mesh2.getMeshDimension());
        self.assertEqual(13,mesh2.getNumberOfCells());
        self.assertEqual(14,revDescIndx.getNbOfElems()); self.assertEqual(14,revDescIndx.getNumberOfTuples());
        self.assertEqual(6,descIndx.getNbOfElems()); self.assertEqual(6,descIndx.getNumberOfTuples());
        self.assertEqual(18,desc.getNbOfElems()); self.assertEqual(18,desc.getNumberOfTuples());
        self.assertEqual(18,revDesc.getNbOfElems()); self.assertEqual(18,revDesc.getNumberOfTuples());
        expected1=[0,1,2,3, 2,4,5, 6,7,4, 8,9,1,10, 11,12,6,9];
        self.assertEqual(expected1,list(desc.getValues()));
        expected2=[0,4,7,10,14,18];
        self.assertEqual(expected2,list(descIndx.getValues()));
        expected3=[0,1,3,5,6,8,9,11,12,13,15,16,17,18];
        self.assertEqual(expected3,list(revDescIndx.getValues()));
        expected4=[0, 0,3, 0,1, 0, 1,2, 1, 2,4, 2, 3, 3,4, 3, 4, 4];
        self.assertEqual(expected4,list(revDesc.getValues()));
        conn=mesh2.getNodalConnectivity();
        connIndex=mesh2.getNodalConnectivityIndex();
        expected5=[0,3,6,9,12,15,18,21,24,27,30,33,36,39];
        self.assertEqual(expected5,list(connIndex.getValues()));
        expected6=[1, 0, 3, 1, 3, 4, 1, 4, 1, 1, 1, 0, 1, 4, 2, 1, 2, 1, 1, 4, 5, 1, 5, 2, 1, 6, 7, 1, 7, 4, 1, 3, 6, 1, 7, 8, 1, 8, 5];
        self.assertEqual(expected6,list(conn.getValues()));
        #
        eltsV=[1,3];
        mesh.convertToPolyTypes(eltsV);
        mesh.checkConsistencyLight();
        #
        desc=DataArrayInt.New();
        descIndx=DataArrayInt.New();
        revDesc=DataArrayInt.New();
        revDescIndx=DataArrayInt.New();
        #
        mesh2=mesh.buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
        mesh2.checkConsistencyLight();
        self.assertEqual(1,mesh2.getMeshDimension());
        self.assertEqual(13,mesh2.getNumberOfCells());
        self.assertEqual(14,revDescIndx.getNbOfElems()); self.assertEqual(14,revDescIndx.getNumberOfTuples());
        self.assertEqual(6,descIndx.getNbOfElems()); self.assertEqual(6,descIndx.getNumberOfTuples());
        self.assertEqual(18,desc.getNbOfElems()); self.assertEqual(18,desc.getNumberOfTuples());
        self.assertEqual(18,revDesc.getNbOfElems()); self.assertEqual(18,revDesc.getNumberOfTuples());
        self.assertEqual(expected1,list(desc.getValues()));
        self.assertEqual(expected2,list(descIndx.getValues()));
        self.assertEqual(expected3,list(revDescIndx.getValues()));
        self.assertEqual(expected4,list(revDesc.getValues()));
        conn=mesh2.getNodalConnectivity();
        connIndex=mesh2.getNodalConnectivityIndex();
        self.assertEqual(expected5,list(connIndex.getValues()));
        self.assertEqual(expected6,list(conn.getValues()));
        pass
    
    def testDescConn3D(self):
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        desc=DataArrayInt.New();
        descIndx=DataArrayInt.New();
        revDesc=DataArrayInt.New();
        revDescIndx=DataArrayInt.New();
        #
        mesh2=mesh.buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
        mesh2.checkConsistencyLight();
        self.assertEqual(2,mesh2.getMeshDimension());
        self.assertEqual(36,mesh2.getNumberOfCells());
        self.assertEqual(37,revDescIndx.getNbOfElems()); self.assertEqual(37,revDescIndx.getNumberOfTuples());
        self.assertEqual(9,descIndx.getNbOfElems()); self.assertEqual(9,descIndx.getNumberOfTuples());
        self.assertEqual(48,desc.getNbOfElems()); self.assertEqual(48,desc.getNumberOfTuples());
        self.assertEqual(48,revDesc.getNbOfElems()); self.assertEqual(48,revDesc.getNumberOfTuples());
        expected1=[0, 6, 12, 18, 24, 30, 36, 42, 48]
        expected2=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 3, 11, 12, 4, 13, 14, 15, 16, 17, 10, 18, 19, 13, 1, 20, 21, 22, 23, 24, 7, 25, 26, 27, 28, 22, 12, 29, 23, 30, 31, 32, 17, 33, 28, 34, 35, 30]
        expected3=[0, 1, 3, 4, 6, 8, 9, 10, 12, 13, 14, 16, 17, 19, 21, 22, 23, 24, 26, 27, 28, 29, 30, 32, 34, 35, 36, 37, 38, 40, 41, 43, 44, 45, 46, 47, 48]
        expected4=[0, 0, 4, 0, 0, 1, 0, 2, 0, 1, 1, 5, 1, 1, 1, 3, 2, 2, 6, 2, 3, 2, 2, 3, 3, 7, 3, 3, 4, 4, 4, 5, 4, 6, 4, 5, 5, 5, 5, 7, 6, 6, 7, 6, 6, 7, 7, 7]
        expected5=[0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180]
        expected6=[4, 0, 1, 4, 3, 4, 9, 12, 13, 10, 4, 0, 9, 10, 1, 4, 1, 10, 13, 4, 4, 4, 13, 12, 3, 4, 3, 12, 9, 0, 4, 1, 2, 5, 4, 4, 10, 13, 14, 11, 4, 1, 10, 11, 2, 4, 2, 11, 14,
                   5, 4, 5, 14, 13, 4, 4, 3, 4, 7, 6, 4, 12, 15, 16, 13, 4, 4, 13, 16, 7, 4, 7, 16, 15, 6, 4, 6, 15, 12, 3, 4, 4, 5, 8, 7, 4, 13, 16, 17, 14, 4, 5, 14, 17, 8, 4, 8,
                   17, 16, 7, 4, 18, 21, 22, 19, 4, 9, 18, 19, 10, 4, 10, 19, 22, 13, 4, 13, 22, 21, 12, 4, 12, 21, 18, 9, 4, 19, 22, 23, 20, 4, 10, 19, 20, 11, 4, 11, 20, 23, 14, 4,
                   14, 23, 22, 13, 4, 21, 24, 25, 22, 4, 13, 22, 25, 16, 4, 16, 25, 24, 15, 4, 15, 24, 21, 12, 4, 22, 25, 26, 23, 4, 14, 23, 26, 17, 4, 17, 26, 25, 16]
        expected7=[4, 0, 1, 4, 3, 4, 9, 12, 13, 10, 4, 0, 9, 10, 1, 4, 1, 10, 13, 4, 4, 4, 13, 12, 3, 4, 3, 12, 9, 0, 5, 1, 2, 5, 4, 5, 10, 13, 14, 11, 5, 1, 10, 11, 2, 5, 2, 11, 14,
                   5, 5, 5, 14, 13, 4, 4, 3, 4, 7, 6, 4, 12, 15, 16, 13, 4, 4, 13, 16, 7, 4, 7, 16, 15, 6, 4, 6, 15, 12, 3, 5, 4, 5, 8, 7, 5, 13, 16, 17, 14, 5, 5, 14, 17, 8, 5, 8,
                   17, 16, 7, 4, 18, 21, 22, 19, 4, 9, 18, 19, 10, 4, 10, 19, 22, 13, 4, 13, 22, 21, 12, 4, 12, 21, 18, 9, 4, 19, 22, 23, 20, 4, 10, 19, 20, 11, 4, 11, 20, 23, 14, 4,
                   14, 23, 22, 13, 4, 21, 24, 25, 22, 4, 13, 22, 25, 16, 4, 16, 25, 24, 15, 4, 15, 24, 21, 12, 4, 22, 25, 26, 23, 4, 14, 23, 26, 17, 4, 17, 26, 25, 16]
        
        self.assertEqual(expected1,list(descIndx.getValues()));
        self.assertEqual(expected2,list(desc.getValues()));
        self.assertEqual(expected3,list(revDescIndx.getValues()));
        self.assertEqual(expected4,list(revDesc.getValues()));
        self.assertEqual(expected5,list(mesh2.getNodalConnectivityIndex().getValues()));
        self.assertEqual(expected6,list(mesh2.getNodalConnectivity().getValues()));
        #
        eltsV=[1,3]
        mesh.convertToPolyTypes(eltsV);
        mesh.checkConsistencyLight();
        desc=DataArrayInt.New();
        descIndx=DataArrayInt.New();
        revDesc=DataArrayInt.New();
        revDescIndx=DataArrayInt.New();
        mesh2=mesh.buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
        mesh2.checkConsistencyLight();
        self.assertEqual(2,mesh2.getMeshDimension());
        self.assertEqual(36,mesh2.getNumberOfCells());
        self.assertEqual(37,revDescIndx.getNbOfElems()); self.assertEqual(37,revDescIndx.getNumberOfTuples());
        self.assertEqual(9,descIndx.getNbOfElems()); self.assertEqual(9,descIndx.getNumberOfTuples());
        self.assertEqual(48,desc.getNbOfElems()); self.assertEqual(48,desc.getNumberOfTuples());
        self.assertEqual(48,revDesc.getNbOfElems()); self.assertEqual(48,revDesc.getNumberOfTuples());
        self.assertEqual(expected1,list(descIndx.getValues()));
        self.assertEqual(expected2,list(desc.getValues()));
        self.assertEqual(expected3,list(revDescIndx.getValues()));
        self.assertEqual(expected4,list(revDesc.getValues()));
        self.assertEqual(expected5,list(mesh2.getNodalConnectivityIndex().getValues()));
        self.assertEqual(expected7,list(mesh2.getNodalConnectivity().getValues()));
        pass

    def testFindBoundaryNodes(self):
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        boundaryNodes=mesh.findBoundaryNodes();
        expected1=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26];
        self.assertEqual(expected1,boundaryNodes.getValues());
        pass

    def testBoundaryMesh(self):
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        mesh2=mesh.buildBoundaryMesh(False);
        self.assertEqual(24,mesh2.getNumberOfCells());
        self.assertEqual(26,mesh2.getNumberOfNodes());
        pass

    def testBuildPartOfMySelf(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        mesh.setName("Toto");
        tab1=[0,4]
        tab2=[0,2,3]
        #
        subMesh=mesh.buildPart(tab1)
        self.assertTrue(isinstance(subMesh,MEDCouplingUMesh))
        subMesh=mesh.buildPartOfMySelf(tab1,True);
        self.assertTrue(isinstance(subMesh,MEDCouplingUMesh))
        name=subMesh.getName();
        self.assertEqual(2,len(mesh.getAllGeoTypes()));
        self.assertEqual(NORM_TRI3,mesh.getAllGeoTypes()[0]);
        self.assertEqual(NORM_QUAD4,mesh.getAllGeoTypes()[1]);
        self.assertEqual(1,len(subMesh.getAllGeoTypes()));
        self.assertEqual(NORM_QUAD4,subMesh.getAllGeoTypes()[0]);
        self.assertEqual(name,"Toto");
        self.assertEqual(2,subMesh.getNumberOfCells());
        subConn=[4,0,3,4,1,4,7,8,5,4];
        subConnIndex=[0,5,10];
        self.assertEqual(10,subMesh.getNodalConnectivity().getNbOfElems());
        self.assertEqual(3,subMesh.getNodalConnectivityIndex().getNbOfElems());
        self.assertEqual(subConn[0:10],list(subMesh.getNodalConnectivity().getValues()));
        self.assertEqual(subConnIndex[0:3],list(subMesh.getNodalConnectivityIndex().getValues()));
        #
        subMesh=mesh.buildPartOfMySelf(tab2[0:3],True);
        self.assertTrue(isinstance(subMesh,MEDCouplingUMesh));
        name=subMesh.getName();
        self.assertEqual(2,len(subMesh.getAllGeoTypes()));
        self.assertEqual(NORM_TRI3,subMesh.getAllGeoTypes()[0]);
        self.assertEqual(NORM_QUAD4,subMesh.getAllGeoTypes()[1]);
        self.assertEqual(name,"Toto");
        self.assertEqual(3,subMesh.getNumberOfCells());
        subConn2=[4,0,3,4,1,3,4,5,2,4,6,7,4,3]
        subConnIndex2=[0,5,9,14]
        self.assertEqual(14,subMesh.getNodalConnectivity().getNbOfElems());
        self.assertEqual(4,subMesh.getNodalConnectivityIndex().getNbOfElems());
        self.assertEqual(subConn2[0:14],list(subMesh.getNodalConnectivity().getValues()));
        self.assertEqual(subConnIndex2[0:4],list(subMesh.getNodalConnectivityIndex().getValues()));
        dd=DataArrayInt.New()
        dd.alloc(3,1)
        dd.iota(0)
        dd.setName("coucou")
        subMesh=subMesh.buildPartOfMySelf(dd,True);
        self.assertEqual("coucou",subMesh.getName());
        pass
    
    def testBuildPartOfMySelfNode(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        tab1=[5,7,8,4]
        subMesh=mesh.buildPartOfMySelfNode(tab1[0:4],True);
        self.assertTrue(isinstance(subMesh,MEDCouplingUMesh))
        self.assertEqual(1,len(subMesh.getAllGeoTypes()));
        self.assertEqual(NORM_QUAD4,subMesh.getAllGeoTypes()[0]);
        self.assertEqual(1,subMesh.getNumberOfCells());
        self.assertEqual(5,subMesh.getNodalConnectivity().getNbOfElems());
        self.assertEqual(2,subMesh.getNodalConnectivityIndex().getNbOfElems());
        subConn=[4,7,8,5,4]
        subConnIndex=[0,5]
        self.assertEqual(subConn[0:5],list(subMesh.getNodalConnectivity().getValues()));
        self.assertEqual(subConnIndex[0:2],list(subMesh.getNodalConnectivityIndex().getValues()));
        #
        ddd=DataArrayInt.New()
        ddd.setValues(tab1[0:2],2,1)
        ddd.setName("ddd")
        subMesh=mesh.buildPartOfMySelfNode(ddd,False);
        self.assertEqual("ddd",subMesh.getName())
        self.assertTrue(isinstance(subMesh,MEDCouplingUMesh))
        self.assertEqual(2,len(subMesh.getAllGeoTypes()));
        self.assertEqual(NORM_TRI3,subMesh.getAllGeoTypes()[0]);
        self.assertEqual(NORM_QUAD4,subMesh.getAllGeoTypes()[1]);
        self.assertEqual(3,subMesh.getNumberOfCells());
        self.assertEqual(14,subMesh.getNodalConnectivity().getNbOfElems());
        self.assertEqual(4,subMesh.getNodalConnectivityIndex().getNbOfElems());
        subConn2=[3,4,5,2,4,6,7,4,3,4,7,8,5,4]
        subConnIndex2=[0,4,9,14]
        self.assertEqual(subConn2[0:14],list(subMesh.getNodalConnectivity().getValues()));
        self.assertEqual(subConnIndex2[0:4],list(subMesh.getNodalConnectivityIndex().getValues()));
        #testing the case where length of tab2 is greater than max number of node per cell.
        tab2=[0,3,2,1,4,5,6]
        subMesh=mesh.buildPartOfMySelfNode(tab2[0:7],True);
        self.assertTrue(isinstance(subMesh,MEDCouplingUMesh))
        self.assertEqual(2,len(subMesh.getAllGeoTypes()));
        self.assertEqual(NORM_TRI3,subMesh.getAllGeoTypes()[0]);
        self.assertEqual(NORM_QUAD4,subMesh.getAllGeoTypes()[1]);
        self.assertEqual(3,subMesh.getNumberOfCells());
        pass
    
    def testZipCoords(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        self.assertEqual(2,len(mesh.getAllGeoTypes()));
        self.assertEqual(2,mesh.getSpaceDimension());
        self.assertEqual(9,mesh.getNumberOfNodes());
        self.assertEqual(5,mesh.getNumberOfCells());
        oldConn=mesh.getNodalConnectivity().getValues()[0:mesh.getNodalConnectivity().getNbOfElems()];
        oldConnIndex=mesh.getNodalConnectivityIndex().getValues()[0:mesh.getNumberOfCells()+1]
        oldCoords=mesh.getCoords();
        mesh.zipCoords();
        self.assertEqual(2,len(mesh.getAllGeoTypes()));
        self.assertEqual(2,mesh.getSpaceDimension());
        self.assertEqual(9,mesh.getNumberOfNodes());
        self.assertEqual(5,mesh.getNumberOfCells());
        self.assertEqual(mesh.getCoords().getValues()[0:2*9],oldCoords.getValues());
        self.assertEqual(list(oldConn),list(mesh.getNodalConnectivity().getValues()));
        self.assertEqual(list(oldConnIndex),list(mesh.getNodalConnectivityIndex().getValues()));
        #
        tab1=[0,4]
        subMesh=mesh.buildPartOfMySelf(tab1,True);
        self.assertTrue(isinstance(subMesh,MEDCouplingUMesh))
        traducer=subMesh.zipCoordsTraducer();
        expectedTraducer=[0, 1, -1, 2, 3, 4, -1, 5, 6]
        self.assertEqual(expectedTraducer,list(traducer.getValues()));
        self.assertEqual(NORM_QUAD4,subMesh.getAllGeoTypes()[0]);
        self.assertEqual(2,subMesh.getNumberOfCells());
        subConn=[4,0,2,3,1,4,5,6,4,3]
        subConnIndex=[0,5,10]
        self.assertEqual(7,subMesh.getNumberOfNodes());
        self.assertEqual(10,subMesh.getNodalConnectivity().getNbOfElems());
        self.assertEqual(3,subMesh.getNodalConnectivityIndex().getNbOfElems());
        self.assertEqual(subConn,list(subMesh.getNodalConnectivity().getValues()));
        self.assertEqual(subConnIndex,list(subMesh.getNodalConnectivityIndex().getValues()));
        #
        subMesh=mesh.buildPartOfMySelf(tab1,False);
        self.assertTrue(isinstance(subMesh,MEDCouplingUMesh))
        self.assertEqual(NORM_QUAD4,subMesh.getAllGeoTypes()[0]);
        self.assertEqual(2,subMesh.getNumberOfCells());
        self.assertEqual(7,subMesh.getNumberOfNodes());
        self.assertEqual(10,subMesh.getNodalConnectivity().getNbOfElems());
        self.assertEqual(3,subMesh.getNodalConnectivityIndex().getNbOfElems());
        self.assertEqual(subConn,list(subMesh.getNodalConnectivity().getValues()));
        self.assertEqual(subConnIndex,list(subMesh.getNodalConnectivityIndex().getValues()));
        pass
    
    def testZipConnectivity(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells1=[2,3,4]
        m3=m2.buildPartOfMySelf(cells1,True);
        self.assertTrue(isinstance(m3,MEDCouplingUMesh))
        m4=MEDCouplingDataForTest.build2DSourceMesh_1();
        m5=MEDCouplingUMesh.MergeUMeshes(m1,m3);
        m6=MEDCouplingUMesh.MergeUMeshes(m5,m4);
        #
        self.assertEqual(10,m6.getNumberOfCells());
        self.assertEqual(22,m6.getNumberOfNodes());
        (arr,areNodesMerged,newNbOfNodes)=m6.mergeNodes(1e-13);
        self.assertTrue(areNodesMerged);
        self.assertEqual(10,m6.getNumberOfCells());
        self.assertEqual(9,m6.getNumberOfNodes());
        #
        arr=m6.zipConnectivityTraducer(0);
        self.assertEqual(7,m6.getNumberOfCells());
        m7=m6.clone(True);
        arr=m6.zipConnectivityTraducer(0);
        self.assertTrue(m7.isEqual(m6,1e-12));
        self.assertEqual(7,m6.getNumberOfCells());
        pass
    
    def testEqualMesh(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_1();
        #
        self.assertTrue(mesh1.isEqual(mesh1,1e-12));
        #
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh2.isEqual(mesh1,1e-12));
        pt=mesh2.getCoords().getValues();
        tmp=pt[1]
        mesh2.getCoords().setIJ(0,1,5.999);
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(not mesh2.isEqual(mesh1,1e-12));
        mesh2.getCoords().setIJ(0,1,tmp);
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh2.isEqual(mesh1,1e-12));
        #
        pt2=mesh1.getNodalConnectivity().getValues();
        mesh1.getNodalConnectivity().setIJ(5,0,int(pt2[5])+1);
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(not mesh2.isEqual(mesh1,1e-12));
        mesh1.getNodalConnectivity().setIJ(5,0,int(pt2[5]));
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh2.isEqual(mesh1,1e-12));
        #
        pt2=mesh1.getNodalConnectivityIndex().getValues();
        mesh1.getNodalConnectivityIndex().setIJ(1,0,int(pt2[1]+1));
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(not mesh2.isEqual(mesh1,1e-12));
        mesh1.getNodalConnectivityIndex().setIJ(1,0,int(pt2[1]));
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh2.isEqual(mesh1,1e-12));
        #
        tmp3=mesh1.getName();
        mesh1.setName("lllll");
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(not mesh2.isEqual(mesh1,1e-12));
        mesh1.setName(tmp3);
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh2.isEqual(mesh1,1e-12));
        #
        tmp3=mesh2.getCoords().getInfoOnComponent(1);
        mesh2.getCoords().setInfoOnComponent(1,"kkkkkk");
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(not mesh2.isEqual(mesh1,1e-12));
        mesh2.getCoords().setInfoOnComponent(1,tmp3);
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh2.isEqual(mesh1,1e-12));
        pass
    
    def testEqualFieldDouble(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_1();
        #
        fieldOnCells1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        fieldOnCells1.setMesh(mesh1);
        fieldOnCells2=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        fieldOnCells2.setMesh(mesh2);
        #
        self.assertTrue(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        #
        fieldOnNodes1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        self.assertTrue(not fieldOnCells1.isEqual(fieldOnNodes1,1e-12,1e-15));
        self.assertTrue(not fieldOnNodes1.isEqual(fieldOnCells1,1e-12,1e-15));
        #
        fieldOnCells2=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        self.assertEqual(fieldOnCells2.getMesh(),None) # to check that convertMesh wrapping do not raise but return Py_None
        self.assertTrue(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        self.assertTrue(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells1.setTime(4.,6,7);
        self.assertTrue(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells2.setTime(4.,6,7);
        self.assertTrue(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells1.setName("Power");
        self.assertTrue(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells2.setName("Power");
        self.assertTrue(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        #
        fieldOnCells1.setMesh(mesh1);
        self.assertTrue(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells2.setMesh(mesh1);
        self.assertTrue(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        arr=DataArrayDouble.New();
        arr.setName("popo");
        arr.setValues(mesh1.getNumberOfCells()*3*[6.],mesh1.getNumberOfCells(),3);
        fieldOnCells1.setArray(arr);
        self.assertTrue(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells2.setArray(arr);
        self.assertTrue(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        #
        arr2=arr.deepCopy();
        fieldOnCells2.setArray(arr2);
        self.assertTrue(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        arr.setIJ(1,2,6.1);
        self.assertTrue(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        arr.setIJ(1,2,6.);
        self.assertTrue(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        arr2.setName("popo2");
        self.assertTrue(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        #
        arr2.setName("popo");
        self.assertTrue(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        #
        arr2.setInfoOnComponent(2,"jjj");
        self.assertTrue(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        arr.setInfoOnComponent(2,"jjj");
        self.assertTrue(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.assertTrue(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        pass

    def testNatureChecking(self):
        field=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        field.setNature(ExtensiveMaximum);
        field.setNature(IntensiveMaximum);
        field.setNature(ExtensiveConservation);
        field=MEDCouplingFieldDouble.New(ON_NODES,NO_TIME);
        field.setNature(IntensiveMaximum);
        self.assertRaises(InterpKernelException,field.setNature,ExtensiveMaximum);
        self.assertRaises(InterpKernelException,field.setNature,ExtensiveConservation);
        pass
      
    def testNatureOperations(self):
        """ Check nature constraints on field operations """
        m = MEDCouplingCMesh()
        m.setCoordsAt(0, DataArrayDouble([1.0,2.0,3.0]))
        m.setCoordsAt(1, DataArrayDouble([1.0,2.0,3.0]))
        m = m.buildUnstructured()
        f1, f2 = MEDCouplingFieldDouble.New(ON_CELLS, NO_TIME), MEDCouplingFieldDouble.New(ON_CELLS, NO_TIME)
        f1.setNature(ExtensiveMaximum)
        f2.setNature(IntensiveMaximum)
        self.assertEqual(ExtensiveMaximum, f1.getNature())
        self.assertEqual(IntensiveMaximum, f2.getNature())
        
        da = DataArrayDouble([1.0,2.0,3.0,4.0])
        f1.setMesh(m); f2.setMesh(m)
        f1.setArray(da); f2.setArray(da.deepCopy())
        # All this should complain about nature:
        self.assertRaises(InterpKernelException, f1.__add__, f2)
        self.assertRaises(InterpKernelException, f1.__iadd__, f2)
        self.assertRaises(InterpKernelException, f1.__sub__, f2)
        self.assertRaises(InterpKernelException, f1.__isub__, f2)
        self.assertRaises(InterpKernelException, f1.__radd__, f2)
        self.assertRaises(InterpKernelException, f1.__rsub__, f2)
        self.assertRaises(InterpKernelException, MEDCouplingFieldDouble.AddFields, f1, f2)
        self.assertRaises(InterpKernelException, MEDCouplingFieldDouble.SubstractFields, f1, f2)
        self.assertRaises(InterpKernelException, MEDCouplingFieldDouble.MaxFields, f1, f2)
        self.assertRaises(InterpKernelException, MEDCouplingFieldDouble.MinFields, f1, f2)
        # Not those ones:
        f3 = MEDCouplingFieldDouble.MultiplyFields(f1,f2)
        self.assertEqual(NoNature, f3.getNature())
        f3 = f1*f2
        self.assertEqual(NoNature, f3.getNature())
        f1Tmp = f1.deepCopy(); f1Tmp.setMesh(m);  f1Tmp *= f2
        self.assertEqual(NoNature, f1Tmp.getNature())
        f3 = MEDCouplingFieldDouble.DivideFields(f1,f2)
        self.assertEqual(NoNature, f3.getNature())
        f3 = f1/f2
        self.assertEqual(NoNature, f3.getNature())
        f1Tmp = f1.deepCopy();  f1Tmp.setMesh(m);  f1Tmp /= f2
        self.assertEqual(NoNature, f1Tmp.getNature())
#         f3 = MEDCouplingFieldDouble.PowFields(f1,f2)
#         self.assertEqual(NoNature, f3.getNature())
        f3 = f1**f2
        self.assertEqual(NoNature, f3.getNature())
        f1Tmp = f1.deepCopy();  f1Tmp.setMesh(m);  f1Tmp **= f2
        self.assertEqual(NoNature, f1Tmp.getNature())
        f3 = MEDCouplingFieldDouble.DotFields(f1,f2)
        self.assertEqual(NoNature, f3.getNature())
        f3 = f1.dot(f2)
        self.assertEqual(NoNature, f3.getNature())
        
        da = DataArrayDouble.Meld([da, da, da])
        f1.setArray(da); f2.setArray(da.deepCopy())
        f3 = MEDCouplingFieldDouble.CrossProductFields(f1,f2)
        self.assertEqual(NoNature, f3.getNature())
        f3 = f1.crossProduct(f2)

    def testBuildSubMeshData(self):
        targetMesh=MEDCouplingDataForTest.build2DTargetMesh_1()
        #check buildSubMesh on field on cells
        fieldCells=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        fieldCells.setMesh(targetMesh);
        elts=[1,2,4]
        ret1,di=fieldCells.buildSubMeshData(elts);
        self.assertTrue(isinstance(ret1,MEDCouplingUMesh))
        self.assertEqual(3,ret1.getNumberOfCells());
        self.assertEqual(9,ret1.getNumberOfNodes());
        self.assertEqual(3,di.getNumberOfTuples());
        self.assertEqual(1,di.getNumberOfComponents());
        toCheck=di.getValues();
        self.assertTrue(elts,toCheck);
        #check buildSubMesh on field on nodes
        fieldNodes=MEDCouplingFieldDouble.New(ON_NODES,NO_TIME);
        fieldNodes.setMesh(targetMesh);
        ret2,di=fieldNodes.buildSubMeshData(elts);
        self.assertTrue(isinstance(ret2,MEDCouplingUMesh))
        self.assertEqual(3,ret2.getNumberOfCells());
        self.assertEqual(6,ret2.getNumberOfNodes());
        self.assertEqual(6,di.getNumberOfTuples());
        self.assertEqual(1,di.getNumberOfComponents());
        toCheck=di.getValues();
        expected=[1,2,4,5,7,8]
        self.assertEqual(expected,list(toCheck));
        pass
    
    def testExtrudedMesh1(self):
        mesh3D,mesh2D=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        ext=MEDCouplingMappedExtrudedMesh.New(mesh3D,mesh2D,1);
        self.assertEqual(18,ext.getNumberOfCells());
        self.assertEqual(60,ext.getNumberOfNodes());
        ids3D=ext.getMesh3DIds();
        ids3DExpected=[5,4,3,2,1,0, 11,10,9,8,7,6, 17,16,15,14,13,12]
        self.assertEqual(18,ids3D.getNumberOfTuples());
        self.assertEqual(1,ids3D.getNumberOfComponents());
        self.assertEqual(ids3DExpected,list(ids3D.getValues()));
        mesh1D=ext.getMesh1D();
        self.assertEqual(4,mesh1D.getNumberOfNodes());
        self.assertEqual(3,mesh1D.getNumberOfCells());
        mesh1DExpected=[0.66666666666666663, 1.4583333333333333, 0, 0.66666666666666663,
                        1.4583333333333333, 1, 0.66666666666666663, 1.4583333333333333,
                        2, 0.66666666666666663, 1.4583333333333333, 3]
        mesh1DCoords=mesh1D.getCoords();
        self.assertEqual(4,mesh1DCoords.getNumberOfTuples());
        self.assertEqual(3,mesh1DCoords.getNumberOfComponents());
        self.assertEqual(mesh1DExpected,mesh1DCoords.getValues());
        conn1D=mesh1D.getNodalConnectivity();
        self.assertEqual(9,conn1D.getNumberOfTuples());
        self.assertEqual(1,conn1D.getNumberOfComponents());
        conn1DExpected=[1,0,1,1,1,2,1,2,3]
        self.assertEqual(conn1DExpected,list(conn1D.getValues()));
        pass

    def testExtrudedMesh3(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m1.changeSpaceDimension(3);
        m2=MEDCouplingDataForTest.buildCU1DMesh_U();
        m2.changeSpaceDimension(3);
        center=[0.,0.,0.]
        vector=[0.,1.,0.]
        m2.rotate(center,vector,-pi/2.);
        m3=m1.buildExtrudedMesh(m2,0);
        #
        m4=MEDCouplingMappedExtrudedMesh.New(m3,m1,0);
        self.assertEqual(15,m4.getNumberOfCells());
        self.assertEqual(5,m4.getMesh2D().getNumberOfCells());
        self.assertEqual(3,m4.getMesh1D().getNumberOfCells());
        m3DIds=m4.getMesh3DIds().getValues();
        self.assertEqual(list(range(15)), list(m3DIds));
        #some random in cells to check that extrusion alg find it correctly
        expected1=[1,3,2,0,6,5,7,10,11,8,12,9,14,13,4]
        m3.renumberCells(expected1,False);
        m4=MEDCouplingMappedExtrudedMesh.New(m3,m1,0);
        self.assertEqual(15,m4.getNumberOfCells());
        self.assertEqual(5,m4.getMesh2D().getNumberOfCells());
        self.assertEqual(3,m4.getMesh1D().getNumberOfCells());
        m3DIds=m4.getMesh3DIds().getValues();
        self.assertEqual(expected1,list(m3DIds));
        #play with polygons and polyedrons
        cells=[2,3]
        m1.convertToPolyTypes(cells);
        m3=m1.buildExtrudedMesh(m2,0);
        self.assertEqual(NORM_HEXA8,m3.getTypeOfCell(0));
        self.assertEqual(NORM_PENTA6,m3.getTypeOfCell(1));
        self.assertEqual(NORM_POLYHED,m3.getTypeOfCell(2));
        self.assertEqual(NORM_POLYHED,m3.getTypeOfCell(3));
        self.assertEqual(NORM_HEXA8,m3.getTypeOfCell(4));
        m3.renumberCells(expected1,False);
        m4=MEDCouplingMappedExtrudedMesh.New(m3,m1,0);
        self.assertEqual(15,m4.getNumberOfCells());
        self.assertEqual(5,m4.getMesh2D().getNumberOfCells());
        self.assertEqual(3,m4.getMesh1D().getNumberOfCells());
        m3DIds=m4.getMesh3DIds().getValues();
        self.assertEqual(expected1,list(m3DIds));
        pass

    def testExtrudedMesh4(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells=[2,4];
        m1.convertToPolyTypes(cells);
        m1.changeSpaceDimension(3);
        m2=MEDCouplingDataForTest.buildCU1DMesh_U();
        m2.changeSpaceDimension(3);
        center=[0.,0.,0.]
        vector=[0.,1.,0.]
        m2.rotate(center,vector,-pi/2.);
        m1.zipCoords()
        m3=m1.buildExtrudedMesh(m2,0);
        expected1=[1,3,2,0,6,5,7,10,11,8,12,9,14,13,4]
        rexpected1=[3, 0, 2, 1, 14, 5, 4, 6, 9, 11, 7, 8, 10, 13, 12]
        m3.renumberCells(expected1,False);
        m4=MEDCouplingMappedExtrudedMesh.New(m3,m1,0);
        self.assertEqual(NORM_HEXA8,m4.getTypeOfCell(0));
        self.assertEqual(NORM_HEXA8,m4.getTypeOfCell(1));
        self.assertEqual(NORM_POLYHED,m4.getTypeOfCell(2));
        self.assertEqual(NORM_PENTA6,m4.getTypeOfCell(7));
        f=m4.getMeasureField(True);
        arr=f.getArray();
        self.assertEqual(15,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        arrPtr=arr.getValues();
        expected2=[0.075,0.0375,0.0375,0.075,0.075,
                   0.1125,0.05625,0.05625,0.1125,0.1125,
                   0.0625,0.03125,0.03125,0.0625,0.0625]
        for i in range(15):
            self.assertAlmostEqual(expected2[rexpected1[i]],arrPtr[i],16);
            pass
        m5=m4.build3DUnstructuredMesh();
        m5.zipCoords()
        self.assertTrue(m5.isEqual(m3,1e-12));
        f=m5.getMeasureField(True);
        f.setMesh(m4)
        self.assertTrue(isinstance(f.getMesh(),MEDCouplingMappedExtrudedMesh))
        arr=f.getArray();
        arrPtr=arr.getValues();
        for i in range(15):
            self.assertAlmostEqual(expected2[rexpected1[i]],arrPtr[i],15);
            pass
        pass

    def testFindCommonNodes(self):
        targetMesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        comm,commI=targetMesh.findCommonNodes(1e-10,-1);
        self.assertEqual(1,commI.getNumberOfTuples());
        self.assertEqual(0,comm.getNumberOfTuples());
        o2n,newNbOfNodes=targetMesh.buildNewNumberingFromCommonNodesFormat(comm,commI);
        self.assertEqual(27,newNbOfNodes);
        self.assertEqual(27,o2n.getNumberOfTuples());
        o2nExp1 = list(range(27))
        self.assertEqual(o2nExp1,list(o2n.getValues()));
        #
        targetMesh=MEDCouplingDataForTest.build3DTargetMeshMergeNode_1();
        self.assertEqual(31,targetMesh.getNumberOfNodes());
        comm,commI=targetMesh.findCommonNodes(1e-10);# testing default parameter
        self.assertEqual(3,commI.getNumberOfTuples());
        self.assertEqual(6,comm.getNumberOfTuples());
        commExpected=[1,27,28,29,23,30]
        commIExpected=[0,4,6]
        self.assertEqual(commExpected,list(comm.getValues()));
        self.assertEqual(commIExpected,list(commI.getValues()));
        o2n,newNbOfNodes=targetMesh.buildNewNumberingFromCommonNodesFormat(comm,commI);
        self.assertEqual(31,o2n.getNumberOfTuples());
        self.assertEqual(27,newNbOfNodes);
        o2nExp2=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
                 21,22,23,24,25,26,1,1,1,23]
        self.assertEqual(o2nExp2,list(o2n.getValues()));
        #
        targetMesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        time=targetMesh.getTimeOfThis();
        o2n,areNodesMerged,newNbOfNodes=targetMesh.mergeNodes(1e-10);
        targetMesh.updateTime();
        self.assertEqual(time,targetMesh.getTimeOfThis());
        self.assertTrue(not areNodesMerged);
        #
        targetMesh=MEDCouplingDataForTest.build3DTargetMeshMergeNode_1();
        time=targetMesh.getTimeOfThis();
        o2n,areNodesMerged,newNbOfNodes=targetMesh.mergeNodes(1e-10);
        targetMesh.updateTime();
        self.assertTrue(time!=targetMesh.getTimeOfThis());
        self.assertTrue(areNodesMerged);
        connExp=[18,0,1,4,3,9,10,13,12, 18,1,2,5,4,10,11,14,13, 18,3,4,7,6,12,13,16,15,
                 18,4,5,8,7,13,14,17,16,
                 18,9,10,13,12,18,19,22,21, 18,10,11,14,13,19,20,23,22, 18,12,13,16,15,21,22,25,24,
                 18,13,14,17,16,22,23,26,25]
        self.assertEqual(72,targetMesh.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(connExp,list(targetMesh.getNodalConnectivity().getValues()));
        self.assertEqual(27,targetMesh.getCoords().getNumberOfTuples());
        coordsExp=[ 0., 0., 0., 50., 0., 0. , 200., 0., 0.  , 0., 50., 0., 50., 50., 0. ,
                    200., 50., 0.,   0., 200., 0., 50., 200., 0. , 200., 200., 0. ,
                    0., 0., 50., 50., 0., 50. , 200., 0., 50.  , 0., 50., 50., 50.,
                    50., 50. , 200., 50., 50.,   0., 200., 50., 50., 200., 50. ,
                    200., 200., 50. , 0., 0., 200., 50., 0., 200. , 200., 0., 200.  
                    , 0., 50., 200., 50., 50., 200. , 200., 50., 200., 
                    0., 200., 200., 50., 200., 200. , 200., 200., 200. ]
        self.assertEqual(coordsExp,targetMesh.getCoords().getValues());
        # 2D
        targetMesh=MEDCouplingDataForTest.build2DTargetMeshMergeNode_1();
        self.assertEqual(18,targetMesh.getNumberOfNodes());
        time=targetMesh.getTimeOfThis();
        o2n,areNodesMerged,newNbOfNodes=targetMesh.mergeNodes(1e-10);
        self.assertTrue(time!=targetMesh.getTimeOfThis());
        self.assertTrue(areNodesMerged);
        self.assertEqual(9,targetMesh.getNumberOfNodes());
        connExp2=[4,0,4,3,1, 3,1,3,2, 3,3,5,2, 4,4,6,7,3, 4,7,8,5,3]
        self.assertEqual(23,targetMesh.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(connExp2,list(targetMesh.getNodalConnectivity().getValues()));
        coordsExp2=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, 0.2,0.2, -0.3,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7]
        self.assertEqual(9,targetMesh.getCoords().getNumberOfTuples());
        self.assertEqual(coordsExp2,targetMesh.getCoords().getValues());
        pass

    def testCheckButterflyCells(self):
        sourceMesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells=sourceMesh.checkButterflyCells();
        self.assertEqual(0,len(cells));
        conn=sourceMesh.getNodalConnectivity()
        tmp=conn.getIJ(15,0)
        conn.setIJ(15,0,conn.getIJ(16,0))
        conn.setIJ(16,0,tmp)
        cells=sourceMesh.checkButterflyCells();
        self.assertEqual(1,len(cells));
        self.assertEqual([3],cells.getValues());
        tmp=conn.getIJ(15,0)
        conn.setIJ(15,0,conn.getIJ(16,0))
        conn.setIJ(16,0,tmp)
        cells=sourceMesh.checkButterflyCells();
        self.assertEqual(0,len(cells));
        # 3D surf
        sourceMesh=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        cells=sourceMesh.checkButterflyCells();
        self.assertEqual(0,len(cells));
        conn=sourceMesh.getNodalConnectivity()
        tmp=conn.getIJ(15,0)
        conn.setIJ(15,0,conn.getIJ(16,0))
        conn.setIJ(16,0,tmp)
        cells=sourceMesh.checkButterflyCells();
        self.assertEqual(1,len(cells));
        self.assertEqual([3],cells.getValues());
        tmp=conn.getIJ(15,0)
        conn.setIJ(15,0,conn.getIJ(16,0))
        conn.setIJ(16,0,tmp)
        cells=sourceMesh.checkButterflyCells();
        self.assertEqual(0,len(cells));
        pass

    def testMergeMesh1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DSourceMesh_1();
        vec=[1.,0.]
        m2.translate(vec);
        m3=m1.mergeMyselfWith(m2);
        self.assertTrue(isinstance(m3,MEDCouplingUMesh));
        m3.checkConsistencyLight();
        m4=MEDCouplingDataForTest.build2DTargetMeshMerged_1();
        self.assertTrue(m3.isEqual(m4,1.e-12));
        da,isMerged,newNbOfNodes=m3.mergeNodes(1.e-12);
        self.assertEqual(11,m3.getNumberOfNodes());
        self.assertTrue(isMerged);
        pass

    def testMergeMeshOnSameCoords1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells = list(range(5));
        m2.convertToPolyTypes(cells);
        m1.tryToShareSameCoords(m2,1e-12);
        m3=MEDCouplingDataForTest.build2DTargetMesh_1();
        m3.tryToShareSameCoords(m2,1e-12);
        meshes=[m1,m2,m3]
        m4=MEDCouplingUMesh.MergeUMeshesOnSameCoords(meshes);
        m4.checkConsistencyLight();
        self.assertEqual(15,m4.getNumberOfCells());
        cells1=[0,1,2,3,4]
        m1_1=m4.buildPartOfMySelf(cells1,True);
        m1_1.setName(m1.getName());
        self.assertTrue(m1.isEqual(m1_1,1e-12));
        cells2=[5,6,7,8,9]
        m2_1=m4.buildPartOfMySelf(cells2,True);
        m2_1.setName(m2.getName());
        self.assertTrue(m2.isEqual(m2_1,1e-12));
        cells3=[10,11,12,13,14]
        m3_1=m4.buildPartOfMySelf(cells3,True);
        m3_1.setName(m3.getName());
        self.assertTrue(m3.isEqual(m3_1,1e-12));
        pass

    def testMergeField1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DSourceMesh_1();
        vec=[1.,0.]
        m2.translate(vec);
        f1=m1.getMeasureField(True);
        f2=m2.getMeasureField(True);
        f3=MEDCouplingFieldDouble.MergeFields(f1,f2);
        f3.checkConsistencyLight();
        m4=MEDCouplingDataForTest.build2DTargetMeshMerged_1();
        self.assertTrue(f3.getMesh().isEqual(m4,1.e-12));
        name=f3.getName();
        self.assertEqual(name,"MeasureOfMesh_");
        self.assertEqual(f3.getTypeOfField(),ON_CELLS);
        self.assertEqual(f3.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f3.getNumberOfComponents());
        self.assertEqual(7,f3.getNumberOfTuples());
        values=[0.25,0.125,0.125,0.25,0.25,0.5,0.5]
        tmp=f3.getArray().getValues();
        self.assertEqual(len(values),len(tmp))
        for i in range(7):
            self.assertTrue(abs(values[i]-tmp[i])<1e-12)
            pass
        pass

    def testFillFromAnalytic(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        m.setTime(3.4,5,6); m.setTimeUnit("us");
        f1=m.fillFromAnalytic(ON_CELLS,1,"x+y");
        self.assertAlmostEqual(3.4,f1.getTime()[0],12) ; self.assertEqual(5,f1.getTime()[1]) ; self.assertEqual(6,f1.getTime()[2])
        self.assertEqual("us",f1.getTimeUnit())
        f1.checkConsistencyLight();                    
        self.assertEqual(f1.getTypeOfField(),ON_CELLS);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        values1=[-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values1),len(tmp))
        for i in range(len(tmp)):
            self.assertTrue(abs(tmp[i]-values1[i])<1.e-12)
            pass
        #
        f1=m.fillFromAnalytic(ON_NODES,1,"x+y");
        f1.checkConsistencyLight();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values2=[-0.6,-0.1,0.4,-0.1,0.4,0.9,0.4,0.9,1.4]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values2),len(tmp))
        for i in range(len(tmp)):
            self.assertTrue(abs(tmp[i]-values2[i])<1.e-12)
            pass
        #
        f1=m.fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+(2*(x+y))*JVec");
        f1.checkConsistencyLight();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(2,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values3=[-0.6,-1.2,-0.1,-0.2,0.4,0.8,-0.1,-0.2,0.4,0.8,0.9,1.8,0.4,0.8,0.9,1.8,1.4,2.8]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values3),len(tmp))
        for i in range(len(tmp)):
            self.assertTrue(abs(tmp[i]-values3[i])<1.e-12)
            pass
        values4=f1.accumulate();
        self.assertEqual(2,len(values4))
        self.assertTrue(abs(3.6-values4[0])<1.e-12);
        self.assertTrue(abs(7.2-values4[1])<1.e-12);
        values4=f1.integral(True);
        self.assertEqual(2,len(values4))
        self.assertTrue(abs(0.5-values4[0])<1.e-12);
        self.assertTrue(abs(1.-values4[1])<1.e-12);
        #
        self.assertRaises(InterpKernelException,m.fillFromAnalytic,ON_NODES,1,"1./(x-0.2)");
        pass

    def testFillFromAnalytic2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=m.fillFromAnalytic(ON_CELLS,1,"y+x");
        f1.checkConsistencyLight();
        self.assertEqual(f1.getTypeOfField(),ON_CELLS);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        values1=[-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values1),len(tmp))
        for i in range(len(values1)):
            self.assertTrue(abs(values1[i]-tmp[i])<1.e-12);
            pass
        #
        f1=m.fillFromAnalytic(ON_NODES,1,"y+2*x");
        f1.checkConsistencyLight();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values2=[-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values2),len(tmp))
        for i in range(len(values2)):
            self.assertTrue(abs(values2[i]-tmp[i])<1.e-12);
            pass
        f1=m.fillFromAnalytic(ON_NODES,1,"2.*x+y");
        f1.checkConsistencyLight();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        tmp=f1.getArray().getValues();
        values2Bis=[-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1]
        self.assertEqual(len(values2Bis),len(tmp))
        for i in range(len(values2Bis)):
            self.assertTrue(abs(values2Bis[i]-tmp[i])<1.e-12);
            pass
        #
        f1=m.fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+2*(x+y)*JVec");
        f1.checkConsistencyLight();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(2,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values3=[-0.6,-1.2,-0.1,-0.2,0.4,0.8,-0.1,-0.2,0.4,0.8,0.9,1.8,0.4,0.8,0.9,1.8,1.4,2.8]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values3),len(tmp))
        for i in range(len(values3)):
            self.assertTrue(abs(values3[i]-tmp[i])<1.e-12);
            pass
        values4=f1.accumulate();
        self.assertTrue(abs(3.6-values4[0])<1.e-12);
        self.assertTrue(abs(7.2-values4[1])<1.e-12);
        values4=f1.integral(True);
        self.assertTrue(abs(0.5-values4[0])<1.e-12);
        self.assertTrue(abs(1.-values4[1])<1.e-12);
        pass

    def testApplyFunc(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=m.fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+(2*(x+y))*JVec");
        f1.checkConsistencyLight();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(2,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        f1.applyFunc(1,"x+y");
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values1=[-1.8,-0.3,1.2,-0.3,1.2,2.7,1.2,2.7,4.2]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values1),len(tmp))
        for i in range(len(tmp)):
            self.assertTrue(abs(tmp[i]-values1[i])<1.e-12)
            pass
        pass

    def testApplyFunc2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=m.fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+2*(x+y)*JVec");
        f1.checkConsistencyLight();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(2,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        #
        f2=f1.clone(True);
        self.assertRaises(InterpKernelException, f2.applyFunc, 1, "a+b+c+d");
        self.assertRaises(InterpKernelException, f2.applyFunc, 1, "a/0");
        self.assertRaises(InterpKernelException, f2.applyFunc, "a/0");
        f2.applyFunc("abs(u)^2.4+2*u");
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(2,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values2=[-0.9065304805418678, -0.85105859001709905, -0.19601892829446504, -0.37898777756476987,
                 0.91090317490482353, 2.1853504664669781, -0.19601892829446504, -0.37898777756476987,
                 0.91090317490482353, 2.1853504664669781, 2.5765725275664879, 7.6987743736515295,
                 0.91090317490482353, 2.1853504664669781, 2.5765725275664879, 7.6987743736515295,
                 5.0423700574830965, 17.435300118916864]
        tmp=f2.getArray().getValues();
        self.assertEqual(len(tmp),len(values2))
        for i in range(len(tmp)):
            self.assertTrue(abs(tmp[i]-values2[i])<1.e-12)
            pass
        #
        f1.applyFunc(1,"x+y");
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values1=[-1.8,-0.3,1.2,-0.3,1.2,2.7,1.2,2.7,4.2]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(tmp),len(values1))
        for i in range(len(tmp)):
            self.assertTrue(abs(tmp[i]-values1[i])<1.e-12)
            pass
        pass

    def testOperationsOnFields(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=m.fillFromAnalytic(ON_NODES,1,"x+y");
        f2=m.fillFromAnalytic(ON_NODES,1,"x+y");
        f1.checkConsistencyLight();
        f2.checkConsistencyLight();
        f3=f1+f2;
        f3.checkConsistencyLight();
        self.assertEqual(f3.getTypeOfField(),ON_NODES);
        self.assertEqual(f3.getTimeDiscretization(),ONE_TIME);
        values1=[-1.2,-0.2,0.8,-0.2,0.8,1.8,0.8,1.8,2.8]
        tmp=f3.getArray().getValues();
        self.assertEqual(len(values1),len(tmp))
        for i in range(len(tmp)):
            self.assertTrue(abs(tmp[i]-values1[i])<1.e-12)
            pass
        #
        f3=f1*f2;
        f3.checkConsistencyLight();
        self.assertEqual(f3.getTypeOfField(),ON_NODES);
        self.assertEqual(f3.getTimeDiscretization(),ONE_TIME);
        values2=[0.36,0.01,0.16,0.01,0.16,0.81,0.16,0.81,1.96]
        tmp=f3.getArray().getValues();
        self.assertEqual(len(values2),len(tmp))
        for i in range(len(tmp)):
            self.assertTrue(abs(tmp[i]-values2[i])<1.e-12)
            pass
        #
        f3=f1+f2;
        f4=f1-f3;
        f4.checkConsistencyLight();
        self.assertEqual(f4.getTypeOfField(),ON_NODES);
        self.assertEqual(f4.getTimeDiscretization(),ONE_TIME);
        values3=[0.6,0.1,-0.4,0.1,-0.4,-0.9,-0.4,-0.9,-1.4]
        tmp=f4.getArray().getValues();
        self.assertEqual(len(values3),len(tmp))
        for i in range(len(tmp)):
            self.assertTrue(abs(tmp[i]-values3[i])<1.e-12)
            pass
        #
        f3=f1+f2;
        f4=f3/f2;
        f4.checkConsistencyLight();
        self.assertEqual(f4.getTypeOfField(),ON_NODES);
        self.assertEqual(f4.getTimeDiscretization(),ONE_TIME);
        tmp=f4.getArray().getValues();
        for i in range(len(tmp)):
            self.assertTrue(abs(tmp[i]-2.)<1.e-12)
            pass
        #
        f4=f2.buildNewTimeReprFromThis(NO_TIME,False);
        f4.checkConsistencyLight();
        self.assertEqual(f4.getTypeOfField(),ON_NODES);
        self.assertEqual(f4.getTimeDiscretization(),NO_TIME);
        self.assertRaises(InterpKernelException,f1.__add__,f4);
        f5=f4.buildNewTimeReprFromThis(ONE_TIME,False);
        self.assertEqual(f5.getTypeOfField(),ON_NODES);
        self.assertEqual(f5.getTimeDiscretization(),ONE_TIME);
        f3=f1+f5;
        tmp=f3.getArray().getValues();
        values4=[-1.2,-0.2,0.8,-0.2,0.8,1.8,0.8,1.8,2.8]
        self.assertEqual(len(values3),len(tmp))
        for i in range(len(tmp)):
            self.assertTrue(abs(tmp[i]-values4[i])<1.e-12)
            pass
        #
        f4=f2.buildNewTimeReprFromThis(NO_TIME,True);
        f4.checkConsistencyLight();
        self.assertEqual(f4.getTypeOfField(),ON_NODES);
        self.assertEqual(f4.getTimeDiscretization(),NO_TIME);
        self.assertRaises(InterpKernelException,f1.__add__,f4);
        f5=f4.buildNewTimeReprFromThis(ONE_TIME,True);
        self.assertEqual(f5.getTypeOfField(),ON_NODES);
        self.assertEqual(f5.getTimeDiscretization(),ONE_TIME);
        f3=f1+f5;
        tmp=f3.getArray().getValues();
        values5=[-1.2,-0.2,0.8,-0.2,0.8,1.8,0.8,1.8,2.8]
        self.assertEqual(len(values5),len(tmp))
        for i in range(len(tmp)):
            self.assertTrue(abs(tmp[i]-values5[i])<1.e-12)
            pass
        pass

    def testOperationsOnFields2(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m.setTime(3.4,5,6); m.setTimeUnit("us");
        f1=m.fillFromAnalytic(ON_NODES,1,"x+y+z");
        f2=m.fillFromAnalytic(ON_NODES,1,"a*a+b+c*c");
        f3=f1/f2;
        f3.checkConsistencyLight();
        self.assertEqual(f3.getTypeOfField(),ON_NODES);
        self.assertEqual(f3.getTimeDiscretization(),ONE_TIME);
        expected1=[-2.4999999999999991, 1.2162162162162162, 0.77868852459016391,
                   0.7407407407407407, 1.129032258064516, 0.81632653061224492,
                   0.86538461538461531, 1.0919540229885056, 0.84302325581395343]
        self.assertEqual(1,f3.getNumberOfComponents());
        self.assertEqual(9,f3.getNumberOfTuples());
        val=f3.getArray().getValues();
        for i in range(9):
            self.assertTrue(abs(expected1[i]-val[i])<1.e-12);
        #
        f1=m.buildOrthogonalField();
        self.assertAlmostEqual(3.4,f1.getTime()[0],12) ; self.assertEqual(5,f1.getTime()[1]) ; self.assertEqual(6,f1.getTime()[2])
        self.assertEqual("us",f1.getTimeUnit())
        f2=m.fillFromAnalytic(ON_CELLS,1,"x");
        f3=f1*f2;
        expected2=[-0.035355339059327376,0.,0.035355339059327376, 0.2592724864350674,0.,-0.2592724864350674, 0.37712361663282529,0.,-0.37712361663282529, -0.035355339059327376,0.,0.035355339059327376, 0.31819805153394637,0.,-0.31819805153394637]
        val=f3.getArray().getValues();
        for i in range(15):
            self.assertTrue(abs(expected2[i]-val[i])<1.e-12);
            pass
        #
        f3=f2*f1;
        val=f3.getArray().getValues();
        for i in range(15):
            self.assertTrue(abs(expected2[i]-val[i])<1.e-12);
            pass
        pass

    def testOperationsOnFields3(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        f1=m.fillFromAnalytic(ON_NODES,1,"x+y+z");
        f2=m.fillFromAnalytic(ON_NODES,1,"a*a+b+c*c");
        f1/=f2
        f1.checkConsistencyLight();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        expected1=[-2.4999999999999991, 1.2162162162162162, 0.77868852459016391,
                   0.7407407407407407, 1.129032258064516, 0.81632653061224492,
                   0.86538461538461531, 1.0919540229885056, 0.84302325581395343]
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        val=f1.getArray().getValues();
        for i in range(9):
            self.assertTrue(abs(expected1[i]-val[i])<1.e-12);
            pass
        #
        f1=m.buildOrthogonalField();
        f2=m.fillFromAnalytic(ON_CELLS,1,"x");
        f1*=f2
        expected2=[-0.035355339059327376,0.,0.035355339059327376, 0.2592724864350674,0.,-0.2592724864350674, 0.37712361663282529,0.,-0.37712361663282529, -0.035355339059327376,0.,0.035355339059327376, 0.31819805153394637,0.,-0.31819805153394637]
        val=f1.getArray().getValues();
        for i in range(15):
            self.assertTrue(abs(expected2[i]-val[i])<1.e-12);
            pass
        #
        f1=m.buildOrthogonalField();
        # to avoid valgrind leaks
        # self.assertRaises(InterpKernelException,f2.__imul__,f1);
        pass

    def testOperationsOnFields4(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        nbOfCells=m.getNumberOfCells();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,CONST_ON_TIME_INTERVAL);
        f1.setMesh(m);
        array=DataArrayDouble.New();
        f1.setArray(array);
        self.assertRaises(InterpKernelException,f1.setEndArray,array);
        self.assertRaises(InterpKernelException,f1.getEndArray);
        arr1=[0.,10.,20.,1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.]
        arr2=[5.,15.,25.,6.,16.,26.,7.,17.,27.,8.,18.,28.,9.,19.,29.]
        array.setValues(arr1,nbOfCells,3);
        f1.setStartTime(2.,0,0);
        f1.setEndTime(3.,0,0);
        f1.checkConsistencyLight();
        pos=[0.3,-0.2]
        res=f1.getValueOn(pos);
        self.assertTrue(abs(arr1[3]-res[0])<1.e-12);
        self.assertTrue(abs(arr1[4]-res[1])<1.e-12);
        self.assertTrue(abs(arr1[5]-res[2])<1.e-12);
        res=None
        res=f1.getValueOn(pos,2.2);
        self.assertTrue(abs(arr1[3]-res[0])<1.e-12);
        self.assertTrue(abs(arr1[4]-res[1])<1.e-12);
        self.assertTrue(abs(arr1[5]-res[2])<1.e-12);
        res=None
        self.assertRaises(InterpKernelException,f1.getValueOn,pos,3.2)
        f2=MEDCouplingFieldDouble.New(ON_CELLS,LINEAR_TIME);
        f2.setMesh(m);
        f2.setArray(f1.getArray());
        f2.setStartTime(2.,3,0);
        f2.setEndTime(4.,13,0);
        self.assertRaises(InterpKernelException,f2.checkConsistencyLight)
        array2=DataArrayDouble.New();
        array2.setValues(arr2,nbOfCells,3);
        f2.setEndArray(array2);
        f2.checkConsistencyLight();
        #
        res=None
        res=f2.getValueOn(pos,3.21);
        self.assertTrue(abs(4.025-res[0])<1.e-12);
        self.assertTrue(abs(14.025-res[1])<1.e-12);
        self.assertTrue(abs(24.025-res[2])<1.e-12);
        f3=f2.clone(True);
        self.assertTrue(f2.isEqual(f3,1e-12,1e-12));
        f3.getEndArray().setIJ(0,0,5.001);
        self.assertTrue(not f2.isEqual(f3,1e-12,1e-12));
        self.assertTrue(f2.isEqual(f3,1e-12,1e-2));
        f3.setStartTime(2.1,3,0);
        self.assertTrue(not f2.isEqual(f3,1e-12,1e-2));
        f3.setStartTime(2.,3,0);
        self.assertTrue(f2.isEqual(f3,1e-12,1e-2));
        f3.setStartTime(2.,4,0);
        self.assertTrue(not f2.isEqual(f3,1e-12,1e-2));
        f3.setStartTime(2.,3,1);
        self.assertTrue(not f2.isEqual(f3,1e-12,1e-2));
        f3.setStartTime(2.,3,0);
        self.assertTrue(f2.isEqual(f3,1e-12,1e-2));
        f3.setEndTime(4.1,13,0);
        self.assertTrue(not f2.isEqual(f3,1e-12,1e-2));
        f3.setEndTime(4.,13,0);
        self.assertTrue(f2.isEqual(f3,1e-12,1e-2));
        f3.setEndTime(4.,14,0);
        self.assertTrue(not f2.isEqual(f3,1e-12,1e-2));
        f3.setEndTime(4.,13,1);
        self.assertTrue(not f2.isEqual(f3,1e-12,1e-2));
        f3.setEndTime(4.,13,0);
        self.assertTrue(f2.isEqual(f3,1e-12,1e-2));
        f4=f2+f2
        res=None
        res=f4.getValueOn(pos,3.21);
        self.assertTrue(abs(8.05-res[0])<1.e-12);
        self.assertTrue(abs(28.05-res[1])<1.e-12);
        self.assertTrue(abs(48.05-res[2])<1.e-12);
        f4+=f2;
        res=None
        res=f4.getValueOn(pos,3.21);
        self.assertTrue(abs(12.075-res[0])<1.e-12);
        self.assertTrue(abs(42.075-res[1])<1.e-12);
        self.assertTrue(abs(72.075-res[2])<1.e-12);
        pass
    
    def testMergeNodesOnField(self):
        targetMesh=MEDCouplingDataForTest.build3DTargetMeshMergeNode_1();
        f1=targetMesh.fillFromAnalytic(ON_NODES,1,"x+y+z");
        f1.mergeNodes(1e-10);
        #
        targetMesh=MEDCouplingDataForTest.build3DTargetMeshMergeNode_1();
        f1=targetMesh.fillFromAnalytic(ON_NODES,1,"x+y+z");
        tmp=f1.getArray()
        tmp.setIJ(0,0,1000.);
        f1.mergeNodes(1e-10);
        #
        targetMesh=MEDCouplingDataForTest.build3DTargetMeshMergeNode_1();
        f1=targetMesh.fillFromAnalytic(ON_NODES,1,"x+y+z");
        tmp=f1.getArray()
        tmp.setIJ(1,0,1000.);
        self.assertRaises(InterpKernelException,f1.mergeNodes,1.e-10)
        pass

    def testCheckConsecutiveCellTypes(self):
        sourceMesh=MEDCouplingDataForTest.build2DSourceMesh_1();
        targetMesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        self.assertTrue(sourceMesh.checkConsecutiveCellTypes());
        order1=[NORM_TRI3,NORM_QUAD4]
        order2=[NORM_QUAD4,NORM_TRI3]
        self.assertTrue(not targetMesh.checkConsecutiveCellTypes());
        self.assertTrue(not targetMesh.checkConsecutiveCellTypesAndOrder(order1));
        self.assertTrue(not targetMesh.checkConsecutiveCellTypesAndOrder(order2));
        da=targetMesh.getRenumArrForConsecutiveCellTypesSpec(order1);
        self.assertEqual(5,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        expected1=[2,0,1,3,4]
        self.assertTrue(expected1==list(da.getValues()));
        da=targetMesh.getRenumArrForConsecutiveCellTypesSpec(order2);
        self.assertEqual(5,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        expected2=[0,3,4,1,2]
        self.assertTrue(expected2==list(da.getValues()));
        renumber1=[4,0,1,2,3]
        targetMesh.renumberCells(renumber1,False);
        self.assertTrue(targetMesh.checkConsecutiveCellTypes());
        self.assertTrue(targetMesh.checkConsecutiveCellTypesAndOrder(order1));
        self.assertTrue(not targetMesh.checkConsecutiveCellTypesAndOrder(order2));
        pass

    def testRearrange2ConsecutiveCellTypes(self):
        m1_1=MEDCouplingDataForTest.build2DSourceMesh_1();
        m2_1=MEDCouplingDataForTest.build2DTargetMesh_1();
        arr1=m1_1.rearrange2ConsecutiveCellTypes();
        m1_2=MEDCouplingDataForTest.build2DSourceMesh_1();
        self.assertTrue(m1_2.isEqual(m1_1,1e-12));
        expected1=[0,1]
        self.assertEqual(2,arr1.getNumberOfTuples());
        self.assertEqual(1,arr1.getNumberOfComponents());
        self.assertEqual(expected1,arr1.getValues());
        expected2=[0,3,4,1,2]
        arr1=m2_1.rearrange2ConsecutiveCellTypes();
        self.assertEqual(5,arr1.getNumberOfTuples());
        self.assertEqual(1,arr1.getNumberOfComponents());
        self.assertEqual(expected2,list(arr1.getValues()));
        m2_2=MEDCouplingDataForTest.build2DTargetMesh_1();
        self.assertEqual(5,arr1.getNumberOfTuples());
        self.assertEqual(1,arr1.getNumberOfComponents());
        self.assertEqual(expected2,list(arr1.getValues()));
        self.assertTrue(not m2_2.isEqual(m2_1,1e-12));
        m2_2.renumberCells(expected2,False);
        self.assertTrue(m2_2.isEqual(m2_1,1e-12));
        pass

    def testSplitByType(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        v=m1.splitByType();
        self.assertEqual(3,len(v));
        m2=MEDCouplingUMesh.MergeUMeshesOnSameCoords(v);
        m2.setName(m1.getName());
        self.assertTrue(m1.isEqual(m2,1.e-12));
        pass

    def testFuseUMeshesOnSameCoords(self):
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells1=[2,3,4]
        m3=m2.buildPartOfMySelf(cells1,True);
        self.assertTrue(isinstance(m3,MEDCouplingUMesh))
        cells2=[1,2,4]
        m4=m2.buildPartOfMySelf(cells2,True);
        self.assertTrue(isinstance(m4,MEDCouplingUMesh))
        cells3=[1,2]
        m5=m2.buildPartOfMySelf(cells3,True);
        self.assertTrue(isinstance(m5,MEDCouplingUMesh))
        meshes=[m3,m4,m5]
        #
        m7,corr=MEDCouplingUMesh.FuseUMeshesOnSameCoords(meshes,0);
        self.assertEqual(4,m7.getNumberOfCells());
        self.assertEqual(3,len(corr));
        expectedVals1=[3,3,2]
        expectedVals2=[[0,1,2],[3,0,2],[3,0]]
        for i in range(3):
            arr=corr[i];
            self.assertEqual(1,arr.getNumberOfComponents());
            nbOfVals=expectedVals1[i];
            self.assertEqual(nbOfVals,arr.getNumberOfTuples());
            vals=arr.getValues();
            self.assertEqual(expectedVals2[i],list(vals));
            pass
        arr2,fidsOfGroups=DataArrayInt.MakePartition(corr,m7.getNumberOfCells());
        fidExp=[5,1,3,4]
        fidsGrp=[[1,3,5],[3,4,5],[4,5]]
        self.assertEqual(3,len(fidsOfGroups));
        self.assertEqual(1,arr2.getNumberOfComponents());
        self.assertEqual(4,arr2.getNumberOfTuples());
        self.assertEqual(fidExp,list(arr2.getValues()));
        for i in range(3):
            nbOfVals=expectedVals1[i];
            self.assertEqual(list(fidsOfGroups[i]),fidsGrp[i]);
            pass
        pass

    def testFuseUMeshesOnSameCoords2(self):
        m1,m2=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        part1=[2,3,6,4,10]
        m3=m1.buildPartOfMySelf(part1,True);
        part2=[5,6,4,7]
        m4=m1.buildPartOfMySelf(part2,True);
        meshes=[m1,m3,m3,m4]
        m5,corr=MEDCouplingUMesh.FuseUMeshesOnSameCoords(meshes,0);
        self.assertEqual(18,m5.getNumberOfCells());
        exp2=[
            [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],
            [2,3,6,4,10],
            [2,3,6,4,10],
            [5,6,4,7]]
        i=0;
        for it in corr:
            self.assertEqual(exp2[i],list(it.getValues()));
            i+=1
            pass
        pass

    def testBuildOrthogonalField(self):
        targetMesh=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        field=targetMesh.buildOrthogonalField();
        expected=[0.70710678118654746,0.,-0.70710678118654746]
        self.assertEqual(5,field.getNumberOfTuples());
        self.assertEqual(3,field.getNumberOfComponents());
        vals=field.getArray().getValues();
        for i in range(15):
            self.assertTrue(abs(expected[i%3]-vals[i])<1e-12);
        # testing
        targetCoords=[0.,0.,0.,0.5,0.,0.5,1.,0.,1.,0.,1.,0.]
        targetConn=[0,1,2,3]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.setMeshDimension(2);
        targetMesh.allocateCells(1);
        targetMesh.insertNextCell(NORM_QUAD4,targetConn[0:4])
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,4,3);
        targetMesh.setCoords(myCoords);
        field=targetMesh.buildOrthogonalField();
        self.assertEqual(1,field.getNumberOfTuples());
        self.assertEqual(3,field.getNumberOfComponents());
        vals=field.getArray().getValues();
        self.assertTrue(abs(-0.70710678118654746-vals[0])<1e-12);
        self.assertTrue(abs(0.-vals[1])<1e-12);
        self.assertTrue(abs(0.70710678118654746-vals[2])<1e-12);
        pass

    def testGetCellsContainingPoint(self):
        targetMesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        pos=[0.,0.,0.4,0.4,0.,0.4,0.1,0.1,0.25,0.,0.65,0.]
        #2D basic
        t1,t2=targetMesh.getCellsContainingPoints(pos,6,1e-12);
        self.assertEqual(6,t1.getNumberOfTuples());
        self.assertEqual(1,t1.getNumberOfComponents());
        self.assertEqual(7,t2.getNumberOfTuples());
        self.assertEqual(1,t2.getNumberOfComponents());
        expectedValues1=[0,4,3,0,1,2]
        expectedValues2=[0,1,2,3,4,5,6]
        self.assertEqual(list(t1.getValues()),expectedValues1);
        self.assertEqual(list(t2.getValues()),expectedValues2);
        #2D with no help of bounding box.
        center=[0.2,0.2]
        MEDCouplingPointSet.Rotate2DAlg(center,0.78539816339744830962,6,pos);
        targetMesh.rotate(center,0.78539816339744830962);
        t1=None
        t2=None
        t1,t2=targetMesh.getCellsContainingPoints(pos,1e-12);
        self.assertEqual(6,t1.getNumberOfTuples());
        self.assertEqual(7,t2.getNumberOfTuples());
        self.assertEqual(list(t1.getValues()),expectedValues1);
        self.assertEqual(list(t2.getValues()),expectedValues2);
        t1,t2=targetMesh.getCellsContainingPoints(DataArrayDouble.New(pos,6,2),1e-12);
        self.assertEqual(6,t1.getNumberOfTuples());
        self.assertEqual(7,t2.getNumberOfTuples());
        self.assertEqual(list(t1.getValues()),expectedValues1);
        self.assertEqual(list(t2.getValues()),expectedValues2);
        self.assertRaises(InterpKernelException,targetMesh.getCellsContainingPoints,DataArrayDouble.New(pos,4,3),1e-12);
        #2D outside
        pos1bis=[-0.3303300858899107,-0.11819805153394641]
        self.assertEqual(-1,targetMesh.getCellContainingPoint(pos1bis,1e-12));
        #test limits 2D
        targetMesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        pos2=[0.2,-0.05]
        t1=None
        t1=targetMesh.getCellsContainingPoint(pos2,1e-12)
        self.assertEqual(2,len(t1));
        expectedValues3=[0,1]
        self.assertEqual(list(t1.getValues()),expectedValues3);
        pos3=[0.2,0.2]
        t1=None
        t1=targetMesh.getCellsContainingPoint(pos3,1e-12);
        self.assertEqual(5,len(t1));
        expectedValues4=[0,1,2,3,4]
        self.assertEqual(list(t1.getValues()),expectedValues4);
        self.assertEqual(0,targetMesh.getCellContainingPoint(pos3,1e-12));
        #3D
        targetMesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        pos4=[25.,25.,25.]
        self.assertEqual(0,targetMesh.getCellContainingPoint(pos4,1e-12));
        pos5=[50.,50.,50.]
        t1=None
        t1=targetMesh.getCellsContainingPoint(pos5,1e-12);
        self.assertEqual(8,len(t1));
        expectedValues5=[0,1,2,3,4,5,6,7]
        self.assertEqual(list(t1.getValues()),expectedValues5);
        pos6=[0., 50., 0.]
        t1=None
        t1=targetMesh.getCellsContainingPoint(pos6,1e-12);
        self.assertEqual(2,len(t1));
        expectedValues6=[0,2]
        self.assertEqual(list(t1.getValues()),expectedValues6);
        #3D outside
        pos7=[-1.0,-1.0,0.]
        self.assertEqual(-1,targetMesh.getCellContainingPoint(pos7,1e-12));
        #3D outside 2
        center2=[0.,0.,0.]
        vec2=[0.,-1.,0.]
        targetMesh.rotate(center2,vec2,0.78539816339744830962);
        pos8=[-25.,25.,12.]
        self.assertEqual(-1,targetMesh.getCellContainingPoint(pos8,1e-12));
        pass

    def testGetValueOn1(self):
        targetMesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        fieldOnCells=MEDCouplingFieldDouble.New(ON_CELLS);
        nbOfCells=targetMesh.getNumberOfCells();
        fieldOnCells.setMesh(targetMesh);
        array=DataArrayDouble.New();
        tmp=2*nbOfCells*[None]
        for i in range(nbOfCells):
            tmp[2*i]=7.+float(i);
            tmp[2*i+1]=17.+float(i)
            pass
        array.setValues(tmp,nbOfCells,2);
        fieldOnCells.setArray(array);
        #
        pos1=[0.25,0.]
        res=fieldOnCells.getValueOn(pos1);
        self.assertEqual(2,len(res))
        self.assertTrue(abs(8.-res[0])<1e-12);
        self.assertTrue(abs(18.-res[1])<1e-12);
        #
        #
        targetMesh=MEDCouplingDataForTest.build2DSourceMesh_1();
        fieldOnNodes=MEDCouplingFieldDouble.New(ON_NODES);
        nbOfNodes=targetMesh.getNumberOfNodes();
        fieldOnNodes.setMesh(targetMesh);
        array=DataArrayDouble.New();
        tmp=2*nbOfNodes*[None]
        for i in range(nbOfNodes):
            tmp[2*i]=17.+float(i);
            tmp[2*i+1]=27.+float(i)
            pass
        array.setValues(tmp,nbOfNodes,2);
        fieldOnNodes.setArray(array);
        #
        pos2=[-0.13333333333333333,-0.13333333333333333]
        res=None
        res=fieldOnNodes.getValueOn(pos2);
        self.assertEqual(2,len(res))
        self.assertTrue(abs(17.5-res[0])<1e-12);
        self.assertTrue(abs(27.5-res[1])<1e-12);
        pos3=[0.033333333333333326,0.36666666666666664]
        res=None
        res=fieldOnNodes.getValueOn(pos3);
        self.assertEqual(2,len(res))
        self.assertTrue(abs(18.666666666666667-res[0])<1e-12);
        self.assertTrue(abs(28.666666666666667-res[1])<1e-12);
        pass

    def testCMesh0(self):
        mesh=MEDCouplingCMesh.New();
        meshEmpty=mesh.clone(True);
        self.assertTrue(meshEmpty.isEqual(mesh, 1e-12));
        
        coordsX=DataArrayDouble.New();
        arrX=[ -1., 1., 2., 4. ]
        coordsX.setValues(arrX, 4, 1);
        coordsY=DataArrayDouble.New();
        arrY=[ -2., 2., 4., 8. ]
        coordsY.setValues(arrY, 4, 1);
        coordsZ=DataArrayDouble.New();
        arrZ=[ -3., 3., 6., 12. ]
        coordsZ.setValues(arrZ, 4, 1);
        mesh.setCoords(coordsX, coordsY, coordsZ);
        #
        fieldOnNodes=mesh.fillFromAnalytic(ON_NODES, 1, "x+y/2.+z/3.");
        self.assertEqual(1, fieldOnNodes.getNumberOfComponents());
        self.assertEqual(64, fieldOnNodes.getNumberOfTuples());
        expected1=[-3., -1., 0., 2., -1., 1., 2., 4., 0., 2., 3., 5., 2., 4., 5., 7., -1., 1., 2.,
                    4., 1., 3., 4., 6., 2., 4., 5., 7., 4., 6., 7., 9., 0., 2., 3., 5., 2., 4., 5.,
                    7., 3., 5., 6., 8., 5., 7., 8., 10., 2., 4., 5.,
                    7., 4., 6., 7., 9., 5., 7., 8., 10., 7., 9., 10., 12.];
        
        val=fieldOnNodes.getArray().getValues();
        for i in range(64):
          self.assertAlmostEqual(expected1[i], val[i], 12)
        res=fieldOnNodes.getValueOnPos(1, 3, 2);
        self.assertAlmostEqual(7., res[0], 12);
        #
        fieldOnCells=mesh.fillFromAnalytic(ON_CELLS, 1, "x+y/2.+z/3.");
        self.assertEqual(1, fieldOnCells.getNumberOfComponents());
        self.assertEqual(27, fieldOnCells.getNumberOfTuples());
        val=fieldOnCells.getArray().getValues();
        expected2=[0, 1.5, 3, 1.5, 3, 4.5, 3, 4.5, 6, 1.5, 3, 4.5, 3, 4.5,
                    6, 4.5, 6, 7.5, 3, 4.5, 6, 4.5, 6, 7.5, 6, 7.5, 9];
        for i in range(27):
          self.assertAlmostEqual(expected2[i], val[i], 12);
        #res=fieldOnCells.getValueOnPos(1,2,1);
        #self.assertAlmostEqual(6.,res,12);
        #
        meshDeepCopy=mesh.deepCopy();
        meshClone=mesh.clone(False);
        
        meshEmpty.copyTinyStringsFrom(mesh);
        #no data in meshEmpty, expected False
        self.assertTrue(not meshEmpty.isEqual(mesh, 1e-12));
        
        self.assertTrue(meshDeepCopy.isEqual(mesh, 1e-12));
        meshDeepCopy.copyTinyStringsFrom(mesh);
        self.assertTrue(meshDeepCopy.isEqual(mesh, 1e-12));
        self.assertTrue(meshClone.isEqual(mesh, 1e-12));
        
        self.assertEqual(CARTESIAN, mesh.getType());
        self.assertEqual(CARTESIAN, meshEmpty.getType());
        self.assertEqual(CARTESIAN, meshDeepCopy.getType());
        self.assertEqual(CARTESIAN, meshClone.getType());
        pass

    def testCMesh1(self):
        mesh1=MEDCouplingCMesh.New();
        coordsX1=DataArrayDouble.New();
        arrX1=[ -1., 1., 2., 4. ]
        coordsX1.setValues(arrX1, 4, 1);
        coordsY1=DataArrayDouble.New();
        arrY1=[ -2., 2., 4., 8. ]
        coordsY1.setValues(arrY1, 4, 1);
        coordsZ1=DataArrayDouble.New();
        arrZ1=[ -3., 3., 6., 12. ]
        coordsZ1.setValues(arrZ1, 4, 1);
        mesh1.setCoords(coordsX1, coordsY1, coordsZ1);
        
        mesh2=MEDCouplingCMesh.New();
        coordsX2=DataArrayDouble.New();
        arrX2=[ -1., 1., 2., 4. ]
        coordsX2.setValues(arrX2, 4, 1);
        coordsY2=DataArrayDouble.New();
        arrY2=[ -2., 2., 4., 8. ]
        coordsY2.setValues(arrY2, 4, 1);
        coordsZ2=DataArrayDouble.New();
        arrZ2=[ -3., 3., 6., 12.+1e-6 ]
        coordsZ2.setValues(arrZ2, 4, 1);
        mesh2.setCoords(coordsX2, coordsY2, coordsZ2);
        
        mesh3=MEDCouplingCMesh.New();
        coordsX3=DataArrayDouble.New();
        arrX3=[-1.]
        coordsX3.setValues(arrX3, 1, 1);
        coordsY3=DataArrayDouble.New();
        arrY3=[-2.]
        coordsY3.setValues(arrY3, 1, 1);
        coordsZ3=DataArrayDouble.New();
        arrZ3=[-3.]
        coordsZ3.setValues(arrZ3, 1, 1);
        mesh3.setCoords(coordsX3, coordsY3, coordsZ3);
        
        self.assertEqual(3, mesh1.getSpaceDimension());
        self.assertEqual(3, mesh1.getMeshDimension());
        
        self.assertTrue(not mesh1.isEqual(mesh2, 1e-12));
        self.assertTrue(not mesh2.isEqual(mesh1, 1e-12));
        self.assertTrue(not mesh2.isEqualWithoutConsideringStr(mesh1, 1e-12));
        self.assertTrue(mesh1.isEqual(mesh2, 1e-5));
        self.assertTrue(not mesh1.isEqual(mesh2, 1e-7));
        
        self.assertRaises(InterpKernelException, mesh3.checkConsistency, 1e-12);
        mesh1.checkConsistency(1e-12);
        self.assertEqual(NORM_HEXA8, mesh1.getTypeOfCell(1));
        
        self.assertEqual(NORM_HEXA8, mesh1.getAllGeoTypes()[0]);
        self.assertEqual(27, mesh1.getNumberOfCellsWithType(NORM_HEXA8));
        self.assertRaises(InterpKernelException, mesh1.getNumberOfCellsWithType, NORM_QUAD4);
        
        coo=mesh1.getCoordinatesOfNode(0);
        self.assertEqual(3, len(coo));
        self.assertAlmostEqual(-1., coo[0], 14);
        self.assertAlmostEqual(-2., coo[1], 14);
        self.assertAlmostEqual(-3., coo[2], 14);
        coo=mesh1.getCoordinatesOfNode(63);
        self.assertEqual(3, len(coo));
        self.assertAlmostEqual(4., coo[0], 14);
        self.assertAlmostEqual(8., coo[1], 14);
        self.assertAlmostEqual(12., coo[2], 14);
        
        a=str(mesh1)
        repr=mesh1.simpleRepr();
        repr=mesh1.advancedRepr();
        self.assertTrue("Cartesian" in repr);
        self.assertTrue("Number of components : 1" in repr);
        self.assertTrue("Number of tuples : 4" in repr);
        self.assertTrue("Z Array :" in repr);
        pass

    def testCMesh2(self):
        mesh1=MEDCouplingCMesh.New();
        coordsX1=DataArrayDouble.New();
        arrX1=[ -1., 1., 2., 4. ]
        coordsX1.setValues(arrX1, 4, 1);
        coordsY1=DataArrayDouble.New();
        arrY1=[ -2., 2., 4., 8. ]
        coordsY1.setValues(arrY1, 4, 1);
        coordsZ1=DataArrayDouble.New();
        arrZ1=[ -3., 3., 6., 12. ]
        coordsZ1.setValues(arrZ1, 4, 1);
        mesh1.setCoords(coordsX1, coordsY1, coordsZ1);
        
        dis=mesh1.getDistributionOfTypes();
        self.assertEqual(1, len(dis));
        self.assertEqual(NORM_HEXA8, dis[0][0]);
        self.assertEqual(27, dis[0][1]);
        self.assertEqual(-1, dis[0][2]);
        
        idsPerType=[]
        self.assertTrue(not mesh1.checkTypeConsistencyAndContig(dis, idsPerType));
        dis[0][0]=NORM_QUAD4;
        self.assertRaises(InterpKernelException, mesh1.checkTypeConsistencyAndContig, dis, idsPerType);
        dis[0][0]=NORM_HEXA8;
        dis[0][2]=0;
        ids=DataArrayInt.New();
        ids.alloc(10, 1);
        ids.fillWithValue(23);
        idsPerType=[ids];
        check=mesh1.checkTypeConsistencyAndContig(dis, idsPerType);
        self.assertTrue(check);
        self.assertTrue(check.isEqual(ids));
        
        code, idsInPflPerType, pfls=mesh1.splitProfilePerType(ids);
        self.assertEqual(1, len(code));
        self.assertEqual(NORM_HEXA8, code[0][0]);
        self.assertEqual(10, code[0][1]);
        self.assertEqual(0, code[0][2]);
        self.assertEqual(1, len(idsInPflPerType));
        self.assertEqual(1, len(pfls));
        self.assertTrue(idsInPflPerType[0].isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9])));
        self.assertTrue(pfls[0].isEqual(ids));
        
        cells1=[0, 1, 25, 26]
        partMesh1=mesh1.buildPart(cells1)
        self.assertTrue(isinstance(partMesh1,MEDCouplingMesh))
        self.assertEqual(4, partMesh1.getNumberOfCellsWithType(NORM_HEXA8));
        self.assertEqual(64, mesh1.getNumberOfNodes());
        self.assertEqual(64, partMesh1.getNumberOfNodes());
        
        cells2=[25, 26]
        partMesh2, arr1=mesh1.buildPartAndReduceNodes(cells2)
        self.assertTrue(isinstance(partMesh2,MEDCouplingCMesh))
        self.assertEqual(2,partMesh2.getNumberOfCellsWithType(NORM_HEXA8));
        self.assertEqual(12,partMesh2.getNumberOfNodes());
        
        cells3=[2, 3]
        partMesh3, arr2=partMesh1.buildPartAndReduceNodes(cells3)
        self.assertTrue(isinstance(partMesh3,MEDCouplingUMesh))
        self.assertEqual(2, partMesh3.getNumberOfCellsWithType(NORM_HEXA8));
        self.assertEqual(12, partMesh3.getNumberOfNodes());
        
        self.assertRaises(InterpKernelException, mesh1.simplexize, 0);
        self.assertRaises(InterpKernelException, mesh1.getMeasureFieldOnNode, True);
        
        #double bbox1[6];
        #double bbox2[6];
        bbox1=mesh1.getBoundingBox(); #[(-1.0, 4.0), (-2.0, 8.0), (-3.0, 12.0)]
        bbox2=partMesh1.getBoundingBox();
        self.assertTrue(bbox1==bbox2);
        bbox1=partMesh3.getBoundingBox();
        bbox2=partMesh2.getBoundingBox();
        self.assertTrue(bbox1==bbox2);
        
        self.assertRaises(InterpKernelException, mesh1.buildOrthogonalField);
        mesh2d=MEDCouplingCMesh.New();
        mesh2d.setCoords(coordsX1, coordsY1);
        f1=mesh2d.buildOrthogonalField();
        
        pass

    def testScale(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        pos=[0.2,0.2]
        mesh.scale(pos,0.5);
        expected1=[-0.05,-0.05, 0.2,-0.05, 0.45,-0.05, -0.05,0.2, 0.2,0.2, 0.45,0.2,
                   -0.05,0.45, 0.2,0.45, 0.45,0.45]
        val=mesh.getCoords().getValues();
        self.assertEqual(18,len(val))
        for i in range(18):
            self.assertTrue(abs(expected1[i]-val[i])<1e-12);
            pass
        pass

    def testTryToShareSameCoords(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        self.assertTrue(m1.getCoords().getHiddenCppPointer()!=m2.getCoords().getHiddenCppPointer());
        m1.tryToShareSameCoords(m2,1e-12);
        self.assertTrue(m1.getCoords().getHiddenCppPointer()==m2.getCoords().getHiddenCppPointer());
        m1.tryToShareSameCoords(m2,1e-12);
        self.assertTrue(m1.getCoords().getHiddenCppPointer()==m2.getCoords().getHiddenCppPointer());
        m2.tryToShareSameCoords(m1,1e-12);
        self.assertTrue(m1.getCoords().getHiddenCppPointer()==m2.getCoords().getHiddenCppPointer());
        #
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_2();
        self.assertTrue(m1.getCoords().getHiddenCppPointer()!=m2.getCoords().getHiddenCppPointer());
        m1.tryToShareSameCoords(m2,1e-12);
        self.assertTrue(m1.getCoords().getHiddenCppPointer()==m2.getCoords().getHiddenCppPointer());
        m1.tryToShareSameCoords(m2,1e-12);
        self.assertTrue(m1.getCoords().getHiddenCppPointer()==m2.getCoords().getHiddenCppPointer());
        m2.tryToShareSameCoords(m1,1e-12);
        self.assertTrue(m1.getCoords().getHiddenCppPointer()==m2.getCoords().getHiddenCppPointer());
        #
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DSourceMesh_1();
        self.assertTrue(m1.getCoords().getHiddenCppPointer()!=m2.getCoords().getHiddenCppPointer());
        self.assertRaises(InterpKernelException,m1.tryToShareSameCoords,m2,1e-12)
        pass

    def testFindNodeOnPlane(self):
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        pt=[300.,300.,0.]
        v=[0.,0.,2.]
        n=mesh.findNodesOnPlane(pt,v,1e-12);
        self.assertEqual(9,len(n));
        m3dSurf=mesh.buildFacePartOfMySelfNode(n,True);
        self.assertTrue(isinstance(m3dSurf,MEDCouplingUMesh))
        me=MEDCouplingMappedExtrudedMesh.New(mesh,m3dSurf,0);
        da=me.getMesh3DIds();
        self.assertEqual(8,me.getNumberOfCells());
        expected=[0,1,2,3,4,5,6,7]
        val=da.getValues();
        self.assertEqual(expected,list(val));
        #
        m3dSurf=mesh.buildFacePartOfMySelfNode(n,True);
        self.assertTrue(isinstance(m3dSurf,MEDCouplingUMesh))
        me=MEDCouplingMappedExtrudedMesh.New(mesh,m3dSurf,0);
        da=me.getMesh3DIds();
        self.assertEqual(8,me.getNumberOfCells());
        expected=[0,1,2,3,4,5,6,7]
        val=da.getValues();
        self.assertEqual(expected,list(val));
        pass

    def testRenumberCells(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        self.assertTrue(m.isEqual(m2,0));
        arr=[12,3,25,2,26]
        m.renumberCells(arr,True);
        self.assertTrue(not m.isEqual(m2,0));
        self.assertEqual(NORM_QUAD4,m.getTypeOfCell(0));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(1));
        self.assertEqual(NORM_QUAD4,m.getTypeOfCell(2));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(3));
        self.assertEqual(NORM_QUAD4,m.getTypeOfCell(4));
        arr2=[5,-1,-5,4,8]
        m.renumberCells(arr2,True);
        self.assertTrue(m.isEqual(m2,0));
        pass

    def testChangeSpaceDimension(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        #
        self.assertEqual(3,m1.getSpaceDimension());
        m1.changeSpaceDimension(2);
        self.assertEqual(2,m1.getSpaceDimension());
        m1.setName(m2.getName());
        self.assertTrue(m1.isEqual(m2,1e-12));
        m1.changeSpaceDimension(3);
        self.assertEqual(3,m1.getSpaceDimension());
        expected=[-0.3,-0.3,0., 0.2,-0.3,0., 0.7,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.7,0.2,0., -0.3,0.7,0., 0.2,0.7,0., 0.7,0.7,0.]
        val=m1.getCoords().getValues();
        for i in range(27):
            self.assertTrue(abs(expected[i]-val[i])<1e-14);
            pass
        pass

    def testGaussPointField1(self):
        _a=0.446948490915965;
        _b=0.091576213509771;
        _p1=0.11169079483905;
        _p2=0.0549758718227661;
        refCoo1=[ 0.,0., 1.,0., 0.,1. ]
        gsCoo1=[ 2*_b-1, 1-4*_b, 2*_b-1, 2.07*_b-1, 1-4*_b,
                 2*_b-1, 1-4*_a, 2*_a-1, 2*_a-1, 1-4*_a, 2*_a-1, 2*_a-1 ]
        wg1=[ 4*_p2, 4*_p2, 4*_p2, 4*_p1, 4*_p1, 4*_p1 ]
        _refCoo1=refCoo1
        _gsCoo1=gsCoo1
        _wg1=wg1
        #
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,NO_TIME);
        f.setMesh(m);
        self.assertEqual(5,f.getNumberOfMeshPlacesExpected());
        self.assertEqual(0,f.getNbOfGaussLocalization());
        f.setGaussLocalizationOnType(NORM_TRI3,_refCoo1,_gsCoo1,_wg1);
        f.setGaussLocalizationOnType(NORM_TRI3,_refCoo1,_gsCoo1,_wg1); # not a bug only to check that it works well
        self.assertRaises(InterpKernelException,f.setGaussLocalizationOnType,NORM_QUAD4,_refCoo1,_gsCoo1,_wg1)
        self.assertEqual(1,f.getNbOfGaussLocalization());
        refCoo2=[ 0.,0., 1.,0., 1.,1., 0.,1. ]
        _refCoo2=refCoo2
        _gsCoo1=_gsCoo1[0:4]
        _wg1=_wg1[0:2]
        f.setGaussLocalizationOnType(NORM_QUAD4,_refCoo2,_gsCoo1,_wg1);
        self.assertEqual(2,f.getNbOfGaussLocalization());
        array=DataArrayDouble.New();
        ptr=18*2*[None]
        for i in range(18 * 2):
            ptr[i]=float(i+1)
        array.setValues(ptr,18,2);
        ptr=array.getPointer();
        f.setArray(array);
        f.setName("MyFirstFieldOnGaussPoint");
        f.checkConsistencyLight();
        self.assertAlmostEqual(27.,f.getIJK(2,5,0),14);
        self.assertAlmostEqual(16.,f.getIJK(1,5,1),14);
        #
        f.clearGaussLocalizations();
        self.assertEqual(0,f.getNbOfGaussLocalization());
        self.assertRaises(InterpKernelException,f.checkConsistencyLight);
        ids1=[0,1,3,4]
        self.assertRaises(InterpKernelException,f.setGaussLocalizationOnCells,ids1,_refCoo2,_gsCoo1,_wg1);
        self.assertEqual(0,f.getNbOfGaussLocalization());
        ids2=[0,4]
        f.setGaussLocalizationOnCells(ids2,_refCoo2,_gsCoo1,_wg1);
        self.assertEqual(1,f.getNbOfGaussLocalization());
        self.assertEqual(0,f.getGaussLocalizationIdOfOneCell(0));
        self.assertRaises(InterpKernelException,f.getGaussLocalizationIdOfOneCell,1);
        ids3=[1,2]
        f.setGaussLocalizationOnCells(ids3,_refCoo1,_gsCoo1,_wg1);
        self.assertEqual(2,f.getNbOfGaussLocalization());
        self.assertEqual(0,f.getGaussLocalizationIdOfOneCell(0));
        self.assertEqual(1,f.getGaussLocalizationIdOfOneCell(1));
        self.assertEqual(1,f.getGaussLocalizationIdOfOneCell(2));
        self.assertRaises(InterpKernelException,f.checkConsistencyLight);#<- cell 3 has no localization
        ids4=[3]
        _gsCoo2=_gsCoo1;
        _wg2=_wg1;
        _gsCoo2[0]=0.8888777776666;
        _wg2[0]=0.1234567892377;
        f.setGaussLocalizationOnCells(ids4,_refCoo2,_gsCoo2,_wg2);
        self.assertEqual(3,f.getNbOfGaussLocalization());
        tmpIds=f.getCellIdsHavingGaussLocalization(0);
        self.assertEqual(ids2,list(tmpIds.getValues()));
        self.assertRaises(InterpKernelException,f.checkConsistencyLight);#<- it's always not ok because undelying array not with the good size.
        array2=f.getArray().subArray(0,10);
        f.setArray(array2);
        f.checkConsistencyLight();#<- here it is OK
        f2=f.clone(True);
        self.assertTrue(f.isEqual(f2,1e-14,1e-14));
        gl1=f2.getGaussLocalization(0);
        tmp=gl1.getGaussCoord(1,1);
        self.assertAlmostEqual(2.07*_b-1,tmp,14);
        gl1.setGaussCoord(1,1,0.07);
        self.assertTrue(not f.isEqual(f2,1e-14,1e-14));
        gl1.setGaussCoord(1,1,tmp);
        self.assertTrue(f.isEqual(f2,1e-14,1e-14));
        f2.checkConsistencyLight();
        pass

    def testGaussPointNEField1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_NE,NO_TIME);
        f.setMesh(m);
        self.assertEqual(5,f.getNumberOfMeshPlacesExpected());
        f.setName("MyFirstFieldOnNE");
        f.setDescription("MyDescriptionNE");
        array=DataArrayDouble.New();
        tmp=18*2*[None]
        for i in range(18 * 2):
            tmp[i]=float(i+7)
            pass
        array.setValues(tmp,18,2);
        ptr=array.getPointer();
        f.setArray(array);
        #
        f.checkConsistencyLight();
        f2=f.clone(True);
        self.assertTrue(f.isEqual(f2,1e-14,1e-14));
        self.assertAlmostEqual(21.,f.getIJK(2,0,0),14);
        self.assertAlmostEqual(18.,f.getIJK(1,1,1),14);
        pass

    def testCellOrientation1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        vec=[0.,0.,-1.]
        self.assertRaises(InterpKernelException,m.are2DCellsNotCorrectlyOriented,vec,False);
        m.changeSpaceDimension(3);
        res1=m.are2DCellsNotCorrectlyOriented(vec,False);
        self.assertTrue(len(res1)==0);
        vec[2]=1.;
        res1=m.are2DCellsNotCorrectlyOriented(vec,False);
        self.assertEqual(5,len(res1));
        #
        vec[2]=-1.;
        # connectivity inversion
        conn=m.getNodalConnectivity().getValues();
        tmp=conn[11];
        conn[11]=conn[12];
        conn[12]=tmp;
        m.getNodalConnectivity().setValues(conn,len(conn),1)
        res1=m.are2DCellsNotCorrectlyOriented(vec,False);
        self.assertEqual(1,len(res1));
        self.assertEqual(2,res1.getValues()[0]);
        m.orientCorrectly2DCells(vec,False);
        res1=m.are2DCellsNotCorrectlyOriented(vec,False);
        self.assertTrue(len(res1)==0);
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2.changeSpaceDimension(3);
        self.assertTrue(m.isEqual(m2,1e-12));
        pass

    def testCellOrientation2(self):
        m2,m1=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        res1=m2.arePolyhedronsNotCorrectlyOriented();
        self.assertEqual(6,len(res1));
        m2.orientCorrectlyPolyhedrons();
        res1=m2.arePolyhedronsNotCorrectlyOriented();
        self.assertTrue(len(res1)==0);
        m2.checkConsistencyLight();
        self.assertEqual(18,m2.getNumberOfCells());
        cellIds2=[0,6,12]
        m2.convertToPolyTypes(cellIds2);
        m2.orientCorrectlyPolyhedrons();
        res1=m2.arePolyhedronsNotCorrectlyOriented();
        self.assertTrue(len(res1)==0);
        f2=m2.getMeasureField(False);
        f2Ptr=f2.getArray().getValues();
        #Test to check global reverse in MEDCouplingUMesh::tryToCorrectPolyhedronOrientation
        m3=MEDCouplingDataForTest.build2DTargetMesh_1();
        vec=[0.,0.,1.]
        m3.changeSpaceDimension(3);
        ids2=[0,1,2,3,4]
        m3.convertToPolyTypes(ids2);
        m3.orientCorrectly2DCells(vec,False);
        m4=MEDCouplingDataForTest.buildCU1DMesh_U();
        m4.changeSpaceDimension(3);
        center=[0.,0.,0.]
        vector=[0.,1.,0.]
        m4.rotate(center,vector,-pi/2.);
        m5=m3.buildExtrudedMesh(m4,0);
        res1=m5.arePolyhedronsNotCorrectlyOriented();
        self.assertEqual(15,len(res1));
        m5.orientCorrectlyPolyhedrons();
        res1=m5.arePolyhedronsNotCorrectlyOriented();
        self.assertTrue(len(res1)==0);
        f3=m5.getMeasureField(False);
        self.assertEqual(15,f3.getArray().getNumberOfTuples());
        self.assertEqual(1,f3.getNumberOfComponents());
        f3Ptr=f3.getArray().getValues();
        expected1=[0.075,0.0375,0.0375,0.075,0.075, 0.1125,0.05625,0.05625,0.1125,0.1125, 0.0625,0.03125,0.03125,0.0625,0.0625];
        for i in range(15):
            self.assertTrue(abs(expected1[i]-f3Ptr[i])<1e-12);
            pass
        f4=m5.computeCellCenterOfMass();
        self.assertEqual(15,f4.getNumberOfTuples());
        self.assertEqual(3,f4.getNumberOfComponents());
        f4Ptr=f4.getValues();
        expected2=[-0.05,-0.05,0.15, 0.3666666666666667,-0.13333333333333333,0.15, 0.53333333333333333,0.033333333333333333,0.15, -0.05,0.45,0.15, 0.45,0.45,0.15,-0.05,-0.05,0.525, 0.3666666666666667,-0.13333333333333333,0.525, 0.53333333333333333,0.033333333333333333,0.525, -0.05,0.45,0.525, 0.45,0.45,0.525,-0.05,-0.05,0.875, 0.3666666666666667,-0.13333333333333333,0.875, 0.53333333333333333,0.033333333333333333,0.875, -0.05,0.45,0.875, 0.45,0.45,0.875];
        for i in range(45):
            self.assertTrue(abs(expected2[i]-f4Ptr[i])<1e-12);
            pass
        pass

    def testCellOrientation3(self):
        from cmath import rect  

        c = [rect(1.0, i * pi / 4.0) for i in range(8)]
        coords = [c[-1].real,c[-1].imag,  c[3].real,c[3].imag,
                   c[5].real,c[5].imag,  c[1].real,c[1].imag]
        connec = [0,1,2,3] 
        baseMesh = MEDCouplingUMesh.New("circle", 2)  
        baseMesh.allocateCells(1)
        meshCoords = DataArrayDouble.New(coords, 4, 2)
        baseMesh.setCoords(meshCoords)
        baseMesh.insertNextCell(NORM_QPOLYG, connec)  # a circle
        baseMesh.finishInsertingCells()  
        baseMesh.changeSpaceDimension(3)
        Oz = [0.0, 0.0, -1.0] 
        cell_lst = baseMesh.are2DCellsNotCorrectlyOriented(Oz, False)
        self.assertEqual(cell_lst.getNumberOfTuples(), 0)
        Oz[2] = 1.0
        cell_lst = baseMesh.are2DCellsNotCorrectlyOriented(Oz, False)
        self.assertEqual(cell_lst.getNumberOfTuples(), 1)

    def testPolyhedronBarycenter(self):
        connN=[0,3,2,1, -1, 4,5,6,7, -1, 0,4,7,3, -1, 3,7,6,2, -1, 2,6,5,1, -1, 1,5,4,0];
        coords=[0.,0.,0., 1.,0.,0., 1.,1.,0., 0.,1.,0., 0.,0.,1., 1.,0.,1., 1.,1.,1., 0.,1.,1., 0.5, 0.5, 0.5];
        meshN=MEDCouplingUMesh.New();
        meshN.setName("ForBary");
        meshN.setMeshDimension(3);
        meshN.allocateCells(4);
        meshN.insertNextCell(NORM_POLYHED,29,connN[0:29])
        meshN.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(coords,9,3);
        meshN.setCoords(myCoords);
        meshN.checkConsistencyLight();
        #
        res1=meshN.arePolyhedronsNotCorrectlyOriented();
        meshN.orientCorrectlyPolyhedrons();
        self.assertTrue(len(res1)==0);
        da=meshN.computeCellCenterOfMass();
        self.assertEqual(1,da.getNumberOfTuples());
        self.assertEqual(3,da.getNumberOfComponents());
        daPtr=da.getValues();
        ref=meshN.getCoords().getValues()[24:];
        for i in range(3):
            self.assertTrue(abs(ref[i]-daPtr[i])<1e-12);
            pass
        #
        center=[0.,0.,0.]
        vec=[0.,2.78,0.]
        da=meshN.computeCellCenterOfMass();
        daPtr=da.getValues();
        ref=meshN.getCoords().getValues()[24:];
        for i in range(3):
            self.assertTrue(abs(ref[i]-daPtr[i])<1e-12);
            pass
        #
        meshN.rotate(center,vec,pi/7.);
        meshN.translate(vec);
        da=meshN.computeCellCenterOfMass();
        daPtr=da.getValues();
        ref=meshN.getCoords().getValues()[24:];
        for i in range(3):
            self.assertTrue(abs(ref[i]-daPtr[i])<1e-12);
            pass
        #
        center2=[1.12,3.45,6.78]
        vec2=[4.5,9.3,2.8]
        meshN.rotate(center2,vec2,e);
        meshN.translate(vec2);
        da=meshN.computeCellCenterOfMass();
        daPtr=da.getValues();
        ref=meshN.getCoords().getValues()[24:];
        for i in range(3):
            self.assertTrue(abs(ref[i]-daPtr[i])<1e-10);
            pass
        pass

    def testNormL12Integ1D(self):
        m1=MEDCouplingDataForTest.build1DTargetMesh_3();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(m1);
        array=DataArrayDouble.New();
        arr=[-5.23,15.45,-25.56,6.67,-16.78,26.89,-7.91,17.23,-27.43,8.21,-18.63,28.72]
        array.setValues(arr,m1.getNumberOfCells(),3);
        f1.setArray(array);
        #
        f3=m1.computeCellCenterOfMass();
        self.assertEqual(4,f3.getNumberOfTuples());
        self.assertEqual(1,f3.getNumberOfComponents());
        expected9=[0.75,5.105,0.8,5.155]
        ptr=f3.getValues();
        for i in range(4):
            self.assertTrue(abs(expected9[i]-ptr[i])<1e-12);
            pass
        #
        f2=m1.getMeasureField(False);
        self.assertEqual(4,f2.getArray().getNumberOfTuples());
        self.assertEqual(1,f2.getNumberOfComponents());
        expected1=[0.5,0.21,-0.6,-0.31]
        ptr=f2.getArray().getValues();
        for i in range(4):
            self.assertTrue(abs(expected1[i]-ptr[i])<1e-12);
            pass
        expected2=[0.5,0.21,0.6,0.31]
        f2=m1.getMeasureField(True);
        ptr=f2.getArray().getValues();
        for i in range(4):
            self.assertTrue(abs(expected2[i]-ptr[i])<1e-12);
            pass
        #integral
        self.assertTrue(4,f1.getNumberOfTuples())
        res=f1.integral(False);
        self.assertTrue(3,len(res))
        expected3=[0.9866,-0.3615,0.4217]
        for i in range(3):
            self.assertTrue(abs(expected3[i]-res[i])<1e-12);
            pass
        self.assertTrue(abs(expected3[0]-f1.integral(0,False))<1e-12);
        self.assertTrue(abs(expected3[1]-f1.integral(1,False))<1e-12);
        self.assertTrue(abs(expected3[2]-f1.integral(2,False))<1e-12);
        res=f1.integral(True);
        expected4=[-3.4152,8.7639,-14.6879]
        for i in range(3):
            self.assertTrue(abs(expected4[i]-res[i])<1e-12);
            pass
        #normL1
        res=f1.normL1();
        self.assertTrue(3,len(res))
        expected5=[6.979506172839505, 16.89018518518518, 27.02969135802469]
        for i in range(3):
            self.assertTrue(abs(expected5[i]-res[i])<1e-12);
            pass
        self.assertTrue(abs(expected5[0]-f1.normL1(0))<1e-12);
        self.assertTrue(abs(expected5[1]-f1.normL1(1))<1e-12);
        self.assertTrue(abs(expected5[2]-f1.normL1(2))<1e-12);
        #normL2
        res=f1.normL2();
        self.assertTrue(3,len(res))
        expected7=[7.090910979452395, 16.9275542960123, 27.053271464160858]
        for i in range(3):
            self.assertTrue(abs(expected7[i]-res[i])<1e-9);
            pass
        self.assertTrue(abs(expected7[0]-f1.normL2(0))<1e-9);
        self.assertTrue(abs(expected7[1]-f1.normL2(1))<1e-9);
        self.assertTrue(abs(expected7[2]-f1.normL2(2))<1e-9);
        #buildMeasureField
        f4=f1.buildMeasureField(False);
        self.assertTrue(abs(-0.2-f4.accumulate(0))<1e-12);
        f4=f1.buildMeasureField(True);
        self.assertTrue(abs(1.62-f4.accumulate(0))<1e-12);
        # Testing with 2D Curve
        m1=MEDCouplingDataForTest.build2DCurveTargetMesh_3();
        f2=m1.getMeasureField(False);
        self.assertEqual(4,f2.getArray().getNumberOfTuples());
        self.assertEqual(1,f2.getNumberOfComponents());
        ptr=f2.getArray().getValues();
        for i in range(4):
            self.assertTrue(abs(sqrt(2.)*expected2[i]-ptr[i])<1e-12);
            pass
        f2=m1.getMeasureField(True);
        self.assertEqual(4,f2.getArray().getNumberOfTuples());
        self.assertEqual(1,f2.getNumberOfComponents());
        ptr=f2.getArray().getValues();
        for i in range(4):
            self.assertTrue(abs(expected2[i]*sqrt(2.)-ptr[i])<1e-12);
            pass
        #bary
        f3=m1.computeCellCenterOfMass();
        self.assertEqual(4,f3.getNumberOfTuples());
        self.assertEqual(2,f3.getNumberOfComponents());
        expected10=[0.75,0.75,5.105,5.105,0.8,0.8,5.155,5.155]
        ptr=f3.getValues();
        for i in range(8):
            self.assertTrue(abs(expected10[i]-ptr[i])<1e-12);
            pass
        #
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(m1);
        array=DataArrayDouble.New();
        array.setValues(arr,m1.getNumberOfCells(),3);
        f1.setArray(array);
        res=f1.integral(False);
        for i in range(3):
            self.assertTrue(abs(sqrt(2.)*expected4[i]-res[i])<1e-12);
            pass
        res=f1.integral(True);
        for i in range(3):
            self.assertTrue(abs(sqrt(2.)*expected4[i]-res[i])<1e-12);
            pass
        res=f1.normL1();
        for i in range(3):
            self.assertTrue(abs(expected5[i]-res[i])<1e-12);
            pass
        res=f1.normL2();
        for i in range(3):
            self.assertTrue(abs(expected7[i]-res[i])<1e-12);
            pass
        pass

    def testAreaBary2D(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_3();
        f1=m1.getMeasureField(False);
        self.assertEqual(10,f1.getArray().getNumberOfTuples());
        self.assertEqual(1,f1.getNumberOfComponents());
        expected1=[-0.5,-1,-1.5,-0.5,-1,  0.5,1,1.5,0.5,1]
        ptr=f1.getArray().getValues();
        for i in range(10):
            self.assertTrue(abs(expected1[i]-ptr[i])<1e-12);
            pass
        f1=m1.getMeasureField(True);
        ptr=f1.getArray().getValues();
        for i in range(10):
            self.assertTrue(abs(abs(expected1[i])-ptr[i])<1e-12);
            pass
        f2=m1.computeCellCenterOfMass();
        self.assertEqual(10,f2.getNumberOfTuples());
        self.assertEqual(2,f2.getNumberOfComponents());
        expected2=[0.5,0.3333333333333333,0.5,0.5,0.5,0.77777777777777777,0.5,0.3333333333333333,0.5,0.5,0.5,0.3333333333333333,0.5,0.5,0.5,0.77777777777777777,0.5,0.3333333333333333,0.5,0.5]
        ptr=f2.getValues();
        for i in range(20):
            self.assertTrue(abs(expected2[i]-ptr[i])<1e-12);
            pass
        m1.changeSpaceDimension(3);
        f1=m1.getMeasureField(False);
        self.assertEqual(10,f1.getArray().getNumberOfTuples());
        self.assertEqual(1,f1.getNumberOfComponents());
        ptr=f1.getArray().getValues();
        for i in range(10):
            self.assertTrue(abs(abs(expected1[i])-ptr[i])<1e-12);
            pass
        f2=m1.computeCellCenterOfMass();
        self.assertEqual(10,f2.getNumberOfTuples());
        self.assertEqual(3,f2.getNumberOfComponents());
        ptr=f2.getValues();
        expected3=[0.5,0.3333333333333333,0.,0.5,0.5,0.,0.5,0.77777777777777777,0.,0.5,0.3333333333333333,0.,0.5,0.5,0., 0.5,0.3333333333333333,0.,0.5,0.5,0.,0.5,0.77777777777777777,0.,0.5,0.3333333333333333,0.,0.5,0.5,0.]
        for i in range(30):
            self.assertTrue(abs(expected3[i]-ptr[i])<1e-12);
            pass
        pass

    def testAreaBary3D(self):
        coords=[ 0.241310763507 , 0.0504777305619 , 0.0682283524903 , 0.252501053866 , -0.0625176732937 , 0.137272639894 ,
                 0.152262663601 , 0.241816569527 , 0.133812556197 , 0.18047750211 , -0.0789949051358 , 0.339098173401 ,
                 0.151741971857 , 0.238885278571 , 0.137715037333 , 0.242532155481 , -0.0928169086456 , 0.0678043417367 ,
                 0.240941965335 , -0.015461491464 , 0.0617186345825 , 0.24127650112 , 0.0499427876717 , 0.0679634099148 ,
                 -0.145828917428 , 0.206291632565 , 0.0310071927543 , 0.0125651775307 , 0.266262085828 , 0.105228430543 ,
                 -0.0994066533286 , 0.233224271238 , 0.0572213839567 , -0.0951345338317 , 0.234819509426 , 0.0592126284538 ,
                 0.136580574205 , -0.205486212579 , 0.0572866072014 , 0.0637270784978 , -0.168886355238 , 0.446614057077 ,
                 0.041337157151 , -0.213402568198 , 0.372407095999 , 0.0411601970268 , -0.202387875756 , 0.411334979491 ,
                 -0.108355701857 , 0.193636239335 , 0.204886756738 , 0.00639779029829 , 0.155296981517 , 0.252585892979 ,
                 0.0262473111702 , -0.112919732543 , 0.424286639249 ,-0.224103052733 , -0.139430015438 , -0.0122352295701 ,
                -0.0312760589481 , -0.274272003594 , 0.0323959636568 , -0.166663422532 , -0.217754445175 , 0.00392109070364 ,
                 -0.30586619777 , -0.0475168041091 , -0.0144585228182 , -0.280881480586 , 0.135571293538 , 0.00623923647986 ,
                 -0.25548538234 , 0.156819217766 , 0.0645277879769 , -0.131567009284 , 0.184133752309 , 0.206021802753 ,
                 -0.196204010965 , 0.151602971681 , 0.212974777736 , -0.183713879463 , 0.0802946639531 , 0.260115662599 ,
                 -0.244241178767 , -0.0738873389604 , 0.144590565817 , -0.155804057829 , -0.164892720025 , 0.210613950558 ,
                 -0.170950800428 , -0.215099334026 , 0.00610122860092 , -0.30552634869 , -0.0490020791904 , -0.0132786533145 ,
                 0.271831011884 , 0.15105657296 , 0.0230534827908 , 0.281919192283 , 0.0898544306288 , -0.0625201489143 ,
                 0.260240727276 , -0.0120688706637 , -0.0532316588626 , 0.244947737722 , 0.0197984684293 , 0.0309341209233 ,
                 0.23439631578 , 0.229825279875 , 0.0508520585381 , 0.160921316875 , 0.265078502128 , 0.121716560626 ,
                 -0.315088694175 , 0.0747700471918 , -0.245836615071 , -0.327728781776 , 0.0857114674649 , -0.239431905957 ,
                 -0.308385460634 , 0.145142997084 , -0.149886828433 , 0.0488236045164 , 0.309462801914 , 0.0849169148265 ,
                -0.0244964803395 , 0.33145611751 , -0.0476415818061 , 0.0060567994229 , 0.32418412014 , 0.0367779543812 ,
                 -0.0950221448063 , 0.236675326003 , 0.0572594453983 , 0.248723023186 , 0.0886648784791 , -0.176629430538 ,
                 0.116796984 , 0.256596599567 , -0.292863523603 , 0.118024552914 , 0.229154257843 , -0.34233232501 ,
                 0.217507892549 , -0.0417822335742 , -0.176771782888 , -0.224429321304 , 0.0125595300114 , -0.362064725588 ,
                 0.0937301100955 , -0.0500824832657 , -0.299713548444 , -0.244162220397 , 0.0383853931293 , -0.389856984411 ,
                 -0.0281989366102 , 0.097392811563 , -0.458244577284 , -0.385010847162 , 0.10122766194 , -0.140052859922 ,
                 -0.377936358012 , 0.110875172128 , -0.176207095463 , 0.244483045556 , -0.0991073977045 , 0.0575134372934 ,
                0.262605120167 , -0.100243191645 , -0.0495620806935 , 0.240306880972 , -0.136153701579 , -0.114745281696 ,
                 0.215763176129 , -0.0836766059189 , -0.183249640616 , 0.237870396603 , -0.132449578286 , -0.121598854639 ,
                 -0.0637683083097 , -0.27921020214 , -0.149112321992 , -0.0856211014977 , -0.2973233473 , -0.0446878139589 ,
                 0.104675342288 , -0.0625908305324 , -0.290346256534 , 0.0248264249186 , -0.247797708548 , -0.165830884019 ,
                 0.0719302438309 , -0.178468260473 , -0.211432157345 , 0.142871843159 , -0.208769948542 , 0.0454101128246 ,
                 0.167803379307 , -0.207851396623 , -0.088802726124 , 0.12868717152 , -0.230920439715 , 0.00760508389036 ,
                 -0.0372812069535 , -0.286740286332 , 0.00963701291166 ]
        
        connN = [ #polyhedron 0
            0 , 1 , 3 , 4 , 2 , -1 , 1 , 5 , 6 , 7 , 0 , -1 , 0 , 7 , 8 , 10 , 11 , 9 , 2 , -1 , 1 , 5 , 12 , 14 , 15 , 13 , 3 , -1 , 16 , 9 , 2 , 4 , 17 , -1
            , 4 , 3 , 13 , 18 , 17 , -1 , 5 , 6 , 19 , 21 , 20 , 12 , -1 , 6 , 7 , 8 , 23 , 22 , 19 , -1 , 23 , 24 , 10 , 8 , -1 , 25 , 11 , 9 , 16 , -1
            , 24 , 26 , 25 , 11 , 10 , -1 , 12 , 14 , 20 , -1 , 27 , 28 , 29 , 15 , 13 , 18 , -1 , 14 , 15 , 29 , 30 , 21 , 20 , -1 , 26 , 27 , 18 , 17 , 16 , 25 , -1
            , 22 , 19 , 21 , 30 , 31 , -1 , 22 , 31 , 28 , 27 , 26 , 24 , 23 , -1 , 31 , 30 , 29 , 28,
            # polyhedron 1
            0 , 7 , 8 , 10 , 11 , 9 , 2 , -1 , 32 , 0 , 7 , 35 , 34 , 33 , -1 , 32 , 0 , 2 , 37 , 36 , -1 , 35 , 7 , 8 , 40 , 39 , 38 , -1
            , 2 , 37 , 41 , 9 , -1 , 40 , 8 , 10 , 44 , 43 , 42 , -1 , 41 , 9 , 11 , 44 , 43 , -1 , 44 , 11 , 10 , -1 , 32 , 33 , 45 , 47 , 46 , 36 , -1
            , 33 , 34 , 48 , 45 , -1 , 35 , 34 , 48 , 50 , 49 , 38 , -1 , 41 , 43 , 42 , 46 , 36 , 37 , -1 , 38 , 39 , 51 , 49 , -1
            , 39 , 40 , 42 , 46 , 47 , 52 , 51 , -1 , 45 , 47 , 52 , 50 , 48 , -1 , 52 , 51 , 49 , 50,
            # polyhedron 2
            6 , 7 , 8 , 23 , 22 , 19 , -1 , 6 , 35 , 7 , -1 , 6 , 35 , 38 , 19 , -1 , 35 , 7 , 8 , 40 , 39 , 38 , -1 , 53 , 22 , 19 , 38 , 39 , 54 , -1
            , 23 , 53 , 54 , 40 , 8 , -1 , 53 , 22 , 23 , -1 , 39 , 54 , 40,
            # polyhedron 3
            35 , 34 , 48 , 50 , 49 , 38 , -1 , 6 , 35 , 34 , 56 , 55 , 5 , -1 , 6 , 35 , 38 , 19 , -1 , 34 , 56 , 57 , 59 , 58 , 48 , -1
            , 60 , 61 , 21 , 19 , 38 , 49 , -1 , 62 , 50 , 48 , 58 , -1 , 60 , 63 , 64 , 62 , 50 , 49 , -1 , 5 , 6 , 19 , 21 , 20 , 12 , -1
            , 55 , 5 , 12 , 65 , -1 , 66 , 67 , 65 , 55 , 56 , 57 , -1 , 63 , 66 , 57 , 59 , 64 , -1 , 64 , 62 , 58 , 59 , -1
            , 60 , 63 , 66 , 67 , 68 , 61 , -1 , 61 , 68 , 20 , 21 , -1 , 67 , 68 , 20 , 12 , 65]
        
        barys = [ -0.0165220465527 , -0.0190922868195 , 0.158882733414 ,
                  0.0287618656076 , 0.135874379934 , -0.14601588119 ,
                  -0.147128055553 , 0.0465995097041 , -0.049391174453 ,
                  -0.00142506732317 , -0.0996953090351 , -0.115159183132 ]
        meshN=MEDCouplingUMesh.New();
        meshN.setName("ForBary");
        meshN.setMeshDimension(3);
        meshN.allocateCells(4);
        meshN.insertNextCell(NORM_POLYHED,113,connN);
        meshN.insertNextCell(NORM_POLYHED,99,connN[113:]);
        meshN.insertNextCell(NORM_POLYHED,43,connN[212:]);
        meshN.insertNextCell(NORM_POLYHED,92,connN[255:]);
        meshN.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(coords,69,3);
        meshN.setCoords(myCoords);
        meshN.checkConsistencyLight();
        res1=meshN.arePolyhedronsNotCorrectlyOriented();
        meshN.orientCorrectlyPolyhedrons();
        res1=meshN.arePolyhedronsNotCorrectlyOriented();
        self.assertTrue(len(res1)==0);
        #
        da=meshN.computeCellCenterOfMass();
        self.assertEqual(4,da.getNumberOfTuples());
        self.assertEqual(3,da.getNumberOfComponents());
        daPtr=da.getValues();
        for i in range(12):
            self.assertTrue(abs(barys[i]-daPtr[i])<1e-12);
            pass
        pass

    def testRenumberCellsForFields(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f.setMesh(m);
        arr=DataArrayDouble.New();
        nbOfCells=m.getNumberOfCells();
        values1=[7.,107.,10007.,8.,108.,10008.,9.,109.,10009.,10.,110.,10010.,11.,111.,10011.]
        arr.setValues(values1,nbOfCells,3);
        f.setArray(arr);
        renumber1=[3,1,0,4,2]
        loc=[-0.05,-0.05, 0.55,-0.25, 0.55,0.15, -0.05,0.45, 0.45,0.45]
        for j in range(5):
            res=f.getValueOn(loc[2*j:2*j+2]);
            for i in range(3):
                self.assertTrue(abs(values1[i+3*j]-res[i])<1e-12);
                pass
            pass
        f.renumberCells(renumber1,False);
        ptr=f.getArray().getValues();
        expected1=[9.,109.,10009.,8.,108.,10008.,11.,111.,10011.,7.,107.,10007.,10.,110.,10010.]
        for i in range(15):
            self.assertTrue(abs(expected1[i]-ptr[i])<1e-12);
            pass
        #check that fields remains the same geometrically
        for j in range(5):
            res=f.getValueOn(loc[2*j:2*(j+1)]);
            for i in range(3):
                self.assertTrue(abs(values1[i+3*j]-res[i])<1e-12);
                pass
            pass
        #On gauss
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,NO_TIME);
        f.setMesh(m);
        _a=0.446948490915965;
        _b=0.091576213509771;
        _p1=0.11169079483905;
        _p2=0.0549758718227661;
        refCoo1=[ 0.,0., 1.,0., 0.,1. ]
        gsCoo1=[ 2*_b-1, 1-4*_b, 2*_b-1, 2.07*_b-1, 1-4*_b, 2*_b-1, 1-4*_a, 2*_a-1, 2*_a-1, 1-4*_a, 2*_a-1, 2*_a-1 ];
        wg1=[ 4*_p2, 4*_p2, 4*_p2, 4*_p1, 4*_p1, 4*_p1 ]
        _refCoo1=refCoo1[0:6];
        _gsCoo1=gsCoo1[0:12];
        _wg1=wg1[0:6];
        f.setGaussLocalizationOnType(NORM_TRI3,_refCoo1,_gsCoo1,_wg1);
        refCoo2=[ 0.,0., 1.,0., 1.,1., 0.,1. ]
        _refCoo2=refCoo2[0:8];
        _gsCoo1=_gsCoo1[0:4]
        _wg1=_wg1[0:2]
        f.setGaussLocalizationOnType(NORM_QUAD4,_refCoo2,_gsCoo1,_wg1);
        arr=DataArrayDouble.New();
        values2=[1.,1001.,2.,1002., 11.,1011.,12.,1012.,13.,1013.,14.,1014.,15.,1015.,16.,1016., 21.,1021.,22.,1022.,23.,1023.,24.,1024.,25.,1025.,26.,1026., 31.,1031.,32.,1032., 41.,1041.,42.,1042.]
        arr.setValues(values2,18,2);
        f.setArray(arr);
        f.checkConsistencyLight();
        fCpy=f.clone(True);
        self.assertTrue(f.isEqual(fCpy,1e-12,1e-12));
        f.renumberCells(renumber1,False);
        self.assertTrue(not f.isEqual(fCpy,1e-12,1e-12));
        expected2=[21.,1021.,22.,1022.,23.,1023.,24.,1024.,25.,1025.,26.,1026., 11.,1011.,12.,1012.,13.,1013.,14.,1014.,15.,1015.,16.,1016., 41.,1041.,42.,1042., 1.,1001.,2.,1002., 31.,1031.,32.,1032.]
        ptr=f.getArray().getValues();
        for i in range(36):
            self.assertTrue(abs(expected2[i]-ptr[i])<1e-12);
            pass
        renumber2=[2,1,4,0,3]
        f.renumberCells(renumber2,False);
        self.assertTrue(f.isEqual(fCpy,1e-12,1e-12));
        #GaussNE
        f=MEDCouplingFieldDouble.New(ON_GAUSS_NE,NO_TIME);
        f.setMesh(m);
        arr=DataArrayDouble.New();
        values3=[1.,1001.,2.,1002.,3.,1003.,4.,1004., 11.,1011.,12.,1012.,13.,1013., 21.,1021.,22.,1022.,23.,1023., 31.,1031.,32.,1032.,33.,1033.,34.,1034., 41.,1041.,42.,1042.,43.,1043.,44.,1044.]
        arr.setValues(values3,18,2);
        f.setArray(arr);
        f.checkConsistencyLight();
        fCpy=f.clone(True);
        self.assertTrue(f.isEqual(fCpy,1e-12,1e-12));
        f.renumberCells(renumber1,False);
        self.assertTrue(not f.isEqual(fCpy,1e-12,1e-12));
        expected3=[21.,1021.,22.,1022.,23.,1023.,11.,1011.,12.,1012.,13.,1013.,41.,1041.,42.,1042.,43.,1043.,44.,1044.,1.,1001.,2.,1002.,3.,1003.,4.,1004.,31.,1031.,32.,1032.,33.,1033.,34.,1034.]
        ptr=f.getArray().getValues();
        for i in range(36):
            self.assertTrue(abs(expected3[i]-ptr[i])<1e-12);
            pass
        f.renumberCells(renumber2,False);#perform reverse operation of renumbering to check that the resulting field is equal.
        self.assertTrue(f.isEqual(fCpy,1e-12,1e-12));
        #
        pass

    def testRenumberNodesForFields(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=MEDCouplingFieldDouble.New(ON_NODES,NO_TIME);
        f.setMesh(m);
        self.assertEqual(9,f.getNumberOfMeshPlacesExpected());
        arr=DataArrayDouble.New();
        nbOfNodes=m.getNumberOfNodes();
        values1=[7.,107.,10007.,8.,108.,10008.,9.,109.,10009.,10.,110.,10010.,11.,111.,10011.,12.,112.,10012.,13.,113.,10013.,14.,114.,10014.,15.,115.,10015.]
        arr.setValues(values1,nbOfNodes,3);
        f.setArray(arr);
        f.checkConsistencyLight();
        renumber1=[0,4,1,3,5,2,6,7,8]
        loc=[0.5432,-0.2432, 0.5478,0.1528]
        expected1=[9.0272, 109.0272, 10009.0272, 11.4124,111.4124,10011.4124]
        for j in range(2):
            res=f.getValueOn(loc[2*j:2*j+2]);
            for i in range(3):
                self.assertTrue(abs(expected1[i+3*j]-res[i])<1e-12);
                pass
            pass
        fCpy=f.clone(True);
        self.assertTrue(f.isEqual(fCpy,1e-12,1e-12));
        f.renumberNodes(renumber1);
        self.assertTrue(not f.isEqual(fCpy,1e-12,1e-12));
        for j in range(2):
            res=f.getValueOn(loc[2*j:2*j+2]);
            for i in range(3):
                self.assertTrue(abs(expected1[i+3*j]-res[i])<1e-12);
                pass
            pass
        expected2=[7.,107.,10007.,9.,109.,10009.,12.,112.,10012.,10.,110.,10010.,8.,108.,10008.,11.,111.,10011.,13.,113.,10013.,14.,114.,10014.,15.,115.,10015.]
        for i in range(27):
            self.assertTrue(abs(expected2[i]-f.getArray().getValues()[i])<1e-12);
            pass
        renumber2=[0,2,5,3,1,4,6,7,8]
        f.renumberNodes(renumber2);
        self.assertTrue(f.isEqual(fCpy,1e-12,1e-12));
        pass

    def testConvertQuadraticCellsToLinear(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_3();
        mesh.checkConsistencyLight();
        types=mesh.getAllGeoTypes();
        types.sort()
        self.assertEqual(5,len(types));
        expected1 = [NORM_POLYGON, NORM_TRI3, NORM_QUAD4, NORM_TRI6, NORM_QUAD8]
        expected1.sort()
        self.assertEqual(expected1,types);
        self.assertTrue(mesh.isPresenceOfQuadratic());
        self.assertEqual(62,mesh.getNodalConnectivityArrayLen());
        f1=mesh.getMeasureField(False);
        #
        mesh.convertQuadraticCellsToLinear();
        self.assertTrue(not mesh.isPresenceOfQuadratic());
        #
        mesh.checkConsistencyLight();
        f2=mesh.getMeasureField(False);
        self.assertTrue(f1.getArray().isEqual(f2.getArray(),1e-12));
        self.assertEqual(48,mesh.getNodalConnectivityArrayLen());
        types2=mesh.getAllGeoTypes();
        types2.sort()
        self.assertEqual(3,len(types2));
        expected2=[NORM_POLYGON, NORM_TRI3, NORM_QUAD4]
        expected2.sort()
        self.assertEqual(expected2,types2);
        pass

    def testCheckGeoEquivalWith(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_3();
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_3();
        #First test mesh1
        cellCor,nodeCor=mesh1.checkGeoEquivalWith(mesh1,0,1e-12);#deepEqual
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        cellCor,nodeCor=mesh1.checkGeoEquivalWith(mesh1,1,1e-12);#fastEqual
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        cellCor,nodeCor=mesh1.checkGeoEquivalWith(mesh1,10,1e-12);#deepEqual with geo permutations
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        #Second test mesh1 and mesh2 are 2 different meshes instance
        cellCor,nodeCor=mesh1.checkGeoEquivalWith(mesh2,0,1e-12);#deepEqual
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        cellCor,nodeCor=mesh1.checkGeoEquivalWith(mesh2,1,1e-12);#fastEqual
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        cellCor,nodeCor=mesh1.checkGeoEquivalWith(mesh2,10,1e-12);#deepEqual with geo permutations
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        #Third test : cell permutation by keeping the first the middle and the last as it is.
        renum=[0,2,1,3,4,5,6,8,7,9]
        mesh2.renumberCells(renum,False);
        self.assertRaises(InterpKernelException,mesh1.checkGeoEquivalWith,mesh2,0,1e-12);#deepEqual fails
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        cellCor,nodeCor=mesh1.checkGeoEquivalWith(mesh2,1,1e-12);#fastEqual do not see anything
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        cellCor,nodeCor=mesh1.checkGeoEquivalWith(mesh2,10,1e-12);#deepEqual with geo permutations
        self.assertTrue(cellCor);
        self.assertEqual(10,cellCor.getNumberOfTuples());
        self.assertEqual(1,cellCor.getNumberOfComponents());
        self.assertEqual(renum,list(cellCor.getValues()))
        self.assertTrue(nodeCor==None);
        cellCor=0;
        self.assertTrue(nodeCor==None);
        a,b=mesh1.checkDeepEquivalWith(mesh2,0,1e-12);
        self.assertEqual(renum,list(a.getValues()))
        self.assertTrue(b==None);
        mesh2.setCoords(mesh1.getCoords())
        a=mesh1.checkDeepEquivalOnSameNodesWith(mesh2,0,1e-12);
        self.assertEqual(renum,list(a.getValues()))
        #4th test : cell and node permutation by keeping the first the middle and the last as it is.
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_3();
        renum2=[0,2,1,3,4,5,6,8,7,9,10]
        mesh2.renumberCells(renum,False);
        mesh2.renumberNodes(renum2,11);
        cellCor=None
        nodeCor=None
        self.assertRaises(InterpKernelException,mesh1.checkGeoEquivalWith,mesh2,0,1e-12);#deepEqual fails
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        cellCor,nodeCor=mesh1.checkGeoEquivalWith(mesh2,1,1e-12);#fastEqual do not see anything
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        cellCor,nodeCor=mesh1.checkGeoEquivalWith(mesh2,10,1e-12);#deepEqual with geo permutations
        self.assertTrue(cellCor);
        self.assertEqual(10,cellCor.getNumberOfTuples());
        self.assertEqual(1,cellCor.getNumberOfComponents());
        self.assertEqual(renum,list(cellCor.getValues()))
        self.assertTrue(nodeCor);
        self.assertEqual(11,nodeCor.getNumberOfTuples());
        self.assertEqual(1,nodeCor.getNumberOfComponents());
        self.assertEqual(renum2,list(nodeCor.getValues()))
        cellCor=0;
        nodeCor=0;
        #5th test : modification of the last cell to check fastCheck detection.
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_3();
        renum3=[0,2,1,3,4,5,6,8,9,7]
        mesh2.renumberCells(renum3,False);
        mesh2.renumberNodes(renum2,11);
        cellCor=None
        nodeCor=None
        self.assertRaises(InterpKernelException,mesh1.checkGeoEquivalWith,mesh2,0,1e-12)
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        self.assertRaises(InterpKernelException,mesh1.checkGeoEquivalWith,mesh2,1,1e-12)
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        cellCor,nodeCor=mesh2.checkGeoEquivalWith(mesh1,10,1e-12);#deepEqual with geo permutations
        self.assertTrue(cellCor!=None);
        self.assertEqual(10,cellCor.getNumberOfTuples());
        self.assertEqual(1,cellCor.getNumberOfComponents());
        self.assertEqual(renum3,list(cellCor.getValues()))
        self.assertTrue(nodeCor!=None);
        self.assertEqual(11,nodeCor.getNumberOfTuples());
        self.assertEqual(1,nodeCor.getNumberOfComponents());
        self.assertEqual(renum2,list(nodeCor.getValues()));
        pass

    def testCheckGeoEquivalWith2(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_4();
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_1();
        cellCor,nodeCor=mesh1.checkGeoEquivalWith(mesh2,10,1e-12);
        self.assertEqual(None,cellCor);
        self.assertNotEqual(None,nodeCor);
        expected1=[0, 1, 3, 4, 5, 6, 7, 8, 9]
        for i in range(9):
            self.assertEqual(expected1[i],nodeCor.getIJ(i,0));
            pass
        pass
      
    def testSwig2CheckDeepEquivalWith1(self):
        eps = 1.0e-8
        mcart = MEDCouplingCMesh()
        mcart.setCoordsAt(0, DataArrayDouble([0.0,1.5,2.0]))
        mcart.setCoordsAt(1, DataArrayDouble([1.0,2.5,3.0,4.0]))
        m = mcart.buildUnstructured()
        m2 = m[1:m.getNumberOfCells()]
        self.assertRaises(InterpKernelException, m.checkDeepEquivalWith, m2, 0, eps)
        self.assertRaises(InterpKernelException, m.checkDeepEquivalWith, m2, 1, eps)
        self.assertRaises(InterpKernelException, m.checkDeepEquivalWith, m2, 2, eps)
        pass

    def testSwig2CheckDeepEquivalWith2(self):
        eps = 1.0e-8
        m = MEDCouplingUMesh("tst", 2)
        m.setCoords(DataArrayDouble([], 0,2))
        m.setConnectivity(DataArrayInt([]), DataArrayInt([0]))
        m2 = m.deepCopy()
        m.checkDeepEquivalWith(m2, 0, eps)  # Should not raise!
        pass

    def testCopyTinyStringsFromOnFields(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        nbOfCells=m.getNumberOfCells();
        f=MEDCouplingFieldDouble.New(ON_CELLS,LINEAR_TIME);
        f.setMesh(m);
        self.assertEqual(5,f.getNumberOfMeshPlacesExpected());
        f.setName("a");
        f.setDescription("b");
        a1=DataArrayDouble.New();
        a1.alloc(nbOfCells,2);
        a1.fillWithZero();
        a1.setInfoOnComponent(0,"c");
        a1.setInfoOnComponent(1,"d");
        a2=a1.deepCopy();
        a2.setInfoOnComponent(0,"e");
        a2.setInfoOnComponent(1,"f");
        f.setArray(a1);
        f.setEndArray(a2);
        f.setEndTime(3.,3,4);
        m.setName("g");
        m.getCoords().setInfoOnComponent(0,"h");
        m.getCoords().setInfoOnComponent(1,"i");
        m.getCoords().setInfoOnComponent(2,"j");
        #
        f.checkConsistencyLight();
        f2=f.clone(True);
        self.assertTrue(f2.isEqual(f,1e-12,1e-12));
        f2.setName("smth");
        self.assertTrue(not f2.isEqual(f,1e-12,1e-12));
        f2.copyTinyStringsFrom(f);
        self.assertTrue(f2.isEqual(f,1e-12,1e-12));
        f2.setDescription("GGG");
        self.assertTrue(not f2.isEqual(f,1e-12,1e-12));
        f2.copyTinyStringsFrom(f);
        self.assertTrue(f2.isEqual(f,1e-12,1e-12));
        f2.getArray().setInfoOnComponent(0,"mmmm");
        self.assertTrue(not f2.isEqual(f,1e-12,1e-12));
        f2.copyTinyStringsFrom(f);
        self.assertTrue(f2.isEqual(f,1e-12,1e-12));
        f2.getEndArray().setInfoOnComponent(1,"mmmm");
        self.assertTrue(not f2.isEqual(f,1e-12,1e-12));
        f2.copyTinyStringsFrom(f);
        self.assertTrue(f2.isEqual(f,1e-12,1e-12));
        m2=m.clone(True);
        self.assertTrue(m2.isEqual(m,1e-12));
        m2.setName("123");
        self.assertTrue(not m2.isEqual(m,1e-12));
        m2.copyTinyStringsFrom(m);
        self.assertTrue(m2.isEqual(m,1e-12));
        m2.getCoords().setInfoOnComponent(1,"eee");
        self.assertTrue(not m2.isEqual(m,1e-12));
        m2.copyTinyStringsFrom(m);
        self.assertTrue(m2.isEqual(m,1e-12));
        pass

    def testTryToShareSameCoordsPermute(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        #self.assertTrue(m.getCoords()!=m2.getCoords());
        m.tryToShareSameCoordsPermute(m2,1e-12);
        #self.assertTrue(m.getCoords()==m2.getCoords());
        self.assertTrue(m2.isEqual(m,1e-12));
        renum1=[1,2,0,5,8,7,4,3,6]
        r1=DataArrayInt.New()
        r1.setValues(renum1,len(renum1),1)
        m.renumberNodes(r1,9);
        #self.assertTrue(m.getCoords()!=m2.getCoords());
        self.assertTrue(not m2.isEqual(m,1e-12));
        m.tryToShareSameCoordsPermute(m2,1e-12);
        #self.assertTrue(m.getCoords()==m2.getCoords());
        self.assertTrue(m2.isEqual(m,1e-12));
        pass

    def testTryToShareSameCoordsPermute2(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_4();
        targetCoords=[-0.3,-0.3, 0.2,-0.3, -0.3,0.2, 0.2,0.2 ]
        targetConn=[0,2,3,1]
        m2=MEDCouplingUMesh.New();
        m2.setMeshDimension(2);
        m2.allocateCells(1);
        m2.insertNextCell(NORM_QUAD4,targetConn[0:4])
        m2.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,4,2);
        m2.setCoords(myCoords);
        m2.checkConsistencyLight();
        m1.checkConsistencyLight();
        #
        expected1=[0.25,0.125,0.125,0.25,0.25]
        f1=m1.getMeasureField(False);
        f2=m2.getMeasureField(False);
        self.assertEqual(5,f1.getArray().getNumberOfTuples());
        self.assertEqual(1,f2.getArray().getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(expected1[i],f1.getIJ(i,0),12);
            pass
        self.assertAlmostEqual(expected1[0],f2.getIJ(0,0),12);
        self.assertRaises(InterpKernelException,m1.tryToShareSameCoordsPermute,m2,1e-12);# <- here in this order the sharing is impossible.
        # Let's go for deeper test of tryToShareSameCoordsPermute
        m2.tryToShareSameCoordsPermute(m1,1e-12);
        f1=m1.getMeasureField(False);
        f2=m2.getMeasureField(False);
        self.assertEqual(5,f1.getArray().getNumberOfTuples());
        self.assertEqual(1,f2.getArray().getNumberOfTuples());
        for i in range(5):
            self.assertAlmostEqual(expected1[i],f1.getIJ(i,0),12);
            pass
        self.assertAlmostEqual(expected1[0],f2.getIJ(0,0),12);
        pass

    def testChangeUnderlyingMesh1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_3();
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_3();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr=[7., 107., 8., 108., 9., 109., 10., 110., 11., 111., 12., 112., 13., 113., 14., 114., 15., 115., 16., 116.]
        array.setValues(arr,mesh1.getNumberOfCells(),2);
        f1.setArray(array);
        #
        renum=[0,2,1,3,4,5,6,8,7,9]
        mesh2.renumberCells(renum,False);
        #self.assertTrue(f1.getMesh()==mesh1);
        f1.changeUnderlyingMesh(mesh1,10,1e-12);# nothing done only to check that nothing done.
        #self.assertTrue(f1.getMesh()==mesh1);
        f1.changeUnderlyingMesh(mesh2,10,1e-12);
        #self.assertTrue(f1.getMesh()==mesh2);
        expected1=[7.,107.,9.,109.,8.,108.,10.,110.,11.,111.,12.,112.,13.,113.,15.,115.,14.,114.,16.,116.]
        for i in range(20):
            self.assertAlmostEqual(expected1[i],f1.getArray().getIJ(0,i),12);
            pass
        #
        f1=MEDCouplingFieldDouble.New(ON_NODES,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr2=[7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.,13.,113.,14.,114.,15.,115.,16.,116.,17.,117.]
        array.setValues(arr2,mesh1.getNumberOfNodes(),2);
        f1.setArray(array);
        #
        renum2=[0,2,10,3,4,5,6,8,7,9,1]
        mesh2.renumberNodes(renum2,11);
        #self.assertTrue(f1.getMesh()==mesh1);
        f1.changeUnderlyingMesh(mesh2,10,1e-12);
        #self.assertTrue(f1.getMesh()==mesh2);
        expected2=[7.,107.,17.,117.,8.,108.,10.,110.,11.,111.,12.,112.,13.,113.,15.,115.,14.,114.,16.,116.,9.,109.]
        for i in range(22):
            self.assertAlmostEqual(expected2[i],f1.getArray().getIJ(0,i),12);
            pass
        pass

    def testGetMaxValue1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        nbOfCells=m.getNumberOfCells();
        f=MEDCouplingFieldDouble.New(ON_CELLS,LINEAR_TIME);
        f.setMesh(m);
        a1=DataArrayDouble.New();
        val1=[3.,4.,5.,6.,7.]
        a1.setValues(val1,nbOfCells,1);
        a2=DataArrayDouble.New();
        val2=[0.,1.,2.,8.,7.]
        a2.setValues(val2,nbOfCells,1);
        f.setArray(a1);
        f.setEndArray(a2);
        f.setEndTime(3.,3,4);
        f.checkConsistencyLight();
        #
        self.assertAlmostEqual(8.,f.getMaxValue(),14);
        self.assertAlmostEqual(0.,f.getMinValue(),14);
        self.assertAlmostEqual(5.,f.getAverageValue(),14);
        self.assertAlmostEqual(5.125,f.getWeightedAverageValue(0,True),14);
        a1.setIJ(0,2,9.5);
        self.assertAlmostEqual(9.5,f.getMaxValue(),14);
        self.assertAlmostEqual(0.,f.getMinValue(),14);
        a2.setIJ(0,0,9.);
        self.assertAlmostEqual(9.5,f.getMaxValue(),14);
        self.assertAlmostEqual(1.,f.getMinValue(),14);
        pass

    def testSubstractInPlaceDM1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_3();
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_3();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr=[7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.,13.,113.,14.,114.,15.,115.,16.,116.]
        array.setValues(arr,mesh1.getNumberOfCells(),2);
        f1.setArray(array);
        #
        self.assertEqual(10,f1.getNumberOfTuples());
        self.assertEqual(2,f1.getNumberOfComponents());
        self.assertEqual(20,f1.getNumberOfValues());
        #
        renum=[0,2,3,1,4,5,6,8,7,9]
        mesh2.renumberCells(renum,False);
        #
        f2=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f2.setMesh(mesh2);
        array=DataArrayDouble.New();
        arr2=[7.1,107.1,10.1,110.1,8.1,108.1,9.1,109.1,11.1,111.1,12.1,112.1,13.1,113.1,15.1,115.1,14.1,114.1,16.1,116.1]
        array.setValues(arr2,mesh2.getNumberOfCells(),2);
        f2.setArray(array);
        #
        f1.substractInPlaceDM(f2,10,1e-12);
        f1.applyFunc(1,"abs(x+y+0.2)");
        self.assertAlmostEqual(0.,f1.getMaxValue(),13);
        pass

    def testDotCrossProduct1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_3();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setTime(2.3,5,6);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[7.,107.,207.,8.,108.,208.,9.,109.,209.,10.,110.,210.,11.,111.,211.,12.,112.,212.,13.,113.,213.,14.,114.,214.,15.,115.,215.,16.,116.,216.]
        array.setValues(arr1,mesh1.getNumberOfCells(),3);
        f1.setArray(array);
        f2=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f2.setTime(7.8,4,5);
        f2.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr2=[1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.]
        array.setValues(arr2,mesh1.getNumberOfCells(),3);
        f2.setArray(array);
        #
        f3=f1.dot(f2);
        expected1=[842.,1820.,2816.,3830.,4862.,5912.,6980.,8066.,9170.,10292.]
        for i in range(10):
            self.assertAlmostEqual(expected1[i],f3.getIJ(i,0),9);
            pass
        #
        f4=f1.crossProduct(f2);
        expected2=[-93., 186., -93., -392., 784., -392., -691., 1382., -691., -990., 1980., -990., -1289., 2578., -1289., -1588., 3176., -1588., -1887., 3774., -1887., -2186., 4372., -2186., -2485., 4970., -2485., -2784., 5568., -2784.]
        for i in range(30):
            self.assertAlmostEqual(expected2[i],f4.getIJ(0,i),9);
            pass
        pass

    pass

if __name__ == '__main__':
    unittest.main()
