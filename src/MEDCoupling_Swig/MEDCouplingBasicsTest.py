#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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
import rlcompleter,readline # this line has to be here, to ensure a usability of MEDCoupling/MEDLoader. B4 removing it please notify to anthony.geay@cea.fr

class MEDCouplingBasicsTest(unittest.TestCase):
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
        arr2=arr1.substr(3);
        self.assertEqual(4,arr2.getNumberOfTuples());
        self.assertEqual(2,arr2.getNumberOfComponents());
        self.assertEqual(arr1Ref[6:],list(arr2.getValues()));
        arr3=arr1.substr(2,5);
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
        for i in xrange(14):
            self.assertTrue(abs(arr4Ref[i]-tmp[i])<1e-14);
            pass
        arr5=arr4.substr(3);
        self.assertEqual(4,arr5.getNumberOfTuples());
        self.assertEqual(2,arr5.getNumberOfComponents());
        tmp=arr5.getValues()
        for i in xrange(8):
            self.assertTrue(abs(arr4Ref[6+i]-tmp[i])<1e-14);
            pass
        arr6=arr4.substr(2,5);
        self.assertEqual(3,arr6.getNumberOfTuples());
        self.assertEqual(2,arr6.getNumberOfComponents());
        tmp=arr6.getValues()
        for i in xrange(6):
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
        mesh.checkCoherency();
        self.assertTrue(mesh.getAllGeoTypes()==[4])
        myFalseConn=DataArrayInt.New()
        myFalseConn.setValues(tab4,6,4)
        self.assertTrue(myFalseConn.getIJ(1,1)==3)
        #
        field=MEDCouplingFieldDouble.New(ON_CELLS)
        field.setMesh(mesh)
        field.setNature(Integral)
        myCoords=DataArrayDouble.New()
        sampleTab=[]
        for i in range(nbOfCells*9):
            sampleTab.append(float(i))
        myCoords.setValues(sampleTab,nbOfCells,9);
        field.setArray(myCoords)
        self.assertTrue(3==mesh.getSpaceDimension())
        field.checkCoherency()
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
        self.assertRaises(InterpKernelException,targetMesh.checkCoherency);
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
        meshM1D.checkCoherency();
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
        fieldOnCells.checkCoherency();
        pass
    
    def testDeepCopy(self):
        array=DataArrayDouble.New();
        array.setValues(5*3*[7.],5,3);
        self.assertEqual(array.getIJ(3,2),7.);
        array2=array.deepCpy();
        self.assertEqual(array2.getIJ(3,2),7.)
        #
        array3=DataArrayInt.New();
        array3.setValues(5*3*[17],5,3);
        self.assertEqual(array3.getIJ(3,2),17);
        array4=array3.deepCpy();
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
        mesh.checkCoherency();
        self.assertEqual(5,mesh.getNumberOfCells());
        self.assertEqual(23,mesh.getNodalConnectivity().getNumberOfTuples());
        expected1=[4, 0, 3, 4, 1, 5, 1, 4, 2, 3, 4, 5, 2, 5, 6, 7, 4, 3, 4, 7, 8, 5, 4]
        self.assertEqual(expected1,list(mesh.getNodalConnectivity().getValues()));
        #
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        mesh.convertToPolyTypes(elts);
        mesh.checkCoherency();
        self.assertEqual(8,mesh.getNumberOfCells());
        self.assertEqual(114,mesh.getNodalConnectivity().getNumberOfTuples());
        mesh.convertToPolyTypes(elts);
        mesh.checkCoherency();
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
        mesh2.checkCoherency();
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
        mesh.checkCoherency();
        #
        desc=DataArrayInt.New();
        descIndx=DataArrayInt.New();
        revDesc=DataArrayInt.New();
        revDescIndx=DataArrayInt.New();
        #
        mesh2=mesh.buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
        mesh2.checkCoherency();
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
        mesh2.checkCoherency();
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
        mesh.checkCoherency();
        desc=DataArrayInt.New();
        descIndx=DataArrayInt.New();
        revDesc=DataArrayInt.New();
        revDescIndx=DataArrayInt.New();
        mesh2=mesh.buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
        mesh2.checkCoherency();
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
        arr2=arr.deepCpy();
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
        field.setNature(Integral);
        field.setNature(ConservativeVolumic);
        field.setNature(IntegralGlobConstraint);
        field=MEDCouplingFieldDouble.New(ON_NODES,NO_TIME);
        field.setNature(ConservativeVolumic);
        self.assertRaises(InterpKernelException,field.setNature,Integral);
        self.assertRaises(InterpKernelException,field.setNature,IntegralGlobConstraint);
        pass

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
        ext=MEDCouplingExtrudedMesh.New(mesh3D,mesh2D,1);
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
        m4=MEDCouplingExtrudedMesh.New(m3,m1,0);
        self.assertEqual(15,m4.getNumberOfCells());
        self.assertEqual(5,m4.getMesh2D().getNumberOfCells());
        self.assertEqual(3,m4.getMesh1D().getNumberOfCells());
        m3DIds=m4.getMesh3DIds().getValues();
        self.assertEqual(range(15),list(m3DIds));
        #some random in cells to check that extrusion alg find it correctly
        expected1=[1,3,2,0,6,5,7,10,11,8,12,9,14,13,4]
        m3.renumberCells(expected1,False);
        m4=MEDCouplingExtrudedMesh.New(m3,m1,0);
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
        m4=MEDCouplingExtrudedMesh.New(m3,m1,0);
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
        m3=m1.buildExtrudedMesh(m2,0);
        expected1=[1,3,2,0,6,5,7,10,11,8,12,9,14,13,4]
        rexpected1=[3, 0, 2, 1, 14, 5, 4, 6, 9, 11, 7, 8, 10, 13, 12]
        m3.renumberCells(expected1,False);
        m4=MEDCouplingExtrudedMesh.New(m3,m1,0);
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
        for i in xrange(15):
            self.assertAlmostEqual(expected2[rexpected1[i]],arrPtr[i],16);
            pass
        m5=m4.build3DUnstructuredMesh();
        self.assertTrue(m5.isEqual(m3,1e-12));
        f=m5.getMeasureField(True);
        f.setMesh(m4)
        self.assertTrue(isinstance(f.getMesh(),MEDCouplingExtrudedMesh))
        arr=f.getArray();
        arrPtr=arr.getValues();
        for i in xrange(15):
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
        o2nExp1=range(27)
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
        m3.checkCoherency();
        m4=MEDCouplingDataForTest.build2DTargetMeshMerged_1();
        self.assertTrue(m3.isEqual(m4,1.e-12));
        da,isMerged,newNbOfNodes=m3.mergeNodes(1.e-12);
        self.assertEqual(11,m3.getNumberOfNodes());
        self.assertTrue(isMerged);
        pass

    def testMergeMeshOnSameCoords1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells=range(5);
        m2.convertToPolyTypes(cells);
        m1.tryToShareSameCoords(m2,1e-12);
        m3=MEDCouplingDataForTest.build2DTargetMesh_1();
        m3.tryToShareSameCoords(m2,1e-12);
        meshes=[m1,m2,m3]
        m4=MEDCouplingUMesh.MergeUMeshesOnSameCoords(meshes);
        m4.checkCoherency();
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
        f3.checkCoherency();
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
        for i in xrange(7):
            self.assertTrue(abs(values[i]-tmp[i])<1e-12)
            pass
        pass

    def testFillFromAnalytic(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        m.setTime(3.4,5,6); m.setTimeUnit("us");
        f1=m.fillFromAnalytic(ON_CELLS,1,"x+y");
        self.assertAlmostEqual(3.4,f1.getTime()[0],12) ; self.assertEqual(5,f1.getTime()[1]) ; self.assertEqual(6,f1.getTime()[2])
        self.assertEqual("us",f1.getTimeUnit())
        f1.checkCoherency();                    
        self.assertEqual(f1.getTypeOfField(),ON_CELLS);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        values1=[-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values1),len(tmp))
        for i in xrange(len(tmp)):
            self.assertTrue(abs(tmp[i]-values1[i])<1.e-12)
            pass
        #
        f1=m.fillFromAnalytic(ON_NODES,1,"x+y");
        f1.checkCoherency();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values2=[-0.6,-0.1,0.4,-0.1,0.4,0.9,0.4,0.9,1.4]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values2),len(tmp))
        for i in xrange(len(tmp)):
            self.assertTrue(abs(tmp[i]-values2[i])<1.e-12)
            pass
        #
        f1=m.fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+(2*(x+y))*JVec");
        f1.checkCoherency();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(2,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values3=[-0.6,-1.2,-0.1,-0.2,0.4,0.8,-0.1,-0.2,0.4,0.8,0.9,1.8,0.4,0.8,0.9,1.8,1.4,2.8]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values3),len(tmp))
        for i in xrange(len(tmp)):
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
        f1.checkCoherency();
        self.assertEqual(f1.getTypeOfField(),ON_CELLS);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        values1=[-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values1),len(tmp))
        for i in xrange(len(values1)):
            self.assertTrue(abs(values1[i]-tmp[i])<1.e-12);
            pass
        #
        f1=m.fillFromAnalytic(ON_NODES,1,"y+2*x");
        f1.checkCoherency();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values2=[-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values2),len(tmp))
        for i in xrange(len(values2)):
            self.assertTrue(abs(values2[i]-tmp[i])<1.e-12);
            pass
        f1=m.fillFromAnalytic(ON_NODES,1,"2.*x+y");
        f1.checkCoherency();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        tmp=f1.getArray().getValues();
        values2Bis=[-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1]
        self.assertEqual(len(values2Bis),len(tmp))
        for i in xrange(len(values2Bis)):
            self.assertTrue(abs(values2Bis[i]-tmp[i])<1.e-12);
            pass
        #
        f1=m.fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+2*(x+y)*JVec");
        f1.checkCoherency();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(2,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values3=[-0.6,-1.2,-0.1,-0.2,0.4,0.8,-0.1,-0.2,0.4,0.8,0.9,1.8,0.4,0.8,0.9,1.8,1.4,2.8]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values3),len(tmp))
        for i in xrange(len(values3)):
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
        f1.checkCoherency();
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
        for i in xrange(len(tmp)):
            self.assertTrue(abs(tmp[i]-values1[i])<1.e-12)
            pass
        pass

    def testApplyFunc2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=m.fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+2*(x+y)*JVec");
        f1.checkCoherency();
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
        for i in xrange(len(tmp)):
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
        for i in xrange(len(tmp)):
            self.assertTrue(abs(tmp[i]-values1[i])<1.e-12)
            pass
        pass

    def testOperationsOnFields(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=m.fillFromAnalytic(ON_NODES,1,"x+y");
        f2=m.fillFromAnalytic(ON_NODES,1,"x+y");
        f1.checkCoherency();
        f2.checkCoherency();
        f3=f1+f2;
        f3.checkCoherency();
        self.assertEqual(f3.getTypeOfField(),ON_NODES);
        self.assertEqual(f3.getTimeDiscretization(),ONE_TIME);
        values1=[-1.2,-0.2,0.8,-0.2,0.8,1.8,0.8,1.8,2.8]
        tmp=f3.getArray().getValues();
        self.assertEqual(len(values1),len(tmp))
        for i in xrange(len(tmp)):
            self.assertTrue(abs(tmp[i]-values1[i])<1.e-12)
            pass
        #
        f3=f1*f2;
        f3.checkCoherency();
        self.assertEqual(f3.getTypeOfField(),ON_NODES);
        self.assertEqual(f3.getTimeDiscretization(),ONE_TIME);
        values2=[0.36,0.01,0.16,0.01,0.16,0.81,0.16,0.81,1.96]
        tmp=f3.getArray().getValues();
        self.assertEqual(len(values2),len(tmp))
        for i in xrange(len(tmp)):
            self.assertTrue(abs(tmp[i]-values2[i])<1.e-12)
            pass
        #
        f3=f1+f2;
        f4=f1-f3;
        f4.checkCoherency();
        self.assertEqual(f4.getTypeOfField(),ON_NODES);
        self.assertEqual(f4.getTimeDiscretization(),ONE_TIME);
        values3=[0.6,0.1,-0.4,0.1,-0.4,-0.9,-0.4,-0.9,-1.4]
        tmp=f4.getArray().getValues();
        self.assertEqual(len(values3),len(tmp))
        for i in xrange(len(tmp)):
            self.assertTrue(abs(tmp[i]-values3[i])<1.e-12)
            pass
        #
        f3=f1+f2;
        f4=f3/f2;
        f4.checkCoherency();
        self.assertEqual(f4.getTypeOfField(),ON_NODES);
        self.assertEqual(f4.getTimeDiscretization(),ONE_TIME);
        tmp=f4.getArray().getValues();
        for i in xrange(len(tmp)):
            self.assertTrue(abs(tmp[i]-2.)<1.e-12)
            pass
        #
        f4=f2.buildNewTimeReprFromThis(NO_TIME,False);
        f4.checkCoherency();
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
        for i in xrange(len(tmp)):
            self.assertTrue(abs(tmp[i]-values4[i])<1.e-12)
            pass
        #
        f4=f2.buildNewTimeReprFromThis(NO_TIME,True);
        f4.checkCoherency();
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
        for i in xrange(len(tmp)):
            self.assertTrue(abs(tmp[i]-values5[i])<1.e-12)
            pass
        pass

    def testOperationsOnFields2(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m.setTime(3.4,5,6); m.setTimeUnit("us");
        f1=m.fillFromAnalytic(ON_NODES,1,"x+y+z");
        f2=m.fillFromAnalytic(ON_NODES,1,"a*a+b+c*c");
        f3=f1/f2;
        f3.checkCoherency();
        self.assertEqual(f3.getTypeOfField(),ON_NODES);
        self.assertEqual(f3.getTimeDiscretization(),ONE_TIME);
        expected1=[-2.4999999999999991, 1.2162162162162162, 0.77868852459016391,
                   0.7407407407407407, 1.129032258064516, 0.81632653061224492,
                   0.86538461538461531, 1.0919540229885056, 0.84302325581395343]
        self.assertEqual(1,f3.getNumberOfComponents());
        self.assertEqual(9,f3.getNumberOfTuples());
        val=f3.getArray().getValues();
        for i in xrange(9):
            self.assertTrue(abs(expected1[i]-val[i])<1.e-12);
        #
        f1=m.buildOrthogonalField();
        self.assertAlmostEqual(3.4,f1.getTime()[0],12) ; self.assertEqual(5,f1.getTime()[1]) ; self.assertEqual(6,f1.getTime()[2])
        self.assertEqual("us",f1.getTimeUnit())
        f2=m.fillFromAnalytic(ON_CELLS,1,"x");
        f3=f1*f2;
        expected2=[-0.035355339059327376,0.,0.035355339059327376, 0.2592724864350674,0.,-0.2592724864350674, 0.37712361663282529,0.,-0.37712361663282529, -0.035355339059327376,0.,0.035355339059327376, 0.31819805153394637,0.,-0.31819805153394637]
        val=f3.getArray().getValues();
        for i in xrange(15):
            self.assertTrue(abs(expected2[i]-val[i])<1.e-12);
            pass
        #
        f3=f2*f1;
        val=f3.getArray().getValues();
        for i in xrange(15):
            self.assertTrue(abs(expected2[i]-val[i])<1.e-12);
            pass
        pass

    def testOperationsOnFields3(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        f1=m.fillFromAnalytic(ON_NODES,1,"x+y+z");
        f2=m.fillFromAnalytic(ON_NODES,1,"a*a+b+c*c");
        f1/=f2
        f1.checkCoherency();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        expected1=[-2.4999999999999991, 1.2162162162162162, 0.77868852459016391,
                   0.7407407407407407, 1.129032258064516, 0.81632653061224492,
                   0.86538461538461531, 1.0919540229885056, 0.84302325581395343]
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        val=f1.getArray().getValues();
        for i in xrange(9):
            self.assertTrue(abs(expected1[i]-val[i])<1.e-12);
            pass
        #
        f1=m.buildOrthogonalField();
        f2=m.fillFromAnalytic(ON_CELLS,1,"x");
        f1*=f2
        expected2=[-0.035355339059327376,0.,0.035355339059327376, 0.2592724864350674,0.,-0.2592724864350674, 0.37712361663282529,0.,-0.37712361663282529, -0.035355339059327376,0.,0.035355339059327376, 0.31819805153394637,0.,-0.31819805153394637]
        val=f1.getArray().getValues();
        for i in xrange(15):
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
        f1.checkCoherency();
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
        self.assertRaises(InterpKernelException,f2.checkCoherency)
        array2=DataArrayDouble.New();
        array2.setValues(arr2,nbOfCells,3);
        f2.setEndArray(array2);
        f2.checkCoherency();
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
        for i in xrange(3):
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
        for i in xrange(3):
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
        for i in xrange(15):
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
        for i in xrange(nbOfCells):
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
        for i in xrange(nbOfNodes):
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
        for i in xrange(64):
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
        for i in xrange(27):
          self.assertAlmostEqual(expected2[i], val[i], 12);
        #res=fieldOnCells.getValueOnPos(1,2,1);
        #self.assertAlmostEqual(6.,res,12);
        #
        meshDeepCopy=mesh.deepCpy();
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
        
        self.assertRaises(InterpKernelException, mesh3.checkCoherency1, 1e-12);
        mesh1.checkCoherency2(1e-12);
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
        for i in xrange(18):
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
        me=MEDCouplingExtrudedMesh.New(mesh,m3dSurf,0);
        da=me.getMesh3DIds();
        self.assertEqual(8,me.getNumberOfCells());
        expected=[0,1,2,3,4,5,6,7]
        val=da.getValues();
        self.assertEqual(expected,list(val));
        #
        m3dSurf=mesh.buildFacePartOfMySelfNode(n,True);
        self.assertTrue(isinstance(m3dSurf,MEDCouplingUMesh))
        me=MEDCouplingExtrudedMesh.New(mesh,m3dSurf,0);
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
        for i in xrange(27):
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
        for i in xrange(18*2):
            ptr[i]=float(i+1)
        array.setValues(ptr,18,2);
        ptr=array.getPointer();
        f.setArray(array);
        f.setName("MyFirstFieldOnGaussPoint");
        f.checkCoherency();
        self.assertAlmostEqual(27.,f.getIJK(2,5,0),14);
        self.assertAlmostEqual(16.,f.getIJK(1,5,1),14);
        #
        f.clearGaussLocalizations();
        self.assertEqual(0,f.getNbOfGaussLocalization());
        self.assertRaises(InterpKernelException,f.checkCoherency);
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
        self.assertRaises(InterpKernelException,f.checkCoherency);#<- cell 3 has no localization
        ids4=[3]
        _gsCoo2=_gsCoo1;
        _wg2=_wg1;
        _gsCoo2[0]=0.8888777776666;
        _wg2[0]=0.1234567892377;
        f.setGaussLocalizationOnCells(ids4,_refCoo2,_gsCoo2,_wg2);
        self.assertEqual(3,f.getNbOfGaussLocalization());
        tmpIds=f.getCellIdsHavingGaussLocalization(0);
        self.assertEqual(ids2,list(tmpIds.getValues()));
        self.assertRaises(InterpKernelException,f.checkCoherency);#<- it's always not ok because undelying array not with the good size.
        array2=f.getArray().substr(0,10);
        f.setArray(array2);
        f.checkCoherency();#<- here it is OK
        f2=f.clone(True);
        self.assertTrue(f.isEqual(f2,1e-14,1e-14));
        gl1=f2.getGaussLocalization(0);
        tmp=gl1.getGaussCoord(1,1);
        self.assertAlmostEqual(2.07*_b-1,tmp,14);
        gl1.setGaussCoord(1,1,0.07);
        self.assertTrue(not f.isEqual(f2,1e-14,1e-14));
        gl1.setGaussCoord(1,1,tmp);
        self.assertTrue(f.isEqual(f2,1e-14,1e-14));
        f2.checkCoherency();
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
        for i in xrange(18*2):
            tmp[i]=float(i+7)
            pass
        array.setValues(tmp,18,2);
        ptr=array.getPointer();
        f.setArray(array);
        #
        f.checkCoherency();
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
        m2.checkCoherency();
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
        for i in xrange(15):
            self.assertTrue(abs(expected1[i]-f3Ptr[i])<1e-12);
            pass
        f4=m5.getBarycenterAndOwner();
        self.assertEqual(15,f4.getNumberOfTuples());
        self.assertEqual(3,f4.getNumberOfComponents());
        f4Ptr=f4.getValues();
        expected2=[-0.05,-0.05,0.15, 0.3666666666666667,-0.13333333333333333,0.15, 0.53333333333333333,0.033333333333333333,0.15, -0.05,0.45,0.15, 0.45,0.45,0.15,-0.05,-0.05,0.525, 0.3666666666666667,-0.13333333333333333,0.525, 0.53333333333333333,0.033333333333333333,0.525, -0.05,0.45,0.525, 0.45,0.45,0.525,-0.05,-0.05,0.875, 0.3666666666666667,-0.13333333333333333,0.875, 0.53333333333333333,0.033333333333333333,0.875, -0.05,0.45,0.875, 0.45,0.45,0.875];
        for i in xrange(45):
            self.assertTrue(abs(expected2[i]-f4Ptr[i])<1e-12);
            pass
        pass

    def testCellOrientation3(self):
        from cmath import rect  

        c = [rect(1.0, i*pi/4.0) for i in range(8)]
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
        meshN.checkCoherency();
        #
        res1=meshN.arePolyhedronsNotCorrectlyOriented();
        meshN.orientCorrectlyPolyhedrons();
        self.assertTrue(len(res1)==0);
        da=meshN.getBarycenterAndOwner();
        self.assertEqual(1,da.getNumberOfTuples());
        self.assertEqual(3,da.getNumberOfComponents());
        daPtr=da.getValues();
        ref=meshN.getCoords().getValues()[24:];
        for i in xrange(3):
            self.assertTrue(abs(ref[i]-daPtr[i])<1e-12);
            pass
        #
        center=[0.,0.,0.]
        vec=[0.,2.78,0.]
        da=meshN.getBarycenterAndOwner();
        daPtr=da.getValues();
        ref=meshN.getCoords().getValues()[24:];
        for i in xrange(3):
            self.assertTrue(abs(ref[i]-daPtr[i])<1e-12);
            pass
        #
        meshN.rotate(center,vec,pi/7.);
        meshN.translate(vec);
        da=meshN.getBarycenterAndOwner();
        daPtr=da.getValues();
        ref=meshN.getCoords().getValues()[24:];
        for i in xrange(3):
            self.assertTrue(abs(ref[i]-daPtr[i])<1e-12);
            pass
        #
        center2=[1.12,3.45,6.78]
        vec2=[4.5,9.3,2.8]
        meshN.rotate(center2,vec2,e);
        meshN.translate(vec2);
        da=meshN.getBarycenterAndOwner();
        daPtr=da.getValues();
        ref=meshN.getCoords().getValues()[24:];
        for i in xrange(3):
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
        f3=m1.getBarycenterAndOwner();
        self.assertEqual(4,f3.getNumberOfTuples());
        self.assertEqual(1,f3.getNumberOfComponents());
        expected9=[0.75,5.105,0.8,5.155]
        ptr=f3.getValues();
        for i in xrange(4):
            self.assertTrue(abs(expected9[i]-ptr[i])<1e-12);
            pass
        #
        f2=m1.getMeasureField(False);
        self.assertEqual(4,f2.getArray().getNumberOfTuples());
        self.assertEqual(1,f2.getNumberOfComponents());
        expected1=[0.5,0.21,-0.6,-0.31]
        ptr=f2.getArray().getValues();
        for i in xrange(4):
            self.assertTrue(abs(expected1[i]-ptr[i])<1e-12);
            pass
        expected2=[0.5,0.21,0.6,0.31]
        f2=m1.getMeasureField(True);
        ptr=f2.getArray().getValues();
        for i in xrange(4):
            self.assertTrue(abs(expected2[i]-ptr[i])<1e-12);
            pass
        #integral
        self.assertTrue(4,f1.getNumberOfTuples())
        res=f1.integral(False);
        self.assertTrue(3,len(res))
        expected3=[0.9866,-0.3615,0.4217]
        for i in xrange(3):
            self.assertTrue(abs(expected3[i]-res[i])<1e-12);
            pass
        self.assertTrue(abs(expected3[0]-f1.integral(0,False))<1e-12);
        self.assertTrue(abs(expected3[1]-f1.integral(1,False))<1e-12);
        self.assertTrue(abs(expected3[2]-f1.integral(2,False))<1e-12);
        res=f1.integral(True);
        expected4=[-3.4152,8.7639,-14.6879]
        for i in xrange(3):
            self.assertTrue(abs(expected4[i]-res[i])<1e-12);
            pass
        #normL1
        res=f1.normL1();
        self.assertTrue(3,len(res))
        expected5=[6.979506172839505, 16.89018518518518, 27.02969135802469]
        for i in xrange(3):
            self.assertTrue(abs(expected5[i]-res[i])<1e-12);
            pass
        self.assertTrue(abs(expected5[0]-f1.normL1(0))<1e-12);
        self.assertTrue(abs(expected5[1]-f1.normL1(1))<1e-12);
        self.assertTrue(abs(expected5[2]-f1.normL1(2))<1e-12);
        #normL2
        res=f1.normL2();
        self.assertTrue(3,len(res))
        expected7=[7.090910979452395, 16.9275542960123, 27.053271464160858]
        for i in xrange(3):
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
        for i in xrange(4):
            self.assertTrue(abs(sqrt(2.)*expected2[i]-ptr[i])<1e-12);
            pass
        f2=m1.getMeasureField(True);
        self.assertEqual(4,f2.getArray().getNumberOfTuples());
        self.assertEqual(1,f2.getNumberOfComponents());
        ptr=f2.getArray().getValues();
        for i in xrange(4):
            self.assertTrue(abs(expected2[i]*sqrt(2.)-ptr[i])<1e-12);
            pass
        #bary
        f3=m1.getBarycenterAndOwner();
        self.assertEqual(4,f3.getNumberOfTuples());
        self.assertEqual(2,f3.getNumberOfComponents());
        expected10=[0.75,0.75,5.105,5.105,0.8,0.8,5.155,5.155]
        ptr=f3.getValues();
        for i in xrange(8):
            self.assertTrue(abs(expected10[i]-ptr[i])<1e-12);
            pass
        #
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(m1);
        array=DataArrayDouble.New();
        array.setValues(arr,m1.getNumberOfCells(),3);
        f1.setArray(array);
        res=f1.integral(False);
        for i in xrange(3):
            self.assertTrue(abs(sqrt(2.)*expected4[i]-res[i])<1e-12);
            pass
        res=f1.integral(True);
        for i in xrange(3):
            self.assertTrue(abs(sqrt(2.)*expected4[i]-res[i])<1e-12);
            pass
        res=f1.normL1();
        for i in xrange(3):
            self.assertTrue(abs(expected5[i]-res[i])<1e-12);
            pass
        res=f1.normL2();
        for i in xrange(3):
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
        for i in xrange(10):
            self.assertTrue(abs(expected1[i]-ptr[i])<1e-12);
            pass
        f1=m1.getMeasureField(True);
        ptr=f1.getArray().getValues();
        for i in xrange(10):
            self.assertTrue(abs(abs(expected1[i])-ptr[i])<1e-12);
            pass
        f2=m1.getBarycenterAndOwner();
        self.assertEqual(10,f2.getNumberOfTuples());
        self.assertEqual(2,f2.getNumberOfComponents());
        expected2=[0.5,0.3333333333333333,0.5,0.5,0.5,0.77777777777777777,0.5,0.3333333333333333,0.5,0.5,0.5,0.3333333333333333,0.5,0.5,0.5,0.77777777777777777,0.5,0.3333333333333333,0.5,0.5]
        ptr=f2.getValues();
        for i in xrange(20):
            self.assertTrue(abs(expected2[i]-ptr[i])<1e-12);
            pass
        m1.changeSpaceDimension(3);
        f1=m1.getMeasureField(False);
        self.assertEqual(10,f1.getArray().getNumberOfTuples());
        self.assertEqual(1,f1.getNumberOfComponents());
        ptr=f1.getArray().getValues();
        for i in xrange(10):
            self.assertTrue(abs(abs(expected1[i])-ptr[i])<1e-12);
            pass
        f2=m1.getBarycenterAndOwner();
        self.assertEqual(10,f2.getNumberOfTuples());
        self.assertEqual(3,f2.getNumberOfComponents());
        ptr=f2.getValues();
        expected3=[0.5,0.3333333333333333,0.,0.5,0.5,0.,0.5,0.77777777777777777,0.,0.5,0.3333333333333333,0.,0.5,0.5,0., 0.5,0.3333333333333333,0.,0.5,0.5,0.,0.5,0.77777777777777777,0.,0.5,0.3333333333333333,0.,0.5,0.5,0.]
        for i in xrange(30):
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
        meshN.checkCoherency();
        res1=meshN.arePolyhedronsNotCorrectlyOriented();
        meshN.orientCorrectlyPolyhedrons();
        res1=meshN.arePolyhedronsNotCorrectlyOriented();
        self.assertTrue(len(res1)==0);
        #
        da=meshN.getBarycenterAndOwner();
        self.assertEqual(4,da.getNumberOfTuples());
        self.assertEqual(3,da.getNumberOfComponents());
        daPtr=da.getValues();
        for i in xrange(12):
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
        for j in xrange(5):
            res=f.getValueOn(loc[2*j:2*j+2]);
            for i in xrange(3):
                self.assertTrue(abs(values1[i+3*j]-res[i])<1e-12);
                pass
            pass
        f.renumberCells(renumber1,False);
        ptr=f.getArray().getValues();
        expected1=[9.,109.,10009.,8.,108.,10008.,11.,111.,10011.,7.,107.,10007.,10.,110.,10010.]
        for i in xrange(15):
            self.assertTrue(abs(expected1[i]-ptr[i])<1e-12);
            pass
        #check that fields remains the same geometrically
        for j in xrange(5):
            res=f.getValueOn(loc[2*j:2*(j+1)]);
            for i in xrange(3):
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
        f.checkCoherency();
        fCpy=f.clone(True);
        self.assertTrue(f.isEqual(fCpy,1e-12,1e-12));
        f.renumberCells(renumber1,False);
        self.assertTrue(not f.isEqual(fCpy,1e-12,1e-12));
        expected2=[21.,1021.,22.,1022.,23.,1023.,24.,1024.,25.,1025.,26.,1026., 11.,1011.,12.,1012.,13.,1013.,14.,1014.,15.,1015.,16.,1016., 41.,1041.,42.,1042., 1.,1001.,2.,1002., 31.,1031.,32.,1032.]
        ptr=f.getArray().getValues();
        for i in xrange(36):
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
        f.checkCoherency();
        fCpy=f.clone(True);
        self.assertTrue(f.isEqual(fCpy,1e-12,1e-12));
        f.renumberCells(renumber1,False);
        self.assertTrue(not f.isEqual(fCpy,1e-12,1e-12));
        expected3=[21.,1021.,22.,1022.,23.,1023.,11.,1011.,12.,1012.,13.,1013.,41.,1041.,42.,1042.,43.,1043.,44.,1044.,1.,1001.,2.,1002.,3.,1003.,4.,1004.,31.,1031.,32.,1032.,33.,1033.,34.,1034.]
        ptr=f.getArray().getValues();
        for i in xrange(36):
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
        f.checkCoherency();
        renumber1=[0,4,1,3,5,2,6,7,8]
        loc=[0.5432,-0.2432, 0.5478,0.1528]
        expected1=[9.0272, 109.0272, 10009.0272, 11.4124,111.4124,10011.4124]
        for j in xrange(2):
            res=f.getValueOn(loc[2*j:2*j+2]);
            for i in xrange(3):
                self.assertTrue(abs(expected1[i+3*j]-res[i])<1e-12);
                pass
            pass
        fCpy=f.clone(True);
        self.assertTrue(f.isEqual(fCpy,1e-12,1e-12));
        f.renumberNodes(renumber1);
        self.assertTrue(not f.isEqual(fCpy,1e-12,1e-12));
        for j in xrange(2):
            res=f.getValueOn(loc[2*j:2*j+2]);
            for i in xrange(3):
                self.assertTrue(abs(expected1[i+3*j]-res[i])<1e-12);
                pass
            pass
        expected2=[7.,107.,10007.,9.,109.,10009.,12.,112.,10012.,10.,110.,10010.,8.,108.,10008.,11.,111.,10011.,13.,113.,10013.,14.,114.,10014.,15.,115.,10015.]
        for i in xrange(27):
            self.assertTrue(abs(expected2[i]-f.getArray().getValues()[i])<1e-12);
            pass
        renumber2=[0,2,5,3,1,4,6,7,8]
        f.renumberNodes(renumber2);
        self.assertTrue(f.isEqual(fCpy,1e-12,1e-12));
        pass

    def testConvertQuadraticCellsToLinear(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_3();
        mesh.checkCoherency();
        types=mesh.getAllGeoTypes();
        types.sort()
        self.assertEqual(5,len(types));
        expected1=[NORM_POLYGON, NORM_TRI3, NORM_QUAD4, NORM_TRI6, NORM_QUAD8]
        expected1.sort()
        self.assertEqual(expected1,types);
        self.assertTrue(mesh.isPresenceOfQuadratic());
        self.assertEqual(62,mesh.getMeshLength());
        f1=mesh.getMeasureField(False);
        #
        mesh.convertQuadraticCellsToLinear();
        self.assertTrue(not mesh.isPresenceOfQuadratic());
        #
        mesh.checkCoherency();
        f2=mesh.getMeasureField(False);
        self.assertTrue(f1.getArray().isEqual(f2.getArray(),1e-12));
        self.assertEqual(48,mesh.getMeshLength());
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
        for i in xrange(9):
            self.assertEqual(expected1[i],nodeCor.getIJ(i,0));
            pass
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
        a2=a1.deepCpy();
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
        f.checkCoherency();
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
        m2.checkCoherency();
        m1.checkCoherency();
        #
        expected1=[0.25,0.125,0.125,0.25,0.25]
        f1=m1.getMeasureField(False);
        f2=m2.getMeasureField(False);
        self.assertEqual(5,f1.getArray().getNumberOfTuples());
        self.assertEqual(1,f2.getArray().getNumberOfTuples());
        for i in xrange(5):
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
        for i in xrange(5):
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
        for i in xrange(20):
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
        for i in xrange(22):
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
        f.checkCoherency();
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
        for i in xrange(10):
            self.assertAlmostEqual(expected1[i],f3.getIJ(i,0),9);
            pass
        #
        f4=f1.crossProduct(f2);
        expected2=[-93., 186., -93., -392., 784., -392., -691., 1382., -691., -990., 1980., -990., -1289., 2578., -1289., -1588., 3176., -1588., -1887., 3774., -1887., -2186., 4372., -2186., -2485., 4970., -2485., -2784., 5568., -2784.]
        for i in xrange(30):
            self.assertAlmostEqual(expected2[i],f4.getIJ(0,i),9);
            pass
        pass

    def testMinMaxFields1(self):
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
        arr2=[6.,108.,206.,9.,107.,209.,8.,110.,208.,11.,109.,211.,10.,112.,210.,13.,111.,213.,12.,114.,212.,15.,113.,215.,14.,116.,214.,17.,115.,217.]
        array.setValues(arr2,mesh1.getNumberOfCells(),3);
        f2.setArray(array);
        #
        f3=f1.max(f2);
        expected1=[7.,108.,207.,9.,108.,209.,9.,110.,209.,11.,110.,211.,11.,112.,211.,13.,112.,213.,13.,114.,213.,15.,114.,215.,15.,116.,215.,17.,116.,217.]
        for i in xrange(30):
            self.assertAlmostEqual(expected1[i],f3.getIJ(0,i),9);
            pass
        #
        f4=f1.min(f2);
        expected2=[6.,107.,206.,8.,107.,208.,8.,109.,208.,10.,109.,210.,10.,111.,210.,12.,111.,212.,12.,113.,212.,14.,113.,214.,14.,115.,214.,16.,115.,216.]
        for i in xrange(30):
            self.assertAlmostEqual(expected2[i],f4.getIJ(0,i),9);
            pass
        #
        pass

    def testApplyLin1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_3();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,LINEAR_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr=[7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.,13.,113.,14.,114.,15.,115.,16.,116.]
        array.setValues(arr,mesh1.getNumberOfCells(),2);
        f1.setArray(array);
        #
        f1.applyLin(2.,3.,0);
        expected1=[17.,107.,19.,108.,21.,109.,23.,110.,25.,111.,27.,112.,29.,113.,31.,114.,33.,115.,35.,116.]
        for i in xrange(20):
            self.assertAlmostEqual(expected1[i],f1.getIJ(0,i),9);
            pass
        #
        arr2=[2.,102.,3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.]
        array=DataArrayDouble.New();
        array.setValues(arr2,mesh1.getNumberOfCells(),2);
        f1.setEndArray(array);
        #
        f1.applyLin(4.,5.,1);
        #
        expected2=[17.,433.,19.,437.,21.,441.,23.,445.,25.,449.,27.,453.,29.,457.,31.,461.,33.,465.,35.,469.]
        for i in xrange(20):
            self.assertAlmostEqual(expected2[i],f1.getIJ(0,i),9);
            pass
        expected3=[2.,413.,3.,417.,4.,421.,5.,425.,6.,429.,7.,433.,8.,437.,9.,441.,10.,445.,11.,449.]
        for i in xrange(20):
            self.assertAlmostEqual(expected3[i],f1.getEndArray().getIJ(0,i),9);
            pass
        #
        pass

    def testGetIdsInRange1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_3();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setTime(2.3,5,6);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[2.,8.,6.,5.,11.,7.,9.,3.,10.,4.]
        array.setValues(arr1,mesh1.getNumberOfCells(),1);
        f1.setArray(array);
        #
        f1.checkCoherency();
        da=f1.getIdsInRange(2.9,7.1);
        self.failUnlessEqual(5,da.getNbOfElems());
        expected1=[2,3,5,7,9]
        self.failUnlessEqual(expected1,list(da.getValues()));
        da=f1.getIdsInRange(8.,12.);
        self.failUnlessEqual(4,da.getNbOfElems());
        expected2=[1,4,6,8]
        self.failUnlessEqual(expected2,list(da.getValues()));
        #
        pass

    def testBuildSubPart1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setTime(2.3,5,6);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[3.,103.,4.,104.,5.,105.,6.,106.,7.,107.]
        array.setValues(arr1,mesh1.getNumberOfCells(),2);
        f1.setArray(array);
        #
        part1=[2,1,4]
        f2=f1[part1];
        f2.zipCoords()
        self.failUnlessEqual(3,f2.getNumberOfTuples());
        self.failUnlessEqual(2,f2.getNumberOfComponents());
        expected1=[5.,105.,4.,104.,7.,107.]
        for i in xrange(6):
            self.assertAlmostEqual(f2.getIJ(0,i),expected1[i],12);
            pass
        self.failUnlessEqual(3,f2.getMesh().getNumberOfCells());
        self.failUnlessEqual(6,f2.getMesh().getNumberOfNodes());
        self.failUnlessEqual(2,f2.getMesh().getSpaceDimension());
        self.failUnlessEqual(2,f2.getMesh().getMeshDimension());
        m2C=f2.getMesh();
        self.failUnlessEqual(13,m2C.getMeshLength());
        expected2=[0.2, -0.3, 0.7, -0.3, 0.2, 0.2, 0.7, 0.2, 0.2, 0.7, 0.7, 0.7]
        for i in xrange(12):
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12);
            pass
        expected3=[3,2,3,1,3,0,2,1,4,4,5,3,2]
        self.failUnlessEqual(expected3,list(m2C.getNodalConnectivity().getValues()));
        expected4=[0,4,8,13]
        self.failUnlessEqual(expected4,list(m2C.getNodalConnectivityIndex().getValues()));
        # Test with field on nodes.
        f1=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
        f1.setTime(2.3,5,6);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr2=[3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.]
        array.setValues(arr2,mesh1.getNumberOfNodes(),2);
        f1.setArray(array);
        part2=[1,2]
        f2=f1.buildSubPart(part2);
        self.failUnlessEqual(4,f2.getNumberOfTuples());
        self.failUnlessEqual(2,f2.getNumberOfComponents());
        expected5=[4.,104.,5.,105.,7.,107.,8.,108.]
        for i in xrange(8):
            self.assertAlmostEqual(f2.getIJ(0,i),expected5[i],12);
            pass
        self.failUnlessEqual(2,f2.getMesh().getNumberOfCells());
        self.failUnlessEqual(4,f2.getMesh().getNumberOfNodes());
        self.failUnlessEqual(2,f2.getMesh().getSpaceDimension());
        self.failUnlessEqual(2,f2.getMesh().getMeshDimension());
        m2C=f2.getMesh();
        self.failUnlessEqual(8,m2C.getMeshLength());
        for i in xrange(8):#8 is not an error
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12);
            pass
        self.failUnlessEqual(expected3[:4],list(m2C.getNodalConnectivity().getValues())[4:]);
        self.failUnlessEqual(expected3[4:8],list(m2C.getNodalConnectivity().getValues())[:4]);
        self.failUnlessEqual(expected4[:3],list(m2C.getNodalConnectivityIndex().getValues()));
        #idem previous because nodes of cell#4 are not fully present in part3
        part3=[1,2]
        arrr=DataArrayInt.New();
        arrr.setValues(part3,2,1);
        f2=f1.buildSubPart(arrr);
        self.failUnlessEqual(4,f2.getNumberOfTuples());
        self.failUnlessEqual(2,f2.getNumberOfComponents());
        for i in xrange(8):
            self.assertAlmostEqual(f2.getIJ(0,i),expected5[i],12);
            pass
        self.failUnlessEqual(2,f2.getMesh().getNumberOfCells());
        self.failUnlessEqual(4,f2.getMesh().getNumberOfNodes());
        self.failUnlessEqual(2,f2.getMesh().getSpaceDimension());
        self.failUnlessEqual(2,f2.getMesh().getMeshDimension());
        m2C=f2.getMesh();
        self.failUnlessEqual(8,m2C.getMeshLength());
        for i in xrange(8):#8 is not an error
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12);
            pass
        self.failUnlessEqual(expected3[:4],list(m2C.getNodalConnectivity().getValues())[4:8]);
        self.failUnlessEqual(expected3[4:8],list(m2C.getNodalConnectivity().getValues())[:4]);
        self.failUnlessEqual(expected4[:3],list(m2C.getNodalConnectivityIndex().getValues()));
        #
        part4=[1,2,4]
        f2=f1.buildSubPart(part4);
        self.failUnlessEqual(6,f2.getNumberOfTuples());
        self.failUnlessEqual(2,f2.getNumberOfComponents());
        expected6=[4.,104.,5.,105.,7.,107.,8.,108.,10.,110.,11.,111.]
        for i in xrange(12):
            self.assertAlmostEqual(f2.getIJ(0,i),expected6[i],12);
            pass
        self.failUnlessEqual(3,f2.getMesh().getNumberOfCells());
        self.failUnlessEqual(6,f2.getMesh().getNumberOfNodes());
        self.failUnlessEqual(2,f2.getMesh().getSpaceDimension());
        self.failUnlessEqual(2,f2.getMesh().getMeshDimension());
        m2C=f2.getMesh();
        self.failUnlessEqual(13,m2C.getMeshLength());
        for i in xrange(12):
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12);
            pass
        self.failUnlessEqual(expected3[0:4],list(m2C.getNodalConnectivity().getValues())[4:8]);
        self.failUnlessEqual(expected3[4:8],list(m2C.getNodalConnectivity().getValues())[0:4]);
        self.failUnlessEqual(expected3[8:13],list(m2C.getNodalConnectivity().getValues())[8:13]);
        self.failUnlessEqual(expected4,list(m2C.getNodalConnectivityIndex().getValues()));
        pass

    def testDoublyContractedProduct1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5]
        array.setValues(arr1,mesh1.getNumberOfCells(),6);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f2=f1.doublyContractedProduct();
        f2.checkCoherency();
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(3906.56,f2.getIJ(i,0),9);
            pass
        #
        pass

    def testDeterminant1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,CONST_ON_TIME_INTERVAL);
        f1.setTime(2.3,5,6);
        f1.setEndTime(3.8,7,3);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5]
        array.setValues(arr1,mesh1.getNumberOfCells(),4);
        f1.setArray(array);
        #4 components
        f1.checkCoherency();
        f2=f1.determinant();
        f2.checkCoherency();
        self.assertEqual(CONST_ON_TIME_INTERVAL,f2.getTimeDiscretization());
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfValues());
        for i in xrange(5):
            self.assertAlmostEqual(-2.42,f2.getIJ(i,0),13);
            pass
        #6 components multi arrays with end array not defined
        f1=MEDCouplingFieldDouble.New(ON_NODES,LINEAR_TIME);
        f1.setTime(2.3,5,6);
        f1.setEndTime(3.8,7,3);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr2=[1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7,
              1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7]
        array.setValues(arr2,mesh1.getNumberOfNodes(),6);
        f1.setArray(array);
        self.assertRaises(InterpKernelException,f1.checkCoherency);#no end array specified !
        #
        f2=f1.determinant();
        self.assertEqual(LINEAR_TIME,f2.getTimeDiscretization());
        self.assertEqual(1,f2.getArray().getNumberOfComponents());
        self.assertEqual(9,f2.getNumberOfTuples());
        for i in xrange(9):
            self.assertAlmostEqual(137.335,f2.getIJ(i,0),10);
            pass
        #6 components multi arrays with end array defined
        array=DataArrayDouble.New();
        arr3=[7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5,
              7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5]
        array.setValues(arr3,mesh1.getNumberOfNodes(),6);
        f1.setEndArray(array);
        f1.checkCoherency();
        f2=f1.determinant();
        f2.checkCoherency();
        self.assertEqual(LINEAR_TIME,f2.getTimeDiscretization());
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(9,f2.getNumberOfTuples());
        time2,it,order=f2.getTime()
        self.assertAlmostEqual(2.3,time2,12);
        self.assertEqual(5,it);
        self.assertEqual(6,order);
        time2,it,order=f2.getEndTime()
        self.assertAlmostEqual(3.8,time2,12);
        self.assertEqual(7,it);
        self.assertEqual(3,order);
        for i in xrange(9):
            self.assertAlmostEqual(137.335,f2.getIJ(i,0),10);
            self.assertAlmostEqual(1289.685,f2.getEndArray().getIJ(i,0),9);
            pass
        #9 components
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setTime(7.8,10,2);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr4=[1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1]
        array.setValues(arr4,mesh1.getNumberOfCells(),9);
        f1.setArray(array);
        #
        f1.checkCoherency();
        f2=f1.determinant();
        f2.checkCoherency();
        self.assertEqual(ONE_TIME,f2.getTimeDiscretization());
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        time2,it,order=f2.getTime()
        self.assertAlmostEqual(7.8,time2,12);
        self.assertEqual(10,it);
        self.assertEqual(2,order);
        for i in xrange(5):
            self.assertAlmostEqual(3.267,f2.getIJ(i,0),13);
            pass
        pass

    def testEigenValues1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7]
        array.setValues(arr1,mesh1.getNumberOfCells(),6);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f2=f1.eigenValues();
        f2.checkCoherency();
        self.assertEqual(3,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        expected1=[13.638813677891717,-4.502313844635971,-2.2364998332557486]
        for i in xrange(5):
            self.assertAlmostEqual(expected1[0],f2.getIJ(i,0),13);
            self.assertAlmostEqual(expected1[1],f2.getIJ(i,1),13);
            self.assertAlmostEqual(expected1[2],f2.getIJ(i,2),13);
            pass
        pass

    def testEigenVectors1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7]
        array.setValues(arr1,mesh1.getNumberOfCells(),6);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f2=f1.eigenVectors();
        f2.checkCoherency();
        self.assertEqual(9,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        expected1=[0.5424262364180696, 0.5351201064614425, 0.6476266283176001,#eigenvect 0
                   0.7381111277307373, 0.06458838384003074, -0.6715804522117897,#eigenvect 1
                   -0.4012053603397987, 0.8423032781211455, -0.3599436712889738#eigenvect 2
                   ]
        for i in xrange(5):
            self.assertAlmostEqual(expected1[0],f2.getIJ(i,0),13);
            self.assertAlmostEqual(expected1[1],f2.getIJ(i,1),13);
            self.assertAlmostEqual(expected1[2],f2.getIJ(i,2),13);
            self.assertAlmostEqual(expected1[3],f2.getIJ(i,3),13);
            self.assertAlmostEqual(expected1[4],f2.getIJ(i,4),13);
            self.assertAlmostEqual(expected1[5],f2.getIJ(i,5),13);
            self.assertAlmostEqual(expected1[6],f2.getIJ(i,6),13);
            self.assertAlmostEqual(expected1[7],f2.getIJ(i,7),13);
            self.assertAlmostEqual(expected1[8],f2.getIJ(i,8),13);
            pass
        #
        pass

    def testInverse1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1]
        array.setValues(arr1,mesh1.getNumberOfCells(),9);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f2=f1.inverse();
        f2.checkCoherency();
        self.assertEqual(9,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        expected1=[-2.6538108356290113, 2.855831037649208, -1.1111111111111067, 3.461891643709813, -4.775022956841121, 2.2222222222222143, -1.1111111111111054, 2.222222222222214, -1.1111111111111072]
        for i in xrange(5):
            self.assertAlmostEqual(expected1[0],f2.getIJ(i,0),13);
            self.assertAlmostEqual(expected1[1],f2.getIJ(i,1),13);
            self.assertAlmostEqual(expected1[2],f2.getIJ(i,2),13);
            self.assertAlmostEqual(expected1[3],f2.getIJ(i,3),13);
            self.assertAlmostEqual(expected1[4],f2.getIJ(i,4),13);
            self.assertAlmostEqual(expected1[5],f2.getIJ(i,5),13);
            self.assertAlmostEqual(expected1[6],f2.getIJ(i,6),13);
            self.assertAlmostEqual(expected1[7],f2.getIJ(i,7),13);
            self.assertAlmostEqual(expected1[8],f2.getIJ(i,8),13);
            pass
        #
        array=DataArrayDouble.New();
        arr3=[7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5]
        array.setValues(arr3,mesh1.getNumberOfCells(),6);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f2=f1.inverse();
        f2.checkCoherency();
        self.assertEqual(6,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        expected3=[-0.3617705098531818, -0.8678630828458127, -0.026843764174972983, 0.5539957431465833, 0.13133439560823013, -0.05301294502145887]
        for i in xrange(5):
            self.assertAlmostEqual(expected3[0],f2.getIJ(i,0),13);
            self.assertAlmostEqual(expected3[1],f2.getIJ(i,1),13);
            self.assertAlmostEqual(expected3[2],f2.getIJ(i,2),13);
            self.assertAlmostEqual(expected3[3],f2.getIJ(i,3),13);
            self.assertAlmostEqual(expected3[4],f2.getIJ(i,4),13);
            self.assertAlmostEqual(expected3[5],f2.getIJ(i,5),13);
            pass
        #
        array=DataArrayDouble.New();
        arr2=[1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5]
        array.setValues(arr2,mesh1.getNumberOfCells(),4);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f2=f1.inverse();
        f2.checkCoherency();
        self.assertEqual(4,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        expected2=[-1.8595041322314059, 0.9504132231404963, 1.404958677685951, -0.49586776859504156]
        for i in xrange(5):
            self.assertAlmostEqual(expected2[0],f2.getIJ(i,0),13);
            self.assertAlmostEqual(expected2[1],f2.getIJ(i,1),13);
            self.assertAlmostEqual(expected2[2],f2.getIJ(i,2),13);
            self.assertAlmostEqual(expected2[3],f2.getIJ(i,3),13);
            pass
        #
        pass

    def testTrace1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1]
        array.setValues(arr1,mesh1.getNumberOfCells(),9);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f2=f1.trace();
        f2.checkCoherency();
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(15.9,f2.getIJ(i,0),13);
            pass
        #
        array=DataArrayDouble.New();
        arr3=[7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5]
        array.setValues(arr3,mesh1.getNumberOfCells(),6);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f2=f1.trace();
        f2.checkCoherency();
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(25.8,f2.getIJ(i,0),13);
            pass
        #
        array=DataArrayDouble.New();
        arr2=[1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5]
        array.setValues(arr2,mesh1.getNumberOfCells(),4);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f2=f1.trace();
        f2.checkCoherency();
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(5.7,f2.getIJ(i,0),13);
            pass
        #
        pass

    def testDeviator1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7]
        array.setValues(arr1,mesh1.getNumberOfCells(),6);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f2=f1.deviator();
        f2.checkCoherency();
        self.assertEqual(6,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        expected1=[-1.1,0.,1.1,4.5,5.6,6.7]
        for i in xrange(5):
            self.assertAlmostEqual(expected1[0],f2.getIJ(i,0),13);
            self.assertAlmostEqual(expected1[1],f2.getIJ(i,1),13);
            self.assertAlmostEqual(expected1[2],f2.getIJ(i,2),13);
            self.assertAlmostEqual(expected1[3],f2.getIJ(i,3),13);
            self.assertAlmostEqual(expected1[4],f2.getIJ(i,4),13);
            self.assertAlmostEqual(expected1[5],f2.getIJ(i,5),13);
            pass
        #
        pass

    def testMagnitude1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6, 1.2,2.3,3.4,4.5,5.6, 1.2,2.3,3.4,4.5,5.6, 1.2,2.3,3.4,4.5,5.6, 1.2,2.3,3.4,4.5,5.6]
        array.setValues(arr1,mesh1.getNumberOfCells(),5);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f2=f1.magnitude();
        f2.checkCoherency();
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(8.3606219864313918,f2.getIJ(i,0),13);
            pass
        #
        pass

    def testMaxPerTuple1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6, 1.2,3.4,4.5,5.6,2.3, 3.4,4.5,5.6,1.2,2.3, 5.6,1.2,2.3,3.4,4.5, 4.5,5.6,1.2,2.3,3.4]
        array.setValues(arr1,mesh1.getNumberOfCells(),5);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f2=f1.maxPerTuple();
        f2.checkCoherency();
        self.assertEqual(1,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(5.6,f2.getIJ(i,0),13);
            pass
        #
        d2,d2I=array.maxPerTupleWithCompoId()
        self.assertEqual(1,d2.getNumberOfComponents());
        self.assertEqual(5,d2.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(5.6,d2.getIJ(i,0),13);
            pass
        self.assertTrue(d2I.isEqual(DataArrayInt([4,3,2,0,1])))
        pass

    def testChangeNbOfComponents(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6, 1.2,3.4,4.5,5.6,2.3, 3.4,4.5,5.6,1.2,2.3, 5.6,1.2,2.3,3.4,4.5, 4.5,5.6,1.2,2.3,3.4]
        array.setValues(arr1,mesh1.getNumberOfCells(),5);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f1.changeNbOfComponents(3,7.77);
        f1.checkCoherency();
        self.assertEqual(3,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        expected1=[1.2,2.3,3.4, 1.2,3.4,4.5, 3.4,4.5,5.6, 5.6,1.2,2.3, 4.5,5.6,1.2]
        for i in xrange(15):
            self.assertAlmostEqual(expected1[i],f1.getIJ(0,i),13);
            pass
        f1.changeNbOfComponents(4,7.77);
        f1.checkCoherency();
        self.assertEqual(4,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        expected2=[1.2,2.3,3.4,7.77, 1.2,3.4,4.5,7.77, 3.4,4.5,5.6,7.77, 5.6,1.2,2.3,7.77, 4.5,5.6,1.2,7.77]
        for i in xrange(20):
            self.assertAlmostEqual(expected2[i],f1.getIJ(0,i),13);
            pass
        #
        pass

    def testSortPerTuple1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setMesh(mesh1);
        array=DataArrayDouble.New();
        arr1=[1.2,2.3,3.4,4.5,5.6, 1.2,3.4,4.5,5.6,2.3, 3.4,4.5,5.6,1.2,2.3, 5.6,1.2,2.3,3.4,4.5, 4.5,5.6,1.2,2.3,3.4]
        array.setValues(arr1,mesh1.getNumberOfCells(),5);
        f1.setArray(array);
        f1.checkCoherency();
        #
        f1.sortPerTuple(True);
        f1.checkCoherency();
        self.assertEqual(5,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(arr1[0],f1.getIJ(i,0),13);
            self.assertAlmostEqual(arr1[1],f1.getIJ(i,1),13);
            self.assertAlmostEqual(arr1[2],f1.getIJ(i,2),13);
            self.assertAlmostEqual(arr1[3],f1.getIJ(i,3),13);
            self.assertAlmostEqual(arr1[4],f1.getIJ(i,4),13);
            pass
        #
        f1.sortPerTuple(False);
        f1.checkCoherency();
        self.assertEqual(5,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(arr1[4],f1.getIJ(i,0),13);
            self.assertAlmostEqual(arr1[3],f1.getIJ(i,1),13);
            self.assertAlmostEqual(arr1[2],f1.getIJ(i,2),13);
            self.assertAlmostEqual(arr1[1],f1.getIJ(i,3),13);
            self.assertAlmostEqual(arr1[0],f1.getIJ(i,4),13);
            pass
        #
        pass

    def testIsEqualWithoutConsideringStr1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_1();
        #
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.setName("rr");
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        da1,da2=mesh1.checkGeoEquivalWith(mesh2,2,1e-12);
        self.assertRaises(InterpKernelException,mesh1.checkGeoEquivalWith,mesh2,0,1e-12);
        mesh2.setName("");
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.getCoords().setInfoOnComponent(0,"tty");
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.getCoords().setInfoOnComponent(0,"");
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.getCoords().setInfoOnComponent(1,"tty");
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.getCoords().setInfoOnComponent(1,"");
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        tmp=mesh2.getCoords().getIJ(0,3);
        mesh2.getCoords().setIJ(0,3,9999.);
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(not mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.getCoords().setIJ(0,3,tmp);
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        tmp2=mesh2.getNodalConnectivity().getIJ(0,4);
        mesh2.getNodalConnectivity().setIJ(0,4,0);
        self.assertTrue(not mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(not mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        mesh2.getNodalConnectivity().setIJ(0,4,tmp2);
        self.assertTrue(mesh1.isEqual(mesh2,1e-12));
        self.assertTrue(mesh1.isEqualWithoutConsideringStr(mesh2,1e-12));
        #
        f1=mesh1.getMeasureField(True);
        f2=mesh2.getMeasureField(True);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        f2.setName("ftest");
        self.assertTrue(not f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        f1.setName("ftest");
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        #
        f2.getArray().setInfoOnComponent(0,"eee");
        self.assertTrue(not f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        f2.getArray().setInfoOnComponent(0,"");
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        #
        f2.getArray().setIJ(1,0,0.123);
        self.assertTrue(not f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(not f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        f2.getArray().setIJ(1,0,0.125);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        self.assertTrue(f1.isEqualWithoutConsideringStr(f2,1e-12,1e-12));
        #
        pass
    
    def testGetNodeIdsOfCell1(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        li=mesh1.getNodeIdsOfCell(1)
        expected1=[1, 4, 2]
        self.assertEqual(expected1,list(li))
        li=mesh1.getCoordinatesOfNode(4)
        self.assertEqual(2,len(li))
        self.assertAlmostEqual(0.2,li[0],13);
        self.assertAlmostEqual(0.2,li[1],13);
        li=mesh1.getCoords().getValuesAsTuple()
        self.assertEqual(9,len(li))
        li2=mesh1.getNodalConnectivityIndex().getValuesAsTuple()
        self.assertEqual(6,len(li2))
        pass

    def testGetEdgeRatioField1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m1.setTime(3.4,5,6); m1.setTimeUnit("us");
        f1=m1.getEdgeRatioField();
        self.assertAlmostEqual(3.4,f1.getTime()[0],12) ; self.assertEqual(5,f1.getTime()[1]) ; self.assertEqual(6,f1.getTime()[2])
        self.assertEqual("us",f1.getTimeUnit())
        self.assertEqual(m1.getNumberOfCells(),f1.getNumberOfTuples());
        self.assertEqual(5,f1.getNumberOfTuples());
        self.assertEqual(1,f1.getNumberOfComponents());
        expected1=[1.,1.4142135623730951, 1.4142135623730951,1.,1.]
        for i in xrange(5):
            self.assertAlmostEqual(expected1[i],f1.getIJ(i,0),14);
            pass
        #
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        f1=m1.getEdgeRatioField();
        self.assertEqual(m1.getNumberOfCells(),f1.getNumberOfTuples());
        self.assertEqual(5,f1.getNumberOfTuples());
        self.assertEqual(1,f1.getNumberOfComponents());
        expected2=[1.4142135623730951, 1.7320508075688772, 1.7320508075688772, 1.4142135623730951, 1.4142135623730951]
        for i in xrange(5):
            self.assertAlmostEqual(expected2[i],f1.getIJ(i,0),14);
            pass
        pass

    def testFillFromAnalytic3(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
        self.assertRaises(InterpKernelException,f1.fillFromAnalytic,1,"y+x");
        f1.setMesh(m)
        f1.setName("myField");
        f1.fillFromAnalytic(1,"y+x");
        f1.checkCoherency();
        self.assertEqual(f1.getName(),"myField");
        self.assertEqual(f1.getTypeOfField(),ON_CELLS);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        values1=[-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values1),len(tmp))
        for i in xrange(len(values1)):
            self.assertTrue(abs(values1[i]-tmp[i])<1.e-12);
            pass
        #
        f1=MEDCouplingFieldDouble.New(ON_NODES,CONST_ON_TIME_INTERVAL)
        f1.setMesh(m)
        f1.fillFromAnalytic(1,"y+2*x");
        f1.setEndTime(1.2,3,4);
        f1.checkCoherency();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),CONST_ON_TIME_INTERVAL);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values2=[-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values2),len(tmp))
        for i in xrange(len(values2)):
            self.assertTrue(abs(values2[i]-tmp[i])<1.e-12);
            pass
        f1=MEDCouplingFieldDouble.New(ON_NODES,LINEAR_TIME);
        f1.setMesh(m)
        f1.fillFromAnalytic(1,"2.*x+y");
        f1.setEndTime(1.2,3,4);
        f1.checkCoherency();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),LINEAR_TIME);
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        tmp=f1.getArray().getValues();
        values2Bis=[-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1]
        self.assertEqual(len(values2Bis),len(tmp))
        for i in xrange(len(values2Bis)):
            self.assertTrue(abs(values2Bis[i]-tmp[i])<1.e-12);
            pass
        tmp=f1.getEndArray().getValues();
        self.assertEqual(len(values2Bis),len(tmp))
        for i in xrange(len(values2Bis)):
            self.assertTrue(abs(values2Bis[i]-tmp[i])<1.e-12);
            pass
        #
        f1=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
        f1.setMesh(m)
        f1.fillFromAnalytic(2,"(x+y)*IVec+2*(x+y)*JVec");
        f1.checkCoherency();
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),ONE_TIME);
        self.assertEqual(2,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        values3=[-0.6,-1.2,-0.1,-0.2,0.4,0.8,-0.1,-0.2,0.4,0.8,0.9,1.8,0.4,0.8,0.9,1.8,1.4,2.8]
        tmp=f1.getArray().getValues();
        self.assertEqual(len(values3),len(tmp))
        for i in xrange(len(values3)):
            self.assertTrue(abs(values3[i]-tmp[i])<1.e-12);
            pass
        values4=f1.accumulate();
        self.assertTrue(abs(3.6-values4[0])<1.e-12);
        self.assertTrue(abs(7.2-values4[1])<1.e-12);
        values4=f1.integral(True);
        self.assertTrue(abs(0.5-values4[0])<1.e-12);
        self.assertTrue(abs(1.-values4[1])<1.e-12);
        #
        f1=MEDCouplingFieldDouble.New(ON_NODES,NO_TIME);
        f1.setMesh(m);
        self.assertRaises(InterpKernelException,f1.fillFromAnalytic,1,"1./(x-0.2)");
        pass

    def testFieldDoubleOpEqual1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        self.assertRaises(InterpKernelException,f1.assign,0.07);
        f1.setMesh(m);
        f1.assign(0.07);
        f1.checkCoherency();
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(0.07,f1.getIJ(i,0),16);
            pass
        f1.assign(0.09);
        f1.checkCoherency();
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(0.09,f1.getIJ(i,0),16);
            pass
        #
        f1=MEDCouplingFieldDouble.New(ON_NODES,LINEAR_TIME);
        f1.setEndTime(4.5,2,3);
        f1.setMesh(m);
        f1.assign(0.08);
        f1.checkCoherency();
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        for i in xrange(9):
            self.assertAlmostEqual(0.08,f1.getIJ(i,0),16);
            pass
        self.assertEqual(1,f1.getEndArray().getNumberOfComponents());
        self.assertEqual(9,f1.getEndArray().getNumberOfTuples());
        for i in xrange(9):
            self.assertAlmostEqual(0.08,f1.getEndArray().getIJ(i,0),16);
            pass
        pass

    def testAreaBary3D2(self):
        coordsForHexa8=[-75.45749305371, 180.95495078401, 39.515472018008,
                        -9.755591679144, 23.394927935279, 5.108794294848,
                        14.337630157832, 61.705351002702, 160.42422501908,
                        -27.273893776752, 167.567731083961, 192.830034145464,
                        99.857193154796,264.499264735586,-8.287335493412,
                        144.939882761126,156.38626563134,-31.896173894226,
                        161.34096835726,182.4654895809,73.832387065572,
                        132.680430393685,255.37973247196,96.15235602819];
        volHexa8=3258520.29637466;
        baryHexa8=[43.925705821778, 155.31893955289, 65.874418109644]
        
        coordsForPenta6=[-68.199829618726,178.938498373416,62.608505919588,
                         8.461744647847,76.653979804423,165.00018874933,
                         -27.273893776752,167.567731083961,192.830034145464,
                         106.586501038965,262.629609408327,13.124533008813,
                         155.465082847275,197.414118382622,78.408350795821,
                         132.680430393685,255.37973247196,96.15235602819];
        volPenta6=944849.868507338;
        baryPenta6=[39.631002313543,182.692711783428,106.98540473964]
        
        coordsForPyra5=[132.680430393685,255.37973247196,96.15235602819,
                        -27.273893776752,167.567731083961,192.830034145464,
                        8.461744647847,76.653979804423,165.00018874933,
                        155.465082847275,197.414118382622,78.408350795821,
                        -68.199829618726,178.938498373416,62.608505919588];
        volPyra5=756943.92980254;
        baryPyra5=[29.204294116618,172.540129749156,118.01035951483]
        mesh=MEDCouplingUMesh.New("Bary3D2",3);
        coo=DataArrayDouble.New();
        tmp=coordsForHexa8+coordsForPenta6+coordsForPyra5
        coo.setValues(tmp,19,3);
        mesh.setCoords(coo);
        #
        tmpConn=[0,1,2,3,4,5,6,7]
        mesh.allocateCells(3);
        self.assertRaises(InterpKernelException,mesh.insertNextCell,NORM_HEXA8,9,tmpConn[0:8])
        mesh.insertNextCell(NORM_HEXA8,tmpConn[0:8])
        mesh.insertNextCell(NORM_PENTA6,6,[i+8 for i in tmpConn])
        mesh.insertNextCell(NORM_PYRA5,5,[i+14 for i in tmpConn])
        mesh.finishInsertingCells();
        mesh.checkCoherency();
        mesh.mergeNodes(1e-7)
        self.assertEqual(12,mesh.getNumberOfNodes());
        vols=mesh.getMeasureField(True);
        self.assertEqual(3,vols.getNumberOfTuples());
        self.assertEqual(1,vols.getNumberOfComponents());
        self.assertAlmostEqual(volHexa8,vols.getIJ(0,0),6);
        self.assertAlmostEqual(volPenta6,vols.getIJ(1,0),7);
        self.assertAlmostEqual(volPyra5,vols.getIJ(2,0),7);
        bary=mesh.getBarycenterAndOwner();
        self.assertEqual(3,bary.getNumberOfTuples());
        self.assertEqual(3,bary.getNumberOfComponents());
        self.assertAlmostEqual(baryHexa8[0],bary.getIJ(0,0),11);
        self.assertAlmostEqual(baryHexa8[1],bary.getIJ(0,1),11);
        self.assertAlmostEqual(baryHexa8[2],bary.getIJ(0,2),11);
        self.assertAlmostEqual(baryPenta6[0],bary.getIJ(1,0),11);
        self.assertAlmostEqual(baryPenta6[1],bary.getIJ(1,1),11);
        self.assertAlmostEqual(baryPenta6[2],bary.getIJ(1,2),11);
        self.assertAlmostEqual(baryPyra5[0],bary.getIJ(2,0),11);
        self.assertAlmostEqual(baryPyra5[1],bary.getIJ(2,1),11);
        self.assertAlmostEqual(baryPyra5[2],bary.getIJ(2,2),11);
        pass

    def testGetMeasureFieldCMesh1(self):
        m=MEDCouplingCMesh.New();
        da=DataArrayDouble.New();
        discX=[2.3,3.4,5.8,10.2]
        discY=[12.3,23.4,45.8]
        discZ=[-0.7,1.2,1.25,2.13,2.67]
        da.setValues(discX,4,1);
        m.setCoordsAt(0,da);
        m.checkCoherency();
        self.assertEqual(4,m.getNumberOfNodes());
        self.assertEqual(3,m.getNumberOfCells());
        self.assertEqual(1,m.getSpaceDimension());
        f=m.getMeasureField(True);
        self.assertEqual(3,f.getNumberOfTuples());
        self.assertEqual(1,f.getNumberOfComponents());
        expected1=[1.1,2.4,4.4]
        for i in xrange(3):
            self.assertAlmostEqual(expected1[i],f.getIJ(i,0),12);
            pass
        coords=m.getCoordinatesAndOwner();
        self.assertEqual(4,coords.getNumberOfTuples());
        self.assertEqual(1,coords.getNumberOfComponents());
        for i in xrange(4):
            self.assertAlmostEqual(discX[i],coords.getIJ(i,0),12);
            pass
        coords=m.getBarycenterAndOwner();
        self.assertEqual(3,coords.getNumberOfTuples());
        self.assertEqual(1,coords.getNumberOfComponents());
        expected1_3=[2.85,4.6,8.]
        for i in xrange(3):
            self.assertAlmostEqual(expected1_3[i],coords.getIJ(i,0),12);
            pass
        #
        da=DataArrayDouble.New();
        da.setValues(discY,3,1);
        m.setCoordsAt(1,da);
        m.checkCoherency();
        self.assertEqual(12,m.getNumberOfNodes());
        self.assertEqual(6,m.getNumberOfCells());
        self.assertEqual(2,m.getSpaceDimension());
        f=m.getMeasureField(True);
        self.assertEqual(6,f.getNumberOfTuples());
        self.assertEqual(1,f.getNumberOfComponents());
        expected2=[12.21,26.64,48.84,24.64,53.76,98.56]
        for i in xrange(6):
            self.assertAlmostEqual(expected2[i],f.getIJ(i,0),12);
            pass
        coords=m.getCoordinatesAndOwner();
        self.assertEqual(12,coords.getNumberOfTuples());
        self.assertEqual(2,coords.getNumberOfComponents());
        expected2_2=[2.3,12.3,3.4,12.3,5.8,12.3,10.2,12.3, 2.3,23.4,3.4,23.4,5.8,23.4,10.2,23.4, 2.3,45.8,3.4,45.8,5.8,45.8,10.2,45.8]
        for i in xrange(24):
            self.assertAlmostEqual(expected2_2[i],coords.getIJ(0,i),12);
            pass
        coords=m.getBarycenterAndOwner();
        self.assertEqual(6,coords.getNumberOfTuples());
        self.assertEqual(2,coords.getNumberOfComponents());
        expected2_3=[2.85,17.85,4.6,17.85,8.,17.85, 2.85,34.6,4.6,34.6,8.,34.6]
        for i in xrange(12):
            self.assertAlmostEqual(expected2_3[i],coords.getIJ(0,i),12);
            pass
        #
        da=DataArrayDouble.New();
        da.setValues(discZ,5,1);
        m.setCoordsAt(2,da);
        m.checkCoherency();
        self.assertEqual(60,m.getNumberOfNodes());
        self.assertEqual(24,m.getNumberOfCells());
        self.assertEqual(3,m.getSpaceDimension());
        f=m.getMeasureField(True);
        self.assertEqual(24,f.getNumberOfTuples());
        self.assertEqual(1,f.getNumberOfComponents());
        expected3=[23.199, 50.616, 92.796, 46.816, 102.144, 187.264, 0.6105, 1.332, 2.442, 1.232, 2.688, 4.928, 10.7448, 23.4432, 42.9792, 21.6832, 47.3088, 86.7328, 6.5934, 14.3856, 26.3736, 13.3056, 29.0304, 53.2224]
        for i in xrange(24):
            self.assertAlmostEqual(expected3[i],f.getIJ(i,0),12);
            pass
        coords=m.getCoordinatesAndOwner();
        self.assertEqual(60,coords.getNumberOfTuples());
        self.assertEqual(3,coords.getNumberOfComponents());
        expected3_2=[
            2.3,12.3,-0.7, 3.4,12.3,-0.7, 5.8,12.3,-0.7, 10.2,12.3,-0.7, 2.3,23.4,-0.7, 3.4,23.4,-0.7, 5.8,23.4,-0.7, 10.2,23.4,-0.7, 2.3,45.8,-0.7, 3.4,45.8,-0.7, 5.8,45.8,-0.7, 10.2,45.8,-0.7,
            2.3,12.3,1.2, 3.4,12.3,1.2, 5.8,12.3,1.2, 10.2,12.3,1.2, 2.3,23.4,1.2, 3.4,23.4,1.2, 5.8,23.4,1.2, 10.2,23.4,1.2, 2.3,45.8,1.2, 3.4,45.8,1.2, 5.8,45.8,1.2, 10.2,45.8,1.2,
            2.3,12.3,1.25, 3.4,12.3,1.25, 5.8,12.3,1.25, 10.2,12.3,1.25, 2.3,23.4,1.25, 3.4,23.4,1.25, 5.8,23.4,1.25, 10.2,23.4,1.25, 2.3,45.8,1.25, 3.4,45.8,1.25, 5.8,45.8,1.25, 10.2,45.8,1.25,
            2.3,12.3,2.13, 3.4,12.3,2.13, 5.8,12.3,2.13, 10.2,12.3,2.13, 2.3,23.4,2.13, 3.4,23.4,2.13, 5.8,23.4,2.13, 10.2,23.4,2.13, 2.3,45.8,2.13, 3.4,45.8,2.13, 5.8,45.8,2.13, 10.2,45.8,2.13,
            2.3,12.3,2.67, 3.4,12.3,2.67, 5.8,12.3,2.67, 10.2,12.3,2.67, 2.3,23.4,2.67, 3.4,23.4,2.67, 5.8,23.4,2.67, 10.2,23.4,2.67, 2.3,45.8,2.67, 3.4,45.8,2.67, 5.8,45.8,2.67, 10.2,45.8,2.67];
        for i in xrange(180):
            self.assertAlmostEqual(expected3_2[i],coords.getIJ(0,i),12);
            pass
        coords=m.getBarycenterAndOwner();
        self.assertEqual(24,coords.getNumberOfTuples());
        self.assertEqual(3,coords.getNumberOfComponents());
        expected3_3=[
            2.85,17.85,0.25,4.6,17.85,0.25,8.,17.85,0.25, 2.85,34.6,0.25,4.6,34.6,0.25,8.,34.6,0.25,
            2.85,17.85,1.225,4.6,17.85,1.225,8.,17.85,1.225, 2.85,34.6,1.225,4.6,34.6,1.225,8.,34.6,1.225,
            2.85,17.85,1.69,4.6,17.85,1.69,8.,17.85,1.69, 2.85,34.6,1.69,4.6,34.6,1.69,8.,34.6,1.69,
            2.85,17.85,2.4,4.6,17.85,2.4,8.,17.85,2.4, 2.85,34.6,2.4,4.6,34.6,2.4,8.,34.6,2.4];
        for i in xrange(72):
            self.assertAlmostEqual(expected3_3[i],coords.getIJ(0,i),12);
            pass
        pass

    def testFieldDoubleZipCoords1(self):
        m=MEDCouplingDataForTest.build2DTargetMeshMergeNode_1();
        f=m.fillFromAnalytic(ON_NODES,2,"x*2.");
        f.getArray().setInfoOnComponent(0,"titi");
        f.getArray().setInfoOnComponent(1,"tutu");
        f.checkCoherency();
        self.assertEqual(18,f.getNumberOfTuples());
        self.assertEqual(2,f.getNumberOfComponents());
        expected1=[-0.6, -0.6, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, 0.4, 0.4]
        for i in xrange(36):
            self.assertAlmostEqual(expected1[i],f.getIJ(0,i),12);
            pass
        self.assertTrue(f.zipCoords());
        f.checkCoherency();
        expected2=[-0.6, -0.6, 1.4, 1.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, 0.4, 0.4]
        for i in xrange(30):
            self.assertAlmostEqual(expected2[i],f.getIJ(0,i),12);
            pass
        self.assertTrue(not f.zipCoords());
        f.checkCoherency();
        for i in xrange(30):
            self.assertAlmostEqual(expected2[i],f.getIJ(0,i),12);
            pass
        self.assertTrue(f.getArray().getInfoOnComponent(0)=="titi");
        self.assertTrue(f.getArray().getInfoOnComponent(1)=="tutu");
        pass

    def testFieldDoubleZipConnectivity1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells1=[2,3,4]
        m3_1=m2.buildPartOfMySelf(cells1,True);
        m3=m3_1;
        m4=MEDCouplingDataForTest.build2DSourceMesh_1();
        m5=MEDCouplingUMesh.MergeUMeshes(m1,m3);
        m6=MEDCouplingUMesh.MergeUMeshes(m5,m4);
        #
        self.assertEqual(10,m6.getNumberOfCells());
        self.assertEqual(22,m6.getNumberOfNodes());
        arr,areNodesMerged,newNbOfNodes=m6.mergeNodes(1e-13);
        self.assertEqual(9,m6.getNumberOfNodes());
        f=m6.fillFromAnalytic(ON_CELLS,2,"x");
        f2=m6.fillFromAnalytic(ON_NODES,2,"x");
        self.assertEqual(10,f.getNumberOfTuples());
        self.assertEqual(2,f.getNumberOfComponents());
        expected1=[-0.05, -0.05, 0.3666666666666667, 0.3666666666666667, 0.53333333333333321, 0.53333333333333321,
                   -0.05, -0.05, 0.45, 0.45, 0.53333333333333321, 0.53333333333333321, -0.05, -0.05, 0.45, 0.45,
                   0.36666666666666659, 0.36666666666666659, 0.033333333333333326, 0.033333333333333326];
        for i in xrange(20):
            self.assertAlmostEqual(expected1[i],f.getIJ(0,i),12);
            pass
        f.getArray().setInfoOnComponent(0,"titi");
        f.getArray().setInfoOnComponent(1,"tutu");
        f.checkCoherency();
        self.assertTrue(f.zipConnectivity(0));
        expected2=[-0.05, -0.05, 0.3666666666666667, 0.3666666666666667, 0.53333333333333321, 0.53333333333333321,
                   -0.05, -0.05, 0.45, 0.45, 0.36666666666666659, 0.36666666666666659, 0.033333333333333326, 0.033333333333333326];
        self.assertEqual(7,f.getNumberOfTuples());
        self.assertEqual(2,f.getNumberOfComponents());
        for i in xrange(14):
            self.assertAlmostEqual(expected2[i],f.getIJ(0,i),12);
            pass
        self.assertTrue(f.getArray().getInfoOnComponent(0)=="titi");
        self.assertTrue(f.getArray().getInfoOnComponent(1)=="tutu");
        self.assertTrue(not f.zipConnectivity(0));
        #
        expected3=[-0.3, -0.3, 0.2, 0.2, 0.7, 0.7, -0.3, -0.3, 0.2, 0.2, 0.7, 0.7,
                   -0.3, -0.3, 0.2, 0.2, 0.7, 0.7];
        self.assertEqual(9,f2.getNumberOfTuples());
        self.assertEqual(2,f2.getNumberOfComponents());
        for i in xrange(18):
            self.assertAlmostEqual(expected3[i],f2.getIJ(0,i),12);
            pass
        self.assertTrue(f2.zipConnectivity(0));
        self.assertEqual(9,f2.getNumberOfTuples());
        self.assertEqual(2,f2.getNumberOfComponents());
        for i in xrange(18):
            self.assertAlmostEqual(expected3[i],f2.getIJ(0,i),12);
            pass
        pass

    def testDaDoubleRenumber1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        a.setInfoOnComponent(0,"toto");
        a.setInfoOnComponent(1,"tata");
        #
        arr2=[3,1,0,6,5,4,2]
        b=a.renumber(arr2);
        self.assertEqual(7,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[3.1, 13.1, 2.1, 12.1, 7.1, 17.1, 1.1, 11.1, 6.1, 16.1, 5.1, 15.1, 4.1, 14.1]
        for i in xrange(14):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.setInfoOnComponent(0,"toto");
        c.setInfoOnComponent(1,"tata");
        d=c.renumber(arr2);
        self.assertEqual(7,d.getNumberOfTuples());
        self.assertEqual(2,d.getNumberOfComponents());
        self.assertTrue(d.getInfoOnComponent(0)=="toto");
        self.assertTrue(d.getInfoOnComponent(1)=="tata");
        expected2=[3, 13, 2, 12, 7, 17, 1, 11, 6, 16, 5, 15, 4, 14]
        for i in xrange(14):
            self.assertEqual(expected2[i],d.getIJ(0,i));
            pass
        pass

    def testDaDoubleRenumberAndReduce1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        a.setInfoOnComponent(0,"toto");
        a.setInfoOnComponent(1,"tata");
        #
        arr2=[2,-1,1,-1,0,4,3]
        b=a.renumberAndReduce(arr2,5);
        self.assertEqual(5,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[5.1,15.1,3.1,13.1,1.1,11.1,7.1,17.1,6.1,16.1]
        for i in xrange(10):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.setInfoOnComponent(0,"toto");
        c.setInfoOnComponent(1,"tata");
        d=c.renumberAndReduce(arr2,5);
        self.assertEqual(5,d.getNumberOfTuples());
        self.assertEqual(2,d.getNumberOfComponents());
        self.assertTrue(d.getInfoOnComponent(0)=="toto");
        self.assertTrue(d.getInfoOnComponent(1)=="tata");
        expected2=[5,15,3,13,1,11,7,17,6,16]
        for i in xrange(10):
            self.assertEqual(expected2[i],d.getIJ(0,i));
            pass
        pass

    def testDaDoubleRenumberInPlace1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        #
        arr2=[3,1,0,6,5,4,2]
        a.renumberInPlace(arr2);
        self.assertEqual(7,a.getNumberOfTuples());
        self.assertEqual(2,a.getNumberOfComponents());
        expected1=[3.1, 13.1, 2.1, 12.1, 7.1, 17.1, 1.1, 11.1, 6.1, 16.1, 5.1, 15.1, 4.1, 14.1]
        for i in xrange(14):
            self.assertAlmostEqual(expected1[i],a.getIJ(0,i),14);
            pass
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.renumberInPlace(arr2);
        self.assertEqual(7,c.getNumberOfTuples());
        self.assertEqual(2,c.getNumberOfComponents());
        expected2=[3, 13, 2, 12, 7, 17, 1, 11, 6, 16, 5, 15, 4, 14]
        for i in xrange(14):
            self.assertEqual(expected2[i],c.getIJ(0,i));
            pass
        pass

    def testDaDoubleRenumberR1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        a.setInfoOnComponent(0,"toto");
        a.setInfoOnComponent(1,"tata");
        #
        arr2=[3,1,0,6,5,4,2]
        b=a.renumberR(arr2);
        self.assertEqual(7,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[4.1, 14.1, 2.1, 12.1, 1.1, 11.1, 7.1, 17.1, 6.1, 16.1, 5.1, 15.1, 3.1, 13.1]
        for i in xrange(14):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.setInfoOnComponent(0,"toto");
        c.setInfoOnComponent(1,"tata");
        d=c.renumberR(arr2);
        self.assertEqual(7,d.getNumberOfTuples());
        self.assertEqual(2,d.getNumberOfComponents());
        self.assertTrue(d.getInfoOnComponent(0)=="toto");
        self.assertTrue(d.getInfoOnComponent(1)=="tata");
        expected2=[4, 14, 2, 12, 1, 11, 7, 17, 6, 16, 5, 15, 3, 13]
        for i in xrange(14):
            self.assertEqual(expected2[i],d.getIJ(0,i));
            pass
        pass

    def testDaDoubleRenumberInPlaceR1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        #
        arr2=[3,1,0,6,5,4,2]
        a.renumberInPlaceR(arr2);
        self.assertEqual(7,a.getNumberOfTuples());
        self.assertEqual(2,a.getNumberOfComponents());
        expected1=[4.1, 14.1, 2.1, 12.1, 1.1, 11.1, 7.1, 17.1, 6.1, 16.1, 5.1, 15.1, 3.1, 13.1]
        for i in xrange(14):
            self.assertAlmostEqual(expected1[i],a.getIJ(0,i),14);
            pass
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.renumberInPlaceR(arr2);
        self.assertEqual(7,c.getNumberOfTuples());
        self.assertEqual(2,c.getNumberOfComponents());
        expected2=[4, 14, 2, 12, 1, 11, 7, 17, 6, 16, 5, 15, 3, 13]
        for i in xrange(14):
            self.assertEqual(expected2[i],c.getIJ(0,i));
            pass
        pass

    def testDaDoubleSelectByTupleId1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        a.setInfoOnComponent(0,"toto");
        a.setInfoOnComponent(1,"tata");
        #
        arr2=[4,2,0,6,5]
        b=a.selectByTupleId(arr2);
        self.assertEqual(5,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[5.1,15.1,3.1,13.1,1.1,11.1,7.1,17.1,6.1,16.1]
        for i in xrange(10):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.setInfoOnComponent(0,"toto");
        c.setInfoOnComponent(1,"tata");
        d=c.selectByTupleId(arr2);
        self.assertEqual(5,d.getNumberOfTuples());
        self.assertEqual(2,d.getNumberOfComponents());
        self.assertTrue(d.getInfoOnComponent(0)=="toto");
        self.assertTrue(d.getInfoOnComponent(1)=="tata");
        expected2=[5,15,3,13,1,11,7,17,6,16]
        for i in xrange(10):
            self.assertEqual(expected2[i],d.getIJ(0,i));
            pass
        pass

    def testDaDoubleGetMinMaxValues1(self):
        a=DataArrayDouble.New();
        arr1=[2.34,4.56,-6.77,4.55,4.56,2.24,2.34,1.02,4.56]
        a.setValues(arr1,9,1);
        m,where=a.getMaxValue();
        self.assertEqual(1,where);
        self.assertAlmostEqual(4.56,m,12);
        m,ws=a.getMaxValue2();
        self.assertAlmostEqual(4.56,m,12);
        self.assertEqual(3,ws.getNumberOfTuples());
        self.assertEqual(1,ws.getNumberOfComponents());
        expected1=[1,4,8]
        for i in xrange(3):
            self.assertEqual(expected1[i],ws.getIJ(i,0));
            pass
        a=DataArrayDouble.New();
        arr2=[-2.34,-4.56,6.77,-4.55,-4.56,-2.24,-2.34,-1.02,-4.56]
        a.setValues(arr2,9,1);
        m,where=a.getMinValue();
        self.assertEqual(1,where);
        self.assertAlmostEqual(-4.56,m,12);
        m,ws=a.getMinValue2();
        self.assertAlmostEqual(-4.56,m,12);
        self.assertEqual(3,ws.getNumberOfTuples());
        self.assertEqual(1,ws.getNumberOfComponents());
        for i in xrange(3):
            self.assertEqual(expected1[i],ws.getIJ(i,0));
            pass
        pass

    def testFieldDoubleGetMinMaxValues2(self):
        m2,m1=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        self.assertEqual(18,m2.getNumberOfCells());
        arr1=[8.71,4.53,-12.41,8.71,-8.71,8.7099,4.55,8.71,5.55,6.77,-1e-200,4.55,8.7099,0.,1.23,0.,2.22,8.71]
        f=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        a=DataArrayDouble.New();
        a.setValues(arr1,18,1);
        f.setArray(a);
        f.setMesh(m2);
        #
        f.checkCoherency();
        m=f.getMaxValue();
        self.assertAlmostEqual(8.71,m,12);
        m,ws=f.getMaxValue2();
        self.assertAlmostEqual(8.71,m,12);
        self.assertEqual(4,ws.getNumberOfTuples());
        self.assertEqual(1,ws.getNumberOfComponents());
        expected1=[0,3,7,17]
        for i in xrange(4):
            self.assertEqual(expected1[i],ws.getIJ(i,0));
            pass
        #
        arr2=[-8.71,-4.53,12.41,-8.71,8.71,-8.7099,-4.55,-8.71,-5.55,-6.77,1e-200,-4.55,-8.7099,0.,-1.23,0.,-2.22,-8.71]
        a.setValues(arr2,18,1);
        f.checkCoherency();
        m=f.getMinValue();
        self.assertAlmostEqual(-8.71,m,12);
        m,ws=f.getMinValue2();
        self.assertAlmostEqual(-8.71,m,12);
        self.assertEqual(4,ws.getNumberOfTuples());
        self.assertEqual(1,ws.getNumberOfComponents());
        for i in xrange(4):
            self.assertEqual(expected1[i],ws.getIJ(i,0));
            pass
        pass

    def testBuildUnstructuredCMesh1(self):
        m=MEDCouplingCMesh.New();
        da=DataArrayDouble.New();
        discX=[2.3,3.4,5.8,10.2]
        discY=[12.3,23.4,45.8]
        discZ=[-0.7,1.2,1.25,2.13,2.67]
        da.setValues(discX,4,1);
        m.setCoordsAt(0,da);
        m.checkCoherency();
        self.assertEqual(0,m.getCellContainingPoint([2.4],1e-12));
        self.assertEqual(1,m.getCellContainingPoint([3.7],1e-12));
        self.assertEqual(2,m.getCellContainingPoint([5.9],1e-12));
        self.assertEqual(-1,m.getCellContainingPoint([10.3],1e-12));
        self.assertEqual(-1,m.getCellContainingPoint([1.3],1e-12));
        #
        m2=m.buildUnstructured();
        m2.checkCoherency();
        f1=m.getMeasureField(False);
        f2=m2.getMeasureField(False);
        self.assertTrue(isinstance(f1.getMesh(),MEDCouplingCMesh))
        self.assertEqual(f1.getNumberOfTuples(),3);
        self.assertEqual(f2.getNumberOfTuples(),3);
        self.assertEqual(1,m2.getMeshDimension());
        self.assertEqual(1,m2.getSpaceDimension());
        for i in xrange(3):
            self.assertAlmostEqual(f1.getIJ(i,0),f2.getIJ(i,0),10);
            pass
        da=DataArrayDouble.New();
        da.setValues(discY,3,1);
        m.setCoordsAt(1,da);
        #
        m2=m.buildUnstructured();
        m2.checkCoherency();
        f1=m.getMeasureField(False);
        f2=m2.getMeasureField(False);
        self.assertEqual(f1.getNumberOfTuples(),6);
        self.assertEqual(f2.getNumberOfTuples(),6);
        self.assertEqual(2,m2.getMeshDimension());
        self.assertEqual(2,m2.getSpaceDimension());
        for i in xrange(6):
            self.assertAlmostEqual(f1.getIJ(i,0),f2.getIJ(i,0),10);
            pass
        #
        da=DataArrayDouble.New();
        da.setValues(discZ,5,1);
        m.setCoordsAt(2,da);
        m2=m.buildUnstructured();
        m2.checkCoherency();
        f1=m.getMeasureField(False);
        f2=m2.getMeasureField(False);
        self.assertEqual(f1.getNumberOfTuples(),24);
        self.assertEqual(f2.getNumberOfTuples(),24);
        self.assertEqual(3,m2.getMeshDimension());
        self.assertEqual(3,m2.getSpaceDimension());
        for i in xrange(24):
            self.assertAlmostEqual(f1.getIJ(i,0),f2.getIJ(i,0),10);
            pass
        #
        pos1=[5.,30.,2.]
        self.assertEqual(16,m.getCellContainingPoint(pos1,1e-12));
        #
        elems=m2.getCellsInBoundingBox([3.5,6.,12.2,25.,0.,1.5],1e-7)
        self.assertEqual([1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17],elems.getValues())
        #
        pt=[2.4,12.7,-3.4]
        m.scale(pt,3.7);
        m3=m.buildUnstructured();
        m2.scale(pt,3.7);
        self.assertTrue(m3.isEqual(m2,1e-12));
        pass

    def testDataArrayIntInvertO2NNO21(self):
        arr1=[2,0,4,1,5,3]
        da=DataArrayInt.New();
        da.setValues(arr1,6,1);
        da2=da.invertArrayO2N2N2O(6);
        self.assertEqual(6,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        expected1=[1,3,0,5,2,4]
        for i in xrange(6):
            self.assertEqual(expected1[i],da2.getIJ(i,0));
            pass
        da3=da2.invertArrayN2O2O2N(6);
        for i in xrange(6):
            self.assertEqual(arr1[i],da3.getIJ(i,0));
            pass
        #
        arr2=[3,-1,5,4,-1,0,-1,1,2,-1]
        da=DataArrayInt.New();
        da.setValues(arr2,10,1);
        da2=da.invertArrayO2N2N2O(6);
        self.assertEqual(6,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        expected2=[5,7,8,0,3,2]
        for i in xrange(6):
            self.assertEqual(expected2[i],da2.getIJ(i,0));
            pass
        da3=da2.invertArrayN2O2O2N(10);
        for i in xrange(10):
            self.assertEqual(arr2[i],da3.getIJ(i,0));
            pass
        pass
    
    def testKeepSetSelectedComponent1(self):
        arr1=[1.,2.,3.,4., 11.,12.,13.,14., 21.,22.,23.,24., 31.,32.,33.,34., 41.,42.,43.,44.]
        a1=DataArrayDouble.New();
        a1.setValues(arr1,5,4);
        expp=[21.,22.,23.,24.]
        self.assertEqual(4,len(a1.getTuple(2)));
        for i in xrange(4):
            self.assertAlmostEqual(expp[i],a1.getTuple(2)[i],12)
            pass
        a1.setInfoOnComponent(0,"aaaa");
        a1.setInfoOnComponent(1,"bbbb");
        a1.setInfoOnComponent(2,"cccc");
        a1.setInfoOnComponent(3,"dddd");
        arr2V=[1,2,1,2,0,0]
        a2=a1.keepSelectedComponents(arr2V);
        self.assertEqual(6,a2.getNumberOfComponents());
        self.assertEqual(5,a2.getNumberOfTuples());
        self.assertTrue(a2.getInfoOnComponent(0)=="bbbb");
        self.assertTrue(a2.getInfoOnComponent(1)=="cccc");
        self.assertTrue(a2.getInfoOnComponent(2)=="bbbb");
        self.assertTrue(a2.getInfoOnComponent(3)=="cccc");
        self.assertTrue(a2.getInfoOnComponent(4)=="aaaa");
        self.assertTrue(a2.getInfoOnComponent(5)=="aaaa");
        expected1=[2.,3.,2.,3.,1.,1., 12.,13.,12.,13.,11.,11., 22.,23.,22.,23.,21.,21., 32.,33.,32.,33.,31.,31., 42.,43.,42.,43.,41.,41.]
        for i in xrange(30):
            self.assertAlmostEqual(expected1[i],a2.getIJ(0,i),14);
            pass
        a3=a1.convertToIntArr();
        self.assertEqual([21,22,23,24],a3.getTuple(2))
        a4=a3.keepSelectedComponents(arr2V);
        self.assertEqual(6,a4.getNumberOfComponents());
        self.assertEqual(5,a4.getNumberOfTuples());
        self.assertTrue(a4.getInfoOnComponent(0)=="bbbb");
        self.assertTrue(a4.getInfoOnComponent(1)=="cccc");
        self.assertTrue(a4.getInfoOnComponent(2)=="bbbb");
        self.assertTrue(a4.getInfoOnComponent(3)=="cccc");
        self.assertTrue(a4.getInfoOnComponent(4)=="aaaa");
        self.assertTrue(a4.getInfoOnComponent(5)=="aaaa");
        for i in xrange(30):
            self.assertEqual(int(expected1[i]),a4.getIJ(0,i));
            pass
        # setSelectedComponents
        arr3V=[3,2]
        a5=a1.keepSelectedComponents(arr3V);
        a5.setInfoOnComponent(0,"eeee");
        a5.setInfoOnComponent(1,"ffff");
        arr4V=[1,2]
        a2.setSelectedComponents(a5,arr4V);
        self.assertEqual(6,a2.getNumberOfComponents());
        self.assertEqual(5,a2.getNumberOfTuples());
        self.assertTrue(a2.getInfoOnComponent(0)=="bbbb");
        self.assertTrue(a2.getInfoOnComponent(1)=="eeee");
        self.assertTrue(a2.getInfoOnComponent(2)=="ffff");
        self.assertTrue(a2.getInfoOnComponent(3)=="cccc");
        self.assertTrue(a2.getInfoOnComponent(4)=="aaaa");
        self.assertTrue(a2.getInfoOnComponent(5)=="aaaa");
        expected2=[2.,4.,3.,3.,1.,1., 12.,14.,13.,13.,11.,11., 22.,24.,23.,23.,21.,21., 32.,34.,33.,33.,31.,31., 42.,44.,43.,43.,41.,41.]
        for i in xrange(30):
            self.assertAlmostEqual(expected2[i],a2.getIJ(0,i),14);
            pass
        a6=a5.convertToIntArr();
        a6.setInfoOnComponent(0,"eeee");
        a6.setInfoOnComponent(1,"ffff");
        a4.setSelectedComponents(a6,arr4V);
        self.assertEqual(6,a4.getNumberOfComponents());
        self.assertEqual(5,a4.getNumberOfTuples());
        self.assertTrue(a4.getInfoOnComponent(0)=="bbbb");
        self.assertTrue(a4.getInfoOnComponent(1)=="eeee");
        self.assertTrue(a4.getInfoOnComponent(2)=="ffff");
        self.assertTrue(a4.getInfoOnComponent(3)=="cccc");
        self.assertTrue(a4.getInfoOnComponent(4)=="aaaa");
        self.assertTrue(a4.getInfoOnComponent(5)=="aaaa");
        for i in xrange(30):
            self.assertEqual(int(expected2[i]),a4.getIJ(0,i));
            pass
        # test of throw
        arr5V=[2,3,6]
        arr6V=[2,7,5]
        arr7V=[2,1,4,6]
        self.assertRaises(InterpKernelException,a2.keepSelectedComponents,arr5V);
        self.assertRaises(InterpKernelException,a2.keepSelectedComponents,arr6V);
        self.assertRaises(InterpKernelException,a2.setSelectedComponents,a1,arr7V);
        arr7V=arr7V[0:3]
        self.assertRaises(InterpKernelException,a2.setSelectedComponents,a1,arr7V);
        #
        pass

    def testKeepSetSelectedComponent2(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        arr1=[1.,2.,3.,4., 11.,12.,13.,14., 21.,22.,23.,24., 31.,32.,33.,34., 41.,42.,43.,44.]
        a1=DataArrayDouble.New();
        a1.setValues(arr1,5,4);
        a1.setInfoOnComponent(0,"aaaa");
        a1.setInfoOnComponent(1,"bbbb");
        a1.setInfoOnComponent(2,"cccc");
        a1.setInfoOnComponent(3,"dddd");
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setTime(2.3,4,5);
        f1.setMesh(m1);
        f1.setName("f1");
        f1.setArray(a1);
        f1.checkCoherency();
        #
        arr2V=[1,2,1,2,0,0]
        f2=f1.keepSelectedComponents(arr2V);
        self.assertTrue(f2.getTimeDiscretization()==ONE_TIME);
        t,dt,it=f2.getTime()
        self.assertAlmostEqual(2.3,t,13);
        self.assertEqual(4,dt);
        self.assertEqual(5,it);
        f2.checkCoherency();
        self.assertEqual(6,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        self.assertTrue(f2.getArray().getInfoOnComponent(0)=="bbbb");
        self.assertTrue(f2.getArray().getInfoOnComponent(1)=="cccc");
        self.assertTrue(f2.getArray().getInfoOnComponent(2)=="bbbb");
        self.assertTrue(f2.getArray().getInfoOnComponent(3)=="cccc");
        self.assertTrue(f2.getArray().getInfoOnComponent(4)=="aaaa");
        self.assertTrue(f2.getArray().getInfoOnComponent(5)=="aaaa");
        expected1=[2.,3.,2.,3.,1.,1., 12.,13.,12.,13.,11.,11., 22.,23.,22.,23.,21.,21., 32.,33.,32.,33.,31.,31., 42.,43.,42.,43.,41.,41.]
        for i in xrange(30):
            self.assertAlmostEqual(expected1[i],f2.getIJ(0,i),14);
            pass
        #setSelectedComponents
        arr3V=[3,2]
        f5=f1.keepSelectedComponents(arr3V);
        f5.setTime(6.7,8,9);
        f5.getArray().setInfoOnComponent(0,"eeee");
        f5.getArray().setInfoOnComponent(1,"ffff");
        f5.checkCoherency();
        arr4V=[1,2]
        f2.setSelectedComponents(f5,arr4V);
        self.assertEqual(6,f2.getNumberOfComponents());
        self.assertEqual(5,f2.getNumberOfTuples());
        f2.checkCoherency();
        t,dt,it=f2.getTime()
        self.assertAlmostEqual(2.3,t,13);
        self.assertEqual(4,dt);
        self.assertEqual(5,it);
        self.assertTrue(f2.getArray().getInfoOnComponent(0)=="bbbb");
        self.assertTrue(f2.getArray().getInfoOnComponent(1)=="eeee");
        self.assertTrue(f2.getArray().getInfoOnComponent(2)=="ffff");
        self.assertTrue(f2.getArray().getInfoOnComponent(3)=="cccc");
        self.assertTrue(f2.getArray().getInfoOnComponent(4)=="aaaa");
        self.assertTrue(f2.getArray().getInfoOnComponent(5)=="aaaa");
        expected2=[2.,4.,3.,3.,1.,1., 12.,14.,13.,13.,11.,11., 22.,24.,23.,23.,21.,21., 32.,34.,33.,33.,31.,31., 42.,44.,43.,43.,41.,41.]
        for i in xrange(30):
            self.assertAlmostEqual(expected2[i],f2.getIJ(0,i),14);
            pass
        #
        pass
    
    def testElementaryDAThrowAndSpecialCases(self):
        da=DataArrayInt.New();
        self.assertRaises(InterpKernelException, da.checkAllocated);
        self.assertRaises(InterpKernelException, da.fillWithValue, 1);
        self.assertRaises(InterpKernelException, da.iota, 1);
        da.alloc(7,1);
        da.fillWithValue(11); #11,11,11,11...
        da.iota(10); #10,11,12,13...
        
        db=DataArrayInt.New();
        db.alloc(7,2);
        
        dbl2=DataArrayDouble.New();
        dbl2.alloc(7,2);
        self.assertRaises(InterpKernelException, dbl2.isUniform, 10., 1e-15);
        self.assertRaises(InterpKernelException, dbl2.sort);
        self.assertRaises(InterpKernelException, dbl2.iota, 10.);
        
        dbl=DataArrayDouble.New();
        #DataArrayDouble not allocated yet
        self.assertRaises(InterpKernelException, dbl.iota, 10.);
        self.assertRaises(InterpKernelException, dbl.isUniform, 10., 1e-15);
        self.assertRaises(InterpKernelException, dbl.sort);
        self.assertRaises(InterpKernelException, dbl.fromNoInterlace);
        self.assertRaises(InterpKernelException, dbl.toNoInterlace);
        
        dbl.alloc(7,1);
        dbl.iota(10.);
        self.assertTrue(not dbl.isUniform(10.,1e-15));
        dbl.sort();
        self.assertTrue(dbl.isMonotonic(True, .99));
        self.assertTrue(dbl.isMonotonic(True, -.99));
        self.assertTrue(not dbl.isMonotonic(True, 1.1));
        self.assertTrue(not dbl.isMonotonic(True, -1.1));
        dbl.reverse();
        self.assertTrue(dbl.isMonotonic(False, .99));
        self.assertTrue(not dbl.isMonotonic(False, 1.1));
        self.assertTrue(not dbl.isMonotonic(False, -1.1));
        
        dc=DataArrayInt.New();
        dc.alloc(14,1);
        
        dd=DataArrayDouble.New();
        self.assertRaises(InterpKernelException, dd.checkAllocated);
        self.assertRaises(InterpKernelException, dd.fillWithValue, 1.);
        self.assertRaises(InterpKernelException, dd.iota, 1.);
        self.assertTrue(not ((dd.repr().find("No data"))==-1));
        
        dd.alloc(0,1); #Allocated but nbOfElements==0!
        self.assertTrue(not ((dd.repr().find("Number of tuples : 0"))==-1));
        self.assertTrue(not ((dd.repr().find("Empty Data"))==-1));
        dd.fillWithValue(11); #?!...ok
        dd.iota(10); #?!...ok
        self.assertTrue(dd.isMonotonic(True, 1.));  #nothing is monotonic
        self.assertTrue(dd.isMonotonic(False, 1.));
        
        self.assertRaises(InterpKernelException, db.copyStringInfoFrom, da);
        self.assertRaises(InterpKernelException, db.copyStringInfoFrom, da);
        cIds=[2,2]
        self.assertRaises(InterpKernelException, da.copyPartOfStringInfoFrom, db, cIds);
        cIds[0]=1;
        cIds[0]=-1;
        self.assertRaises(InterpKernelException, da.copyPartOfStringInfoFrom, db, cIds);
        
        info=["infoOfOneComponent"]*2;
        self.assertRaises(InterpKernelException, da.setInfoOnComponents, info);
        self.assertRaises(InterpKernelException, da.setInfoOnComponent, 1, info[0]);
        db.setInfoOnComponents(info);
        
        self.assertRaises(InterpKernelException, da.getInfoOnComponent, -1);
        self.assertRaises(InterpKernelException, da.getInfoOnComponent, 2);
        self.assertTrue(db.getInfoOnComponent(1)==db.getInfoOnComponent(0));
        self.assertRaises(InterpKernelException, db.getVarOnComponent, -1);
        self.assertRaises(InterpKernelException, db.getVarOnComponent, 2);
        self.assertRaises(InterpKernelException, db.getUnitOnComponent, -1);
        self.assertRaises(InterpKernelException, db.getUnitOnComponent, 2);
        
        self.assertTrue(da.GetVarNameFromInfo("varname unit ")=="varname unit ");
        self.assertTrue(da.GetVarNameFromInfo("varname]unit[")=="varname]unit[");
        self.assertTrue(da.GetVarNameFromInfo("[unit]")=="");
        self.assertTrue(da.GetVarNameFromInfo("varname [unit]")=="varname");
        
        self.assertTrue(da.GetUnitFromInfo("varname unit ")=="");
        self.assertTrue(da.GetUnitFromInfo("varname]unit[")=="");
        self.assertTrue(da.GetUnitFromInfo("[unit]")=="unit");
        self.assertTrue(da.GetUnitFromInfo("varname [unit]")=="unit");
        
        self.assertRaises(InterpKernelException, da.checkNbOfTuplesAndComp, db, "theMessageInThrow");
        self.assertRaises(InterpKernelException, da.checkNbOfTuplesAndComp, dc, "theMessageInThrow");
        self.assertRaises(InterpKernelException, db.checkNbOfTuplesAndComp, dc, "theMessageInThrow");
        
        self.assertRaises(InterpKernelException, da.checkNbOfTuplesAndComp, 7, 2, "theMessageInThrow");
        da.checkNbOfTuplesAndComp(7,1,"theMessageInThrow");
        
        self.assertRaises(InterpKernelException, db.checkNbOfElems, 7*2+1, "theMessageInThrow");
        db.checkNbOfElems(7*2,"theMessageInThrow");
        
        self.assertRaises(InterpKernelException, db.GetNumberOfItemGivenBES, 10, 9, 1, "theMessageInThrow");
        self.assertRaises(InterpKernelException, db.GetNumberOfItemGivenBES, 0, 1, -1, "theMessageInThrow");
        self.assertEqual(10,db.GetNumberOfItemGivenBES(0,10,1,"theMessageInThrow"));
        self.assertEqual(5,db.GetNumberOfItemGivenBES(0,10,2,"theMessageInThrow"));
        self.assertEqual(6,db.GetNumberOfItemGivenBES(0,11,2,"theMessageInThrow"));
        
        self.assertTrue(not ((da.repr().find("Number of components : 1"))==-1));
        self.assertTrue(not ((dd.repr().find("Number of components : 1"))==-1));
        self.assertTrue(not ((dbl.repr().find("Number of components : 1"))==-1));
        
        self.assertTrue(not ((da.reprZip().find("Number of components : 1"))==-1));
        self.assertTrue(not ((dd.reprZip().find("Number of components : 1"))==-1));
        self.assertTrue(not ((dbl.reprZip().find("Number of components : 1"))==-1));
        
        self.assertRaises(InterpKernelException, dbl.selectByTupleId2, 0, 1, -1);
        self.assertRaises(InterpKernelException, dbl.substr, -1, 1);
        self.assertRaises(InterpKernelException, dbl.substr, 8, 1);
        self.assertRaises(InterpKernelException, dbl.substr, 0, 8);
        self.assertRaises(InterpKernelException, dbl.meldWith, dd);
        
        self.assertRaises(InterpKernelException, dbl.setPartOfValuesAdv, dbl2, da); #dbl dbl2 not have the same number of components
        self.assertRaises(InterpKernelException, dbl.setPartOfValuesAdv, dd, da);  #da tuple selector DataArrayInt instance not have exactly 2 components
        
        dbl3=DataArrayDouble.New();
        dbl3.alloc(6,2);
        dbl3.fillWithValue(11.);
        #bad number of components
        self.assertRaises(InterpKernelException, dbl3.getMaxValue);
        self.assertRaises(InterpKernelException, dd.getMaxValue);
        self.assertRaises(InterpKernelException, dbl3.getMinValue);
        self.assertRaises(InterpKernelException, dd.getMinValue);
        self.assertRaises(InterpKernelException, dbl3.getAverageValue);
        self.assertRaises(InterpKernelException, dd.getAverageValue);
        self.assertRaises(InterpKernelException, dd.accumulate, 100);
        self.assertRaises(InterpKernelException, dbl.fromPolarToCart);
        self.assertRaises(InterpKernelException, dbl3.fromCylToCart);
        self.assertRaises(InterpKernelException, dbl3.fromSpherToCart);
        self.assertRaises(InterpKernelException, dbl3.doublyContractedProduct);
        self.assertRaises(InterpKernelException, dbl3.determinant);
        self.assertRaises(InterpKernelException, dbl3.eigenValues);
        self.assertRaises(InterpKernelException, dbl3.eigenVectors);
        self.assertRaises(InterpKernelException, dbl3.inverse);
        self.assertRaises(InterpKernelException, dbl3.trace);
        self.assertRaises(InterpKernelException, dbl3.deviator);
        
        dbl3.setIJ(5,1,12.);
        self.assertTrue(dbl3.getMaxValueInArray()==12.);
        self.assertTrue(dbl3.getMinValueInArray()==11.);
        
        db.fillWithValue(100); #bad Ids
        self.assertRaises(InterpKernelException, dbl3.setPartOfValuesAdv, dbl2, db);
        db.fillWithValue(-1); #bad Ids
        self.assertRaises(InterpKernelException, dbl3.setPartOfValuesAdv, dbl2, db);
        db.fillWithValue(6); #bad Ids for dbl3
        self.assertRaises(InterpKernelException, dbl3.setPartOfValuesAdv, dbl2, db);
        
        dbl3.checkNoNullValues();
        dbl3.setIJ(5,0,0.);
        self.assertRaises(InterpKernelException, dbl3.checkNoNullValues);
        self.assertRaises(InterpKernelException, dbl3.applyInv, 1.);  #div by zero
        self.assertRaises(InterpKernelException, dbl2.getIdsInRange, 1., 2.);
        a=[]
        self.assertRaises(InterpKernelException, DataArrayDouble_Aggregate, a);
        self.assertRaises(InterpKernelException, DataArrayDouble_Meld, a);
        
        a=[dbl2,dbl]; #Nb of components mismatch
        self.assertRaises(InterpKernelException, DataArrayDouble_Aggregate, a);
        
        self.assertRaises(InterpKernelException, DataArrayDouble_Dot, dbl2, dbl);
        
        self.assertRaises(InterpKernelException, DataArrayDouble_CrossProduct, dbl2, dbl); #Nb of components mismatch
        self.assertRaises(InterpKernelException, DataArrayDouble_CrossProduct, dbl2, dbl2); #Nb of components must be equal to 3
        dbl4=DataArrayDouble.New();
        dbl4.alloc(6,3);
        dbl5=DataArrayDouble.New();
        dbl5.alloc(7,3);
        self.assertRaises(InterpKernelException, DataArrayDouble_CrossProduct, dbl4, dbl5); #Nb of tuples mismatch
        
        a[0]=dbl4; #Nb of tuple mismatch
        a[1]=dbl5; #Nb of tuple mismatch
        self.assertRaises(InterpKernelException, DataArrayDouble_Meld, a);
        self.assertRaises(InterpKernelException, DataArrayDouble_Dot, dbl4, dbl5);
        pass

    def testDAIGetIdsEqual1(self):
        tab1=[5,-2,-4,-2,3,2,-2];
        da=DataArrayInt.New();
        da.setValues(tab1,7,1);
        da2=da.getIdsEqual(-2);
        self.assertEqual(3,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        expected1=[1,3,6];
        self.assertEqual(expected1,da2.getValues());
        pass

    def testDAIGetIdsEqualList1(self):
        tab1=[5,-2,-4,-2,3,2,-2];
        da=DataArrayInt.New();
        da.setValues(tab1,7,1);
        da2=da.getIdsEqualList([3,-2,0]);
        self.assertEqual(4,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        expected1=[1,3,4,6];
        self.assertEqual(expected1,da2.getValues());
        pass

    def testDAFromNoInterlace1(self):
        tab1=[1,11,21,31,41,2,12,22,32,42,3,13,23,33,43]
        da=DataArrayInt.New();
        da.setValues(tab1,5,3);
        da2=da.fromNoInterlace();
        expected1=[1,2,3,11,12,13,21,22,23,31,32,33,41,42,43]
        self.assertEqual(5,da2.getNumberOfTuples());
        self.assertEqual(3,da2.getNumberOfComponents());# it's not a bug. Avoid to have 1 million components !
        self.assertEqual(expected1,da2.getValues());
        da3=da.convertToDblArr();
        da4=da3.fromNoInterlace();
        self.assertEqual(5,da4.getNumberOfTuples());
        self.assertEqual(3,da4.getNumberOfComponents());# it's not a bug. Avoid to have 1 million components !
        for i in xrange(15):
            self.assertAlmostEqual(expected1[i],da4.getIJ(0,i),14);
            pass
        pass
    
    def testDAToNoInterlace1(self):
        tab1=[1,2,3,11,12,13,21,22,23,31,32,33,41,42,43]
        da=DataArrayInt.New();
        da.setValues(tab1,5,3);
        da2=da.toNoInterlace();
        expected1=[1,11,21,31,41,2,12,22,32,42,3,13,23,33,43]
        self.assertEqual(5,da2.getNumberOfTuples());
        self.assertEqual(3,da2.getNumberOfComponents());# it's not a bug. Avoid to have 1 million components !
        self.assertEqual(expected1,da2.getValues());
        da3=da.convertToDblArr();
        da4=da3.toNoInterlace();
        self.assertEqual(5,da4.getNumberOfTuples());
        self.assertEqual(3,da4.getNumberOfComponents());# it's not a bug. Avoid to have 1 million components !
        for i in xrange(15):
            self.assertAlmostEqual(expected1[i],da4.getIJ(0,i),14);
            pass
        pass
    
    def testDAIsUniform1(self):
        tab1=[1,1,1,1,1]
        da=DataArrayInt.New();
        da.setValues(tab1,5,1);
        self.assertTrue(da.isUniform(1));
        da.setIJ(2,0,2);
        self.assertTrue(not da.isUniform(1));
        da.setIJ(2,0,1);
        self.assertTrue(da.isUniform(1));
        da2=da.convertToDblArr();
        self.assertTrue(da2.isUniform(1.,1.e-12));
        da2.setIJ(1,0,1.+1.e-13);
        self.assertTrue(da2.isUniform(1.,1.e-12));
        da2.setIJ(1,0,1.+1.e-11);
        self.assertTrue(not da2.isUniform(1.,1.e-12));
        pass
    
    def testDADFromPolarToCart1(self):
        tab1=[2.,0.2,2.5,0.7]
        da=DataArrayDouble.New();
        da.setValues(tab1,2,2);
        da2=da.fromPolarToCart();
        expected1=[1.9601331556824833,0.39733866159012243, 1.9121054682112213,1.6105442180942275]
        for i in xrange(4):
            self.assertAlmostEqual(expected1[i],da2.getIJ(0,i),13);
            pass
        pass
    
    def testDADFromCylToCart1(self):
        tab1=[2.,0.2,4.,2.5,0.7,9.]
        da=DataArrayDouble.New();
        da.setValues(tab1,2,3);
        da2=da.fromCylToCart();
        expected1=[1.9601331556824833,0.39733866159012243,4., 1.9121054682112213,1.6105442180942275,9.]
        for i in xrange(6):
            self.assertAlmostEqual(expected1[i],da2.getIJ(0,i),13);
            pass
        pass
    
    def testDADFromSpherToCart1(self):
        tab1=[2.,0.2,0.3,2.5,0.7,0.8]
        da=DataArrayDouble.New();
        da.setValues(tab1,2,3);
        da2=da.fromSpherToCart();
        expected1=[0.37959212195737485,0.11742160338765303,1.9601331556824833, 1.1220769624465328,1.1553337045129035,1.9121054682112213]
        for i in xrange(6):
            self.assertAlmostEqual(expected1[i],da2.getIJ(0,i),13);
            pass
        pass

    def testUnPolyze1(self):
        elts=[0,1,2,3,4,5,6,7]
        eltsV=elts;
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        mesh.convertToPolyTypes(eltsV);
        mesh.unPolyze();
        mesh2=MEDCouplingDataForTest.build3DTargetMesh_1();
        mesh.checkCoherency();
        self.assertTrue(mesh.isEqual(mesh2,1e-12));
        mesh.convertToPolyTypes(eltsV);
        self.assertTrue(not mesh.isEqual(mesh2,1e-12));
        mesh.getNodalConnectivity().setIJ(0,6,10);
        mesh.getNodalConnectivity().setIJ(0,7,9);
        mesh.getNodalConnectivity().setIJ(0,8,12);
        mesh.getNodalConnectivity().setIJ(0,9,13);
        mesh.unPolyze();
        self.assertTrue(mesh.isEqual(mesh2,1e-12));
        mesh.convertToPolyTypes(eltsV);
        mesh.getNodalConnectivity().setIJ(0,6,12);
        mesh.getNodalConnectivity().setIJ(0,7,13);
        mesh.getNodalConnectivity().setIJ(0,8,10);
        mesh.getNodalConnectivity().setIJ(0,9,9);
        mesh.unPolyze();
        self.assertTrue(mesh.isEqual(mesh2,1e-12));
        mesh.convertToPolyTypes(eltsV);
        mesh.getNodalConnectivity().setIJ(0,6,12);
        mesh.getNodalConnectivity().setIJ(0,7,10);
        mesh.getNodalConnectivity().setIJ(0,8,13);
        mesh.getNodalConnectivity().setIJ(0,9,9);
        mesh.unPolyze();
        self.assertTrue(not mesh.isEqual(mesh2,1e-12));
        # Test for 2D mesh
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_1();
        eltsV=eltsV[:5];
        mesh.convertToPolyTypes(eltsV);
        self.assertTrue(not mesh.isEqual(mesh2,1e-12));
        mesh.unPolyze();
        self.assertTrue(mesh.isEqual(mesh2,1e-12));
        pass

    def testConvertDegeneratedCells1(self):
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        conn=[0,1,3,3,9,10,12,12, 0,1,3,4,9,9,9,9, 1,1,1,1,10,12,9,10, 10,11,12,9,1,1,1,1]
        mesh.allocateCells(4);
        mesh.insertNextCell(NORM_HEXA8,8,conn[0:8])
        mesh.insertNextCell(NORM_HEXA8,8,conn[8:16])
        mesh.insertNextCell(NORM_HEXA8,8,conn[16:24])
        mesh.insertNextCell(NORM_HEXA8,8,conn[24:32])
        mesh.finishInsertingCells();
        mesh.checkCoherency();
        self.assertEqual(4,mesh.getNumberOfCells());
        self.assertEqual(NORM_HEXA8,mesh.getTypeOfCell(0));
        self.assertEqual(NORM_HEXA8,mesh.getTypeOfCell(1));
        self.assertEqual(NORM_HEXA8,mesh.getTypeOfCell(2));
        self.assertEqual(NORM_HEXA8,mesh.getTypeOfCell(3));
        f1=mesh.getMeasureField(True);
        mesh.convertDegeneratedCells();
        mesh.checkCoherency();
        f2=mesh.getMeasureField(True);
        self.assertEqual(4,mesh.getNumberOfCells());
        self.assertEqual(NORM_PENTA6,mesh.getTypeOfCell(0));
        self.assertEqual(NORM_PYRA5,mesh.getTypeOfCell(1));
        self.assertEqual(NORM_TETRA4,mesh.getTypeOfCell(2));
        self.assertEqual(NORM_PYRA5,mesh.getTypeOfCell(3));
        for i in xrange(4):
            self.assertAlmostEqual(f1.getArray().getIJ(0,i),f2.getArray().getIJ(0,i),5);
            pass
        pass

    def testGetNodeIdsNearPoints1(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        coords=mesh.getCoords();
        tmp=DataArrayDouble.New();
        vals=[0.2,0.2,0.1,0.2,0.2,0.2]
        tmp.setValues(vals,3,2);
        tmp2=DataArrayDouble.Aggregate(coords,tmp);
        mesh.setCoords(tmp2);
        pts=[0.2,0.2,0.1,0.3,-0.3,0.7]
        c=mesh.getNodeIdsNearPoint(pts[:2],1e-7);
        self.assertEqual([4,9,11],c.getValues());
        c,cI=mesh.getNodeIdsNearPoints(pts,3,1e-7);
        self.assertEqual([0,3,3,4],cI.getValues());
        self.assertEqual([4,9,11,6],c.getValues());
        c,cI=mesh.getNodeIdsNearPoints(pts,1e-7);
        self.assertEqual([0,3,3,4],cI.getValues());
        self.assertEqual([4,9,11,6],c.getValues());
        c,cI=mesh.getNodeIdsNearPoints(DataArrayDouble.New(pts,3,2),1e-7);
        self.assertEqual([0,3,3,4],cI.getValues());
        self.assertEqual([4,9,11,6],c.getValues());
        self.assertRaises(InterpKernelException,mesh.getNodeIdsNearPoints,DataArrayDouble.New(pts,2,3),1e-7);
        pass

    def testFieldCopyTinyAttrFrom1(self):
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setName("f1");
        f1.setTimeTolerance(1.e-5);
        f1.setDescription("f1Desc");
        f1.setTime(1.23,4,5);
        f2=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f2.setName("f2");
        f2.setDescription("f2Desc");
        f2.setTime(6.78,9,10);
        f2.setTimeTolerance(4.556e-12);
        #
        f1.copyTinyAttrFrom(f2);
        self.assertAlmostEqual(4.556e-12,f1.getTimeTolerance(),24);
        t,dt,it=f1.getTime()
        self.assertAlmostEqual(6.78,t,12);
        self.assertEqual(9,dt);
        self.assertEqual(10,it);
        self.assertTrue(f1.getName()=="f1");#name unchanged
        self.assertTrue(f1.getDescription()=="f1Desc");#description unchanged
        #
        f1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f1.setName("f1");
        f1.setTimeTolerance(1.e-5);
        f1.setDescription("f1Desc");
        f2=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f2.setName("f2");
        f2.setDescription("f2Desc");
        f2.setTimeTolerance(4.556e-12);
        #
        f1.copyTinyAttrFrom(f2);
        self.assertAlmostEqual(4.556e-12,f1.getTimeTolerance(),24);
        self.assertTrue(f1.getName()=="f1");#name unchanged
        self.assertTrue(f1.getDescription()=="f1Desc");#description unchanged
        #
        f1=MEDCouplingFieldDouble.New(ON_CELLS,CONST_ON_TIME_INTERVAL);
        f1.setName("f1");
        f1.setTimeTolerance(1.e-5);
        f1.setDescription("f1Desc");
        f1.setTime(1.23,4,5);
        f1.setEndTime(5.43,2,1);
        f2=MEDCouplingFieldDouble.New(ON_CELLS,CONST_ON_TIME_INTERVAL);
        f2.setName("f2");
        f2.setDescription("f2Desc");
        f2.setTimeTolerance(4.556e-12);
        f2.setTime(6.78,9,10);
        f2.setEndTime(10.98,7,6);
        #
        f1.copyTinyAttrFrom(f2);
        self.assertAlmostEqual(4.556e-12,f1.getTimeTolerance(),24);
        self.assertTrue(f1.getName()=="f1");#name unchanged
        self.assertTrue(f1.getDescription()=="f1Desc");#description unchanged
        t,dt,it=f1.getTime()
        self.assertAlmostEqual(6.78,t,12);
        self.assertEqual(9,dt);
        self.assertEqual(10,it);
        t,dt,it=f1.getEndTime()
        self.assertAlmostEqual(10.98,t,12);
        self.assertEqual(7,dt);
        self.assertEqual(6,it);
        #
        f1=MEDCouplingFieldDouble.New(ON_CELLS,LINEAR_TIME);
        f1.setName("f1");
        f1.setTimeTolerance(1.e-5);
        f1.setDescription("f1Desc");
        f1.setTime(1.23,4,5);
        f1.setEndTime(5.43,2,1);
        f2=MEDCouplingFieldDouble.New(ON_CELLS,LINEAR_TIME);
        f2.setName("f2");
        f2.setDescription("f2Desc");
        f2.setTimeTolerance(4.556e-12);
        f2.setTime(6.78,9,10);
        f2.setEndTime(10.98,7,6);
        #
        f1.copyTinyAttrFrom(f2);
        self.assertAlmostEqual(4.556e-12,f1.getTimeTolerance(),24);
        self.assertTrue(f1.getName()=="f1");#name unchanged
        self.assertTrue(f1.getDescription()=="f1Desc");#description unchanged
        t,dt,it=f1.getTime()
        self.assertAlmostEqual(6.78,t,12);
        self.assertEqual(9,dt);
        self.assertEqual(10,it);
        t,dt,it=f1.getEndTime()
        self.assertAlmostEqual(10.98,t,12);
        self.assertEqual(7,dt);
        self.assertEqual(6,it);
        pass

    def testExtrudedMesh5(self):
        coo1=[0.,1.,2.,3.5]
        a=DataArrayDouble.New();
        a.setValues(coo1,4,1);
        b=MEDCouplingCMesh.New();
        b.setCoordsAt(0,a);
        c=b.buildUnstructured();
        self.assertEqual(1,c.getSpaceDimension());
        c.changeSpaceDimension(2);
        #
        d=DataArrayDouble.New();
        d.alloc(13,1);
        d.iota();
        e=MEDCouplingCMesh.New();
        e.setCoordsAt(0,d);
        f=e.buildUnstructured();
        g=f.getCoords().applyFunc(2,"3.5*IVec+x/6*3.14159265359*JVec");
        self.assertRaises(InterpKernelException,f.getCoords().applyFunc,2,"3.5*IVec+x/6*3.14159265359*KVec"); # KVec refers to component #2 and there is only 2 components !
        h=g.fromPolarToCart();
        f.setCoords(h);
        i=c.buildExtrudedMesh(f,1);
        self.assertEqual(52,i.getNumberOfNodes());
        tmp,tmp2,tmp3=i.mergeNodes(1e-9);
        self.assertTrue(tmp2);
        self.assertEqual(37,tmp3);
        i.convertDegeneratedCells();
        i.checkCoherency();
        self.assertEqual(36,i.getNumberOfCells());
        self.assertEqual(37,i.getNumberOfNodes());
        self.assertEqual(12,i.getNumberOfCellsWithType(NORM_TRI3));
        self.assertEqual(24,i.getNumberOfCellsWithType(NORM_QUAD4));
        expected1=[0.25,0.75,2.0625]
        j=i.getMeasureField(True);
        for ii in xrange(12):
            for k in xrange(3):
                self.assertAlmostEqual(expected1[k],j.getIJ(0,ii*3+k),10);
                pass
            pass
        expected2=[0.62200846792814113, 0.16666666666681595, 1.4513530918323276, 0.38888888888923495, 2.6293994326053212, 0.7045454545460802, 0.45534180126145435, 0.45534180126150181, 1.0624642029433926, 1.0624642029435025, 1.9248539780597826, 1.9248539780599816, 0.16666666666661334, 0.62200846792815856, 0.38888888888876294, 1.4513530918323678, 0.70454545454522521, 2.629399432605394, -0.16666666666674007, 0.62200846792812436, -0.38888888888906142, 1.4513530918322881, -0.70454545454576778, 2.6293994326052488, -0.45534180126154766, 0.45534180126140844, -1.0624642029436118, 1.0624642029432834, -1.9248539780601803, 1.9248539780595841, -0.62200846792817499, 0.1666666666665495, -1.451353091832408, 0.388888888888613, -2.6293994326054668, 0.70454545454495332, -0.62200846792810593, -0.16666666666680507, -1.451353091832247, -0.38888888888921297, -2.6293994326051746, -0.70454545454604123, -0.45534180126135926, -0.45534180126159562, -1.0624642029431723, -1.0624642029437235, -1.9248539780593836, -1.9248539780603811, -0.1666666666664828, -0.62200846792819242, -0.38888888888846079, -1.4513530918324489, -0.70454545454467987, -2.6293994326055397, 0.16666666666687083, -0.62200846792808862, 0.38888888888936374, -1.4513530918322073, 0.70454545454631357, -2.6293994326051022, 0.45534180126164348, -0.45534180126131207, 1.0624642029438327, -1.0624642029430627, 1.9248539780605791, -1.9248539780591853, 0.62200846792821063, -0.16666666666641802, 1.4513530918324888, -0.38888888888831086, 2.6293994326056125, -0.70454545454440853]
        m=i.getBarycenterAndOwner();
        for i in xrange(72):
            self.assertAlmostEqual(expected2[i],m.getIJ(0,i),10);
            pass
        #
        pass

    def testExtrudedMesh6(self):
        coo1=[0.,1.,2.,3.5]
        a=DataArrayDouble.New();
        a.setValues(coo1,4,1);
        b=MEDCouplingCMesh.New();
        b.setCoordsAt(0,a);
        c=b.buildUnstructured();
        self.assertEqual(1,c.getSpaceDimension());
        c.changeSpaceDimension(2);
        #
        d=DataArrayDouble.New();
        d.alloc(5);
        d.iota();
        e=MEDCouplingCMesh.New();
        e.setCoordsAt(0,d);
        f=e.buildUnstructured();
        d2=f.getCoords().applyFunc("x*x/2");
        f.setCoords(d2);
        f.changeSpaceDimension(2);
        #
        center=[0.,0.]
        f.rotate(center,None,pi/3);
        g=c.buildExtrudedMesh(f,0);
        g.checkCoherency();
        expected1=[ 0.4330127018922193, 0.4330127018922193, 0.649519052838329, 1.2990381056766578, 1.299038105676658, 1.948557158514987, 2.1650635094610955, 2.1650635094610964, 3.2475952641916446, 3.031088913245533, 3.0310889132455352, 4.546633369868303 ]
        f1=g.getMeasureField(True);
        for i in xrange(12):
            self.assertAlmostEqual(expected1[i],f1.getIJ(0,i),12);
            pass
        expected2=[0.625, 0.21650635094610962, 1.625, 0.21650635094610959, 2.8750000000000004, 0.21650635094610965, 1.1250000000000002, 1.0825317547305482, 2.125, 1.0825317547305482, 3.3750000000000004, 1.0825317547305484, 2.125, 2.8145825622994254, 3.125, 2.8145825622994254, 4.375, 2.8145825622994254, 3.6250000000000009, 5.4126587736527414, 4.625, 5.4126587736527414, 5.875, 5.4126587736527414]
        f2=g.getBarycenterAndOwner();
        for i in xrange(24):
            self.assertAlmostEqual(expected2[i],f2.getIJ(0,i),12);
            pass
        pass

    def testExtrudedMesh7(self):
        coo1=[0.,1.,2.,3.5]
        a=DataArrayDouble.New();
        a.setValues(coo1,4,1);
        b=MEDCouplingCMesh.New();
        b.setCoordsAt(0,a);
        c=b.buildUnstructured();
        self.assertEqual(1,c.getSpaceDimension());
        c.changeSpaceDimension(2);
        #
        d=DataArrayDouble.New();
        d.alloc(13,1);
        d.iota();
        e=MEDCouplingCMesh.New();
        e.setCoordsAt(0,d);
        f=e.buildUnstructured();
        g=f.getCoords().applyFunc(2,"3.5*IVec+x/6*3.14159265359*JVec");
        h=g.fromPolarToCart();
        f.setCoords(h);
        i=c.buildExtrudedMesh(f,1);
        self.assertEqual(52,i.getNumberOfNodes());
        tmp,tmp2,tmp3=i.mergeNodes(1e-9);
        self.assertTrue(tmp2);
        self.assertEqual(37,tmp3);
        i.convertDegeneratedCells();
        vec1=[10.,0]
        i.translate(vec1);
        g2=h.applyFunc(3,"13.5/3.5*x*IVec+0*JVec+13.5/3.5*y*KVec");
        f.setCoords(g2);
        i.changeSpaceDimension(3);
        i3=i.buildExtrudedMesh(f,1);
        f2=i3.getMeasureField(True);
        tmp,tmp2,tmp3=i.mergeNodes(1e-9);
        self.assertTrue(tmp2);
        self.assertEqual(444,tmp3);
        expected1=[1.327751058489274, 4.2942574094314701, 13.024068164857139, 1.3069177251569044, 4.1484240761012954, 12.297505664866796, 1.270833333332571, 3.8958333333309674, 11.039062499993179, 1.2291666666659207, 3.6041666666644425, 9.585937499993932, 1.1930822748415895, 3.3515759238941376, 8.3274943351204556, 1.1722489415082769, 3.2057425905609289, 7.6009318351210622, 1.1722489415082862, 3.2057425905609884, 7.6009318351213713, 1.1930822748416161, 3.3515759238943001, 8.3274943351212727, 1.2291666666659564, 3.6041666666646734, 9.5859374999950777, 1.2708333333326081, 3.8958333333311868, 11.039062499994293, 1.3069177251569224, 4.1484240761014384, 12.297505664867627, 1.3277510584902354, 4.2942574094346071, 13.024068164866796]
        for ii in xrange(12):
            for jj in xrange(36):
                self.assertAlmostEqual(expected1[jj],f2.getIJ(0,ii*36+jj),9);
                pass
        #
        pass

    def testSimplexize1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m.convertToPolyTypes([3]);
        da=m.simplexize(0);
        self.assertEqual(7,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        expected2=[0,0,1,2,3,4,4]
        for i in xrange(7):
            self.assertEqual(expected2[i],da.getIJ(i,0));
            pass
        m.checkCoherency();
        self.assertEqual(7,m.getNumberOfCells());
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(0));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(1));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(2));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(3));
        self.assertEqual(NORM_POLYGON,m.getTypeOfCell(4));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(5));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(6));
        expected1=[0.125,0.125,0.125,0.125,0.25,0.125,0.125]
        f=m.getMeasureField(False);
        for i in xrange(7):
            self.assertAlmostEqual(expected1[i]*sqrt(2.),f.getIJ(i,0),10);
            pass
        types=m.getAllGeoTypes();
        self.assertEqual([NORM_TRI3,NORM_POLYGON],types);
        #
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m.convertToPolyTypes([3]);
        da=m.simplexize(1);
        self.assertEqual(7,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        for i in xrange(7):
            self.assertEqual(expected2[i],da.getIJ(i,0));
            pass
        m.checkCoherency();
        types=m.getAllGeoTypes();
        self.assertEqual([NORM_TRI3,NORM_POLYGON],types);
        self.assertEqual(7,m.getNumberOfCells());
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(0));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(1));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(2));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(3));
        self.assertEqual(NORM_POLYGON,m.getTypeOfCell(4));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(5));
        self.assertEqual(NORM_TRI3,m.getTypeOfCell(6));
        f=m.getMeasureField(False);
        for i in xrange(7):
            self.assertAlmostEqual(expected1[i]*sqrt(2.),f.getIJ(i,0),10);
            pass
        pass

    def testSimplexize2(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m.convertToPolyTypes([3]);
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setMesh(m);
        arr=DataArrayDouble.New();
        arr1=[10.,110.,20.,120.,30.,130.,40.,140.,50.,150.]
        arr.setValues(arr1,5,2);
        f1.setArray(arr);
        #
        f1.checkCoherency();
        self.assertTrue(f1.simplexize(0));
        f1.checkCoherency();
        expected1=[10.,110.,10.,110.,20.,120.,30.,130.,40.,140.,50.,150.,50.,150.]
        for i in xrange(14):
            self.assertAlmostEqual(expected1[i],f1.getIJ(0,i),10);
            pass
        self.assertTrue(not f1.simplexize(0));
        for i in xrange(14):
            self.assertAlmostEqual(expected1[i],f1.getIJ(0,i),10);
            pass
        #
        pass

    def testDAMeld1(self):
        da1=DataArrayDouble.New();
        da1.alloc(7,2);
        da2=DataArrayDouble.New();
        da2.alloc(7,1);
        #
        da1.fillWithValue(7.);
        da2.iota(0.);
        da3=da2.applyFunc(3,"10*x*IVec+100*x*JVec+1000*x*KVec");
        #
        da1.setInfoOnComponent(0,"c0da1");
        da1.setInfoOnComponent(1,"c1da1");
        da3.setInfoOnComponent(0,"c0da3");
        da3.setInfoOnComponent(1,"c1da3");
        da3.setInfoOnComponent(2,"c2da3");
        #
        da1C=da1.deepCpy();
        da1.meldWith(da3);
        self.assertEqual(5,da1.getNumberOfComponents());
        self.assertEqual(7,da1.getNumberOfTuples());
        self.assertTrue(da1.getInfoOnComponent(0)=="c0da1");
        self.assertTrue(da1.getInfoOnComponent(1)=="c1da1");
        self.assertTrue(da1.getInfoOnComponent(2)=="c0da3");
        self.assertTrue(da1.getInfoOnComponent(3)=="c1da3");
        self.assertTrue(da1.getInfoOnComponent(4)=="c2da3");
        #
        expected1=[7.,7.,0.,0.,0., 7.,7.,10.,100.,1000., 7.,7.,20.,200.,2000., 7.,7.,30.,300.,3000., 7.,7.,40.,400.,4000.,7.,7.,50.,500.,5000.,7.,7.,60.,600.,6000.]
        for i in xrange(35):
            self.assertAlmostEqual(expected1[i],da1.getIJ(0,i),10);
            pass
        #
        dai1=da1C.convertToIntArr();
        dai3=da3.convertToIntArr();
        dai1.meldWith(dai3);
        self.assertEqual(5,dai1.getNumberOfComponents());
        self.assertEqual(7,dai1.getNumberOfTuples());
        self.assertTrue(dai1.getInfoOnComponent(0)=="c0da1");
        self.assertTrue(dai1.getInfoOnComponent(1)=="c1da1");
        self.assertTrue(dai1.getInfoOnComponent(2)=="c0da3");
        self.assertTrue(dai1.getInfoOnComponent(3)=="c1da3");
        self.assertTrue(dai1.getInfoOnComponent(4)=="c2da3");
        for i in xrange(35):
            self.assertEqual(int(expected1[i]),dai1.getIJ(0,i));
            pass
        # test of static method DataArrayDouble::meld
        da4=DataArrayDouble.Meld(da1C,da3);
        tmp=DataArrayDouble.Meld([da1C,da3]);
        self.assertTrue(da4.isEqual(tmp,1e-10))
        self.assertEqual(5,da4.getNumberOfComponents());
        self.assertEqual(7,da4.getNumberOfTuples());
        self.assertTrue(da4.getInfoOnComponent(0)=="c0da1");
        self.assertTrue(da4.getInfoOnComponent(1)=="c1da1");
        self.assertTrue(da4.getInfoOnComponent(2)=="c0da3");
        self.assertTrue(da4.getInfoOnComponent(3)=="c1da3");
        self.assertTrue(da4.getInfoOnComponent(4)=="c2da3");
        for i in xrange(35):
            self.assertAlmostEqual(expected1[i],da4.getIJ(0,i),10);
            pass
        # test of static method DataArrayInt::meld
        dai1=da1C.convertToIntArr();
        dai4=DataArrayInt.Meld(dai1,dai3);
        tmp=DataArrayInt.Meld([dai1,dai3]);
        self.assertTrue(dai4.isEqual(tmp))
        self.assertEqual(5,dai4.getNumberOfComponents());
        self.assertEqual(7,dai4.getNumberOfTuples());
        self.assertTrue(dai4.getInfoOnComponent(0)=="c0da1");
        self.assertTrue(dai4.getInfoOnComponent(1)=="c1da1");
        self.assertTrue(dai4.getInfoOnComponent(2)=="c0da3");
        self.assertTrue(dai4.getInfoOnComponent(3)=="c1da3");
        self.assertTrue(dai4.getInfoOnComponent(4)=="c2da3");
        for i in xrange(35):
            self.assertEqual(int(expected1[i]),dai4.getIJ(0,i));
            pass
        pass

    def testFieldMeld1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setMesh(m);
        da1=DataArrayDouble.New();
        arr1=[12.,23.,34.,45.,56.]
        da1.setValues(arr1,5,1);
        da1.setInfoOnComponent(0,"aaa");
        f1.setArray(da1);
        f1.setTime(3.4,2,1);
        f1.checkCoherency();
        #
        f2=f1.deepCpy();
        f2.setMesh(f1.getMesh());
        f2.checkCoherency();
        f2.changeNbOfComponents(2,5.);
        f2.assign(5.);
        f2.getArray().setInfoOnComponent(0,"bbb");
        f2.getArray().setInfoOnComponent(1,"ccc");
        f2.checkCoherency();
        #
        f3=MEDCouplingFieldDouble.MeldFields(f2,f1);
        f3.checkCoherency();
        self.assertEqual(5,f3.getNumberOfTuples());
        self.assertEqual(3,f3.getNumberOfComponents());
        self.assertTrue(f3.getArray().getInfoOnComponent(0)=="bbb");
        self.assertTrue(f3.getArray().getInfoOnComponent(1)=="ccc");
        self.assertTrue(f3.getArray().getInfoOnComponent(2)=="aaa");
        expected1=[5.,5.,12.,5.,5.,23.,5.,5.,34.,5.,5.,45.,5.,5.,56.]
        for i in xrange(15):
            self.assertAlmostEqual(expected1[i],f3.getIJ(0,i),12);
            pass
        time,dt,it=f3.getTime();
        self.assertAlmostEqual(3.4,time,14);
        self.assertEqual(2,dt);
        self.assertEqual(1,it);
        #
        f4=f2.buildNewTimeReprFromThis(NO_TIME,False);
        f5=f1.buildNewTimeReprFromThis(NO_TIME,False);
        f6=MEDCouplingFieldDouble.MeldFields(f4,f5);
        f6.checkCoherency();
        self.assertEqual(5,f6.getNumberOfTuples());
        self.assertEqual(3,f6.getNumberOfComponents());
        self.assertTrue(f6.getArray().getInfoOnComponent(0)=="bbb");
        self.assertTrue(f6.getArray().getInfoOnComponent(1)=="ccc");
        self.assertTrue(f6.getArray().getInfoOnComponent(2)=="aaa");
        for i in xrange(15):
            self.assertAlmostEqual(expected1[i],f6.getIJ(0,i),12);
            pass
        #
        pass

    def testMergeNodes2(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        vec=[0.002,0.]
        m2.translate(vec);
        #
        m3=MEDCouplingUMesh.MergeUMeshes([m1,m2]);
        da,b,newNbOfNodes=m3.mergeNodes2(0.01);
        self.assertEqual(9,m3.getNumberOfNodes());
        expected1=[-0.299,-0.3, 0.201,-0.3, 0.701,-0.3, -0.299,0.2, 0.201,0.2, 0.701,0.2, -0.299,0.7, 0.201,0.7, 0.701,0.7]
        for i in xrange(18):
            self.assertAlmostEqual(expected1[i],m3.getCoords().getIJ(0,i),13);
            pass
        #
        pass

    def testMergeField2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setMesh(m);
        arr=DataArrayDouble.New();
        arr.alloc(5,2);
        arr.fillWithValue(2.);
        f1.setArray(arr);
        f2=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f2.setMesh(m);
        arr=DataArrayDouble.New();
        arr.alloc(5,2);
        arr.fillWithValue(5.);
        f2.setArray(arr);
        f3=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f3.setMesh(m);
        arr=DataArrayDouble.New();
        arr.alloc(5,2);
        arr.fillWithValue(7.);
        f3.setArray(arr);
        #
        f4=MEDCouplingFieldDouble.MergeFields([f1,f2,f3]);
        self.assertEqual(15,f4.getMesh().getNumberOfCells());
        expected1=[2.,2.,2.,2.,2.,2.,2.,2.,2.,2., 5.,5.,5.,5.,5.,5.,5.,5.,5.,5., 7.,7.,7.,7.,7.,7.,7.,7.,7.,7.]
        for i in xrange(30):
            self.assertAlmostEqual(expected1[i],f4.getIJ(0,i),13);
            pass
        #
        pass

    def testDAIBuildComplement1(self):
        a=DataArrayInt.New();
        tab=[3,1,7,8]
        a.setValues(tab,4,1);
        b=a.buildComplement(12);
        self.assertEqual(8,b.getNumberOfTuples());
        self.assertEqual(1,b.getNumberOfComponents());
        expected1=[0,2,4,5,6,9,10,11]
        for i in xrange(8):
            self.assertEqual(expected1[i],b.getIJ(0,i));
            pass
        pass

    def testDAIBuildUnion1(self):
        a=DataArrayInt.New();
        tab1=[3,1,7,8]
        a.setValues(tab1,4,1);
        c=DataArrayInt.New();
        tab2=[5,3,0,18,8]
        c.setValues(tab2,5,1);
        b=a.buildUnion(c);
        self.assertEqual(7,b.getNumberOfTuples());
        self.assertEqual(1,b.getNumberOfComponents());
        expected1=[0,1,3,5,7,8,18]
        for i in xrange(7):
            self.assertEqual(expected1[i],b.getIJ(0,i));
            pass
        b=DataArrayInt.BuildUnion([a,c]);
        self.assertEqual(7,b.getNumberOfTuples());
        self.assertEqual(1,b.getNumberOfComponents());
        expected1=[0,1,3,5,7,8,18]
        for i in xrange(7):
            self.assertEqual(expected1[i],b.getIJ(0,i));
            pass
        pass

    def testDAIBuildIntersection1(self):
        a=DataArrayInt.New();
        tab1=[3,1,7,8]
        a.setValues(tab1,4,1);
        c=DataArrayInt.New();
        tab2=[5,3,0,18,8]
        c.setValues(tab2,5,1);
        b=a.buildIntersection(c);
        self.assertEqual(2,b.getNumberOfTuples());
        self.assertEqual(1,b.getNumberOfComponents());
        expected1=[3,8]
        for i in xrange(2):
            self.assertEqual(expected1[i],b.getIJ(0,i));
            pass
        b=DataArrayInt.BuildIntersection([a,c]);
        self.assertEqual(2,b.getNumberOfTuples());
        self.assertEqual(1,b.getNumberOfComponents());
        expected1=[3,8]
        for i in xrange(2):
            self.assertEqual(expected1[i],b.getIJ(0,i));
            pass
        pass

    def testDAIDeltaShiftIndex1(self):
        a=DataArrayInt.New();
        tab=[1,3,6,7,7,9,15]
        a.setValues(tab,7,1);
        b=a.deltaShiftIndex();
        self.assertEqual(6,b.getNumberOfTuples());
        self.assertEqual(1,b.getNumberOfComponents());
        expected1=[2,3,1,0,2,6]
        for i in xrange(6):
            self.assertEqual(expected1[i],b.getIJ(0,i));
            pass
        pass

    def testDaDoubleSelectByTupleIdSafe1(self):
        a=DataArrayDouble.New();
        arr1=[1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1]
        a.setValues(arr1,7,2);
        a.setInfoOnComponent(0,"toto");
        a.setInfoOnComponent(1,"tata");
        #
        arr2=[4,2,0,6,5]
        b=a.selectByTupleIdSafe(arr2);
        self.assertEqual(5,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[5.1,15.1,3.1,13.1,1.1,11.1,7.1,17.1,6.1,16.1]
        for i in xrange(10):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        arr4=[4,-1,0,6,5]
        self.assertRaises(InterpKernelException,a.selectByTupleIdSafe,arr4);
        arr5=[4,2,0,6,7]
        self.assertRaises(InterpKernelException,a.selectByTupleIdSafe,arr5);
        #
        c=DataArrayInt.New();
        arr3=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        c.setValues(arr3,7,2);
        c.setInfoOnComponent(0,"toto");
        c.setInfoOnComponent(1,"tata");
        d=c.selectByTupleIdSafe(arr2);
        self.assertEqual(5,d.getNumberOfTuples());
        self.assertEqual(2,d.getNumberOfComponents());
        self.assertTrue(d.getInfoOnComponent(0)=="toto");
        self.assertTrue(d.getInfoOnComponent(1)=="tata");
        expected2=[5,15,3,13,1,11,7,17,6,16]
        for i in xrange(10):
            self.assertEqual(expected2[i],d.getIJ(0,i));
            pass
        self.assertRaises(InterpKernelException,c.selectByTupleIdSafe,arr4);
        self.assertRaises(InterpKernelException,c.selectByTupleIdSafe,arr5);
        pass

    def testAreCellsIncludedIn1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        pt=[1,3]
        m2=m.buildPartOfMySelf(pt,True);
        ret,tmp=m.areCellsIncludedIn(m2,0)
        self.assertTrue(ret);
        self.assertEqual(2,tmp.getNumberOfTuples());
        self.assertEqual(1,tmp.getNumberOfComponents());
        self.assertEqual(pt[0],tmp.getIJ(0,0));
        self.assertEqual(pt[1],tmp.getIJ(0,1));
        ret,tmp=m2.areCellsIncludedIn(m,0)
        self.assertTrue(not ret);
        m3=MEDCouplingUMesh.MergeUMeshesOnSameCoords(m,m2)
        c,cI=m3.findCommonCells(2,m.getNumberOfCells())
        self.assertTrue(c.isEqual(DataArrayInt([1,5,3,6])))
        self.assertTrue(cI.isEqual(DataArrayInt([0,2,4])))
        pass

    def testSwigErrorProtection1(self):
        m=MEDCouplingDataForTest.build3DTargetMesh_1();
        m.rotate([0.,0.,0.],[0.3,0.6,1.2],0.37)
        m.rotate([0.,0.,0.],[0.3,6,1.2],0.37)
        self.assertRaises(InterpKernelException,m.rotate,[0.,0.,0.],(0.3,6,"1.2"),0.37)
        self.assertRaises(InterpKernelException,m.rotate,[0.,"0.",0.],[0.3,0.6,1.2],0.37)
        self.assertRaises(InterpKernelException,m.rotate,[0.,0.,0.],[0.3,'0.6',1.2],0.37)
        m2=m.buildPartOfMySelf([2,5],True)
        m3=m.buildPartOfMySelf((2,5),True)
        self.assertTrue(m2.isEqual(m3,1e-12))
        self.assertRaises(InterpKernelException,m.buildPartOfMySelf,[2,5.],True)
        da1=m.getCoords().keepSelectedComponents([1])
        da2=m.getCoords().keepSelectedComponents((1,))
        self.assertTrue(da1.isEqual(da2,1e-12))
        self.assertRaises(InterpKernelException,m.getCoords().keepSelectedComponents,["1"])
        pass

    def testDAIBuildSubstraction1(self):
        a=DataArrayInt.New()
        aa=[2,3,6,8,9]
        a.setValues(aa,5,1)
        b=DataArrayInt.New()
        bb=[1,3,5,9,11]
        b.setValues(bb,5,1)
        self.assertEqual([2,6,8],a.buildSubstraction(b).getValues())
        pass

    def testBuildOrthogonalField2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        d1=DataArrayInt.New();
        d2=DataArrayInt.New();
        d3=DataArrayInt.New();
        d4=DataArrayInt.New();
        m1=m.buildDescendingConnectivity(d1,d2,d3,d4);
        #
        f1=m1.buildOrthogonalField();
        da1=f1.getArray();
        self.assertEqual(2,da1.getNumberOfComponents());
        self.assertEqual(13,da1.getNumberOfTuples());
        #
        expected1=[-1.,0.,0.,1.,1.,0.,0.,-1.,0.707106781186548,0.707106781186548,0.,-1.,0.,1.,1.,0.,0.,1.,1.,0.,-1.,0.,0.,1.,1.,0.];
        for i in xrange(26):
            self.assertAlmostEqual(expected1[i],da1.getIJ(0,i),14);
            pass
        pass

    def testSwigErrorProtection2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        coo=m.getCoords()
        c=m.getNodalConnectivity()
        ci=m.getNodalConnectivityIndex()
        del m
        self.assertEqual(2,coo.getNumberOfComponents());
        self.assertEqual(6,ci.getNumberOfTuples());
        self.assertEqual(23,c.getNumberOfTuples());
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=m.getMeasureField(True)
        c=f.getArray()
        del f
        self.assertEqual(1,c.getNumberOfComponents());
        m=MEDCouplingCMesh.New()
        x=DataArrayDouble.New()
        x.setValues([1.,2.,4.],3,1)
        m.setCoordsAt(0,x)
        del x
        xx=m.getCoordsAt(0)
        del m
        self.assertEqual(3,xx.getNumberOfTuples());
        #
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=m.getMeasureField(True)
        m2=f.getMesh()
        del m
        del f
        self.assertEqual(5,m2.getNumberOfCells());
        pass

    def testUMInsertNextCell1(self):
        targetCoords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ]
        targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        targetMesh=MEDCouplingUMesh.New();
        targetMesh.allocateCells(5);
        self.assertRaises(InterpKernelException,targetMesh.insertNextCell,NORM_QUAD4,4,targetConn[0:4])
        targetMesh.setMeshDimension(2);
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
        self.assertRaises(InterpKernelException,targetMesh.insertNextCell,NORM_TETRA4,4,targetConn[0:4])
        self.assertRaises(InterpKernelException,targetMesh.insertNextCell,NORM_SEG2,2,targetConn[0:2])
        self.assertRaises(InterpKernelException,targetMesh.insertNextCell,NORM_POINT1,1,targetConn[0:1])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[4:7])
        targetMesh.insertNextCell(NORM_TRI3,3,targetConn[7:10])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[10:14])
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[14:18])
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,9,2);
        targetMesh.setCoords(myCoords);
        targetMesh.checkCoherency();
        pass

    def testFieldOperatorDivDiffComp1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        m1,d0,d1,d2,d3=m.buildDescendingConnectivity();
        #
        f1=m1.buildOrthogonalField();
        arr1=[2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.]
        arr=DataArrayDouble.New();
        arr.setValues(arr1,13,1);
        f2=MEDCouplingFieldDouble.New(ON_CELLS);
        f2.setArray(arr);
        f2.setMesh(m1);
        f2.checkCoherency();
        #
        f3=f1/f2;
        self.assertRaises(InterpKernelException,f2.__div__,f1)
        f3.checkCoherency();
        f1/=f2;
        #self.assertRaises(InterpKernelException,f2.__idiv__,f1) # mem leaks
        self.assertTrue(f1.isEqual(f3,1e-10,1e-10));
        expected1=[-0.5, 0.0, 0.0, 0.33333333333333331, 0.25, 0.0, 0.0, -0.20000000000000001, 0.117851130197758, 0.117851130197758, 0.0, -0.14285714285714285, 0.0, 0.125, 0.1111111111111111, 0.0, 0.0, 0.10000000000000001, 0.090909090909090912, 0.0, -0.083333333333333329, 0.0, 0.0, 0.076923076923076927, 0.071428571428571425, 0.0]
        for i in xrange(26):
            self.assertAlmostEqual(expected1[i],f3.getIJ(0,i),10);
            pass
        pass

    def testDARearrange1(self):
        da1=DataArrayInt.New();
        da1.alloc(12,1);
        da1.iota(0);
        #
        self.assertEqual(12,da1.getNbOfElems());
        self.assertEqual(1,da1.getNumberOfComponents());
        self.assertEqual(12,da1.getNumberOfTuples());
        da1.rearrange(4);
        self.assertEqual(12,da1.getNbOfElems());
        self.assertEqual(4,da1.getNumberOfComponents());
        self.assertEqual(3,da1.getNumberOfTuples());
        for i in xrange(12):
            self.assertEqual(i,da1.getIJ(0,i));
        #
        da1.rearrange(6);
        self.assertEqual(12,da1.getNbOfElems());
        self.assertEqual(6,da1.getNumberOfComponents());
        self.assertEqual(2,da1.getNumberOfTuples());
        for i in xrange(12):
            self.assertEqual(i,da1.getIJ(0,i));
        #
        self.assertRaises(InterpKernelException,da1.rearrange,7);
        #
        da1.rearrange(12);
        self.assertEqual(12,da1.getNbOfElems());
        self.assertEqual(12,da1.getNumberOfComponents());
        self.assertEqual(1,da1.getNumberOfTuples());
        for i in xrange(12):
            self.assertEqual(i,da1.getIJ(0,i));
        #
        da1.rearrange(3);
        self.assertEqual(12,da1.getNbOfElems());
        self.assertEqual(3,da1.getNumberOfComponents());
        self.assertEqual(4,da1.getNumberOfTuples());
        for i in xrange(12):
            self.assertEqual(i,da1.getIJ(0,i));
        #double
        da2=da1.convertToDblArr();
        st=da2.getHiddenCppPointer()
        #
        self.assertEqual(12,da2.getNbOfElems());
        self.assertEqual(3,da2.getNumberOfComponents());
        self.assertEqual(4,da2.getNumberOfTuples());
        da2.rearrange(4);
        self.assertEqual(12,da2.getNbOfElems());
        self.assertEqual(4,da2.getNumberOfComponents());
        self.assertEqual(3,da2.getNumberOfTuples());
        for i in xrange(12):
            self.assertAlmostEqual(float(i),da2.getIJ(0,i),14);
        #
        da2.rearrange(6);
        self.assertEqual(12,da2.getNbOfElems());
        self.assertEqual(6,da2.getNumberOfComponents());
        self.assertEqual(2,da2.getNumberOfTuples());
        for i in xrange(12):
            self.assertAlmostEqual(float(i),da2.getIJ(0,i),14);
        #
        self.assertRaises(InterpKernelException,da2.rearrange,7);
        #
        da2.rearrange(1);
        self.assertEqual(st,da2.getHiddenCppPointer())
        self.assertEqual(12,da2.getNbOfElems());
        self.assertEqual(1,da2.getNumberOfComponents());
        self.assertEqual(12,da2.getNumberOfTuples());
        for i in xrange(12):
            self.assertAlmostEqual(float(i),da2.getIJ(0,i),14);
        #
        da2.rearrange(3);
        self.assertEqual(12,da2.getNbOfElems());
        self.assertEqual(3,da2.getNumberOfComponents());
        self.assertEqual(4,da2.getNumberOfTuples());
        for i in xrange(12):
            self.assertAlmostEqual(float(i),da2.getIJ(0,i),14);
        pass

    def testDARearrange2(self):
        da1=DataArrayInt.New();
        arr=[1,2,3,2,2,3,5,1,5,5,2,2]
        da1.setValues(arr,4,3);
        s=da1.getDifferentValues();
        expected1=DataArrayInt([1,2,3,5])
        self.assertTrue(expected1.isEqual(s));
        pass

    def testSwigErrorProtection3(self):
        da=DataArrayInt.New()
        da.setValues([1,2,3,4,0,0,0,0,0,0,0,0],4,3)
        self.assertEqual([1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        da=DataArrayInt.New()
        da.setValues(((1,2,3),(4,4,3),(0,0,0),(0,0,0)),4,3)
        self.assertEqual([1, 2, 3, 4, 4, 3, 0, 0, 0, 0, 0, 0],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        da.setValues((10*[1]+290*[2])[:12],4,3)
        self.assertEqual(10*[1]+[2,2],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        #
        da=DataArrayDouble.New()
        da.setValues([1,2,3.,4,0,0,0,0,0,0,0,0],4,3)
        self.assertEqual([1., 2., 3., 4., 0., 0., 0., 0., 0., 0., 0., 0.],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        da=DataArrayDouble.New()
        da.setValues(((1,2,3),(4.,4,3),(0,0,0),(0,0,0)),4,3)
        self.assertEqual([1., 2., 3., 4., 4., 3., 0., 0., 0., 0., 0., 0.],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        da.setValues((10*[1]+290*[2])[:12],4,3)
        self.assertEqual(10*[1.]+[2.,2.],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        pass

    def testDAIBuildPermutationArr1(self):
        a=DataArrayInt.New()
        a.setValues([4,5,6,7,8],5,1)
        b=DataArrayInt.New()
        b.setValues([5,4,8,6,7],5,1)
        c=a.buildPermutationArr(b)
        self.assertEqual([1,0,4,2,3],c.getValues())
        self.assertTrue(a.isEqualWithoutConsideringStrAndOrder(b))
        b.setIJ(0,0,9)
        self.assertTrue(not a.isEqualWithoutConsideringStrAndOrder(b))
        self.assertRaises(InterpKernelException,a.buildPermutationArr,b)
        a.setIJ(3,0,4)
        b.setIJ(0,0,5)
        b.setIJ(4,0,4)#a==[4,5,6,4,8] and b==[5,4,8,6,4]
        self.assertTrue(a.isEqualWithoutConsideringStrAndOrder(b))
        c=a.buildPermutationArr(b)
        self.assertEqual([1,3,4,2,3],c.getValues())
        d=b.convertToDblArr()
        expect3=[4,4,5,6,8]
        b.sort()
        self.assertEqual(expect3,b.getValues())
        d.sort()
        self.assertEqual(5,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        for i in xrange(5):
            self.assertAlmostEqual(float(expect3[i]),d.getIJ(i,0),14);
            pass
        pass

    def testAreCellsIncludedIn2(self):
        myName="Vitoo";
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=m.buildPartOfMySelf([],True);
        self.assertEqual(0,m2.getNumberOfCells());
        self.assertEqual(3,m2.getSpaceDimension());
        self.assertEqual(2,m2.getMeshDimension());
        m2.setName(myName);
        test,tmp=m.areCellsIncludedIn(m2,0)
        self.assertTrue(test);
        self.assertEqual(myName,tmp.getName());
        self.assertEqual(0,tmp.getNumberOfTuples())
        self.assertEqual(1,tmp.getNumberOfComponents())
        pass

    def testUMeshGetPartBarycenterAndOwner1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        part1=[1,0,4];
        part=DataArrayInt.New();
        part.setValues(part1,3,1);
        b=m1.getPartBarycenterAndOwner(part);
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertEqual(3,b.getNumberOfTuples());
        expected1=[0.36666666666666665,-0.13333333333333333,-0.05,-0.05,0.45,0.45];
        for i in xrange(6):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        pass

    def testUMeshGetPartMeasureField1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        part1=[1,0,4];
        part=DataArrayInt.New();
        part.setValues(part1,3,1);
        b=m1.getPartMeasureField(True,part);
        self.assertEqual(1,b.getNumberOfComponents());
        self.assertEqual(3,b.getNumberOfTuples());
        expected1=[0.125,0.25,0.25];
        for i in xrange(3):
            self.assertAlmostEqual(expected1[i],b.getIJ(0,i),14);
            pass
        pass

    def testUMeshBuildPartOrthogonalField1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m1.changeSpaceDimension(3);
        part1=[1,0,4];
        part=DataArrayInt.New();
        part.setValues(part1,3,1);
        b=m1.buildPartOrthogonalField(part);
        self.assertEqual(3,b.getArray().getNumberOfComponents());
        self.assertEqual(3,b.getArray().getNumberOfTuples());
        expected1=[0.,0.,-1.,0.,0.,-1.,0.,0.,-1.];
        for i in xrange(9):
            self.assertAlmostEqual(expected1[i],b.getArray().getIJ(0,i),14);
            pass
        pass

    def testUMeshGetTypesOfPart1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        part1=[0,3,4];
        p1=DataArrayInt.New()
        p1.setValues(part1,3,1)
        s=m1.getTypesOfPart(p1);
        self.assertEqual([NORM_QUAD4],s);
        part2=[2,2,2,1];
        p2=DataArrayInt.New()
        p2.setValues(part2,4,1)
        s=m1.getTypesOfPart(p2);
        self.assertEqual([NORM_TRI3],s);
        part3=[3,2,1];
        p3=DataArrayInt.New()
        p3.setValues(part3,3,1)
        s=m1.getTypesOfPart(p3);
        self.assertEqual(s,[NORM_TRI3,NORM_QUAD4]);
        pass

    def testUMeshKeepCellIdsByType1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        part1=[0,3,4]
        p1=DataArrayInt.New()
        p1.setValues(part1,3,1)
        p1.setName("p1")
        a=m1.keepCellIdsByType(NORM_TRI3,p1);
        self.assertEqual("p1",a.getName())
        self.assertEqual(1,a.getNumberOfComponents());
        self.assertEqual(0,a.getNumberOfTuples());
        #
        part2=[3,2,0,2,4]
        p2=DataArrayInt.New()
        p2.setValues(part2,5,1)
        p2.setName("p2")
        a=m1.keepCellIdsByType(NORM_TRI3,p2);
        self.assertEqual("p2",a.getName())
        self.assertEqual(1,a.getNumberOfComponents());
        self.assertEqual(2,a.getNumberOfTuples());
        self.assertEqual(2,a.getIJ(0,0));
        self.assertEqual(2,a.getIJ(1,0));
        #
        a=m1.keepCellIdsByType(NORM_QUAD4,p2);
        self.assertEqual("p2",a.getName())
        self.assertEqual(1,a.getNumberOfComponents());
        self.assertEqual(3,a.getNumberOfTuples());
        self.assertEqual(3,a.getIJ(0,0));
        self.assertEqual(0,a.getIJ(1,0));
        self.assertEqual(4,a.getIJ(2,0));
        pass
    
    def testSwigErrorDaIntSelectByTupleId1(self):
        a=DataArrayInt.New();
        arr1=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        a.setValues(arr1,7,2);
        a.setInfoOnComponent(0,"toto");
        a.setInfoOnComponent(1,"tata");
        #
        arr2=[4,2,0,6,5]
        b=a.selectByTupleId(arr2);
        self.assertEqual(5,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[5,15,3,13,1,11,7,17,6,16]
        self.assertEqual(expected1,b.getValues())
        #
        a2=DataArrayInt.New()
        a2.setValues(arr2,5,1)
        b=a.selectByTupleId(a2);
        self.assertEqual(5,b.getNumberOfTuples());
        self.assertEqual(2,b.getNumberOfComponents());
        self.assertTrue(b.getInfoOnComponent(0)=="toto");
        self.assertTrue(b.getInfoOnComponent(1)=="tata");
        expected1=[5,15,3,13,1,11,7,17,6,16]
        self.assertEqual(expected1,b.getValues())
        pass

    def testSwigErrorRenum(self):
        da=DataArrayDouble.New()
        da.setValues([7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.,13.,113.,14.,114.,15.,115.,16.,116.],10,2)
        d=DataArrayInt.New()
        d.setValues([0,2,3,1,4,5,6,8,7,9],10,1)
        da.renumberInPlace(d)
        da.renumber(d)
        pass

    def testSwigGetItem1(self):
        da=DataArrayInt.New()
        da.alloc(16,3)
        da.rearrange(1)
        da.iota(7)
        da.rearrange(3)
        da.setInfoOnComponent(0,"X [m]")
        da.setInfoOnComponent(1,"Y [m]")
        da.setInfoOnComponent(2,"Z [km]")
        da2=da[5:-1]
        self.assertEqual([22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51],da2.getValues())
        da2=da[4]
        self.assertEqual([19, 20, 21],da2.getValues())
        try:
            da2=da[4:17]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        da2=da[5:-2,2]
        self.assertEqual([24, 27, 30, 33, 36, 39, 42, 45, 48],da2.getValues())
        da2=da[5:8,:]
        self.assertEqual([22, 23, 24, 25, 26, 27, 28, 29, 30],da2.getValues())
        da2=da[:]
        self.assertTrue(da2.isEqual(da))
        da2=da[:,:]
        self.assertTrue(da2.isEqual(da))
        try:
            da2=da[:,:,:]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        self.assertTrue(da[5:8,-2].isEqualWithoutConsideringStr(DataArrayInt([23,26,29])))
        da2=da[5:8,:-2]
        self.assertEqual([22, 25, 28],da2.getValues())
        try:
            da2=da[5:-18,2]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        da2=da[5:5,2]
        self.assertEqual([],da2.getValues())
        pass

    def testSwigGetItem2(self):
        da=DataArrayDouble.New()
        da.alloc(16,3)
        da.rearrange(1)
        da.iota(7)
        da.rearrange(3)
        da.setInfoOnComponent(0,"X [m]")
        da.setInfoOnComponent(1,"Y [m]")
        da.setInfoOnComponent(2,"Z [km]")
        da2=da[5:-1]
        self.assertEqual([22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41., 42., 43., 44., 45., 46., 47., 48., 49., 50., 51.],da2.getValues())
        da2=da[4]
        self.assertEqual([19., 20., 21],da2.getValues())
        try:
            da2=da[4:17]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        da2=da[5:-2,2]
        self.assertEqual([24., 27., 30., 33., 36., 39., 42., 45., 48.],da2.getValues())
        da2=da[5:8,:]
        self.assertEqual([22., 23., 24., 25., 26., 27., 28., 29., 30.],da2.getValues())
        da2=da[:]
        self.assertTrue(da2.isEqual(da,1e-12))
        da2=da[:,:]
        self.assertTrue(da2.isEqual(da,1e-12))
        try:
            da2=da[:,:,:]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        self.assertTrue(da[5:8,-2].isEqualWithoutConsideringStr(DataArrayDouble([23.,26.,29.]),1e-12))
        da2=da[5:8,:-2]
        self.assertEqual([22., 25., 28.],da2.getValues())
        try:
            da2=da[5:-18,2]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        da2=da[5:5,2]
        self.assertEqual([],da2.getValues())
        pass

    def testSwigSetItem1(self):
        da=DataArrayInt.New()
        da.alloc(20,1)
        da.iota(7)
        da.rearrange(5)
        da.setInfoOnComponent(0,"X [m]") ; da.setInfoOnComponent(1,"Y [km]") ; da.setInfoOnComponent(2,"Y [m]")
        da.setInfoOnComponent(3,"Z [W]") ; da.setInfoOnComponent(4,"ZZ [km]") ; 
        da[:,2]=3
        self.assertEqual([7, 8, 3, 10, 11, 12, 13, 3, 15, 16, 17, 18, 3, 20, 21, 22, 23, 3, 25, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[2]=3
        self.assertEqual([7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 3, 3, 3, 3, 3, 22, 23, 24, 25, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[[0,3]]=-1
        self.assertEqual([-1, -1, -1, -1, -1, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, -1, -1, -1, -1, -1],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[:,[1,3,4]]=-3
        self.assertEqual([7, -3, 9, -3, -3, 12, -3, 14, -3, -3, 17, -3, 19, -3, -3, 22, -3, 24, -3, -3],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da2=DataArrayInt.New() ; da2.setValues([0,2,3],3,1)
        da[da2]=-7
        self.assertEqual([-7, -7, -7, -7, -7, 12, 13, 14, 15, 16, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,-2:]=-7
        self.assertEqual([7, 8, 9, -7, -7, 12, 13, 14, 15, 16, 17, 18, 19, -7, -7, 22, 23, 24, -7, -7],da.getValues())
        # Let's test with DAI right hand side
        da1=DataArrayInt.New()
        da1.setValues([25,26,27,125,126,127],2,3)
        #
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[-2:,1:4]=da1
        self.assertEqual([7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 25, 26, 27, 21, 22, 125, 126, 127, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[1:,3]=[225,226,227]
        self.assertEqual([7, 8, 9, 10, 11, 12, 13, 14, 225, 16, 17, 18, 19, 226, 21, 22, 23, 24, 227, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[1,2:]=[225,226,227]
        self.assertEqual([7, 8, 9, 10, 11, 12, 13, 225, 226, 227, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,-2:]=[88,99,1010,1111,1212,1313]
        self.assertEqual([7, 8, 9, 88, 99, 12, 13, 14, 15, 16, 17, 18, 19, 1010, 1111, 22, 23, 24, 1212, 1313],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da3=DataArrayInt.New(); da3.setValues([88,99,1010,1111,1212,1313],3,2)
        da[da2,-2:]=da3
        self.assertEqual([7, 8, 9, 88, 99, 12, 13, 14, 15, 16, 17, 18, 19, 1010, 1111, 22, 23, 24, 1212, 1313],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,[0,2]]=da3
        self.assertEqual([88, 8, 99, 10, 11, 12, 13, 14, 15, 16, 1010, 18, 1111, 20, 21, 1212, 23, 1313, 25, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,0:3:2]=da3
        self.assertEqual([88, 8, 99, 10, 11, 12, 13, 14, 15, 16, 1010, 18, 1111, 20, 21, 1212, 23, 1313, 25, 26],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,0:3:2]=-8
        self.assertEqual([-8, 8, -8, 10, 11, 12, 13, 14, 15, 16, -8, 18, -8, 20, 21, -8, 23, -8, 25, 26],da.getValues())
        pass

    def testSwigSetItem2(self):
        da=DataArrayDouble.New()
        da.alloc(20,1)
        da.iota(7)
        da.rearrange(5)
        da.setInfoOnComponent(0,"X [m]") ; da.setInfoOnComponent(1,"Y [km]") ; da.setInfoOnComponent(2,"Y [m]")
        da.setInfoOnComponent(3,"Z [W]") ; da.setInfoOnComponent(4,"ZZ [km]") ; 
        da[:,2]=3.
        self.assertEqual([7., 8., 3., 10., 11., 12., 13., 3., 15., 16., 17., 18., 3., 20., 21., 22., 23., 3., 25., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[2]=3.
        self.assertEqual([7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 3., 3., 3., 3., 3., 22., 23., 24., 25., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[[0,3]]=-1.
        self.assertEqual([-1., -1., -1., -1., -1., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., -1., -1., -1., -1., -1.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[:,[1,3,4]]=-3.
        self.assertEqual([7., -3., 9., -3., -3., 12., -3., 14., -3., -3., 17., -3., 19., -3., -3., 22., -3., 24., -3., -3.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da2=DataArrayInt.New() ; da2.setValues([0,2,3],3,1)
        da[da2]=-7.
        self.assertEqual([-7., -7., -7., -7., -7., 12., 13., 14., 15., 16., -7., -7., -7., -7., -7., -7., -7., -7., -7., -7.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,-2:]=-7
        self.assertEqual([7., 8., 9., -7., -7., 12., 13., 14., 15., 16., 17., 18., 19., -7., -7., 22., 23., 24., -7., -7.],da.getValues())
        # Let's test with DAI right hand side
        da1=DataArrayDouble.New()
        da1.setValues([25,26,27,125,126,127],2,3)
        #
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[-2:,1:4]=da1
        self.assertEqual([7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 25., 26., 27., 21., 22., 125., 126., 127., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[1:,3]=[225.,226.,227.]
        self.assertEqual([7., 8., 9., 10., 11., 12., 13., 14., 225., 16., 17., 18., 19., 226., 21., 22., 23., 24., 227., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[1,2:]=[225,226,227]
        self.assertEqual([7., 8., 9., 10., 11., 12., 13., 225., 226., 227., 17., 18., 19., 20., 21., 22., 23., 24., 25., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,-2:]=[88,99,1010,1111,1212,1313]
        self.assertEqual([7., 8., 9., 88., 99., 12., 13., 14., 15., 16., 17., 18., 19., 1010., 1111., 22., 23., 24., 1212., 1313.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da3=DataArrayDouble.New(); da3.setValues([88,99,1010,1111,1212,1313],3,2)
        da[da2,-2:]=da3
        self.assertEqual([7., 8., 9., 88., 99., 12., 13., 14., 15., 16., 17., 18., 19., 1010., 1111., 22., 23., 24., 1212., 1313.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,[0,2]]=da3
        self.assertEqual([88., 8., 99., 10., 11., 12., 13., 14., 15., 16., 1010., 18., 1111., 20., 21., 1212., 23., 1313., 25., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,0:3:2]=da3
        self.assertEqual([88., 8., 99., 10., 11., 12., 13., 14., 15., 16., 1010., 18., 1111., 20., 21., 1212., 23., 1313., 25., 26.],da.getValues())
        da.rearrange(1) ; da.iota(7) ; da.rearrange(5)
        da[da2,0:3:2]=-8.
        self.assertEqual([-8., 8., -8., 10., 11., 12., 13., 14., 15., 16., -8., 18., -8., 20., 21., -8., 23., -8., 25., 26.],da.getValues())
        pass

    def testSwigDADOp(self):
        da=DataArrayDouble.New()
        da.alloc(12,1)
        da.iota(7.)
        da1=DataArrayDouble.New()
        da1.alloc(12,1)
        da1.iota(8.)
        da2=da+da1
        self.assertEqual([15., 17., 19., 21., 23., 25., 27., 29., 31., 33., 35., 37.],da2.getValues())
        da2=da+3
        da3=3+da
        self.assertTrue(da2.isEqual(da3,1e-12))
        da2=da-1.
        self.assertEqual([6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0],da2.getValues())
        da2=1-da
        self.assertEqual([-6.0, -7.0, -8.0, -9.0, -10.0, -11.0, -12.0, -13.0, -14.0, -15.0, -16.0, -17.0],da2.getValues())
        da2=da*3
        self.assertEqual([21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0],da2.getValues())
        da2=3.*da
        self.assertEqual([21.0, 24.0, 27.0, 30.0, 33.0, 36.0, 39.0, 42.0, 45.0, 48.0, 51.0, 54.0],da2.getValues())
        da2=da*da1
        self.assertEqual([56.0, 72.0, 90.0, 110.0, 132.0, 156.0, 182.0, 210.0, 240.0, 272.0, 306.0, 342.0],da2.getValues())
        da2=da/4.
        self.assertEqual([1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5],da2.getValues())
        da3=4./da
        da4=da3*da2
        self.assertTrue(da4.isUniform(1.,1e-12))
        st1=da.getHiddenCppPointer()
        da+=1
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertTrue(da.isEqual(da1,1e-12))
        da-=8
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual(range(12),da.getValues())
        da+=da1
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0],da.getValues())
        da*=0.5
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0],da.getValues())
        da*=da1
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([32.0, 45.0, 60.0, 77.0, 96.0, 117.0, 140.0, 165.0, 192.0, 221.0, 252.0, 285.0],da.getValues())
        da/=da1
        self.assertEqual(st1,st2)
        self.assertEqual([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0],da.getValues())
        da/=2
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5],da.getValues())
        da.rearrange(3)
        da5=DataArrayDouble.New()
        da5.setValues([5.,4.,3.,2.],4,1)
        da*=da5 # it works with unmathing number of compo
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([10.0, 12.5, 15.0, 14.0, 16.0, 18.0, 15.0, 16.5, 18.0, 13.0, 14.0, 15.0],da.getValues())
        #
        da.alloc(30,1)
        da.iota(7.)
        da.rearrange(3)
        ids=DataArrayInt.New()
        ids.setValues([3,4,7],3,1)
        da[ids,:]=[5.,8.,9.]
        self.assertEqual([7.,8.,9.,10.,11.,12.,13.,14.,15.,5.,8.,9.,5.,8.,9.,22.,23.,24.,25.,26.,27.,5.,8.,9.,31.,32.,33.,34.,35.,36.0],da.getValues())
        #
        da.rearrange(1) ; da.iota(7) ; da.rearrange(3)
        da[ids,[1,2]]=[5,8]
        self.assertEqual([7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,5.,8.,19.,5.,8.,22.,23.,24.,25.,26.,27.,28.,5.,8.,31.,32.,33.,34.,35.,36.],da.getValues())
        pass

    def testSwigDAIOp(self):
        da=DataArrayInt.New()
        da.alloc(12,1)
        da.iota(7)
        da1=DataArrayInt.New()
        da1.alloc(12,1)
        da1.iota(8)
        da2=da+da1
        self.assertEqual([15,17,19,21,23,25,27,29,31,33,35,37],da2.getValues())
        da2=da+3
        da3=3+da
        self.assertTrue(da2.isEqual(da3))
        da2=da-1
        self.assertEqual([6,7,8,9,10,11,12,13,14,15,16,17],da2.getValues())
        da2=1-da
        self.assertEqual([-6,-7,-8,-9,-10,-11,-12,-13,-14,-15,-16,-17],da2.getValues())
        da2=da*3
        self.assertEqual([21,24,27,30,33,36,39,42,45,48,51,54.0],da2.getValues())
        da2=3*da
        self.assertEqual([21,24,27,30,33,36,39,42,45,48,51,54.0],da2.getValues())
        da2=da*da1
        self.assertEqual([56,72,90,110,132,156,182,210,240,272,306,342.0],da2.getValues())
        da2=da/4
        self.assertEqual([1,2,2,2,2,3,3,3,3,4,4,4],da2.getValues())
        da3=4/da
        da4=da3*da2
        self.assertTrue(da4.isUniform(0))
        st1=da.getHiddenCppPointer()
        da+=1
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertTrue(da.isEqual(da1))
        da-=8
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual(range(12),da.getValues())
        da+=da1
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([8,10,12,14,16,18,20,22,24,26,28,30],da.getValues())
        da/=2
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([4,5,6,7,8,9,10,11,12,13,14,15],da.getValues())
        da*=da1
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([32,45,60,77,96,117,140,165,192,221,252,285],da.getValues())
        da/=da1
        self.assertEqual(st1,st2)
        self.assertEqual([4,5,6,7,8,9,10,11,12,13,14,15],da.getValues())
        da/=2
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([2,2, 3,3, 4,4, 5,5, 6,6, 7,7],da.getValues())
        da.rearrange(3)
        da5=DataArrayInt.New()
        da5.setValues([5,4,3,2],4,1)
        da*=da5 # it works with unmathing number of compo
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([10,10, 15,12,16,16,15,15, 18,12,14,14],da.getValues())
        da%=6
        st2=da.getHiddenCppPointer()
        self.assertEqual(st1,st2)
        self.assertEqual([4,4,3,0,4,4,3,3,0,0,2,2],da.getValues())
        #
        da.alloc(30,1)
        da.iota(7)
        da.rearrange(3)
        ids=DataArrayInt.New()
        ids.setValues([3,4,7],3,1)
        da[ids,:]=[5,8,9]
        self.assertEqual([7,8,9,10,11,12,13,14,15,5,8,9,5,8,9,22,23,24,25,26,27,5,8,9,31,32,33,34,35,36],da.getValues())
        #
        da.rearrange(1) ; da.iota(7) ; da.rearrange(3)
        da[ids,[1,2]]=[5,8]
        self.assertEqual([7,8,9,10,11,12,13,14,15,16,5,8,19,5,8,22,23,24,25,26,27,28,5,8,31,32,33,34,35,36],da.getValues())
        pass

    def testSwigDAIOp2(self):
        da=DataArrayInt.New()
        st=da.getHiddenCppPointer()
        da.alloc(10,3)
        da.rearrange(1)
        da.iota(0)
        da.rearrange(3)
        da[:,1]+=4
        da[-2:,2]+=10
        da[-2:,2]+=10
        da[:,2]+=da[:,0]
        da[da[0],:]=7
        self.assertEqual(st,da.getHiddenCppPointer())
        self.assertEqual(da.getValues(),[7,7,7,3,8,8,7,7,7,9,14,20,12,17,26,7,7,7,18,23,38,21,26,44,24,29,70,27,32,76])
        pass

    def testSwigDAIOp3(self):
        da=DataArrayInt.New()
        self.assertRaises(InterpKernelException,da.__len__)
        self.assertRaises(InterpKernelException,da.__int__)
        for elt in da:
            self.assertTrue(False)
            pass
        da.alloc(12,3)
        da.rearrange(1) ; da.fillWithZero()
        l1=list(da)
        self.assertEqual(36,len(da));
        da.rearrange(3)
        tmp=da[0]
        self.assertRaises(InterpKernelException,tmp.__int__)
        self.assertEqual(12,len(da));
        l=list(da)
        for elt in enumerate(l):
            elt[1][2]=elt[0]
            pass
        ref=[0,0,0,0,0,1,0,0,2,0,0,3,0,0,4,0,0,5,0,0,6,0,0,7,0,0,8,0,0,9,0,0,10,0,0,11]
        self.assertEqual(ref,da.getValues());
        da.rearrange(1)
        l=[int(elt) for elt in l1]
        self.assertEqual(ref,da.getValues());
        self.assertEqual(11,int(da[-1:]))
        pass

    def testSwigDADOp3(self):
        da=DataArrayDouble.New()
        self.assertRaises(InterpKernelException,da.__len__)
        self.assertRaises(InterpKernelException,da.__float__)
        for elt in da:
            self.assertTrue(False)
            pass
        da.alloc(12,3)
        da.rearrange(1) ; da.fillWithZero()
        l1=list(da)
        self.assertEqual(36,len(da));
        da.rearrange(3)
        tmp=da[0]
        self.assertRaises(InterpKernelException,tmp.__float__)
        self.assertEqual(12,len(da));
        l=list(da)
        for elt in enumerate(l):
            elt[1][2]=elt[0]
            pass
        ref=[0.,0.,0.,0.,0.,1.,0.,0.,2.,0.,0.,3.,0.,0.,4.,0.,0.,5.,0.,0.,6.,0.,0.,7.,0.,0.,8.,0.,0.,9.,0.,0.,10.,0.,0.,11.]
        self.assertEqual(ref,da.getValues());
        da.rearrange(1)
        l=[float(elt) for elt in l1]
        self.assertEqual(ref,da.getValues());
        self.assertEqual(11.,float(da[-1:]))
        pass

    def testSwigDataArrayIntIterator1(self):
        da=DataArrayInt.New()
        da.alloc(12,1)
        da.iota(2)
        da.rearrange(3)
        # __getitem__ testing
        li=[]
        for it in da:
            li+=it[1:]
            pass
        self.assertEqual([3, 4, 6, 7, 9, 10, 12, 13],li)
        li=[]
        for it in da:
            li+=[it[-1]]
            pass
        self.assertEqual([4, 7, 10, 13],li)
        li=[]
        for it in da:
            li+=it[[2,1,0]]
            pass
        self.assertEqual([4, 3, 2, 7, 6, 5, 10, 9, 8, 13, 12, 11],li)
        # __setitem__ testing
        da3=da.deepCpy()
        da2=DataArrayInt.New()
        da2.alloc(12,1)
        da2.iota(2002)
        da2.rearrange(3)
        it2=da2.__iter__()
        i=0
        for it in da:
            pt=it2.next()
            it[:]=pt
            pass
        self.assertTrue(da.isEqual(da2))
        da=da3
        da3=da.deepCpy()
        #
        for it in da:
            it[:]=5
            pass
        da.rearrange(1)
        self.assertTrue(da.isUniform(5))
        da=da3
        da3=da.deepCpy()
        #
        for it in da:
            it[:]=[8,9,12]
            pass
        self.assertEqual([8, 9, 12, 8, 9, 12, 8, 9, 12, 8, 9, 12],da.getValues())
        da=da3
        da3=da.deepCpy()
        #
        for it in da:
            it[2]=[7]
            pass
        self.assertEqual([2, 3, 7, 5, 6, 7, 8, 9, 7, 11, 12, 7],da.getValues())
        pass

    def testSwigDataArrayDoubleIterator1(self):
        da=DataArrayDouble.New()
        da.alloc(12,1)
        da.iota(2)
        da.rearrange(3)
        # __getitem__ testing
        li=[]
        for it in da:
            li+=it[1:]
            pass
        self.assertEqual([3, 4, 6, 7, 9, 10, 12, 13],li)
        li=[]
        for it in da:
            li+=[it[-1]]
            pass
        self.assertEqual([4, 7, 10, 13],li)
        li=[]
        for it in da:
            li+=it[[2,1,0]]
            pass
        self.assertEqual([4, 3, 2, 7, 6, 5, 10, 9, 8, 13, 12, 11],li)
        # __setitem__ testing
        da3=da.deepCpy()
        da2=DataArrayDouble.New()
        da2.alloc(12,1)
        da2.iota(2002)
        da2.rearrange(3)
        it2=da2.__iter__()
        i=0
        for it in da:
            pt=it2.next()
            it[:]=pt
            pass
        self.assertTrue(da.isEqual(da2,1e-12))
        da=da3
        da3=da.deepCpy()
        #
        for it in da:
            it[:]=5
            pass
        da.rearrange(1)
        self.assertTrue(da.isUniform(5,1e-12))
        da=da3
        da3=da.deepCpy()
        #
        for it in da:
            it[:]=[8,9,12]
            pass
        self.assertEqual([8, 9, 12, 8, 9, 12, 8, 9, 12, 8, 9, 12],da.getValues())
        da=da3
        da3=da.deepCpy()
        #
        for it in da:
            it[2]=[7]
            pass
        self.assertEqual([2, 3, 7, 5, 6, 7, 8, 9, 7, 11, 12, 7],da.getValues())
        pass

    def testSwigUMeshIterator1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        li1=[]
        li2=[]
        for cell in m:
            li1+=cell.getAllConn()[1:]
            li2+=[cell.getType()]
            pass
        self.assertEqual(li1,[0, 3, 4, 1, 1, 4, 2, 4, 5, 2, 6, 7, 4, 3, 7, 8, 5, 4])
        self.assertEqual(li2,[4, 3, 3, 4, 4])
        pass

    def testSwigUMeshIterator2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        self.assertRaises(InterpKernelException,m.cellsByType);
        m.rearrange2ConsecutiveCellTypes()
        li1=[]
        li2=[]
        li3=[]
        for cellsByType in m.cellsByType():
            li1.append(cellsByType.getType())
            li2.append(cellsByType.getNumberOfElems())
            temp=[]
            for cell in cellsByType:
                t=[None,None]
                t[0]=cell.getType()
                t[1]=cell.getAllConn()[1:]
                temp.append(t)
                pass
            li3.append(temp)
            pass
        self.assertEqual(li1,[4, 3])
        self.assertEqual(li2,[3, 2])
        self.assertEqual(li3,[[[4, (0, 3, 4, 1)], [4, (6, 7, 4, 3)], [4, (7, 8, 5, 4)]], [[3, (1, 4, 2)], [3, (4, 5, 2)]]])
        pass

    def testDAIAggregateMulti1(self):
        a=DataArrayInt.New()
        a.setValues(range(4),2,2)
        a.setName("aa")
        b=DataArrayInt.New()
        b.setValues(range(6),3,2)
        c=DataArrayInt.Aggregate([a,b])
        self.assertEqual(range(4)+range(6),c.getValues())
        self.assertEqual("aa",c.getName())
        self.assertEqual(5,c.getNumberOfTuples())
        self.assertEqual(2,c.getNumberOfComponents())
        pass

    def testMergeUMeshes2(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m3=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        #
        vec1=[0,2,3]
        m2_2=m2.buildPartOfMySelf(vec1,False);
        vec2=[1,1]
        m3_2=m3.buildPartOfMySelf(vec2,False);
        #
        ms=[m1,m2_2,m3_2];
        #
        self.assertRaises(InterpKernelException,MEDCouplingUMesh.MergeUMeshes,ms+[None]);
        self.assertRaises(InterpKernelException,MEDCouplingUMesh.MergeUMeshes,ms+[3.4])
        m4=MEDCouplingUMesh.MergeUMeshes(ms);
        m4.checkCoherency();
        self.assertEqual(10,m4.getNumberOfCells());
        self.assertEqual(20,m4.getNumberOfNodes());
        self.assertEqual(45,m4.getMeshLength());
        m4bis=MEDCouplingMesh.MergeMeshes(ms);
        self.assertTrue(m4.isEqual(m4bis,1e-12))
        del m4bis
        #
        vec3=[0,1,2,3,4]
        m4_1=m4.buildPartOfMySelf(vec3,False);
        m4_1.setName(m1.getName());
        self.assertTrue(m4_1.isEqual(m1,1e-12));
        #
        vec4=[5,6,7]
        m4_2=m4.buildPartOfMySelf(vec4,False);
        cellCor,nodeCor=m4_2.checkGeoEquivalWith(m2_2,10,1e-12);
        #
        vec5=[8,9]
        m4_3=m4.buildPartOfMySelf(vec5,False);
        self.assertEqual(2,m4_3.getNumberOfCells());
        self.assertEqual(3,m4_3.getNumberOfNodes());
        m3_2.zipCoords();
        m4_3.setName(m3_2.getName());
        self.assertTrue(m4_3.isEqual(m3_2,1e-12));
        #
        pass

    def testBuild0DMeshFromCoords1(self):
        sourceCoords=[-0.3,-0.3,0., 0.7,-0.3,0., -0.3,0.7,0., 0.7,0.7,0.]
        coo=DataArrayDouble.New();
        coo.setValues(sourceCoords,4,3);
        coo.setName("My0D");
        m=MEDCouplingUMesh.Build0DMeshFromCoords(coo);
        m.checkCoherency();
        self.assertEqual(4,m.getNumberOfNodes());
        self.assertEqual(4,m.getNumberOfCells());
        self.assertEqual(3,m.getSpaceDimension());
        self.assertEqual(0,m.getMeshDimension());
        types1=m.getAllGeoTypes();
        self.assertEqual([NORM_POINT1],types1);
        for i in xrange(4):
            conn=m.getNodeIdsOfCell(i);
            self.assertEqual([i],conn);
            self.assertTrue(NORM_POINT1==m.getTypeOfCell(i));
            pass
        self.assertEqual(m.getName(),"My0D");
        pass

    def testDescriptionInMeshTimeUnit1(self):
        text1="totoTTEDD";
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        m.setDescription(text1);
        self.assertEqual(m.getDescription(),text1);
        m2=m.deepCpy();
        self.assertTrue(m.isEqual(m2,1e-12));
        self.assertEqual(m2.getDescription(),text1);
        m2.setDescription("ggg");
        self.assertTrue(not m.isEqual(m2,1e-12));
        #
        f=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f.setTimeUnit(text1);
        self.assertEqual(f.getTimeUnit(),text1);
        f2=f.deepCpy();
        self.assertEqual(f2.getTimeUnit(),text1);
        #
        pass

    def testMultiFields1(self):
        mfs=MEDCouplingDataForTest.buildMultiFields_1();
        ms=mfs.getMeshes();
        dms,refs=mfs.getDifferentMeshes()
        das=mfs.getArrays();
        das2,refs2=mfs.getDifferentArrays()
        self.assertEqual(5,len(mfs.getFields()))
        self.assertEqual(1,len(mfs.getFields()[0].getArrays()));
        self.assertEqual(2,len(mfs.getFields()[1].getArrays()));
        self.assertEqual(1,len(mfs.getFields()[2].getArrays()));
        self.assertEqual(1,len(mfs.getFields()[3].getArrays()));
        self.assertEqual(1,len(mfs.getFields()[4].getArrays()));
        self.assertEqual(5,len(ms));
        self.assertEqual(2,len(dms));
        self.assertEqual(6,len(das));
        self.assertEqual(5,len(das2));
        mfs2=mfs.deepCpy();
        self.assertTrue(mfs.isEqual(mfs2,1e-12,1e-12))
        pass

    def testFieldOverTime1(self):
        fs=MEDCouplingDataForTest.buildMultiFields_2();
        self.assertRaises(InterpKernelException,MEDCouplingFieldOverTime.New,fs);
        f4bis=fs[4].buildNewTimeReprFromThis(ONE_TIME,False);
        fs[4]=f4bis;
        self.assertRaises(InterpKernelException,MEDCouplingFieldOverTime.New,fs);
        f4bis.setTime(2.7,20,21);
        fot=MEDCouplingFieldOverTime.New(fs);
        dt=fot.getDefinitionTimeZone();
        hs=dt.getHotSpotsTime();
        self.assertEqual(6,len(hs));
        expected1=[0.2,0.7,1.2,1.35,1.7,2.7]
        for i in xrange(6):
            self.assertAlmostEqual(expected1[i],hs[i],12);
            pass
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeRight(0.2);
        self.assertEqual(0,meshId);
        self.assertEqual(0,arrId);
        self.assertEqual(0,arrIdInField);
        self.assertEqual(0,fieldId);
        #
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeRight(0.7);
        self.assertEqual(0,meshId);
        self.assertEqual(1,arrId);
        self.assertEqual(0,arrIdInField);
        self.assertEqual(1,fieldId);
        #
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeLeft(1.2);#**** WARNING left here
        self.assertEqual(0,meshId);
        self.assertEqual(2,arrId);
        self.assertEqual(1,arrIdInField);
        self.assertEqual(1,fieldId);
        #
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeRight(1.2);#**** WARNING right again here
        self.assertEqual(1,meshId);
        self.assertEqual(3,arrId);
        self.assertEqual(0,arrIdInField);
        self.assertEqual(2,fieldId);
        #
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeRight(1.35);
        self.assertEqual(1,meshId);
        self.assertEqual(3,arrId);
        self.assertEqual(0,arrIdInField);
        self.assertEqual(2,fieldId);
        #
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeRight(1.7);
        self.assertEqual(0,meshId);
        self.assertEqual(3,arrId);
        self.assertEqual(0,arrIdInField);
        self.assertEqual(3,fieldId);
        #
        meshId,arrId,arrIdInField,fieldId=dt.getIdsOnTimeRight(2.7);
        self.assertEqual(1,meshId);
        self.assertEqual(4,arrId);
        self.assertEqual(0,arrIdInField);
        self.assertEqual(4,fieldId);
        #
        dt2=MEDCouplingDefinitionTime();
        self.assertTrue(not dt2.isEqual(dt));
        dt2.assign(dt);
        dt2.assign(dt);#to check memory management
        self.assertTrue(dt2.isEqual(dt));
        #
        dt3=MEDCouplingDefinitionTime();
        #
        pass

    def testDAICheckAndPreparePermutation1(self):
        vals1=[9,10,0,6,4,11,3,7];
        expect1=[5,6,0,3,2,7,1,4];
        vals2=[9,10,0,6,10,11,3,7];
        da=DataArrayInt.New();
        da.setValues(vals1,8,1);
        da2=da.checkAndPreparePermutation();
        self.assertEqual(8,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        for i in xrange(8):
            self.assertEqual(expect1[i],da2.getIJ(i,0));
            pass
        #
        da=DataArrayInt.New();
        da.alloc(8,1);
        da.iota(0);
        da2=da.checkAndPreparePermutation();
        self.assertEqual(8,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        self.assertTrue(da2.isIdentity());
        #
        da=DataArrayInt.New();
        da.alloc(8,1);
        da.setValues(vals2,8,1);
        self.assertRaises(InterpKernelException,da.checkAndPreparePermutation);
        pass

    def testDAIChangeSurjectiveFormat1(self):
        vals1=[0,3,2,3,2,2,1,2]
        expected1=[0,1,2,6,8]
        expected2=[0,  6,  2,4,5,7,  1,3]
        da=DataArrayInt.New();
        da.setValues(vals1,8,1);
        #
        da2,da2I=da.changeSurjectiveFormat(4);
        self.assertEqual(5,da2I.getNumberOfTuples());
        self.assertEqual(8,da2.getNumberOfTuples());
        self.assertEqual(expected1,da2I.getValues());
        self.assertEqual(expected2,da2.getValues());
        #
        self.assertRaises(InterpKernelException,da.changeSurjectiveFormat,3);
        #
        pass

    def testUMeshGetCellIdsLyingOnNodes1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        nodeIds1=[1,2,3,4,6]
        nodeIds2=[6,7]
        da=m.getCellIdsLyingOnNodes(nodeIds1,True);
        self.assertEqual(1,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        self.assertEqual(1,da.getIJ(0,0));
        da2=DataArrayInt.New()
        da2.setValues(nodeIds2,2,1)
        da=m.getCellIdsLyingOnNodes(da2,False);
        self.assertEqual(2,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        self.assertEqual(3,da.getIJ(0,0));
        self.assertEqual(4,da.getIJ(1,0));
        pass

    def testUMeshFindCellIdsOnBoundary1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        da5=m.findCellIdsOnBoundary();
        self.assertEqual(5,da5.getNumberOfTuples());
        self.assertTrue(da5.isIdentity());
        pass

    def testMeshSetTime1(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        #
        self.assertTrue(m1.isEqual(m2,1e-12));
        m1.setTime(3.14,6,7);
        tmp3,tmp1,tmp2=m1.getTime();
        self.assertEqual(6,tmp1);
        self.assertEqual(7,tmp2);
        self.assertAlmostEqual(3.14,tmp3,12);
        self.assertTrue(not m1.isEqual(m2,1e-12));
        m2.setTime(3.14,6,7);
        self.assertTrue(m1.isEqual(m2,1e-12));
        m1.setTimeUnit("ms");
        self.assertTrue(m1.getTimeUnit()=="ms");
        m1.setTimeUnit("us");
        self.assertTrue(m1.getTimeUnit()=="us");
        self.assertTrue(not m1.isEqual(m2,1e-12));
        m2.setTimeUnit("us");
        self.assertTrue(m1.isEqual(m2,1e-12));
        m2.setTime(3.14,6,8);
        self.assertTrue(not m1.isEqual(m2,1e-12));
        m2.setTime(3.14,7,7);
        self.assertTrue(not m1.isEqual(m2,1e-12));
        m2.setTime(3.15,6,7);
        self.assertTrue(not m1.isEqual(m2,1e-12));
        #
        m1.setTime(10.34,55,12);
        m3=m1.deepCpy();
        self.assertTrue(m1.isEqual(m3,1e-12));
        tmp3,tmp1,tmp2=m3.getTime();
        self.assertEqual(55,tmp1);
        self.assertEqual(12,tmp2);
        self.assertAlmostEqual(10.34,tmp3,12);
        #
        # testing CMesh
        coo1=[0.,1.,2.,3.5]
        a=DataArrayDouble.New();
        a.setValues(coo1,4,1);
        b=MEDCouplingCMesh.New();
        b.setCoordsAt(0,a);
        #
        b.setTime(5.67,8,100);
        tmp3,tmp1,tmp2=b.getTime();
        self.assertEqual(8,tmp1);
        self.assertEqual(100,tmp2);
        self.assertAlmostEqual(5.67,tmp3,12);
        c=b.deepCpy();
        self.assertTrue(c.isEqual(b,1e-12));
        tmp3,tmp1,tmp2=c.getTime();
        self.assertEqual(8,tmp1);
        self.assertEqual(100,tmp2);
        self.assertAlmostEqual(5.67,tmp3,12);
        pass

    def testApplyFuncTwo1(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setMesh(m1);
        #
        vals=[1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.,5.,15.,25.]
        da=DataArrayDouble.New();
        da.setValues(vals,5,3);
        f1.setArray(da);
        #
        self.assertRaises(InterpKernelException,da.applyFunc2,1,"y+z");
        da.setInfoOnComponent(0,"x [m]");
        da.setInfoOnComponent(1,"y [mm]");
        da.setInfoOnComponent(2,"z [km]");
        
        self.assertRaises(InterpKernelException, da.applyFunc2, 1, "x+y+zz+zzz");
        self.assertRaises(InterpKernelException, da.applyFunc2, 1, "toto(x+y)");
        self.assertRaises(InterpKernelException, da.applyFunc2, 1, "x/0");
        
        da2=da.applyFunc2(1,"y+z");
        self.assertEqual(1,da2.getNumberOfComponents());
        self.assertEqual(5,da2.getNumberOfTuples());
        expected1=[32.,34.,36.,38.,40.]
        for i in xrange(5):
            self.assertAlmostEqual(expected1[i],da2.getIJ(0,i),12);
            pass
        da2=da.applyFunc(1,"y+z");
        expected2=[12.,14.,16.,18.,20.]
        for i in xrange(5):
            self.assertAlmostEqual(expected2[i],da2.getIJ(0,i),12);
            pass
        #
        self.assertEqual(3,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        f1.applyFunc2(1,"y+z");
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(expected1[i],f1.getArray().getIJ(0,i),12);
            pass
        #
        pass

    def testApplyFuncThree1(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setMesh(m1);
        #
        vals=[1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.,5.,15.,25.]
        da=DataArrayDouble.New();
        da.setValues(vals,5,3);
        f1.setArray(da);
        #
        vs=3*[None];
        vs[0]="x"; vs[1]="Y"; vs[2]="z";
        self.assertRaises(InterpKernelException, da.applyFunc3, 1, vs, "y+z");
        self.assertRaises(InterpKernelException, da.applyFunc3, 1, vs, "x+Y+z+zz+zzz");
        self.assertRaises(InterpKernelException, da.applyFunc3, 1, vs, "x/0");
        vs[1]="y";
        da2=da.applyFunc3(1,vs,"y+z");
        expected1=[32.,34.,36.,38.,40.]
        for i in xrange(5):
            self.assertAlmostEqual(expected1[i],da2.getIJ(0,i),12);
            pass
        self.assertRaises(InterpKernelException, da.applyFunc3, 1, ["x","y","z","a"],"x+a")
        f1.setArray(da);
        self.assertEqual(3,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        f1.applyFunc3(1,vs,"y+z");
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(5,f1.getNumberOfTuples());
        for i in xrange(5):
            self.assertAlmostEqual(expected1[i],f1.getArray().getIJ(0,i),12);
            pass
        pass

    def testFillFromAnalyticTwo1(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m1.setTime(3.4,5,6); m1.setTimeUnit("us");
        self.assertRaises(InterpKernelException,m1.fillFromAnalytic2,ON_NODES,1,"y+z");
        m1.getCoords().setInfoOnComponent(0,"x [m]");
        m1.getCoords().setInfoOnComponent(1,"y");
        m1.getCoords().setInfoOnComponent(2,"z");
        f1=m1.fillFromAnalytic2(ON_NODES,1,"y+z");
        self.assertAlmostEqual(3.4,f1.getTime()[0],12) ; self.assertEqual(5,f1.getTime()[1]) ; self.assertEqual(6,f1.getTime()[2])
        self.assertEqual("us",f1.getTimeUnit())
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        expected1=[0.2, 0.7, 1.2, 0.7, 1.2, 1.7, 1.2, 1.7, 2.2]
        for i in xrange(9):
            self.assertAlmostEqual(expected1[i],f1.getArray().getIJ(0,i),12);
            pass
        pass

    def testFillFromAnalyticThree1(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m1.setTime(3.4,5,6); m1.setTimeUnit("us");
        vs=3*[None];
        vs[0]="x"; vs[1]="Y"; vs[2]="z";
        self.assertRaises(InterpKernelException,m1.fillFromAnalytic3,ON_NODES,1,vs,"y+z");
        vs[1]="y";
        f1=m1.fillFromAnalytic3(ON_NODES,1,vs,"y+z");
        self.assertAlmostEqual(3.4,f1.getTime()[0],12) ; self.assertEqual(5,f1.getTime()[1]) ; self.assertEqual(6,f1.getTime()[2])
        self.assertEqual("us",f1.getTimeUnit())
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        expected1=[0.2, 0.7, 1.2, 0.7, 1.2, 1.7, 1.2, 1.7, 2.2]
        for i in xrange(9):
            self.assertAlmostEqual(expected1[i],f1.getArray().getIJ(0,i),12);
            pass
        pass

    def testDAUnitVar1(self):
        da=DataArrayDouble.New();
        da.alloc(1,3);
        da.setInfoOnComponent(0,"XPS [m]");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="XPS");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2=="m");
        #
        da.setInfoOnComponent(0,"XPS         [m]");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="XPS");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2=="m");
        #
        da.setInfoOnComponent(0,"XPP         [m]");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="XPP");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2=="m");
        #
        da.setInfoOnComponent(0,"XPP kdep  kefer   [ m  ]");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="XPP kdep  kefer");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2==" m  ");
        #
        da.setInfoOnComponent(0,"     XPP k[  dep  k]efer   [ m^ 2/s^3*kJ  ]");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="     XPP k[  dep  k]efer");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2==" m^ 2/s^3*kJ  ");
        #
        da.setInfoOnComponent(0,"     XPP kefer   ");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="     XPP kefer   ");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2=="");
        #
        da.setInfoOnComponent(0,"temperature( bof)");
        st1=da.getVarOnComponent(0);
        self.assertTrue(st1=="temperature( bof)");
        st2=da.getUnitOnComponent(0);
        self.assertTrue(st2=="");
        #
        da.setInfoOnComponent(0,"kkk [m]");
        da.setInfoOnComponent(1,"ppp   [m^2/kJ]");
        da.setInfoOnComponent(2,"abcde   [MW/s]");
        #
        vs=da.getVarsOnComponent();
        self.assertEqual(3,len(vs));
        self.assertTrue(vs[0]=="kkk");
        self.assertTrue(vs[1]=="ppp");
        self.assertTrue(vs[2]=="abcde");
        vs=da.getUnitsOnComponent();
        self.assertEqual(3,len(vs));
        self.assertTrue(vs[0]=="m");
        self.assertTrue(vs[1]=="m^2/kJ");
        self.assertTrue(vs[2]=="MW/s");
        pass

    def testGaussCoordinates1(self):
        #Testing 1D cell types
        m1=MEDCouplingDataForTest.build1DMultiTypes_1();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME);
        f.setMesh(m1);
        wg1=[0.3];
        gsCoo1=[0.2];
        refCoo1=[-1.0,1.0];
        f.setGaussLocalizationOnType(NORM_SEG2,refCoo1,gsCoo1,wg1);
        wg2=wg1;
        gsCoo2=[0.2];
        refCoo2=[-1.0,1.0,0.0];
        f.setGaussLocalizationOnType(NORM_SEG3,refCoo2,gsCoo2,wg2);
        #
        resToTest=f.getLocalizationOfDiscr();
        self.assertEqual(3,resToTest.getNumberOfComponents());
        self.assertEqual(2,resToTest.getNumberOfTuples());
        expected1=[0.6,0.6,0.6, 0.6,0.6,0.6]
        for i in xrange(6):
            self.assertAlmostEqual(expected1[i],resToTest.getIJ(0,i),14);
            pass
        #
        #Testing 2D cell types
        m2=MEDCouplingDataForTest.build2DMultiTypes_1();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME);
        f.setMesh(m2);
        wg3=[0.3,0.3];
        tria3CooGauss=[ 0.1, 0.8, 0.2, 0.7 ]
        gsCoo3=tria3CooGauss
        tria3CooRef=[ 0.0, 0.0, 1.0 , 0.0, 0.0, 1.0 ]
        refCoo3=tria3CooRef;
        f.setGaussLocalizationOnType(NORM_TRI3,refCoo3,gsCoo3,wg3);
        wg4=[0.3,0.3,0.3];
        tria6CooGauss=[ 0.3, 0.2, 0.2, 0.1, 0.2, 0.4 ]
        gsCoo4=tria6CooGauss;
        tria6CooRef=[0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5]
        refCoo4=tria6CooRef;
        f.setGaussLocalizationOnType(NORM_TRI6,refCoo4,gsCoo4,wg4);
        wg5=[0.3,0.3,0.3,0.3];
        quad4CooGauss=[ 0.3, 0.2, 0.2, 0.1, 0.2, 0.4, 0.15, 0.27 ]
        gsCoo5=quad4CooGauss;
        quad4CooRef=[-1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0]
        refCoo5=quad4CooRef;
        f.setGaussLocalizationOnType(NORM_QUAD4,refCoo5,gsCoo5,wg5);
        wg6=[0.3,0.3,0.3,0.3];
        quad8CooGauss=[ 0.34, 0.16, 0.21, 0.3, 0.23, 0.4, 0.14, 0.37 ]
        gsCoo6=quad8CooGauss;
        quad8CooRef=[ -1.0, -1.0, 1.0, -1.0, 1.0,  1.0, -1.0,  1.0, 0.0, -1.0, 1.0,  0.0, 0.0,  1.0, -1.0,  0.0]
        refCoo6=quad8CooRef;
        f.setGaussLocalizationOnType(NORM_QUAD8,refCoo6,gsCoo6,wg6);
        #
        resToTest=f.getLocalizationOfDiscr();
        self.assertEqual(3,resToTest.getNumberOfComponents());
        self.assertEqual(13,resToTest.getNumberOfTuples());#2+3+4+4 gauss points for resp TRI3,TRI6,QUAD4,QUAD8
        expected2=[5.1,1.55,0.0, 4.7,1.65,0.0,
                   2.32,1.52,0.0, 1.6,1.32,0.0, 3.52,1.26,0.0,#TRI6
                   2.6,1.6,0.0, 2.4,1.8,0.0, 2.4,1.2,0.0, 2.3,1.46,0.0,#QUAD4
                   2.32,2.68,0.0, 2.6,2.42,0.0, 2.8,2.46,0.0, 2.74,2.28,0.0 ];#QUAD8
        for i in xrange(39):
            self.assertAlmostEqual(expected2[i],resToTest.getIJ(0,i),14);
            pass
        #
        #Testing 3D cell types
        m3=MEDCouplingDataForTest.build3DMultiTypes_1();
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME);
        f.setMesh(m3);
        #
        wg7=[0.3];
        tetra4CooGauss=[0.34, 0.16, 0.21]
        gsCoo7=tetra4CooGauss;
        tetra4CooRef=[0.0,1.0,0.0, 0.0,0.0,1.0, 0.0,0.0,0.0, 1.0,0.0,0.0]
        refCoo7=tetra4CooRef;
        f.setGaussLocalizationOnType(NORM_TETRA4,refCoo7,gsCoo7,wg7);
        wg8=[0.3];
        tetra10CooGauss=[0.2, 0.3, 0.1]
        gsCoo8=tetra10CooGauss;
        tetra10CooRef=[0.0,1.0,0.0, 0.0,0.0,0.0, 0.0,0.0,1.0, 1.0,0.0,0.0, 0.0,0.5,0.0, 0.0,0.0,0.5, 0.0,0.5,0.5, 0.5,0.5,0.0, 0.5,0.0,0.0, 0.5,0.0,0.5]
        refCoo8=tetra10CooRef;
        f.setGaussLocalizationOnType(NORM_TETRA10,refCoo8,gsCoo8,wg8);
        wg9=[0.3];
        pyra5CooGauss=[0.2, 0.3, 0.1]
        gsCoo9=pyra5CooGauss;
        pyra5CooRef=[1.0,0.0,0.0, 0.0,1.0,0.0, -1.0,0.0,0.0, 0.0,-1.0,0.0, 0.0,0.0,1.0]
        refCoo9=pyra5CooRef;
        f.setGaussLocalizationOnType(NORM_PYRA5,refCoo9,gsCoo9,wg9);
        wg10=[0.3];
        pyra13CooGauss=[0.1, 0.2, 0.7]
        gsCoo10=pyra13CooGauss;
        pyra13CooRef=[1.0,0.0,0.0, 0.0,1.0,0.0,-1.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.5,0.5,0.0,-0.5,0.5,0.0,-0.5,-0.5,0.0,0.5,-0.5,0.0,0.5,0.0,0.5,0.0,0.5,0.5,-0.5,0.0,0.5,0.0,-0.5,0.5]
        refCoo10=pyra13CooRef;
        f.setGaussLocalizationOnType(NORM_PYRA13,refCoo10,gsCoo10,wg10);
        wg11=[0.3];
        penta6CooGauss=[0.2, 0.3, 0.1]
        gsCoo11=penta6CooGauss;
        penta6CooRef=[-1.0,1.0,0.0,-1.0,-0.0,1.0,-1.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,0.0,0.0]
        refCoo11=penta6CooRef;
        f.setGaussLocalizationOnType(NORM_PENTA6,refCoo11,gsCoo11,wg11);
        wg12=[0.3];
        penta15CooGauss=[0.2, 0.3,0.15]
        gsCoo12=penta15CooGauss;
        penta15CooRef=[-1.0,1.0,0.0,-1.0,0.0,1.0,-1.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,0.0,0.0,-1.0,0.5,0.5,-1.0,0.0,0.5,-1.0,0.5,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.5,0.5,1.0,0.0, 0.5,1.0,0.5,0.0]
        refCoo12=penta15CooRef;
        f.setGaussLocalizationOnType(NORM_PENTA15,refCoo12,gsCoo12,wg12);
        wg13=[0.3];
        hexa8CooGauss=[0.2,0.3,0.15]
        gsCoo13=hexa8CooGauss;
        hexa8CooRef=[-1.0,-1.0,-1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,1.0,1.0]
        refCoo13=hexa8CooRef;
        f.setGaussLocalizationOnType(NORM_HEXA8,refCoo13,gsCoo13,wg13);
        wg14=[0.3];
        hexa20CooGauss=[0.11,0.3,0.55]
        gsCoo14=hexa20CooGauss;
        hexa20CooRef=[-1.0,-1.0,-1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,1.0,1.0,0.0,-1.0,-1.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,0.0,-1.0,-1.0,-1.0,0.0,1.0,-1.0,0.0,1.0,1.0,0.0,-1.0,1.0,0.0,0.0,-1.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,-1.0,0.0,1.0]
        refCoo14=hexa20CooRef;
        f.setGaussLocalizationOnType(NORM_HEXA20,refCoo14,gsCoo14,wg14);
        #
        resToTest=f.getLocalizationOfDiscr();
        self.assertEqual(3,resToTest.getNumberOfComponents());
        self.assertEqual(8,resToTest.getNumberOfTuples());#2+3+4+4 gauss points for resp TRI3,TRI6,QUAD4,QUAD8
        expected3=[1.312,3.15,1.02, 0.56,3.3,0.6, 2.18,1.1,0.2, 1.18,1.54,0.98, 1.56,0.3,3.6, 1.613,0.801,4.374, 2.6,2.4,2.3, 2.31232,2.3933985,1.553255]
        for i in xrange(24):
            self.assertAlmostEqual(expected3[i],resToTest.getIJ(0,i),14);
            pass
        #
        pass

    def testP2Localization1(self):
        m=MEDCouplingUMesh.New("testP2",2);
        coords=[0.,2.,3.5,0.,4.5,1.5,1.2,0.32,3.4,1.,2.1,2.4]
        conn=[0,1,2,3,4,5]
        coo=DataArrayDouble.New();
        coo.setValues(coords,6,2);
        m.setCoords(coo);
        m.allocateCells(1);
        m.insertNextCell(NORM_TRI6,6,conn[0:6])
        m.finishInsertingCells();
        #
        f=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
        f.setMesh(m);
        da=DataArrayDouble.New();
        vals1=[1.2,2.3,3.4, 2.2,3.3,4.4, 3.2,4.3,5.4, 4.2,5.3,6.4, 5.2,6.3,7.4, 6.2,7.3,8.4]
        da.setValues(vals1,6,3);
        f.setArray(da);
        #
        loc=[2.27,1.3]
        locs=f.getValueOnMulti(loc);
        expected1=[6.0921164547752236, 7.1921164547752232, 8.2921164547752255]
        for i in xrange(3):
            self.assertAlmostEqual(expected1[i],locs.getIJ(0,i),12);
            pass
        pass

    def testP2Localization2(self):
        m=MEDCouplingUMesh.New("testP2_2",3);
        coords=[0.33312787792955395, -0.35155740179580952, -0.03567564825034563, 1.307146326477638, -0.57234557776250305, -0.08608044208272235, 0.5551834466499993, 0.62324964668794192, -0.014638951108536295, 0.37761817224442129, -0.38324019806913578, 0.96283164472856886, 0.79494856035658679, -0.40628057809270046, 0.0021004190225864614, 1.023740446371799, 0.07665912970471335, -0.072889657161871096, 0.54564584619517376, 0.11132872093429744, 0.039647326652013051, 0.27164784387819052, -0.42018012100866675, 0.46563376500745146, 0.89501965094896418, -0.56148455362735061, 0.43337469695473035, 0.49118025152924394, 0.093884938060727313, 0.47216346905220891]
        conn=[0,1,2,3,4,5,6,7,8,9]
        coo=DataArrayDouble.New();
        coo.setValues(coords,10,3);
        m.setCoords(coo);
        m.allocateCells(1);
        m.insertNextCell(NORM_TETRA10,10,conn[0:10])
        m.finishInsertingCells();
        #
        f=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
        f.setMesh(m);
        da=DataArrayDouble.New();
        vals1=[1.1,2.1,3.1,4.1,5.2,6.2,7.2,8.2,9.2,10.2]
        da.setValues(vals1,10,1);
        f.setArray(da);
        #
        loc=[0.64637931739890486, -0.16185896817550552, 0.22678966365273748]
        locs=f.getValueOnMulti(loc);
        expected1=[10.0844021968047]
        for i in xrange(1):
            self.assertAlmostEqual(expected1[i],locs.getIJ(0,i),12);
            pass
        pass

    def testGetValueOn2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f.setMesh(m);
        arr=DataArrayDouble.New();
        nbOfCells=m.getNumberOfCells();
        f.setArray(arr);
        values1=[7.,107.,10007.,8.,108.,10008.,9.,109.,10009.,10.,110.,10010.,11.,111.,10011.]
        arr.setValues(values1,nbOfCells,3);
        loc=[-0.05,-0.05, 0.55,-0.25, 0.55,0.15, -0.05,0.45, 0.45,0.45]
        f.checkCoherency();
        locs=f.getValueOnMulti(loc);
        self.assertEqual(5,locs.getNumberOfTuples());
        self.assertEqual(3,locs.getNumberOfComponents());
        for j in xrange(15):
            self.assertAlmostEqual(values1[j],locs.getIJ(0,j),12);
            pass
        # Testing ON_NODES
        f=MEDCouplingFieldDouble.New(ON_NODES,NO_TIME);
        f.setMesh(m);
        arr=DataArrayDouble.New();
        nbOfNodes=m.getNumberOfNodes();
        f.setArray(arr);
        values2=[7.,107.,10007.,8.,108.,10008.,9.,109.,10009.,10.,110.,10010.,11.,111.,10011.,12.,112.,10012.,13.,113.,10013.,14.,114.,10014.,15.,115.,10015.]
        arr.setValues(values2,nbOfNodes,3);
        loc2=[0.5432,-0.2432, 0.5478,0.1528, 0.5432,-0.2432, 0.5432,-0.2432]
        expected2=[9.0272, 109.0272, 10009.0272, 11.4124,111.4124,10011.4124, 9.0272, 109.0272, 10009.0272, 9.0272, 109.0272, 10009.0272]
        f.checkCoherency();
        loc3=DataArrayDouble.New()
        loc3.setValues(loc2,4,2);
        locs=f.getValueOnMulti(loc3);
        self.assertEqual(4,locs.getNumberOfTuples());
        self.assertEqual(3,locs.getNumberOfComponents());
        for i in xrange(12):
            self.assertAlmostEqual(expected2[i],locs.getIJ(0,i),12);
            pass
        #
        pass

    def testDAIGetIdsNotEqual1(self):
        d=DataArrayInt.New();
        vals1=[2,3,5,6,8,5,5,6,1,-5]
        d.setValues(vals1,10,1);
        d2=d.getIdsNotEqual(5);
        self.assertEqual(7,d2.getNumberOfTuples());
        self.assertEqual(1,d2.getNumberOfComponents());
        expected1=[0,1,3,4,7,8,9]
        for i in xrange(7):
            self.assertEqual(expected1[i],d2.getIJ(0,i));
            pass
        d.rearrange(2);
        self.assertRaises(InterpKernelException,d.getIdsNotEqual,5);
        vals2=[-4,5,6]
        vals3=vals2;
        d.rearrange(1);
        d3=d.getIdsNotEqualList(vals3);
        self.assertEqual(5,d3.getNumberOfTuples());
        self.assertEqual(1,d3.getNumberOfComponents());
        expected2=[0,1,4,8,9]
        for i in xrange(5):
            self.assertEqual(expected2[i],d3.getIJ(0,i));
            pass
        pass

    def testDAIComputeOffsets1(self):
        d=DataArrayInt.New();
        vals1=[3,5,1,2,0,8]
        expected1=[0,3,8,9,11,11]
        d.setValues(vals1,6,1);
        d.computeOffsets();
        self.assertEqual(6,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        for i in xrange(6):
            self.assertEqual(expected1[i],d.getIJ(0,i));
            pass
        pass

    def testUMeshHexagonPrism1(self):
        coords=[0.8660254037844386, 0.5, 0.0, 0.0, 1.0, 0.0, -0.8660254037844386, 0.5, 0.0, -0.8660254037844386, -0.5, 0.0, 0.0, -1.0, 0.0, 0.8660254037844386, -0.5, 0.0,
                0.8660254037844386, 0.5, 2.0, 0.0, 1.0, 2.0, -0.8660254037844386, 0.5, 2.0, -0.8660254037844386, -0.5, 2.0, 0.0, -1.0, 2.0, 0.8660254037844386, -0.5, 2.0];
        conn=[1,2,3,4,5,0,7,8,9,10,11,6]
        mesh=MEDCouplingUMesh.New("MyFirstHexagonalPrism",3);
        coo=DataArrayDouble.New();
        coo.setValues(coords,12,3);
        mesh.setCoords(coo);
        mesh.allocateCells(1);
        mesh.insertNextCell(NORM_HEXGP12,12,conn[0:12])
        mesh.finishInsertingCells();
        #
        mesh.checkCoherency();
        vols=mesh.getMeasureField(False);
        self.assertEqual(1,vols.getNumberOfTuples());
        self.assertEqual(1,vols.getNumberOfComponents());
        self.assertAlmostEqual(-5.196152422706632,vols.getIJ(0,0),12);
        bary=mesh.getBarycenterAndOwner();
        self.assertEqual(1,bary.getNumberOfTuples());
        self.assertEqual(3,bary.getNumberOfComponents());
        expected1=[0.,0.,1.]
        for i in xrange(3):
            self.assertAlmostEqual(expected1[i],bary.getIJ(0,i),12);
            pass
        d1=DataArrayInt.New();
        d2=DataArrayInt.New();
        d3=DataArrayInt.New();
        d4=DataArrayInt.New();
        m2=mesh.buildDescendingConnectivity(d1,d2,d3,d4);
        self.assertEqual(8,m2.getNumberOfCells());
        expected4=[[1,2,3,4,5,0],[7,6,11,10,9,8],[1,7,8,2],[2,8,9,3],[3,9,10,4],[4,10,11,5],[5,11,6,0],[0,6,7,1]];
        expected2=[NORM_POLYGON, NORM_POLYGON, NORM_QUAD4, NORM_QUAD4, NORM_QUAD4, NORM_QUAD4, NORM_QUAD4, NORM_QUAD4];
        expected3=[6,6,4,4,4,4,4,4]
        for i in xrange(8):
            self.assertTrue(m2.getTypeOfCell(i)==expected2[i]);
            v=m2.getNodeIdsOfCell(i);
            self.assertTrue(len(v)==expected3[i]);
            self.assertEqual(expected4[i],v);
        #
        mesh.convertAllToPoly();
        self.assertTrue(NORM_POLYHED==mesh.getTypeOfCell(0));
        mesh.unPolyze();
        self.assertTrue(NORM_HEXGP12==mesh.getTypeOfCell(0));
        self.assertEqual(13,mesh.getMeshLength());
        #
        pass

    def testDADCheckIsMonotonic(self):
        da=DataArrayDouble.New();
        da.setValues([-1.,1.01,2.03,6.],2,2);
        self.assertRaises(InterpKernelException,da.isMonotonic,True,1e-12);
        da.rearrange(1);
        self.assertTrue(da.isMonotonic(True,1e-12));
        da.checkMonotonic(True,1e-12);
        da.setIJ(2,0,6.1);
        self.assertTrue(not da.isMonotonic(True,1e-12));
        self.assertRaises(InterpKernelException,da.checkMonotonic,True,1e-12);
        da.setIJ(2,0,5.99);
        self.assertTrue(da.isMonotonic(True,1e-12));
        self.assertTrue(not da.isMonotonic(True,1e-1));
        pass

    def testCheckCoherencyDeeper1(self):
        m=MEDCouplingDataForTest.build3DSourceMesh_1();
        m.checkCoherency();
        m.checkCoherency1();
        m.getNodalConnectivity().setIJ(8,0,-1);
        m.checkCoherency();
        self.assertRaises(InterpKernelException,m.checkCoherency1);
        m.getNodalConnectivity().setIJ(8,0,-6);
        m.checkCoherency();
        self.assertRaises(InterpKernelException,m.checkCoherency1);
        m.getNodalConnectivity().setIJ(8,0,9);#9>=NbOfNodes
        m.checkCoherency();
        self.assertRaises(InterpKernelException,m.checkCoherency1);
        m.getNodalConnectivity().setIJ(8,0,8);#OK
        m.checkCoherency();
        m.checkCoherency1();
        elts=[1,5]
        m.convertToPolyTypes(elts);
        m.checkCoherency();
        m.checkCoherency1();
        m.getNodalConnectivity().setIJ(2,0,9);#9>=NbOfNodes
        m.checkCoherency();
        self.assertRaises(InterpKernelException,m.checkCoherency1);
        m.getNodalConnectivity().setIJ(2,0,-3);
        m.checkCoherency();
        self.assertRaises(InterpKernelException,m.checkCoherency1);
        m.getNodalConnectivity().setIJ(2,0,-1);
        m.checkCoherency();
        self.assertRaises(InterpKernelException,m.checkCoherency1);#Throw because cell#0 is not a polyhedron
        m.getNodalConnectivity().setIJ(2,0,4);
        m.checkCoherency();
        m.checkCoherency1();
        m.getNodalConnectivity().setIJ(7,0,-1);
        m.checkCoherency();
        m.checkCoherency1();#OK because we are in polyhedron connec
        m.getNodalConnectivity().setIJ(36,0,14);
        m.checkCoherency();
        self.assertRaises(InterpKernelException,m.checkCoherency1);#Throw beacause now cell 5 is a TETRA4 (14) so mimatch of number index and static type.
        pass

    def testUnPolyze2(self):
        m=MEDCouplingUMesh.New("jjj",3);
        coo=DataArrayDouble.New();
        coo.alloc(4,3);
        coo.rearrange(1);
        coo.iota(0);
        coo.rearrange(3);
        m.setCoords(coo);
        m.allocateCells(2);
        m.insertNextCell(NORM_TETRA4,4,[0,1,2,3]);
        m.insertNextCell(NORM_TETRA4,4,[0,1,2,3]);
        m.finishInsertingCells();
        m2=MEDCouplingUMesh.MergeUMeshesOnSameCoords(4*[m]);
        m2.convertToPolyTypes([2]);
        m2.unPolyze();
        self.assertEqual(NORM_TETRA4,m2.getTypeOfCell(2));
        self.assertEqual(40,m2.getMeshLength());
        temp2=m2.getNodeIdsOfCell(2);
        self.assertEqual(temp2,[0,1,2,3]);
        m2.checkCoherency1();
        m3=m2.deepCpy();
        m2.unPolyze();
        self.assertTrue(m3.isEqual(m2,1e-12));
        pass

    def testDACpyFrom1(self):
        d=DataArrayDouble.New();
        d.alloc(12,1);
        d.iota(14.);
        d.rearrange(3);
        d.setName("Toto");
        d.setInfoOnComponent(0,"X [m]");
        d.setInfoOnComponent(1,"Y [m]");
        d.setInfoOnComponent(2,"Z [m]");
        #
        d1=DataArrayDouble.New();
        self.assertTrue(not d.isEqual(d1,1e-12));
        d1.cpyFrom(d);
        self.assertTrue(d.isEqual(d1,1e-12));
        d1.cpyFrom(d);
        self.assertTrue(d.isEqual(d1,1e-12));
        d1.rearrange(2);
        self.assertTrue(not d.isEqual(d1,1e-12));
        d1.cpyFrom(d);
        self.assertTrue(d.isEqual(d1,1e-12));
        #
        d2=d.convertToIntArr();
        d4=DataArrayInt.New();
        self.assertTrue(not d2.isEqual(d4));
        d4.cpyFrom(d2);
        self.assertTrue(d2.isEqual(d4));
        d4.cpyFrom(d2);
        self.assertTrue(d2.isEqual(d4));
        d4.rearrange(2);
        self.assertTrue(not d2.isEqual(d4));
        d4.cpyFrom(d2);
        self.assertTrue(d2.isEqual(d4));
        pass

    def testDAITransformWithIndArr1(self):
        tab1=[17,18,22,19]
        tab2=[0,1,1,3,3,0,1,3,2,2,3,0]
        expected=[17,18,18,19,19,17,18,19,22,22,19,17]
        d=DataArrayInt.New();
        d.setValues(tab1,4,1);
        d1=DataArrayInt.New();
        d1.setValues(tab2,12,1);
        d2=d1[:]
        #
        d1.transformWithIndArr(d);
        self.assertEqual(12,d1.getNumberOfTuples());
        self.assertEqual(1,d1.getNumberOfComponents());
        for i in xrange(12):
            self.assertEqual(expected[i],d1.getIJ(i,0));
            pass
        #
        d1=d2
        d1.transformWithIndArr(tab1)
        self.assertEqual(12,d1.getNumberOfTuples());
        self.assertEqual(1,d1.getNumberOfComponents());
        for i in xrange(12):
            self.assertEqual(expected[i],d1.getIJ(i,0));
            pass
        pass

    def testDAIBuildPermArrPerLevel1(self):
        arr=[2,0,1,1,0,1,2,0,1,1,0,0]
        expected1=[10,0,5,6,1,7,11,2,8,9,3,4]
        da=DataArrayInt.New();
        da.setValues(arr,12,1);
        da2=da.buildPermArrPerLevel();
        self.assertEqual(12,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        for i in xrange(12):
            self.assertEqual(expected1[i],da2.getIJ(i,0));
            pass
        pass

    def testDAIOperations1(self):
        arr1=[-1,-2,4,7,3,2,6,6,4,3,0,1]
        da=DataArrayInt.New();
        da.setValues(arr1,4,3);
        da1=DataArrayInt.New();
        da1.alloc(12,1);
        da1.iota(2);
        self.assertRaises(InterpKernelException,DataArrayInt.Add,da,da1);#not same number of tuples/Components
        da1.rearrange(3);
        da2=DataArrayInt.Add(da,da1);
        self.assertEqual(4,da2.getNumberOfTuples());
        self.assertEqual(3,da2.getNumberOfComponents());
        expected1=[1,1,8,12,9,9,14,15,14,14,12,14]
        for i in xrange(12):
            self.assertEqual(expected1[i],da2.getIJ(0,i));
            pass
        da1.substractEqual(da);
        expected2=[3,5,0,-2,3,5,2,3,6,8,12,12]
        for i in xrange(12):
            self.assertEqual(expected2[i],da1.getIJ(0,i));
            pass
        da1.rearrange(1); da1.iota(2); da1.rearrange(3);
        da1.addEqual(da);
        for i in xrange(12):
            self.assertEqual(expected1[i],da1.getIJ(0,i));
            pass
        da1.rearrange(1); da1.iota(2); da1.rearrange(3);
        da2=DataArrayInt.Multiply(da,da1);
        self.assertEqual(4,da2.getNumberOfTuples());
        self.assertEqual(3,da2.getNumberOfComponents());
        expected3=[-2,-6,16,35,18,14,48,54,40,33,0,13]
        for i in xrange(12):
            self.assertEqual(expected3[i],da2.getIJ(0,i));
            pass
        da.divideEqual(da1);
        self.assertEqual(4,da.getNumberOfTuples());
        self.assertEqual(3,da.getNumberOfComponents());
        expected4=[0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        for i in xrange(12):
            self.assertEqual(expected4[i],da.getIJ(0,i));
            pass
        da.setValues(arr1,4,3);
        da1.multiplyEqual(da);
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        for i in xrange(12):
            self.assertEqual(expected3[i],da1.getIJ(0,i));
            pass
        da1.rearrange(1); da1.iota(2); da1.rearrange(3);
        da2=DataArrayInt.Divide(da,da1);
        self.assertEqual(4,da2.getNumberOfTuples());
        self.assertEqual(3,da2.getNumberOfComponents());
        for i in xrange(12):
            self.assertEqual(expected4[i],da2.getIJ(0,i));
            pass
        da1.applyInv(321);
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        expected5=[160,107,80,64,53,45,40,35,32,29,26,24]
        for i in xrange(12):
            self.assertEqual(expected5[i],da1.getIJ(0,i));
            pass
        da1.applyDivideBy(2);
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        expected6=[80,53,40,32,26,22,20,17,16,14,13,12]
        for i in xrange(12):
            self.assertEqual(expected6[i],da1.getIJ(0,i));
            pass
        expected7=[3,4,5,4,5,1,6,3,2,0,6,5]
        da1.applyModulus(7);
        for i in xrange(12):
            self.assertEqual(expected7[i],da1.getIJ(0,i));
            pass
        da1.applyLin(1,1);
        expected8=[3,3,3,3,3,1,3,3,0,0,3,3]
        da1.applyRModulus(3);
        for i in xrange(12):
            self.assertEqual(expected8[i],da1.getIJ(0,i));
            pass
        pass

    def testEmulateMEDMEMBDC1(self):
        m,m1=MEDCouplingDataForTest.buildPointe_1();
        m2,da1,da2,da3,da4,da5,da0=m.emulateMEDMEMBDC(m1)
        expected0=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,36,37,32,33,34,35,38,39,40,41,42,43,44,45,46]
        expected1=[1,32,29,23,41,36]
        self.assertEqual(47,da0.getNumberOfTuples());
        self.assertEqual(1,da0.getNumberOfComponents());
        for i in xrange(47):
            self.assertEqual(expected0[i],da0.getIJ(0,i));
            pass
        self.assertEqual(6,da5.getNumberOfTuples());
        self.assertEqual(1,da5.getNumberOfComponents());
        for i in xrange(6):
            self.assertEqual(expected1[i],da5.getIJ(0,i));
            pass
        expected2=[0,1,2,3,4,0,5,6,7,4,8,9,1,7,10,11,12,13,14,5,15,16,17,8,18,19,20,10,21,22,23,2,13,24,25,21,16,26,27,12,19,28,29,15,22,30,31,18,36,26,28,30,24,37,32,33,34,35,38,36,39,40,41,42,37,38,43,44,45,46]
        self.assertEqual(70,da1.getNumberOfTuples());
        self.assertEqual(1,da1.getNumberOfComponents());
        for i in xrange(70):
            self.assertEqual(expected2[i],da1.getIJ(0,i));
            pass
        expected3=[0,4,8,12,16,20,24,28,32,36,40,44,48,53,58,64,70]
        self.assertEqual(17,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        for i in xrange(17):
            self.assertEqual(expected3[i],da2.getIJ(0,i));
            pass
        expected4=[0,2,4,6,7,9,11,12,14,16,17,19,20,22,24,25,27,29,30,32,34,35,37,39,40,42,43,45,46,48,49,51,52,53,54,55,56,58,60,62,63,64,65,66,67,68,69,70]
        #expected4=[0,2,4,6,7,9,11,12,14,16,17,19,20,22,24,25,27,29,30,32,34,35,37,39,40,42,43,45,46,48,49,51,52,54,56,57,58,59,60,62,63,64,65,66,67,68,69,70];
        self.assertEqual(48,da4.getNumberOfTuples());
        self.assertEqual(1,da4.getNumberOfComponents());
        for i in xrange(48):
            self.assertEqual(expected4[i],da4.getIJ(0,i));
            pass
        expected5=[0,1,0,3,0,7,0,1,2,1,4,1,2,3,2,5,2,3,6,3,4,9,4,8,4,5,10,5,9,5,6,11,6,10,6,7,8,7,11,7,8,12,8,9,12,9,10,12,10,11,12,11,13,13,13,13,12,14,13,15,14,15,14,14,14,14,15,15,15,15]
        self.assertEqual(70,da3.getNumberOfTuples());
        self.assertEqual(1,da3.getNumberOfComponents());
        for i in xrange(70):
            self.assertEqual(expected5[i],da3.getIJ(0,i));
            pass
        pass

    def testGetLevArrPerCellTypes1(self):
        m,m1=MEDCouplingDataForTest.buildPointe_1();
        m1,d0,d1,d2,d3=m.buildDescendingConnectivity();
        order=[NORM_TRI3,NORM_QUAD4];
        da0,da1=m1.getLevArrPerCellTypes(order);
        expected0=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1]
        expected1=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,36,37,32,33,34,35,38,39,40,41,42,43,44,45,46]
        self.assertEqual(47,da0.getNumberOfTuples());
        self.assertEqual(1,da0.getNumberOfComponents());
        for i in xrange(47):
            self.assertEqual(expected0[i],da0.getIJ(0,i));
            pass
        self.assertEqual(2,da1.getNumberOfTuples());
        self.assertEqual(1,da1.getNumberOfComponents());
        self.assertEqual(36,da1.getIJ(0,0));#36 TRI3
        self.assertEqual(11,da1.getIJ(1,0));#11 QUAD4
        #
        da2=da0.buildPermArrPerLevel();
        #
        self.assertEqual(47,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        for i in xrange(47):
            self.assertEqual(expected1[i],da2.getIJ(0,i));
            pass
        pass

    def testSortCellsInMEDFileFrmt1(self):
        m,m1=MEDCouplingDataForTest.buildPointe_1();
        m2=m.deepCpy()
        da=DataArrayInt.New()
        da.setValues([0,1,2,14,3,12,4,5,15,6,7,8,9,10,11,13],16,1)
        daa=da.invertArrayN2O2O2N(16)
        m.renumberCells(daa,False)
        da2=m.sortCellsInMEDFileFrmt()
        self.assertEqual(da2.getValues(),[0,1,2,14,3,12,4,5,15,6,7,8,9,10,11,13])
        self.assertTrue(m.isEqual(m2,1e-12))
        self.assertTrue(da.isEqual(da2))
        pass

    def testBuildPartAndReduceNodes1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        arr=[1,0]
        m2,da=m.buildPartAndReduceNodes(arr);
        self.assertEqual(5,m2.getNumberOfNodes());
        self.assertEqual(2,m2.getNumberOfCells());
        f=m2.getMeasureField(True);
        self.assertAlmostEqual(0.125,f.getArray().getIJ(0,0),12);
        self.assertAlmostEqual(0.25,f.getArray().getIJ(1,0),12);
        #
        arr2=DataArrayInt.New()
        arr2.setValues(arr,2,1)
        m2,da=m.buildPartAndReduceNodes(arr2);
        self.assertEqual(5,m2.getNumberOfNodes());
        self.assertEqual(2,m2.getNumberOfCells());
        f=m2.getMeasureField(True);
        self.assertAlmostEqual(0.125,f.getArray().getIJ(0,0),12);
        self.assertAlmostEqual(0.25,f.getArray().getIJ(1,0),12);
        pass

    def testDAITransformWithIndArrR1(self):
        tab1=[2,4,5,3,6,7]
        tab2=[-1,-1,0,1,2,3,4,5,-1,-1,-1,-1]
        expected=[0,3,1,2,4,5]
        d=DataArrayInt.New();
        d.setValues(tab1,6,1);
        d1=DataArrayInt.New();
        d1.setValues(tab2,12,1);
        d2=d1[:]
        #
        d3=d.transformWithIndArrR(d1);
        self.assertEqual(6,d3.getNumberOfTuples());
        self.assertEqual(1,d3.getNumberOfComponents());
        for i in xrange(6):
            self.assertEqual(expected[i],d3.getIJ(i,0));
            pass
        #
        d1=d2
        d3=d.transformWithIndArrR(tab2)
        self.assertEqual(6,d3.getNumberOfTuples());
        self.assertEqual(1,d3.getNumberOfComponents());
        for i in xrange(6):
            self.assertEqual(expected[i],d3.getIJ(i,0));
            pass
        pass

    def testDAISplitByValueRange1(self):
        val1=[6,5,0,3,2,7,8,1,4]
        val2=[0,4,9]
        d=DataArrayInt.New();
        d.setValues(val1,9,1);
        e,f,g=d.splitByValueRange(val2);
        self.assertEqual(9,e.getNumberOfTuples());
        self.assertEqual(1,e.getNumberOfComponents());
        self.assertEqual(9,f.getNumberOfTuples());
        self.assertEqual(1,f.getNumberOfComponents());
        self.assertEqual(2,g.getNumberOfTuples());
        self.assertEqual(1,g.getNumberOfComponents());
        #
        expected1=[1,1,0,0,0,1,1,0,1]
        expected2=[2,1,0,3,2,3,4,1,0]
        for i in xrange(9):
            self.assertEqual(expected1[i],e.getIJ(i,0));
            self.assertEqual(expected2[i],f.getIJ(i,0));
            pass
        self.assertEqual(0,g.getIJ(0,0));
        self.assertEqual(1,g.getIJ(1,0));
        #
        d.setIJ(6,0,9);
        self.assertRaises(InterpKernelException,d.splitByValueRange,val2);
        pass

    def testUMeshSplitProfilePerType1(self):
        val0=[2,0,1,3,4]
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        m.renumberCells(val0,False);
        #
        val1=[0,2,3]
        d=DataArrayInt.New();
        d.setValues(val1,3,1);
        d.setName("sup")
        code,idsInPflPerType,pfls=m.splitProfilePerType(d);
        self.assertEqual(2,len(code));
        self.assertEqual(2,len(idsInPflPerType));
        expected1=[[3,1,0], [4,2,1]]
        self.assertEqual(expected1,code)
        self.assertEqual(2,len(idsInPflPerType));
        self.assertEqual(1,idsInPflPerType[0].getNumberOfTuples());
        self.assertEqual(0,idsInPflPerType[0].getIJ(0,0));
        self.assertEqual(2,idsInPflPerType[1].getNumberOfTuples());
        self.assertEqual(1,idsInPflPerType[1].getIJ(0,0));
        self.assertEqual(2,idsInPflPerType[1].getIJ(1,0));
        #
        self.assertEqual(2,len(pfls));
        self.assertEqual("sup",pfls[0].getName())
        self.assertEqual(1,pfls[0].getNumberOfTuples());
        self.assertEqual(0,pfls[0].getIJ(0,0));
        self.assertEqual("sup",pfls[1].getName())
        self.assertEqual(2,pfls[1].getNumberOfTuples());
        self.assertEqual(0,pfls[1].getIJ(0,0));
        self.assertEqual(1,pfls[1].getIJ(1,0));
        #
        val2=[0,2,3,4]
        d=DataArrayInt.New();
        d.setValues(val2,4,1);
        code,idsInPflPerType,pfls=m.splitProfilePerType(d);
        self.assertEqual(2,len(code));
        self.assertEqual(2,len(idsInPflPerType));
        expected2=[[3,1,0], [4,3,-1]]
        self.assertEqual(expected2,code);
        self.assertEqual(2,len(idsInPflPerType));
        self.assertEqual(1,idsInPflPerType[0].getNumberOfTuples());
        self.assertEqual(0,idsInPflPerType[0].getIJ(0,0));
        self.assertEqual(3,idsInPflPerType[1].getNumberOfTuples());
        self.assertEqual(1,idsInPflPerType[1].getIJ(0,0));
        self.assertEqual(2,idsInPflPerType[1].getIJ(1,0));
        self.assertEqual(3,idsInPflPerType[1].getIJ(2,0));
        #
        self.assertEqual(1,len(pfls));
        self.assertEqual(1,pfls[0].getNumberOfTuples());
        self.assertEqual(0,pfls[0].getIJ(0,0));
        #
        val3=[1,0,2]
        d=DataArrayInt.New();
        d.setValues(val3,3,1);
        code,idsInPflPerType,pfls=m.splitProfilePerType(d);
        self.assertEqual(2,len(code));
        self.assertEqual(2,len(idsInPflPerType));
        expected3=[[3,2,0], [4,1,1]]
        self.assertEqual(expected3,code);
        self.assertEqual(2,len(idsInPflPerType));
        self.assertEqual(2,idsInPflPerType[0].getNumberOfTuples());
        self.assertEqual(0,idsInPflPerType[0].getIJ(0,0));
        self.assertEqual(1,idsInPflPerType[0].getIJ(1,0));
        self.assertEqual(1,idsInPflPerType[1].getNumberOfTuples());
        self.assertEqual(2,idsInPflPerType[1].getIJ(0,0));
        #
        self.assertEqual(2,len(pfls));
        self.assertEqual(2,pfls[0].getNumberOfTuples());
        self.assertEqual(1,pfls[0].getIJ(0,0));
        self.assertEqual(0,pfls[0].getIJ(1,0));
        self.assertEqual(0,pfls[1].getIJ(0,0));
        #
        val4=[3,4]
        d=DataArrayInt.New();
        d.setValues(val4,2,1);
        code,idsInPflPerType,pfls=m.splitProfilePerType(d);
        self.assertEqual(1,len(code));
        self.assertEqual(1,len(idsInPflPerType));
        expected4=[[4,2,0]]
        self.assertEqual(expected4,code);
        self.assertEqual(1,len(idsInPflPerType));
        self.assertEqual(2,idsInPflPerType[0].getNumberOfTuples());
        self.assertEqual(0,idsInPflPerType[0].getIJ(0,0));
        self.assertEqual(1,idsInPflPerType[0].getIJ(1,0));
        #
        self.assertEqual(1,len(pfls));
        self.assertEqual(2,pfls[0].getNumberOfTuples());
        self.assertEqual(1,pfls[0].getIJ(0,0));
        self.assertEqual(2,pfls[0].getIJ(1,0));
        pass

    def testDAIBuildExplicitArrByRanges1(self):
        d=DataArrayInt.New();
        vals1=[0,2,3]
        d.setValues(vals1,3,1);
        e=DataArrayInt.New();
        vals2=[0,3,6,10,14,20]
        e.setValues(vals2,6,1);
        #
        f=d.buildExplicitArrByRanges(e);
        self.assertEqual(11,f.getNumberOfTuples());
        self.assertEqual(1,f.getNumberOfComponents());
        expected1=[0,1,2,6,7,8,9,10,11,12,13]
        for i in xrange(11):
            self.assertEqual(expected1[i],f.getIJ(i,0));
            pass
        pass

    def testDAIComputeOffsets2(self):
        d=DataArrayInt.New();
        vals1=[3,5,1,2,0,8]
        expected1=[0,3,8,9,11,11,19]
        d.setValues(vals1,6,1);
        d.computeOffsets2();
        self.assertEqual(7,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        for i in xrange(7):
            self.assertEqual(expected1[i],d.getIJ(0,i));
            pass
        pass

    def testMergeField3(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        m.getCoords().setInfoOnComponent(0,"x [m]");
        m.getCoords().setInfoOnComponent(1,"z [km]");
        m.setName("m");
        m.setDescription("desc");
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setName("f1");
        f1.setMesh(m);
        arr=DataArrayDouble.New();
        arr.alloc(5,2);
        arr.setInfoOnComponent(0,"X [m]");
        arr.setInfoOnComponent(1,"YY [mm]");
        arr.fillWithValue(2.);
        f1.setArray(arr);
        #
        f2=MEDCouplingFieldDouble.MergeFields([f1]);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        #
        pass
    
    def testGetDistributionOfTypes1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        tab1=[2,0,1,3,4]
        self.assertRaises(InterpKernelException,m.getDistributionOfTypes);
        m.renumberCells(tab1,False);
        code=m.getDistributionOfTypes();
        self.assertEqual(2,len(code));
        self.assertEqual(3,code[0][0]);
        self.assertEqual(2,code[0][1]);
        self.assertEqual(-1,code[0][2]);
        self.assertEqual(4,code[1][0]);
        self.assertEqual(3,code[1][1]);
        self.assertEqual(-1,code[1][2]);
        pass

    def testNorm2_1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f.setMesh(m);
        #
        d=DataArrayDouble.New();
        tab=[1.2,1.3,2.2,2.3,3.2,3.3,4.2,4.3,5.2,5.3]
        d.setValues(tab,5,2);
        f.setArray(d);
        f.checkCoherency();
        #
        self.assertAlmostEqual(11.209371079592289,f.norm2(),14);
        #
        pass

    def testNormMax1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f.setMesh(m);
        #
        d=DataArrayDouble.New();
        tab=[2.3,-1.2,6.3,-7.8,2.9,7.7,2.1,0.,3.6,-7.6]
        d.setValues(tab,5,2);
        f.setArray(d);
        f.checkCoherency();
        #
        self.assertAlmostEqual(7.8,f.normMax(),14);
        #
        pass

    def testFindAndCorrectBadOriented3DExtrudedCells1(self):
        coords=[0.0011180339887498999, -0.0011755705045849499, 0.0, -0.0012331070204200001, -0.0011755705045849499, 0.0, -0.00067557050458494599, -0.00145964954842536, 0.0, -0.00050000000000000001, -0.00086602540378443902, 0.0, 0.00140211303259031, -0.00061803398874989504, 0.0, 0.00086602540378443902, -0.00050000000000000001, 0.0, 0.001, 0.0, 0.0, 0.00034561537182258202, 0.000269164072574575, 0.0, 0.0, 0.001, 0.0, -0.00050000000000000001, 0.00086602540378443902, 0.0, -0.000269164072574575, 0.00034561537182258202, 0.0, -0.001, 0.0, 0.0, -0.00086602540378443902, -0.00050000000000000001, 0.0, -0.00034561537182258202, -0.000269164072574575, 0.0, 0.0, -0.001, 0.0, 0.00050000000000000001, -0.00086602540378443902, 0.0, 0.000269164072574575, -0.00034561537182258202, 0.0, 0.0015, -6.01853107621011e-36, 0.0, 0.00056049747291484397, -0.00145964954842536, 0.0, 0.0011180339887498999, -0.0011755705045849499, 0.00050000000000000001, -0.0012331070204200001, -0.0011755705045849499, 0.00050000000000000001, -0.00067557050458494599, -0.00145964954842536, 0.00050000000000000001, -0.00050000000000000001, -0.00086602540378443902, 0.00050000000000000001, 0.00140211303259031, -0.00061803398874989504, 0.00050000000000000001, 0.00086602540378443902, -0.00050000000000000001, 0.00050000000000000001, 0.001, 0.0, 0.00050000000000000001, 0.00034561537182258202, 0.000269164072574575, 0.00050000000000000001, 0.0, 0.001, 0.00050000000000000001, -0.00050000000000000001, 0.00086602540378443902, 0.00050000000000000001, -0.000269164072574575, 0.00034561537182258202, 0.00050000000000000001, -0.001, 0.0, 0.00050000000000000001, -0.00086602540378443902, -0.00050000000000000001, 0.00050000000000000001, -0.00034561537182258202, -0.000269164072574575, 0.00050000000000000001, 0.0, -0.001, 0.00050000000000000001, 0.00050000000000000001, -0.00086602540378443902, 0.00050000000000000001, 0.000269164072574575, -0.00034561537182258202, 0.00050000000000000001, 0.0015, -6.01853107621011e-36, 0.00050000000000000001, 0.00056049747291484397, -0.00145964954842536, 0.00050000000000000001];
        conn=[2, 1, 3, 21, 20, 22, 4, 0, 5, 23, 19, 24, 8, 9, 10, 27, 28, 29, 11, 12, 13, 30, 31, 32, 0, 18, 15, 5, 19, 37, 34, 24, 6, 17, 4, 5, 25, 36, 23, 24, 3, 14, 16, 13, 22, 33, 35, 32, 13, 16, 7, 10, 32, 35, 26, 29]
        connExp=[16, 2, 1, 3, 21, 20, 22, 16, 4, 0, 5, 23, 19, 24, 16, 8, 10, 9, 27, 29, 28, 16, 11, 13, 12, 30, 32, 31, 18, 0, 18, 15, 5, 19, 37, 34, 24,18, 6, 17, 4, 5, 25, 36, 23, 24, 18, 3, 13, 16, 14, 22, 32, 35, 33, 18, 13, 10, 7, 16, 32, 29, 26, 35]
        invalidCells=[2,3,6,7]
        m=MEDCouplingUMesh.New("Example",3);
        coo=DataArrayDouble.New();
        coo.setValues(coords,38,3);
        m.setCoords(coo);
        m.allocateCells(8);
        m.insertNextCell(NORM_PENTA6,6,conn[0:6])
        m.insertNextCell(NORM_PENTA6,6,conn[6:12])
        m.insertNextCell(NORM_PENTA6,6,conn[12:18])
        m.insertNextCell(NORM_PENTA6,6,conn[18:24])
        m.insertNextCell(NORM_HEXA8,8,conn[24:32])
        m.insertNextCell(NORM_HEXA8,8,conn[32:40])
        m.insertNextCell(NORM_HEXA8,8,conn[40:48])
        m.insertNextCell(NORM_HEXA8,8,conn[48:56])
        m.finishInsertingCells();
        #
        v=m.findAndCorrectBadOriented3DExtrudedCells();
        self.assertEqual(4,len(v));
        self.assertEqual(v.getValues(),invalidCells);
        self.assertEqual(connExp,m.getNodalConnectivity().getValues());
        self.assertTrue(m.findAndCorrectBadOriented3DExtrudedCells().empty())
        #
        pass

    def testConvertExtrudedPolyhedra1(self):
        conn=[1,2,3,4, 5,6,7,8,9,10,11,12, 13,14,15,16, 17,18,19,20,21,22, 23,24,25,26,27,28, 29,30,31,32,33,34,35,36,37,38, 39,40,41,42,43,44,45,46, 47,48,49,50,51,52,53,54,55,56,57,58, 59,60,61,62,63,64,65,66,67,68,69,70,71,72]
        m=MEDCouplingUMesh.New("Example",3);
        coo=DataArrayDouble.New();
        coo.alloc(73,3);
        coo.rearrange(1); coo.iota(0); coo.rearrange(3);
        m.setCoords(coo);
        m.allocateCells(9);
        m.insertNextCell(NORM_TETRA4,4,conn[0:4])
        m.insertNextCell(NORM_HEXA8,8,conn[4:12])
        m.insertNextCell(NORM_TETRA4,4,conn[12:16])
        m.insertNextCell(NORM_POLYHED,6,conn[16:22])
        m.insertNextCell(NORM_PENTA6,6,conn[22:28])
        m.insertNextCell(NORM_POLYHED,10,conn[28:38])
        m.insertNextCell(NORM_HEXA8,8,conn[38:46])
        m.insertNextCell(NORM_HEXGP12,12,conn[46:58])
        m.insertNextCell(NORM_POLYHED,14,conn[58:72])
        m.finishInsertingCells();
        #
        m.convertExtrudedPolyhedra();
        da=m.getNodalConnectivity();
        dai=m.getNodalConnectivityIndex();
        self.assertEqual(10,dai.getNbOfElems());
        self.assertEqual(159,da.getNbOfElems());
        #
        expected1=[14,1,2,3,4,18,5,6,7,8,9,10,11,12,14,13,14,15,16,31,17,18,19,-1,20,22,21,-1,17,20,21,18,-1,18,21,22,19,-1,19,22,20,17,16,23,24,25,26,27,28,31,29,30,31,32,33,-1,34,38,37,36,35,-1,29,34,35,30,-1,30,35,36,31,-1,31,36,37,32,-1,32,37,38,33,-1,33,38,34,29,18,39,40,41,42,43,44,45,46,22,47,48,49,50,51,52,53,54,55,56,57,58,31,59,60,61,62,63,64,65,-1,66,72,71,70,69,68,67,-1,59,66,67,60,-1,60,67,68,61,-1,61,68,69,62,-1,62,69,70,63,-1,63,70,71,64,-1,64,71,72,65,-1,65,72,66,59];
        expected2=[0,5,14,19,42,49,86,95,108,159]
        self.assertEqual(expected1,da.getValues());
        self.assertEqual(expected2,dai.getValues());
        m.checkCoherency2()
        pass

    def testNonRegressionCopyTinyStrings(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f1=m.getMeasureField(True)
        f1.getArray().setInfoOnComponent(0,"P [N/m^2]")
        bary=m.getBarycenterAndOwner()
        f2=f1.buildNewTimeReprFromThis(NO_TIME,False)
        f2.setArray(bary)
        self.assertRaises(InterpKernelException,f1.copyTinyAttrFrom,f2)
        pass

    def testDaDSetPartOfValuesAdv1(self):
        tab1=[3.,4.,5., 13.,14.,15., 23.,24.,25., 33.,34.,35., 43.,44.,45., 53.,54.,55.]
        tab2=[6.,7.,8., 16.,17.,18., 26.,27.,28.]
        tab3=[4,1, 2,2, 3,0]
        a=DataArrayDouble.New();
        a.setValues(tab1,6,3);
        b=DataArrayDouble.New();
        b.setValues(tab2,3,3);
        c=DataArrayInt.New();
        c.setValues(tab3,3,2);
        #
        a.setPartOfValuesAdv(b,c);
        expected1=[3.,4.,5., 13.,14.,15., 26.,27.,28., 6.,7.,8., 16.,17.,18., 53.,54.,55.]
        self.assertEqual(expected1,a.getValues());
        pass

    def testUMeshBuildSetInstanceFromThis1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=m.buildSetInstanceFromThis(3);
        self.assertTrue(m.isEqual(m2,1e-12));
        #
        m=MEDCouplingUMesh.New("toto",2);
        m2=m.buildSetInstanceFromThis(3);
        self.assertEqual(0,m2.getNumberOfNodes());
        self.assertEqual(0,m2.getNumberOfCells());
        pass

    def testUMeshMergeMeshesCVW1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=MEDCouplingUMesh.New("toto",2);
        m3=MEDCouplingUMesh.MergeUMeshes([m,m2]);
        m3.setName(m.getName());
        self.assertTrue(m.isEqual(m3,1e-12));
        pass
    
    def testChangeUnderlyingMeshWithCMesh1(self):
        mesh=MEDCouplingCMesh.New();
        coordsX=DataArrayDouble.New();
        arrX=[ -1., 1., 2., 4. ]
        coordsX.setValues(arrX,4,1);
        coordsY=DataArrayDouble.New();
        arrY=[ -2., 2., 4., 8. ]
        coordsY.setValues(arrY,4,1);
        coordsZ=DataArrayDouble.New();
        arrZ=[ -3., 3., 6., 12. ]
        coordsZ.setValues(arrZ,4,1);
        mesh.setCoords(coordsX,coordsY,coordsZ);
        f=mesh.getMeasureField(True)
        mesh2=mesh.deepCpy()
        for myId in [0,1,2,10,11,12,20,21,22]:
            f=mesh.getMeasureField(True)
            f.changeUnderlyingMesh(mesh2,myId,1e-12);
            pass
        mesh2.setName("uuuu")
        for myId in [1,2,10,11,12,20,21,22]:
            f=mesh.getMeasureField(True)
            f.changeUnderlyingMesh(mesh2,myId,1e-12);
            pass
        pass

    def testDADFindCommonTuples1(self):
        da=DataArrayDouble.New();
        # nbOftuples=1
        array1=[2.3,1.2,1.3,2.3,2.301,0.8]
        da.setValues(array1,6,1)
        c,cI=da.findCommonTuples(1e-2);
        expected1=[0,3,4]
        expected2=[0,3]
        self.assertEqual(3,c.getNbOfElems());
        self.assertEqual(2,cI.getNbOfElems());
        self.assertEqual(expected1,c.getValues())
        self.assertEqual(expected2,cI.getValues())
        c,cI=da.findCommonTuples(2e-1)
        expected3=[0,3,4,1,2]
        expected4=[0,3,5]
        self.assertEqual(5,c.getNbOfElems());
        self.assertEqual(3,cI.getNbOfElems());
        self.assertEqual(expected3,c.getValues())
        self.assertEqual(expected4,cI.getValues())
        # nbOftuples=2
        array2=[2.3,2.3,1.2,1.2,1.3,1.3,2.3,2.3,2.301,2.301,0.8,0.8]
        da.setValues(array2,6,2)
        c,cI=da.findCommonTuples(1e-2);
        self.assertEqual(3,c.getNbOfElems());
        self.assertEqual(2,cI.getNbOfElems());
        self.assertEqual(expected1,c.getValues())
        self.assertEqual(expected2,cI.getValues())
        c,cI=da.findCommonTuples(2e-1)
        self.assertEqual(5,c.getNbOfElems());
        self.assertEqual(3,cI.getNbOfElems());
        self.assertEqual(expected3,c.getValues())
        self.assertEqual(expected4,cI.getValues())
        # nbOftuples=3
        array3=[2.3,2.3,2.3,1.2,1.2,1.2,1.3,1.3,1.3,2.3,2.3,2.3,2.301,2.301,2.301,0.8,0.8,0.8]
        da.setValues(array3,6,3)
        c,cI=da.findCommonTuples(1e-2);
        self.assertEqual(3,c.getNbOfElems());
        self.assertEqual(2,cI.getNbOfElems());
        self.assertEqual(expected1,c.getValues())
        self.assertEqual(expected2,cI.getValues())
        c,cI=da.findCommonTuples(2e-1)
        self.assertEqual(5,c.getNbOfElems());
        self.assertEqual(3,cI.getNbOfElems());
        self.assertEqual(expected3,c.getValues())
        self.assertEqual(expected4,cI.getValues())
        # nbOftuples=1, no common groups
        array11=[2.3,1.2,1.3,2.4,2.5,0.8]
        da.setValues(array11,6,1)
        c,cI=da.findCommonTuples(1e-2);
        self.assertEqual(0,c.getNbOfElems());
        self.assertEqual(1,cI.getNbOfElems());
        self.assertEqual([0],cI.getValues())
        
        array12=[0.]*(6*5)
        da.setValues(array12,6,5) #bad NumberOfComponents
        self.assertRaises(InterpKernelException, da.findCommonTuples, 1e-2);
        pass

    def testDABack1(self):
        da=DataArrayDouble.New();
        array1=[2.3,1.2,1.3,2.3,2.301,0.8]
        da.setValues(array1,6,1);
        self.assertAlmostEqual(0.8,da.back(),14);
        da.rearrange(2);
        self.assertRaises(InterpKernelException,da.back);
        da.alloc(0,1);
        self.assertRaises(InterpKernelException,da.back);
        #
        da=DataArrayInt.New();
        array2=[4,7,8,2]
        da.setValues(array2,4,1);
        self.assertEqual(2,da.back());
        da.rearrange(2);
        self.assertRaises(InterpKernelException,da.back);
        da.alloc(0,1);
        self.assertRaises(InterpKernelException,da.back);
        pass

    def testDADGetDifferentValues1(self):
        da=DataArrayDouble.New();
        array1=[2.3,1.2,1.3,2.3,2.301,0.8]
        da.setValues(array1,6,1)
        #
        expected1=[2.301,1.2,1.3,0.8]
        dv=da.getDifferentValues(1e-2);
        self.assertEqual(4,dv.getNbOfElems());
        for i in xrange(4):
            self.assertAlmostEqual(expected1[i],dv.getIJ(i,0),14);
            pass
        #
        dv=da.getDifferentValues(2e-1);
        expected2=[2.301,1.3,0.8]
        self.assertEqual(3,dv.getNbOfElems());
        for i in xrange(3):
            self.assertAlmostEqual(expected2[i],dv.getIJ(i,0),14);
            pass
        pass

    def testDAIBuildOld2NewArrayFromSurjectiveFormat2(self):
        arr=[0,3, 5,7,9]
        arrI=[0,2,5]
        a=DataArrayInt.New();
        a.setValues(arr,5,1);
        b=DataArrayInt.New();
        b.setValues(arrI,3,1);
        ret,newNbTuple=DataArrayInt.BuildOld2NewArrayFromSurjectiveFormat2(10,a,b);
        expected=[0,1,2,0,3,4,5,4,6,4]
        self.assertEqual(10,ret.getNbOfElems());
        self.assertEqual(7,newNbTuple);
        self.assertEqual(1,ret.getNumberOfComponents());
        self.assertEqual(expected,ret.getValues());
        self.assertRaises(InterpKernelException,DataArrayInt.BuildOld2NewArrayFromSurjectiveFormat2,9,a,b);
        pass

    def testDADIReverse1(self):
        arr=[0,3,5,7,9,2]
        a=DataArrayInt.New();
        a.setValues(arr,6,1);
        self.assertEqual(2,a.back());
        a.reverse();
        for i in xrange(6):
            self.assertEqual(arr[5-i],a.getIJ(i,0));
            pass
        a.setValues(arr[:-1],5,1);
        a.reverse();
        for i in xrange(5):
            self.assertEqual(arr[4-i],a.getIJ(i,0));
            pass
        #
        arr2=[0.,3.,5.,7.,9.,2.]
        b=DataArrayDouble.New();
        b.setValues(arr2,6,1);
        b.reverse();
        for i in xrange(6):
            self.assertAlmostEqual(arr2[5-i],b.getIJ(i,0),14);
            pass
        b.setValues(arr2[:5],5,1);
        self.assertAlmostEqual(9.,b.back(),14)
        b.reverse();
        for i in xrange(5):
            self.assertAlmostEqual(arr2[4-i],b.getIJ(i,0),14);
            pass
        pass

    def testGetNodeIdsInUse1(self):
        m0=MEDCouplingDataForTest.build2DTargetMesh_1();
        CellIds=[1,2]
        m1=m0.buildPartOfMySelf(CellIds,True);
        arr,newNbOfNodes=m1.getNodeIdsInUse();
        expected=[-1,0,1,-1,2,3,-1,-1,-1]
        self.assertEqual(4,newNbOfNodes);
        self.assertEqual(9,arr.getNbOfElems());
        self.assertEqual(expected,arr.getValues());
        arr2=arr.invertArrayO2N2N2O(newNbOfNodes);
        self.assertEqual(4,arr2.getNbOfElems());
        expected2=[1,2,4,5]
        self.assertEqual(expected2,arr2.getValues());
        pass

    def testBuildDescendingConnec2(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        #
        mesh2,desc,descIndx,revDesc,revDescIndx=mesh.buildDescendingConnectivity2();
        mesh2.checkCoherency();
        self.assertEqual(1,mesh2.getMeshDimension());
        self.assertEqual(13,mesh2.getNumberOfCells());
        self.assertEqual(14,revDescIndx.getNbOfElems()); self.assertEqual(14,revDescIndx.getNumberOfTuples());
        self.assertEqual(6,descIndx.getNbOfElems()); self.assertEqual(6,descIndx.getNumberOfTuples());
        self.assertEqual(18,desc.getNbOfElems()); self.assertEqual(18,desc.getNumberOfTuples());
        self.assertEqual(18,revDesc.getNbOfElems()); self.assertEqual(18,revDesc.getNumberOfTuples());
        expected1=[1,2,3,4,-3,5,6, 7,8,-5,9,10,-2,11, 12,13,-7,-10]
        self.assertEqual(expected1,desc.getValues());
        expected2=[0,4,7,10,14,18]
        self.assertEqual(expected2,descIndx.getValues());
        expected3=[0,1,3,5,6,8,9,11,12,13,15,16,17,18]
        self.assertEqual(expected3,revDescIndx.getValues());
        expected4=[0, 0,3, 0,1, 0, 1,2, 1, 2,4, 2, 3, 3,4, 3, 4, 4]
        self.assertEqual(expected4,revDesc.getValues());
        conn=mesh2.getNodalConnectivity();
        connIndex=mesh2.getNodalConnectivityIndex();
        expected5=[0,3,6,9,12,15,18,21,24,27,30,33,36,39]
        self.assertEqual(expected5,connIndex.getValues());
        expected6=[1, 0, 3, 1, 3, 4, 1, 4, 1, 1, 1, 0, 1, 4, 2, 1, 2, 1, 1, 4, 5, 1, 5, 2, 1, 6, 7, 1, 7, 4, 1, 3, 6, 1, 7, 8, 1, 8, 5]
        self.assertEqual(expected6,conn.getValues());
        pass

    def testIntersect2DMeshesTmp1(self):
        m1c=MEDCouplingCMesh.New();
        coordsX=DataArrayDouble.New();
        arrX=[ -1., 1., 2., 4. ]
        coordsX.setValues(arrX,4,1);
        m1c.setCoordsAt(0,coordsX);
        coordsY=DataArrayDouble.New();
        arrY=[ -2., 2., 4., 8. ]
        coordsY.setValues(arrY,4,1);
        m1c.setCoordsAt(1,coordsY);
        m1=m1c.buildUnstructured()
        m1bis=m1.buildPartOfMySelf([3,4,5],False)
        m2=m1.deepCpy()
        m2=m2.buildPartOfMySelf([0,1,2],False)
        m2.translate([0.5,0.5])
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1bis,m2,1e-10)
        expected1=[0,0,1,1,1,2,2,2]
        expected2=[0,-1,0,1,-1,1,2,-1]
        self.assertEqual(8,d1.getNumberOfTuples());
        self.assertEqual(8,d2.getNumberOfTuples());
        self.assertEqual(8,m3.getNumberOfCells());
        self.assertEqual(22,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[5,17,1,16,12,5,16,0,4,5,17,12,5,18,1,17,13,5,19,2,18,13,5,17,5,6,19,13,5,20,2,19,14,5,21,3,20,14,5,19,6,7,21,14]
        expected4=[0,5,12,17,22,28,33,38,44]
        expected5=[-1.0,2.0,1.0,2.0,2.0,2.0,4.0,2.0,-1.0,4.0,1.0,4.0,2.0,4.0,4.0,4.0,-0.5,-1.5,1.5,-1.5,2.5,-1.5,4.5,-1.5,-0.5,2.5,1.5,2.5,2.5,2.5,4.5,2.5,-0.5,2.0,1.0,2.5,1.5,2.0,2.0,2.5,2.5,2.0,4.0,2.5]
        self.assertEqual(44,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(9,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in xrange(44):
            self.assertAlmostEqual(expected5[i],m3.getCoords().getIJ(0,i),12);
            pass
        pass

    def testFindNodesOnLine1(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        pt=[-0.3,-0.3]
        pt2=[0.,0.,0.]
        pt3=[-0.3,0.,0.]
        vec=[0.,1.]
        vec2=[1.,0.,0.]
        vec3=[0.,1.,1.]
        expected1=[0,3,6]
        res=mesh.findNodesOnLine(pt,vec,1e-12);
        self.assertEqual(3,len(res));
        self.assertEqual(expected1,res.getValues());
        #
        mesh.changeSpaceDimension(3);
        mesh.rotate(pt2,vec2,pi/4.);
        res=mesh.findNodesOnLine(pt3,vec3,1e-12);
        self.assertEqual(3,len(res));
        self.assertEqual(expected1,res.getValues());
        pass

    def testIntersect2DMeshesTmp2(self):
        m1c=MEDCouplingCMesh.New();
        coordsX1=DataArrayDouble.New();
        arrX1=[ 0., 1., 1.5, 2. ]
        coordsX1.setValues(arrX1,4,1);
        m1c.setCoordsAt(0,coordsX1);
        coordsY1=DataArrayDouble.New();
        arrY1=[ 0., 1.5, 3.]
        coordsY1.setValues(arrY1,3,1);
        m1c.setCoordsAt(1,coordsY1);
        m1=m1c.buildUnstructured();
        m2c=MEDCouplingCMesh.New();
        coordsX2=DataArrayDouble.New();
        arrX2=[ 0., 1., 2. ]
        coordsX2.setValues(arrX2,3,1);
        m2c.setCoordsAt(0,coordsX2);
        coordsY2=DataArrayDouble.New();
        arrY2=[ 0., 1., 3.]
        coordsY2.setValues(arrY2,3,1);
        m2c.setCoordsAt(1,coordsY2);
        m2=m2c.buildUnstructured();
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1,m2,1e-10)
        #
        expected1=[0,0,1,1,2,2,3,4,5]
        expected2=[0,2,1,3,1,3,2,3,3]
        self.assertEqual(9,d1.getNumberOfTuples());
        self.assertEqual(9,d2.getNumberOfTuples());
        self.assertEqual(9,m3.getNumberOfCells());
        self.assertEqual(22,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[5,16,13,12,15,5,15,4,5,16,5,21,2,13,16,5,16,5,6,21,5,17,14,2,21,5,21,6,7,17,5,4,18,19,5,5,5,19,10,6,5,6,10,20,7]
        expected4=[0,5,10,15,20,25,30,35,40,45]
        expected5=[0.0,0.0,1.0,0.0,1.5,0.0,2.0,0.0,0.0,1.5,1.0,1.5,1.5,1.5,2.0,1.5,0.0,3.0,1.0,3.0,1.5,3.0,2.0,3.0,0.0,0.0,1.0,0.0,2.0,0.0,0.0,1.0,1.0,1.0,2.0,1.0,0.0,3.0,1.0,3.0,2.0,3.0,1.5,1.0]
        self.assertEqual(45,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(10,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in xrange(44):
            self.assertAlmostEqual(expected5[i],m3.getCoords().getIJ(0,i),12);
            pass
        pass
    
    def testBuildPartOfMySelfSafe1(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1()
        self.assertRaises(InterpKernelException,mesh.buildPartOfMySelf,[0,-1,4,2],True)
        self.assertRaises(InterpKernelException,mesh.buildPartOfMySelf,[0,4,5,4],True)
        pass

    def testIntersect2DMeshesTmp3(self):
        m1Coords=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1,0.,-1.5,0.5,0.,1.25,0.,0.70710678118654757,0.70710678118654757,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.70710678118654757,0.70710678118654757,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.70710678118654757,-0.70710678118654757,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.70710678118654757,-0.70710678118654757,1.0606601717798214,-1.0606601717798214];
        m1Conn=[0,3,1,13,11,9, 3,4,2,1,14,12,10,11, 5,3,0,15,13,17, 6,4,3,5,16,14,15,18, 5,0,7,17,21,19, 6,5,7,8,18,19,22,20, 0,1,7,9,23,21, 1,2,8,7,10,24,22,23];
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(2);
        m1.allocateCells(8);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[0:6]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[6:14]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[14:20]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[20:28]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[28:34]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[34:42]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[42:48]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[48:56]);
        m1.finishInsertingCells();
        myCoords1=DataArrayDouble.New();
        myCoords1.setValues(m1Coords,25,2);
        m1.setCoords(myCoords1);
        #
        m2Coords=[0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1,-1.1,-1.,0.,-1.,1.1,-1,1.7,-1.]
        m2Conn=[0,3,2,1, 1,2,5,4, 7,6,3,0, 8,9,6,7, 7,0,12,11, 8,7,11,10, 0,1,13,12, 1,4,14,13]
        m2=MEDCouplingUMesh.New();
        m2.setMeshDimension(2);
        m2.allocateCells(8);
        for i in xrange(8):
            m2.insertNextCell(NORM_QUAD4,4,m2Conn[4*i:4*(i+1)])
            pass
        m2.finishInsertingCells();
        myCoords2=DataArrayDouble.New();
        myCoords2.setValues(m2Coords,15,2);
        m2.setCoords(myCoords2);
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1,m2,1e-10)
        m3.unPolyze()
        #
        expected1=[0,1,1,1,2,3,3,3,4,5,5,5,6,7,7,7]
        expected2=[0,0,1,-1,2,2,3,-1,4,4,5,-1,6,6,7,-1]
        self.assertEqual(16,d1.getNumberOfTuples());
        self.assertEqual(16,d2.getNumberOfTuples());
        self.assertEqual(16,m3.getNumberOfCells());
        self.assertEqual(104,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[6,28,1,25,44,45,46,8,26,1,28,27,47,48,49,50,8,40,2,26,27,51,52,53,54,8,28,4,40,27,55,56,57,58,6,28,25,5,59,60,61,8,28,5,32,31,62,63,64,65,8,32,6,41,31,66,67,68,69,8,41,4,28,31,70,71,72,73,6,25,37,5,74,75,76,8,32,5,37,36,77,78,79,80,8,42,6,32,36,81,82,83,84,8,37,8,42,36,85,86,87,88,6,1,37,25,89,90,91,8,37,1,26,38,92,93,94,95,8,26,2,43,38,96,97,98,99,8,43,8,37,38,100,101,102,103]
        expected4=[0,7,16,25,34,41,50,59,68,75,84,93,102,109,118,127,136]
        expected5=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1.,-1.1,-1.,0.,-1.,1.1,-1.,1.7,-1.,1.118033988749895,1.,-1.118033988749895,1.,-1.118033988749895,-1.,1.118033988749895,-1.,0.7071067811865477,0.7071067811865476,0.5,0.,0.,0.5,1.05,0.,0.7071067811865475,0.7071067811865477,0.55,1.,1.1,0.5,1.4012585384440737,0.535233134659635,1.3,0.,1.1,0.5,1.1090169943749475,1.,0.,1.25,0.6123724356957946,1.369306393762915,1.1090169943749475,1.,0.55,1.,0.,0.5,-0.5,0.,-0.7071067811865477,0.7071067811865476,-0.7071067811865475,0.7071067811865477,-1.05,0.,-1.1,0.5,-0.55,1.,-1.3,0.,-1.4012585384440737,0.5352331346596344,-1.1090169943749475,1.,-1.1,0.5,-0.6123724356957941,1.3693063937629155,0.,1.25,-0.55,1.,-1.1090169943749475,1.,0.,-0.5,-0.7071067811865475,-0.7071067811865477,-0.5,0.,-1.05,0.,-0.7071067811865478,-0.7071067811865475,-0.55,-1.,-1.1,-0.5,-1.4012585384440734,-0.5352331346596354,-1.3,0.,-1.1,-0.5,-1.1090169943749475,-1.,0.,-1.25,-0.6123724356957945,-1.369306393762915,-1.1090169943749475,-1.,-0.55,-1.,0.7071067811865475,-0.7071067811865477,0.,-0.5,0.5,0.,0.7071067811865477,-0.7071067811865475,1.05,0.,1.1,-0.5,0.55,-1.,1.3,0.,1.4012585384440737,-0.535233134659635,1.1090169943749475,-1.,1.1,-0.5,0.6123724356957946,-1.369306393762915,0.,-1.25,0.55,-1.,1.1090169943749475,-1.0]
        self.assertEqual(136,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(17,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in xrange(208):
            self.assertAlmostEqual(expected5[i],m3.getCoords().getIJ(0,i),12);
            pass
        pass

    def testUMeshTessellate2D1(self):
        m1Coords=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1,0.,-1.5,0.5,0.,1.25,0.,0.70710678118654757,0.70710678118654757,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.70710678118654757,0.70710678118654757,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.70710678118654757,-0.70710678118654757,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.70710678118654757,-0.70710678118654757,1.0606601717798214,-1.0606601717798214];
        m1Conn=[0,3,1,13,11,9, 3,4,2,1,14,12,10,11, 5,3,0,15,13,17, 6,4,3,5,16,14,15,18, 5,0,7,17,21,19, 6,5,7,8,18,19,22,20, 0,1,7,9,23,21, 1,2,8,7,10,24,22,23];
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(2);
        m1.allocateCells(8);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[0:6]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[6:14]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[14:20]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[20:28]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[28:34]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[34:42]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[42:48]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[48:56]);
        m1.finishInsertingCells();
        myCoords1=DataArrayDouble.New();
        myCoords1.setValues(m1Coords,25,2);
        m1.setCoords(myCoords1);
        #
        m11=m1.deepCpy();
        m11.tessellate2D(1.);
        self.assertTrue(m11.getCoords().isEqual(m11.getCoords(),1e-12));
        expected1=[5,0,3,11,1,5,3,4,12,2,1,11,5,5,15,3,0,5,6,16,4,3,15,5,5,5,0,7,19,5,6,5,19,7,8,20,5,0,1,23,7,5,1,2,24,8,7,23]
        expected2=[0,5,12,17,24,29,36,41,48]
        self.assertEqual(48,m11.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(9,m11.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected1,m11.getNodalConnectivity().getValues());
        self.assertEqual(expected2,m11.getNodalConnectivityIndex().getValues());
        #
        m12=m1.deepCpy();
        m12.tessellate2D(0.5);
        self.assertEqual(41,m12.getNumberOfNodes());
        expected3=[5,0,3,25,26,1,5,3,4,27,28,2,1,26,25,5,5,29,30,3,0,5,6,31,32,4,3,30,29,5,5,5,0,7,33,34,5,6,5,34,33,7,8,35,36,5,0,1,37,38,7,5,1,2,39,40,8,7,38,37]
        expected4=[0,6,15,21,30,36,45,51,60]
        expected5=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,0.479425538604203,0.8775825618903728,0.8414709848078964,0.54030230586814,0.7191383079063044,1.3163738428355591,1.2622064772118446,0.8104534588022099,-0.877582561890373,0.4794255386042027,-0.5403023058681399,0.8414709848078964,-1.3163738428355596,0.7191383079063038,-0.8104534588022098,1.2622064772118446,-0.4794255386042031,-0.8775825618903728,-0.8414709848078965,-0.5403023058681399,-0.7191383079063045,-1.3163738428355591,-1.2622064772118449,-0.8104534588022098,0.8775825618903729,-0.47942553860420295,0.54030230586814,-0.8414709848078964,1.3163738428355594,-0.7191383079063043,0.8104534588022099,-1.2622064772118446]
        for i in xrange(82):
            self.assertAlmostEqual(expected5[i],m12.getCoords().getIJ(0,i),12);
            pass
        self.assertEqual(60,m12.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(9,m12.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m12.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m12.getNodalConnectivityIndex().getValues());
        pass

    def testUMeshTessellate2DCurve1(self):
        # A quarter of circle:
        mcoords = [0.4,0.0,   0.0,-0.4,   0.283,-0.283]
        mconnec = [0,1,2]

        m1 = MEDCouplingUMesh.New()
        m1.setMeshDimension(1)
        m1.allocateCells(1)
        m1.insertNextCell(NORM_SEG3, mconnec)

        myCoords = DataArrayDouble.New(mcoords, 3, 2)
        m1.setCoords(myCoords)
        
        m2 = m1.deepCpy()
        m2.tessellate2DCurve(0.1)
        # If the following raises, the test will fail automatically:
        m2.checkCoherency1(0.0) # eps param not used

    def testIntersect2DMeshesTmp4(self):
        m1Coords=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1,0.,-1.5,0.5,0.,1.25,0.,0.70710678118654757,0.70710678118654757,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.70710678118654757,0.70710678118654757,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.70710678118654757,-0.70710678118654757,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.70710678118654757,-0.70710678118654757,1.0606601717798214,-1.0606601717798214];
        m1Conn=[0,3,1,13,11,9, 3,4,2,1,14,12,10,11, 5,3,0,15,13,17, 6,4,3,5,16,14,15,18, 5,0,7,17,21,19, 6,5,7,8,18,19,22,20, 0,1,7,9,23,21, 1,2,8,7,10,24,22,23];
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(2);
        m1.allocateCells(8);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[0:6]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[6:14]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[14:20]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[20:28]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[28:34]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[34:42]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[42:48]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[48:56]);
        m1.finishInsertingCells();
        myCoords1=DataArrayDouble.New();
        myCoords1.setValues(m1Coords,25,2);
        m1.setCoords(myCoords1);
        #
        m2Coords=[0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1,-1.1,-1.,0.,-1.,1.1,-1,1.7,-1.]
        m2Conn=[0,3,2,1, 1,2,5,4, 7,6,3,0, 8,9,6,7, 7,0,12,11, 8,7,11,10, 0,1,13,12, 1,4,14,13]
        m2=MEDCouplingUMesh.New();
        m2.setMeshDimension(2);
        m2.allocateCells(8);
        for i in xrange(8):
            m2.insertNextCell(NORM_QUAD4,4,m2Conn[4*i:4*(i+1)])
            pass
        m2.finishInsertingCells();
        myCoords2=DataArrayDouble.New();
        myCoords2.setValues(m2Coords,15,2);
        m2.setCoords(myCoords2);
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m2,m1,1e-10)
        m3.unPolyze()
        #
        expected1=[0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7]
        expected2=[0,1,1,-1,2,3,3,-1,4,5,5,-1,6,7,7,-1]
        self.assertEqual(16,d1.getNumberOfTuples());
        self.assertEqual(16,d2.getNumberOfTuples());
        self.assertEqual(16,m3.getNumberOfCells());
        self.assertEqual(104,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[6,16,15,18,44,45,46,8,18,2,1,16,47,48,49,50,8,17,1,2,40,51,52,53,54,8,40,5,4,17,55,56,57,58,6,18,15,20,59,60,61,8,20,7,6,18,62,63,64,65,8,41,6,7,21,66,67,68,69,8,21,8,9,41,70,71,72,73,6,20,15,22,74,75,76,8,22,11,7,20,77,78,79,80,8,21,7,11,42,81,82,83,84,8,42,10,8,21,85,86,87,88,6,22,15,16,89,90,91,8,16,1,13,22,92,93,94,95,8,43,13,1,17,96,97,98,99,8,17,4,14,43,100,101,102,103]
        expected4=[0,7,16,25,34,41,50,59,68,75,84,93,102,109,118,127,136]
        expected5=[0.,0.,1.1, 0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1.,-1.1,-1.,0.,-1.,1.1,-1.,1.7,-1.,0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,1.1180339887498951,1.,-1.1180339887498951,1.,-1.1180339887498951,-1.,1.1180339887498951,-1.,0.5,0.,0.,0.5,0.7071067811865477,0.7071067811865476,0.55,1.,1.1,0.5,1.05,0.,0.7071067811865477,0.7071067811865475,1.3,0.,1.1,0.5,1.1090169943749475,1.,1.4012585384440737,0.535233134659635,1.4090169943749475,1.,1.7,0.5,1.6,0.,1.4012585384440737,0.535233134659635,0.,0.5,-0.5,0.,-0.7071067811865477,0.7071067811865476,-1.05,0.,-1.1,0.5,-0.55,1.,-0.7071067811865478,0.7071067811865475,-1.1090169943749475,1.,-1.1,0.5,-1.3,0.,-1.4012585384440737,0.5352331346596344,-1.6,0.,-1.7,0.5,-1.4090169943749475,1.,-1.4012585384440737,0.5352331346596344,-0.5,0.,0.,-0.5,-0.7071067811865475,-0.7071067811865477,-0.55,-1.,-1.1,-0.5,-1.05,0.,-0.7071067811865475,-0.7071067811865477,-1.3,0.,-1.1,-0.5,-1.1090169943749475,-1.,-1.4012585384440734,-0.5352331346596354,-1.4090169943749475,-1.,-1.7,-0.5,-1.6,0.,-1.4012585384440732,-0.5352331346596354,0.,-0.5,0.5,0.,0.7071067811865475,-0.7071067811865477,1.05,0.,1.1,-0.5,0.55,-1.,0.7071067811865475,-0.7071067811865477,1.1090169943749475,-1.,1.1,-0.5,1.3,0.,1.4012585384440737,-0.535233134659635,1.6,0.,1.7,-0.5,1.4090169943749475,-1.,1.4012585384440737,-0.535233134659635]
        self.assertEqual(136,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(17,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in xrange(208):
            self.assertAlmostEqual(expected5[i],m3.getCoords().getIJ(0,i),12);
            pass
        pass

    def testGetCellIdsCrossingPlane1(self):
        mesh3D,mesh2D=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        vec=[-0.07,1.,0.07]
        origin=[1.524,1.4552,1.74768]
        ids1=mesh3D.getCellIdsCrossingPlane(origin,vec,1e-10)
        self.assertEqual([1,3,4,7,9,10,13,15,16],ids1.getValues())
        vec2=[0.,0.,1.]
        ids2=mesh3D.getCellIdsCrossingPlane(origin,vec2,1e-10)
        self.assertEqual([6,7,8,9,10,11],ids2.getValues())
        pass

    def testBuildSlice3D1(self):
        mesh3D,mesh2D=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        vec1=[-0.07,1.,0.07]
        origin1=[1.524,1.4552,1.74768]
        slice1,ids=mesh3D.buildSlice3D(origin1,vec1,1e-10);
        expected1=[1,3,4,7,9,10,13,15,16]
        expected2=[5,42,41,40,43,44,5,42,46,45,41,5,44,43,40,47,48,5,49,42,44,50,5,49,51,46,42,5,50,44,48,52,5,53,49,50,54,5,53,55,51,49,5,54,50,52,56]
        expected3=[0,6,11,17,22,27,32,37,42,47]
        expected4=[1.,1.,0.,1.,1.25,0.,1.,1.5,0.,2.,1.,0.,1.,2.,0.,0.,2.,0.,3.,1.,0.,3.,2.,0.,0.,1.,0.,2.,2.,0.,1.,1.,1.,1.,1.25,1.,1.,1.5,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,2.,2.,1.,1.,1.,2.,1.,1.25,2.,1.,1.5,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,2.,2.,2.,1.,1.,3.,1.,1.25,3.,1.,1.5,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,2.,2.,3.,1.,1.5408576,0.,2.,1.6108576000000001,0.,2.,1.5408576,1.,1.,1.5,0.5836800000000008,1.,1.4708576,1.,3.,1.6808576,0.,3.,1.6108576000000001,1.,0.,1.4708576,0.,0.,1.4008576,1.,2.,1.4708576,2.,1.,1.4008576000000001,2.,3.,1.5408575999999998,2.,0.,1.3308575999999999,2.,2.,1.4008576,3.,1.,1.3308576,3.,3.,1.4708576,3.,0.,1.2608576,3.]
        self.assertEqual(2,slice1.getMeshDimension());
        self.assertEqual(3,slice1.getSpaceDimension());
        self.assertEqual(57,slice1.getNumberOfNodes());
        self.assertEqual(9,slice1.getNumberOfCells());
        self.assertEqual(9,ids.getNumberOfTuples());
        self.assertEqual(47,slice1.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(10,slice1.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected1,ids.getValues());
        self.assertEqual(expected2,slice1.getNodalConnectivity().getValues());
        self.assertEqual(expected3,slice1.getNodalConnectivityIndex().getValues());
        for i in xrange(171):
            self.assertAlmostEqual(expected4[i],slice1.getCoords().getIJ(0,i),12);
            pass
        # 2nd slice based on already existing nodes of mesh3D.
        vec2=[0.,3.,1.]
        origin2=[2.5,1.,3.]
        slice1,ids=mesh3D.buildSlice3D(origin2,vec2,1e-10);
        expected5=[5,50,10,4,51,5,50,52,7,10,5,51,4,5,53,5,54,50,51,55,56,5,54,57,52,50,5,56,55,51,53,58,5,38,59,56,54,43,5,54,57,46,43,5,38,59,56,58,48]
        expected6=[0,5,10,15,21,26,32,38,43,49]
        expected7=[1.,1.,0.,1.,1.25,0.,1.,1.5,0.,2.,1.,0.,1.,2.,0.,0.,2.,0.,3.,1.,0.,3.,2.,0.,0.,1.,0.,1.,3.,0.,2.,2.,0.,2.,3.,0.,1.,1.,1.,1.,1.25,1.,1.,1.5,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,1.,3.,1.,2.,2.,1.,2.,3.,1.,0.,0.,2.,1.,1.,2.,1.,1.25,2.,1.,0.,2.,1.,1.5,2.,2.,0.,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,2.,2.,2.,0.,0.,3.,1.,1.,3.,1.,1.25,3.,1.,0.,3.,1.,1.5,3.,2.,0.,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,2.,2.,3.,2.,1.6666666666666667,1.,1.,1.6666666666666667,1.,3.,1.6666666666666667,1.,0.,1.6666666666666667,1.,2.,1.3333333333333335,2.,1.,1.5,1.5,1.,1.3333333333333333,2.,3.,1.3333333333333335,2.,0.,1.3333333333333335,2.,1.,1.25,2.25]
        self.assertEqual(2,slice1.getMeshDimension());
        self.assertEqual(3,slice1.getSpaceDimension());
        self.assertEqual(60,slice1.getNumberOfNodes());
        self.assertEqual(9,slice1.getNumberOfCells());
        self.assertEqual(9,ids.getNumberOfTuples());
        self.assertEqual(49,slice1.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(10,slice1.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected1,ids.getValues());
        self.assertEqual(expected5,slice1.getNodalConnectivity().getValues());
        self.assertEqual(expected6,slice1.getNodalConnectivityIndex().getValues());
        for i in xrange(180):
            self.assertAlmostEqual(expected7[i],slice1.getCoords().getIJ(0,i),12);
            pass
        # 3rd slice based on shared face of mesh3D.
        vec3=[0.,0.,1.]
        origin3=[2.5,1.,2.]
        slice1,ids=mesh3D.buildSlice3D(origin3,vec3,1e-10);
        expected8=[6,7,8,9,10,11,12,13,14,15,16,17]
        expected9=[5,15,26,16,18,5,16,21,28,22,19,17,5,18,20,21,16,5,21,24,25,28,5,26,16,17,19,22,23,5,22,27,29,28,5,15,26,16,18,5,16,21,28,22,19,17,5,18,20,21,16,5,21,24,25,28,5,26,16,17,19,22,23,5,22,27,29,28]
        expected10=[0,5,12,17,22,29,34,39,46,51,56,63,68]
        expected11=[0.,0.,1.,1.,1.,1.,1.,1.25,1.,1.,0.,1.,1.,1.5,1.,2.,0.,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,1.,3.,1.,2.,2.,1.,2.,3.,1.,0.,0.,2.,1.,1.,2.,1.,1.25,2.,1.,0.,2.,1.,1.5,2.,2.,0.,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,1.,3.,2.,2.,2.,2.,2.,3.,2.,0.,0.,3.,1.,1.,3.,1.,1.25,3.,1.,0.,3.,1.,1.5,3.,2.,0.,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,1.,3.,3.,2.,2.,3.,2.,3.,3.]
        self.assertEqual(2,slice1.getMeshDimension());
        self.assertEqual(3,slice1.getSpaceDimension());
        self.assertEqual(45,slice1.getNumberOfNodes());
        self.assertEqual(12,slice1.getNumberOfCells());
        self.assertEqual(12,ids.getNumberOfTuples());
        self.assertEqual(68,slice1.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(13,slice1.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected8,ids.getValues());
        self.assertEqual(expected9,slice1.getNodalConnectivity().getValues());
        self.assertEqual(expected10,slice1.getNodalConnectivityIndex().getValues());
        for i in xrange(135):
            self.assertAlmostEqual(expected11[i],slice1.getCoords().getIJ(0,i),12);
            pass
        pass

    def testBuildSlice3DSurf1(self):
        mesh3D,mesh2D=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        mesh2D=mesh3D.buildDescendingConnectivity()[0];
        vec1=[-0.07,1.,0.07]
        origin1=[1.524,1.4552,1.74768]
        slice1,ids=mesh2D.buildSlice3DSurf(origin1,vec1,1e-10);
        expected1=[6,8,10,11,13,18,19,21,23,25,26,38,41,43,47,49,52,53,64,67,69,73,75,78,79]
        expected2=[1,40,41,1,42,41,1,40,43,1,44,43,1,42,44,1,45,41,1,42,46,1,46,45,1,47,40,1,47,48,1,44,48,1,49,42,1,44,50,1,49,50,1,49,51,1,51,46,1,48,52,1,50,52,1,53,49,1,50,54,1,53,54,1,53,55,1,55,51,1,52,56,1,54,56]
        expected3=[0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75];
        expected4=[1.,1.,0.,1.,1.25,0.,1.,1.5,0.,2.,1.,0.,1.,2.,0.,0.,2.,0.,3.,1.,0.,3.,2.,0.,0.,1.,0.,2.,2.,0.,1.,1.,1.,1.,1.25,1.,1.,1.5,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,2.,2.,1.,1.,1.,2.,1.,1.25,2.,1.,1.5,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,2.,2.,2.,1.,1.,3.,1.,1.25,3.,1.,1.5,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,2.,2.,3.,1.,1.5408576,0.,2.,1.6108576000000001,0.,2.,1.5408576,1.,1.,1.5,0.5836800000000008,1.,1.4708576,1.,3.,1.6808576,0.,3.,1.6108576000000001,1.,0.,1.4708576,0.,0.,1.4008576,1.,2.,1.4708576,2.,1.,1.4008576000000001,2.,3.,1.5408575999999998,2.,0.,1.3308575999999999,2.,2.,1.4008576,3.,1.,1.3308576,3.,3.,1.4708576,3.,0.,1.2608576,3.]
        self.assertEqual(1,slice1.getMeshDimension());
        self.assertEqual(3,slice1.getSpaceDimension());
        self.assertEqual(57,slice1.getNumberOfNodes());
        self.assertEqual(25,slice1.getNumberOfCells());
        self.assertEqual(25,ids.getNumberOfTuples());
        self.assertEqual(75,slice1.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(26,slice1.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected1,ids.getValues());
        self.assertEqual(expected2,slice1.getNodalConnectivity().getValues());
        self.assertEqual(expected3,slice1.getNodalConnectivityIndex().getValues());
        for i in xrange(171):
            self.assertAlmostEqual(expected4[i],slice1.getCoords().getIJ(0,i),12);
            pass
        #
        vec2=[0.,0.,1.]
        origin2=[2.5,1.,2.]
        slice1,ids=mesh2D.buildSlice3DSurf(origin2,vec2,1e-10);
        expected5=[32,32,32,32,33,34,35,36,37,38,39,40,41,42,43,43,43,43,43,43,44,44,44,44,45,46,47,47,47,47,48,49,50,51,52,53,53,53,53,53,53,54,54,54,54,55,56,57,59,60,61,62,63,64,65,66,67,68,71,72,74,75,76,77,78,81,82,83]
        expected6=[1,15,18,1,18,16,1,16,26,1,26,15,1,26,15,1,16,26,1,18,16,1,15,18,1,16,21,1,21,28,1,22,28,1,19,22,1,17,19,1,16,17,1,16,21,1,21,28,1,28,22,1,22,19,1,19,17,1,17,16,1,16,18,1,18,20,1,20,21,1,21,16,1,20,21,1,18,20,1,28,21,1,21,24,1,24,25,1,25,28,1,25,28,1,24,25,1,21,24,1,23,22,1,26,23,1,26,16,1,16,17,1,17,19,1,19,22,1,22,23,1,23,26,1,22,28,1,28,29,1,29,27,1,27,22,1,27,22,1,29,27,1,28,29,1,26,15,1,16,26,1,18,16,1,15,18,1,16,21,1,21,28,1,22,28,1,19,22,1,17,19,1,16,17,1,20,21,1,18,20,1,25,28,1,24,25,1,21,24,1,23,22,1,26,23,1,27,22,1,29,27,1,28,29]
        expected7=[0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111,114,117,120,123,126,129,132,135,138,141,144,147,150,153,156,159,162,165,168,171,174,177,180,183,186,189,192,195,198,201,204];
        expected8=[0.,0.,1.,1.,1.,1.,1.,1.25, 1.,1.,0.,1.,1.,1.5, 1.,2.,0.,1.,2.,1.,1.,1.,2.,1.,0.,2.,1.,3.,1.,1.,3.,2.,1.,0.,1.,1.,1.,3.,1.,2.,2.,1.,2.,3.,1.,0.,0.,2.,1.,1.,2.,1.,1.25, 2.,1.,0.,2.,1.,1.5, 2.,2.,0.,2.,2.,1.,2.,1.,2.,2.,0.,2.,2.,3.,1.,2.,3.,2.,2.,0.,1.,2.,1.,3.,2.,2.,2.,2.,2.,3.,2.,0.,0.,3.,1.,1.,3.,1.,1.25, 3.,1.,0.,3.,1.,1.5, 3.,2.,0.,3.,2.,1.,3.,1.,2.,3.,0.,2.,3.,3.,1.,3.,3.,2.,3.,0.,1.,3.,1.,3.,3.,2.,2.,3.,2.,3.,3.]
        self.assertEqual(1,slice1.getMeshDimension());
        self.assertEqual(3,slice1.getSpaceDimension());
        self.assertEqual(45,slice1.getNumberOfNodes());
        self.assertEqual(68,slice1.getNumberOfCells());
        self.assertEqual(68,ids.getNumberOfTuples());
        self.assertEqual(204,slice1.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(69,slice1.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected5,ids.getValues());
        self.assertEqual(expected6,slice1.getNodalConnectivity().getValues());
        self.assertEqual(expected7,slice1.getNodalConnectivityIndex().getValues());
        for i in xrange(135):
            self.assertAlmostEqual(expected8[i],slice1.getCoords().getIJ(0,i),12);
            pass
        pass

    def testDataArrayDoubleAdvSetting1(self):
        data1=[1.,11.,2.,12.,3.,13.,4.,14.,5.,15.,6.,16.,7.,17.]
        data2=[8.,38.,9.,39.,0.,30.,11.,41.,12.,42.]
        compsCpp=["comp1","comp2"]
        da=DataArrayDouble.New();
        da.setInfoAndChangeNbOfCompo(compsCpp);
        da.setName("da");
        da.alloc(7,2);
        compsCpp=compsCpp[:-1]
        self.assertRaises(InterpKernelException,da.setInfoAndChangeNbOfCompo,compsCpp);
        da.setValues(data1,7,2)
        #
        p=[(0,3),(3,5),(5,7)]
        tmp=da.selectByTupleRanges(p);
        self.assertTrue(tmp.isEqual(da,1e-14));
        p=[(0,2),(3,4),(5,7)]
        tmp=da.selectByTupleRanges(p);
        expected1=[1.,11.,2.,12.,4.,14.,6.,16.,7.,17.]
        self.assertEqual(5,tmp.getNumberOfTuples());
        self.assertEqual(2,tmp.getNumberOfComponents());
        for i in xrange(10):
            self.assertAlmostEqual(expected1[i],tmp.getIJ(0,i),14);
            pass
        p=[(0,2),(0,2),(5,6)]
        tmp=da.selectByTupleRanges(p);
        expected2=[1.,11.,2.,12.,1.,11.,2.,12.,6.,16.]
        self.assertEqual(5,tmp.getNumberOfTuples());
        self.assertEqual(2,tmp.getNumberOfComponents());
        for i in xrange(10):
            self.assertAlmostEqual(expected2[i],tmp.getIJ(0,i),14);
            pass
        p=[(0,2),(-1,2),(5,6)]
        self.assertRaises(InterpKernelException,da.selectByTupleRanges,p);
        p=[(0,2),(0,2),(5,8)]
        self.assertRaises(InterpKernelException,da.selectByTupleRanges,p);
        #
        da2=DataArrayDouble.New();
        da2.setValues(data2,5,2);
        #
        dac=da.deepCpy();
        dac.setContigPartOfSelectedValues2(1,da2,2,4,1);
        expected3=[1.,11.,0.,30.,11.,41.,4.,14.,5.,15.,6.,16.,7.,17.]
        for i in xrange(14):
            self.assertAlmostEqual(expected3[i],dac.getIJ(0,i),14);
            pass
        #
        dac=da.deepCpy();
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues2,3,da2,0,5,1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues2,0,da2,4,6,1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues2,3,da2,5,0,1);
        dac.setContigPartOfSelectedValues2(3,da2,1,5,1);
        expected4=[1.,11.,2.,12.,3.,13.,9.,39.,0.,30.,11.,41.,12.,42.]
        for i in xrange(14):
            self.assertAlmostEqual(expected4[i],dac.getIJ(0,i),14);
            pass
        #
        ids=DataArrayInt.New();
        ids.alloc(3,1);
        dac=da.deepCpy();
        ids.setIJ(0,0,2); ids.setIJ(1,0,0); ids.setIJ(2,0,4);
        dac.setContigPartOfSelectedValues(2,da2,ids);
        expected5=[1.,11.,2.,12.,0.,30.,8.,38.,12.,42.,6.,16.,7.,17.]
        for i in xrange(14):
            self.assertAlmostEqual(expected5[i],dac.getIJ(0,i),14);
            pass
        #
        dac=da.deepCpy();
        ids.setIJ(0,0,2); ids.setIJ(1,0,5); ids.setIJ(2,0,4);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues,1,da2,ids);
        ids.setIJ(0,0,2); ids.setIJ(1,0,2); ids.setIJ(2,0,-1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues,1,da2,ids);
        ids.setIJ(0,0,2); ids.setIJ(1,0,2); ids.setIJ(2,0,1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues,5,da2,ids);
        #
        ids.setIJ(0,0,2); ids.setIJ(1,0,2); ids.setIJ(2,0,1);
        dac=da.deepCpy();
        dac.setContigPartOfSelectedValues(4,da2,ids);
        expected6=[1.,11.,2.,12.,3.,13.,4.,14.,0.,30.,0.,30.,9.,39.]
        for i in xrange(14):
            self.assertAlmostEqual(expected6[i],dac.getIJ(0,i),14);
            pass
        pass

    def testDataArrayIntAdvSetting1(self):
        data1=[1,11,2,12,3,13,4,14,5,15,6,16,7,17]
        data2=[8,38,9,39,0,30,11,41,12,42]
        compsCpp=["comp1","comp2"]
        da=DataArrayInt.New();
        da.setInfoAndChangeNbOfCompo(compsCpp);
        da.setName("da");
        da.alloc(7,2);
        compsCpp=compsCpp[:-1]
        self.assertRaises(InterpKernelException,da.setInfoAndChangeNbOfCompo,compsCpp);
        da.setValues(data1,7,2)
        #
        p=[(0,3),(3,5),(5,7)]
        tmp=da.selectByTupleRanges(p);
        self.assertTrue(tmp.isEqual(da));
        p=[(0,2),(3,4),(5,7)]
        tmp=da.selectByTupleRanges(p);
        expected1=[1,11,2,12,4,14,6,16,7,17]
        self.assertEqual(5,tmp.getNumberOfTuples());
        self.assertEqual(2,tmp.getNumberOfComponents());
        for i in xrange(10):
            self.assertEqual(expected1[i],tmp.getIJ(0,i));
            pass
        p=[(0,2),(0,2),(5,6)]
        tmp=da.selectByTupleRanges(p);
        expected2=[1,11,2,12,1,11,2,12,6,16]
        self.assertEqual(5,tmp.getNumberOfTuples());
        self.assertEqual(2,tmp.getNumberOfComponents());
        for i in xrange(10):
            self.assertEqual(expected2[i],tmp.getIJ(0,i));
            pass
        p=[(0,2),(-1,2),(5,6)]
        self.assertRaises(InterpKernelException,da.selectByTupleRanges,p);
        p=[(0,2),(0,2),(5,8)]
        self.assertRaises(InterpKernelException,da.selectByTupleRanges,p);
        #
        da2=DataArrayInt.New();
        da2.setValues(data2,5,2);
        #
        dac=da.deepCpy();
        dac.setContigPartOfSelectedValues2(1,da2,2,4,1);
        expected3=[1,11,0,30,11,41,4,14,5,15,6,16,7,17]
        for i in xrange(14):
            self.assertEqual(expected3[i],dac.getIJ(0,i));
            pass
        #
        dac=da.deepCpy();
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues2,3,da2,0,5,1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues2,0,da2,4,6,1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues2,3,da2,5,0,1);
        dac.setContigPartOfSelectedValues2(3,da2,1,5,1);
        expected4=[1,11,2,12,3,13,9,39,0,30,11,41,12,42]
        for i in xrange(14):
            self.assertEqual(expected4[i],dac.getIJ(0,i));
            pass
        #
        ids=DataArrayInt.New();
        ids.alloc(3,1);
        dac=da.deepCpy();
        ids.setIJ(0,0,2); ids.setIJ(1,0,0); ids.setIJ(2,0,4);
        dac.setContigPartOfSelectedValues(2,da2,ids);
        expected5=[1,11,2,12,0,30,8,38,12,42,6,16,7,17]
        for i in xrange(14):
            self.assertEqual(expected5[i],dac.getIJ(0,i));
            pass
        #
        dac=da.deepCpy();
        ids.setIJ(0,0,2); ids.setIJ(1,0,5); ids.setIJ(2,0,4);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues,1,da2,ids);
        ids.setIJ(0,0,2); ids.setIJ(1,0,2); ids.setIJ(2,0,-1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues,1,da2,ids);
        ids.setIJ(0,0,2); ids.setIJ(1,0,2); ids.setIJ(2,0,1);
        self.assertRaises(InterpKernelException,dac.setContigPartOfSelectedValues,5,da2,ids);
        #
        ids.setIJ(0,0,2); ids.setIJ(1,0,2); ids.setIJ(2,0,1);
        dac=da.deepCpy();
        dac.setContigPartOfSelectedValues(4,da2,ids);
        expected6=[1,11,2,12,3,13,4,14,0,30,0,30,9,39]
        for i in xrange(14):
            self.assertEqual(expected6[i],dac.getIJ(0,i));
            pass
        pass

    def testBuildDescendingConnec2Of3DMesh1(self):
        mesh=MEDCouplingDataForTest.build3DSourceMesh_1();
        #
        mesh2,desc,descIndx,revDesc,revDescIndx=mesh.buildDescendingConnectivity2();
        mesh2.checkCoherency();
        self.assertEqual(2,mesh2.getMeshDimension());
        self.assertEqual(30,mesh2.getNumberOfCells());
        self.assertEqual(31,revDescIndx.getNbOfElems()); self.assertEqual(31,revDescIndx.getNumberOfTuples());
        self.assertEqual(13,descIndx.getNbOfElems()); self.assertEqual(13,descIndx.getNumberOfTuples());
        self.assertEqual(48,desc.getNbOfElems()); self.assertEqual(48,desc.getNumberOfTuples());
        self.assertEqual(48,revDesc.getNbOfElems()); self.assertEqual(48,revDesc.getNumberOfTuples());
        expected1=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,-10,15,-5,-13,16,17,-14,18,-4,19,-2,20,21,22,23,24,25,-11,26,-1,-12,-25,-22,27,28,-7,-20,-24,29,-16,-18,30,-8,-28]
        self.assertEqual(expected1,desc.getValues());
        expected2=[0,4,8,12,16,20,24,28,32,36,40,44,48]
        self.assertEqual(expected2,descIndx.getValues());
        expected3=[0,2,4,5,7,9,10,12,14,15,17,19,21,23,25,26,28,29,31,32,34,35,37,38,40,42,43,44,46,47,48]
        self.assertEqual(expected3,revDescIndx.getValues());
        expected4=[0,8,0,6,0,0,5,1,4,1,1,9,1,11,2,2,3,2,7,2,8,3,4,3,5,3,4,10,4,5,11,5,6,10,6,6,9,7,7,10,7,8,8,9,9,11,10,11]
        self.assertEqual(expected4,revDesc.getValues());
        conn=mesh2.getNodalConnectivity();
        connIndex=mesh2.getNodalConnectivityIndex();
        expected5=[0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120]
        self.assertEqual(expected5,connIndex.getValues());
        expected6=[3,8,1,7,3,8,3,1,3,1,3,7,3,7,3,8,3,6,0,8,3,6,2,0,3,0,2,8,3,8,2,6,3,7,4,5,3,7,8,4,3,4,8,5,3,5,8,7,3,6,8,4,3,6,7,8,3,4,7,6,3,8,4,0,3,0,4,6,3,6,3,8,3,7,3,6,3,8,0,1,3,1,0,3,3,3,0,8,3,4,1,5,3,4,8,1,3,1,8,5,3,1,7,5,3,0,2,3,3,3,2,8,3,1,4,0,3,3,2,6]
        self.assertEqual(expected6,conn.getValues());
        pass

    def testAre2DCellsNotCorrectlyOriented1(self):
        m1Coords=[1.,1.,-1.,-1.,-1.,-1.,1.,-1.]
        m1Conn=[0,3,1,2]
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(2);
        m1.allocateCells(1);
        m1.insertNextCell(NORM_QUAD4,4,m1Conn[0:4])
        m1.finishInsertingCells();
        myCoords1=DataArrayDouble.New();
        myCoords1.setValues(m1Coords,4,2);
        m1.setCoords(myCoords1);
        #
        vec1=[0.,0.,1.]
        for i in xrange(18):
            vec2=[3.*cos(pi/9.*i),3.*sin(pi/9.*i)];
            m1Cpy=m1.deepCpy();
            m1Cpy.translate(vec2);
            self.assertRaises(InterpKernelException,m1Cpy.are2DCellsNotCorrectlyOriented,vec1,False);
            m1Cpy.changeSpaceDimension(3);
            res=m1Cpy.are2DCellsNotCorrectlyOriented(vec1,False)
            self.assertEqual([0],res.getValues());
            pass
        pass

    def testDataArrayAbs1(self):
        d1=DataArrayDouble.New();
        val1=[2.,-3.,-5.,6.,-7.,-8.,9.,10.,-11.,-12.,-13.,-15.]
        expected1=[2.,3.,5.,6.,7.,8.,9.,10.,11.,12.,13.,15.]
        d1.setValues(val1,6,2);
        d2=d1.convertToIntArr();
        #
        d1.abs();
        for i in xrange(12):
            self.assertAlmostEqual(expected1[i],d1.getIJ(0,i),14);
            pass
        #
        expected2=[2,3,5,6,7,8,9,10,11,12,13,15]
        d2.abs();
        for i in xrange(12):
            self.assertEqual(expected2[i],d2.getIJ(0,i));
            pass
        #
        pass

    # test on 1D
    def testGetValueOn3(self):
        v=[0.,1.,1.5,2.]
        v2=[0.7,1.25,0.,2.,1.5]
        disp=[5.,50.,500.,6.,60.,600.,7.,70.,700.,8.,80.,800.]
        m=MEDCouplingUMesh.New("myMesh",1)
        nbNodes=len(v)
        nbCells=nbNodes-1
        m.allocateCells(nbCells)
        coords=DataArrayDouble.New() ; coords.setValues(v,nbNodes,1)
        m.setCoords(coords)
        m.insertNextCell(NORM_SEG2,2,[0,1])
        m.insertNextCell(NORM_SEG2,2,[2,1])
        m.insertNextCell(NORM_SEG2,2,[2,3])
        m.finishInsertingCells()
        f=MEDCouplingFieldDouble.New(ON_NODES)
        f.setMesh(m)
        array=DataArrayDouble.New(); array.setValues(disp,m.getNumberOfNodes(),3)
        f.setArray(array)
        arr1=f.getValueOnMulti(v2)
        self.assertEqual(5,arr1.getNumberOfTuples());
        self.assertEqual(3,arr1.getNumberOfComponents());
        expected1=[5.7,57.,570.,6.5,65.,650.,5.,50.,500.,8.,80.,800.,7.,70.,700.]
        for i in xrange(15):
            self.assertAlmostEqual(expected1[i],arr1.getIJ(0,i),14);
            pass
        pass

    def testGetNodeIdsOfCell2(self):
        m1c=MEDCouplingCMesh.New();
        coordsX=DataArrayDouble.New();
        arrX=[ -1., 1., 2., 4., 4.5 ]
        coordsX.setValues(arrX,5,1);
        coordsY=DataArrayDouble.New();
        arrY=[ -2., 2., 4., 8.]
        coordsY.setValues(arrY,4,1);
        coordsZ=DataArrayDouble.New();
        arrZ=[ -2., 2., 4.]
        coordsZ.setValues(arrZ,3,1);
        # test in 1D
        m1c.setCoordsAt(0,coordsX);
        expected1=[[0,1],[1,2],[2,3],[3,4]]
        self.assertEqual(4,m1c.getNumberOfCells())
        for i in xrange(m1c.getNumberOfCells()):
            self.assertEqual(expected1[i],m1c.getNodeIdsOfCell(i))
            pass
        # test in 2D
        m1c.setCoordsAt(1,coordsY);
        self.assertEqual(12,m1c.getNumberOfCells())
        self.assertEqual(20,m1c.getNumberOfNodes())
        expected2=[[0,1,6,5],[1,2,7,6],[2,3,8,7],[3,4,9,8],[5,6,11,10],[6,7,12,11],[7,8,13,12],[8,9,14,13],[10,11,16,15],[11,12,17,16],[12,13,18,17],[13,14,19,18]]
        for i in xrange(m1c.getNumberOfCells()):
            self.assertEqual(expected2[i],m1c.getNodeIdsOfCell(i))
            pass
        # test in 3D
        m1c.setCoordsAt(2,coordsZ);
        self.assertEqual(24,m1c.getNumberOfCells())
        self.assertEqual(60,m1c.getNumberOfNodes())
        expected3=[[0,1,6,5,20,21,26,25],[1,2,7,6,21,22,27,26],[2,3,8,7,22,23,28,27],[3,4,9,8,23,24,29,28],[5,6,11,10,25,26,31,30],[6,7,12,11,26,27,32,31],[7,8,13,12,27,28,33,32],[8,9,14,13,28,29,34,33],[10,11,16,15,30,31,36,35],[11,12,17,16,31,32,37,36],[12,13,18,17,32,33,38,37],[13,14,19,18,33,34,39,38],[20,21,26,25,40,41,46,45],[21,22,27,26,41,42,47,46],[22,23,28,27,42,43,48,47],[23,24,29,28,43,44,49,48],[25,26,31,30,45,46,51,50],[26,27,32,31,46,47,52,51],[27,28,33,32,47,48,53,52],[28,29,34,33,48,49,54,53],[30,31,36,35,50,51,56,55],[31,32,37,36,51,52,57,56],[32,33,38,37,52,53,58,57],[33,34,39,38,53,54,59,58]]
        self.assertEqual(24,m1c.getNumberOfCells())
        for i in xrange(m1c.getNumberOfCells()):
            self.assertEqual(expected3[i],m1c.getNodeIdsOfCell(i))
            pass
        pass
    
    def testSwigDADOp4(self):
        da=DataArrayDouble.New(range(6,30),12,2)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),float(i+6),13)
            pass
        # operator transpose
        da.transpose()
        self.assertEqual(2,da.getNumberOfTuples());
        self.assertEqual(12,da.getNumberOfComponents());
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),float(i+6),13)
            pass
        da.transpose()
        # operator __neg__
        da2=DataArrayDouble.New(12,1)
        da2.iota(0.)
        dabis=-da
        for i in xrange(24):
            self.assertAlmostEqual(dabis.getIJ(0,i),-float(i+6),13)
            pass
        # operator+=
        da+=da2
        expected1=[6.,7.,9.,10.,12.,13.,15.,16.,18.,19.,21.,22.,24.,25.,27.,28.,30.,31.,33.,34.,36.,37.,39.,40.]
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected1[i],13)
            pass
        da=-dabis
        da+=[100.,101.]
        expected2=[106.,108.,108.,110.,110.,112.,112.,114.,114.,116.,116.,118.,118.,120.,120.,122.,122.,124.,124.,126.,126.,128.,128.,130.]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected2[i],13)
            pass
        for pos,elt in enumerate(dabis):
            da[pos]+=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for elt in da:
            li=elt[:]
            self.assertAlmostEqual(li[0],100.,13) ; self.assertAlmostEqual(li[1],101.,13)
            pass
        # operator-=
        da=DataArrayDouble.New(range(6,30),12,2)
        da2=DataArrayDouble.New(range(12),12,1)
        dabis=-da
        da-=da2
        expected1=[6.,7.,7.,8.,8.,9.,9.,10.,10.,11.,11.,12.,12.,13.,13.,14.,14.,15.,15.,16.,16.,17.,17.,18.]
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected1[i],13)
            pass
        da=-dabis
        da-=[100.,101.]
        expected2=[-94.,-94.,-92.,-92.,-90.,-90.,-88.,-88.,-86.,-86.,-84.,-84.,-82.,-82.,-80.,-80.,-78.,-78.,-76.,-76.,-74.,-74.,-72.,-72.]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected2[i],13)
            pass
        for pos,elt in enumerate(dabis):
            da[pos]-=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        expected3=[-88.,-87.,-84.,-83.,-80.,-79.,-76.,-75.,-72.,-71.,-68.,-67.,-64.,-63.,-60.,-59.,-56.,-55.,-52.,-51.,-48.,-47.,-44.,-43.]
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected3[i],13)
            pass
        # operator*=
        da=DataArrayDouble.New(range(6,30),12,2)
        da2=DataArrayDouble.New(range(12),12,1)
        dabis=-da
        da*=da2
        expected1=[0.,0.,8.,9.,20.,22.,36.,39.,56.,60.,80.,85.,108.,114.,140.,147.,176.,184.,216.,225.,260.,270.,308.,319.]
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected1[i],13)
            pass
        da=-dabis
        da*=[100.,101.]
        expected2=[600.,707.,800.,909.,1000.,1111.,1200.,1313.,1400.,1515.,1600.,1717.,1800.,1919.,2000.,2121.,2200.,2323.,2400.,2525.,2600.,2727.,2800.,2929.]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected2[i],13)
            pass
        for pos,elt in enumerate(dabis):
            da[pos]*=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        expected3=[-3600.,-4949.,-6400.,-8181.,-10000.,-12221.,-14400.,-17069.,-19600.,-22725.,-25600.,-29189.,-32400.,-36461.,-40000.,-44541.,-48400.,-53429.,-57600.,-63125.,-67600.,-73629.,-78400.,-84941.0]
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected3[i],13)
            pass
        # operator/=
        da=DataArrayDouble.New(range(6,30),12,2)
        da2=DataArrayDouble.New(range(1,13),12,1)
        dabis=-da
        da/=da2
        expected1=[6.0,7.0,4.0,4.5,3.3333333333333335,3.6666666666666665,3.0,3.25,2.8,3.0,2.6666666666666665,2.8333333333333335,2.5714285714285716,2.7142857142857144,2.5,2.625,2.4444444444444446,2.5555555555555554,2.4,2.5,2.3636363636363638,2.4545454545454546,2.3333333333333335,2.4166666666666665]
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected1[i],13)
            pass
        da=-dabis
        da/=[100.,101.]
        expected2=[0.06,0.06930693069306931,0.08,0.0891089108910891,0.1,0.10891089108910891,0.12,0.12871287128712872,0.14,0.1485148514851485,0.16,0.16831683168316833,0.18,0.18811881188118812,0.2,0.2079207920792079,0.22,0.22772277227722773,0.24,0.24752475247524752,0.26,0.26732673267326734,0.28,0.2871287128712871]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected2[i],13)
            pass
        for pos,elt in enumerate(dabis):
            da[pos]/=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        expected3=[-0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.009900990099009901, -0.01, -0.0099009900990099]
        for i in xrange(24):
            self.assertAlmostEqual(da.getIJ(0,i),expected3[i],13)
            pass
        pass

    def testSwigDAIOp4(self):
        da=DataArrayInt.New(range(6,30),12,2)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in xrange(24):
            self.assertEqual(da.getIJ(0,i),i+6)
            pass
        # operator transpose
        da.transpose()
        self.assertEqual(2,da.getNumberOfTuples());
        self.assertEqual(12,da.getNumberOfComponents());
        for i in xrange(24):
            self.assertEqual(da.getIJ(0,i),i+6)
            pass
        da.transpose()
        # operator __neg__
        da2=DataArrayInt.New(12,1)
        da2.iota(0)
        dabis=-da
        for i in xrange(24):
            self.assertEqual(dabis.getIJ(0,i),-(i+6))
            pass
        # operator+=
        da+=da2
        expected1=[6,7,9,10,12,13,15,16,18,19,21,22,24,25,27,28,30,31,33,34,36,37,39,40]
        for i in xrange(24):
            self.assertEqual(da.getIJ(0,i),expected1[i])
            pass
        da=-dabis
        da+=[100,101]
        expected2=[106,108,108,110,110,112,112,114,114,116,116,118,118,120,120,122,122,124,124,126,126,128,128,130]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in xrange(24):
            self.assertEqual(da.getIJ(0,i),expected2[i])
            pass
        for pos,elt in enumerate(dabis):
            da[pos]+=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for elt in da:
            li=elt[:]
            self.assertEqual(li[0],100) ; self.assertEqual(li[1],101)
            pass
        # operator-=
        da=DataArrayInt.New(range(6,30),12,2)
        da2=DataArrayInt.New(range(12),12,1)
        dabis=-da
        da-=da2
        expected1=[6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18]
        for i in xrange(24):
            self.assertEqual(da.getIJ(0,i),expected1[i])
            pass
        da=-dabis
        da-=[100,101]
        expected2=[-94,-94,-92,-92,-90,-90,-88,-88,-86,-86,-84,-84,-82,-82,-80,-80,-78,-78,-76,-76,-74,-74,-72,-72]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in xrange(24):
            self.assertEqual(da.getIJ(0,i),expected2[i])
            pass
        for pos,elt in enumerate(dabis):
            da[pos]-=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        expected3=[-88,-87,-84,-83,-80,-79,-76,-75,-72,-71,-68,-67,-64,-63,-60,-59,-56,-55,-52,-51,-48,-47,-44,-43]
        for i in xrange(24):
            self.assertEqual(da.getIJ(0,i),expected3[i])
            pass
        # operator*=
        da=DataArrayInt.New(range(6,30),12,2)
        da2=DataArrayInt.New(range(12),12,1)
        dabis=-da
        da*=da2
        expected1=[0,0,8,9,20,22,36,39,56,60,80,85,108,114,140,147,176,184,216,225,260,270,308,319]
        for i in xrange(24):
            self.assertEqual(da.getIJ(0,i),expected1[i])
            pass
        da=-dabis
        da*=[100,101]
        expected2=[600,707,800,909,1000,1111,1200,1313,1400,1515,1600,1717,1800,1919,2000,2121,2200,2323,2400,2525,2600,2727,2800,2929]
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in xrange(24):
            self.assertEqual(da.getIJ(0,i),expected2[i])
            pass
        for pos,elt in enumerate(dabis):
            da[pos]*=elt
            pass
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        expected3=[-3600,-4949,-6400,-8181,-10000,-12221,-14400,-17069,-19600,-22725,-25600,-29189,-32400,-36461,-40000,-44541,-48400,-53429,-57600,-63125,-67600,-73629,-78400,-84941.0]
        for i in xrange(24):
            self.assertEqual(da.getIJ(0,i),expected3[i])
            pass
        # operator/=
        da=DataArrayInt.New(range(6,30),12,2)
        da2=DataArrayInt.New(range(1,13),12,1)
        dabis=-da
        da/=da2
        expected1=[6,7,4,4,3,3,3,3,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
        for i in xrange(24):
            self.assertEqual(da.getIJ(0,i),expected1[i])
            pass
        da=-dabis
        da/=DataArrayInt.New([2,3],1,2)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        expected2=[3,2,4,3,5,3,6,4,7,5,8,5,9,6,10,7,11,7,12,8,13,9,14,9]
        for i in xrange(24):
            self.assertEqual(da.getIJ(0,i),expected2[i])
            pass
        pass

    def testSwigDADOp5(self):
        da=DataArrayDouble.New([5,6,7,8,9,6,7,-2,3,9,8,10])
        da.rearrange(3)
        da2=DataArrayDouble.New([5.,8.,10.,12])
        self.assertEqual(4,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        da3=da+da2
        self.assertEqual(4,da3.getNumberOfTuples());
        self.assertEqual(3,da3.getNumberOfComponents());
        expected1=[10.,11.,12.,16.,17.,14.,17.,8.,13.,21.,20.,22.]
        for i in xrange(12):
            self.assertAlmostEqual(da3.getIJ(0,i),expected1[i],13)
            pass
        da3=da2+da
        self.assertEqual(4,da3.getNumberOfTuples());
        self.assertEqual(3,da3.getNumberOfComponents());
        for i in xrange(12):
            self.assertAlmostEqual(da3.getIJ(0,i),expected1[i],13)
            pass
        # Test new API of classmethod DataArrayDouble.New
        vals=[5,6,7,8,9,6,7,-2,3,9,8,10]
        da=DataArrayDouble.New(vals)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        for i in xrange(12):
            self.assertAlmostEqual(da.getIJ(0,i),vals[i],13)
            pass
        da=DataArrayDouble.New(vals,12)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        for i in xrange(12):
            self.assertAlmostEqual(da.getIJ(0,i),vals[i],13)
            pass
        da=DataArrayDouble.New(vals,1,12)
        self.assertEqual(1,da.getNumberOfTuples());
        self.assertEqual(12,da.getNumberOfComponents());
        for i in xrange(12):
            self.assertAlmostEqual(da.getIJ(0,i),vals[i],13)
            pass
        da=DataArrayDouble.New(vals,6,2)
        self.assertEqual(6,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in xrange(12):
            self.assertAlmostEqual(da.getIJ(0,i),vals[i],13)
            pass
        da=DataArrayDouble.New(vals,4,3)
        self.assertEqual(4,da.getNumberOfTuples());
        self.assertEqual(3,da.getNumberOfComponents());
        for i in xrange(12):
            self.assertAlmostEqual(da.getIJ(0,i),vals[i],13)
            pass
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,11);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,13);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,5,2);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,7,2);
        pass

    def testSwigDADOp6(self):
        da=DataArrayInt.New([5,6,7,8,9,6,7,-2,3,9,8,10])
        da.rearrange(3)
        da2=DataArrayInt.New([5,8,10,12])
        self.assertEqual(4,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        da3=da+da2
        self.assertEqual(4,da3.getNumberOfTuples());
        self.assertEqual(3,da3.getNumberOfComponents());
        expected1=[10,11,12,16,17,14,17,8,13,21,20,22]
        for i in xrange(12):
            self.assertEqual(da3.getIJ(0,i),expected1[i])
            pass
        da3=da2+da
        self.assertEqual(4,da3.getNumberOfTuples());
        self.assertEqual(3,da3.getNumberOfComponents());
        for i in xrange(12):
            self.assertEqual(da3.getIJ(0,i),expected1[i])
            pass
        da3=da+DataArrayInt.New(da2.getValues())
        # Test new API of classmethod DataArrayInt.New
        vals=[5,6,7,8,9,6,7,-2,3,9,8,10]
        da=DataArrayDouble.New(vals)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        for i in xrange(12):
            self.assertEqual(da.getIJ(0,i),vals[i])
            pass
        da=DataArrayDouble.New(vals,12)
        self.assertEqual(12,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        for i in xrange(12):
            self.assertEqual(da.getIJ(0,i),vals[i])
            pass
        da=DataArrayDouble.New(vals,1,12)
        self.assertEqual(1,da.getNumberOfTuples());
        self.assertEqual(12,da.getNumberOfComponents());
        for i in xrange(12):
            self.assertEqual(da.getIJ(0,i),vals[i])
            pass
        da=DataArrayDouble.New(vals,6,2)
        self.assertEqual(6,da.getNumberOfTuples());
        self.assertEqual(2,da.getNumberOfComponents());
        for i in xrange(12):
            self.assertEqual(da.getIJ(0,i),vals[i])
            pass
        da=DataArrayDouble.New(vals,4,3)
        self.assertEqual(4,da.getNumberOfTuples());
        self.assertEqual(3,da.getNumberOfComponents());
        for i in xrange(12):
            self.assertEqual(da.getIJ(0,i),vals[i])
            pass
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,11);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,13);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,5,2);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,vals,7,2);
        pass

    def testSwigDADOp9(self):
        l1=[(1.,2.,3),(4.,5.,6.),(7.,8.,9.),[10.,11.,12.]]
        da1=DataArrayDouble(l1,4,3)
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        da2=DataArrayDouble(12) ; da2.iota(1.) ; da2.rearrange(3)
        self.assertTrue(da2.isEqual(da1,1e-12))
        self.assertRaises(InterpKernelException,DataArrayDouble.New,l1,3,4);
        da3=DataArrayDouble(l1,4)
        self.assertTrue(da3.isEqual(da1,1e-12))
        self.assertRaises(InterpKernelException,DataArrayDouble.New,l1,3);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,l1,5);
        l1=[(1.,2.,3),(4.,(5.),((6.))),(7.,8.,9.),[10.,11.,12.]]
        da1=DataArrayDouble(l1,4,3)
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        da2=DataArrayDouble(12) ; da2.iota(1.) ; da2.rearrange(3)
        self.assertTrue(da2.isEqual(da1,1e-12))
        self.assertRaises(InterpKernelException,DataArrayDouble.New,l1,3,4);
        da3=DataArrayDouble(l1,4)
        self.assertTrue(da3.isEqual(da1,1e-12))
        self.assertRaises(InterpKernelException,DataArrayDouble.New,l1,3);
        self.assertRaises(InterpKernelException,DataArrayDouble.New,l1,5);
        #
        l1=[(1,2,3),(4,5,6),(7,8,9),[10,11,12]]
        da1=DataArrayInt(l1,4,3)
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        da2=DataArrayInt(12) ; da2.iota(1) ; da2.rearrange(3)
        self.assertTrue(da2.isEqual(da1))
        self.assertRaises(InterpKernelException,DataArrayInt.New,l1,3,4);
        da3=DataArrayInt(l1,4)
        self.assertTrue(da3.isEqual(da1))
        self.assertRaises(InterpKernelException,DataArrayInt.New,l1,3);
        self.assertRaises(InterpKernelException,DataArrayInt.New,l1,5);
        l1=[(1,[2],3),(4,[(5)],6),((([7])),8,9),[10,11,12]]
        da1=DataArrayInt(l1,4,3)
        self.assertEqual(4,da1.getNumberOfTuples());
        self.assertEqual(3,da1.getNumberOfComponents());
        da2=DataArrayInt(12) ; da2.iota(1) ; da2.rearrange(3)
        self.assertTrue(da2.isEqual(da1))
        self.assertRaises(InterpKernelException,DataArrayInt.New,l1,3,4);
        da3=DataArrayInt(l1,4)
        self.assertTrue(da3.isEqual(da1))
        self.assertRaises(InterpKernelException,DataArrayInt.New,l1,3);
        self.assertRaises(InterpKernelException,DataArrayInt.New,l1,5);
        pass

    def testRenumberNodesInConn1(self):
        mesh2DCoords=[-0.3,-0.3,0., 0.2,-0.3,0., 0.7,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.7,0.2,0., -0.3,0.7,0., 0.2,0.7,0., 0.7,0.7,0. ]
        mesh2DConn=[1,4,2, 4,5,2, 0,3,4,1, 6,7,4,3, 7,8,5,4]
        mesh2D=MEDCouplingUMesh.New("mesh",2);
        mesh2D.allocateCells(5);
        mesh2D.insertNextCell(NORM_TRI3,3,mesh2DConn[0:3])
        mesh2D.insertNextCell(NORM_TRI3,3,mesh2DConn[3:6])
        mesh2D.insertNextCell(NORM_QUAD4,4,mesh2DConn[6:10])
        mesh2D.insertNextCell(NORM_QUAD4,4,mesh2DConn[10:14])
        mesh2D.insertNextCell(NORM_QUAD4,4,mesh2DConn[14:18])
        mesh2D.finishInsertingCells();
        myCoords=DataArrayDouble.New(mesh2DCoords,9,3);
        mesh2D.setCoords(myCoords);
        mesh2D.checkCoherency();
        #
        mesh3DCoords=[-0.3,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.2,-0.3,0., -0.3,-0.3,1., -0.3,0.2,1., 0.2,0.2,1., 0.2,-0.3,1. ]
        mesh3DConn=[0,1,2,3,4,5,6,7]
        mesh3D=MEDCouplingUMesh.New("mesh",3);
        mesh3D.allocateCells(1);
        mesh3D.insertNextCell(NORM_HEXA8,8,mesh3DConn[:])
        mesh3D.finishInsertingCells();
        myCoords3D=DataArrayDouble.New(mesh3DCoords,8,3);
        mesh3D.setCoords(myCoords3D);
        mesh3D.checkCoherency();
        #
        mesh3D_2=mesh3D.deepCpy();
        mesh2D_2=mesh2D.deepCpy();
        mesh3D_4=mesh3D.deepCpy();
        mesh2D_4=mesh2D.deepCpy();
        oldNbOf3DNodes=mesh3D.getNumberOfNodes();
        renumNodes=DataArrayInt.New();
        renumNodes.alloc(mesh2D.getNumberOfNodes(),1);
        renumNodes.iota(oldNbOf3DNodes);
        coo=DataArrayDouble.Aggregate(mesh3D.getCoords(),mesh2D.getCoords());
        mesh3D.setCoords(coo);
        mesh2D.setCoords(coo);
        mesh2DCpy=mesh2D.deepCpy()
        mesh2D_3=mesh2D.deepCpy();
        mesh2D_3.shiftNodeNumbersInConn(oldNbOf3DNodes);
        mesh2D.renumberNodesInConn(renumNodes);
        mesh2DCpy.renumberNodesInConn(renumNodes.getValues());
        self.assertTrue(mesh2D.isEqual(mesh2DCpy,1e-12))
        self.assertTrue(mesh2D.isEqual(mesh2D_3,1e-12))
        #
        da1,da2=mesh3D.checkGeoEquivalWith(mesh3D_2,10,1e-12);
        self.assertTrue(da1==None);
        self.assertEqual(8,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        expected1=[8,11,12,9,4,5,6,7]
        for i in xrange(8):
            self.assertEqual(expected1[i],da2.getIJ(i,0));
            pass
        #
        da1,da2=mesh2D.checkGeoEquivalWith(mesh2D_2,10,1e-12);
        self.assertTrue(da1==None);
        self.assertEqual(9,da2.getNumberOfTuples());
        self.assertEqual(1,da2.getNumberOfComponents());
        for i in xrange(9):
            self.assertEqual(8+i,da2.getIJ(i,0));
            pass
        #
        mesh2D_5=mesh2D_4.deepCpy();
        mesh2D_5.translate([1.,0.,0.]);
        meshes=[mesh3D_4,mesh2D_4,mesh2D_5];
        MEDCouplingUMesh.PutUMeshesOnSameAggregatedCoords(meshes);
        self.assertTrue(mesh3D_4.getCoords().getHiddenCppPointer()==mesh2D_4.getCoords().getHiddenCppPointer());
        self.assertTrue(mesh2D_4.getCoords().getHiddenCppPointer()==mesh2D_5.getCoords().getHiddenCppPointer());
        mesh3D_4.checkCoherency(); mesh2D_4.checkCoherency(); mesh2D_5.checkCoherency();
        self.assertEqual(26,mesh3D_4.getNumberOfNodes());
        self.assertEqual(3,mesh3D_4.getSpaceDimension());
        self.assertEqual(9,mesh3D_4.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(23,mesh2D_4.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(23,mesh2D_5.getNodalConnectivity().getNumberOfTuples());
        expected2=[18,0,1,2,3,4,5,6,7]
        expected3=[3,9,12,10, 3,12,13,10, 4,8,11,12,9, 4,14,15,12,11, 4,15,16,13,12]
        expected4=[3,18,21,19, 3,21,22,19, 4,17,20,21,18, 4,23,24,21,20, 4,24,25,22,21]
        expected5=[-0.3,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.2,-0.3,0., -0.3,-0.3,1., -0.3,0.2,1., 0.2,0.2,1., 0.2,-0.3,1., -0.3,-0.3,0., 0.2,-0.3,0., 0.7,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.7,0.2,0., -0.3,0.7,0., 0.2,0.7,0., 0.7,0.7,0., 0.7, -0.3, 0.0, 1.2, -0.3, 0.0, 1.7, -0.3, 0.0, 0.7, 0.2, 0.0, 1.2, 0.2, 0.0, 1.7, 0.2, 0.0, 0.7, 0.7, 0.0, 1.2, 0.7, 0.0, 1.7, 0.7, 0.0]
        self.assertEqual(expected2,mesh3D_4.getNodalConnectivity().getValues());
        self.assertEqual(expected3,mesh2D_4.getNodalConnectivity().getValues());
        self.assertEqual(expected4,mesh2D_5.getNodalConnectivity().getValues());
        for i in xrange(78):
            self.assertAlmostEqual(expected5[i],mesh3D_4.getCoords().getIJ(0,i),12);
            pass
        #
        MEDCouplingUMesh.MergeNodesOnUMeshesSharingSameCoords(meshes,1e-12);
        mesh3D_4.checkCoherency(); mesh2D_4.checkCoherency(); mesh2D_5.checkCoherency();
        self.assertTrue(mesh3D_4.getCoords().getHiddenCppPointer()==mesh2D_4.getCoords().getHiddenCppPointer());
        self.assertTrue(mesh2D_4.getCoords().getHiddenCppPointer()==mesh2D_5.getCoords().getHiddenCppPointer());
        self.assertEqual(19,mesh3D_4.getNumberOfNodes());
        self.assertEqual(3,mesh3D_4.getSpaceDimension());
        self.assertEqual(9,mesh3D_4.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(23,mesh2D_4.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(23,mesh2D_5.getNodalConnectivity().getNumberOfTuples());
        expected6=[18,0,1,2,3,4,5,6,7]
        expected7=[3,3,2,8, 3,2,9,8, 4,0,1,2,3, 4,10,11,2,1, 4,11,12,9,2]
        expected8=[3,13,15,14, 3,15,16,14, 4,8,9,15,13, 4,12,17,15,9, 4,17,18,16,15]
        expected9=[-0.3, -0.3, 0., -0.3, 0.2, 0., 0.2, 0.2, 0., 0.2, -0.3, 0., -0.3, -0.3, 1., -0.3, 0.2, 1.,
                    0.2, 0.2, 1., 0.2, -0.3, 1., 0.7, -0.3, 0., 0.7, 0.2, 0., -0.3, 0.7, 0., 0.2, 0.7, 0.,
                    0.7, 0.7, 0., 1.2, -0.3, 0., 1.7, -0.3, 0., 1.2, 0.2, 0., 1.7, 0.2, 0., 1.2, 0.7, 0., 1.7, 0.7, 0.]
        self.assertEqual(expected6,mesh3D_4.getNodalConnectivity().getValues());
        self.assertEqual(expected7,mesh2D_4.getNodalConnectivity().getValues());
        self.assertEqual(expected8,mesh2D_5.getNodalConnectivity().getValues());
        for i in xrange(57):
            self.assertAlmostEqual(expected9[i],mesh3D_4.getCoords().getIJ(0,i),12);
            pass
        #
        pass
    
    def testComputeNeighborsOfCells1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        d1,d2=m.computeNeighborsOfCells();
        self.assertEqual(6,d2.getNumberOfTuples());
        self.assertEqual(10,d1.getNumberOfTuples());
        expected1=[0,2,4,6,8,10]
        expected2=[3,1,0,2,4,1,4,0,2,3]
        self.assertEqual(expected1,d2.getValues());
        self.assertEqual(expected2,d1.getValues());
        pass

    def testCheckButterflyCellsBug1(self):
        mesh2DCoords=[323.85,120.983748908684,317.5,131.982271536747,336.55,120.983748908686,330.2,131.982271536751,323.85,142.98079416481]
        mesh2DConn=[4,1,0,2,3]
        mesh2D=MEDCouplingUMesh.New("mesh",2);
        mesh2D.allocateCells(1);
        mesh2D.insertNextCell(NORM_POLYGON,5,mesh2DConn[0:5])
        mesh2D.finishInsertingCells();
        myCoords=DataArrayDouble.New(mesh2DCoords,5,2);
        mesh2D.setCoords(myCoords);
        mesh2D.checkCoherency();
        #
        v=mesh2D.checkButterflyCells();
        self.assertTrue(v.empty());
        pass

    def testDataArrayIntRange1(self):
        d=DataArrayInt.Range(2,17,7);
        expected1=[2,9,16]
        self.assertEqual(3,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(expected1,d.getValues());
        #
        d=DataArrayInt.Range(2,23,7);
        self.assertEqual(3,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(expected1,d.getValues());
        #
        d=DataArrayInt.Range(2,24,7);
        expected2=[2,9,16,23]
        self.assertEqual(4,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(expected2,d.getValues());
        #
        d=DataArrayInt.Range(24,2,-7);
        expected3=[24,17,10,3]
        self.assertEqual(4,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(expected3,d.getValues());
        #
        d=DataArrayInt.Range(23,2,-7);
        expected4=[23,16,9]
        self.assertEqual(3,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(expected4,d.getValues());
        #
        d=DataArrayInt.Range(23,22,-7);
        self.assertEqual(1,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(23,d.getIJ(0,0));
        #
        d=DataArrayInt.Range(22,23,7);
        self.assertEqual(1,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        self.assertEqual(22,d.getIJ(0,0));
        #
        d=DataArrayInt.Range(22,22,7);
        self.assertEqual(0,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        #
        d=DataArrayInt.Range(22,22,-7);
        self.assertEqual(0,d.getNumberOfTuples());
        self.assertEqual(1,d.getNumberOfComponents());
        #
        self.assertRaises(InterpKernelException,DataArrayInt.Range,22,23,-7);
        self.assertRaises(InterpKernelException,DataArrayInt.Range,23,22,7);
        self.assertRaises(InterpKernelException,DataArrayInt.Range,23,22,0);
        self.assertRaises(InterpKernelException,DataArrayInt.Range,22,23,0);
        pass

    def testSwigUMeshGetItem1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        subMesh=m.buildPartOfMySelf([1,3],True);
        self.assertTrue(isinstance(subMesh,MEDCouplingUMesh))
        m1=m[[1,3]]
        self.assertTrue(isinstance(m1,MEDCouplingUMesh))
        m2=m[(1,3)]
        self.assertTrue(isinstance(m2,MEDCouplingUMesh))
        m3=m[1::2]
        self.assertTrue(isinstance(m3,MEDCouplingUMesh))
        m4=m[DataArrayInt.New([1,3])]
        m5_1=m[1]
        self.assertTrue(isinstance(m5_1,MEDCouplingUMesh))
        m5_2=m[3]
        self.assertTrue(isinstance(m5_2,MEDCouplingUMesh))
        m5=MEDCouplingUMesh.MergeUMeshesOnSameCoords([m5_1,m5_2]);
        m5.setName(subMesh.getName())
        self.assertTrue(isinstance(m4,MEDCouplingUMesh))
        self.assertTrue(subMesh.isEqual(m1,1e-12))
        self.assertTrue(subMesh.isEqual(m2,1e-12))
        self.assertTrue(subMesh.isEqual(m3,1e-12))
        self.assertTrue(subMesh.isEqual(m4,1e-12))
        self.assertTrue(subMesh.isEqual(m5,1e-12))
        self.assertRaises(InterpKernelException,m.buildPartOfMySelf,[1,5],True);
        pass
    
    def testSwigGetItem3(self):
        da=DataArrayInt.New([4,5,6])
        self.assertEqual(5,da[1])
        self.assertEqual(6,da[-1])
        self.assertRaises(InterpKernelException,da.__getitem__,3)
        da=DataArrayInt.New([4,5,6,7,8,9],2,3)
        self.assertEqual(9,da[1,2])
        da=DataArrayDouble.New([4.1,5.2,6.3])
        self.assertAlmostEqual(5.2,da[1],12)
        self.assertAlmostEqual(6.3,da[-1],12)
        self.assertRaises(InterpKernelException,da.__getitem__,3)
        da=DataArrayDouble.New([4.12,5.12,6.12,7.12,8.12,9.12],2,3)
        self.assertAlmostEqual(9.12,da[1,2],12)
        pass

    def testSwigDADISub1(self):
        mesh3D,mesh2D=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        bary=mesh3D.getBarycenterAndOwner()
        bary=bary[:,:2]
        pts=bary.getDifferentValues(1e-12)
        expected=[[0,6,12],[1,7,13],[2,8,14],[3,9,15],[4,10,16],[5,11,17]]
        for pos,pt in enumerate(pts):
            bary2=bary[:,:2]
            bary2[:]-=pt
            norm=bary2.magnitude()
            self.assertEqual(expected[pos],norm.getIdsInRange(-1.,1e-5).getValues())
            pass
        expected2=[[3.,54.],[-141.,180.],[21.,54.],[39.,72.],[-15.,90.],[21.,90.]]
        for pos,pt in enumerate(pts):
            bary2=bary[:,:2]
            bary2[:]+=pt
            self.assertAlmostEqual(expected2[pos][0],bary2.accumulate()[0],12);
            self.assertAlmostEqual(expected2[pos][1],bary2.accumulate()[1],12);
            pass
        expected3=[[-3.,22.5],[45.,337.5],[-9., 22.5],[-15.,67.5],[3.,112.5],[-9.,112.5]]
        for pos,pt in enumerate(pts):
            bary2=bary[:,:2]
            bary2[:]*=pt
            self.assertAlmostEqual(expected3[pos][0],bary2.accumulate()[0],12);
            self.assertAlmostEqual(expected3[pos][1],bary2.accumulate()[1],12);
            pass
        expected4=[[-12.,90.],[0.8,6.],[-4,90.],[-2.4,30.],[12.,18],[-4,18.]]
        for pos,pt in enumerate(pts):
            bary2=bary[:,:2]
            bary2[:]/=pt
            self.assertAlmostEqual(expected4[pos][0],bary2.accumulate()[0],12);
            self.assertAlmostEqual(expected4[pos][1],bary2.accumulate()[1],12);
            pass
        #
        d=DataArrayInt.New([1,2,0,1,0,2],3,2)
        e=DataArrayInt.New([1,11,101,2,12,102,3,13,103,4,14,104],4,3)
        expected5=[[1,11,101,77,77,77,77,77,77,4,14,104],[77,77,77,77,77,77,3,13,103,4,14,104],[77,77,77,2,12,102,77,77,77,4,14,104]]
        expected6=[[1,77,77,2,77,77,3,77,77,4,77,77],[77,77,101,77,77,102,77,77,103,77,77,104],[77,11,77,77,12,77,77,13,77,77,14,77]]
        for pos,tup in enumerate(d):
            f=e[:]
            self.assertTrue(isinstance(f,DataArrayInt))
            f[tup]=77
            self.assertEqual(expected5[pos],f.getValues())
            self.assertEqual(6*[77],f[tup].getValues())
            f=e[:]
            f[:,tup]=77
            self.assertEqual(expected6[pos],f.getValues())
            self.assertEqual(8*[77],f[:,tup].getValues())
            pass
        #
        e=e.convertToDblArr()
        for pos,tup in enumerate(d):
            f=e[:]
            self.assertTrue(isinstance(f,DataArrayDouble))
            f[tup]=77.
            self.assertEqual(expected5[pos],f.convertToIntArr().getValues())
            self.assertEqual(6*[77],f[tup].convertToIntArr().getValues())
            f=e[:]
            f[:,tup]=77.
            self.assertEqual(expected6[pos],f.convertToIntArr().getValues())
            self.assertEqual(8*[77],f[:,tup].convertToIntArr().getValues())
            pass
        pass

    def testDataArrayDoubleGetMinMaxPerComponent1(self):
        values1=[1.,2.,3.,-0.9,2.1,3.,1.3,1.7,3.,1.,1.8,3.]
        d1=DataArrayDouble.New();
        self.assertRaises(InterpKernelException,d1.getMinMaxPerComponent)
        d1=DataArrayDouble.New(values1,4,3);
        res=d1.getMinMaxPerComponent();
        self.assertTrue(isinstance(res,list))
        self.assertEqual(3,len(res))
        for i in xrange(3):
            self.assertTrue(isinstance(res[i],tuple))
            self.assertEqual(2,len(res[i]))
            pass
        expected1=[-0.9,1.3,1.7,2.1,3.,3.]
        for i in xrange(6):
            self.assertAlmostEqual(expected1[i],res[i/2][i%2],14)
            pass
        #
        d1.rearrange(2);
        res=d1.getMinMaxPerComponent();
        self.assertTrue(isinstance(res,list))
        self.assertEqual(2,len(res))
        for i in xrange(2):
            self.assertTrue(isinstance(res[i],tuple))
            self.assertEqual(2,len(res[i]))
            pass
        expected2=[1.,3.,-0.9,3.]
        for i in xrange(4):
            self.assertAlmostEqual(expected2[i],res[i/2][i%2],14)
            pass
        #
        d1.rearrange(1);
        res=d1.getMinMaxPerComponent();
        self.assertTrue(isinstance(res,list))
        self.assertEqual(1,len(res))
        for i in xrange(1):
            self.assertTrue(isinstance(res[i],tuple))
            self.assertEqual(2,len(res[i]))
            pass
        expected3=[-0.9,3.]
        for i in xrange(2):
            self.assertAlmostEqual(expected3[i],res[i/2][i%2],14)
            pass
        pass

    def testDataArrayIntGetHashCode1(self):
        d1=DataArrayInt.New(range(3545))
        d2=DataArrayInt.New(range(3545))
        self.assertEqual(d2.getHashCode(),d1.getHashCode())
        self.assertEqual(232341068,d1.getHashCode())
        d1[886]=6
        self.assertEqual(232340188,d1.getHashCode())
        pass

    def testZipConnectivityPol1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells1=[2,3,4]
        m2_1=m1.buildPartOfMySelf(cells1,True);
        m2=m2_1
        self.assertTrue(isinstance(m2,MEDCouplingUMesh))
        # no permutation policy 0
        isOk,arr=m1.areCellsIncludedIn(m2,0)
        self.assertTrue(isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells1,arr.getValues())
        # no permutation policy 1
        isOk,arr=m1.areCellsIncludedIn(m2,1)
        self.assertTrue(isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells1,arr.getValues())
        # no permutation policy 2
        isOk,arr=m1.areCellsIncludedIn(m2,2)
        self.assertTrue(isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells1,arr.getValues())
        # some modification into m2
        modif1=[2,4,5]
        m2.getNodalConnectivity()[1:4]=modif1
        #policy 0 fails because cell0 in m2 has same orientation be not same connectivity
        expected1=[5,3,4]
        isOk,arr=m1.areCellsIncludedIn(m2,0)
        self.assertTrue(not isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(expected1,arr.getValues())
        #policy 1 succeeds because cell0 in m2 has not exactly the same conn
        isOk,arr=m1.areCellsIncludedIn(m2,1)
        self.assertTrue(isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells1,arr.getValues())
        #policy 2 succeeds because cell0 in m2 has same nodes in connectivity
        isOk,arr=m1.areCellsIncludedIn(m2,2)
        self.assertTrue(isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells1,arr.getValues())
        #some new modification into m2
        modif2=[2,5,4]
        m2.getNodalConnectivity()[1:4]=modif2
        #policy 0 fails because cell0 in m2 has not exactly the same conn
        isOk,arr=m1.areCellsIncludedIn(m2,0)
        self.assertTrue(not isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(expected1,arr.getValues())
        #policy 1 fails too because cell0 in m2 has not same orientation
        isOk,arr=m1.areCellsIncludedIn(m2,1)
        self.assertTrue(not isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(expected1,arr.getValues())
        #policy 2 succeeds because cell0 in m2 has same nodes in connectivity
        isOk,arr=m1.areCellsIncludedIn(m2,2)
        self.assertTrue(isOk);
        self.assertEqual(3,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells1,arr.getValues())
        # Now 1D
        cells2=[3,2]
        m1=MEDCouplingDataForTest.build1DSourceMesh_2();
        m2_1=m1.buildPartOfMySelf(cells2,True);
        m2=m2_1
        self.assertTrue(isinstance(m2,MEDCouplingUMesh))
        # no permutation policy 0
        isOk,arr=m1.areCellsIncludedIn(m2,0)
        self.assertTrue(isOk);
        self.assertEqual(2,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells2,arr.getValues())
        # no permutation policy 1
        isOk,arr=m1.areCellsIncludedIn(m2,1)
        self.assertTrue(isOk);
        self.assertEqual(2,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells2,arr.getValues())
        # no permutation policy 2
        isOk,arr=m1.areCellsIncludedIn(m2,2)
        self.assertTrue(isOk);
        self.assertEqual(2,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells2,arr.getValues())
        # some modification into m2
        modif3=[4,3]
        m2.getNodalConnectivity()[1:3]=modif3
        #policy 0 fails because cell0 in m2 has not exactly the same conn
        expected2=[4,2]
        isOk,arr=m1.areCellsIncludedIn(m2,0)
        self.assertTrue(not isOk);
        self.assertEqual(2,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(expected2,arr.getValues())
        #policy 1 fails too because cell0 in m2 has not same orientation
        isOk,arr=m1.areCellsIncludedIn(m2,1)
        self.assertTrue(not isOk);
        self.assertEqual(2,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(expected2,arr.getValues())
        #policy 2 succeeds because cell0 in m2 has same nodes in connectivity
        isOk,arr=m1.areCellsIncludedIn(m2,2)
        self.assertTrue(isOk);
        self.assertEqual(2,arr.getNumberOfTuples());
        self.assertEqual(1,arr.getNumberOfComponents());
        self.assertEqual(cells2,arr.getValues())
        pass

    def toSeeIfDaIIopsAreOK(self,d):
        d+=5
        d*=6
        d/=3
        d-=2
        d%=7
        pass
        
    def testSwigDAIOp5(self):
        d=DataArrayInt.New([4,5,6,10,3,-1],2,3)
        self.toSeeIfDaIIopsAreOK(d)
        dExp=DataArrayInt.New([2,4,6,0,0,6],2,3)
        self.assertTrue(d.isEqual(dExp));
        pass
    
    def toSeeIfDaDIopsAreOK(self,d):
        d+=5
        d*=6
        d/=3
        d-=2
        pass

    def testSwigDADOp7(self):
        d=DataArrayDouble.New([4.,5.,6.,10.,3.,-1.],2,3)
        self.toSeeIfDaDIopsAreOK(d)
        dExp=DataArrayDouble.New([16.,18.,20.,28.,14.,6.],2,3)
        self.assertTrue(d.isEqual(dExp,1e-14));
        pass

    def testConvexEnvelop2D1(self):
        coords=[7.54758495819e-14,-1.12270326253e-12,8.43143594193,-1.02835845055e-12,4.21571797096,7.30183771609,-4.21571797097,7.30183771609,-8.43143594193,-1.09439981894e-12,-4.21571797097,-7.30183771609,4.21571797097,-7.30183771609,16.8628718839,-1.02835845055e-12,12.6471539129,7.30183771609,8.43143594193,14.6036754322,2.26427548746e-13,14.6036754322,-8.43143594193,14.6036754322,-12.6471539129,7.30183771609,-16.8628718839,-1.39630321727e-12,-12.6471539129,-7.30183771609,-8.43143594193,-14.6036754322,3.7737924791e-14,-14.6036754322,8.43143594193,-14.6036754322,12.6471539129,-7.30183771609,25.2943078258,-1.07553085654e-12,21.0785898548,7.30183771609,16.8628718839,14.6036754322,12.6471539129,21.9055131483,4.21571797096,21.9055131483,-4.21571797097,21.9055131483,-12.6471539129,21.9055131483,-16.8628718839,14.6036754322,-21.0785898548,7.30183771609,-25.2943078258,-1.02835845055e-12,-21.0785898548,-7.30183771609,-16.8628718839,-14.6036754322,-12.6471539129,-21.9055131483,-4.21571797097,-21.9055131483,4.21571797097,-21.9055131483,12.6471539129,-21.9055131483,16.8628718839,-14.6036754322,21.0785898548,-7.30183771609,33.7257437677,-7.45324014622e-13,29.5100257968,7.30183771609,25.2943078258,14.6036754322,21.0785898548,21.9055131483,16.8628718839,29.2073508644,8.43143594193,29.2073508644,-1.20761359331e-12,29.2073508644,-8.43143594193,29.2073508644,-16.8628718839,29.2073508644,-21.0785898548,21.9055131483,-25.2943078258,14.6036754322,-29.5100257968,7.30183771609,-33.7257437677,-7.26455052226e-13,-29.5100257968,-7.30183771609,-25.2943078258,-14.6036754322,-21.0785898548,-21.9055131483,-16.8628718839,-29.2073508644,-8.43143594193,-29.2073508644,4.15117172701e-13,-29.2073508644,8.43143594193,-29.2073508644,16.8628718839,-29.2073508644,21.0785898548,-21.9055131483,25.2943078258,-14.6036754322,29.5100257968,-7.30183771609,42.1571797097,-1.86802727715e-12,37.9414617387,7.30183771609,33.7257437677,14.6036754322,29.5100257968,21.9055131483,25.2943078258,29.2073508644,21.0785898548,36.5091885805,12.6471539129,36.5091885805,4.21571797096,36.5091885805,-4.21571797096,36.5091885805,-12.6471539129,36.5091885805,-21.0785898548,36.5091885805,-25.2943078258,29.2073508644,-29.5100257968,21.9055131483,-33.7257437677,14.6036754322,-37.9414617387,7.30183771609,-42.1571797097,-9.81186044565e-13,-37.9414617387,-7.30183771609,-33.7257437677,-14.6036754322,-29.5100257968,-21.9055131483,-25.2943078258,-29.2073508644,-21.0785898548,-36.5091885805,-12.6471539129,-36.5091885805,-4.21571797097,-36.5091885805,4.21571797097,-36.5091885805,12.6471539129,-36.5091885805,21.0785898548,-36.5091885805,25.2943078258,-29.2073508644,29.5100257968,-21.9055131483,33.7257437677,-14.6036754322,37.9414617387,-7.30183771609,50.5886156516,-6.98151608633e-13,46.3728976806,7.30183771609,42.1571797097,14.6036754322,37.9414617387,21.9055131483,33.7257437677,29.2073508644,29.5100257968,36.5091885805,25.2943078258,43.8110262966,16.8628718839,43.8110262966,8.43143594193,43.8110262966,-1.84915831476e-12,43.8110262966,-8.43143594193,43.8110262966,-16.8628718839,43.8110262966,-25.2943078258,43.8110262966,-29.5100257968,36.5091885805,-33.7257437677,29.2073508644,-37.9414617387,21.9055131483,-42.1571797097,14.6036754322,-46.3728976806,7.30183771609,-50.5886156516,-1.47177906685e-12,-46.3728976806,-7.30183771609,-42.1571797097,-14.6036754322,-37.9414617387,-21.9055131483,-33.7257437677,-29.2073508644,-29.5100257968,-36.5091885805,-25.2943078258,-43.8110262966,-16.8628718839,-43.8110262966,-8.43143594193,-43.8110262966,7.54758495819e-14,-43.8110262966,8.43143594193,-43.8110262966,16.8628718839,-43.8110262966,25.2943078258,-43.8110262966,29.5100257968,-36.5091885805,33.7257437677,-29.2073508644,37.9414617387,-21.9055131483,42.1571797097,-14.6036754322,46.3728976806,-7.30183771609,59.0200515935,-7.9249642061e-13,54.8043336225,7.30183771609,50.5886156516,14.6036754322,46.3728976806,21.9055131483,42.1571797097,29.2073508644,37.9414617387,36.5091885805,33.7257437677,43.8110262966,29.5100257968,51.1128640127,21.0785898548,51.1128640127,12.6471539129,51.1128640127,4.21571797096,51.1128640127,-4.21571797096,51.1128640127,-12.6471539129,51.1128640127,-21.0785898548,51.1128640127,-29.5100257968,51.1128640127,-33.7257437677,43.8110262966,-37.9414617387,36.5091885805,-42.1571797097,29.2073508644,-46.3728976806,21.9055131483,-50.5886156516,14.6036754322,-54.8043336226,7.30183771609,-59.0200515935,-1.31139288649e-12,-54.8043336226,-7.30183771609,-50.5886156516,-14.6036754322,-46.3728976806,-21.9055131483,-42.1571797097,-29.2073508644,-37.9414617387,-36.5091885805,-33.7257437677,-43.8110262966,-29.5100257968,-51.1128640127,-21.0785898548,-51.1128640127,-12.6471539129,-51.1128640127,-4.21571797097,-51.1128640127,4.21571797097,-51.1128640127,12.6471539129,-51.1128640127,21.0785898548,-51.1128640127,29.5100257968,-51.1128640127,33.7257437677,-43.8110262966,37.9414617387,-36.5091885805,42.1571797097,-29.2073508644,46.3728976806,-21.9055131483,50.5886156516,-14.6036754322,54.8043336225,-7.30183771609,67.4514875354,-2.14162723189e-12,63.2357695645,7.30183771609,59.0200515935,14.6036754322,54.8043336226,21.9055131483,50.5886156516,29.2073508644,46.3728976806,36.5091885805,42.1571797097,43.8110262966,37.9414617387,51.1128640127,33.7257437677,58.4147017287,25.2943078258,58.4147017287,16.8628718839,58.4147017287,8.43143594193,58.4147017287,6.79282646237e-13,58.4147017287,-8.43143594193,58.4147017287,-16.8628718839,58.4147017287,-25.2943078258,58.4147017287,-33.7257437677,58.4147017287,-37.9414617387,51.1128640127,-42.1571797097,43.8110262966,-46.3728976806,36.5091885805,-50.5886156516,29.2073508644,-54.8043336226,21.9055131483,-59.0200515935,14.6036754322,-63.2357695645,7.30183771609,-67.4514875354,-1.16044118732e-12,-63.2357695645,-7.30183771609,-59.0200515935,-14.6036754322,-54.8043336226,-21.9055131483,-50.5886156516,-29.2073508644,-46.3728976806,-36.5091885805,-42.1571797097,-43.8110262966,-37.9414617387,-51.1128640127,-33.7257437677,-58.4147017287,-25.2943078258,-58.4147017287,-16.8628718839,-58.4147017287,-8.43143594193,-58.4147017287,-5.66068871864e-14,-58.4147017287,8.43143594193,-58.4147017287,16.8628718839,-58.4147017287,25.2943078258,-58.4147017287,33.7257437677,-58.4147017287,37.9414617387,-51.1128640127,42.1571797097,-43.8110262966,46.3728976806,-36.5091885805,50.5886156516,-29.2073508644,54.8043336226,-21.9055131483,59.0200515935,-14.6036754322,63.2357695645,-7.30183771609,75.8829234774,-2.29257893105e-12,71.6672055064,7.30183771609,67.4514875354,14.6036754322,63.2357695645,21.9055131483,59.0200515935,29.2073508644,54.8043336226,36.5091885805,50.5886156516,43.8110262966,46.3728976806,51.1128640127,42.1571797097,58.4147017287,37.9414617387,65.7165394448,29.5100257968,65.7165394448,21.0785898548,65.7165394448,12.6471539129,65.7165394448,4.21571797097,65.7165394448,-4.21571797096,65.7165394448,-12.6471539129,65.7165394448,-21.0785898548,65.7165394448,-29.5100257968,65.7165394448,-37.9414617387,65.7165394448,-42.1571797097,58.4147017287,-46.3728976806,51.1128640127,-50.5886156516,43.8110262966,-54.8043336226,36.5091885805,-59.0200515935,29.2073508644,-63.2357695645,21.9055131483,-67.4514875354,14.6036754322,-71.6672055064,7.30183771609,-75.8829234774,-1.31139288649e-12,-71.6672055064,-7.30183771609,-67.4514875354,-14.6036754322,-63.2357695645,-21.9055131483,-59.0200515935,-29.2073508644,-54.8043336226,-36.5091885805,-50.5886156516,-43.8110262966,-46.3728976806,-51.1128640127,-42.1571797097,-58.4147017287,-37.9414617387,-65.7165394448,-29.5100257968,-65.7165394448,-21.0785898548,-65.7165394448,-12.6471539129,-65.7165394448,-4.21571797097,-65.7165394448,4.21571797097,-65.7165394448,12.6471539129,-65.7165394448,21.0785898548,-65.7165394448,29.5100257968,-65.7165394448,37.9414617387,-65.7165394448,42.1571797097,-58.4147017287,46.3728976806,-51.1128640127,50.5886156516,-43.8110262966,54.8043336226,-36.5091885805,59.0200515935,-29.2073508644,63.2357695645,-21.9055131483,67.4514875354,-14.6036754322,71.6672055064,-7.30183771609,84.3143594193,-1.49064802924e-12,80.0986414483,7.30183771609,75.8829234774,14.6036754322,71.6672055064,21.9055131483,67.4514875354,29.2073508644,63.2357695645,36.5091885805,59.0200515935,43.8110262966,54.8043336226,51.1128640127,50.5886156516,58.4147017287,46.3728976806,65.7165394448,42.1571797097,73.0183771609,33.7257437677,73.0183771609,25.2943078258,73.0183771609,16.8628718839,73.0183771609,8.43143594193,73.0183771609,2.0755858635e-12,73.0183771609,-8.43143594193,73.0183771609,-16.8628718839,73.0183771609,-25.2943078258,73.0183771609,-33.7257437677,73.0183771609,-42.1571797097,73.0183771609,-46.3728976806,65.7165394448,-50.5886156516,58.4147017287,-54.8043336226,51.1128640127,-59.0200515935,43.8110262966,-63.2357695645,36.5091885805,-67.4514875354,29.2073508644,-71.6672055064,21.9055131483,-75.8829234774,14.6036754322,-80.0986414483,7.30183771609,-84.3143594193,-1.11326878133e-12,-80.0986414483,-7.30183771609,-75.8829234774,-14.6036754322,-71.6672055064,-21.9055131483,-67.4514875354,-29.2073508644,-63.2357695645,-36.5091885805,-59.0200515935,-43.8110262966,-54.8043336226,-51.1128640127,-50.5886156516,-58.4147017287,-46.3728976806,-65.7165394448,-42.1571797097,-73.0183771609,-33.7257437677,-73.0183771609,-25.2943078258,-73.0183771609,-16.8628718839,-73.0183771609,-8.43143594193,-73.0183771609,-5.66068871864e-14,-73.0183771609,8.43143594193,-73.0183771609,16.8628718839,-73.0183771609,25.2943078258,-73.0183771609,33.7257437677,-73.0183771609,42.1571797097,-73.0183771609,46.3728976806,-65.7165394448,50.5886156516,-58.4147017287,54.8043336226,-51.1128640127,59.0200515935,-43.8110262966,63.2357695645,-36.5091885805,67.4514875354,-29.2073508644,71.6672055064,-21.9055131483,75.8829234774,-14.6036754322,80.0986414483,-7.3018377161]
        conn=[0,2,3,4,5,6,1,1,8,2,0,6,18,7,2,9,10,3,0,1,8,3,10,11,12,4,0,2,4,3,12,13,14,5,0,5,0,4,14,15,16,6,6,1,0,5,16,17,18,7,20,8,1,18,36,19,8,21,9,2,1,7,20,9,22,23,10,2,8,21,10,23,24,11,3,2,9,11,24,25,26,12,3,10,12,11,26,27,13,4,3,13,12,27,28,29,14,4,14,4,13,29,30,15,5,15,5,14,30,31,32,16,16,6,5,15,32,33,17,17,18,6,16,33,34,35,18,7,1,6,17,35,36,19,38,20,7,36,60,37,20,39,21,8,7,19,38,21,40,22,9,8,20,39,22,41,42,23,9,21,40,23,42,43,24,10,9,22,24,43,44,25,11,10,23,25,44,45,46,26,11,24,26,25,46,47,27,12,11,27,26,47,48,28,13,12,28,27,48,49,50,29,13,29,13,28,50,51,30,14,30,14,29,51,52,31,15,31,15,30,52,53,54,32,32,16,15,31,54,55,33,33,17,16,32,55,56,34,34,35,17,33,56,57,58,35,36,18,17,34,58,59,36,19,7,18,35,59,60,37,62,38,19,60,90,61,38,63,39,20,19,37,62,39,64,40,21,20,38,63,40,65,41,22,21,39,64,41,66,67,42,22,40,65,42,67,68,43,23,22,41,43,68,69,44,24,23,42,44,69,70,45,25,24,43,45,70,71,72,46,25,44,46,45,72,73,47,26,25,47,46,73,74,48,27,26,48,47,74,75,49,28,27,49,48,75,76,77,50,28,50,28,49,77,78,51,29,51,29,50,78,79,52,30,52,30,51,79,80,53,31,53,31,52,80,81,82,54,54,32,31,53,82,83,55,55,33,32,54,83,84,56,56,34,33,55,84,85,57,57,58,34,56,85,86,87,58,59,35,34,57,87,88,59,60,36,35,58,88,89,60,37,19,36,59,89,90,61,92,62,37,90,126,91,62,93,63,38,37,61,92,63,94,64,39,38,62,93,64,95,65,40,39,63,94,65,96,66,41,40,64,95,66,97,98,67,41,65,96,67,98,99,68,42,41,66,68,99,100,69,43,42,67,69,100,101,70,44,43,68,70,101,102,71,45,44,69,71,102,103,104,72,45,70,72,71,104,105,73,46,45,73,72,105,106,74,47,46,74,73,106,107,75,48,47,75,74,107,108,76,49,48,76,75,108,109,110,77,49,77,49,76,110,111,78,50,78,50,77,111,112,79,51,79,51,78,112,113,80,52,80,52,79,113,114,81,53,81,53,80,114,115,116,82,82,54,53,81,116,117,83,83,55,54,82,117,118,84,84,56,55,83,118,119,85,85,57,56,84,119,120,86,86,87,57,85,120,121,122,87,88,58,57,86,122,123,88,89,59,58,87,123,124,89,90,60,59,88,124,125,90,61,37,60,89,125,126,91,128,92,61,126,168,127,92,129,93,62,61,91,128,93,130,94,63,62,92,129,94,131,95,64,63,93,130,95,132,96,65,64,94,131,96,133,97,66,65,95,132,97,134,135,98,66,96,133,98,135,136,99,67,66,97,99,136,137,100,68,67,98,100,137,138,101,69,68,99,101,138,139,102,70,69,100,102,139,140,103,71,70,101,103,140,141,142,104,71,102,104,103,142,143,105,72,71,105,104,143,144,106,73,72,106,105,144,145,107,74,73,107,106,145,146,108,75,74,108,107,146,147,109,76,75,109,108,147,148,149,110,76,110,76,109,149,150,111,77,111,77,110,150,151,112,78,112,78,111,151,152,113,79,113,79,112,152,153,114,80,114,80,113,153,154,115,81,115,81,114,154,155,156,116,116,82,81,115,156,157,117,117,83,82,116,157,158,118,118,84,83,117,158,159,119,119,85,84,118,159,160,120,120,86,85,119,160,161,121,121,122,86,120,161,162,163,122,123,87,86,121,163,164,123,124,88,87,122,164,165,124,125,89,88,123,165,166,125,126,90,89,124,166,167,126,91,61,90,125,167,168,127,170,128,91,168,216,169,128,171,129,92,91,127,170,129,172,130,93,92,128,171,130,173,131,94,93,129,172,131,174,132,95,94,130,173,132,175,133,96,95,131,174,133,176,134,97,96,132,175,134,177,178,135,97,133,176,135,178,179,136,98,97,134,136,179,180,137,99,98,135,137,180,181,138,100,99,136,138,181,182,139,101,100,137,139,182,183,140,102,101,138,140,183,184,141,103,102,139,141,184,185,186,142,103,140,142,141,186,187,143,104,103,143,142,187,188,144,105,104,144,143,188,189,145,106,105,145,144,189,190,146,107,106,146,145,190,191,147,108,107,147,146,191,192,148,109,108,148,147,192,193,194,149,109,149,109,148,194,195,150,110,150,110,149,195,196,151,111,151,111,150,196,197,152,112,152,112,151,197,198,153,113,153,113,152,198,199,154,114,154,114,153,199,200,155,115,155,115,154,200,201,202,156,156,116,115,155,202,203,157,157,117,116,156,203,204,158,158,118,117,157,204,205,159,159,119,118,158,205,206,160,160,120,119,159,206,207,161,161,121,120,160,207,208,162,162,163,121,161,208,209,210,163,164,122,121,162,210,211,164,165,123,122,163,211,212,165,166,124,123,164,212,213,166,167,125,124,165,213,214,167,168,126,125,166,214,215,168,127,91,126,167,215,216,169,218,170,127,216,270,217,170,219,171,128,127,169,218,171,220,172,129,128,170,219,172,221,173,130,129,171,220,173,222,174,131,130,172,221,174,223,175,132,131,173,222,175,224,176,133,132,174,223,176,225,177,134,133,175,224,177,226,227,178,134,176,225,178,227,228,179,135,134,177,179,228,229,180,136,135,178,180,229,230,181,137,136,179,181,230,231,182,138,137,180,182,231,232,183,139,138,181,183,232,233,184,140,139,182,184,233,234,185,141,140,183,185,234,235,236,186,141,184,186,185,236,237,187,142,141,187,186,237,238,188,143,142,188,187,238,239,189,144,143,189,188,239,240,190,145,144,190,189,240,241,191,146,145,191,190,241,242,192,147,146,192,191,242,243,193,148,147,193,192,243,244,245,194,148,194,148,193,245,246,195,149,195,149,194,246,247,196,150,196,150,195,247,248,197,151,197,151,196,248,249,198,152,198,152,197,249,250,199,153,199,153,198,250,251,200,154,200,154,199,251,252,201,155,201,155,200,252,253,254,202,202,156,155,201,254,255,203,203,157,156,202,255,256,204,204,158,157,203,256,257,205,205,159,158,204,257,258,206,206,160,159,205,258,259,207,207,161,160,206,259,260,208,208,162,161,207,260,261,209,209,210,162,208,261,262,263,210,211,163,162,209,263,264,211,212,164,163,210,264,265,212,213,165,164,211,265,266,213,214,166,165,212,266,267,214,215,167,166,213,267,268,215,216,168,167,214,268,269,216,169,127,168,215,269,270,217,272,218,169,270,330,271,218,273,219,170,169,217,272,219,274,220,171,170,218,273,220,275,221,172,171,219,274,221,276,222,173,172,220,275,222,277,223,174,173,221,276,223,278,224,175,174,222,277,224,279,225,176,175,223,278,225,280,226,177,176,224,279,226,281,282,227,177,225,280,227,282,283,228,178,177,226,228,283,284,229,179,178,227,229,284,285,230,180,179,228,230,285,286,231,181,180,229,231,286,287,232,182,181,230,232,287,288,233,183,182,231,233,288,289,234,184,183,232,234,289,290,235,185,184,233,235,290,291,292,236,185,234,236,235,292,293,237,186,185,237,236,293,294,238,187,186,238,237,294,295,239,188,187,239,238,295,296,240,189,188,240,239,296,297,241,190,189,241,240,297,298,242,191,190,242,241,298,299,243,192,191,243,242,299,300,244,193,192,244,243,300,301,302,245,193,245,193,244,302,303,246,194,246,194,245,303,304,247,195,247,195,246,304,305,248,196,248,196,247,305,306,249,197,249,197,248,306,307,250,198,250,198,249,307,308,251,199,251,199,250,308,309,252,200,252,200,251,309,310,253,201,253,201,252,310,311,312,254,254,202,201,253,312,313,255,255,203,202,254,313,314,256,256,204,203,255,314,315,257,257,205,204,256,315,316,258,258,206,205,257,316,317,259,259,207,206,258,317,318,260,260,208,207,259,318,319,261,261,209,208,260,319,320,262,262,263,209,261,320,321,322,263,264,210,209,262,322,323,264,265,211,210,263,323,324,265,266,212,211,264,324,325,266,267,213,212,265,325,326,267,268,214,213,266,326,327,268,269,215,214,267,327,328,269,270,216,215,268,328,329,270,217,169,216,269,329,330,271,272,217,330,273,218,217,271,274,219,218,272,275,220,219,273,276,221,220,274,277,222,221,275,278,223,222,276,279,224,223,277,280,225,224,278,281,226,225,279,281,282,226,280,283,227,226,281,284,228,227,282,285,229,228,283,286,230,229,284,287,231,230,285,288,232,231,286,289,233,232,287,290,234,233,288,291,235,234,289,291,292,235,290,291,293,236,235,292,294,237,236,293,295,238,237,294,296,239,238,295,297,240,239,296,298,241,240,297,299,242,241,298,300,243,242,299,301,244,243,301,300,302,244,244,301,303,245,245,302,304,246,246,303,305,247,247,304,306,248,248,305,307,249,249,306,308,250,250,307,309,251,251,308,310,252,252,309,311,253,311,253,310,312,254,253,311,313,255,254,312,314,256,255,313,315,257,256,314,316,258,257,315,317,259,258,316,318,260,259,317,319,261,260,318,320,262,261,319,321,321,322,262,320,323,263,262,321,324,264,263,322,325,265,264,323,326,266,265,324,327,267,266,325,328,268,267,326,329,269,268,327,330,270,269,328,271,217,270,329]
        connI=[0,7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140,147,154,161,168,175,182,189,196,203,210,217,224,231,238,245,252,259,266,273,280,287,294,301,308,315,322,329,336,343,350,357,364,371,378,385,392,399,406,413,420,427,434,441,448,455,462,469,476,483,490,497,504,511,518,525,532,539,546,553,560,567,574,581,588,595,602,609,616,623,630,637,644,651,658,665,672,679,686,693,700,707,714,721,728,735,742,749,756,763,770,777,784,791,798,805,812,819,826,833,840,847,854,861,868,875,882,889,896,903,910,917,924,931,938,945,952,959,966,973,980,987,994,1001,1008,1015,1022,1029,1036,1043,1050,1057,1064,1071,1078,1085,1092,1099,1106,1113,1120,1127,1134,1141,1148,1155,1162,1169,1176,1183,1190,1197,1204,1211,1218,1225,1232,1239,1246,1253,1260,1267,1274,1281,1288,1295,1302,1309,1316,1323,1330,1337,1344,1351,1358,1365,1372,1379,1386,1393,1400,1407,1414,1421,1428,1435,1442,1449,1456,1463,1470,1477,1484,1491,1498,1505,1512,1519,1526,1533,1540,1547,1554,1561,1568,1575,1582,1589,1596,1603,1610,1617,1624,1631,1638,1645,1652,1659,1666,1673,1680,1687,1694,1701,1708,1715,1722,1729,1736,1743,1750,1757,1764,1771,1778,1785,1792,1799,1806,1813,1820,1827,1834,1841,1848,1855,1862,1869,1876,1883,1890,1897,1901,1905,1909,1913,1917,1921,1925,1929,1933,1937,1941,1945,1949,1953,1957,1961,1965,1969,1973,1977,1981,1985,1989,1993,1997,2001,2005,2009,2013,2017,2021,2025,2029,2033,2037,2041,2045,2049,2053,2057,2061,2065,2069,2073,2077,2081,2085,2089,2093,2097,2101,2105,2109,2113,2117,2121,2125,2129,2133,2137]
        #
        m=MEDCouplingUMesh.New("convexhull",2);
        m.allocateCells(331);
        for i in xrange(331):
            m.insertNextCell(NORM_POLYGON,conn[connI[i]:connI[i+1]]);
            pass
        m.finishInsertingCells();
        coordsDa=DataArrayDouble.New(coords,331,2);
        m.setCoords(coordsDa);
        m.checkCoherency();
        #
        da=m.convexEnvelop2D();
        m.checkCoherency()
        self.assertEqual(coordsDa.getHiddenCppPointer(),m.getCoords().getHiddenCppPointer())
        daC=da.buildComplement(m.getNumberOfCells());
        expected2=DataArrayInt.New([271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,302,303,304,305,306,307,308,309,310,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330]);
        self.assertTrue(expected2.isEqual(daC));
        #
        vals=m.getMeasureField(ON_CELLS).getArray()
        ref=271*[184.69493088478035]+3*[-61.564976961404426,-92.34746544254946,-92.34746544259811,-92.34746544253488,-92.3474654425349,-92.34746544180479,-92.34746544253493,-92.3474654419026,-92.34746544190256,-92.34746544253491]+2*[61.564976961404426,-92.34746544254946,-92.34746544259811,-92.34746544253488,-92.3474654425349,-92.34746544180479,-92.34746544253493,-92.3474654419026,-92.34746544190256,-92.34746544253491]+[-61.564976961404426,-92.34746544254946,-92.34746544259811,-92.34746544253488,-92.3474654425349,-92.34746544180479,-92.34746544253493,-92.3474654419026,-92.34746544190256,-92.34746544253491]
        vals-=DataArrayDouble.New(ref)
        vals.abs()
        theTest=vals.getIdsInRange(-1.,1e-7)
        self.assertTrue(theTest.isIdentity())
        self.assertEqual(331,len(theTest))
        pass

    def testSwigDAIOp8(self):
        da=DataArrayInt.New([7,5,6,7,8,9,9,10,12,13,47,15])
        self.assertTrue(7 in da)
        self.assertTrue(47 in da)
        self.assertTrue(15 in da)
        self.assertEqual(0,da.index(7))
        self.assertEqual(10,da.index(47))
        self.assertTrue(14 not in da)
        self.assertEqual(5,da.search([9,9]))
        self.assertEqual(-1,da.search([5,8]))
        da.rearrange(2)
        self.assertTrue([47,16] not in da)
        self.assertTrue([5,6] not in da)
        self.assertTrue([6,7] in da)
        self.assertEqual(4,da.index([12,13]))
        pass

    def testDataArraySort1(self):
        arr=DataArrayInt.New();
        self.assertRaises(InterpKernelException,arr.sort,True)
        self.assertRaises(InterpKernelException,arr.sort,False)
        values=[2,1,6,5,4,7]
        arr.alloc(3,2);
        self.assertRaises(InterpKernelException,arr.sort,True)
        self.assertRaises(InterpKernelException,arr.sort,False)
        arr.rearrange(1);
        arr.setValues(values,6,1)
        arr1=arr.deepCpy();
        arr2=arr.deepCpy();
        arr1.sort(True);
        expected1=[1,2,4,5,6,7]
        self.assertEqual(6,arr1.getNumberOfTuples());
        self.assertEqual(1,arr1.getNumberOfComponents());
        self.assertEqual(expected1,arr1.getValues());
        arr2.sort(False);
        expected2=[7,6,5,4,2,1]
        self.assertEqual(6,arr2.getNumberOfTuples());
        self.assertEqual(1,arr2.getNumberOfComponents());
        self.assertTrue(expected2,arr2.getValues());
        #
        ard=DataArrayDouble.New();
        self.assertRaises(InterpKernelException,ard.sort,True)
        self.assertRaises(InterpKernelException,ard.sort,False)
        valuesD=[2.,1.,6.,5.,4.,7.]
        ard.alloc(3,2);
        self.assertRaises(InterpKernelException,ard.sort,True)
        self.assertRaises(InterpKernelException,ard.sort,False)
        ard.rearrange(1);
        ard.setValues(valuesD,6,1)
        ard1=ard.deepCpy();
        ard2=ard.deepCpy();
        ard1.sort(True);
        expected3=[1.,2.,4.,5.,6.,7.]
        self.assertEqual(6,ard1.getNumberOfTuples());
        self.assertEqual(1,ard1.getNumberOfComponents());
        for i in xrange(6):
            self.assertAlmostEqual(expected3[i],ard1.getIJ(i,0),12)
            pass
        ard2.sort(False);
        expected4=[7.,6.,5.,4.,2.,1.]
        self.assertEqual(6,ard2.getNumberOfTuples());
        self.assertEqual(1,ard2.getNumberOfComponents());
        for i in xrange(6):
            self.assertAlmostEqual(expected4[i],ard2.getIJ(i,0),12)
            pass
        pass
    
    def testPartitionBySpreadZone1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        m4=MEDCouplingUMesh.MergeUMeshes([m,m[-3:],m[0:2]]);
        m4.renumberCells([5,2,9,6,4,7,0,1,3,8]);
        #
        v2=m4.partitionBySpreadZone();
        self.assertTrue(3,len(v2));
        self.assertTrue(v2[0].isEqual(DataArrayInt.New([0,1,7])))
        self.assertTrue(v2[1].isEqual(DataArrayInt.New([2,4,5,6,9])))
        self.assertTrue(v2[2].isEqual(DataArrayInt.New([3,8])))
        #
        m5=m4.buildSpreadZonesWithPoly();
        self.assertEqual(3,m5.getNumberOfCells());
        self.assertTrue(m5.getCoords().getHiddenCppPointer()==m4.getCoords().getHiddenCppPointer());
        self.assertEqual([5,15,16,17,14,11,13,12,5,2,1,0,3,6,7,8,5,5,18,21,22,20,19],m5.getNodalConnectivity().getValues())
        self.assertEqual([0,8,17,23],m5.getNodalConnectivityIndex().getValues())
        #
        pass

    def testGiveCellsWithType1(self):
        expected0=[1,2]
        expected1=[0,3,4]
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        da=m.giveCellsWithType(NORM_TRI3);
        self.assertEqual(2,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        self.assertEqual(expected0,da.getValues())
        #
        da=m.giveCellsWithType(NORM_QUAD4);
        self.assertEqual(3,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        self.assertEqual(expected1,da.getValues())
        #
        da=m.giveCellsWithType(NORM_TRI6);
        self.assertEqual(0,da.getNumberOfTuples());
        self.assertEqual(1,da.getNumberOfComponents());
        #
        self.assertRaises(InterpKernelException,m.giveCellsWithType,NORM_SEG2)
        self.assertRaises(InterpKernelException,m.giveCellsWithType,NORM_HEXA8)
        pass

    def testSwigDAOp1(self):
        d=DataArrayDouble.New(5,2)
        d.rearrange(1) ; d.iota(2.) ; d.rearrange(2)
        d.setInfoOnComponents(["X [m]","Y [m]"])
        d.setName("AName")
        #
        d1=d+[8,9]
        self.assertTrue(d1.isEqualWithoutConsideringStr(DataArrayDouble.New([10.0,12.0,12.0,14.0,14.0,16.0,16.0,18.0,18.0,20.0]),1e-12))
        d1bis=DataArrayDouble.New([8,9],1,2)+d
        self.assertTrue(d1bis.isEqual(d1,1e-12))
        d1ter=[8,9]+d
        self.assertTrue(d1ter.isEqual(d1,1e-12))
        #
        d2=d1-[8,9]
        self.assertTrue(d2.isEqual(d,1e-12))
        self.assertRaises(InterpKernelException,d1.__rsub__,[8,9])#[8,9]-d1
        #
        d3=d*[8,9]
        self.assertTrue(d3.isEqualWithoutConsideringStr(DataArrayDouble.New([16.0,27.0,32.0,45.0,48.0,63.0,64.0,81.0,80.0,99.0]),1e-12))
        d3bis=DataArrayDouble.New([8,9],1,2)*d
        self.assertTrue(d3bis.isEqual(d3,1e-12))
        d3ter=[8,9]*d
        self.assertTrue(d3ter.isEqual(d3,1e-12))
        #
        d4=d3/[8,9]
        self.assertTrue(d4.isEqual(d,1e-12))
        #
        d=DataArrayInt.New(5,2)
        d.rearrange(1) ; d.iota(2) ; d.rearrange(2)
        d.setInfoOnComponents(["X [m]","Y [m]"])
        d.setName("AName")
        #
        d1=d+[8,9]
        self.assertEqual(d1.getValues(),[10,12,12,14,14,16,16,18,18,20])
        d1bis=DataArrayInt.New([8,9],1,2)+d
        self.assertTrue(d1bis.isEqual(d1))
        d1ter=[8,9]+d
        self.assertTrue(d1ter.isEqual(d1))
        #
        d2=d1-[8,9]
        self.assertTrue(d2.isEqual(d))
        self.assertRaises(InterpKernelException,d1.__rsub__,[8,9])
        #
        d3=d*[8,9]
        self.assertEqual(d3.getValues(),[16,27,32,45,48,63,64,81,80,99])
        d3bis=DataArrayInt.New([8,9],1,2)*d
        self.assertTrue(d3bis.isEqual(d3))
        d3ter=[8,9]*d
        self.assertTrue(d3ter.isEqual(d3))
        #
        d4=d3/[8,9]
        self.assertTrue(d4.isEqual(d))
        #
        d5=d%[4,5]
        self.assertEqual(d5.getValues(),[2,3,0,0,2,2,0,4,2,1])
        pass

    def testSwigSelectTupleId2DAIBug1(self):
        da=DataArrayInt.New([0,1,2,3,12,13,4,5,6,7,14,15,8,9,10,11,16,17])
        self.assertEqual([2,6,10],da[2::6].getValues())
        self.assertEqual([0,4,8],da[::6].getValues())
        self.assertEqual([5,9],da[7::6].getValues())
        self.assertEqual([5],da[7:-5:6].getValues())
        pass

    def testSwigCpp5Safe1(self):
        m=MEDCouplingUMesh.New("toto",2)
        coords=DataArrayDouble.New([0.,0.,1.,0.,1.,1.,0.,1.],4,2)
        m.setCoords(coords)
        vecs=DataArrayDouble.New([2.,3.,4.,5.,6.,7.],3,2)
        expected1=[[2.,3.,3.,3.,3.,4.,2.,4.0],[4.,5.,5.,5.,5.,6.,4.,6.0],[6.,7.,7.,7.,7.,8.,6.,8.0]]
        for pos,vec in enumerate(vecs):
            m2=m.deepCpy()
            m2.translate(vec)
            self.assertTrue(m2.getCoords().isEqual(DataArrayDouble.New(expected1[pos],4,2),1e-12))
            pass
        for pos,vec in enumerate(vecs):
            m2=m.deepCpy()
            m2.translate(vec.buildDADouble())
            self.assertTrue(m2.getCoords().isEqual(DataArrayDouble.New(expected1[pos],4,2),1e-12))
            pass
        pass
    
    def testSwigBugNonRegressionZipDA(self):
        angles=map(lambda x:pi/3*x,xrange(6))
        radius=3
        #
        dad=DataArrayDouble.New(6, 2)
        dad[:,0]=radius
        dad[:,1]=angles
        #
        dad2=dad.fromPolarToCart()
        dads=[dad2.deepCpy() for elt in 7*[None]]
        #
        translationToPerform=[[0.01,0.02],[3./2.*radius,-radius*sqrt(3.)/2],[3./2.*radius,radius*sqrt(3.)/2],[0.,radius*sqrt(3.)],[-3./2.*radius,radius*sqrt(3.)/2],[-3./2.*radius,-radius*sqrt(3.)/2],[0.,-radius*sqrt(3.)]]
        for d,t in zip(dads,translationToPerform):
            d+=t
            pass
        for elt in dads:
            self.assertTrue(not dad2.isEqual(elt,1e-12))
            pass
        for d,t in zip(dads,translationToPerform):
            d-=t
            pass
        for elt in dads:
            self.assertTrue(dad2.isEqual(elt,1e-12))
            pass
        pass

    def testBuildSlice3D2(self):
        mesh3D,mesh2D=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        vec1=[-0.07,1.,0.07]
        origin1=[1.524,1.4552,1.74768]
        slice1,ids=mesh3D.buildSlice3D(origin1,vec1,1e-10);
        f=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
        f.setTime(4.5,6,7) ; f.setMesh(mesh3D)
        arr=DataArrayDouble(mesh3D.getNumberOfCells(),2)
        arr.rearrange(1) ; arr.iota(2.) ; arr.rearrange(2)
        f.setArray(arr)
        f.checkCoherency()
        expected1=DataArrayInt([1,3,4,7,9,10,13,15,16])
        self.assertTrue(expected1.isEqual(ids))
        arr2=arr[expected1]
        #
        f2=f.extractSlice3D(origin1,vec1,1e-10)
        self.assertTrue(f2.getArray().isEqual(arr2,1e-12));
        self.assertTrue(slice1.isEqual(f2.getMesh(),1e-12))
        self.assertEqual(6,f2.getTime()[1]) ; self.assertEqual(7,f2.getTime()[2])
        self.assertAlmostEqual(4.5,f2.getTime()[0],12);
        pass

    def testComputeTupleIdsToSelectFromCellIds1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_3()
        f=MEDCouplingFieldDouble.New(ON_GAUSS_NE,NO_TIME);
        f.setMesh(m);
        arr=DataArrayDouble(52,2) ; arr.rearrange(1) ; arr.iota(7.) ; arr.rearrange(2)
        f.setArray(arr)
        #
        f2=f.buildSubPart([1,5,9])
        f2.checkCoherency()
        cI=m.computeNbOfNodesPerCell()
        cI.computeOffsets2()
        sel=DataArrayInt([1,5,9])
        res=sel.buildExplicitArrByRanges(cI)
        arr2=arr[res]
        self.assertTrue(arr2.isEqual(DataArrayDouble([13,14,15,16,17,18,19,20,59,60,61,62,63,64,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110],15,2),1e-12))
        self.assertTrue(arr2.isEqual(f2.getArray(),1e-12))
        pass

    def testComputeSkin1(self):
        arrX=DataArrayDouble([2.,3.4,5.6,7.7,8.0]) ; arrY=DataArrayDouble([2.,3.4,5.6,7.7,9.0,14.2])
        cmesh=MEDCouplingCMesh() ; cmesh.setCoordsAt(0,arrX) ; cmesh.setCoordsAt(1,arrY)
        umesh=cmesh.buildUnstructured()
        #
        skin=umesh.computeSkin()
        self.assertEqual(18,skin.getNumberOfCells())
        self.assertEqual(1,skin.getMeshDimension())
        self.assertTrue(skin.getCoords().getHiddenCppPointer()==umesh.getCoords().getHiddenCppPointer())
        self.assertEqual([0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54],skin.getNodalConnectivityIndex().getValues())
        self.assertEqual([1,1,0,1,0,5,1,2,1,1,3,2,1,4,3,1,9,4,1,5,10,1,14,9,1,10,15,1,19,14,1,15,20,1,24,19,1,20,25,1,25,26,1,26,27,1,27,28,1,28,29,1,29,24],skin.getNodalConnectivity().getValues())
        ids=skin.computeFetchedNodeIds()
        self.assertEqual([0,1,2,3,4,5,9,10,14,15,19,20,24,25,26,27,28,29],ids.getValues())
        part=umesh.buildFacePartOfMySelfNode(ids,True)
        part.setName(skin.getName());
        self.assertTrue(part.isEqual(skin,1e-12))
        part2=part[1::2]
        part[::2]=part2
        self.assertTrue(not part.isEqual(skin,1e-12))
        trad=part.zipConnectivityTraducer(0)
        self.assertEqual(9,part.getNumberOfCells())
        self.assertEqual([0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8],trad.getValues())
        pass

    def testUMeshSetPartOfMySelf2(self):
        # resize with explicit ids list
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        self.assertEqual([3,4],m.getAllGeoTypes())
        part=m[[0,3,4]]
        part.simplexize(0)
        part2=part[[1,2,5]]
        m[[0,3,4]]=part2
        self.assertEqual([3,0,4,1,3,1,4,2,3,4,5,2,3,6,7,4,3,7,5,4],m.getNodalConnectivity().getValues())
        self.assertEqual([0,4,8,12,16,20],m.getNodalConnectivityIndex().getValues())
        self.assertEqual([3],m.getAllGeoTypes())
        # no resize with explicit ids list
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        part=m[[0,3]]
        part.convertAllToPoly()
        m[[3,4]]=part
        self.assertEqual([4,0,3,4,1,3,1,4,2,3,4,5,2,5,0,3,4,1,5,6,7,4,3],m.getNodalConnectivity().getValues())
        self.assertEqual([0,5,9,13,18,23],m.getNodalConnectivityIndex().getValues())
        self.assertEqual([3,4,5],m.getAllGeoTypes())
        # resize with range ids
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        part=m[3:]
        m[1:3]=part
        self.assertEqual([4,0,3,4,1,4,6,7,4,3,4,7,8,5,4,4,6,7,4,3,4,7,8,5,4],m.getNodalConnectivity().getValues())
        self.assertEqual([0,5,10,15,20,25],m.getNodalConnectivityIndex().getValues())
        self.assertEqual([4],m.getAllGeoTypes())
        # no resize with range ids
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        part=m[0::3]
        part.convertAllToPoly()
        m[3:]=part
        self.assertEqual([4,0,3,4,1,3,1,4,2,3,4,5,2,5,0,3,4,1,5,6,7,4,3],m.getNodalConnectivity().getValues())
        self.assertEqual([0,5,9,13,18,23],m.getNodalConnectivityIndex().getValues())
        self.assertEqual([3,4,5],m.getAllGeoTypes())
        # no resize with range ids negative direction
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        part=m[3::-3]
        part.convertAllToPoly()
        m[:-3:-1]=part
        self.assertEqual([4,0,3,4,1,3,1,4,2,3,4,5,2,5,0,3,4,1,5,6,7,4,3],m.getNodalConnectivity().getValues())
        self.assertEqual([0,5,9,13,18,23],m.getNodalConnectivityIndex().getValues())
        self.assertEqual([3,4,5],m.getAllGeoTypes())
        pass

    def testUnPolyze3(self):
        coord=[0.0,0.5,-0.5,-0.5,-0.5,-0.5,0.5,-0.5,-0.5,0.0,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5,0.5]
        conn=[1,2,5,4,-1,4,3,0,1,-1,2,0,3,5,-1,0,2,1,-1,4,5,3]
        m=MEDCouplingUMesh.New("a mesh",3);
        m.allocateCells(1);
        m.insertNextCell(NORM_POLYHED,22,conn[0:22])
        m.finishInsertingCells();
        coords=DataArrayDouble(coord,6,3);
        m.setCoords(coords);
        m.checkCoherency();
        #
        vol=m.getMeasureField(ON_CELLS);
        self.assertEqual(1,vol.getArray().getNumberOfTuples());
        self.assertAlmostEqual(0.5,vol.getArray().getIJ(0,0),12)
        #
        m.unPolyze();
        #
        self.assertEqual([NORM_PENTA6],m.getAllGeoTypes())
        self.assertTrue(DataArrayInt([0,7]).isEqual(m.getNodalConnectivityIndex()))
        self.assertTrue(DataArrayInt([16,0,2,1,3,5,4]).isEqual(m.getNodalConnectivity()))
        #
        vol=m.getMeasureField(ON_CELLS);
        self.assertEqual(1,vol.getArray().getNumberOfTuples());
        self.assertAlmostEqual(0.5,vol.getArray().getIJ(0,0),12)
        pass

    def testKrSpatialDiscretization1(self):
        srcPointCoordsX=[0.8401877171547095, 0.7830992237586059, 0.9116473579367843, 0.335222755714889, 0.2777747108031878, 0.4773970518621602, 0.3647844727918433, 0.9522297251747128, 0.6357117279599009, 0.1416025553558034]
        srcFieldValsOnPoints=[2.129892434968836, 2.295320474540621, 1.931948594981134, 2.728013590937196, 2.715603240418478, 2.661778472822935, 2.695696990104364, 1.893710234970982, 2.529628016549284, 2.728432341300668]
        targetPointCoordsX=[-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,-6.93889390391e-17,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45]
        targetFieldValsExpected=[2.975379475824351, 2.95613491917003, 2.936890362515361, 2.917645805861018, 2.898401249206574, 2.879156692552137, 2.859912135897732, 2.840667579243201, 2.821423022588731, 2.802178465934342, 2.78293390927989, 2.763689352625457, 2.744444795971001, 2.725209522098197, 2.709077577124666, 2.706677252549218, 2.727467797847971, 2.713338094723676, 2.671342424824244, 2.664877370146978, 2.653840141412181, 2.619607861392791, 2.569777214476479, 2.513263929794591, 2.450732752808528, 2.368313560985155, 2.250909795670307, 2.098194272085416, 1.954257891732065, 1.895040660973802, 1.865256788315972, 1.835475248687992, 1.80569370905998, 1.775912169431971, 1.746130629803976, 1.716349090175918, 1.686567550547855, 1.656786010919941, 1.627004471291988, 1.597222931663817]
        coeffsExpected=DataArrayDouble.New([52.238272642008695, 26.186513281350948, -173.42106377948534, 324.56733663875184, -104.64968873410248, 34.375030568158316, -256.12372208190425, 105.2292032463934, -16.239907618144965, 7.838025836978943, 2.621910745077291, -0.4902609628247241])
        #
        nbOfInputPoints=10;
        f=MEDCouplingFieldDouble.New(ON_NODES_KR,ONE_TIME);
        srcArrX=DataArrayDouble.New(srcPointCoordsX,nbOfInputPoints,1);
        cmesh=MEDCouplingCMesh.New("aMesh");
        cmesh.setCoordsAt(0,srcArrX);
        umesh=cmesh.buildUnstructured();
        f.setMesh(umesh);
        srcVals=DataArrayDouble.New(srcFieldValsOnPoints,nbOfInputPoints,1);
        f.setArray(srcVals);
        f.checkCoherency();
        #
        res0=f.getValueOn(targetPointCoordsX[:1]);
        self.assertAlmostEqual(targetFieldValsExpected[0],res0[0],10)
        #
        valuesToTest=f.getValueOnMulti(targetPointCoordsX);
        self.assertEqual(40,valuesToTest.getNumberOfTuples());
        self.assertEqual(1,valuesToTest.getNumberOfComponents());
        for i in xrange(40):
            self.assertAlmostEqual(targetFieldValsExpected[i],valuesToTest.getIJ(i,0),10)
            pass
        fd=f.getDiscretization()
        del f
        self.assertTrue(isinstance(fd,MEDCouplingFieldDiscretizationKriging))
        coeffs,isDrift=fd.computeVectorOfCoefficients(umesh,srcVals)
        self.assertEqual(2,isDrift)
        self.assertTrue(coeffsExpected.isEqual(coeffs,1e-8))
        #
        pass

    def testDuplicateEachTupleNTimes1(self):
        d=DataArrayDouble.New([9.,8.,7.,6.],4,1) ; d.setInfoOnComponents(["mass [kg]"]) ; d.setName("aname")
        d2=d.duplicateEachTupleNTimes(3)
        self.assertTrue(d2.isEqualWithoutConsideringStr(DataArrayDouble.New([9.,9.,9.,8.,8.,8.,7.,7.,7.,6.,6.,6.],4*3,1),1e-14))
        self.assertEqual("aname",d2.getName())
        self.assertEqual(["mass [kg]"],d2.getInfoOnComponents())
        #
        d=DataArrayInt.New([9,8,7,6],4,1) ; d.setInfoOnComponents(["mass [kg]"]) ; d.setName("aname")
        d2=d.duplicateEachTupleNTimes(3)
        self.assertTrue(d2.isEqualWithoutConsideringStr(DataArrayInt.New([9,9,9,8,8,8,7,7,7,6,6,6],4*3,1)))
        self.assertEqual("aname",d2.getName())
        self.assertEqual(["mass [kg]"],d2.getInfoOnComponents())
        pass

    def testSwigComputeTupleIdsNearTuples1(self):
        da=DataArrayDouble([5.,6.,-5.,-6.,5.,-6.,-5.,6.,5.,6.],5,2)
        arr,arrI=da.computeTupleIdsNearTuples(DataArrayDouble([5.,-6.,5.,6.,-5.,-6.],3,2),1e-10)
        self.assertEqual([2,0,4,1],arr.getValues())
        self.assertEqual([0,1,3,4],arrI.getValues())
        arr,arrI=da.computeTupleIdsNearTuples([5.,-6.,5.,6.,-5.,-6.],1e-10)
        self.assertEqual([2,0,4,1],arr.getValues())
        self.assertEqual([0,1,3,4],arrI.getValues())
        expected0=[[2],[0,4],[1]]
        expected1=[[0,1],[0,2],[0,1]]
        for pos,it in enumerate(DataArrayDouble([5.,-6.,5.,6.,-5.,-6.],3,2)):
            arr,arrI=da.computeTupleIdsNearTuples(it,1e-10)
            self.assertEqual(expected0[pos],arr.getValues())
            self.assertEqual(expected1[pos],arrI.getValues())
            pass
        pass

    def testSwigDataTupleIOp1(self):
        d=DataArrayDouble(10,1)
        d.iota(7.)
        for elt in d:
            elt+=2.
            pass
        toTest=DataArrayDouble([9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0])
        self.assertTrue(toTest.isEqual(d,1e-12))
        for elt in d:
            elt-=2.
            pass
        toTest=DataArrayDouble([7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0])
        self.assertTrue(toTest.isEqual(d,1e-12))
        for elt in d:
            elt*=2.
            pass
        toTest=DataArrayDouble([14.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0])
        self.assertTrue(toTest.isEqual(d,1e-12))
        for elt in d:
            elt/=2.
            pass
        toTest=DataArrayDouble([7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0])
        self.assertTrue(toTest.isEqual(d,1e-12))
        #
        d=DataArrayInt(10,1)
        d.iota(7)
        for elt in d:
            elt+=2
            pass
        self.assertEqual(d.getValues(),[9,10,11,12,13,14,15,16,17,18])
        for elt in d:
            elt-=2
            pass
        self.assertEqual(d.getValues(),[7,8,9,10,11,12,13,14,15,16])
        for elt in d:
            elt*=2
            pass
        self.assertEqual(d.getValues(),[14,16,18,20,22,24,26,28,30,32])
        for elt in d:
            elt/=2
            pass
        self.assertEqual(d.getValues(),[7,8,9,10,11,12,13,14,15,16])
        for elt in d:
            elt%=3
            pass
        self.assertEqual(d.getValues(),[1,2,0,1,2,0,1,2,0,1])
        pass

    def testIntersect2DMeshesTmp5(self):
        coords=DataArrayDouble.New([41,0,42,0,0,42,0,41,41.5,0,29.698484809834998,29.698484809834994,0,41.5,28.991378028648452,28.991378028648445,-42,0,-41,0,-29.698484809834994,29.698484809834998,-41.5,0,-28.991378028648445,28.991378028648452,0,-42,0,-41,-29.698484809835001,-29.698484809834994,0,-41.5,-28.991378028648455,-28.991378028648445,29.698484809834987,-29.698484809835001,28.991378028648441,-28.991378028648455,43,0,0,43,42.5,0,30.405591591021544,30.40559159102154,0,42.5,-43,0,-30.40559159102154,30.405591591021544,-42.5,0,0,-43,-30.405591591021551,-30.40559159102154,0,-42.5,30.405591591021537,-30.405591591021551,44,0,0,44,43.5,0,31.112698372208094,31.112698372208087,0,43.5,-44,0,-31.112698372208087,31.112698372208094,-43.5,0,0,-44,-31.112698372208097,-31.112698372208087,0,-43.5,31.112698372208083,-31.112698372208097,45,0,0,45,44.5,0,31.81980515339464,31.819805153394636,0,44.5,-45,0,-31.819805153394636,31.81980515339464,-44.5,0,0,-45,-31.819805153394647,-31.819805153394636,0,-44.5,31.819805153394629,-31.819805153394647,47,0,0,47,46,0,33.234018715767739,33.234018715767732,0,46,-47,0,-33.234018715767732,33.234018715767739,-46,0,0,-47,-33.234018715767739,-33.234018715767732,0,-46,33.234018715767725,-33.234018715767739,49,0,0,49,48,0,34.648232278140831,34.648232278140824,0,48,-49,0,-34.648232278140824,34.648232278140831,-48,0,0,-49,-34.648232278140839,-34.648232278140824,0,-48,34.648232278140817,-34.648232278140839,51,0,0,51,50,0,36.062445840513924,36.062445840513924,0,50,-51,0,-36.062445840513924,36.062445840513924,-50,0,0,-51,-36.062445840513931,-36.062445840513924,0,-50,36.062445840513917,-36.062445840513931,53,0,0,53,52,0,37.476659402887023,37.476659402887016,0,52,-53,0,-37.476659402887016,37.476659402887023,-52,0,0,-53,-37.47665940288703,-37.476659402887016,0,-52,37.476659402887009,-37.47665940288703,55,0,0,55,54,0,38.890872965260115,38.890872965260108,0,54,-55,0,-38.890872965260108,38.890872965260115,-54,0,0,-55,-38.890872965260122,-38.890872965260108,0,-54,38.890872965260101,-38.890872965260122,59,0,0,59,57,0,41.719300090006307,41.7193000900063,0,57,-59,0,-41.7193000900063,41.719300090006307,-57,0,0,-59,-41.719300090006314,-41.7193000900063,0,-57,41.719300090006293,-41.719300090006314,63,0,0,63,61,0,44.547727214752499,44.547727214752491,0,61,-63,0,-44.547727214752491,44.547727214752499,-61,0,0,-63,-44.547727214752506,-44.547727214752491,0,-61,44.547727214752484,-44.547727214752506,67,0,0,67,65,0,47.37615433949869,47.376154339498683,0,65,-67,0,-47.376154339498683,47.37615433949869,-65,0,0,-67,-47.376154339498697,-47.376154339498683,0,-65,47.376154339498676,-47.376154339498697,71,0,0,71,69,0,50.204581464244875,50.204581464244868,0,69,-71,0,-50.204581464244868,50.204581464244875,-69,0,0,-71,-50.204581464244889,-50.204581464244868,0,-69,50.20458146424486,-50.204581464244889,75,0,0,75,73,0,53.033008588991066,53.033008588991059,0,73,-75,0,-53.033008588991059,53.033008588991066,-73,0,0,-75,-53.033008588991073,-53.033008588991059,0,-73,53.033008588991052,-53.033008588991073,80,0,0,80,77.5,0,56.568542494923804,56.568542494923797,0,77.5,-80,0,-56.568542494923797,56.568542494923804,-77.5,0,0,-80,-56.568542494923818,-56.568542494923797,0,-77.5,56.56854249492379,-56.568542494923818],188,2)
        conn=DataArrayInt.New([8,0,1,2,3,4,5,6,7,8,3,2,8,9,6,10,11,12,8,9,8,13,14,11,15,16,17,8,14,13,1,0,16,18,4,19,8,1,20,21,2,22,23,24,5,8,2,21,25,8,24,26,27,10,8,8,25,28,13,27,29,30,15,8,13,28,20,1,30,31,22,18,8,20,32,33,21,34,35,36,23,8,21,33,37,25,36,38,39,26,8,25,37,40,28,39,41,42,29,8,28,40,32,20,42,43,34,31,8,32,44,45,33,46,47,48,35,8,33,45,49,37,48,50,51,38,8,37,49,52,40,51,53,54,41,8,40,52,44,32,54,55,46,43,8,44,56,57,45,58,59,60,47,8,45,57,61,49,60,62,63,50,8,49,61,64,52,63,65,66,53,8,52,64,56,44,66,67,58,55,8,56,68,69,57,70,71,72,59,8,57,69,73,61,72,74,75,62,8,61,73,76,64,75,77,78,65,8,64,76,68,56,78,79,70,67,8,68,80,81,69,82,83,84,71,8,69,81,85,73,84,86,87,74,8,73,85,88,76,87,89,90,77,8,76,88,80,68,90,91,82,79,8,80,92,93,81,94,95,96,83,8,81,93,97,85,96,98,99,86,8,85,97,100,88,99,101,102,89,8,88,100,92,80,102,103,94,91,8,92,104,105,93,106,107,108,95,8,93,105,109,97,108,110,111,98,8,97,109,112,100,111,113,114,101,8,100,112,104,92,114,115,106,103,8,104,116,117,105,118,119,120,107,8,105,117,121,109,120,122,123,110,8,109,121,124,112,123,125,126,113,8,112,124,116,104,126,127,118,115,8,116,128,129,117,130,131,132,119,8,117,129,133,121,132,134,135,122,8,121,133,136,124,135,137,138,125,8,124,136,128,116,138,139,130,127,8,128,140,141,129,142,143,144,131,8,129,141,145,133,144,146,147,134,8,133,145,148,136,147,149,150,137,8,136,148,140,128,150,151,142,139,8,140,152,153,141,154,155,156,143,8,141,153,157,145,156,158,159,146,8,145,157,160,148,159,161,162,149,8,148,160,152,140,162,163,154,151,8,152,164,165,153,166,167,168,155,8,153,165,169,157,168,170,171,158,8,157,169,172,160,171,173,174,161,8,160,172,164,152,174,175,166,163,8,164,176,177,165,178,179,180,167,8,165,177,181,169,180,182,183,170,8,169,181,184,172,183,185,186,173,8,172,184,176,164,186,187,178,175],540)
        connI=DataArrayInt.New([0,9,18,27,36,45,54,63,72,81,90,99,108,117,126,135,144,153,162,171,180,189,198,207,216,225,234,243,252,261,270,279,288,297,306,315,324,333,342,351,360,369,378,387,396,405,414,423,432,441,450,459,468,477,486,495,504,513,522,531,540],61)
        #
        m1=MEDCouplingUMesh.New("Fix",2);
        m1.setCoords(coords);
        m1.setConnectivity(conn,connI,True);
        #
        coords=DataArrayDouble([46.5,-2.5,53.5,-2.5,53.5,2.5,46.5,2.5,50,-2.5,53.5,0,50,2.5,46.5,0,60.5,-2.5,60.5,2.5,57,-2.5,60.5,0,57,2.5,53.5,7.5,46.5,7.5,53.5,5,50,7.5,46.5,5,60.5,7.5,60.5,5,57,7.5,-2,47,2,47,2,53,-2,53,0,47,2,50,0,53,-2,50,6,47,6,53,4,47,6,50,4,53,2,59,-2,59,2,56,0,59,-2,56,6,59,6,56,4,59],42,2)
        # connectivity
        conn=DataArrayInt([8,0,1,2,3,4,5,6,7,8,1,8,9,2,10,11,12,5,8,3,2,13,14,6,15,16,17,8,2,9,18,13,12,19,20,15,8,21,22,23,24,25,26,27,28,8,22,29,30,23,31,32,33,26,8,24,23,34,35,27,36,37,38,8,23,30,39,34,33,40,41,36],72);
        conn.setName("");
        connI=DataArrayInt([0,9,18,27,36,45,54,63,72],9)
        m2=MEDCouplingUMesh.New("Mobile",2);
        m2.setCoords(coords);
        m2.setConnectivity(conn,connI,True);
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1,m2,1e-10);
        self.assertEqual(105,m3.getNumberOfCells());
        self.assertEqual(105,d1.getNumberOfTuples());
        self.assertEqual(105,d2.getNumberOfTuples());
        self.assertEqual(704,m3.getNumberOfNodes());
        #
        areaExpected=[-65.18804756198824,-65.18804756198824,-65.18804756198824,-65.18804756198824,-66.75884388878285,-66.75884388878285,-66.7588438887833,-66.75884388878308,-68.32964021557768,-68.32964021557768,-68.32964021557814,-68.32964021557791,-69.9004365423732,-69.9004365423732,-69.90043654237297,-69.90043654237297,-1.194568659706448,-1.0869994447159463,-142.2316939607081,-144.51326206513068,-144.5132620651309,-1.1945686597064424,-143.3186934054243,-5.002264310862817,-10.0261332846393,-3.9727823117092953,-7.290862524642649,-124.504404940456,-3.9727823117093237,-146.82366506060032,-150.79644737231024,-5.002264310862776,-145.79418306144626,-5.00208651738126,-10.054764051268958,-4.001067863263231,-8.027932154428669,-129.99378209314813,-4.001067863263216,-153.07856481622616,-157.0796326794898,-5.0020865173811915,-152.07754616210832,-5.001928880064381,-10.050590216368969,-4.00098721602491,-8.025810856794209,-136.28350081741684,-4.000987216024939,-159.36183077064402,-163.36281798667005,-5.0019288800643285,-158.36088910660442,-1.2991516319851801,-3.702636830195414,-3.7815130030068254,-6.265364371195623,-0.02516260900254963,-0.6553944641345026,-3.975752765070567,-7.368528340442765,-142.57249927881398,-0.02516260900254963,-3.9757527650706095,-165.64508791977525,-169.64600329384803,-1.299151631985167,-3.7026368301953885,-164.6442148316677,-10.00321285677458,-20.08414323176165,-8.001644468035863,-16.042954878437143,-304.0096070742277,-8.00164446803587,-350.1399180412005,-358.1415625092368,-10.003212856774468,-348.13834965246224,-3.794150313030109,-8.65049239704272,-0.02260276689354157,-0.5885167811200915,-370.2185414798688,-0.022602766893559393,-383.2517009710623,-383.2743037379555,-3.7941503130300576,-379.48015342492505,-408.40704496667513,-408.4070449666742,-408.4070449666742,-408.4070449666742,-433.53978619538975,-433.5397861953902,-433.5397861953911,-433.53978619539066,-458.67252742410983,-458.6725274241094,-458.67252742410983,-458.6725274241089,-608.6835766330232,-608.6835766330232,-608.6835766330232,-608.6835766330241]
        expected1=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,16,17,18,19,19,20,20,20,20,20,21,21,22,23,23,24,24,24,24,24,25,25,26,27,27,28,28,28,28,28,29,29,30,31,31,32,32,32,32,32,32,32,32,32,33,33,33,34,35,35,35,36,36,36,36,36,37,37,38,39,39,40,40,40,40,40,41,41,42,43,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59]
        expected2=[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,2,-1,-1,-1,0,-1,0,2,4,5,-1,4,-1,-1,0,-1,0,2,4,5,-1,4,-1,-1,0,-1,0,2,4,5,-1,4,-1,-1,0,-1,0,1,2,3,4,5,6,7,-1,4,6,-1,-1,0,1,-1,1,3,6,7,-1,6,-1,-1,1,-1,1,3,6,7,-1,6,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]
        f3=m3.getMeasureField(ON_CELLS).getArray().getValues();
        for i in xrange(105):
            self.assertAlmostEqual(areaExpected[i],f3[i],10)
            pass
        self.assertEqual(expected1,d1.getValues())
        self.assertEqual(expected2,d2.getValues())
        pass

    def testSwig2Intersect2DMeshesQuadra1(self):
        import cmath
        def createDiagCircle(lX, lY, R, cells=[0,1]):  
            """ A circle in a square box, cut along the diagonal. 
            """    
            c = []
            for i in range(8):
              c.append(cmath.rect(R, i*pi/4))
        
            coords = [0.0,0.0,          c[3].real,c[3].imag,       -lX/2.0, lY/2.0,
                      0.0, lY/2.0,      lX/2.0,lY/2.0,             lX/2.0,0.0,
                      #   6                  7                              8
                      lX/2.0,-lY/2.0,   c[7].real,c[7].imag,       c[1].real,c[1].imag,
                      #   9                  10                            11  
                      c[5].real,c[5].imag,   -lX/2.0,-lY/2.0,      0.0, -lY/2.0,
                      #   12                  13                            14
                      -lX/2.0,0.0,         0.0,0.0,                  0.0, 0.0]
            # Points 13 (reps. 14) are average of points (6,7) (resp (1,2))
            coords[13*2]   = 0.5*(coords[6*2]+coords[7*2])
            coords[13*2+1] = 0.5*(coords[6*2+1]+coords[7*2+1])
            coords[14*2]   = 0.5*(coords[1*2]+coords[2*2])
            coords[14*2+1] = 0.5*(coords[1*2+1]+coords[2*2+1])
            connec  = [1,7,8,0]      # half circle up right
            connec3 = [6,7,1,2,4,13,8,14,3,5]
            
            baseMesh = MEDCouplingUMesh.New("box_circle", 2)  
            baseMesh.allocateCells(2)
            meshCoords = DataArrayDouble.New(coords, len(coords)/2, 2)
            meshCoords.setInfoOnComponents(["X [au]", "Y [au]"])
            baseMesh.setCoords(meshCoords)
            
            if 0 in cells:
              baseMesh.insertNextCell(NORM_QPOLYG, connec)  
            if 1 in cells: 
              baseMesh.insertNextCell(NORM_QPOLYG, connec3) 
            baseMesh.finishInsertingCells()  
            baseMesh.checkCoherency() 
            return baseMesh 
        
        eps = 1.0e-7
        m1 = createDiagCircle(1.0, 1.0, 0.5*0.90, cells=[0,1])  
        m2 = createDiagCircle(1.0, 1.0, 0.5*0.95, cells=[0])
        m3, _, _= MEDCouplingUMesh.Intersect2DMeshes(m1, m2, eps)
        m3.mergeNodes(eps)
        m3.convertDegeneratedCells()
        m3.zipCoords()        
        m4 = m3.deepCpy()
        m5, _, _ = MEDCouplingUMesh.Intersect2DMeshes(m3, m4, eps)
        m5.mergeNodes(eps)
        # Check coordinates:
        self.assertTrue(m3.getCoords().isEqual(m5.getCoords(), eps))

    def testIntersect2DMeshesTmp7(self):
        eps = 1.0e-8
        coords = [-0.5,-0.5,   -0.5, 0.5, 0.5, 0.5,    0.5,-0.5]
        connec = range(4)
        m1 = MEDCouplingUMesh.New("box", 2)  
        m1.allocateCells(1)
        meshCoords = DataArrayDouble.New(coords, len(coords)/2, 2)
        m1.setCoords(meshCoords)
        m1.insertNextCell(NORM_POLYGON, connec)
        m1.finishInsertingCells()  
     
        m2 = MEDCouplingDataForTest.buildCircle(0.25, 0.2, 0.4)
        # Was looping indefinitly:
        m_intersec, resToM1, resToM2 = MEDCouplingUMesh.Intersect2DMeshes(m1, m2, eps)
        m_intersec.zipCoords()
        coo_tgt = DataArrayDouble([-0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, -0.5, -0.03284271247461901, 0.4828427124746191, 
          -0.014575131106459124, 0.5000000000000001, 0.5, -0.11224989991991996, 0.24271243444677046, 0.5, 0.5, 0.19387505004004, 
          -0.04799910280454185, -0.06682678787499614, -0.023843325638122054, 0.4915644577163915, 0.5, -0.30612494995996, 0.0, -0.5, 
          -0.5, 0.0, -0.25728756555322957, 0.5, -0.023843325638122026, 0.49156445771639157, -0.04799910280454181, -0.06682678787499613], 17 ,2)
        conn_tgt = [32, 5, 2, 6, 4, 7, 8, 9, 10, 32, 6, 3, 0, 1, 5, 4, 11, 12, 13, 14, 15, 16]
        connI_tgt = [0, 9, 22]
        res1_tgt  = [0, 0]
        res2_tgt = [0, -1]
        self.assert_(coo_tgt.isEqualWithoutConsideringStr(m_intersec.getCoords(), 1e-12))
        self.assertEqual(conn_tgt, m_intersec.getNodalConnectivity().getValues())
        self.assertEqual(connI_tgt, m_intersec.getNodalConnectivityIndex().getValues())
        self.assertEqual(res1_tgt, resToM1.getValues())
        self.assertEqual(res2_tgt, resToM2.getValues())
        
    def testDAIBuildUnique1(self):
        d=DataArrayInt([1,2,2,3,3,3,3,4,5,5,7,7,7,19])
        e=d.buildUnique()
        self.assertTrue(e.isEqual(DataArrayInt([1,2,3,4,5,7,19])))
        pass

    def testDAIPartitionByDifferentValues1(self):
        d=DataArrayInt([1,0,1,2,0,2,2,-3,2])
        expected=[[-3,[7]],[0,[1,4]],[1,[0,2]],[2,[3,5,6,8]]]
        for i,elt in enumerate(zip(*d.partitionByDifferentValues())):
            self.assertEqual(expected[i][0],elt[1])
            self.assertEqual(expected[i][1],elt[0].getValues())
            pass
        pass

    def testFieldGaussMultiDiscPerType1(self):
        coords=DataArrayDouble([0.,0.,0.,1.,1.,1.,1.,0.,0.,0.5,0.5,1.,1.,0.5,0.5,0.],8,2)
        mQ8=MEDCouplingUMesh("",2) ; mQ8.setCoords(coords)
        mQ8.allocateCells(1)
        mQ8.insertNextCell(NORM_QUAD8,range(8))
        mQ8.finishInsertingCells()
        mQ4=MEDCouplingUMesh("",2) ; mQ4.setCoords(coords)
        mQ4.allocateCells(1)
        mQ4.insertNextCell(NORM_QUAD4,range(4))
        mQ4.finishInsertingCells()
        mT3=MEDCouplingUMesh("",2) ; mT3.setCoords(coords)
        mT3.allocateCells(1)
        mT3.insertNextCell(NORM_TRI3,range(3))
        mT3.finishInsertingCells()
        
        tr=[[0.,0.],[2.,0.], [0.,2.],[2.,2.],[4.,2.],[6.,2.],[8.,2.],[10.,2.],[12.,2.],[0.,4.],[2.,4.],[4.,4.],[6.,4.],[8.,4.],[10.,4.],[12.,4.],[14.,4.],[16.,4.],[18.,4.],[20.,4.],[22.,4.]]
        ms=2*[mQ4]+7*[mQ8]+11*[mT3]
        ms[:]=(elt.deepCpy() for elt in ms)
        for m,t in zip(ms,tr):
            d=m.getCoords() ; d+= t
            pass
        m=MEDCouplingUMesh.MergeUMeshes(ms)
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,NO_TIME)
        f.setMesh(m)
        # throw because cell 0,1 are QUAD4 and cell 3 is QUAD8
        self.assertRaises(InterpKernelException,f.setGaussLocalizationOnCells,[0,1,3],[0.,0.,1.,0.,1.,1.,0.,1.],[0.3,0.3,0.7,0.7],[0.8,0.2])
        f.setGaussLocalizationOnCells([0,1],[0.,0.,1.,0.,1.,1.,0.,1.],[0.3,0.3,0.7,0.7],[0.8,0.2])
        f.setGaussLocalizationOnCells([3,2,5],[0.,0.,1.,0.,1.,1.,0.,1.,0.5,0.,1.,0.5,0.5,1.,0.,0.5],[0.3,0.3,0.7,0.7,0.9,0.9],[0.8,0.05,0.15])
        f.setGaussLocalizationOnCells([4,6,8,7],[0.,0.,1.,0.,1.,1.,0.,1.,0.5,0.,1.,0.5,0.5,1.,0.,0.5],[0.3,0.3,0.7,0.7,0.9,0.9,-0.1,0.3],[0.7,0.05,0.15,0.1])
        f.setGaussLocalizationOnCells([9,10,11,12,13],[0.,0.,1.,0.,1.,1.],[0.4,0.4],[1.])
        f.setGaussLocalizationOnCells([14,15,16,17,18,19],[0.,0.,1.,0.,1.,1.],[0.4,0.4,0.14,0.16],[0.22,0.78])
        self.assertEqual(46,f.getNumberOfTuplesExpected())
        vals=DataArrayDouble.New(46*3,1) ; vals.iota(7.7) ; vals.rearrange(3)
        f.setArray(vals)
        f.checkCoherency()
        #f.getLocalizationOfDiscr()
        self.assertRaises(InterpKernelException,f.getGaussLocalizationIdOfOneType,NORM_QUAD8) #throw because several loc
        self.assertEqual([1,2],f.getGaussLocalizationIdsOfOneType(NORM_QUAD8))
        self.assertEqual([0,0,1,1,2,1,2,2,2,3,3,3,3,3,4,4,4,4,4,4],f.getDiscretization().getArrayOfDiscIds().getValues())
        fc=f[[1,2,3,8]]
        fc.checkCoherency()
        self.assertTrue(DataArrayDouble([13.7,14.7,15.7,16.7,17.7,18.7,19.7,20.7,21.7,22.7,23.7,24.7,25.7,26.7,27.7,28.7,29.7,30.7,31.7,32.7,33.7,34.7,35.7,36.7,82.7,83.7,84.7,85.7,86.7,87.7,88.7,89.7,90.7,91.7,92.7,93.7],12,3).isEqual(fc.getArray(),1e-10))
        fc.renumberCells([3,2,0,1])
        self.assertTrue(DataArrayDouble([28.7, 29.7, 30.7, 31.7, 32.7, 33.7, 34.7, 35.7, 36.7, 82.7, 83.7, 84.7, 85.7, 86.7, 87.7, 88.7, 89.7, 90.7, 91.7, 92.7, 93.7, 19.7, 20.7, 21.7, 22.7, 23.7, 24.7, 25.7, 26.7, 27.7, 13.7, 14.7, 15.7, 16.7, 17.7, 18.7],12,3).isEqual(fc.getArray(),1e-10))
        fc.getArray()
        pass

    def testSwigRotate(self):
        d=DataArrayDouble([1.,2.,3.,4.,6.,5.],2,3)
        MEDCouplingPointSet.Rotate3DAlg([0.,0.,0.],[0.,1.,0.],1.5707963267948966,d)
        self.assertTrue(d.isEqual(DataArrayDouble([3.,2.,-1.,5.,6.,-4.],2,3),1e-12))
        d=DataArrayDouble([1.,2.,3.,4.,6.,5.],3,2)
        MEDCouplingPointSet.Rotate2DAlg([0.,0.],1.5707963267948966,d)
        self.assertTrue(d.isEqual(DataArrayDouble([-2.,1.,-4.,3.,-5.,6.],3,2),1e-12))
        pass

    def testSwigCMeshProtection(self):
        cm=MEDCouplingCMesh()
        self.assertRaises(InterpKernelException,cm.setCoordsAt,0,DataArrayDouble([4.,4.5,6.,7.],2,2))
        self.assertRaises(InterpKernelException,cm.setCoords,DataArrayDouble([4.,4.5,6.,7.],2,2))
        pass

    def testSwigCellsInBoundingBox1(self):
        m3D=MEDCouplingDataForTest.build3DExtrudedUMesh_1()[0]
        self.assertTrue(m3D.getCellsInBoundingBox([(0,3),(0,3),(0,1)],-1e-12).isEqual(DataArrayInt([0,1,2,3,4,5])))
        self.assertRaises(InterpKernelException,m3D.getCellsInBoundingBox,[(0,3,0),(3,0,1)],-1e-12)
        pass

    def testDAICheckMonotonic1(self):
        data1=[-1,0,2,2,4,5]
        data2=[6,2,0,-8,-9,-56]
        data3=[-1,0,3,2,4,6]
        data4=[7,5,2,3,0,-6]
        d=DataArrayInt.New(data1);
        self.assertTrue(d.isMonotonic(True));
        self.assertTrue(not d.isMonotonic(False));
        d.checkMonotonic(True);
        self.assertRaises(InterpKernelException,d.checkMonotonic,False)
        d=DataArrayInt.New(data2);
        self.assertTrue(d.isMonotonic(False));
        self.assertTrue(not d.isMonotonic(True));
        d.checkMonotonic(False);
        self.assertRaises(InterpKernelException,d.checkMonotonic,True)
        d=DataArrayInt.New(data3);
        self.assertTrue(not d.isMonotonic(False));
        self.assertTrue(not d.isMonotonic(True));
        self.assertRaises(InterpKernelException,d.checkMonotonic,True)
        self.assertRaises(InterpKernelException,d.checkMonotonic,False)
        d=DataArrayInt.New(data4);
        self.assertTrue(not d.isMonotonic(False));
        self.assertTrue(not d.isMonotonic(True));
        self.assertRaises(InterpKernelException,d.checkMonotonic,True)
        self.assertRaises(InterpKernelException,d.checkMonotonic,False)
        d=DataArrayInt.New(0,1)
        self.assertTrue(d.isMonotonic(True));
        self.assertTrue(d.isMonotonic(False));
        d.checkMonotonic(True);
        d.checkMonotonic(False);
        d=DataArrayInt.New(data4,3,2);#throw because nbComp!=1
        self.assertRaises(InterpKernelException,d.isMonotonic,True)
        self.assertRaises(InterpKernelException,d.isMonotonic,False)
        self.assertRaises(InterpKernelException,d.checkMonotonic,True)
        self.assertRaises(InterpKernelException,d.checkMonotonic,False)
        pass

    def testSwigDASetItemOnEmpty1(self):
        d=DataArrayInt(0,1)
        isThrow=False
        try:
            d[0:1000:2]=4
        except InterpKernelException as e:
            isThrow=True
            pass
        self.assertTrue(isThrow)
        d[:]=4
        d[::2]=5
        #
        d=DataArrayDouble(0,1)
        isThrow=False
        try:
            d[0:1000:2]=4
        except InterpKernelException as e:
            isThrow=True
            pass
        self.assertTrue(isThrow)
        d[:]=4
        d[::2]=5
        d=DataArrayInt([],0,1)
        d2=DataArrayInt(0)
        self.assertTrue(d2.isEqual(d))
        d=DataArrayDouble([],0,1)
        d2=DataArrayDouble(0)
        self.assertTrue(d2.isEqual(d,1e-12))
        pass

    def testSwigDAITransformWithIndArr1(self):
        arr=DataArrayInt([0,4,5,1])
        d=DataArrayInt([7,8,9,10])
        self.assertRaises(InterpKernelException,arr.transformWithIndArr,d)
        pass

    def testIntersect2DMeshesTmp6(self):
        # coordinates
        coords=DataArrayDouble.New([2.7554552980815448e-15,45,-45,5.5109105961630896e-15,-31.819805153394636,31.81980515339464,2.8779199779962799e-15,47,2.8166876380389124e-15,46,-47,5.7558399559925599e-15,-33.234018715767732,33.234018715767739,-46,5.6333752760778247e-15],8,2);
        # connectivity
        conn=DataArrayInt.New([8,0,3,5,1,4,6,7,2])
        connI=DataArrayInt.New([0,9]);
        m1=MEDCouplingUMesh.New("Fixe",2);
        m1.setCoords(coords);
        m1.setConnectivity(conn,connI,True);
        #
        coords=DataArrayDouble.New([-7.3800475508445391,41.854329503018846,-3.7041190667754655,42.338274668899189,-3.7041190667754655,45.338274668899189,-7.3800475508445382,44.854329503018839,-5.5473631693521845,42.136406608386956,-3.7041190667754655,43.838274668899189,-5.5420833088100014,45.09630208595901,-7.3800475508445382,43.354329503018839,-3.7041190667754651,52.338274668899189,-7.3800475508445382,51.854329503018839,-3.7041190667754655,48.838274668899189,-5.5420833088100014,52.09630208595901,-7.3800475508445382,48.354329503018839],13,2);
        # connectivity
        conn=DataArrayInt.New([8,0,1,2,3,4,5,6,7,8,3,2,8,9,6,10,11,12]);
        connI=DataArrayInt.New([0,9,18]);
        #
        m2=MEDCouplingUMesh.New("Mobile",2);
        m2.setCoords(coords);
        m2.setConnectivity(conn,connI,True);
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1,m2,1e-10);
        self.assertTrue(d1.isEqual(DataArrayInt([0,0,0,0])));
        self.assertTrue(d2.isEqual(DataArrayInt([0,1,-1,-1])));
        self.assertEqual(4,m3.getNumberOfCells());
        self.assertEqual(4,d1.getNumberOfTuples());
        self.assertEqual(4,d2.getNumberOfTuples());
        self.assertEqual(43,m3.getNumberOfNodes());
        dI,areMerged,newNbOfNodes=m3.mergeNodes(1e-12)
        self.assertEqual(35,m3.getNumberOfNodes());
        m3.zipCoords();
        self.assertEqual(23,m3.getNumberOfNodes());
        #
        f=m3.getMeasureField(True);
        valuesExpected=DataArrayDouble([1.6603638692585716,5.747555728471923,129.68907101754394,7.4162714498559694])
        self.assertTrue(f.getArray().isEqual(valuesExpected,1e-12))
        pass

    def testDAPushBack(self):
        d=DataArrayDouble(0,1)
        for i in xrange(8):
            d.pushBackSilent(i)
            pass
        self.assertEqual(d.getNumberOfTuples(),8)
        self.assertEqual(d.getNbOfElemAllocated(),8)
        d.pushBackSilent(4.44)
        self.assertEqual(d.getNumberOfTuples(),9)
        self.assertEqual(d.getNbOfElemAllocated(),16)
        self.assertTrue(d.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,4.44]),1e-12))
        e=d.deepCpy()
        self.assertEqual(e.getNumberOfTuples(),9)
        self.assertEqual(e.getNbOfElemAllocated(),9)
        self.assertTrue(e.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,4.44]),1e-12))
        self.assertAlmostEqual(d.popBackSilent(),4.44,12)
        self.assertEqual(d.getNumberOfTuples(),8)
        self.assertEqual(d.getNbOfElemAllocated(),16)
        self.assertTrue(d.isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.]),1e-12))
        f=DataArrayDouble()
        f.reserve(1000)
        f.pushBackSilent(4.)
        self.assertTrue(f.isEqual(DataArrayDouble([4.]),1e-12))
        self.assertEqual(f.getNumberOfTuples(),1)
        self.assertEqual(f.getNbOfElemAllocated(),1000)
        ff=f[:]
        self.assertTrue(ff.isEqual(DataArrayDouble([4.]),1e-12))
        self.assertEqual(ff.getNumberOfTuples(),1)
        self.assertEqual(ff.getNbOfElemAllocated(),1)
        d=DataArrayDouble()
        d.pushBackSilent(4.44)
        d.pushBackSilent(5.55)
        d.pushBackSilent(6.66)
        self.assertTrue(d.isEqual(DataArrayDouble([4.44,5.55,6.66]),1e-12))
        #
        d=DataArrayInt(0,1)
        for i in xrange(8):
            d.pushBackSilent(i)
            pass
        self.assertEqual(d.getNumberOfTuples(),8)
        self.assertEqual(d.getNbOfElemAllocated(),8)
        d.pushBackSilent(444)
        self.assertEqual(d.getNumberOfTuples(),9)
        self.assertEqual(d.getNbOfElemAllocated(),16)
        self.assertTrue(d.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,444])))
        e=d.deepCpy()
        self.assertEqual(e.getNumberOfTuples(),9)
        self.assertEqual(e.getNbOfElemAllocated(),9)
        self.assertTrue(e.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,444])))
        self.assertEqual(d.popBackSilent(),444)
        self.assertEqual(d.getNumberOfTuples(),8)
        self.assertEqual(d.getNbOfElemAllocated(),16)
        self.assertTrue(d.isEqual(DataArrayInt([0,1,2,3,4,5,6,7])))
        f=DataArrayInt()
        f.reserve(1000)
        f.pushBackSilent(4)
        self.assertTrue(f.isEqual(DataArrayInt([4])))
        self.assertEqual(f.getNumberOfTuples(),1)
        self.assertEqual(f.getNbOfElemAllocated(),1000)
        ff=f[:]
        self.assertTrue(ff.isEqual(DataArrayInt([4])))
        self.assertEqual(ff.getNumberOfTuples(),1)
        self.assertEqual(ff.getNbOfElemAllocated(),1)
        d=DataArrayInt()
        d.pushBackSilent(444)
        d.pushBackSilent(555)
        d.pushBackSilent(666)
        self.assertTrue(d.isEqual(DataArrayInt([444,555,666])))
        #
        d=DataArrayInt()
        d.alloc(10,1)
        d.setInfoOnComponent(0,"ABC")
        d.setName("dEf")
        d.iota(7)
        e=DataArrayInt([7,8,9,10,11,12,13,14,15,16]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e))
        self.assertEqual(10,d.getNbOfElemAllocated())
        d.pushBackSilent(55)
        e=DataArrayInt([7,8,9,10,11,12,13,14,15,16,55]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e))
        self.assertEqual(20,d.getNbOfElemAllocated())
        d.reserve(4)
        e=DataArrayInt([7,8,9,10]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e))
        self.assertEqual(4,d.getNbOfElemAllocated())
        d.pushBackSilent(5)
        e=DataArrayInt([7,8,9,10,5]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e))
        self.assertEqual(8,d.getNbOfElemAllocated())
        self.assertEqual(5,d.popBackSilent())
        e=DataArrayInt([7,8,9,10]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e))
        self.assertEqual(8,d.getNbOfElemAllocated())
        self.assertRaises(OverflowError,d.reserve,-1)
        e=DataArrayInt([7,8,9,10]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e))
        self.assertEqual(8,d.getNbOfElemAllocated())
        d.reserve(0)
        e=DataArrayInt([]) ; e.setInfoOnComponent(0,"ABC") ; e.setName("dEf") ; self.assertTrue(d.isEqual(e))
        self.assertEqual(0,d.getNbOfElemAllocated())
        #
        d=DataArrayDouble()
        d.alloc(10,1)
        d.setInfoOnComponent(0,"ABC")
        d.setName("dEf")
        d.iota(7)
        e=DataArrayDouble([7,8,9,10,11,12,13,14,15,16]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(10,d.getNbOfElemAllocated())
        d.pushBackSilent(55)
        e=DataArrayDouble([7,8,9,10,11,12,13,14,15,16,55]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(20,d.getNbOfElemAllocated())
        d.reserve(4)
        e=DataArrayDouble([7,8,9,10]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(4,d.getNbOfElemAllocated())
        d.pushBackSilent(5)
        e=DataArrayDouble([7,8,9,10,5]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(8,d.getNbOfElemAllocated())
        self.assertEqual(5.,d.popBackSilent())
        e=DataArrayDouble([7,8,9,10]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(8,d.getNbOfElemAllocated())
        self.assertRaises(OverflowError,d.reserve,-1)
        e=DataArrayDouble([7,8,9,10]) ; e.copyStringInfoFrom(d) ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(8,d.getNbOfElemAllocated())
        d.reserve(0)
        e=DataArrayDouble([]) ; e.setInfoOnComponent(0,"ABC") ; e.setName("dEf") ; self.assertTrue(d.isEqual(e,1e-14))
        self.assertEqual(0,d.getNbOfElemAllocated())
        pass

    def testDAIBuildSubstractionOptimized1(self):
        da1=DataArrayInt.New([1,3,5,6,7,9,13])
        da2=DataArrayInt.New([3,5,9])
        da3=DataArrayInt.New([1,3,5])
        da4=DataArrayInt.New([1,3,5,6,7,9,13])
        #
        a=da1.buildSubstractionOptimized(da2);
        self.assertTrue(a.isEqual(DataArrayInt([1,6,7,13])));
        #
        a=da1.buildSubstractionOptimized(da3);
        self.assertTrue(a.isEqual(DataArrayInt([6,7,9,13])));
        #
        a=da1.buildSubstractionOptimized(da4);
        self.assertTrue(a.isEqual(DataArrayInt([])));
        pass

    def testDAIIsStrictlyMonotonic1(self):
        da1=DataArrayInt.New([1,3,5,6,7,9,13])
        self.assertTrue(da1.isStrictlyMonotonic(True));
        da1.checkStrictlyMonotonic(True);
        self.assertTrue(da1.isMonotonic(True));
        da1.checkMonotonic(True);
        self.assertTrue(not da1.isStrictlyMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,False)
        self.assertTrue(not da1.isMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,False)
        #
        da1=DataArrayInt.New([1,3,5,6,6,9,13])
        self.assertTrue(not da1.isStrictlyMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,True)
        self.assertTrue(da1.isMonotonic(True));
        da1.checkMonotonic(True);
        self.assertTrue(not da1.isStrictlyMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,False)
        self.assertTrue(not da1.isMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,False)
        #
        da1=DataArrayInt.New([1,3,5,6,5,9,13])
        self.assertTrue(not da1.isStrictlyMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,True)
        self.assertTrue(not da1.isMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,True)
        self.assertTrue(not da1.isStrictlyMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,False)
        self.assertTrue(not da1.isMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,False)
        #
        da1=DataArrayInt.New([13,9,7,6,5,3,1])
        self.assertTrue(not da1.isStrictlyMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,True)
        self.assertTrue(not da1.isMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,True)
        self.assertTrue(da1.isStrictlyMonotonic(False));
        da1.checkStrictlyMonotonic(False);
        self.assertTrue(da1.isMonotonic(False));
        da1.checkMonotonic(False);
        #
        da1=DataArrayInt.New([13,9,6,6,5,3,1])
        self.assertTrue(not da1.isStrictlyMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,True)
        self.assertTrue(not da1.isMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,True)
        self.assertTrue(not da1.isStrictlyMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,False)
        self.assertTrue(da1.isMonotonic(False));
        da1.checkMonotonic(False);
        #
        da1=DataArrayInt.New([13,9,5,6,5,3,1])
        self.assertTrue(not da1.isStrictlyMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,True)
        self.assertTrue(not da1.isMonotonic(True));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,True)
        self.assertTrue(not da1.isStrictlyMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkStrictlyMonotonic,False)
        self.assertTrue(not da1.isMonotonic(False));
        self.assertRaises(InterpKernelException,da1.checkMonotonic,False)
        #
        da1=DataArrayInt.New([])
        self.assertTrue(da1.isStrictlyMonotonic(True));
        da1.checkStrictlyMonotonic(True);
        self.assertTrue(da1.isMonotonic(True));
        da1.checkMonotonic(True);
        self.assertTrue(da1.isStrictlyMonotonic(False));
        da1.checkStrictlyMonotonic(False);
        self.assertTrue(da1.isMonotonic(False));
        da1.checkMonotonic(False);
        #
        da1=DataArrayInt.New([13])
        self.assertTrue(da1.isStrictlyMonotonic(True));
        da1.checkStrictlyMonotonic(True);
        self.assertTrue(da1.isMonotonic(True));
        da1.checkMonotonic(True);
        self.assertTrue(da1.isStrictlyMonotonic(False));
        da1.checkStrictlyMonotonic(False);
        self.assertTrue(da1.isMonotonic(False));
        da1.checkMonotonic(False);
        pass

    def testFindAndCorrectBadOriented3DCells1(self):
        nbOfDisc=20
        vects=([0,0,-1],[0.3,0.7,0.2],[-0.3,0.7,0.2],[-0.3,-0.7,0.2])
        #
        m0=MEDCouplingUMesh("m",3) ; m0.allocateCells(0); m0.insertNextCell(NORM_TETRA4,[0,1,2,3]); #Well oriented
        m1=MEDCouplingUMesh("m",3) ; m1.allocateCells(0); m1.insertNextCell(NORM_PYRA5,[0,1,2,3,4]); #Well oriented
        m2=MEDCouplingUMesh("m",3) ; m2.allocateCells(0); m2.insertNextCell(NORM_PENTA6,[0,1,2,3,4,5]); #Well oriented 
        m3=MEDCouplingUMesh("m",3) ; m3.allocateCells(0); m3.insertNextCell(NORM_HEXA8,[0,1,2,3,4,5,6,7]); #Well oriented
        m4=MEDCouplingUMesh("m",3) ; m4.allocateCells(0)
        self.assertRaises(InterpKernelException,m4.insertNextCell,NORM_HEXGP12,[0,1,2,3,4,5,6,7,8,9,10,11,12]);
        m4.insertNextCell(NORM_HEXGP12,[0,1,2,3,4,5,6,7,8,9,10,11]); #Well oriented
        c0=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,0.,1.],4,3) ; m0.setCoords(c0)
        c1=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,1.,0.,1.,0.,0.,0.,0.,1.],5,3) ; m1.setCoords(c1)
        c2=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,0.,0., 0.,0.,1.,0.,1.,1.,1.,0.,1.],6,3) ; m2.setCoords(c2)
        c3=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,1.,0.,1.,0.,0.,0.,0.,1.,0.,1.,1.,1.,1.,1.,1.,0.,1.],8,3) ; m3.setCoords(c3)
        c4=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,1.,0.,1.,0.,0.,0.8,0.,0.,0.45,0.,0.,   0.,0.,1.,0.,1.,1.,1.,1.,1.,1.,0.,1.,0.8,0.,1.,0.45,0.,1.],12,3) ; m4.setCoords(c4)
        m=MEDCouplingMesh.MergeMeshes([m0,m1,m2,m3,m4])
        expected1=DataArrayDouble([0.16666666666666666,0.3333333333333333,0.5,1.,1.])
        for v in vects:
            for i in xrange(nbOfDisc):
                mm=m.deepCpy()
                mm.rotate([0.,0.,0.],[0.3,0.7,0.2],float(i)/float(nbOfDisc)*2*pi)
                mm2=mm.deepCpy()
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                self.assertTrue(mm.findAndCorrectBadOriented3DCells().empty())
                self.assertTrue(mm.isEqual(mm2,1e-14))
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                mm.convertAllToPoly()
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                pass
            pass
        #
        mOK=m.deepCpy()
        m0=MEDCouplingUMesh("m",3) ; m0.allocateCells(0); m0.insertNextCell(NORM_TETRA4,[0,2,1,3]); #Not well oriented
        m1=MEDCouplingUMesh("m",3) ; m1.allocateCells(0); m1.insertNextCell(NORM_PYRA5,[0,1,2,3,4]); #Well oriented 
        m2=MEDCouplingUMesh("m",3) ; m2.allocateCells(0); m2.insertNextCell(NORM_PENTA6,[0,1,2,3,4,5]); #Well oriented 
        m3=MEDCouplingUMesh("m",3) ; m3.allocateCells(0); m3.insertNextCell(NORM_HEXA8,[0,3,2,1,4,7,6,5]); #Not well oriented
        m4=MEDCouplingUMesh("m",3) ; m4.allocateCells(0); m4.insertNextCell(NORM_HEXGP12,[0,5,4,3,2,1,6,11,10,9,8,7]); #Not well oriented
        m0.setCoords(c0) ; m1.setCoords(c1) ; m2.setCoords(c2) ; m3.setCoords(c3) ; m4.setCoords(c4)
        m=MEDCouplingMesh.MergeMeshes([m0,m1,m2,m3,m4])
        expected2=DataArrayDouble([-0.16666666666666666,0.3333333333333333,0.5,-1.,-1.])
        for v in vects:
            for i in xrange(nbOfDisc):
                mm=m.deepCpy()
                mm.rotate([0.,0.,0.],[0.3,0.7,0.2],float(i)/float(nbOfDisc)*2*pi)
                mm2=mm.deepCpy() ; mm3=mm.deepCpy() ; mm3.convertAllToPoly()
                self.assertTrue(mm3.getMeasureField(False).getArray().isEqual(expected2,1e-14))
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected2,1e-14))
                self.assertTrue(mm.findAndCorrectBadOriented3DCells().isEqual(DataArrayInt([0,3,4])))
                mOK.setCoords(mm.getCoords())
                self.assertTrue(mm.isEqual(mOK,1e-14))
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                mmm=mm.deepCpy()
                self.assertTrue(mmm.findAndCorrectBadOriented3DCells().empty())
                mm.convertAllToPoly()
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                pass
            pass
        #
        m0=MEDCouplingUMesh("m",3) ; m0.allocateCells(0); m0.insertNextCell(NORM_TETRA4,[0,1,2,3]); #Well oriented
        m1=MEDCouplingUMesh("m",3) ; m1.allocateCells(0); m1.insertNextCell(NORM_PYRA5,[0,3,2,1,4]); #Not well oriented 
        m2=MEDCouplingUMesh("m",3) ; m2.allocateCells(0); m2.insertNextCell(NORM_PENTA6,[0,2,1,3,5,4]); #Not well oriented 
        m3=MEDCouplingUMesh("m",3) ; m3.allocateCells(0); m3.insertNextCell(NORM_HEXA8,[0,1,2,3,4,5,6,7]); #Well oriented
        m4=MEDCouplingUMesh("m",3) ; m4.allocateCells(0); m4.insertNextCell(NORM_HEXGP12,range(12)); #Well oriented
        m0.setCoords(c0) ; m1.setCoords(c1) ; m2.setCoords(c2) ; m3.setCoords(c3) ; m4.setCoords(c4)
        m=MEDCouplingMesh.MergeMeshes([m0,m1,m2,m3,m4])
        expected3=DataArrayDouble([0.16666666666666666,-0.3333333333333333,-0.5,1.,1.])
        for v in vects:
            for i in xrange(nbOfDisc):
                mm=m.deepCpy()
                mm.rotate([0.,0.,0.],[0.3,0.7,0.2],float(i)/float(nbOfDisc)*2*pi)
                mm2=mm.deepCpy() ; mm3=mm.deepCpy() ; mm3.convertAllToPoly()
                self.assertTrue(mm3.getMeasureField(False).getArray().isEqual(expected3,1e-14))
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected3,1e-14))
                self.assertTrue(mm.findAndCorrectBadOriented3DCells().isEqual(DataArrayInt([1,2])))
                mOK.setCoords(mm.getCoords())
                self.assertTrue(mm.isEqual(mOK,1e-14))
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                mmm=mm.deepCpy()
                self.assertTrue(mmm.findAndCorrectBadOriented3DCells().empty())
                mm.convertAllToPoly()
                self.assertTrue(mm.getMeasureField(False).getArray().isEqual(expected1,1e-14))
                pass
            pass
        pass

    def testSwig2CellOrientation1(self):
        coords=DataArrayDouble([-0.21606,-0.10803,0.29999999999999999,-0.21606,-0.10803,0.37700000000000006,0,-0.10803,0.29999999999999999,0,-0.10803,0.37700000000000006,0,0.10803,0.29999999999999999,0,0.10803,0.37700000000000006,-0.21606,0.10803,0.29999999999999999,-0.21606,0.10803,0.37700000000000006,0,0.03601,0.29999999999999999,0,0.03601,0.37700000000000006,0,-0.03601,0.29999999999999999,0,-0.03601,0.37700000000000006],12,3)
        conn=[[0,2,10,8,4,6],[1,3,11,9,5,7],[0,1,3,2],[2,3,11,10],[10,11,9,8],[8,9,5,4],[4,5,7,6],[6,7,1,0]]
        for i in xrange(256):
            mesh=MEDCouplingUMesh("FluidMesh_1",3);
            mesh.allocateCells(0)
            conn2=[elt[:] for elt in conn]
            code=bin(i)[2:] ; code='0'*(8-len(code))+code
            for face,rev in zip(conn2,code):
                if bool(int(rev)):
                    face.reverse()
                    pass
                pass
            conn3=[elt+[-1] for elt in conn2]
            conn3=sum(conn3,[])[:-1]
            mesh.insertNextCell(NORM_POLYHED,conn3)
            mesh.setCoords(coords)
            mesh.orientCorrectlyPolyhedrons()
            self.assertTrue(mesh.getBarycenterAndOwner().isEqual(DataArrayDouble([-0.10803,0.,0.3385],1,3),1e-12))
            pass
        pass

    def testSwig2CheckConsecutiveCellTypesForMEDFileFrmt1(self):
        m1=MEDCouplingUMesh("",2) ; m1.allocateCells(0)
        m1.insertNextCell(NORM_QUAD4,[0,1,2,3])
        m1.insertNextCell(NORM_TRI3,[0,1,2])
        d=DataArrayDouble(4,3) ; d[:]=0.
        m1.setCoords(d)
        self.assertTrue(m1.checkConsecutiveCellTypes())
        self.assertTrue(not m1.checkConsecutiveCellTypesForMEDFileFrmt())
        m1.renumberCells([1,0])
        self.assertTrue(m1.checkConsecutiveCellTypes())
        self.assertTrue(m1.checkConsecutiveCellTypesForMEDFileFrmt())
        pass

    def testSwig2DAAccumulate1(self):
        d=DataArrayInt(10) ; d.iota(0)
        self.assertEqual([45],d.accumulate())
        self.assertEqual(45,d.accumulate(0))
        d=DataArrayInt(30) ; d.iota(0) ; d.rearrange(3)
        self.assertEqual([135,145,155],d.accumulate())
        self.assertEqual(135,d.accumulate(0))
        self.assertEqual(145,d.accumulate(1))
        self.assertEqual(155,d.accumulate(2))
        d=DataArrayDouble(10) ; d.iota(0.)
        self.assertEqual([45.],d.accumulate())
        self.assertEqual(45.,d.accumulate(0))
        d=DataArrayDouble(30) ; d.iota(0) ; d.rearrange(3)
        self.assertEqual([135.,145.,155.],d.accumulate())
        self.assertEqual(135.,d.accumulate(0))
        self.assertEqual(145.,d.accumulate(1))
        self.assertEqual(155.,d.accumulate(2))
        pass

    def testSwig2UMeshDistanceToMesh1(self):
        m=MEDCouplingUMesh("toto",2)
        coords=DataArrayDouble([2.3,3.4,5.6,6.5,-4.3,3.2,-9.8,7.6,-5.4],3,3)
        m.setCoords(coords)
        m.allocateCells(0)
        m.insertNextCell(NORM_TRI3,[0,1,2])
        a,b=m.distanceToPoint([-0.335,2.27,1.21])
        self.assertEqual(0,b)
        self.assertAlmostEqual(0.022360988100374124,a,14);
        a,b=m.distanceToPoint(DataArrayDouble([-0.335,2.27,1.21],1,3))
        self.assertEqual(0,b)
        self.assertAlmostEqual(0.022360988100374124,a,14);
        a,b=coords.distanceToTuple([-0.335,2.27,1.21])
        self.assertAlmostEqual(5.243302871282566,a,14)
        self.assertEqual(0,b)
        #
        m=MEDCouplingUMesh("toto",2)
        coords=DataArrayDouble([0.,0.,0., 8.,0.,0., 8.,8.,0., 0.,8.,0.],4,3)
        m.setCoords(coords)
        m.allocateCells(0)
        m.insertNextCell(NORM_QUAD4,[0,1,2,3])
        m.checkCoherency2()
        self.assertEqual([4,0,1,2,3],m.getNodalConnectivity().getValues())
        a,b=m.distanceToPoint([5.,2.,0.1])
        self.assertAlmostEqual(0.1,a,14) ; self.assertEqual(0,b)
        a,b=m.distanceToPoint([5.,-2.,4.])
        self.assertAlmostEqual(sqrt(2*2+4*4),a,14) ; self.assertEqual(0,b)
        m.allocateCells(0)
        m.insertNextCell(NORM_POLYGON,[0,1,2,3])
        m.checkCoherency2()
        self.assertEqual([5,0,1,2,3],m.getNodalConnectivity().getValues())
        a,b=m.distanceToPoint([11.,3.,4.])
        self.assertAlmostEqual(sqrt(3*3+4*4),a,14) ; self.assertEqual(0,b)
        a,b=m.distanceToPoint([4.,12.,5.])
        self.assertAlmostEqual(sqrt(4*4+5*5),a,14) ; self.assertEqual(0,b)
        d=DataArrayDouble([-1.2,3.,2.],1,3)
        for elt in d:
            a,b=m.distanceToPoint(d)
            self.assertAlmostEqual(sqrt(1.2*1.2+2*2),a,14) ; self.assertEqual(0,b)
            pass
        #
        m=MEDCouplingUMesh("toto",1)
        coords=DataArrayDouble([0.,0.,4.,0.,0.,4.],3,2) ; m.setCoords(coords)
        m.allocateCells(0) ; m.insertNextCell(NORM_SEG2,[0,1]) ; m.insertNextCell(NORM_SEG2,[1,2])
        a,b=m.distanceToPoint([-0.1,4.1])
        self.assertAlmostEqual(0.14142135623730925,a,14)  # b==1 self.assertEqual(2,c)
        a,b=m.distanceToPoint([0.,3.9])
        self.assertAlmostEqual(0.07071067811865482,a,14) ; self.assertEqual(1,b) # self.assertEqual(2,c)
        pass

    def testSwig2NonRegressionPartitionBySpreadZone1(self):
        m=MEDCouplingCMesh()
        arr=DataArrayDouble(6) ; arr.iota(0.)
        m.setCoords(arr,arr,arr)
        m=m.buildUnstructured()
        mPart=m[50,80,85,87,92,122]
        zones=mPart.partitionBySpreadZone()
        self.assertEqual(4,len(zones))
        self.assertTrue(zones[0].isEqual(DataArrayInt([0])))
        self.assertTrue(zones[1].isEqual(DataArrayInt([1,2])))
        self.assertTrue(zones[2].isEqual(DataArrayInt([3,4])))
        self.assertTrue(zones[3].isEqual(DataArrayInt([5])))
        #
        n,ni=m.computeNeighborsOfCells()
        a,b=MEDCouplingUMesh.ComputeSpreadZoneGraduallyFromSeed(0,n,ni)
        self.assertEqual(13,b) ; self.assertEqual(125,len(a)) ; self.assertTrue(a.isIdentity())
        a,b=MEDCouplingUMesh.ComputeSpreadZoneGraduallyFromSeed([1],n,ni)
        self.assertEqual(12,b) ; self.assertEqual(125,len(a)) ; self.assertTrue(a.isIdentity())
        a,b=MEDCouplingUMesh.ComputeSpreadZoneGraduallyFromSeed((2,),n,ni)
        self.assertEqual(11,b) ; self.assertEqual(125,len(a)) ; self.assertTrue(a.isIdentity())
        a,b=MEDCouplingUMesh.ComputeSpreadZoneGraduallyFromSeed(DataArrayInt([3]),n,ni)
        self.assertEqual(12,b) ; self.assertEqual(125,len(a)) ; self.assertTrue(a.isIdentity())
        pass

    def testSwigUMeshInsertNextCell1(self):
        m=MEDCouplingUMesh("toto",2)
        #
        coords=DataArrayDouble([0.,0.,1.,1.,1.,0.]) ; m.setCoords(coords)
        da=DataArrayInt([0,1,2])
        m.allocateCells(0)
        for i in xrange(5):
            m.insertNextCell(NORM_TRI3,da)
            pass
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8,12,16,20])))
        #
        da=DataArrayInt([0,1,2,3])
        m.allocateCells(0)
        for i in xrange(5):
            m.insertNextCell(NORM_TRI3,3,da)
            pass
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8,12,16,20])))
        #
        da=DataArrayInt([0,1])
        m.allocateCells(0)
        self.assertRaises(InterpKernelException,m.insertNextCell,NORM_TRI3,3,da)
        #
        da=DataArrayInt([0,1,2,0,1,3,0,1,4,0,1,5,0,1,6],5,3)
        m.allocateCells(0)
        for t in da:
            m.insertNextCell(NORM_TRI3,t)
            pass
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([3,0,1,2,3,0,1,3,3,0,1,4,3,0,1,5,3,0,1,6])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8,12,16,20])))
        self.assertRaises(InterpKernelException,m.insertNextCell,NORM_TRI3,None)
        pass

    def testSwigCurveLinearMesh1(self):
        m=MEDCouplingCurveLinearMesh("toto")
        m.setNodeGridStructure([2,3])
        coords=DataArrayDouble([0.,0., 2.,0., 0.,1., 1.9,1.1, 0.3,1.9, 2.2,2.1],6,2)
        m.setCoords(coords)
        m.checkCoherency()
        m0=m.deepCpy()
        self.assertTrue(m0.isEqual(m,1e-12))
        m.getCoords().setInfoOnComponents(["X [m]","Y [m]"])
        self.assertTrue(not m0.isEqual(m,1e-12))
        m0=m.deepCpy()
        self.assertTrue(m0.isEqual(m,1e-12))
        self.assertEqual(m.getNodeGridStructure(),(2,3))
        pass

    def testSimplexize3(self):
        m=MEDCouplingUMesh("toto",3)
        m.allocateCells(0)
        m.insertNextCell(NORM_TETRA4,[0,1,2,3])
        self.assertEqual([NORM_TETRA4],m.getAllGeoTypesSorted())
        m.insertNextCell(NORM_HEXA8,[4,5,6,7,8,9,10,11])
        self.assertEqual([NORM_TETRA4,NORM_HEXA8],m.getAllGeoTypesSorted())
        m.insertNextCell(NORM_HEXA8,[12,13,14,15,16,17,18,19])
        self.assertEqual([NORM_TETRA4,NORM_HEXA8],m.getAllGeoTypesSorted())
        m.insertNextCell(NORM_TETRA4,[20,21,22,23])
        self.assertEqual([NORM_TETRA4,NORM_HEXA8,NORM_TETRA4],m.getAllGeoTypesSorted())
        c1=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,0.,1.],4,3)
        c2=DataArrayDouble([0.,0.,0.,0.,1.,0.,1.,1.,0.,1.,0.,0., 0.,0.,1.,0.,1.,1.,1.,1.,1.,1.,0.,1.],8,3) ; c2+=[2.,0.,0.]
        c3=c2+[2.,0.,0.]
        c4=c1+[6.,0.,0.]
        c=DataArrayDouble.Aggregate([c1,c2,c3,c4])
        m.setCoords(c)
        m.checkCoherency2()
        #
        m1=m.deepCpy()
        d1=m1.simplexize(PLANAR_FACE_5)
        m1.checkCoherency2()
        vol1=m1.getMeasureField(ON_CELLS).getArray()
        self.assertTrue(vol1.isEqual(DataArrayDouble([1./6, 1./6, 1./6,1./6, 1./6, 1./3,1./6, 1./6, 1./6, 1./6, 1./3, 1./6]),1e-12))
        self.assertEqual(m1.getNodalConnectivity().getValues(),[14,0,1,2,3,14,4,9,5,6,14,4,8,9,11,14,4,7,11,6,14,9,11,10,6,14,4,9,6,11,14,12,17,13,14,14,12,16,17,19,14,12,15,19,14,14,17,19,18,14,14,12,17,14,19,14,20,21,22,23])
        self.assertEqual(m1.getNodalConnectivityIndex().getValues(),[0,5,10,15,20,25,30,35,40,45,50,55,60])
        self.assertTrue(d1.isEqual(DataArrayInt([0,1,1,1,1,1,2,2,2,2,2,3])))
        #
        m2=m.deepCpy()
        d2=m2.simplexize(PLANAR_FACE_6)
        m2.checkCoherency2()
        vol2=m2.getMeasureField(ON_CELLS).getArray()
        self.assertTrue(vol2.isEqual(DataArrayDouble([1./6, 1./6, 1./6,1./6, 1./6, 1./6,1./6,1./6, 1./6, 1./6, 1./6, 1./6,1./6,1./6]),1e-12))
        self.assertEqual(m2.getNodalConnectivity().getValues(),[14,0,1,2,3,14,4,9,5,10,14,4,5,6,10,14,4,8,9,10,14,4,11,8,10,14,4,6,7,10,14,4,7,11,10,14,12,17,13,18,14,12,13,14,18,14,12,16,17,18,14,12,19,16,18,14,12,14,15,18,14,12,15,19,18,14,20,21,22,23])
        self.assertEqual(m2.getNodalConnectivityIndex().getValues(),[0,5,10,15,20,25,30,35,40,45,50,55,60,65,70])
        self.assertTrue(d2.isEqual(DataArrayInt([0,1,1,1,1,1,1,2,2,2,2,2,2,3])))
        pass

    def testSwig2CurveLinearMesh2(self):
        c=MEDCouplingCMesh()
        #2D
        arr1=DataArrayDouble([0,1,3,7])
        arr2=DataArrayDouble([0,1,1.5])
        c.setCoords(arr1,arr2)
        u=c.buildUnstructured()
        coo=u.getCoords()
        cl=MEDCouplingCurveLinearMesh()
        cl.setCoords(coo)
        cl.setNodeGridStructure([4,3])
        cl.checkCoherency2()
        li1=[1.,2.,4.,0.5,1.,2.]
        self.assertTrue(cl.getMeasureField(False).getArray().isEqual(DataArrayDouble(li1),1e-14))
        self.assertTrue(u.getMeasureField(False).getArray().isEqual(DataArrayDouble(li1),1e-14))
        li1_1=[0.5,0.5,2.,0.5,5.,0.5,0.5,1.25,2.,1.25,5.,1.25]
        self.assertTrue(cl.getBarycenterAndOwner().isEqual(DataArrayDouble(li1_1,6,2),1e-14))
        self.assertTrue(u.getBarycenterAndOwner().isEqual(DataArrayDouble(li1_1,6,2),1e-14))
        #3D
        c.setCoords(arr1,arr2,arr2)
        u=c.buildUnstructured()
        coo=u.getCoords()
        cl=MEDCouplingCurveLinearMesh()
        cl.setCoords(coo)
        cl.setNodeGridStructure([4,3,3])
        cl.checkCoherency2()
        li2=[1.,2.,4.,0.5, 1.,2.,0.5,1.,2.,0.25,0.5,1.]
        li2_1=[0.5,0.5,0.5,2.,0.5,0.5,5.,0.5,0.5,0.5,1.25,0.5,2.,1.25,0.5,5.,1.25,0.5,0.5,0.5,1.25,2.,0.5,1.25,5.,0.5,1.25,0.5,1.25,1.25,2.,1.25,1.25,5.,1.25,1.25]
        self.assertTrue(cl.getMeasureField(False).getArray().isEqual(DataArrayDouble(li2),1e-14))
        self.assertTrue(u.getMeasureField(False).getArray().isEqual(DataArrayDouble(li2),1e-14))
        self.assertTrue(cl.getBarycenterAndOwner().isEqual(DataArrayDouble(li2_1,12,3),1e-14))
        self.assertTrue(u.getBarycenterAndOwner().isEqual(DataArrayDouble(li2_1,12,3),1e-14))
        #1D spaceDim 1
        coo=DataArrayDouble(5) ; coo.iota(0.)
        coo=coo*coo
        cl.setCoords(coo)
        cl.setNodeGridStructure([5])
        cl.checkCoherency2()
        li3=[1.,3.,5.,7.]
        li3_1=[0.5,2.5,6.5,12.5]
        self.assertTrue(cl.getMeasureField(False).getArray().isEqual(DataArrayDouble(li3),1e-14))
        self.assertTrue(cl.buildUnstructured().getMeasureField(False).getArray().isEqual(DataArrayDouble(li3),1e-14))
        self.assertTrue(cl.getBarycenterAndOwner().isEqual(DataArrayDouble(li3_1),1e-14))
        self.assertTrue(cl.buildUnstructured().getBarycenterAndOwner().isEqual(DataArrayDouble(li3_1),1e-14))
        #1D spaceDim 2
        coo=DataArrayDouble.Meld(coo,coo)
        cl.setCoords(coo)
        cl.checkCoherency2()
        li4=[sqrt(2.)*elt for elt in [1.,3.,5.,7.]]
        li4_1=[0.5,0.5,2.5,2.5,6.5,6.5,12.5,12.5]
        self.assertEqual(2,cl.getSpaceDimension())
        self.assertEqual(1,cl.getMeshDimension())
        self.assertEqual(4,cl.getNumberOfCells())
        self.assertEqual(5,cl.getNumberOfNodes())
        self.assertTrue(cl.getMeasureField(False).getArray().isEqual(DataArrayDouble(li4),1e-14))
        self.assertTrue(cl.buildUnstructured().getMeasureField(False).getArray().isEqual(DataArrayDouble(li4),1e-14))
        self.assertTrue(cl.getBarycenterAndOwner().isEqual(DataArrayDouble(li4_1,4,2),1e-14))
        self.assertTrue(cl.buildUnstructured().getBarycenterAndOwner().isEqual(DataArrayDouble(li4_1,4,2),1e-14))
        pass

    def testSwig2CurveLinearMeshNonRegression1(self):
        coords=DataArrayDouble([0.0, 0.0, 0.10000000149011612, 0.6000000238418579, 0.10000000149011612, 0.30000001192092896, 1.100000023841858, 0.10000000149011612, 0.20000000298023224, 0.10000000149011612, 0.6000000238418579, 0.20000000298023224, 0.699999988079071, 0.6000000238418579, 0.10000000149011612, 1.2000000476837158, 0.6000000238418579, 0.30000001192092896, 0.10000000149011612, 1.100000023841858, 0.30000001192092896, 0.5, 1.100000023841858, 0.20000000298023224, 1.0, 1.2000000476837158, 0.10000000149011612, 0.0, 0.10000000149011612, 0.5, 0.5, 0.10000000149011612, 0.6000000238418579, 1.2000000476837158, 0.10000000149011612, 0.699999988079071, 0.10000000149011612, 0.6000000238418579, 0.699999988079071, 0.6000000238418579, 0.6000000238418579, 0.5, 1.100000023841858, 0.6000000238418579, 0.6000000238418579, 0.10000000149011612, 1.0, 0.6000000238418579, 0.699999988079071, 1.2000000476837158, 0.699999988079071, 0.8999999761581421, 1.0, 0.5, 0.10000000149011612, 0.10000000149011612, 1.2000000476837158, 0.699999988079071, 0.10000000149011612, 1.0, 1.0, 0.10000000149011612, 1.100000023841858, 0.10000000149011612, 0.6000000238418579, 1.100000023841858, 0.6000000238418579, 0.6000000238418579, 1.100000023841858, 1.100000023841858, 0.6000000238418579, 1.2000000476837158, 0.10000000149011612, 1.2000000476837158, 1.0, 0.5, 1.100000023841858, 1.2000000476837158, 1.2000000476837158, 1.100000023841858, 1.0],27,3)
        m=MEDCouplingCurveLinearMesh("toto")
        m.setCoords(coords)
        m.setNodeGridStructure([3,3,3])
        #
        vol=m.getMeasureField(False).getArray()
        self.assertTrue(vol.isEqual(DataArrayDouble([0.11450000709295281, 0.10583334351579375,0.11149999939029423,0.08866666863113633, 0.1404166805123294,0.1250000135352219,0.1270833433481557,0.13258334288001067]),1e-12))
        self.assertTrue(vol.isEqual(m.buildUnstructured().getMeasureField(False).getArray(),1e-12))
        #
        self.assertTrue(m.getBarycenterAndOwner().isEqual(m.buildUnstructured().getBarycenterAndOwner(),1e-12))
        pass

    def testSwig2NonRegressionDASetSelectedComponents1(self):
        da=DataArrayDouble.New([1.,2.,3.,4.,5.,6.],3,2)
        dv=DataArrayDouble.New();
        dv.alloc(4,4)
        dv.fillWithZero()
        # da has less tuples than dv
        dv.setSelectedComponents(da,[1,0])
        #
        self.assertTrue(dv.isEqual(DataArrayDouble([2.,1.,0.,0.,4.,3.,0.,0.,6.,5.,0.,0.,0.,0.,0.,0.],4,4),1e-14))
        #
        da=DataArrayInt.New([1,2,3,4,5,6],3,2)
        dv=DataArrayInt.New();
        dv.alloc(4,4)
        dv.fillWithZero()
        # da has less tuples than dv
        dv.setSelectedComponents(da,[1,0])
        #
        self.assertTrue(dv.isEqual(DataArrayInt([2,1,0,0,4,3,0,0,6,5,0,0,0,0,0,0],4,4)))
        pass

    def testSwigSetItem3(self):
        # 1-2 
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[3]=[1,2]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,0,0,0,0,1,2,0,0,0,0],6,2),1e-14))
        # 2-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[5,3,2]]=[1,2]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,0,0,1,2,1,2,0,0,1,2],6,2),1e-14))
        # 3-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[:]=[1,2]
        self.assertTrue(d.isEqual(DataArrayDouble([1,2,1,2,1,2,1,2,1,2,1,2],6,2),1e-14))
        # 4-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[DataArrayInt([0,3,4])]=[1,2]
        self.assertTrue(d.isEqual(DataArrayDouble([1,2,0,0,0,0,1,2,1,2,0,0],6,2),1e-14))
        # 5-2
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[5,1]=[7]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,7],6,2),1e-14))
        # 6-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[3,5],1]=[7]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,0,0,0,0,0,7,0,0,0,7],6,2),1e-14))
        # 7-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[:-1:2,1]=[7]
        self.assertTrue(d.isEqual(DataArrayDouble([0,7,0,0,0,7,0,0,0,7,0,0],6,2),1e-14))
        # 8-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[DataArrayInt([0,3,4]),1]=[7]
        self.assertTrue(d.isEqual(DataArrayDouble([0,7,0,0,0,0,0,7,0,7,0,0],6,2),1e-14))
        # 9-2
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[3,[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,0,0,0,0,8,7,0,0,0,0],6,2),1e-14))
        # 10-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[1,3,4],[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,8,7,0,0,8,7,8,7,0,0],6,2),1e-14))
        # 11-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[1::2,[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,8,7,0,0,8,7,0,0,8,7],6,2),1e-14))
        # 12-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[DataArrayInt([1,4]),[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,8,7,0,0,0,0,8,7,0,0],6,2),1e-14))
        # 13-2
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[1,:-1]=[9]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,9,0,0,0,0,0,0,0,0,0],6,2),1e-14))
        # 14-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[1,4,5],:]=[7,8]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,7,8,0,0,0,0,7,8,7,8],6,2),1e-14))
        # 15-2 false
        d=DataArrayDouble([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[1::2,:]=[3,9]
        self.assertTrue(d.isEqual(DataArrayDouble([0,0,3,9,0,0,3,9,0,0,3,9],6,2),1e-14))
        # 1-2 
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[3]=[1,2]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,0,0,1,2,0,0,0,0],6,2)))
        # 2-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[5,3,2]]=[1,2]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,1,2,1,2,0,0,1,2],6,2)))
        # 3-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[:]=[1,2]
        self.assertTrue(d.isEqual(DataArrayInt([1,2,1,2,1,2,1,2,1,2,1,2],6,2)))
        # 4-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[DataArrayInt([0,3,4])]=[1,2]
        self.assertTrue(d.isEqual(DataArrayInt([1,2,0,0,0,0,1,2,1,2,0,0],6,2)))
        # 5-2
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[5,1]=[7]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,7],6,2)))
        # 6-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[3,5],1]=[7]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,0,0,0,7,0,0,0,7],6,2)))
        # 7-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[:-1:2,1]=[7]
        self.assertTrue(d.isEqual(DataArrayInt([0,7,0,0,0,7,0,0,0,7,0,0],6,2)))
        # 8-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[DataArrayInt([0,3,4]),1]=[7]
        self.assertTrue(d.isEqual(DataArrayInt([0,7,0,0,0,0,0,7,0,7,0,0],6,2)))
        # 9-2
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[3,[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,0,0,8,7,0,0,0,0],6,2)))
        # 10-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[1,3,4],[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,8,7,0,0,8,7,8,7,0,0],6,2)))
        # 11-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[1::2,[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,8,7,0,0,8,7,0,0,8,7],6,2)))
        # 12-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[DataArrayInt([1,4]),[1,0]]=[7,8]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,8,7,0,0,0,0,8,7,0,0],6,2)))
        # 13-2
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[1,:-1]=[9]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,9,0,0,0,0,0,0,0,0,0],6,2)))
        # 14-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[[1,4,5],:]=[7,8]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,7,8,0,0,0,0,7,8,7,8],6,2)))
        # 15-2 false
        d=DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0],6,2)
        d[1::2,:]=[3,9]
        self.assertTrue(d.isEqual(DataArrayInt([0,0,3,9,0,0,3,9,0,0,3,9],6,2)))
        pass

    def testSwig2ConvertLinearCellsToQuadratic1(self):
        coordsExp=DataArrayDouble([-0.3,-0.3,0.2,-0.3,0.7,-0.3,-0.3,0.2,0.2,0.2,0.7,0.2,-0.3,0.7,0.2,0.7,0.7,0.7,-0.3,-0.05,-0.05,0.2,0.2,-0.05,-0.05,-0.3,0.45,-0.05,0.45,-0.3,0.45,0.2,0.7,-0.05,-0.05,0.7,0.2,0.45,-0.3,0.45,0.45,0.7,0.7,0.45],22,2)
        # 2D
        m2D=MEDCouplingDataForTest.build2DTargetMesh_1()
        m2D.convertLinearCellsToQuadratic(0)
        m2D.checkCoherency1()
        self.assertEqual(m2D.getNodalConnectivity().getValues(),[8,0,3,4,1,9,10,11,12,6,1,4,2,11,13,14,6,4,5,2,15,16,13,8,6,7,4,3,17,18,10,19,8,7,8,5,4,20,21,15,18])
        self.assertEqual(m2D.getNodalConnectivityIndex().getValues(),[0,9,16,23,32,41])
        self.assertTrue(m2D.getCoords().isEqual(coordsExp,1e-14))
        # 1D
        m1D=MEDCouplingDataForTest.build2DTargetMesh_1().buildDescendingConnectivity()[0]
        m1D.convertLinearCellsToQuadratic(0)
        m1D.checkCoherency1()
        self.assertEqual(m1D.getNodalConnectivity().getValues(),[2,0,3,9,2,3,4,10,2,4,1,11,2,1,0,12,2,4,2,13,2,2,1,14,2,4,5,15,2,5,2,16,2,6,7,17,2,7,4,18,2,3,6,19,2,7,8,20,2,8,5,21])
        self.assertEqual(m1D.getNodalConnectivityIndex().getValues(),[0,4,8,12,16,20,24,28,32,36,40,44,48,52])
        self.assertTrue(m1D.getCoords().isEqual(coordsExp,1e-14))
        # 3D
        m2D=MEDCouplingDataForTest.build2DTargetMesh_1()
        m2D.changeSpaceDimension(3)
        arr=DataArrayDouble(4);  arr.iota(0) ; z=MEDCouplingCMesh() ; z.setCoords(arr)
        m1D=z.buildUnstructured() ; m1D.setCoords(arr.changeNbOfComponents(3,0.))
        m1D.getCoords()[:]=m1D.getCoords()[:,[1,2,0]]
        cooTmp=m2D.getCoords()[:]
        m3D=m2D.buildExtrudedMesh(m1D,0)
        m3D.convertLinearCellsToQuadratic(0)
        m3D.checkCoherency1()
        # check of new m3D content
        coordsExp2=[coordsExp.changeNbOfComponents(3,i) for i in xrange(4)]
        coordsExp3=[DataArrayDouble.Meld(cooTmp[:,[0,1]],cooTmp[:,2]+(0.5+float(i))) for i in xrange(3)]
        coordsExp4=DataArrayDouble.Aggregate([coordsExp2[0],coordsExp3[0],coordsExp2[1],coordsExp3[1],coordsExp2[2],coordsExp3[2],coordsExp2[3]])
        c=DataArrayDouble.Aggregate(m3D.getCoords(),coordsExp4)
        self.assertEqual(len(coordsExp4),115)
        self.assertEqual(len(m3D.getCoords()),115)
        a,b=c.findCommonTuples(1e-14)
        self.assertEqual(len(b),len(coordsExp4)+1)
        e,f=DataArrayInt.BuildOld2NewArrayFromSurjectiveFormat2(2*115,a,b)
        self.assertEqual(f,115)
        self.assertTrue(e.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,0,1,2,3,4,5,6,7,8,36,37,38,39,48,49,53,54,58,59,60,66,67,44,47,52,45,46,57,64,65,70,9,10,11,12,13,14,15,16,17,40,41,42,43,50,51,55,56,61,62,63,68,69,75,78,81,76,77,84,88,89,92,18,19,20,21,22,23,24,25,26,71,72,73,74,79,80,82,83,85,86,87,90,91,97,100,103,98,99,106,110,111,114,27,28,29,30,31,32,33,34,35,93,94,95,96,101,102,104,105,107,108,109,112,113])))
        self.assertTrue(DataArrayInt([30,0,3,4,1,9,12,13,10,36,37,38,39,40,41,42,43,44,45,46,47,25,1,4,2,10,13,11,38,48,49,42,50,51,47,46,52,25,4,5,2,13,14,11,53,54,48,55,56,50,46,57,52,30,6,7,4,3,15,16,13,12,58,59,37,60,61,62,41,63,64,65,46,45,30,7,8,5,4,16,17,14,13,66,67,53,59,68,69,55,62,65,70,57,46,30,9,12,13,10,18,21,22,19,40,41,42,43,71,72,73,74,75,76,77,78,25,10,13,11,19,22,20,42,50,51,73,79,80,78,77,81,25,13,14,11,22,23,20,55,56,50,82,83,79,77,84,81,30,15,16,13,12,24,25,22,21,61,62,41,63,85,86,72,87,88,89,77,76,30,16,17,14,13,25,26,23,22,68,69,55,62,90,91,82,86,89,92,84,77,30,18,21,22,19,27,30,31,28,71,72,73,74,93,94,95,96,97,98,99,100,25,19,22,20,28,31,29,73,79,80,95,101,102,100,99,103,25,22,23,20,31,32,29,82,83,79,104,105,101,99,106,103,30,24,25,22,21,33,34,31,30,85,86,72,87,107,108,94,109,110,111,99,98,30,25,26,23,22,34,35,32,31,90,91,82,86,112,113,104,108,111,114,106,99]).isEqual(m3D.getNodalConnectivity()))
        self.assertTrue(DataArrayInt([0,21,37,53,74,95,116,132,148,169,190,211,227,243,264,285]).isEqual(m3D.getNodalConnectivityIndex()))
        # testing explode3DMeshTo1D
        m3DSlice0=m3D[:5]
        m3DSlice0.zipCoords()
        a,b,c,d,e=m3DSlice0.explode3DMeshTo1D()
        self.assertTrue(b.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,2,12,13,6,14,15,11,10,16,17,18,12,19,20,14,10,21,16,22,23,1,24,25,26,5,27,28,29,10,9,30,31,17,23,32,33,19,26,29,34,21,10])))
        self.assertTrue(c.isEqual(DataArrayInt([0,12,21,30,42,54])))
        self.assertTrue(d.isEqual(DataArrayInt([0,0,3,0,1,0,0,0,3,0,1,0,0,0,3,0,1,2,3,4,0,1,1,2,1,1,2,1,1,2,2,4,2,2,4,2,2,4,3,3,4,3,3,3,4,3,3,3,4,4,4,4,4,4])))
        self.assertTrue(e.isEqual(DataArrayInt([0,1,3,5,6,7,9,11,12,13,15,20,22,24,25,27,28,30,32,33,35,36,38,39,41,42,43,45,46,47,49,50,51,52,53,54])))
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([2,0,3,18,2,3,4,19,2,4,1,20,2,1,0,21,2,9,12,22,2,12,13,23,2,13,10,24,2,10,9,25,2,0,9,26,2,3,12,27,2,4,13,28,2,1,10,29,2,4,2,30,2,2,1,31,2,13,11,32,2,11,10,33,2,2,11,34,2,4,5,35,2,5,2,36,2,13,14,37,2,14,11,38,2,5,14,39,2,6,7,40,2,7,4,41,2,3,6,42,2,15,16,43,2,16,13,44,2,12,15,45,2,6,15,46,2,7,16,47,2,7,8,48,2,8,5,49,2,16,17,50,2,17,14,51,2,8,17,52])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140])))
        self.assertTrue(a.getCoords().isEqual(DataArrayDouble([-0.3,-0.3,0.0,0.2,-0.3,0.0,0.7,-0.3,0.0,-0.3,0.2,0.0,0.2,0.2,0.0,0.7,0.2,0.0,-0.3,0.7,0.0,0.2,0.7,0.0,0.7,0.7,0.0,-0.3,-0.3,1.0,0.2,-0.3,1.0,0.7,-0.3,1.0,-0.3,0.2,1.0,0.2,0.2,1.0,0.7,0.2,1.0,-0.3,0.7,1.0,0.2,0.7,1.0,0.7,0.7,1.0,-0.3,-0.05,0.0,-0.05,0.2,0.0,0.2,-0.05,0.0,-0.05,-0.3,0.0,-0.3,-0.05,1.0,-0.05,0.2,1.0,0.2,-0.05,1.0,-0.05,-0.3,1.0,-0.3,-0.3,0.5,-0.3,0.2,0.5,0.2,0.2,0.5,0.2,-0.3,0.5,0.45,-0.05,0.0,0.45,-0.3,0.0,0.45,-0.05,1.0,0.45,-0.3,1.0,0.7,-0.3,0.5,0.45,0.2,0.0,0.7,-0.05,0.0,0.45,0.2,1.0,0.7,-0.05,1.0,0.7,0.2,0.5,-0.05,0.7,0.0,0.2,0.45,0.0,-0.3,0.45,0.0,-0.05,0.7,1.0,0.2,0.45,1.0,-0.3,0.45,1.0,-0.3,0.7,0.5,0.2,0.7,0.5,0.45,0.7,0.0,0.7,0.45,0.0,0.45,0.7,1.0,0.7,0.45,1.0,0.7,0.7,0.5],53,3),1e-14))
        pass

    def testSwig2DataArrayPushBackValsSilent1(self):
        d=DataArrayDouble()
        d.pushBackValsSilent([4,5,6])
        self.assertTrue(d.isEqual(DataArrayDouble([4.,5.,6.]),1e-14))
        e=DataArrayDouble([1,2,3],1,3)
        for t in e: d.pushBackValsSilent(t)
        self.assertTrue(d.isEqual(DataArrayDouble([4.,5.,6.,1.,2.,3.]),1e-14))
        d.pushBackValsSilent(DataArrayDouble([9,10.]))
        self.assertTrue(d.isEqual(DataArrayDouble([4.,5.,6.,1.,2.,3.,9.,10.]),1e-14))
        d.pushBackValsSilent(DataArrayDouble(0,1))
        self.assertTrue(d.isEqual(DataArrayDouble([4.,5.,6.,1.,2.,3.,9.,10.]),1e-14))
        e=DataArrayDouble([1,2,3],3,1)
        for t in e: d.pushBackValsSilent(t)
        self.assertTrue(d.isEqual(DataArrayDouble([4.,5.,6.,1.,2.,3.,9.,10.,1.,2.,3.]),1e-14))
        d.pushBackValsSilent(77)
        self.assertTrue(d.isEqual(DataArrayDouble([4.,5.,6.,1.,2.,3.,9.,10.,1.,2.,3.,77.]),1e-14))
        #
        d=DataArrayInt()
        d.pushBackValsSilent([4,5,6])
        self.assertTrue(d.isEqual(DataArrayInt([4,5,6])))
        e=DataArrayInt([1,2,3],1,3)
        for t in e: d.pushBackValsSilent(t)
        self.assertTrue(d.isEqual(DataArrayInt([4,5,6,1,2,3])))
        d.pushBackValsSilent(DataArrayInt([9,10]))
        self.assertTrue(d.isEqual(DataArrayInt([4,5,6,1,2,3,9,10])))
        d.pushBackValsSilent(DataArrayInt(0,1))
        self.assertTrue(d.isEqual(DataArrayInt([4,5,6,1,2,3,9,10])))
        e=DataArrayInt([1,2,3],3,1)
        for t in e: d.pushBackValsSilent(t)
        self.assertTrue(d.isEqual(DataArrayInt([4,5,6,1,2,3,9,10,1,2,3])))
        d.pushBackValsSilent(77)
        self.assertTrue(d.isEqual(DataArrayInt([4,5,6,1,2,3,9,10,1,2,3,77])))
        pass

    def testSwig2ConvertLinearCellsToQuadratic2(self):
        m2D=MEDCouplingDataForTest.build2DTargetMesh_1()
        ret=m2D.convertLinearCellsToQuadratic(1)
        self.assertTrue(ret.isIdentity())
        self.assertEqual(5,len(ret))
        m2D.checkCoherency1()
        coordsExp=DataArrayDouble([-0.3,-0.3,0.2,-0.3,0.7,-0.3,-0.3,0.2,0.2,0.2,0.7,0.2,-0.3,0.7,0.2,0.7,0.7,0.7,-0.3,-0.05,-0.05,0.2,0.2,-0.05,-0.05,-0.3,0.45,-0.05,0.45,-0.3,0.45,0.2,0.7,-0.05,-0.05,0.7,0.2,0.45,-0.3,0.45,0.45,0.7,0.7,0.45,-0.05,-0.05,0.3666666666666667,-0.1333333333333333,0.5333333333333332,0.03333333333333334,-0.05,0.45,0.45,0.45],27,2)
        self.assertTrue(m2D.getCoords().isEqual(coordsExp,1e-14))
        self.assertTrue(m2D.getNodalConnectivity().isEqual(DataArrayInt([9,0,3,4,1,9,10,11,12,22,7,1,4,2,11,13,14,23,7,4,5,2,15,16,13,24,9,6,7,4,3,17,18,10,19,25,9,7,8,5,4,20,21,15,18,26])))
        self.assertTrue(m2D.getNodalConnectivityIndex().isEqual(DataArrayInt([0,10,18,26,36,46])))
        #
        m2D=MEDCouplingDataForTest.build2DTargetMesh_1()[(0,3)] ; m2D.zipCoords()
        m2D.changeSpaceDimension(3)
        arr=DataArrayDouble(3);  arr.iota(0) ; z=MEDCouplingCMesh() ; z.setCoords(arr)
        m1D=z.buildUnstructured() ; m1D.setCoords(arr.changeNbOfComponents(3,0.))
        m1D.getCoords()[:]=m1D.getCoords()[:,[1,2,0]]
        cooTmp=m2D.getCoords()[:]
        m3D=m2D.buildExtrudedMesh(m1D,0)
        ret=m3D.convertLinearCellsToQuadratic(1)
        self.assertTrue(ret.isIdentity())
        self.assertEqual(4,len(ret))
        m3D.checkCoherency1()
        coordsExp2=DataArrayDouble([-0.3,-0.3,0.0,0.2,-0.3,0.0,-0.3,0.2,0.0,0.2,0.2,0.0,-0.3,0.7,0.0,0.2,0.7,0.0,-0.3,-0.3,1.0,0.2,-0.3,1.0,-0.3,0.2,1.0,0.2,0.2,1.0,-0.3,0.7,1.0,0.2,0.7,1.0,-0.3,-0.3,2.0,0.2,-0.3,2.0,-0.3,0.2,2.0,0.2,0.2,2.0,-0.3,0.7,2.0,0.2,0.7,2.0,-0.3,-0.05,0.0,-0.05,0.2,0.0,0.2,-0.05,0.0,-0.05,-0.3,0.0,-0.3,-0.05,1.0,-0.05,0.2,1.0,0.2,-0.05,1.0,-0.05,-0.3,1.0,-0.3,-0.3,0.5,-0.3,0.2,0.5,0.2,0.2,0.5,0.2,-0.3,0.5,-0.05,0.7,0.0,0.2,0.45,0.0,-0.3,0.45,0.0,-0.05,0.7,1.0,0.2,0.45,1.0,-0.3,0.45,1.0,-0.3,0.7,0.5,0.2,0.7,0.5,-0.3,-0.05,2.0,-0.05,0.2,2.0,0.2,-0.05,2.0,-0.05,-0.3,2.0,-0.3,-0.3,1.5,-0.3,0.2,1.5,0.2,0.2,1.5,0.2,-0.3,1.5,-0.05,0.7,2.0,0.2,0.45,2.0,-0.3,0.45,2.0,-0.3,0.7,1.5,0.2,0.7,1.5,-0.05,-0.05,0.0,-0.3,-0.05,0.5,-0.05,0.2,0.5,0.2,-0.05,0.5,-0.05,-0.3,0.5,-0.05,-0.05,1.0,-0.05,0.45,0.0,-0.05,0.7,0.5,0.2,0.45,0.5,-0.3,0.45,0.5,-0.05,0.45,1.0,-0.3,-0.05,1.5,-0.05,0.2,1.5,0.2,-0.05,1.5,-0.05,-0.3,1.5,-0.05,-0.05,2.0,-0.05,0.7,1.5,0.2,0.45,1.5,-0.3,0.45,1.5,-0.05,0.45,2.0,-0.05,-0.05,0.5,-0.05,0.45,0.5,-0.05,-0.05,1.5,-0.05,0.45,1.5],75,3)
        self.assertTrue(m3D.getCoords().isEqual(coordsExp2,1e-14))
        self.assertTrue(m3D.getNodalConnectivity().isEqual(DataArrayInt([27,0,2,3,1,6,8,9,7,18,19,20,21,22,23,24,25,26,27,28,29,51,52,53,54,55,56,71,27,4,5,3,2,10,11,9,8,30,31,19,32,33,34,23,35,36,37,28,27,57,58,59,53,60,61,72,27,6,8,9,7,12,14,15,13,22,23,24,25,38,39,40,41,42,43,44,45,56,62,63,64,65,66,73,27,10,11,9,8,16,17,15,14,33,34,23,35,46,47,39,48,49,50,44,43,61,67,68,63,69,70,74])))
        self.assertTrue(m3D.getNodalConnectivityIndex().isEqual(DataArrayInt([0,28,56,84,112])))
        pass

    def testSwig2GaussNEIntegral1(self):
        m2D=MEDCouplingDataForTest.build2DTargetMesh_1()
        m0=m2D[0] ; m0.zipCoords()
        m1=m2D[[1,2]] ; m1.zipCoords()
        m2=m2D[[3,4]] ; m2.zipCoords()
        m0.convertLinearCellsToQuadratic(1)
        m1.convertLinearCellsToQuadratic(0)
        m2.convertLinearCellsToQuadratic(1)
        m=MEDCouplingUMesh.MergeUMeshes([m0,m1,m2])
        m.mergeNodes(1e-12)
        f=MEDCouplingFieldDouble(ON_GAUSS_NE)
        f.setMesh(m)
        arr=DataArrayDouble([1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9,
                             11.1,12.2,13.3,14.4,15.5,16.6,
                             21.1,22.2,23.3,24.4,25.5,26.6,
                             31.1,32.2,33.3,34.4,35.5,36.6,37.7,38.8,39.9,
                             41.1,42.2,43.3,44.4,45.5,46.6,47.7,48.8,49.9])
        arr2=DataArrayDouble(len(arr),2)
        arr2[:,0]=arr ; arr2[:,1]=arr+100
        f.setArray(arr2)
        f.checkCoherency()
        res=f.integral(False)
        # a=25./81 ; b=40./81 ; c=64./81
        # p1=0.11169079483905 ; p2=0.0549758718227661
        # 1st compo
        # c0=(a*(1.1+2.2+3.3+4.4)+b*(5.5+6.6+7.7+8.8)+c*9.9)*0.25/3.9999999999999978 ; c0=1.5837962962962973
        # c1=(p2*(11.1+12.2+13.3)+p1*(14.4+15.5+16.6))*0.125/0.4999999999854482 ; c1=1.8014347172346943
        # c2=(p2*(21.1+22.2+23.3)+p1*(24.4+25.5+26.6))*0.125/0.4999999999854482 ; c2=3.0514347172346943
        # c3=(a*(31.1+32.2+33.3+34.4)+b*(35.5+36.6+37.7+38.8)+c*39.9)*0.25/3.9999999999999978 ; c3=9.0837962962963
        # c4=(a*(41.1+42.2+43.3+44.4)+b*(45.5+46.6+47.7+48.8)+c*49.9)*0.25/3.9999999999999978 ; c4=11.583796296296303
        # c0+c1+c2+c3+c4=27.104258323358287
        integExp0=27.104258323358287
        self.assertAlmostEqual(res[0],integExp0,13)
        # 2nd compo
        # c0=(a*(101.1+102.2+103.3+104.4)+b*(105.5+106.6+107.7+108.8)+c*109.9)*0.25/3.9999999999999978 ; c0=26.58379629629631
        # c1=(p2*(111.1+112.2+113.3)+p1*(114.4+115.5+116.6))*0.125/0.4999999999854482 ; c1=14.301434717234699
        # c2=(p2*(121.1+122.2+123.3)+p1*(124.4+125.5+126.6))*0.125/0.4999999999854482 ; c2=15.5514347172347
        # c3=(a*(131.1+132.2+133.3+134.4)+b*(135.5+136.6+137.7+138.8)+c*139.9)*0.25/3.9999999999999978 ; c3=34.08379629629631
        # c4=(a*(141.1+142.2+143.3+144.4)+b*(145.5+146.6+147.7+148.8)+c*149.9)*0.25/3.9999999999999978 ; c4=36.58379629629632
        # c0+c1+c2+c3+c4=127.10425832335835
        integExp1=127.10425832335835
        self.assertAlmostEqual(res[1],integExp1,12)
        meas=f.getDiscretization().getMeasureField(f.getMesh(),False)
        intPerTuple=meas*f
        res2=intPerTuple.accumulate()
        self.assertAlmostEqual(res2[0],integExp0,13)
        self.assertAlmostEqual(res2[1],integExp1,12)
        #
        meas2=f.buildMeasureField(False)
        intPerTuple=meas2*f
        res3=intPerTuple.accumulate()
        self.assertAlmostEqual(res3[0],integExp0,13)
        self.assertAlmostEqual(res3[1],integExp1,12)
        #
        res4=f.getWeightedAverageValue(False) # res4==res2 because sum of area of mesh is equal to 1
        self.assertAlmostEqual(res4[0],integExp0,13)
        self.assertAlmostEqual(res4[1],integExp1,12)
        #
        m.scale([0,0],2.)
        #
        res5=f.getWeightedAverageValue() # res4==res4 because weighted average is not sensitive to the scaling
        self.assertAlmostEqual(res5[0],integExp0,13)
        self.assertAlmostEqual(res5[1],integExp1,12)
        meas3=f.buildMeasureField(False)
        delta=4*meas2.getArray()-meas3.getArray()
        delta.abs()
        self.assertTrue(delta.isUniform(0.,1e-16))
        res6=f.integral(False)
        self.assertAlmostEqual(res6[0],4.*integExp0,12)
        self.assertAlmostEqual(res6[1],4.*integExp1,11)
        pass

    def testSwig2SlowDADFindClosestTupleId(self):
        nbPts=[10,]
        for nbPt in nbPts:
            d=DataArrayDouble(nbPt) ; d.iota() ; d*=1./(nbPt-1)
            c=MEDCouplingCMesh() ; c.setCoords(d,d) ; m=c.buildUnstructured() ; pts=m.getCoords() ; del m
            #
            d0=DataArrayDouble((nbPt-1)*(nbPt-1)) ; d0.iota() ; d0*=(3./((nbPt-1)*(nbPt-1))) ; d0=d0.applyFunc("exp(x)-1")
            d1=DataArrayDouble((nbPt-1)*(nbPt-1)) ; d1.iota()
            d2=DataArrayDouble.Meld(d0,d1) ; d2=d2.fromPolarToCart() ; d2+=[0.32,0.73]
            ids=pts.findClosestTupleId(d2)
            #print "Start of costly computation"
            idsExpected=DataArrayInt(len(d2))
            tmp=1e300
            for i,elt in enumerate(d2):
                l,m=(pts-elt).magnitude().getMinValue()
                idsExpected.setIJSilent(i,0,m)
                if l<tmp:
                    tmp=l ; tmp1=m ; tmp2=i
                    pass
                pass
            #print "End of costly computation"
            self.assertTrue(idsExpected.isEqual(ids))
            a,b,c=pts.minimalDistanceTo(d2)
            self.assertEqual(tmp,a)
            self.assertEqual(tmp1,b)
            self.assertEqual(tmp2,c)
            #
            l=[d2[:,i] for i in [0,1]]
            for elt in l: elt.reverse()
            d2i=DataArrayDouble.Meld(l)
            ids1=pts.findClosestTupleId(d2i)
            idsExpectedI=idsExpected.deepCpy() ; idsExpectedI.reverse()
            self.assertTrue(idsExpectedI.isEqual(ids1))
            #
            l=[pts[:,i] for i in [0,1]]
            for elt in l: elt.reverse()
            ptsi=DataArrayDouble.Meld(l)
            ids2=ptsi.findClosestTupleId(d2)
            idsExpected2=nbPt*nbPt-1-ids
            self.assertTrue(idsExpected2.isEqual(ids2))
            #
            ids3=ptsi.findClosestTupleId(d2i)
            idsExpected3=idsExpected2.deepCpy() ; idsExpected3.reverse()
            self.assertTrue(idsExpected3.isEqual(ids3))
            pass

    def testSwig2DataArrayAsciiChar1(self):
        alpha=DataArrayInt(26) ; alpha.iota(ord("A"))
        d=DataArrayAsciiChar(alpha.getValues(),2,13)
        d.setInfoOnComponents(["c%i"%(v) for v in xrange(13)])
        self.assertEqual('ABCDEFGHIJKLM',d.getTuple(0))
        self.assertEqual('NOPQRSTUVWXYZ',d.getTuple(1))
        self.assertEqual(2,d.getNumberOfTuples())
        self.assertEqual(26,d.getNbOfElems())
        self.assertEqual(13,d.getNumberOfComponents())
        dd=d.deepCpy()
        self.assertTrue(d.isEqual(dd))
        dd.setIJ(0,3,'d')
        self.assertTrue(not d.isEqual(dd))
        d.setIJ(0,3,ord('d'))
        self.assertTrue(d.isEqual(dd))
        d.rearrange(1)
        d.reserve(20)
        self.assertEqual(20,d.getNumberOfTuples())
        self.assertEqual(20,d.getNbOfElems())
        self.assertEqual(1,d.getNumberOfComponents())
        #
        d0=DataArrayAsciiChar([ord('a')],1,1)
        self.assertEqual('a',d0.asciiCharValue())
        self.assertTrue(not d0.empty())
        d0=DataArrayAsciiChar(0,3)
        self.assertTrue(d0.empty())
        d.pushBackSilent("U") ; d.pushBackSilent("V") ; d.pushBackSilent("W")
        self.assertEqual("W",d.popBackSilent())
        d.rearrange(2)
        self.assertEqual(['AB','Cd','EF','GH','IJ','KL','MN','OP','QR','ST','UV'],d.toStrList())
        d.fillWithZero()
        self.assertEqual(11*[''],d.toStrList())
        d.fillWithValue('T')
        self.assertEqual(11*["TT"],d.toStrList())
        d.rearrange(1)
        self.assertTrue(d.isUniform("T"))
        d.rearrange(2)
        #
        dd.rearrange(2)
        dd2=dd.deepCpy()
        dd.renumberInPlace([3,1,2,4,0,11,10,9,8,7,5,12,6])
        self.assertEqual(dd.toStrList(),['IJ','Cd','EF','AB','GH','UV','YZ','ST','QR','OP','MN','KL','WX'])
        dd.renumberInPlaceR([3,1,2,4,0,11,10,9,8,7,5,12,6])
        self.assertEqual(['AB','Cd','EF','GH','IJ','KL','MN','OP','QR','ST','UV','WX','YZ'],dd.toStrList())
        e=dd.renumber([3,1,2,4,0,11,10,9,8,7,5,12,6])
        self.assertEqual(e.toStrList(),['IJ','Cd','EF','AB','GH','UV','YZ','ST','QR','OP','MN','KL','WX'])
        e=dd.renumberR([3,1,2,4,0,11,10,9,8,7,5,12,6])
        self.assertEqual(e.toStrList(),['GH','Cd','EF','IJ','AB','WX','UV','ST','QR','OP','KL','YZ','MN'])
        e=dd.renumberAndReduce([1,1,1,1,1,1,1,2,0,0,0,0,0],3)
        self.assertEqual(['YZ','MN','OP'],e.toStrList())
        self.assertEqual(['GH','IJ'],dd.selectByTupleIdSafe([3,4]).toStrList())
        self.assertEqual(['AB','GH','MN','ST','YZ'],dd.selectByTupleId2(0,13,3).toStrList())
        dd3=dd.changeNbOfComponents(3,"G")
        self.assertEqual(['ABG','CdG','EFG','GHG','IJG','KLG','MNG','OPG','QRG','STG','UVG','WXG','YZG'],dd3.toStrList())
        dd3.rearrange(1) ; self.assertEqual("G",dd3.back()) ; dd3.rearrange(3)
        self.assertTrue(dd3.changeNbOfComponents(2,"\0").isEqual(dd))
        self.assertEqual(len(dd),13)
        d=DataArrayAsciiChar(13,2) ; d.fillWithValue('Y')
        dd3.meldWith(d)
        self.assertEqual(['ABGYY','CdGYY','EFGYY','GHGYY','IJGYY','KLGYY','MNGYY','OPGYY','QRGYY','STGYY','UVGYY','WXGYY','YZGYY'],dd3.toStrList())
        self.assertEqual("d",dd3.getIJ(0,6))
        self.assertRaises(InterpKernelException,dd3.getIJSafe,0,6)
        self.assertEqual("d",dd3.getIJSafe(1,1))
        dd3.rearrange(1)
        e=dd3.getIdsEqual("Y")
        self.assertTrue(e.isEqual(DataArrayInt([3,4,8,9,13,14,18,19,23,24,28,29,33,34,38,39,43,44,48,49,53,54,58,59,60,63,64])))
        e=dd3.getIdsNotEqual("Y")
        self.assertTrue(e.isEqual(DataArrayInt([0,1,2,5,6,7,10,11,12,15,16,17,20,21,22,25,26,27,30,31,32,35,36,37,40,41,42,45,46,47,50,51,52,55,56,57,61,62])))
        self.assertEqual(("d",6),dd3.getMaxValue())
        self.assertEqual(("A",0),dd3.getMinValue())
        self.assertEqual(26,dd3.search("LGYYM"))
        self.assertEqual(-1,dd3.search("LGYYN"))
        dd3.rearrange(5)
        self.assertEqual(7,dd3.locateTuple("OPGYY"))
        self.assertTrue("OPGYY" in dd3)
        self.assertEqual(7,dd3.index("OPGYY"))
        self.assertEqual(-1,dd3.locateTuple("OPGYP"))
        dd3.rearrange(1)
        self.assertEqual(2,dd3.locateValue("OPGYY"))
        self.assertTrue(dd3.presenceOfValue("OPGYY"))
        self.assertTrue("O" in dd3)
        self.assertTrue(not dd3.presenceOfValue("z"))
        self.assertTrue("z" not in dd3)
        dd3.rearrange(5)
        l=list(dd3)
        self.assertEqual([e.buildDAAsciiChar().toStrList()[0] for e in list(dd3)],dd3.toStrList())
        dd3.reAlloc(5)
        dd4=DataArrayChar.Aggregate(dd3,dd3)
        self.assertEqual(['ABGYY','CdGYY','EFGYY','GHGYY','IJGYY','ABGYY','CdGYY','EFGYY','GHGYY','IJGYY'],dd4.toStrList())
        dd5=DataArrayChar.Aggregate([dd4,dd3,dd4])
        self.assertEqual(['ABGYY','CdGYY','EFGYY','GHGYY','IJGYY','ABGYY','CdGYY','EFGYY','GHGYY','IJGYY','ABGYY','CdGYY','EFGYY','GHGYY','IJGYY','ABGYY','CdGYY','EFGYY','GHGYY','IJGYY','ABGYY','CdGYY','EFGYY','GHGYY','IJGYY'],dd5.toStrList())
        # getitem,__iter__,__setitem__
        a=list(dd3)
        self.assertEqual("ABGYY",str(a[0]))
        dd4=dd3[::2]
        self.assertEqual(['ABGYY','EFGYY','IJGYY'],dd4.toStrList())
        dd4=dd3[(3,2,1)]
        self.assertEqual(['GHGYY','EFGYY','CdGYY'],dd4.toStrList())
        dd4=dd3[:]
        dd4[::2]=["12","345","67890"]
        self.assertEqual(['12   ','CdGYY','345  ','GHGYY','67890'],dd4.toStrList())
        dd4=dd3[:]
        dd4[[1,2]]=" "
        self.assertEqual(['ABGYY','     ','     ','GHGYY','IJGYY'],dd4.toStrList())
        dd4=dd3[:]
        dd4[4]='12345'
        self.assertEqual(['ABGYY','CdGYY','EFGYY','GHGYY','12345'],dd4.toStrList())
        dd4[0]=dd4[1]
        self.assertEqual(['CdGYY','CdGYY','EFGYY','GHGYY','12345'],dd4.toStrList())
        dd4=DataArrayAsciiChar(["abc","de","fghi"])
        self.assertEqual(['abc ','de  ','fghi'],dd4.toStrList())
        dd4=DataArrayAsciiChar(["abc","de","fghi"],"t")
        self.assertEqual(['abct','dett','fghi'],dd4.toStrList())
        pass

    def testSwig2GaussNELocalizationOfDiscValues(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_3()[[0,1,3,4,5,6,8,9]] # suppression of polygons
        f=MEDCouplingFieldDouble(ON_GAUSS_NE)
        f.setMesh(m)
        loc=f.getLocalizationOfDiscr()
        self.assertEqual(42,len(loc))
        self.assertTrue(loc.isEqual(DataArrayDouble([0.,0.,1.,0.,0.5,1.,0.,0.,1.,0.,1.,1.,0.,1.,0.,0.,1.,0.,0.5,1.,0.5,0.,0.75,0.5,0.25,0.5,0.,0.,1.,0.,1.,1.,0.,1.,0.5,0.,1.,0.5,0.5,1.,0.,0.5,0.,0.,0.5,1.,1.,0.,0.,0.,0.,1.,1.,1.,1.,0.,0.,0.,0.5,1.,1.,0.,0.25,0.5,0.75,0.5,0.5,0.,0.,0.,0.,1.,1.,1.,1.,0.,0.,0.5,0.5,1.,1.,0.5,0.5,0.],42,2),1e-13))
        m.changeSpaceDimension(3)
        m.getCoords()[:,2]=7.
        loc=f.getLocalizationOfDiscr()
        self.assertEqual(42,len(loc))
        self.assertTrue(loc.isEqual(DataArrayDouble([0.,0.,7.,1.,0.,7.,0.5,1.,7.,0.,0.,7.,1.,0.,7.,1.,1.,7.,0.,1.,7.,0.,0.,7.,1.,0.,7.,0.5,1.,7.,0.5,0.,7.,0.75,0.5,7.,0.25,0.5,7.,0.,0.,7.,1.,0.,7.,1.,1.,7.,0.,1.,7.,0.5,0.,7.,1.,0.5,7.,0.5,1.,7.,0.,0.5,7.,0.,0.,7.,0.5,1.,7.,1.,0.,7.,0.,0.,7.,0.,1.,7.,1.,1.,7.,1.,0.,7.,0.,0.,7.,0.5,1.,7.,1.,0.,7.,0.25,0.5,7.,0.75,0.5,7.,0.5,0.,7.,0.,0.,7.,0.,1.,7.,1.,1.,7.,1.,0.,7.,0.,0.5,7.,0.5,1.,7.,1.,0.5,7.,0.5,0.,7.],42,3),1e-13))
        pass

    def testSwig2GaussMeasureAndIntegral(self):
        ft=MEDCouplingDataForTest.buildFieldOnGauss_1()
        mea=ft.buildMeasureField(False)
        mea.checkCoherency()
        self.assertTrue(mea.getArray().isEqual(DataArrayDouble([-0.08504076274779823,-0.06378057206084897,-0.08504076274779869,-0.10630095343474463,-0.12756114412169625,-0.10630095343474734,-0.0637805720608491,-0.0850407627477968,-0.1063009534347449,-0.0850407627477994,-0.10630095343474809,-0.1275611441216954,-0.037205333702161475,-0.037205333702161475,-0.037205333702161475,-0.037205333702161475,-0.047835429045636084,-0.047835429045636084,-0.047835429045636084,-0.047835429045636084,-0.05846552438911087,-0.05846552438911087,-0.05846552438911087,-0.05846552438911087,-0.037205333702161725,-0.037205333702161725,-0.037205333702161725,-0.037205333702161725,-0.047835429045635834,-0.047835429045635834,-0.047835429045635834,-0.047835429045635834,-0.05846552438911058,-0.05846552438911058,-0.05846552438911058,-0.05846552438911058,-0.03879154890291829,-0.03879154890291829,-0.03879154890291829,-0.04120270848015563,-0.04120270848015563,-0.04120270848015563,-0.03393028948486933,-0.03393028948486933,-0.03393028948486933,-0.03151955746491709,-0.03151955746491709,-0.03151955746491709,-0.02424752187358276,-0.02424752187358276,-0.02424752187358276,-0.026657914642918758,-0.026657914642918758,-0.026657914642918758,-0.04120270848015456,-0.04120270848015456,-0.04120270848015456,-0.03879154890291757,-0.03879154890291757,-0.03879154890291757,-0.031519557464916595,-0.031519557464916595,-0.031519557464916595,-0.03393028948487046,-0.03393028948487046,-0.03393028948487046,-0.0266579146429191,-0.0266579146429191,-0.0266579146429191,-0.024247521873582645,-0.024247521873582645,-0.024247521873582645,-0.01851718920904466,-0.01851718920904466,-0.01851718920904466,-0.01851718920904466,-0.029627502734471456,-0.029627502734471456,-0.029627502734471456,-0.029627502734471456,-0.04740400437515433,-0.015150427534672922,-0.015150427534672922,-0.015150427534672922,-0.015150427534672922,-0.024240684055476674,-0.024240684055476674,-0.024240684055476674,-0.024240684055476674,-0.038785094488762675,-0.011783665860301345,-0.011783665860301345,-0.011783665860301345,-0.011783665860301345,-0.018853865376482152,-0.018853865376482152,-0.018853865376482152,-0.018853865376482152,-0.030166184602371443,-0.018517189209044892,-0.018517189209044892,-0.018517189209044892,-0.018517189209044892,-0.029627502734471827,-0.029627502734471827,-0.029627502734471827,-0.029627502734471827,-0.04740400437515492,-0.015150427534672776,-0.015150427534672776,-0.015150427534672776,-0.015150427534672776,-0.02424068405547644,-0.02424068405547644,-0.02424068405547644,-0.02424068405547644,-0.03878509448876231,-0.011783665860301277,-0.011783665860301277,-0.011783665860301277,-0.011783665860301277,-0.01885386537648204,-0.01885386537648204,-0.01885386537648204,-0.01885386537648204,-0.030166184602371266]),1e-14))
        f=MEDCouplingFieldDouble(ft)
        arr=DataArrayDouble(126,2)
        arr[:,0]=range(126)
        arr[:,1]=range(126)
        arr[:,1]+=1000
        f.setArray(arr)
        f.checkCoherency()
        self.assertTrue(DataArrayDouble(f.integral(False)).isEqual(DataArrayDouble([-211.66121638700983,-4863.9563007698835]),1e-11))
        self.assertTrue(DataArrayDouble(f.getWeightedAverageValue()).isEqual(DataArrayDouble([45.496085813113595,1045.496085813114]),1e-11))
        self.assertTrue(DataArrayDouble(f.normL1()).isEqual(DataArrayDouble([45.49608581311362,1045.496085813114]),1e-11))
        self.assertTrue(DataArrayDouble(f.normL2()).isEqual(DataArrayDouble([58.16846378340898,1046.1241521947334]),1e-11))
        pass

    def testSwig2FieldDiscretizationComputeMeshRestrictionFromTupleIds1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_GAUSS_NE)
        f.setMesh(m)
        a,b=f.getDiscretization().computeMeshRestrictionFromTupleIds(f.getMesh(),[3,4,5,6,8,9,10,14,15,16,17])
        self.assertTrue(a.isEqual(DataArrayInt([1,4])))
        self.assertTrue(b.isEqual(DataArrayInt([4,5,6,14,15,16,17])))
        a,b=f.getDiscretization().computeMeshRestrictionFromTupleIds(f.getMesh(),DataArrayInt([0,1,2,3,5,7,8,9,10,11,12,18]))
        self.assertTrue(a.isEqual(DataArrayInt([0,2])))
        self.assertTrue(b.isEqual(DataArrayInt([0,1,2,3,7,8,9])))
        #
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        a,b=f.getDiscretization().computeMeshRestrictionFromTupleIds(f.getMesh(),[3,4])
        self.assertTrue(a.isEqual(DataArrayInt([3,4])))
        self.assertTrue(b.isEqual(DataArrayInt([3,4])))
        #
        f=MEDCouplingFieldDouble(ON_NODES)
        f.setMesh(m)
        a,b=f.getDiscretization().computeMeshRestrictionFromTupleIds(f.getMesh(),[1,2,3,4])
        self.assertTrue(a.isEqual(DataArrayInt([1])))
        self.assertTrue(b.isEqual(DataArrayInt([1,2,4])))
        #
        f=MEDCouplingDataForTest.buildFieldOnGauss_1()
        a,b=f.getDiscretization().computeMeshRestrictionFromTupleIds(f.getMesh(),[0,11,12,13,14,15,17,18,19,36,37,38,115,117,118,119,120,121,122,123,124,125])
        self.assertTrue(a.isEqual(DataArrayInt([0,11,12,18,35])))
        self.assertTrue(b.isEqual(DataArrayInt([0,11,12,13,14,15,36,37,38,117,118,119,120,121,122,123,124,125])))
        #
        d=DataArrayInt([0,3,7,9,15,18])
        e=DataArrayInt([0,1,2,3,7,8,15,16,17])
        a,b=d.searchRangesInListOfIds(e)
        self.assertTrue(a.isEqual(DataArrayInt([0,2,4])))
        self.assertTrue(b.isEqual(DataArrayInt([0,1,2,7,8,15,16,17])))
        pass
    
    def testSwig2BigMem(self):
        if MEDCouplingSizeOfVoidStar()==64:
            d=DataArrayAsciiChar(223456789,16)
            self.assertTrue(d.getNumberOfTuples(),223456789)
            self.assertTrue(d.getNumberOfComponents(),16)
            d.setIJ(223456788,5,"r")
            self.assertTrue(d.getIJ(223456788,5),'r')
            d[223456787]="1234567890123456"
            self.assertTrue(d[223456787],'1234567890123456')
            self.assertRaises(InterpKernelException,d.rearrange,1)# fails because it would lead to nb of tuples > 2147483647
            pass
        pass

    def testSwig2DAReverseMultiCompo1(self):
        d=DataArrayDouble(6,2)
        d[:,0]=range(6)
        d[:,1]=range(10,16)
        d.reverse()
        self.assertTrue(d.isEqual(DataArrayDouble([5.,15.,4.,14.,3.,13.,2.,12.,1.,11.,0.,10.],6,2),1e-14))
        d=DataArrayDouble(7,2)
        d[:,0]=range(7)
        d[:,1]=range(10,17)
        d.reverse()
        self.assertTrue(d.isEqual(DataArrayDouble([6.,16.,5.,15.,4.,14.,3.,13.,2.,12.,1.,11.,0.,10.],7,2),1e-14))
        #
        d=DataArrayInt(6,2)
        d[:,0]=range(6)
        d[:,1]=range(10,16)
        d.reverse()
        self.assertTrue(d.isEqual(DataArrayInt([5,15,4,14,3,13,2,12,1,11,0,10],6,2)))
        d=DataArrayInt(7,2)
        d[:,0]=range(7)
        d[:,1]=range(10,17)
        d.reverse()
        self.assertTrue(d.isEqual(DataArrayInt([6,16,5,15,4,14,3,13,2,12,1,11,0,10],7,2)))
        pass

    def testSwigDAPow1(self):
        d=DataArrayInt(10)
        d.iota(0)
        d1=d.deepCpy()
        d.setIJ(2,0,-2)
        self.assertTrue((d**2).isEqual(DataArrayInt([0,1,4,9,16,25,36,49,64,81])))
        self.assertTrue((d**3).isEqual(DataArrayInt([0,1,-8,27,64,125,216,343,512,729])))
        for elt in [d]:
            elt**=2
            pass
        self.assertTrue(d.isEqual(DataArrayInt([0,1,4,9,16,25,36,49,64,81])))
        self.assertTrue((d1[:4]**d1[:4]).isEqual(DataArrayInt([1,1,4,27])))
        self.assertTrue((3**d1[:4]).isEqual(DataArrayInt([1,3,9,27])))
        d2=d1[:4]
        d2**=d2
        self.assertTrue(d2.isEqual(DataArrayInt([1,1,4,27])))
        self.assertRaises(InterpKernelException,d2.__pow__,-1)#non supporting negative pow in DataArrayInt.__pow__
        self.assertRaises(InterpKernelException,d2.__ipow__,-1)#non supporting negative pow in DataArrayInt.__pow__
        #
        d=DataArrayDouble(10)
        d.iota(0)
        d1=d.deepCpy()
        d.setIJ(2,0,-2.)
        self.assertTrue((d**2).isEqual(DataArrayDouble([0,1,4,9,16,25,36,49,64,81]),1e-12))
        self.assertTrue((d**3).isEqual(DataArrayDouble([0,1,-8,27,64,125,216,343,512,729]),1e-12))
        self.assertRaises(InterpKernelException,d.__pow__,3.1)#3.1 is double not integer -> not supporting negative values in d
        for elt in [d]:
            elt**=2
            pass
        self.assertTrue(d.isEqual(DataArrayDouble([0,1,4,9,16,25,36,49,64,81]),1e-12))
        self.assertTrue((d1[:4]**d1[:4]).isEqual(DataArrayDouble([1,1,4,27]),1e-12))
        self.assertTrue((3**d1[:4]).isEqual(DataArrayDouble([1,3,9,27]),1e-12))
        d2=d1[:4]
        d2**=d2
        self.assertTrue(d2.isEqual(DataArrayDouble([1,1,4,27]),1e-12))
        d2**=-0.5
        self.assertTrue(d2.isEqual(DataArrayDouble([1,1,1./2,1./sqrt(27.)]),1e-14))
        d3=-1./d1[1:5]
        self.assertTrue((3**d3).isEqual(DataArrayDouble([0.3333333333333333,0.5773502691896257,0.6933612743506348,0.7598356856515925]),1e-14))
        d4=d3.deepCpy() ; d4.abs()
        self.assertTrue((d4**d3).isEqual(DataArrayDouble([1.,sqrt(2.),1.4422495703074083,sqrt(2.)]),1e-14))
        d4**=d3
        self.assertTrue(d4.isEqual(DataArrayDouble([1.,sqrt(2.),1.4422495703074083,sqrt(2.)]),1e-14))
        pass
    
    def testSwig2Baryenter3DForCellsWithVolumeZero1(self):
        coo=DataArrayDouble([0.,0.,0.,1.,0.,0.,0.,1.,0.],3,3)
        m2=MEDCouplingUMesh("mesh",2)
        m2.allocateCells(0)
        m2.insertNextCell(NORM_POLYGON,[0,1,2])
        m2.setCoords(coo)
        m2.checkCoherency1()
        #
        coo2=DataArrayDouble([0.,0.,0.,0.,0.,0.,0.,0.,2.],3,3)
        m1=MEDCouplingUMesh("mesh",1)
        m1.allocateCells(0)
        m1.insertNextCell(NORM_SEG2,[0,1])
        m1.insertNextCell(NORM_SEG2,[1,2])
        m1.setCoords(coo2)
        m1.checkCoherency1()
        #
        m3=m2.buildExtrudedMesh(m1,0)
        m3.insertNextCell(NORM_POLYHED,[3,4,5,-1,8,7,6,-1,4,3,6,7,-1,5,4,7,8,-1,5,4,-1,3,5,8,6])# addition of face #4 with null surface
        self.assertTrue(m3.getBarycenterAndOwner().isEqual(DataArrayDouble([0.3333333333333333,0.3333333333333333,0.,0.3333333333333333,0.3333333333333333,1.,0.3333333333333333,0.3333333333333333,1.],3,3),1e-13))
        m4,a,b,c,d=m3.buildDescendingConnectivity()
        self.assertTrue(m4.getBarycenterAndOwner().isEqual(DataArrayDouble([0.3333333333333333,0.3333333333333333,0.,0.3333333333333333,0.3333333333333333,0.,0.5,0.,0.,0.5,0.5,0.,0.,0.5,0.,0.3333333333333333,0.3333333333333333,2.,0.5,0.,1.,0.5,0.5,1.,0.,0.5,1.,0.5,0.5,0.],10,3),1e-13))
        pass

    def testSwigRepr1(self):
        d=DataArrayDouble()
        self.assertTrue(len(d.__repr__())<120)
        d.alloc(1000,0) ; self.assertTrue(len(d.__repr__())<100)
        for i in xrange(100):
            d.alloc(i,1) ; d.iota(1.1234567890123456) ; d*=1e123
            self.assertTrue(len(d.__repr__())<500)
            pass
        for i in xrange(50):
            d.alloc(i,2) ; d.rearrange(1) ; d.iota(1.1234567890123456) ; d.rearrange(2) ; d*=1e123
            self.assertTrue(len(d.__repr__())<500)
            pass
        d.alloc(4000,1) ; d.iota() ; self.assertTrue(len(d.__repr__())<500)
        for i in xrange(2,4):
            d.alloc(362880,1) ; d.iota() ; d.rearrange(i) ; self.assertTrue(len(d.__repr__())<500)
            pass
        d.alloc(0,9)
        self.assertTrue(len(d.__repr__())<120)
        #
        d=DataArrayInt()
        self.assertTrue(len(d.__repr__())<100)
        d.alloc(1000,0) ; self.assertTrue(len(d.__repr__())<100)
        for i in xrange(100):
            d.alloc(i,1) ; d.iota(123456789)
            self.assertTrue(len(d.__repr__())<500)
            pass
        for i in xrange(50):
            d.alloc(i,2) ; d.rearrange(1) ; d.iota(123456789) ; d.rearrange(2)
            self.assertTrue(len(d.__repr__())<500)
            pass
        d.alloc(4000,1) ; d.iota() ; self.assertTrue(len(d.__repr__())<500)
        for i in xrange(2,10):
            d.alloc(362880,1) ; d.iota() ; d.rearrange(i) ; self.assertTrue(len(d.__repr__())<500)
            pass
        d.alloc(0,9)
        self.assertTrue(len(d.__repr__())<100)
        #
        d=DataArrayAsciiChar()
        d.alloc(1000,0) ; self.assertTrue(len(d.__repr__())<100)
        d.alloc(2,16) ; d[:]='1234567890ABCDEF'
        self.assertTrue(len(d.__repr__())<500)
        d.alloc(2000,16) ; d[:]='1234567890ABCDEF'
        self.assertTrue(len(d.__repr__())<500)
        d.alloc(0,16) ; d[:]='1234567890ABCDEF'
        self.assertTrue(len(d.__repr__())<120)
        #
        d=DataArrayByte()
        self.assertTrue(len(d.__repr__())<100)
        d.alloc(1000,0) ; self.assertTrue(len(d.__repr__())<100)
        d.alloc(0,16) ; self.assertTrue(len(d.__repr__())<100)
        d.alloc(5,1) ; d.fillWithValue(127)
        self.assertTrue(len(d.__repr__())<200)
        d.alloc(1000,1) ; d.fillWithValue(127)
        self.assertTrue(len(d.__repr__())<500)
        d.alloc(1000,3) ; d.fillWithValue(127)
        self.assertTrue(len(d.__repr__())<500)
        pass
    
    def testSwig2MeshComputeIsoBarycenterOfNodesPerCell1(self):
        coo=DataArrayDouble([26.17509821414239,5.0374,200.,26.175098214142388,-5.0374,200.,17.450065476094927,20.1496,200.,8.725032738047464,25.187,200.,43.62516369023732,5.0374,200.,34.90013095218986,10.0748,200.,34.900130952189855,-10.0748,200.,43.625163690237315,-5.0374,200.,26.175098214142402,25.187,200.,26.175098214142395,35.2618,200.,17.45006547609493,40.2992,200.,8.725032738047469,35.2618,200.,26.17509821414239,5.0374,200.,26.175098214142388,-5.0374,200.,17.450065476094927,20.1496,200.,8.725032738047464,25.187,200.,43.62516369023732,5.0374,200.,34.90013095218986,10.0748,200.,34.900130952189855,-10.0748,200.,43.625163690237315,-5.0374,200.,26.175098214142402,25.187,200.,26.175098214142395,35.2618,200.,17.45006547609493,40.2992,200.,8.725032738047469,35.2618,200.],24,3)
        m=MEDCouplingUMesh.New("toto",3)
        m.allocateCells(0)
        m.insertNextCell(NORM_POLYHED,[4,5,0,1,6,7,-1,19,18,13,12,17,16,-1,5,4,16,17,-1,0,5,17,12,-1,1,0,12,13,-1,6,1,13,18,-1,7,6,18,19,-1,4,7,19,16])
        m.insertNextCell(NORM_POLYHED,[9,10,11,3,2,8,-1,20,14,15,23,22,21,-1,10,9,21,22,-1,11,10,22,23,-1,3,11,23,15,-1,2,3,15,14,-1,8,2,14,20,-1,9,8,20,21])
        m.setCoords(coo)
        m.checkCoherency1()
        #
        dReference=DataArrayDouble([(34.900130952189848,0.,200),(17.450065476094931,30.2244,200.)])
        self.assertTrue(m.computeIsoBarycenterOfNodesPerCell().isEqual(dReference,1e-12))
        m.getNodalConnectivity().setIJ(87,0,24)
        self.assertRaises(InterpKernelException,m.computeIsoBarycenterOfNodesPerCell)
        m.getNodalConnectivity().setIJ(87,0,-2)
        self.assertRaises(InterpKernelException,m.computeIsoBarycenterOfNodesPerCell)
        m.getNodalConnectivity().setIJ(87,0,21)# put again 21 as at the beginning
        #
        self.assertTrue(m.unPolyze())
        self.assertEqual([NORM_HEXGP12],m.getAllGeoTypes())
        self.assertTrue(m.computeIsoBarycenterOfNodesPerCell().isEqual(dReference,1e-12))
        m.getNodalConnectivity().setIJ(25,0,24)
        self.assertRaises(InterpKernelException,m.computeIsoBarycenterOfNodesPerCell)
        m.getNodalConnectivity().setIJ(25,0,-1)
        self.assertRaises(InterpKernelException,m.computeIsoBarycenterOfNodesPerCell)
        pass

    def testSwig2NonRegressionBugDescHexa20(self):
        coo=DataArrayDouble([0.,0.,0.,1.23,0.,0.,0.615,0.,0.,0.,2.1,0.,0.615,2.1,0.,1.23,2.1,0.,1.23,1.05,0.,0.,1.05,0.,0.,0.,2.16,1.23,0.,2.16,1.23,2.1,2.16,0.,2.1,2.16,0.,0.,4.32,0.615,0.,4.32,1.23,0.,4.32,1.23,1.05,4.32,1.23,2.1,4.32,0.615,2.1,4.32,0.,2.1,4.32,0.,1.05,4.32],20,3)
        m=MEDCouplingUMesh('mesh',3)
        m.allocateCells(0)
        m.insertNextCell(NORM_HEXA20,[0,3,5,1,12,18,16,14,7,4,6,2,19,17,15,13,8,11,10,9])
        m.setCoords(coo)
        m.checkCoherency1()
        #
        a,b,c,d,e=m.buildDescendingConnectivity()
        m2=MEDCouplingUMesh('mesh',2)
        m2.allocateCells(0)
        m2.setCoords(coo)
        conn2=[[0,3,5,1,7,4,6,2],[12,14,16,18,13,15,17,19],[0,12,18,3,8,19,11,7],[3,18,16,5,11,17,10,4],[5,16,14,1,10,15,9,6],[1,14,12,0,9,13,8,2]]
        for i in xrange(6):
            m2.insertNextCell(NORM_QUAD8,conn2[i])
            pass
        self.assertTrue(m2.isEqual(a,1e-12))
        self.assertTrue(b.isEqual(DataArrayInt([0,1,2,3,4,5])))
        self.assertTrue(c.isEqual(DataArrayInt([0,6])))
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,0,0])))
        self.assertTrue(e.isEqual(DataArrayInt([0,1,2,3,4,5,6])))
        #
        m.convertQuadraticCellsToLinear() ; m.zipCoords()
        m.convertLinearCellsToQuadratic(1)
        #
        coo2=DataArrayDouble([0.,0.,0.,1.23,0.,0.,0.,2.1,0.,1.23,2.1,0.,0.,0.,4.32,1.23,0.,4.32,1.23,2.1,4.32,0.,2.1,4.32,0.,1.05,0.,0.615,2.1,0.,1.23,1.05,0.,0.615,0.,0.,0.,1.05,4.32,0.615,2.1,4.32,1.23,1.05,4.32,0.615,0.,4.32,0.,0.,2.16,0.,2.1,2.16,1.23,2.1,2.16,1.23,0.,2.16,0.615,1.05,0.,0.,1.05,2.16,0.615,2.1,2.16,1.23,1.05,2.16,0.615,0.,2.16,0.615,1.05,4.32,0.615,1.05,2.16],27,3)
        m3=MEDCouplingUMesh("mesh",3)
        m3.allocateCells(1)
        m3.insertNextCell(NORM_HEXA27,[0,2,3,1,4,7,6,5,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26])
        m3.setCoords(coo2)
        self.assertTrue(m3.isEqual(m,1e-12))
        #
        a,b,c,d,e=m.buildDescendingConnectivity()
        conn4=[[0,2,3,1,8,9,10,11,20],[4,5,6,7,15,14,13,12,25],[0,4,7,2,16,12,17,8,21],[2,7,6,3,17,13,18,9,22],[3,6,5,1,18,14,19,10,23],[1,5,4,0,19,15,16,11,24]]
        m4=MEDCouplingUMesh("mesh",2)
        m4.allocateCells(0)
        for i in xrange(6):
            m4.insertNextCell(NORM_QUAD9,conn4[i])
            pass
        m4.setCoords(coo2)
        self.assertTrue(m4.isEqual(a,1e-12))
        self.assertTrue(b.isEqual(DataArrayInt([0,1,2,3,4,5])))
        self.assertTrue(c.isEqual(DataArrayInt([0,6])))
        self.assertTrue(d.isEqual(DataArrayInt([0,0,0,0,0,0])))
        self.assertTrue(e.isEqual(DataArrayInt([0,1,2,3,4,5,6])))
        pass
    
    def testSwigAdvGauss(self):
        f=MEDCouplingFieldTemplate(ON_GAUSS_PT)
        f.setDiscretization(None)
        f.__repr__() ; f.__str__()
        #
        f=MEDCouplingFieldTemplate(ON_GAUSS_PT)
        d=f.getDiscretization()
        i=DataArrayInt() ; i.alloc(10,1) ; i.iota(1)
        d.setArrayOfDiscIds(i)
        f.__repr__() ; f.__str__()
        i2=d.getArrayOfDiscIds()
        self.assertEqual(i.__repr__(),i2.__repr__())
        #
        f=MEDCouplingFieldDouble(ON_GAUSS_PT)
        f.setDiscretization(None)
        f.__repr__() ; f.__str__()
        #
        f=MEDCouplingFieldDouble(ON_GAUSS_PT)
        d=f.getDiscretization()
        i=DataArrayInt() ; i.alloc(10,1) ; i.iota(1)
        d.setArrayOfDiscIds(i)
        f.__repr__() ; f.__str__()
        #
        gl=MEDCouplingGaussLocalization(NORM_SEG2,[0,1],[0.5],[1.])
        gl.setWeights([3.])
        gl.__repr__() ; gl.__str__()
        gl=MEDCouplingGaussLocalization(NORM_ERROR)
        gl.setWeights([3.])
        gl.__repr__() ; gl.__str__()
        pass

    def testSwig2NonRegressionBugSubstractInPlaceDM(self):
        m0=MEDCouplingCMesh()
        arr=DataArrayDouble(5,1) ; arr.iota(0.)
        m0.setCoords(arr,arr)
        m0=m0.buildUnstructured()
        m00=m0[::2] ; m00.simplexize(0) ; m01=m0[1::2]
        m0=MEDCouplingUMesh.MergeUMeshes([m00,m01])
        m0.getCoords()[:]*=1/4.
        m0.setName("mesh")
        #
        NodeField=MEDCouplingFieldDouble(ON_NODES,ONE_TIME) ; NodeField.setTime(5.6,5,6) ; NodeField.setMesh(m0)
        NodeField.setName("NodeField")
        NodeField.fillFromAnalytic(1,"exp(-((x-1)*(x-1)+(y-1)*(y-1)))") ; NodeField.getArray().setInfoOnComponent(0,"powernode [W]")
        proc0=m0.getCellsInBoundingBox([(0.,0.4),(0.,0.4)],1e-10)
        proc1=proc0.buildComplement(m0.getNumberOfCells())
        #
        NodeField0=NodeField[proc0] ; NodeField0.getMesh().setName(m0.getName())
        NodeField1=NodeField[proc1] ; NodeField1.getMesh().setName(m0.getName())
        #
        NodeField_read=MEDCouplingFieldDouble.MergeFields([NodeField0,NodeField1])
        NodeField_read.mergeNodes(1e-10)
        NodeFieldCpy=NodeField.deepCpy()
        NodeFieldCpy.mergeNodes(1e-10)
        NodeField.checkCoherency()
        self.assertTrue(not NodeField.getArray().isUniform(0.,1e-12))
        NodeField.substractInPlaceDM(NodeField_read,10,1e-12)
        self.assertTrue(NodeField.getArray().isUniform(0.,1e-12))
        pass

    def testSwigFieldOperationOpen1(self):
        ## MEDCouplingFieldDouble.__add__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(5) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__add__,2)
        self.assertRaises(InterpKernelException,f.__add__,range(5))
        self.assertRaises(InterpKernelException,f.__add__,arr)
        self.assertRaises(InterpKernelException,f.__add__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__add__,2)
        self.assertRaises(InterpKernelException,f.__add__,range(5))
        self.assertRaises(InterpKernelException,f.__add__,arr)
        self.assertRaises(InterpKernelException,f.__add__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(5) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        ff=f+2
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(2,9),(3,10),(4,11),(5,12),(6,13)]),1e-12))
        ff=f+arr
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,7),(2,10),(4,13),(6,16),(8,19)]),1e-12))
        self.assertRaises(InterpKernelException,f.__add__,f2)
        f2.setArray(arr)
        ff=f+f2
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,7),(2,10),(4,13),(6,16),(8,19)]),1e-12))
        ff=f+[5,8]
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(5,15),(6,16),(7,17),(8,18),(9,19)]),1e-12))
        ### MEDCouplingFieldDouble.__sub__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(5) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__sub__,2)
        self.assertRaises(InterpKernelException,f.__sub__,range(5))
        self.assertRaises(InterpKernelException,f.__sub__,arr)
        self.assertRaises(InterpKernelException,f.__sub__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__sub__,2)
        self.assertRaises(InterpKernelException,f.__sub__,range(5))
        self.assertRaises(InterpKernelException,f.__sub__,arr)
        self.assertRaises(InterpKernelException,f.__sub__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(5) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        ff=f-2
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(-2,5),(-1,6),(0,7),(1,8),(2,9)]),1e-12))
        ff=f-arr
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,7),(0,6),(0,5),(0,4),(0,3)]),1e-12))
        self.assertRaises(InterpKernelException,f.__sub__,f2)
        f2.setArray(arr)
        ff=f-f2
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,7),(0,6),(0,5),(0,4),(0,3)]),1e-12))
        ff=f-[5,8]
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(-5,-1),(-4,0),(-3,1),(-2,2),(-1,3)]),1e-12))
        ### MEDCouplingFieldDouble.__mul__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(5) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__mul__,2)
        self.assertRaises(InterpKernelException,f.__mul__,range(5))
        self.assertRaises(InterpKernelException,f.__mul__,arr)
        self.assertRaises(InterpKernelException,f.__mul__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__mul__,2)
        self.assertRaises(InterpKernelException,f.__mul__,range(5))
        self.assertRaises(InterpKernelException,f.__mul__,arr)
        self.assertRaises(InterpKernelException,f.__mul__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(5) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        ff=f*2
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,14),(2,16),(4,18),(6,20),(8,22)]),1e-12))
        ff=f*arr
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,0),(1,16),(4,36),(9,60),(16,88)]),1e-12))
        self.assertRaises(InterpKernelException,f.__mul__,f2)
        f2.setArray(arr)
        ff=f*f2
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,0),(1,16),(4,36),(9,60),(16,88)]),1e-12))
        ff=f*[5,8]
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,56),(5,64),(10,72),(15,80),(20,88)]),1e-12))
        ### MEDCouplingFieldDouble.__div__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(1,6) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__div__,2)
        self.assertRaises(InterpKernelException,f.__div__,range(5))
        self.assertRaises(InterpKernelException,f.__div__,arr)
        self.assertRaises(InterpKernelException,f.__div__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__div__,2)
        self.assertRaises(InterpKernelException,f.__div__,range(5))
        self.assertRaises(InterpKernelException,f.__div__,arr)
        self.assertRaises(InterpKernelException,f.__div__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(5) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        self.assertRaises(InterpKernelException,f.__div__,0)
        ff=f/2
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,3.5),(0.5,4),(1,4.5),(1.5,5),(2,5.5)]),1e-12))
        ff=f/arr
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,3.5),(0.5,2),(0.6666666666666666,1.5),(0.75,1.25),(0.8,1.1)]),1e-12))
        self.assertRaises(InterpKernelException,f.__div__,f2)
        f2.setArray(arr)
        ff=f/f2
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,3.5),(0.5,2),(0.6666666666666666,1.5),(0.75,1.25),(0.8,1.1)]),1e-12))
        ff=f/[5,8]
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,0.875),(0.2,1),(0.4,1.125),(0.6,1.25),(0.8,1.375)]),1e-12))
        ### MEDCouplingFieldDouble.__pow__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5)
        arr[:]=[1,1,3,2,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__div__,2)
        self.assertRaises(InterpKernelException,f.__div__,range(5))
        self.assertRaises(InterpKernelException,f.__div__,arr)
        self.assertRaises(InterpKernelException,f.__div__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__div__,2)
        self.assertRaises(InterpKernelException,f.__div__,range(5))
        self.assertRaises(InterpKernelException,f.__div__,arr)
        self.assertRaises(InterpKernelException,f.__div__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,1)
        f.getArray()[:]=range(2,7)
        ff=f**2
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([4,9,16,25,36]),1e-12))
        ff=f**arr
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([2,3,64,25,1]),1e-12))
        f2.setArray(arr)
        ff=f**f2
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([2,3,64,25,1]),1e-12))
        ## MEDCouplingFieldDouble.__iadd__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(5) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__iadd__,2)
        self.assertRaises(InterpKernelException,f.__iadd__,range(5))
        self.assertRaises(InterpKernelException,f.__iadd__,arr)
        self.assertRaises(InterpKernelException,f.__iadd__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__iadd__,2)
        self.assertRaises(InterpKernelException,f.__iadd__,range(5))
        self.assertRaises(InterpKernelException,f.__iadd__,arr)
        self.assertRaises(InterpKernelException,f.__iadd__,f2)
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(5) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        f.checkCoherency()
        f+=2
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(2,9),(3,10),(4,11),(5,12),(6,13)]),1e-12))
        f+=arr
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(2,9),(4,12),(6,15),(8,18),(10,21)]),1e-12))
        f2.setArray(arr)
        f+=f2
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(2,9),(5,14),(8,19),(11,24),(14,29)]),1e-12))
        f+=[0.1,0.2]
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(2.1,9.2),(5.1,14.2),(8.1,19.2),(11.1,24.2),(14.1,29.2)]),1e-12))
        ## MEDCouplingFieldDouble.__isub__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(5) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__isub__,2)
        self.assertRaises(InterpKernelException,f.__isub__,range(5))
        self.assertRaises(InterpKernelException,f.__isub__,arr)
        self.assertRaises(InterpKernelException,f.__isub__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__isub__,2)
        self.assertRaises(InterpKernelException,f.__isub__,range(5))
        self.assertRaises(InterpKernelException,f.__isub__,arr)
        self.assertRaises(InterpKernelException,f.__isub__,f2)
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(5) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        f.checkCoherency()
        f-=2
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(-2,5),(-1,6),(0,7),(1,8),(2,9)]),1e-12))
        f-=arr
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(-2,5),(-2,4),(-2,3),(-2,2),(-2,1)]),1e-12))
        f2.setArray(arr)
        f-=f2
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(-2,5),(-3,2),(-4,-1),(-5,-4),(-6,-7)]),1e-12))
        f-=[0.1,0.2]
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(-2.1,4.8),(-3.1,1.8),(-4.1,-1.2),(-5.1,-4.2),(-6.1,-7.2)]),1e-12))
        ## MEDCouplingFieldDouble.__imul__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(5) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__imul__,2)
        self.assertRaises(InterpKernelException,f.__imul__,range(5))
        self.assertRaises(InterpKernelException,f.__imul__,arr)
        self.assertRaises(InterpKernelException,f.__imul__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__imul__,2)
        self.assertRaises(InterpKernelException,f.__imul__,range(5))
        self.assertRaises(InterpKernelException,f.__imul__,arr)
        self.assertRaises(InterpKernelException,f.__imul__,f2)
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(5) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        f.checkCoherency()
        f*=2
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,14),(2,16),(4,18),(6,20),(8,22)]),1e-12))
        f*=arr
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,0),(2,32),(8,72),(18,120),(32,176)]),1e-12))
        f2.setArray(arr)
        f*=f2
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,0),(2,64),(16,288),(54,720),(128,1408)]),1e-12))
        f*=[0.1,0.2]
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,0),(0.2,12.8),(1.6,57.6),(5.4,144),(12.8,281.6)]),1e-12))
        ## MEDCouplingFieldDouble.__idiv__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(1,6) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__idiv__,2)
        self.assertRaises(InterpKernelException,f.__idiv__,range(5))
        self.assertRaises(InterpKernelException,f.__idiv__,arr)
        self.assertRaises(InterpKernelException,f.__idiv__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__idiv__,2)
        self.assertRaises(InterpKernelException,f.__idiv__,range(5))
        self.assertRaises(InterpKernelException,f.__idiv__,arr)
        self.assertRaises(InterpKernelException,f.__idiv__,f2)
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(5) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        f.checkCoherency()
        f/=2
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,3.5),(0.5,4),(1,4.5),(1.5,5),(2,5.5)]),1e-12))
        f/=arr
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,1.75),(0.25,1),(0.3333333333333333,0.75),(0.375,0.625),(0.4,0.55)]),1e-12))
        f2.setArray(arr)
        f/=f2
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,0.875),(0.125,0.25),(0.1111111111111111,0.125),(0.09375,0.078125),(0.08,0.055)]),1e-12))
        f/=[0.1,0.2]
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,4.375),(1.25,1.25),(1.1111111111111111,0.625),(0.9375,0.390625),(0.8,0.275)]),1e-12))
        ## MEDCouplingFieldDouble.__ipow__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(1,6) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__ipow__,2)
        self.assertRaises(InterpKernelException,f.__ipow__,range(5))
        self.assertRaises(InterpKernelException,f.__ipow__,arr)
        self.assertRaises(InterpKernelException,f.__ipow__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__ipow__,2)
        self.assertRaises(InterpKernelException,f.__ipow__,range(5))
        self.assertRaises(InterpKernelException,f.__ipow__,arr)
        self.assertRaises(InterpKernelException,f.__ipow__,f2)
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(5) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        f.checkCoherency()
        f**=2
        f.checkCoherency()
        self.assertTrue(f.getArray().isEqual(DataArrayDouble([(0,49),(1,64),(4,81),(9,100),(16,121)]),1e-12))
         ## MEDCouplingFieldDouble.__radd__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(5) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__radd__,2)
        self.assertRaises(InterpKernelException,f.__radd__,range(5))
        self.assertRaises(InterpKernelException,f.__radd__,arr)
        self.assertRaises(InterpKernelException,f.__radd__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__radd__,2)
        self.assertRaises(InterpKernelException,f.__radd__,range(5))
        self.assertRaises(InterpKernelException,f.__radd__,arr)
        self.assertRaises(InterpKernelException,f.__radd__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(5) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        ff=2+f
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(2,9),(3,10),(4,11),(5,12),(6,13)]),1e-12))
        ff=arr+f
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,7),(2,10),(4,13),(6,16),(8,19)]),1e-12))
        self.assertRaises(InterpKernelException,f.__radd__,f2)
        ff=[5,8]+f
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(5,15),(6,16),(7,17),(8,18),(9,19)]),1e-12))
        ### MEDCouplingFieldDouble.__rsub__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(5) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__rsub__,2)
        self.assertRaises(InterpKernelException,f.__rsub__,range(5))
        self.assertRaises(InterpKernelException,f.__rsub__,arr)
        self.assertRaises(InterpKernelException,f.__rsub__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__rsub__,2)
        self.assertRaises(InterpKernelException,f.__rsub__,range(5))
        self.assertRaises(InterpKernelException,f.__rsub__,arr)
        self.assertRaises(InterpKernelException,f.__rsub__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(5) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        ff=2-f
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(2,-5),(1,-6),(0,-7),(-1,-8),(-2,-9)]),1e-12))
        ff=arr-f
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,-7),(0,-6),(0,-5),(0,-4),(0,-3)]),1e-12))
        self.assertRaises(InterpKernelException,f.__rsub__,f2)
        ### MEDCouplingFieldDouble.__rmul__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(5) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__rmul__,2)
        self.assertRaises(InterpKernelException,f.__rmul__,range(5))
        self.assertRaises(InterpKernelException,f.__rmul__,arr)
        self.assertRaises(InterpKernelException,f.__rmul__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__rmul__,2)
        self.assertRaises(InterpKernelException,f.__rmul__,range(5))
        self.assertRaises(InterpKernelException,f.__rmul__,arr)
        self.assertRaises(InterpKernelException,f.__rmul__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(5) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        ff=2*f
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,14),(2,16),(4,18),(6,20),(8,22)]),1e-12))
        ff=arr*f
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,0),(1,16),(4,36),(9,60),(16,88)]),1e-12))
        self.assertRaises(InterpKernelException,f.__rmul__,f2)
        ff=f*[5,8]
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(0,56),(5,64),(10,72),(15,80),(20,88)]),1e-12))
        ### MEDCouplingFieldDouble.__rdiv__
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2)
        arr[:,0]=range(1,6) ; arr[:,1]=2*arr[:,0]
        f2=f.clone(True)
        self.assertRaises(InterpKernelException,f.__rdiv__,2)
        self.assertRaises(InterpKernelException,f.__rdiv__,range(5))
        self.assertRaises(InterpKernelException,f.__rdiv__,arr)
        self.assertRaises(InterpKernelException,f.__rdiv__,f2)
        f.setArray(DataArrayDouble())
        self.assertRaises(InterpKernelException,f.__rdiv__,2)
        self.assertRaises(InterpKernelException,f.__rdiv__,range(5))
        self.assertRaises(InterpKernelException,f.__rdiv__,arr)
        self.assertRaises(InterpKernelException,f.__rdiv__,f2)
        self.assertRaises(InterpKernelException,f.__getitem__,(slice(None),0))
        f.getArray().alloc(5,2)
        f.getArray()[:,0]=range(1,6) ; f.getArray()[:,1]=f.getArray()[:,0]+7
        ff=2/f
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(2,0.25),(1,0.22222222222222221),(0.66666666666666663,0.20000000000000001),(0.5,0.18181818181818182),(0.40000000000000002,0.16666666666666666)]),1e-12))
        ff=arr/f
        ff.checkCoherency()
        self.assertTrue(ff.getArray().isEqual(DataArrayDouble([(1,0.25),(1,0.44444444444444442),(1,0.59999999999999998),(1,0.72727272727272729),(1,0.83333333333333337)]),1e-12))
        self.assertRaises(InterpKernelException,f.__rdiv__,f2)
        pass
    
    def testSwig2FieldDoubleBuildSubPartRange1(self):
        #ON_CELLS
        m=MEDCouplingDataForTest.build2DTargetMesh_1()
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(m)
        arr=DataArrayDouble(5,2) ; arr[:,0]=range(7,12) ; arr[:,1]=100+arr[:,0]
        f.setArray(arr)
        f.checkCoherency()
        ff=f[1:-1:2]
        ff.checkCoherency()
        self.assertTrue((m.buildPartOfMySelf([1,3],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(2,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[1,3]],1e-12))
        #
        a,b=f.buildSubMeshDataRange(2,5,1)
        self.assertTrue(m.buildPartOfMySelf([2,3,4],True).isEqual(a,1e-12))
        self.assertEqual(b,slice(2,5,1))
        ff=f[2:]
        ff.checkCoherency()
        self.assertTrue((m.buildPartOfMySelf([2,3,4],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[2,3,4]],1e-12))
        #
        ff=f[-2:0:-1]
        ff.checkCoherency()
        self.assertTrue((m.buildPartOfMySelf([3,2,1],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[3,2,1]],1e-12))
        self.assertTrue(f[-2:0:-1,1].getArray().isEqual(arr[[3,2,1],1],1e-12))
        #ON_NODES
        f=MEDCouplingFieldDouble(ON_NODES)
        f.setMesh(m)
        arr=DataArrayDouble(9,2) ; arr[:,0]=range(7,16) ; arr[:,1]=100+arr[:,0]
        f.setArray(arr)
        f.checkCoherency()
        ff=f[1:-1:2]
        ff.checkCoherency()
        self.assertTrue((m.buildPartOfMySelf([1,3],False)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(6,ff.getMesh().getNumberOfNodes())
        self.assertTrue(2,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[1,2,3,4,6,7]],1e-12))
        #
        m2=m.buildPartRange(2,5,1)
        self.assertTrue(m.buildPartOfMySelf([2,3,4],True).isEqual(m2,1e-12))
        m2,b=m.buildPartRangeAndReduceNodes(2,5,1)
        self.assertTrue(m.buildPartOfMySelf([2,3,4],False).isEqual(m2,1e-12))
        self.assertTrue(b.isEqual(DataArrayInt([-1,-1,0,1,2,3,4,5,6])))
        a,b=f.buildSubMeshDataRange(2,5,1)
        self.assertTrue(m.buildPartOfMySelf([2,3,4],False).isEqual(a,1e-12))
        self.assertTrue(b.isEqual(DataArrayInt([2,3,4,5,6,7,8])))
        ff=f[2:]
        ff.checkCoherency()
        self.assertTrue((m.buildPartOfMySelf([2,3,4],False)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(7,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[2,3,4,5,6,7,8]],1e-12))
        #
        ff=f[-2:0:-1]
        ff.checkCoherency()
        self.assertTrue((m.buildPartOfMySelf([3,2,1],False)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(7,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[1,2,3,4,5,6,7]],1e-12))
        self.assertTrue(f[-2:0:-1,1].getArray().isEqual(arr[[1,2,3,4,5,6,7],1],1e-12))
        #ON_GAUSS_NE
        f=MEDCouplingFieldDouble(ON_GAUSS_NE)
        f.setMesh(m)
        arr=DataArrayDouble(18,2) ; arr[:,0]=range(7,25) ; arr[:,1]=100+arr[:,0]
        f.setArray(arr)
        f.checkCoherency()
        ff=f[1:-1:2]
        ff.checkCoherency()
        self.assertTrue((m.buildPartOfMySelf([1,3],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(2,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[4,5,6,10,11,12,13]],1e-12))
        #
        a,b=f.buildSubMeshDataRange(2,5,1)
        self.assertTrue(m.buildPartOfMySelf([2,3,4],True).isEqual(a,1e-12))
        self.assertEqual(b,slice(7,18,1))
        ff=f[2:]
        ff.checkCoherency()
        self.assertTrue((m.buildPartOfMySelf([2,3,4],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[7,8,9,10,11,12,13,14,15,16,17]],1e-12))
        #
        ff=f[-2:0:-1]
        ff.checkCoherency()
        self.assertTrue((m.buildPartOfMySelf([3,2,1],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[10,11,12,13,7,8,9,4,5,6]],1e-12))
        self.assertTrue(f[-2:0:-1,1].getArray().isEqual(arr[[10,11,12,13,7,8,9,4,5,6],1],1e-12))
        #ON_GAUSS_PT
        f=MEDCouplingFieldDouble(ON_GAUSS_PT)
        f.setMesh(m)
        f.setGaussLocalizationOnCells([0,4],[0,0,1,0,1,1,1,0],[1.1,1.1,2.2,2.2],[0.2,0.8]);
        f.setGaussLocalizationOnCells([3],[0,0,1,0,1,1,1,0],[1.1,1.1,2.2,2.2,3.,3.],[0.2,0.4,0.4]);
        f.setGaussLocalizationOnCells([1],[0,0,1,0,1,0],[1.1,1.1,2.2,2.2,3.,3.,4.,4.],[0.1,0.1,0.4,0.4]);
        f.setGaussLocalizationOnCells([2],[0,0,1,0,1,0],[1.1,1.1,2.2,2.2,3.,3.,4.,4.,5.,5.],[0.1,0.1,0.4,0.3,0.1]);
        arr=DataArrayDouble(16,2) ; arr[:,0]=range(7,23) ; arr[:,1]=100+arr[:,0]
        f.setArray(arr)
        f.checkCoherency()
        ff=f[1:-1:2]
        ff.checkCoherency()
        self.assertTrue((m.buildPartOfMySelf([1,3],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(2,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[2,3,4,5,11,12,13]],1e-12))
        #
        a,b=f.buildSubMeshDataRange(2,5,1)
        self.assertTrue(m.buildPartOfMySelf([2,3,4],True).isEqual(a,1e-12))
        self.assertEqual(b,slice(6,16,1))
        ff=f[2:]
        ff.checkCoherency()
        self.assertTrue((m.buildPartOfMySelf([2,3,4],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[6,7,8,9,10,11,12,13,14,15]],1e-12))
        #
        ff=f[-2:0:-1]
        ff.checkCoherency()
        self.assertTrue((m.buildPartOfMySelf([3,2,1],True)).isEqual(ff.getMesh(),1e-12))
        self.assertTrue(9,ff.getMesh().getNumberOfNodes())
        self.assertTrue(3,ff.getMesh().getNumberOfCells())
        self.assertTrue(ff.getArray().isEqual(arr[[11,12,13,6,7,8,9,10,2,3,4,5]],1e-12))
        self.assertTrue(f[-2:0:-1,0].getArray().isEqual(arr[[11,12,13,6,7,8,9,10,2,3,4,5],0],1e-12))
        pass

    def testSwig2FieldDoubleApplyFuncBug1(self):
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setMesh(MEDCouplingDataForTest.build2DTargetMesh_1())
        f.applyFunc(3,700.)
        f.checkCoherency()
        self.assertEqual(3,f.getArray().getNumberOfComponents())
        f.getArray().rearrange(1)
        self.assertTrue(f.getArray().isUniform(700.,1e-10))
        f.getArray().rearrange(3)
        f.checkCoherency()
        f.applyFunc(4,800.)
        f.checkCoherency()
        self.assertEqual(4,f.getArray().getNumberOfComponents())
        f.getArray().rearrange(1)
        self.assertTrue(f.getArray().isUniform(800.,1e-10))
        f.getArray().rearrange(4)
        f.checkCoherency()
        pass

    def testSwig2ComputeTupleIdsNearTupleBug1(self):
        coords=[1.1,0.0, 1.1,0.0 ];
        coordsArr=DataArrayDouble(coords,2,2);
        mesh=MEDCouplingUMesh();
        mesh.setCoords(coordsArr);
        points=[1.1, 0.002]
        c,cI=mesh.getNodeIdsNearPoints(points,0.00185);
        self.assertTrue(c.isEqual(DataArrayInt([])))
        self.assertTrue(cI.isEqual(DataArrayInt([0,0])))
        c,cI=mesh.getNodeIdsNearPoints(points,0.00200000000000001);
        self.assertTrue(c.isEqual(DataArrayInt([0,1])))
        self.assertTrue(cI.isEqual(DataArrayInt([0,2])))
        pass

    def testSwig2NonRegressionBugChangeUnderlyingWithZeroCells(self):
        coords1=[0.,1.,2.,3.]
        coords2=[2.,1.,0.,3.] #0 <==> #2
        # mesh 1
        mesh1=MEDCouplingUMesh.New();
        coordsArr=DataArrayDouble.New(coords1,4,1);
        mesh1.setCoords(coordsArr);
        mesh1.setMeshDimension(0);
        mesh1.allocateCells(0);
        mesh1.finishInsertingCells();
        # mesh 2
        mesh2=mesh1.deepCpy();
        coordsArr=DataArrayDouble.New(coords2,4,1);
        mesh2.setCoords(coordsArr);
        field = mesh1.fillFromAnalytic(ON_NODES,1,"x")
        field.checkCoherency()
        levOfCheck = 10
        field.changeUnderlyingMesh( mesh2, levOfCheck, 1e-13, 0 )
        self.assertTrue( field.getArray().getValues() == coords2 )
        pass

    def testSwig2UMeshDistanceToMesh2(self):
        sz=5
        m=MEDCouplingCMesh()
        arr=DataArrayDouble(sz+1) ; arr.iota() ; arr/=sz
        m.setCoords(arr,arr,arr)
        m=m.buildUnstructured()
        m1=m.computeSkin()
        m1.zipCoords()
        c=m1.getCoords()[:]
        d=2*(c-[0.5,0.5,0.5])+[0.5,0.5,0.5]
        time_deb = datetime.now()
        #print "go.."
        a,b=m1.distanceToPoints(d)
        #print 'time spent in distanceToPoints %s ' %str(datetime.now() - time_deb)
        time_deb = datetime.now()
        a1=DataArrayDouble(len(d))
        b1=DataArrayInt(len(d))
        m1s=[m1[i] for i in xrange(m1.getNumberOfCells())]
        for j,pt in enumerate(d):
            eter=1e308
            fter=-1
            for i,miter in enumerate(m1s):
                e,f=miter.distanceToPoint(pt)
                self.assertEqual(0,f)
                if e<eter:
                    eter=e ; fter=i
                    pass
                pass
            a1[j]=eter
            b1[j]=fter
            pass
        #print 'time spent in naive distanceToPoints  %s ' %str(datetime.now() - time_deb)
        self.assertTrue(a.isEqual(a1,1e-12))
        self.assertTrue(b.isEqual(b1))
        self.assertTrue(a.isEqual(DataArrayDouble([0.8660254037844386,0.714142842854285,0.7071067811865476,0.7071067811865476,0.714142842854285,0.8660254037844386,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706632,0.714142842854285,0.7071067811865475,0.5099019513592785,0.5,0.5,0.5099019513592785,0.7071067811865476,0.7071067811865475,0.5099019513592785,0.5,0.5,0.5099019513592785,0.7071067811865476,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706632,0.714142842854285,0.8660254037844386,0.714142842854285,0.7071067811865476,0.7071067811865476,0.714142842854285,0.8660254037844386,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706631,0.714142842854285,0.5196152422706631,0.5196152422706632,0.5099019513592784,0.5099019513592785,0.5099019513592784,0.5099019513592785,0.5196152422706631,0.5196152422706632,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706632,0.714142842854285,0.7071067811865475,0.5099019513592785,0.5,0.5,0.5099019513592784,0.7071067811865475,0.5099019513592784,0.5099019513592785,0.5,0.5,0.5,0.5,0.5099019513592785,0.5099019513592785,0.7071067811865476,0.5099019513592785,0.5,0.5,0.5099019513592785,0.7071067811865476,0.7071067811865475,0.5099019513592785,0.5,0.5,0.5099019513592784,0.7071067811865475,0.5099019513592784,0.5099019513592785,0.5,0.5,0.5,0.5,0.5099019513592785,0.5099019513592785,0.7071067811865476,0.5099019513592785,0.5,0.5,0.5099019513592785,0.7071067811865476,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706631,0.714142842854285,0.5196152422706631,0.5196152422706632,0.5099019513592784,0.5099019513592785,0.5099019513592784,0.5099019513592785,0.5196152422706631,0.5196152422706632,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706632,0.714142842854285,0.8660254037844386,0.714142842854285,0.7071067811865476,0.7071067811865476,0.714142842854285,0.8660254037844386,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706632,0.714142842854285,0.7071067811865475,0.5099019513592785,0.5,0.5,0.5099019513592785,0.7071067811865476,0.7071067811865475,0.5099019513592785,0.5,0.5,0.5099019513592785,0.7071067811865476,0.714142842854285,0.5196152422706632,0.5099019513592785,0.5099019513592785,0.5196152422706632,0.714142842854285,0.8660254037844386,0.714142842854285,0.7071067811865476,0.7071067811865476,0.714142842854285,0.8660254037844386]),1e-12))
        self.assertTrue(b.isEqual(DataArrayInt([0,0,3,7,9,9,0,0,3,7,9,9,12,12,14,16,17,17,26,26,28,30,31,31,33,33,36,40,42,42,33,33,36,40,42,42,0,0,3,7,11,9,0,9,12,17,26,31,33,42,33,33,36,40,42,42,45,45,47,49,51,51,45,50,52,53,56,57,58,63,58,58,60,62,63,63,85,85,87,89,91,91,85,90,92,93,96,97,98,103,98,98,100,102,103,103,105,105,108,112,116,114,105,114,117,122,131,136,138,147,138,138,141,145,147,147,105,105,108,112,114,114,105,105,108,112,114,114,117,117,119,121,122,122,131,131,133,135,136,136,138,138,141,145,147,147,138,138,141,145,147,147])))
        pass

    def testSwig2NonRegressionBugDistance1(self):
        pt=DataArrayDouble([(8.8452994616207476,3.1547005383792515,3.1547005383792515)])
        coo=DataArrayDouble([(8,0,0),(8,0,8),(8,8,8),(8,8,0),(16,0,0),(16,0,8),(16,8,8),(16,8,0),(8,0,4),(8,4,8),(8,8,4),(8,4,0),(16,0,4),(16,4,8),(16,8,4),(16,4,0),(12,0,0),(12,0,8),(12,8,8),(12,8,0),(8,4,4),(16,4,4),(12,0,4),(12,4,8),(12,8,4),(12,4,0)])
        conn=DataArrayInt([4,15,21,12,4,16,25,15,12,22,16,4,0,8,20,11,16,0,11,25,22,8,0,16,15,7,14,21,15,25,19,7,7,19,24,14,11,20,10,3,25,11,3,19,19,3,10,24,12,21,13,5,13,23,17,5,5,17,22,12,8,1,9,20,23,9,1,17,17,1,8,22,21,14,6,13,14,24,18,6 ,6,18,23,13,20,9,2,10,24,10,2,18,18,2,9,23])
        m=MEDCouplingUMesh("mesh",2)
        m.setCoords(coo)
        m.allocateCells()
        for i in xrange(24):
            m.insertNextCell(NORM_QUAD4,conn[4*i:4*i+4])
            pass
        m.checkCoherency2()
        m0=m[3] ; m0.zipCoords()
        expectedDist=0.8452994616207476
        a,b=m0.distanceToPoint(pt)
        self.assertAlmostEqual(expectedDist,a,14)
        self.assertEqual(0,b)
        #
        a,b=m.distanceToPoint(pt)
        self.assertAlmostEqual(expectedDist,a,14)
        self.assertEqual(3,b)
        #
        fd=MEDCouplingFieldDiscretization.New(ON_CELLS)
        self.assertEqual(24,fd.getNumberOfTuples(m))
        fd=MEDCouplingFieldDiscretization.New(ON_NODES)
        self.assertEqual(26,fd.getNumberOfTuples(m))
        pass

    def testSwig2AreaBarySeg3Quad8Tri6QPolyg(self):
        #QUAD8 representing a circle of center zeBary and radius zeRadius
        zeBary=[5,6]
        zeRadius=3
        d=DataArrayDouble(8,2)
        d[:,0]=zeRadius
        d[:,1]=[87,-100,-170,110,5,-130,175,95] # angle in degree
        d[:,1]*=pi/180. # angle in radian
        d=d.fromPolarToCart()
        d+=zeBary
        m=MEDCouplingUMesh("quad8",2) ; m.allocateCells() ; m.insertNextCell(NORM_QUAD8,range(8)) ; m.setCoords(d)
        self.assertTrue(m.getBarycenterAndOwner().isEqual(DataArrayDouble(zeBary,1,2),1e-13))
        self.assertAlmostEqual(float(m.getMeasureField(False).getArray()),pi*zeRadius*zeRadius,12)
        tri32D=m.buildDescendingConnectivity()[0][0] ; tri32D.zipCoords()
        # spaceDim=3 QUAD8 becomes QUAD4 ... for the moment
        m.setCoords(m.getCoords().changeNbOfComponents(3,0.))
        m2=m.deepCpy()
        m2.convertQuadraticCellsToLinear()
        self.assertAlmostEqual(float(m.getMeasureField(False).getArray()),float(m2.getMeasureField(False).getArray()),12)
        self.assertTrue(m.getBarycenterAndOwner().isEqual(m2.getBarycenterAndOwner(),1e-13))
        #TRI6 representing a circle of center zeBary and radius zeRadius
        zeBary=[5,6]
        zeRadius=3
        d=DataArrayDouble(6,2)
        d[:,0]=zeRadius
        d[:,1]=[87,-100,110,5,175,95] # angle in degree
        d[:,1]*=pi/180. # angle in radian
        d=d.fromPolarToCart()
        d+=zeBary
        m=MEDCouplingUMesh("tri6",2) ; m.allocateCells() ; m.insertNextCell(NORM_TRI6,range(6)) ; m.setCoords(d)
        self.assertTrue(m.getBarycenterAndOwner().isEqual(DataArrayDouble(zeBary,1,2),1e-13))
        self.assertAlmostEqual(float(m.getMeasureField(False).getArray()),pi*zeRadius*zeRadius,12)
        # spaceDim=3 TRI6 becomes TRI3 ... for the moment
        m.setCoords(m.getCoords().changeNbOfComponents(3,0.))
        m2=m.deepCpy()
        m2.convertQuadraticCellsToLinear()
        self.assertAlmostEqual(float(m.getMeasureField(False).getArray()),float(m2.getMeasureField(False).getArray()),12)
        self.assertTrue(m.getBarycenterAndOwner().isEqual(m2.getBarycenterAndOwner(),1e-13))
        # QPOLYG representing a circle of center zeBary and radius zeRadius
        zeBary=[5,6]
        zeRadius=3
        d=DataArrayDouble(10,2)
        d[:,0]=zeRadius
        d[:,1]=[87,-80,-100,-170,110,5,-90,-130,175,95] # angle in degree
        d[:,1]*=pi/180. # angle in radian
        d=d.fromPolarToCart()
        d+=zeBary
        m=MEDCouplingUMesh("qpolyg",2) ; m.allocateCells() ; m.insertNextCell(NORM_QPOLYG,range(10)) ; m.setCoords(d)
        self.assertTrue(m.getBarycenterAndOwner().isEqual(DataArrayDouble(zeBary,1,2),1e-13))
        self.assertAlmostEqual(float(m.getMeasureField(False).getArray()),pi*zeRadius*zeRadius,12)
        # spaceDim=3 QPOLYG becomes POLYG ... for the moment
        m.setCoords(m.getCoords().changeNbOfComponents(3,0.))
        m2=m.deepCpy()
        m2.convertQuadraticCellsToLinear() ; m2.checkCoherency2()
        self.assertTrue(m2.getAllGeoTypes()==[NORM_POLYGON] and m2.getNodalConnectivity().getValues()==[5,0,1,2,3,4])
        self.assertAlmostEqual(float(m.getMeasureField(False).getArray()),float(m2.getMeasureField(False).getArray()),12)
        self.assertTrue(m.getBarycenterAndOwner().isEqual(m2.getBarycenterAndOwner(),1e-13))
        # TRI3
        self.assertAlmostEqual(float(tri32D.getMeasureField(False).getArray()),(87+100)*pi/180*zeRadius,13)
        exp=DataArrayDouble(1,2) ; exp[:,0]=3 ; exp[:,1]=(87-100)/2. ; exp[:,1]*=pi/180. ;  exp=exp.fromPolarToCart() ; exp+=DataArrayDouble([5,6],1,2)
        self.assertTrue(tri32D.getBarycenterAndOwner().isEqual(exp,1e-12))
        # spaceDim=3 TRI3 becomes TRI2 ... for the moment
        tri32D.changeSpaceDimension(3)
        tri2=tri32D.deepCpy() ; tri2.convertQuadraticCellsToLinear()
        self.assertAlmostEqual(float(tri32D.getMeasureField(False).getArray()),float(tri2.getMeasureField(False).getArray()),13)
        self.assertTrue(tri32D.getBarycenterAndOwner().isEqual(tri2.getBarycenterAndOwner(),1e-12))
        tri32D.changeSpaceDimension(1)
        self.assertAlmostEqual(float(tri32D.getMeasureField(False).getArray()),-0.67795240172962323,12)
        pass

    # this bug 5/6/2013 is swig specific
    def testSwigNonRegressionBugRotate3D1(self):
        m=MEDCouplingUMesh.New()
        dataArray=DataArrayDouble.New(100,3)
        dataArray[:]=0.
        dataArray[0]=[0.,1,3]
        m.setCoords(dataArray[0])
        m1=m.deepCpy()
        m.rotate([0.,0.,3.],[1.,0.,0.],0.5*pi)
        self.assertTrue(m.getCoords().isEqual(DataArrayDouble([0.,0.,4.],1,3),1e-15))
        #
        d1=DataArrayDouble([0.,0.,3.],1,3) ; d2=DataArrayDouble([1.,0.,0.],1,3)
        pts=[[0.,0.,3.],[(0.,0.,3.)],DataArrayDouble([0.,0.,3.],1,3),list(d1)[0]]
        vec=[[1.,0.,0.],[(1.,0.,0.)],DataArrayDouble([1.,0.,0.],1,3),list(d2)[0]]
        for p in pts:
            for v in vec:
                m2=m1.deepCpy()
                m2.rotate(p,v,0.5*pi)
                self.assertTrue(m2.getCoords().isEqual(DataArrayDouble([0.,0.,4.],1,3),1e-15))
                pass
        pass

    def testSwig2DataArrayCount1(self):
        d=DataArrayInt([])
        self.assertEqual(0,d.getNumberOfTuples())
        self.assertEqual(1,d.getNumberOfComponents())
        self.assertEqual(0,d.count(0))
        self.assertEqual(0,d.count(1))
        self.assertEqual(0,d.count(-1))
        d=DataArrayInt([2,1,-2,-3,2,0,0,7,2,-2,3,0])
        self.assertEqual(12,d.getNumberOfTuples())
        self.assertEqual(1,d.getNumberOfComponents())
        self.assertEqual(3,d.count(0))
        self.assertEqual(1,d.count(1))
        self.assertEqual(0,d.count(-1))
        self.assertEqual(2,d.count(-2))
        self.assertEqual(3,d.count(2))
        e=d.getDifferentValues()
        f=DataArrayInt()
        for it in e:
            f.pushBackSilent(d.count(int(it)))
            pass
        self.assertEqual(12,f.accumulate()[0])
        #
        eps=1e-12
        d=DataArrayDouble([])
        self.assertEqual(0,d.getNumberOfTuples())
        self.assertEqual(1,d.getNumberOfComponents())
        self.assertEqual(0,d.count(0,eps))
        self.assertEqual(0,d.count(1,eps))
        self.assertEqual(0,d.count(-1,eps))
        d=DataArrayDouble([2,1,-2,-3,2,0,eps/10,7,2+eps/10,-2,3,0])
        self.assertEqual(12,d.getNumberOfTuples())
        self.assertEqual(1,d.getNumberOfComponents())
        self.assertEqual(3,d.count(0,eps))
        self.assertEqual(1,d.count(1,eps))
        self.assertEqual(0,d.count(-1,eps))
        self.assertEqual(2,d.count(-2,eps))
        self.assertEqual(3,d.count(2,eps))
        self.assertEqual(3,d.count(2,eps))
        self.assertEqual(2,d.count(2,eps/100))
        e=d.getDifferentValues(eps)
        f=DataArrayInt()
        for it in e:
            f.pushBackSilent(d.count(float(it),eps))
            pass
        self.assertEqual(12,f.accumulate()[0])
        pass

    def testSwig2DataArrayGetSlice1(self):
        s=slice(2,18,1)
        self.assertEqual(DataArray.GetNumberOfItemGivenBESRelative(s),16)
        self.assertEqual(DataArray.GetNumberOfItemGivenBES(s),16)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(2,6,1))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(6,10,1))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(10,14,1))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(14,18,1))
        #
        s=slice(2,18,2)
        self.assertEqual(DataArray.GetNumberOfItemGivenBESRelative(s),8)
        self.assertEqual(DataArray.GetNumberOfItemGivenBES(s),8)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(2,6,2))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(6,10,2))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(10,14,2))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(14,18,2))
        #
        s=slice(1,18,1)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(1,5,1))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(5,9,1))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(9,13,1))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(13,18,1))# 18 not 17
        #
        s=slice(1,18,2)
        self.assertEqual(DataArray.GetNumberOfItemGivenBESRelative(s),9)
        self.assertEqual(DataArray.GetNumberOfItemGivenBES(s),9)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(1,5,2))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(5,9,2))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(9,13,2))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(13,18,2))# 18 not 17
        #
        s=slice(18,2,-1)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(18,14,-1))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(14,10,-1))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(10,6,-1))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(6,2,-1))
        #
        s=slice(18,2,-2)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(18,14,-2))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(14,10,-2))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(10,6,-2))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(6,2,-2))
        #
        s=slice(18,1,-1)
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(18,14,-1))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(14,10,-1))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(10,6,-1))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(6,1,-1))# 1 not 2
        #
        s=slice(18,1,-2)
        self.assertEqual(DataArray.GetNumberOfItemGivenBESRelative(s),9)
        self.assertRaises(InterpKernelException,DataArray.GetNumberOfItemGivenBES,s)
        self.assertEqual(sum([DataArray.GetNumberOfItemGivenBESRelative(DataArray.GetSlice(s,i,4)) for i in xrange(4)]),DataArray.GetNumberOfItemGivenBESRelative(s))
        self.assertEqual(DataArray.GetSlice(s,0,4),slice(18,14,-2))
        self.assertEqual(DataArray.GetSlice(s,1,4),slice(14,10,-2))
        self.assertEqual(DataArray.GetSlice(s,2,4),slice(10,6,-2))
        self.assertEqual(DataArray.GetSlice(s,3,4),slice(6,1,-2))# 1 not 2
        self.assertRaises(InterpKernelException,DataArray.GetSlice,slice(0,None,2),0,4)
        #
        d=DataArrayInt.Range(0,18,1)
        s=slice(2,None,1)
        self.assertEqual(d.getNumberOfItemGivenBES(s),16)
        self.assertEqual(d.getNumberOfItemGivenBESRelative(s),16)
        self.assertEqual(d.getSlice(s,0,4),slice(2,6,1))
        self.assertEqual(d.getSlice(s,1,4),slice(6,10,1))
        self.assertEqual(d.getSlice(s,2,4),slice(10,14,1))
        self.assertEqual(d.getSlice(s,3,4),slice(14,18,1))
        #
        d=DataArrayInt.Range(0,18,1)
        s=slice(2,-2,1)
        self.assertEqual(d.getSlice(s,0,4),slice(2,5,1))
        self.assertEqual(d.getSlice(s,1,4),slice(5,8,1))
        self.assertEqual(d.getSlice(s,2,4),slice(8,11,1))
        self.assertEqual(d.getSlice(s,3,4),slice(11,16,1))
        #
        d=DataArrayInt.Range(0,18,1)
        s=slice(None,None,1)
        self.assertEqual(d.getSlice(s,0,4),slice(0,4,1))
        self.assertEqual(d.getSlice(s,1,4),slice(4,8,1))
        self.assertEqual(d.getSlice(s,2,4),slice(8,12,1))
        self.assertEqual(d.getSlice(s,3,4),slice(12,18,1))
        #
        d=DataArrayInt.Range(0,18,1)
        s=slice(None,2,-2)
        self.assertRaises(InterpKernelException,d.getNumberOfItemGivenBES,s)
        self.assertEqual(d.getNumberOfItemGivenBESRelative(s),8)
        self.assertEqual(d.getSlice(s,0,4),slice(17,13,-2))
        self.assertEqual(d.getSlice(s,1,4),slice(13,9,-2))
        self.assertEqual(d.getSlice(s,2,4),slice(9,5,-2))
        self.assertEqual(d.getSlice(s,3,4),slice(5,2,-2))
        pass

    def testSwig2AccumulatePerChunk1(self):
        arr=DataArrayDouble(11) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
        m=m.buildUnstructured()
        m0=m[::2] ; ids0=m0.simplexize(0) ; m1=m[1::2]
        m=MEDCouplingUMesh.MergeUMeshesOnSameCoords(m0,m1) ; m.setName("mesh")
        m.checkConsecutiveCellTypesForMEDFileFrmt()
        #
        formula="7-sqrt((x-5.)*(x-5.)+(y-5.)*(y-5.))"
        f=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f.setMesh(m)
        f.fillFromAnalytic(1,formula)
        f.setName("Field1") ; f.setTime(1.1,1,-1)
        f.checkCoherency()
        #
        arr=f.getArray()
        arr2=DataArrayDouble(len(arr),2) ; arr2[:,0]=arr
        arr2=DataArrayDouble(len(arr),2) ; arr2[:,0]=arr ; arr2[:,1]=2*arr
        f.setArray(arr2)
        f.checkCoherency()
        # here the compact code to obviously put field on cell to nodes
        rn,rni=f.getMesh().getReverseNodalConnectivity()
        arr2=f.getArray()[rn]
        arr4=arr2.accumulatePerChunck(rni)
        nbOfCellsSharingNodes=rni.deltaShiftIndex()
        arr4/=nbOfCellsSharingNodes.convertToDblArr()
        #
        maxNbCSN=nbOfCellsSharingNodes.getMaxValue()[0]
        arr3=DataArrayDouble(f.getMesh().getNumberOfNodes(),f.getArray().getNumberOfComponents()) ; arr3[:]=0.
        for i in xrange(1,maxNbCSN+1):
            ids=nbOfCellsSharingNodes.getIdsEqual(i)
            if len(ids)==0:
                continue
            for j in range(i):
                rni2=rni[ids] ; rni2+=j
                arr3[ids]+=arr2[rni2]
                pass
            arr3[ids]/=i
            pass
        fNode=MEDCouplingFieldDouble(ON_NODES,ONE_TIME) ; fNode.setMesh(m)
        fNode.setName("Field1Node") ; fNode.setTime(1.1,1,-1)
        fNode.setArray(arr3) ; fNode.checkCoherency()
        self.assertTrue(arr3.isEqual(arr4,1e-12))
        #
        d=DataArrayInt.Range(0,20,1)
        self.assertTrue(d.accumulatePerChunck([2,4,12]).isEqual(DataArrayInt([5,60])))
        #
        a=DataArrayDouble(12) ; a.iota() ; a.rearrange(3)
        b=DataArrayDouble(12) ; b.iota(20) ; b.rearrange(3)
        ids=DataArrayInt([])
        self.assertEqual(len(a[ids]),0)
        self.assertEqual(len(b[ids]),0)
        a2=a.deepCpy() ;  a2[ids]+=b[ids] ; self.assertTrue(a2.isEqual(a,1e-15))
        a2=a.deepCpy() ;  a2[ids]*=b[ids] ; self.assertTrue(a2.isEqual(a,1e-15))
        a2=a.deepCpy() ;  a2[ids]/=b[ids] ; self.assertTrue(a2.isEqual(a,1e-15))
        a2=a.deepCpy() ;  a2[ids]-=b[ids] ; self.assertTrue(a2.isEqual(a,1e-15))
        pass

    def testSwig2CheckAndPreparePermutation1(self):
        a=DataArrayInt([10003,9999999,5,67])
        self.assertTrue(a.checkAndPreparePermutation().isEqual(DataArrayInt([2,3,0,1])))
        a=DataArrayInt([10003,-9999999,5,67])
        self.assertTrue(a.checkAndPreparePermutation().isEqual(DataArrayInt([3,0,1,2])))
        a=DataArrayInt([])
        self.assertTrue(a.checkAndPreparePermutation().isEqual(DataArrayInt([])))
        a=DataArrayInt([])
        a.iota();
        self.assertTrue(a.isEqual(DataArrayInt([])))
        pass

    def testSwig21SGTUMesh1(self):
        m=MEDCoupling1GTUMesh.New("m",NORM_PENTA6)
        m.__repr__() ; m.__str__()
        self.assertTrue(isinstance(m,MEDCoupling1SGTUMesh))
        m.setCoords(DataArrayDouble(20,3))
        m.allocateCells()
        m.__repr__() ; m.__str__()
        m.insertNextCell([0,1,2,5,7,2])
        self.assertEqual(1,m.getNumberOfCells())
        self.assertTrue(DataArrayInt([6]).isEqual(m.computeNbOfNodesPerCell()))
        self.assertTrue(DataArrayInt([5]).isEqual(m.computeNbOfFacesPerCell()))
        m.__repr__() ; m.__str__()
        m.checkCoherency()
        m.checkCoherency2()
        #
        cm=MEDCouplingCMesh() ; cm.setName("m")
        arr0=DataArrayDouble(6) ; arr0.iota()
        arr1=DataArrayDouble([0,1])
        cm.setCoords(arr0,arr1,arr1) ; um=cm.buildUnstructured()
        #
        m=MEDCoupling1SGTUMesh("m",NORM_QUAD4)
        mem_m=m.getHeapMemorySize()
        m.allocateCells(5)
        self.assertIn(m.getHeapMemorySize()-mem_m,xrange(5*4*4,5*4*4+32))
        self.assertEqual(m.getNodalConnectivity().getNbOfElemAllocated(),20)
        m.setCoords(um.getCoords())
        m.insertNextCell([1,0,6,7])
        self.assertEqual(1,m.getNumberOfCells())
        m.insertNextCell([2,1,7,8])
        m.insertNextCell([3,2,8,9])
        m.insertNextCell([4,3,9,10])
        m.insertNextCell([5,4,10,11])
        self.assertEqual(5,m.getNumberOfCells())
        self.assertRaises(InterpKernelException,m.insertNextCell,[0,6,7])
        self.assertRaises(InterpKernelException,m.insertNextCell,[0,6,7,1,2])
        self.assertEqual(m.getNodalConnectivity().getNbOfElemAllocated(),20)
        f=m.getMeasureField(ON_CELLS)
        self.assertEqual(f.getMesh().getHiddenCppPointer(),m.getHiddenCppPointer())
        self.assertTrue(f.getArray().isUniform(1,1e-14))
        self.assertEqual(m.getType(),10)
        self.assertEqual(m.getCellModelEnum(),NORM_QUAD4)
        mo=MEDCoupling1SGTUMesh("m",NORM_QUAD4) ; mo.setCoords(m.getCoords())
        mo.setNodalConnectivity(DataArrayInt([1,0,6,7,2,1,7,8,3,2,8,9,4,3,9,10,5,4,10,11]))
        self.assertTrue(m.isEqual(mo,1e-12))
        #
        mo2=MEDCoupling1SGTUMesh.Merge1SGTUMeshesOnSameCoords([m[[0,1]],m[[2]],m[[3,4]]])
        mo2.setName(m.getName())
        self.assertTrue(m.isEqual(mo2,1e-12))
        #
        mp0=m[[0]] ; mp0.zipCoords() ; mp1=m[2] ; mp1.zipCoords() ; mp2=m[4] ; mp2.zipCoords()
        mo3=MEDCoupling1SGTUMesh.Merge1SGTUMeshes([mp0,mp1,mp2])
        self.assertTrue(isinstance(mo3,MEDCoupling1SGTUMesh))
        mo3.setName(m.getName())
        m_ref=m[(0,2,4)] ; m_ref.zipCoords()
        m_ref.tryToShareSameCoordsPermute(mo3,1e-12)
        self.assertTrue(m_ref.isEqual(mo3,1e-12))
        #
        m1=um.buildDescendingConnectivity()[0]
        ids=m1.getCellIdsFullyIncludedInNodeIds(DataArrayInt.Range(0,12,1))
        m1=m1[ids]
        m1c=m1.convertIntoSingleGeoTypeMesh()
        self.assertTrue(isinstance(m1c,MEDCoupling1SGTUMesh))
        self.assertEqual(m1c.getCoords().getHiddenCppPointer(),m.getCoords().getHiddenCppPointer())
        m1c.checkCoherency2()
        self.assertTrue(m1c.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,2,1,7,8,3,2,8,9,4,3,9,10,5,4,10,11])))
        self.assertEqual(20,m1c.getNodalConnectivityLength())
        self.assertTrue(m.isEqual(m1c,1e-12))
        m.getNodalConnectivity().setIJ(1,0,1)
        self.assertTrue(not m.isEqual(m1c,1e-12))
        m.getNodalConnectivity().setIJ(1,0,0)
        self.assertTrue(m.isEqual(m1c,1e-12))
        m1c.setCoords(m.getCoords().deepCpy())
        self.assertTrue(m.isEqual(m1c,1e-12))
        m1c.getCoords().setIJ(0,1,0.1)
        self.assertTrue(not m.isEqual(m1c,1e-12))
        m1c.getCoords().setIJ(0,1,0)
        self.assertTrue(m.isEqual(m1c,1e-12))
        m1c.getCoords().setInfoOnComponent(1,"X")
        self.assertTrue(not m.isEqual(m1c,1e-12) and m.isEqualWithoutConsideringStr(m1c,1e-12))
        m.getCoords().setInfoOnComponent(1,"X")
        self.assertTrue(m.isEqual(m1c,1e-12) and m.isEqualWithoutConsideringStr(m1c,1e-12))
        m.setName("m2")
        self.assertTrue(not m.isEqual(m1c,1e-12) and m.isEqualWithoutConsideringStr(m1c,1e-12))
        #
        m.checkCoherency() ; m.checkCoherency1() ; m.checkCoherency2()
        self.assertEqual(m.getMeshDimension(),2)
        self.assertTrue(m.giveCellsWithType(NORM_QUAD4).isEqual(DataArrayInt([0,1,2,3,4])))
        self.assertTrue(m.giveCellsWithType(NORM_TRI3).isEqual(DataArrayInt([])))
        self.assertEqual(m.getNumberOfCellsWithType(NORM_QUAD4),5)
        self.assertEqual(m.getNumberOfCellsWithType(NORM_TRI3),0)
        self.assertEqual(m.getTypeOfCell(3),NORM_QUAD4)
        self.assertRaises(InterpKernelException,m.getTypeOfCell,5)
        self.assertEqual(m.getAllGeoTypes(),[NORM_QUAD4])
        self.assertEqual(m.getDistributionOfTypes(),[[NORM_QUAD4,5,-1]])
        ##
        pfl1=DataArrayInt([1,3,4])
        a,b,c=m.splitProfilePerType(pfl1)
        d,e,f=m.buildUnstructured().splitProfilePerType(pfl1)
        self.assertTrue(a==[[4,3,0]] and len(b)==1 and b[0].isEqual(DataArrayInt([0,1,2])) and len(c)==1 and c[0].getHiddenCppPointer()==pfl1.getHiddenCppPointer())
        self.assertTrue(a==d and len(b)==1 and b[0].isEqual(e[0]) and len(c)==1 and c[0].isEqual(f[0]))
        #
        pfl2=DataArrayInt([0,1,2,3])
        a,b,c=m.splitProfilePerType(pfl2)
        d,e,f=m.buildUnstructured().splitProfilePerType(pfl2)
        self.assertTrue(a==[[4,4,0]] and len(b)==1 and b[0].isEqual(DataArrayInt([0,1,2,3])) and len(c)==1 and c[0].getHiddenCppPointer()==pfl2.getHiddenCppPointer())
        self.assertTrue(a==d and len(b)==1 and b[0].isEqual(e[0]) and len(c)==1 and c[0].isEqual(f[0]))
        #
        pfl3=DataArrayInt([0,1,2,3,4])
        a,b,c=m.splitProfilePerType(pfl3)
        d,e,f=m.buildUnstructured().splitProfilePerType(pfl3)
        self.assertTrue(a==[[4,5,-1]] and len(b)==1 and b[0].isEqual(DataArrayInt([0,1,2,3,4])) and c==[])
        self.assertTrue(a==d and len(b)==1 and b[0].isEqual(e[0]) and c==[])
        #
        invalidPfl=DataArrayInt([1,2,3,4,5])
        self.assertRaises(InterpKernelException,m.splitProfilePerType,invalidPfl)
        self.assertRaises(InterpKernelException,m.buildUnstructured().splitProfilePerType,invalidPfl)
        ##
        pfl1=DataArrayInt([1,2,3])
        a=m.checkTypeConsistencyAndContig([NORM_QUAD4,3,0],[pfl1])
        b=m.buildUnstructured().checkTypeConsistencyAndContig([NORM_QUAD4,3,0],[pfl1])
        self.assertTrue(a.isEqual(b) and pfl1.getHiddenCppPointer(),a.getHiddenCppPointer())
        #
        pfl2=DataArrayInt([0,1,2,3])
        a=m.checkTypeConsistencyAndContig([NORM_QUAD4,4,0],[pfl2])
        b=m.buildUnstructured().checkTypeConsistencyAndContig([NORM_QUAD4,4,0],[pfl2])
        self.assertTrue(a.isEqual(b) and pfl2.getHiddenCppPointer()==a.getHiddenCppPointer())
        #
        pfl3=DataArrayInt([0,1,2,3,4])
        a=m.checkTypeConsistencyAndContig([NORM_QUAD4,4,0],[pfl3])
        b=m.buildUnstructured().checkTypeConsistencyAndContig([NORM_QUAD4,5,0],[pfl3])
        self.assertTrue(a.isEqual(b) and pfl3.getHiddenCppPointer()==a.getHiddenCppPointer())
        #
        invalidPfl=DataArrayInt([1,2,3,4,5])
        self.assertRaises(InterpKernelException,m.checkTypeConsistencyAndContig,[NORM_QUAD4,5,0],[invalidPfl])
        self.assertRaises(InterpKernelException,m.buildUnstructured().checkTypeConsistencyAndContig,[NORM_QUAD4,5,0],[invalidPfl])
        ##
        self.assertTrue(DataArrayInt([4,4,4,4,4]).isEqual(m.computeNbOfNodesPerCell()))
        ##
        self.assertEqual(m.getNodeIdsOfCell(1),[2,1,7,8])
        ##
        self.assertTrue(m.computeIsoBarycenterOfNodesPerCell().isEqual(DataArrayDouble([(0.5,0.5,0),(1.5,0.5,0),(2.5,0.5,0),(3.5,0.5,0),(4.5,0.5,0)]),1e-13))
        ##
        ref=m.getCoords().getHiddenCppPointer()
        mcpy=m.deepCpy() ; mcpy.insertNextCell([1,0,6,7])
        c=m.getNodalConnectivity().deepCpy()
        o2n=DataArrayInt([2,0,1,4,3])
        m.renumberCells(o2n,False)
        c.rearrange(4) ; c.renumberInPlace(o2n) ; c.rearrange(1)
        self.assertTrue(c.isEqual(m.getNodalConnectivity()))
        self.assertEqual(ref,m.getCoords().getHiddenCppPointer())
        m2=mcpy.mergeMyselfWith(m)
        self.assertTrue(isinstance(m2,MEDCoupling1SGTUMesh))
        self.assertEqual(11,m2.getNumberOfCells())
        self.assertEqual(48,m2.getNumberOfNodes())
        self.assertTrue(m2.getCoords().isEqual(DataArrayDouble.Aggregate([m.getCoords(),m.getCoords()]),1e-12))
        self.assertTrue(m2.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,2,1,7,8,3,2,8,9,4,3,9,10,5,4,10,11,1,0,6,7,26,25,31,32,27,26,32,33,25,24,30,31,29,28,34,35,28,27,33,34])))
        ##
        mu=m.buildUnstructured()
        mu.checkCoherency2()
        self.assertEqual(mu.getCoords().getHiddenCppPointer(),m.getCoords().getHiddenCppPointer())
        self.assertEqual(2,mu.getMeshDimension())
        self.assertEqual([NORM_QUAD4],mu.getAllGeoTypes())
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([4,2,1,7,8,4,3,2,8,9,4,1,0,6,7,4,5,4,10,11,4,4,3,9,10])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25])))
        ##
        for typ in [0,1]:
            mcpy2=m.deepCpy() ; umcpy2=mcpy2.buildUnstructured()
            ids=mcpy2.simplexize(typ) ; ids2=umcpy2.simplexize(typ)
            self.assertTrue(ids.isEqual(ids2))
            mcpy3=umcpy2.convertIntoSingleGeoTypeMesh()
            self.assertTrue(mcpy2.isEqual(mcpy3,1e-14))
            pass
        um1=um.convertIntoSingleGeoTypeMesh()
        self.assertEqual(8,um1.getNumberOfNodesPerCell())
        for typ in [PLANAR_FACE_5,PLANAR_FACE_6]:
            mcpy2=um1.deepCpy() ; umcpy2=mcpy2.buildUnstructured()
            ids=mcpy2.simplexize(typ) ; ids2=umcpy2.simplexize(typ)
            self.assertTrue(ids.isEqual(ids2))
            mcpy3=umcpy2.convertIntoSingleGeoTypeMesh()
            self.assertTrue(mcpy2.isEqual(mcpy3,1e-14))
            pass
        ##
        self.assertRaises(InterpKernelException,mcpy.mergeMyselfWithOnSameCoords,m)
        mcpy.tryToShareSameCoords(m,1e-14)
        m3=mcpy.mergeMyselfWithOnSameCoords(m)
        self.assertTrue(isinstance(m3,MEDCoupling1SGTUMesh))
        self.assertEqual(11,m3.getNumberOfCells())
        self.assertEqual(24,m3.getNumberOfNodes())
        self.assertEqual(m3.getCoords().getHiddenCppPointer(),mcpy.getCoords().getHiddenCppPointer())
        self.assertTrue(m3.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,2,1,7,8,3,2,8,9,4,3,9,10,5,4,10,11,1,0,6,7,2,1,7,8,3,2,8,9,1,0,6,7,5,4,10,11,4,3,9,10])))
        ##
        ref=mcpy.getCoords().deepCpy()
        c3=mcpy.getNodalConnectivity()[:]
        mcpy.getNodalConnectivity().setIJ(int(c3.getIdsEqual(11)),0,24)
        c2=DataArrayDouble.Aggregate([mcpy.getCoords(),mcpy.getCoords()[11:]])
        mcpy.setCoords(c2)
        mcpy.checkCoherency2()
        a,b=mcpy.getNodeIdsInUse()
        self.assertEqual(12,b)
        self.assertTrue(a.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])))
        ids=mcpy.zipCoordsTraducer()
        self.assertTrue(ids.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1])))
        self.assertTrue(mcpy.getCoords().isEqual(ref[:12],1e-12))
        self.assertTrue(mcpy.getNodalConnectivity().isEqual(c3))
        mcpy.checkCoherency2()
        ##
        m4=mcpy[DataArrayInt([0,3,4])]
        m5=mcpy.buildPartOfMySelfKeepCoords(DataArrayInt([0,3,4]))
        self.assertTrue(isinstance(m4,MEDCoupling1SGTUMesh))
        self.assertTrue(m4.isEqual(m5,-1e-14))# < 0 not a bug it proves that coordinates pointer are equal
        self.assertTrue(m4.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,4,3,9,10,5,4,10,11])))
        m6=mcpy[::2]
        self.assertTrue(isinstance(m6,MEDCoupling1SGTUMesh))
        self.assertTrue(m6.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,3,2,8,9,5,4,10,11])))
        ##
        mcpy.setCoords(DataArrayDouble.Aggregate([mcpy.getCoords(),mcpy.getCoords()]))
        mcpy.checkCoherency2()
        ##
        mcppy=mcpy.deepCpyConnectivityOnly()
        self.assertTrue(mcppy.isEqual(mcpy,1e-12))
        self.assertTrue(mcppy.getCoords().getHiddenCppPointer()==mcpy.getCoords().getHiddenCppPointer())
        self.assertTrue(mcppy.getNodalConnectivity().isEqual(mcpy.getNodalConnectivity()))
        self.assertTrue(mcppy.getNodalConnectivity().getHiddenCppPointer()!=mcpy.getNodalConnectivity().getHiddenCppPointer())
        ##
        a,b=mcpy.getReverseNodalConnectivity()
        self.assertTrue(a.isEqual(DataArrayInt([0,5,0,1,5,1,2,2,3,3,4,4,0,5,0,1,5,1,2,2,3,3,4,4])))
        self.assertTrue(b.isEqual(DataArrayInt([0,2,5,7,9,11,12,14,17,19,21,23,24,24,24,24,24,24,24,24,24,24,24,24,24])))
        self.assertTrue(mcpy.fillCellIdsToKeepFromNodeIds([0,1,6,7],False).isEqual(DataArrayInt([0,1,5])))
        self.assertTrue(mcpy.fillCellIdsToKeepFromNodeIds([0,1,6,7],True).isEqual(DataArrayInt([0,5])))
        self.assertTrue(mcpy.getCellsInBoundingBox([(0,1),(0,1),(0,1)],1e-12).isEqual(DataArrayInt([0,1,5])))
        f=mcpy.buildOrthogonalField()
        self.assertEqual(f.getMesh().getHiddenCppPointer(),mcpy.getHiddenCppPointer())
        self.assertTrue(f.getArray().isEqual(DataArrayDouble(6*[(0,0,-1)]),1e-12))
        mcpy.changeSpaceDimension(2)
        self.assertEqual(1,mcpy.getCellContainingPoint([1.5,0.5],1e-12))
        ##
        self.assertTrue(mcpy.fillCellIdsToKeepFromNodeIds(DataArrayInt([6,7]),False).isEqual(DataArrayInt([0,1,5])))
        ##
        mcpy2=mcpy.deepCpy()
        self.assertEqual([None,None],mcpy.checkGeoEquivalWith(mcpy2,1,1e-12))#fast equal
        mcpy.checkFastEquivalWith(mcpy2,1e-12)
        mcpy2.renumberCells([0,2,4,3,1,5])
        mcpy.checkFastEquivalWith(mcpy2,1e-12)
        self.assertEqual([None,None],mcpy.checkGeoEquivalWith(mcpy2,1,1e-12))#fast equal
        mcpy2.renumberCells([0,2,4,3,1,5])
        mcpy2.renumberCells([1,3,5,0,2,4])
        self.assertRaises(InterpKernelException,mcpy.checkFastEquivalWith,mcpy2,1e-12)
        self.assertRaises(InterpKernelException,mcpy.checkGeoEquivalWith,mcpy2,1,1e-12)#fast equal
        pass

    def testSwig21DGTUMesh1(self):
        a0=DataArrayInt([0,2,3,5,6,8])
        a1=DataArrayInt([0,4,7,11,14,18,21,25])
        a2=DataArrayInt([0,1,4,5])
        self.assertTrue(DataArrayInt.AggregateIndexes([a0,a1,a2]).isEqual(DataArrayInt([0,2,3,5,6,8,12,15,19,22,26,29,33,34,37,38])))
        self.assertEqual(a1[3:].front(),11)
        self.assertEqual(a1[4:].convertToDblArr().front(),14.)
        a1c=DataArrayInt([5,7,1,2, 8,11,0, 5,6,3,12, 1,5,2, 13,12,11,7, 6,1,0, 20,21,19,17])
        d,e=MEDCouplingUMesh.ExtractFromIndexedArrays2(1,5,2,a1c,a1)
        self.assertTrue(d.isEqual(DataArrayInt([8,11,0,1,5,2])))
        self.assertTrue(e.isEqual(DataArrayInt([0,3,6])))
        #
        m=MEDCouplingDataForTest.build2DTargetMesh_1()[0,3,4]
        ref=DataArrayInt([0,3,4,1,6,7,4,3,7,8,5,4])
        self.assertTrue(m.convertNodalConnectivityToStaticGeoTypeMesh().isEqual(ref))
        d,e=m.convertNodalConnectivityToDynamicGeoTypeMesh()
        self.assertTrue(d.isEqual(ref))
        self.assertTrue(e.isEqual(DataArrayInt.Range(0,13,4)))
        self.assertTrue(m.fillCellIdsToKeepFromNodeIds(DataArrayInt([6,7]),False).isEqual(DataArrayInt([1,2])))
        #
        m=MEDCoupling1GTUMesh.New("m",NORM_POLYHED)
        self.assertTrue(isinstance(m,MEDCoupling1DGTUMesh))
        m.__repr__() ; m.__str__()
        m.setCoords(DataArrayDouble(20,3))
        m.allocateCells()
        m.__repr__() ; m.__str__()
        m.insertNextCell([0,1,2,5,7,2,-1,1,3])
        self.assertEqual(1,m.getNumberOfCells())
        self.assertTrue(DataArrayInt([8]).isEqual(m.computeNbOfNodesPerCell()))
        self.assertTrue(DataArrayInt([2]).isEqual(m.computeNbOfFacesPerCell()))
        m.__repr__() ; m.__str__()
        m.checkCoherency()
        m.checkCoherency2()
        #
        cm=MEDCouplingCMesh() ; cm.setName("m")
        arr0=DataArrayDouble(6) ; arr0.iota()
        arr1=DataArrayDouble([0,1])
        cm.setCoords(arr0,arr1,arr1) ; um=cm.buildUnstructured() ; um.convertAllToPoly()
        um2=um.deepCpyConnectivityOnly()
        self.assertTrue(um2.isEqual(um,1e-12))
        self.assertEqual(um2.getCoords().getHiddenCppPointer(),um.getCoords().getHiddenCppPointer())
        self.assertTrue(um2.getNodalConnectivity().isEqual(um.getNodalConnectivity()))
        self.assertTrue(um2.getNodalConnectivity().getHiddenCppPointer()!=um.getNodalConnectivity().getHiddenCppPointer())
        self.assertTrue(um2.getNodalConnectivityIndex().isEqual(um.getNodalConnectivityIndex()))
        self.assertTrue(um2.getNodalConnectivityIndex().getHiddenCppPointer()!=um.getNodalConnectivityIndex().getHiddenCppPointer())
        #
        self.assertRaises(InterpKernelException,MEDCoupling1SGTUMesh.New,"m",NORM_POLYHED)
        m=MEDCoupling1DGTUMesh("m",NORM_POLYHED)
        m.allocateCells(5)
        self.assertEqual(15,m.getNodalConnectivity().getNbOfElemAllocated())
        self.assertEqual(6,m.getNodalConnectivityIndex().getNbOfElemAllocated())
        m.setCoords(um.getCoords())
        m.insertNextCell([1,0,6,7,-1,7,6,1])
        self.assertEqual(1,m.getNumberOfCells())
        m.insertNextCell([2,1,7,8,-1,2,1,-1,8,-1,7])
        m.insertNextCell([3,2,8,9])
        m.insertNextCell([4,3,9,10,-1,5,3,9])
        m.insertNextCell([5,4,10,11,-1,11,10,-1,5])
        m.checkCoherency()
        m.checkCoherency2()
        self.assertEqual(5,m.getNumberOfCells())
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,8,19,23,31,40])))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,-1,7,6,1,2,1,7,8,-1,2,1,-1,8,-1,7,3,2,8,9,4,3,9,10,-1,5,3,9,5,4,10,11,-1,11,10,-1,5])))
        #
        m4=m.deepCpy()
        self.assertTrue(m.isEqual(m4,1e-12))
        m4.getNodalConnectivity().setIJ(2,0,5)
        self.assertTrue(not m.isEqual(m4,1e-12))
        m4.getNodalConnectivity().setIJ(2,0,6)
        self.assertTrue(m.isEqual(m4,1e-12))
        m4.getNodalConnectivityIndex().setIJ(2,0,21)
        self.assertTrue(not m.isEqual(m4,1e-12))
        m4.getNodalConnectivityIndex().setIJ(2,0,19)
        self.assertTrue(m.isEqual(m4,1e-12))
        m4.getCoords().setIJ(10,1,1.1)
        self.assertTrue(not m.isEqual(m4,1e-12))
        m4.getCoords().setIJ(10,1,1.)
        self.assertTrue(m.isEqual(m4,1e-12))
        m4.getNodalConnectivity().pushBackSilent(7)
        self.assertTrue(not m.isEqual(m4,1e-12))
        self.assertEqual(7,m4.getNodalConnectivity().popBackSilent())
        self.assertTrue(m.isEqual(m4,1e-12))
        m4.setName("m4")
        self.assertTrue(not m.isEqual(m4,1e-12))
        m4.setName("m")
        self.assertTrue(m.isEqual(m4,1e-12))
        #
        self.assertEqual(6,m.getNodalConnectivityIndex().getNbOfElemAllocated())
        self.assertEqual(60,m.getNodalConnectivity().getNbOfElemAllocated())
        self.assertTrue(m.computeNbOfNodesPerCell().isEqual(DataArrayInt([7,8,4,7,7])))
        self.assertTrue(m.computeNbOfFacesPerCell().isEqual(DataArrayInt([2,4,1,2,3])))
        self.assertEqual(m.getNodeIdsOfCell(1),[2,1,7,8,-1,2,1,-1,8,-1,7])
        f=m.computeIsoBarycenterOfNodesPerCell()
        self.assertTrue(DataArrayDouble([(0.5714285714285714,0.5714285714285714,0),(1.5,0.5,0),(2.5,0.5,0),(3.5714285714285712,0.42857142857142855,0),(4.5714285714285712,0.5714285714285714,0)]).isEqual(f,1e-14))
        mu0=m.buildUnstructured()
        o2n=[1,2,0,4,3]
        m2=m.deepCpy()
        m3=m.deepCpyConnectivityOnly()
        self.assertTrue(m3.isEqual(m,1e-12))
        self.assertEqual(m3.getCoords().getHiddenCppPointer(),m.getCoords().getHiddenCppPointer())
        self.assertTrue(m3.getNodalConnectivity().getHiddenCppPointer()!=m.getNodalConnectivity().getHiddenCppPointer())
        self.assertTrue(m3.getNodalConnectivity().isEqual(m.getNodalConnectivity()))
        self.assertTrue(m3.getNodalConnectivityIndex().getHiddenCppPointer()!=m.getNodalConnectivityIndex().getHiddenCppPointer())
        self.assertTrue(m3.getNodalConnectivityIndex().isEqual(m.getNodalConnectivityIndex()))
        m.renumberCells(o2n)
        mu0.renumberCells(o2n)
        self.assertTrue(mu0.isEqual(m.buildUnstructured(),1e-12))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,12,23,32,40])))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([3,2,8,9,1,0,6,7,-1,7,6,1,2,1,7,8,-1,2,1,-1,8,-1,7,5,4,10,11,-1,11,10,-1,5,4,3,9,10,-1,5,3,9])))
        #
        mcpy0=m.buildUnstructured()
        self.assertTrue(isinstance(mcpy0,MEDCouplingUMesh))
        self.assertTrue(mcpy0.getNodalConnectivity().isEqual(DataArrayInt([31,3,2,8,9,31,1,0,6,7,-1,7,6,1,31,2,1,7,8,-1,2,1,-1,8,-1,7,31,5,4,10,11,-1,11,10,-1,5,31,4,3,9,10,-1,5,3,9])))
        self.assertTrue(mcpy0.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,14,26,36,45])))
        self.assertEqual(mcpy0.getAllGeoTypes(),[NORM_POLYHED])
        mcpy0.checkCoherency()
        mcpy0.checkCoherency2()
        mcpy1=mcpy0.convertIntoSingleGeoTypeMesh()
        self.assertTrue(mcpy1.isEqual(m,1e-12))
        #
        m_mrg=MEDCoupling1DGTUMesh.Merge1DGTUMeshes([m2,m,m2])
        self.assertTrue(m_mrg.getNodalConnectivityIndex().isEqual(DataArrayInt([0,8,19,23,31,40,44,52,63,72,80,88,99,103,111,120])))
        self.assertTrue(m_mrg.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,-1,7,6,1,2,1,7,8,-1,2,1,-1,8,-1,7,3,2,8,9,4,3,9,10,-1,5,3,9,5,4,10,11,-1,11,10,-1,5,27,26,32,33,25,24,30,31,-1,31,30,25,26,25,31,32,-1,26,25,-1,32,-1,31,29,28,34,35,-1,35,34,-1,29,28,27,33,34,-1,29,27,33,49,48,54,55,-1,55,54,49,50,49,55,56,-1,50,49,-1,56,-1,55,51,50,56,57,52,51,57,58,-1,53,51,57,53,52,58,59,-1,59,58,-1,53])))
        m_mrg2=MEDCoupling1DGTUMesh.Merge1DGTUMeshesOnSameCoords([m3,m,m3])
        self.assertTrue(m_mrg2.getNodalConnectivityIndex().isEqual(DataArrayInt([0,8,19,23,31,40,44,52,63,72,80,88,99,103,111,120])))
        self.assertTrue(m_mrg2.getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,-1,7,6,1,2,1,7,8,-1,2,1,-1,8,-1,7,3,2,8,9,4,3,9,10,-1,5,3,9,5,4,10,11,-1,11,10,-1,5,3,2,8,9,1,0,6,7,-1,7,6,1,2,1,7,8,-1,2,1,-1,8,-1,7,5,4,10,11,-1,11,10,-1,5,4,3,9,10,-1,5,3,9,1,0,6,7,-1,7,6,1,2,1,7,8,-1,2,1,-1,8,-1,7,3,2,8,9,4,3,9,10,-1,5,3,9,5,4,10,11,-1,11,10,-1,5])))
        a,b=m_mrg2.getReverseNodalConnectivity()
        self.assertTrue(b.isEqual(DataArrayInt([0,3,15,24,33,39,48,54,66,75,84,93,99,99,99,99,99,99,99,99,99,99,99,99,99])))
        self.assertTrue(a.isEqual(DataArrayInt([0,6,10,0,0,1,1,6,6,7,7,10,10,11,11,1,1,2,5,7,7,11,11,12,2,3,3,5,9,9,12,13,13,3,4,8,9,13,14,3,4,4,8,8,9,13,14,14,0,0,6,6,10,10,0,0,1,1,6,6,7,7,10,10,11,11,1,1,2,5,7,7,11,11,12,2,3,3,5,9,9,12,13,13,3,4,4,8,8,9,13,14,14,4,4,8,8,14,14])))
        self.assertTrue(m_mrg2.fillCellIdsToKeepFromNodeIds([7],False).isEqual(DataArrayInt([0,1,6,7,10,11])))
        self.assertTrue(m_mrg2.fillCellIdsToKeepFromNodeIds([0,1,6,7],True).isEqual(DataArrayInt([0,6,10])))
        #
        self.assertTrue(m_mrg2.isPacked())
        self.assertEqual(120,m_mrg2.getNodalConnectivityIndex().popBackSilent())
        self.assertEqual(m_mrg2.getNumberOfCells(),14)
        m_mrg2.checkCoherency2()
        self.assertTrue(not m_mrg2.isPacked())
        m_mrg4,b=m_mrg2.copyWithNodalConnectivityPacked()
        self.assertTrue(not b)
        m_mrg4.checkCoherency2()
        self.assertEqual(m_mrg4.getNumberOfCells(),14)
        self.assertTrue(m_mrg4.getNodalConnectivityIndex().isEqual(m_mrg2.getNodalConnectivityIndex()))
        self.assertEqual(len(m_mrg4.getNodalConnectivity()),111)
        self.assertEqual(len(m_mrg2.getNodalConnectivity()),120)
        self.assertTrue(m_mrg4.getNodalConnectivity().isEqual(m_mrg2.getNodalConnectivity()[:111]))
        #
        m0=m_mrg2[:5]
        m1=m_mrg2[[5,6,7,8,9]]
        m2=m_mrg2[10:]
        self.assertTrue(m1.isEqualWithoutConsideringStr(m,1e-12))
        a,b=m.checkGeoEquivalWith(m0,12,1e-12)
        self.assertTrue(a.isEqual(DataArrayInt(o2n)))
        self.assertTrue(b is None)
        pass

    def testSwig2DADAreIncludedInMe1(self):
        a=DataArrayDouble(30) ; a.iota() ; a.rearrange(3)
        p=DataArrayInt([5,2,1,9])
        b,c=a.areIncludedInMe(a[p],1e-12)
        self.assertTrue(b)
        self.assertTrue(c.isEqual(p))
        d=a[p]
        d.setIJ(3,1,28.1)
        b,c=a.areIncludedInMe(d,1e-12)
        self.assertTrue(not b)
        self.assertTrue(c.isEqual(DataArrayInt([5,2,1,10])))
        pass

    def testSwig2DADesallocate1(self):
        d=DataArrayDouble([(1,2),(6,7),(6,8)]) ; d.setInfoOnComponents(["aa","bbb"])
        self.assertTrue(d.isAllocated())
        d.checkAllocated()
        self.assertEqual(d.getInfoOnComponents(),["aa","bbb"])
        ref=d.getHeapMemorySize()
        d.desallocate()
        self.assertEqual(ref-d.getHeapMemorySize(),6*8)
        self.assertTrue(not d.isAllocated())
        self.assertEqual(d.getInfoOnComponents(),["aa","bbb"])
        self.assertRaises(InterpKernelException,d.checkAllocated)
        #
        d=DataArrayInt([(1,2),(6,7),(6,8)]) ; d.setInfoOnComponents(["aa","bbb"])
        self.assertTrue(d.isAllocated())
        d.checkAllocated()
        self.assertEqual(d.getInfoOnComponents(),["aa","bbb"])
        ref=d.getHeapMemorySize()
        d.desallocate()
        self.assertEqual(ref-d.getHeapMemorySize(),6*4)
        self.assertTrue(not d.isAllocated())
        self.assertEqual(d.getInfoOnComponents(),["aa","bbb"])
        self.assertRaises(InterpKernelException,d.checkAllocated)
        pass

    def testSwig2IsPartStructured1(self):
        #dim 1
        d10=DataArrayInt([2,3,4,5,6,7,8,9,10,11])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d10,[13])
        self.assertTrue(a) ; self.assertEqual(b,[(2,12)])
        d11=DataArrayInt([2,3,4,5,6,7,8,10,9,11])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d11,[13])
        self.assertTrue(not a)
        self.assertRaises(InterpKernelException,MEDCouplingStructuredMesh.IsPartStructured,d10,[11])
        #dim 2
        st=[10,4]
        d20=DataArrayInt([1,2,3,4,11,12,13,14,21,22,23,24])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d20,st)
        self.assertTrue(a) ; self.assertEqual(b,[(1,5),(0,3)])
        d20=DataArrayInt([1,2,3,4,12,11,13,14,21,22,23,24])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d20,st)
        self.assertTrue(not a)
        d20=DataArrayInt([1,2,3,4,11,12,13,15,21,22,23,24])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d20,st)
        self.assertTrue(not a)
        d21=DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d21,st)
        self.assertTrue(a) ; self.assertEqual(b,[(0,10),(0,4)])
        d22=DataArrayInt([1,2,3,4,11,12,13,14,21,22,23,24,31,32,33,34,41,42,43,44])
        self.assertRaises(InterpKernelException,MEDCouplingStructuredMesh.IsPartStructured,d22,st)
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d22,[10,5])
        self.assertTrue(a) ; self.assertEqual(b,[(1,5),(0,5)])
        #dim 3
        d30=DataArrayInt([11,12,13,14,21,22,23,24,51,52,53,54,61,62,63,64])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d30,[10,4,2])
        self.assertTrue(a) ; self.assertEqual(b,[(1,5),(1,3),(0,2)])
        d31=DataArrayInt([11,12,13,14,21,22,24,23,51,52,53,54,61,62,63,64])
        a,b=MEDCouplingStructuredMesh.IsPartStructured(d31,[10,4,2])
        self.assertTrue(not a)
        self.assertRaises(InterpKernelException,MEDCouplingStructuredMesh.IsPartStructured,d30,[10,4,1])
        pass

    def testSwig2PartStructured1(self):
        c=MEDCouplingCMesh() ; c.setName("toto")
        arr0=DataArrayDouble(10); arr0.iota()
        arr1=DataArrayDouble(4) ; arr1.iota(3)
        c.setCoords(arr0,arr1)
        self.assertEqual(c.getNodeGridStructure(),(10,4))
        self.assertEqual(c.getCellGridStructure(),(9,3))
        d20=DataArrayInt([1,2,3,4,10,11,12,13,19,20,21,22])
        self.assertEqual(27,c.getNumberOfCells())
        self.assertEqual(40,c.getNumberOfNodes())
        self.assertEqual(2,c.getMeshDimension())
        c.checkCoherency()
        #
        arr2=MEDCouplingStructuredMesh.BuildExplicitIdsFrom([9,3],[(1,5),(0,3)])
        self.assertTrue(arr2.isEqual(DataArrayInt([1,2,3,4,10,11,12,13,19,20,21,22])))
        # CMesh
        c2=c.buildStructuredSubPart([(1,5),(0,3)])
        c2.checkCoherency()
        self.assertTrue(isinstance(c2,MEDCouplingCMesh))
        self.assertEqual(12,c2.getNumberOfCells())
        self.assertEqual(20,c2.getNumberOfNodes())
        self.assertEqual(2,c2.getMeshDimension())
        self.assertEqual("toto",c2.getName())
        self.assertTrue(c2.getCoordsAt(0).isEqual(DataArrayDouble([1.,2.,3.,4.,5.]),1e-12))
        self.assertTrue(c2.getCoordsAt(1).isEqual(DataArrayDouble([3.,4.,5.,6.]),1e-12))
        #
        a,b=c.buildPartAndReduceNodes(d20)
        a.checkCoherency()
        exp2=DataArrayInt([-1,0,1,2,3,4,-1,-1,-1,-1,-1,5,6,7,8,9,-1,-1,-1,-1,-1,10,11,12,13,14,-1,-1,-1,-1,-1,15,16,17,18,19,-1,-1,-1,-1])
        self.assertTrue(exp2.isEqual(b))
        self.assertTrue(isinstance(a,MEDCouplingCMesh))
        self.assertTrue(a.buildUnstructured().isEqual(c.buildUnstructured().buildPartAndReduceNodes(d20)[0],1e-12))
        # CurveLinearMesh
        c2=MEDCouplingCurveLinearMesh() ; c2.setName("toto")
        c2.setCoords(c.buildUnstructured().getCoords())
        c2.setNodeGridStructure([10,4])
        c2.checkCoherency()
        a,b=c2.buildPartAndReduceNodes(d20)
        a.checkCoherency()
        self.assertTrue(exp2.isEqual(b))
        self.assertTrue(isinstance(a,MEDCouplingCurveLinearMesh))
        self.assertTrue(a.buildUnstructured().isEqual(c2.buildUnstructured().buildPartAndReduceNodes(d20)[0],1e-12))
        pass

    def testSwig2FindPermutationFromFirstToSecond1(self):
        ids1=DataArrayInt([3,1,103,4,6,10,-7,205])
        ids2=DataArrayInt([-7,1,205,10,6,3,103,4])
        ids3=DataArrayInt.FindPermutationFromFirstToSecond(ids1,ids2)
        self.assertTrue(ids3.isEqual(DataArrayInt([5,1,6,7,4,3,0,2])))
        ids2ToTest=ids1.renumber(ids3)
        self.assertTrue(ids2ToTest.isEqual(ids2))
        self.assertRaises(InterpKernelException,DataArrayInt.FindPermutationFromFirstToSecond,DataArrayInt([3,1,103]),DataArrayInt([1,103]))
        self.assertRaises(InterpKernelException,DataArrayInt.FindPermutationFromFirstToSecond,DataArrayInt([3,1,103]),DataArrayInt([1,103,2]))
        self.assertRaises(InterpKernelException,DataArrayInt.FindPermutationFromFirstToSecond,DataArrayInt([3,1,103]),DataArrayInt([1,103,1]))
        self.assertTrue(DataArrayInt.FindPermutationFromFirstToSecond(DataArrayInt([]),DataArrayInt([])).empty())
        pass

    def testSwig2BugStructuredMeshGetNodeIdsOfCell1(self):
        m=MEDCouplingCMesh("mesh")
        coordsX=DataArrayDouble([0,1.1,2.2,3.3,4.4]) ; coordsX.setInfoOnComponents(["XX [m]"])
        coordsY=DataArrayDouble([0,1.7,3.4]) ; coordsY.setInfoOnComponents(["YYY [km]"])
        m.setCoords(coordsX,coordsY)
        self.assertEqual([2,3,8,7],m.getNodeIdsOfCell(2))
        self.assertEqual([3,4,9,8],m.getNodeIdsOfCell(3))
        self.assertEqual([7,8,13,12],m.getNodeIdsOfCell(6))
        self.assertEqual([8,9,14,13],m.getNodeIdsOfCell(7))
        pass

    def testSwig2ThrowOnDAIInvertN2O2ON2(self):
        p1=DataArrayInt([3,5,8])
        p2=DataArrayInt([0,3,4,5,6,7,8,9,10])
        p1.transformWithIndArr(p2.invertArrayN2O2O2N(11))
        self.assertTrue(p1.isEqual(DataArrayInt([1,3,6])))
        self.assertTrue(p2.invertArrayN2O2O2N(11).isEqual(DataArrayInt([0,-1,-1,1,2,3,4,5,6,7,8])))
        self.assertRaises(InterpKernelException,p2.invertArrayN2O2O2N,10)
        pass

    def testSwig2ComputeEffectiveNbOfNodesPerCell1(self):
        coords=DataArrayDouble([ 0.241310763507 , 0.0504777305619 , 0.0682283524903 , 0.252501053866 , -0.0625176732937 , 0.137272639894 ,
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
                 -0.0372812069535 , -0.286740286332 , 0.00963701291166 ], 69, 3)
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
        meshN=MEDCouplingUMesh.New()
        meshN.setName("ForBary")
        meshN.setMeshDimension(3) ; meshN.setCoords(coords)
        meshN.allocateCells(4)
        meshN.insertNextCell(NORM_POLYHED,113,connN);
        meshN.insertNextCell(NORM_POLYHED,99,connN[113:])
        meshN.insertNextCell(NORM_POLYHED,43,connN[212:])
        meshN.insertNextCell(NORM_POLYHED,92,connN[255:])
        d=meshN.computeEffectiveNbOfNodesPerCell()
        e=meshN.computeNbOfNodesPerCell()
        self.assertTrue(d.isEqual(DataArrayInt([32,28,12,26])))
        self.assertTrue(e.isEqual(DataArrayInt([96,84,36,78])))
        m0=MEDCoupling1DGTUMesh(meshN)
        c=MEDCouplingCMesh()
        arr=DataArrayDouble(3) ; arr.iota(10)
        c.setCoords(arr,arr,arr)
        m10=c.buildUnstructured()
        m11=c.build1SGTUnstructured()
        m12=MEDCoupling1SGTUMesh.New(m10)
        self.assertTrue(m12.isEqual(m11,1e-12))
        m12.setCoords(m0.getCoords()) # m12 is not OK geometrically but the aim of the test is only connectivity values
        m3=MEDCoupling1GTUMesh.AggregateOnSameCoordsToUMesh([m12,m0])
        m3.checkCoherency()
        self.assertEqual(m3.getCoords().getHiddenCppPointer(),m12.getCoords().getHiddenCppPointer())
        self.assertTrue(m3.getNodalConnectivity().isEqual(DataArrayInt([18,1,0,3,4,10,9,12,13,18,2,1,4,5,11,10,13,14,18,4,3,6,7,13,12,15,16,18,5,4,7,8,14,13,16,17,18,10,9,12,13,19,18,21,22,18,11,10,13,14,20,19,22,23,18,13,12,15,16,22,21,24,25,18,14,13,16,17,23,22,25,26,31,0,1,3,4,2,-1,1,5,6,7,0,-1,0,7,8,10,11,9,2,-1,1,5,12,14,15,13,3,-1,16,9,2,4,17,-1,4,3,13,18,17,-1,5,6,19,21,20,12,-1,6,7,8,23,22,19,-1,23,24,10,8,-1,25,11,9,16,-1,24,26,25,11,10,-1,12,14,20,-1,27,28,29,15,13,18,-1,14,15,29,30,21,20,-1,26,27,18,17,16,25,-1,22,19,21,30,31,-1,22,31,28,27,26,24,23,-1,31,30,29,28,31,0,7,8,10,11,9,2,-1,32,0,7,35,34,33,-1,32,0,2,37,36,-1,35,7,8,40,39,38,-1,2,37,41,9,-1,40,8,10,44,43,42,-1,41,9,11,44,43,-1,44,11,10,-1,32,33,45,47,46,36,-1,33,34,48,45,-1,35,34,48,50,49,38,-1,41,43,42,46,36,37,-1,38,39,51,49,-1,39,40,42,46,47,52,51,-1,45,47,52,50,48,-1,52,51,49,50,31,6,7,8,23,22,19,-1,6,35,7,-1,6,35,38,19,-1,35,7,8,40,39,38,-1,53,22,19,38,39,54,-1,23,53,54,40,8,-1,53,22,23,-1,39,54,40,31,35,34,48,50,49,38,-1,6,35,34,56,55,5,-1,6,35,38,19,-1,34,56,57,59,58,48,-1,60,61,21,19,38,49,-1,62,50,48,58,-1,60,63,64,62,50,49,-1,5,6,19,21,20,12,-1,55,5,12,65,-1,66,67,65,55,56,57,-1,63,66,57,59,64,-1,64,62,58,59,-1,60,63,66,67,68,61,-1,61,68,20,21,-1,67,68,20,12,65])))
        self.assertTrue(m3.getNodalConnectivityIndex().isEqual(DataArrayInt([0,9,18,27,36,45,54,63,72,186,286,330,423])))
        pass

    def testSwig2Tetrahedrize1(self):
        d=DataArrayInt([0,3,6,10,14,20])
        d2=d.buildExplicitArrOfSliceOnScaledArr(slice(0,5,2))
        self.assertTrue(d2.isEqual(DataArrayInt([0,0,0, 2,2,2,2, 4,4,4,4,4,4])))
        m=MEDCouplingUMesh("Penta6",3)
        m.setCoords(DataArrayDouble([0,0,0,0,1,0,1,0,0,0,0,2,0,1,2,1,0,2],6,3)) ; m.getCoords().setInfoOnComponents(["X","YY","ZZZ"])
        m.allocateCells()
        m.insertNextCell(NORM_PENTA6,[1,2,0,4,5,3])
        st=m.getCoords().getHiddenCppPointer()
        c,a,b=m.tetrahedrize(PLANAR_FACE_5)
        c.checkCoherency2()
        self.assertTrue(a.isEqual(DataArrayInt([0,0,0])))
        self.assertEqual(0,b)
        self.assertEqual(m.getCoords().getHiddenCppPointer(),c.getCoords().getHiddenCppPointer())
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([1,2,0,4,4,3,5,0,5,0,2,4])))
        del m,c
        #
        m2=MEDCouplingUMesh("octa12",3)
        coords=DataArrayDouble([1.,0.,0.,0.5,0.8660254037844386,0.,-0.5,0.8660254037844387,0.,-1.,1.2246467991473532e-16,0.,-0.5,-0.8660254037844384,0.,0.5,-0.866025403784439,0.,1.,0.,2.,0.5,0.8660254037844386,2.,-0.5,0.8660254037844387,2.,-1.,1.2246467991473532e-16,2.,-0.5,-0.8660254037844384,2.,0.5,-0.866025403784439,2.0],12,3)
        m2.setCoords(coords)
        m2.allocateCells()
        m2.insertNextCell(NORM_HEXGP12,[3,2,1,0,5,4,9,8,7,6,11,10])
        c,a,b=m2.tetrahedrize(PLANAR_FACE_5)
        c.checkCoherency2()
        self.assertTrue(a.isEqual(DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0])))
        self.assertEqual(0,b)
        self.assertEqual(c.getCoords().getHiddenCppPointer(),coords.getHiddenCppPointer())
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,2,4,9,9,10,8,4,8,4,2,9,2,5,4,8,8,10,11,4,11,4,5,8,2,1,5,8,8,11,7,5,7,5,1,8,1,0,5,7,7,11,6,5,6,5,0,7])))
        del m2,coords,c
        #
        coords=DataArrayDouble([0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,2.,1.,0.,2.,1.,1.,2.,0.,1.,2.],8,3) ; coords.setInfoOnComponents(["X","YY","ZZZ"])
        m3=MEDCouplingUMesh("hexa8",3)
        m3.setCoords(coords)
        m3.allocateCells(0)
        m3.insertNextCell(NORM_HEXA8,[3,2,1,0,7,6,5,4])
        st=m3.getCoords().getHiddenCppPointer()
        c,a,b=m3.tetrahedrize(PLANAR_FACE_5)
        c.checkCoherency2()
        a.isEqual(DataArrayInt([0,0,0,0,0]))
        self.assertEqual(0,b)
        self.assertEqual(m3.getCoords().getHiddenCppPointer(),coords.getHiddenCppPointer())
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,6,2,1,3,7,6,4,3,0,4,1,6,4,5,1,3,6,1,4])))
        #
        m4=MEDCouplingUMesh("hexa8",3)
        m4.setCoords(coords)
        m4.allocateCells(0)
        m4.insertNextCell(NORM_HEXA8,[3,2,1,0,7,6,5,4])
        c,a,b=m4.tetrahedrize(PLANAR_FACE_6)
        c.checkCoherency2()
        a.isEqual(DataArrayInt([0,0,0,0,0,0]))
        self.assertEqual(0,b)
        self.assertEqual(c.getCoords().getHiddenCppPointer(),coords.getHiddenCppPointer())
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,6,2,5,3,2,1,5,3,7,6,5,3,4,7,5,3,1,0,5,3,0,4,5])))
        #
        m4=MEDCouplingUMesh("hexa8",3)
        m4.setCoords(coords)
        m4.allocateCells(0)
        m4.insertNextCell(NORM_HEXA8,[3,2,1,0,7,6,5,4])
        st=m4.getCoords().getHiddenCppPointer()
        c,a,b=m4.tetrahedrize(GENERAL_24)
        c.checkCoherency2()
        a.isEqual(DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
        self.assertEqual(7,b)
        self.assertTrue(c.getCoords().getHiddenCppPointer()!=coords.getHiddenCppPointer())
        self.assertTrue(c.getCoords()[:8].isEqual(coords,0))
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,7,8,14,7,6,8,14,6,2,8,14,2,3,8,14,3,2,9,14,2,1,9,14,1,0,9,14,0,3,9,14,3,0,10,14,0,4,10,14,4,7,10,14,7,3,10,14,2,6,11,14,6,5,11,14,5,1,11,14,1,2,11,14,7,4,12,14,4,5,12,14,5,6,12,14,6,7,12,14,1,5,13,14,5,4,13,14,4,0,13,14,0,1,13,14])))
        m4CoordsExp=DataArrayDouble([0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,2.,1.,0.,2.,1.,1.,2.,0.,1.,2.,0.5,1.,1.,0.5,0.5,0.,0.,0.5,1.,1.,0.5,1.,0.5,0.5,2.,0.5,0.,1.,0.5,0.5,1.],15,3)
        m4CoordsExp.setInfoOnComponents(["X","YY","ZZZ"])
        self.assertTrue(c.getCoords().isEqual(m4CoordsExp,1e-12))
        self.assertAlmostEqual(2.,c.getMeasureField(False).accumulate()[0],12)
        #
        m6=MEDCouplingUMesh("hexa8",3)
        m6.setCoords(coords)
        m6.allocateCells(0)
        m6.insertNextCell(NORM_HEXA8,[3,2,1,0,7,6,5,4])
        st=m6.getCoords().getHiddenCppPointer()
        c,a,b=m6.tetrahedrize(GENERAL_48)
        c.checkCoherency2()
        a.isEqual(DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
        self.assertEqual(19,b)
        self.assertTrue(c.getCoords().getHiddenCppPointer()!=coords.getHiddenCppPointer())
        self.assertTrue(c.getCoords()[:8].isEqual(coords,0))
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,20,8,26,3,8,21,26,3,9,20,26,3,22,9,26,3,21,12,26,3,12,22,26,8,10,2,23,8,2,13,23,8,20,10,23,8,26,20,23,8,13,21,23,8,21,26,23,12,26,21,25,12,21,16,25,12,22,26,25,12,17,22,25,12,16,0,25,12,0,17,25,21,23,13,18,21,13,1,18,21,26,23,18,21,25,26,18,21,1,16,18,21,16,25,18,9,11,20,24,9,20,26,24,9,7,11,24,9,14,7,24,9,26,22,24,9,22,14,24,20,6,10,15,20,10,23,15,20,11,6,15,20,24,11,15,20,23,26,15,20,26,24,15,22,24,26,19,22,26,25,19,22,14,24,19,22,4,14,19,22,25,17,19,22,17,4,19,26,15,23,5,26,23,18,5,26,24,15,5,26,19,24,5,26,18,25,5,26,25,19,5])))
        m6CoordsExp=DataArrayDouble([0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,2.,1.,0.,2.,1.,1.,2.,0.,1.,2.,0.5,1.,0.,0.,1.,1.,1.,1.,1.,0.5,1.,2.,0.,0.5,0.,1.,0.5,0.,0.,0.5,2.,1.,0.5,2.,0.5,0.,0.,0.,0.,1.,1.,0.,1.,0.5,0.,2.,0.5,1.,1.,0.5,0.5,0.,0.,0.5,1.,1.,0.5,1.,0.5,0.5,2.,0.5,0.,1.,0.5,0.5,1.],27,3)
        m6CoordsExp.setInfoOnComponents(["X","YY","ZZZ"])
        self.assertTrue(c.getCoords().isEqual(m6CoordsExp,1e-12))
        self.assertAlmostEqual(2.,c.getMeasureField(False).accumulate()[0],12)
        #
        m7=MEDCouplingUMesh("polyhed",3)
        coords=DataArrayDouble([1.,0.,0.,0.5,0.8660254037844386,0.,-0.5,0.8660254037844387,0.,-1.,0.,0.,-0.5,-0.8660254037844384,0.,0.5,-0.866025403784439,0.,1.,0.,2.,0.5,0.8660254037844386,2.,-0.5,0.8660254037844387,2.,-1.,0.,2.,-0.5,-0.8660254037844384,2.,0.5,-0.866025403784439,2.0],12,3) ; coords.setInfoOnComponents(["X","YY","ZZZ"])
        m7.setCoords(coords)
        m7.allocateCells()
        m7.insertNextCell(NORM_POLYHED,[3,2,1,0,5,4,-1,9,10,11,6,7,8,-1,3,9,8,2,-1,2,8,7,1,-1,1,7,6,0,-1,0,6,11,5,-1,5,11,10,4,-1,4,10,9,3])
        c,a,b=m7.tetrahedrize(PLANAR_FACE_5)
        c.checkCoherency2()
        self.assertTrue(a.isEqual(DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])))
        self.assertEqual(9,b)
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,2,12,20,2,1,12,20,1,0,12,20,0,5,12,20,5,4,12,20,4,3,12,20,9,10,13,20,10,11,13,20,11,6,13,20,6,7,13,20,7,8,13,20,8,9,13,20,3,9,14,20,9,8,14,20,8,2,14,20,2,3,14,20,2,8,15,20,8,7,15,20,7,1,15,20,1,2,15,20,1,7,16,20,7,6,16,20,6,0,16,20,0,1,16,20,0,6,17,20,6,11,17,20,11,5,17,20,5,0,17,20,5,11,18,20,11,10,18,20,10,4,18,20,4,5,18,20,4,10,19,20,10,9,19,20,9,3,19,20,3,4,19,20])))
        self.assertAlmostEqual(5.196152422706635,c.getMeasureField(False).accumulate()[0],12)
        m7CoordsExp=DataArrayDouble([1.0,0.0,0.0,0.5,0.8660254037844386,0.0,-0.5,0.8660254037844387,0.0,-1.0,0.,0.0,-0.5,-0.8660254037844384,0.0,0.5,-0.866025403784439,0.0,1.0,0.0,2.0,0.5,0.8660254037844386,2.0,-0.5,0.8660254037844387,2.0,-1.0,0.,2.0,-0.5,-0.8660254037844384,2.0,0.5,-0.866025403784439,2.0,0.0,0.0,0.0,0.0,0.,2.0,-0.75,0.4330127018922194,1.0,0.0,0.8660254037844386,1.0,0.75,0.4330127018922193,1.0,0.75,-0.4330127018922195,1.0,0.0,-0.8660254037844387,1.0,-0.75,-0.4330127018922191,1.0,0.0,0.,1.0],21,3)
        m7CoordsExp.setInfoOnComponents(["X","YY","ZZZ"])
        self.assertTrue(c.getCoords().isEqual(m7CoordsExp,1e-12))
        del m7,coords,c
        #
        coords=DataArrayDouble([0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,2.,1.,0.,2.,1.,1.,2.,0.,1.,2.],8,3) ; coords.setInfoOnComponents(["X","YY","ZZZ"])
        m8=MEDCouplingUMesh("pyra5",3)
        m8.setCoords(coords)
        m8.allocateCells(0)
        m8.insertNextCell(NORM_PYRA5,[3,2,1,0,7])
        st=m8.getCoords().getHiddenCppPointer()
        c,a,b=m8.tetrahedrize(PLANAR_FACE_5)
        self.assertEqual(m8.getCoords().getHiddenCppPointer(),coords.getHiddenCppPointer())
        c.checkCoherency2()
        self.assertTrue(a.isEqual(DataArrayInt([0,0])))
        self.assertEqual(0,b)
        self.assertTrue(c.getNodalConnectivity().isEqual(DataArrayInt([3,2,1,7,3,1,0,7])))
        self.assertAlmostEqual(0.6666666666666667,c.getMeasureField(False).accumulate()[0],12)
        pass

    def testDualMesh3D1(self):
        arr=DataArrayDouble(2) ; arr.iota()
        c=MEDCouplingCMesh() ; c.setCoords(arr,arr,arr)
        m=c.buildUnstructured()
        t=m.tetrahedrize(PLANAR_FACE_5)[0]
        d=t.computeDualMesh()
        self.assertTrue(d.getNodalConnectivityIndex().isEqual(DataArrayInt([0,29,118,207,236,325,354,383,472])))
        self.assertTrue(d.getNodalConnectivity().isEqual(DataArrayInt([26,11,42,8,-1,25,8,42,10,-1,29,10,42,11,-1,0,26,8,25,-1,0,25,10,29,-1,0,29,11,26,24,9,42,8,-1,26,8,42,11,-1,27,11,42,9,-1,1,24,8,26,-1,1,26,11,27,-1,30,13,43,12,-1,24,12,43,15,-1,32,15,43,13,-1,1,30,12,24,-1,1,32,13,30,-1,35,17,44,16,-1,32,16,44,19,-1,27,19,44,17,-1,1,35,16,32,-1,1,27,17,35,-1,24,15,46,9,-1,27,9,46,19,-1,32,19,46,15,27,9,42,11,-1,29,11,42,10,-1,28,10,42,9,-1,2,29,10,28,-1,2,27,11,29,-1,27,17,44,19,-1,38,19,44,18,-1,37,18,44,17,-1,2,37,17,27,-1,2,38,18,37,-1,28,21,45,23,-1,41,23,45,22,-1,38,22,45,21,-1,2,41,22,38,-1,2,28,23,41,-1,27,19,46,9,-1,28,9,46,21,-1,38,21,46,19,35,16,44,17,-1,36,18,44,16,-1,37,17,44,18,-1,3,36,16,35,-1,3,35,17,37,-1,3,37,18,36,24,8,42,9,-1,25,10,42,8,-1,28,9,42,10,-1,4,25,8,24,-1,4,28,10,25,-1,24,15,43,12,-1,31,12,43,14,-1,34,14,43,15,-1,4,24,12,31,-1,4,31,14,34,-1,34,21,45,20,-1,40,20,45,23,-1,28,23,45,21,-1,4,34,20,40,-1,4,40,23,28,-1,24,9,46,15,-1,28,21,46,9,-1,34,15,46,21,30,12,43,13,-1,31,14,43,12,-1,33,13,43,14,-1,5,31,12,30,-1,5,30,13,33,-1,5,33,14,31,40,23,45,20,-1,39,20,45,22,-1,41,22,45,23,-1,6,40,20,39,-1,6,39,22,41,-1,6,41,23,40,32,13,43,15,-1,34,15,43,14,-1,33,14,43,13,-1,7,33,13,32,-1,7,34,14,33,-1,32,19,44,16,-1,36,16,44,18,-1,38,18,44,19,-1,7,32,16,36,-1,7,36,18,38,-1,34,20,45,21,-1,39,22,45,20,-1,38,21,45,22,-1,7,39,20,34,-1,7,38,22,39,-1,32,15,46,19,-1,38,19,46,21,-1,34,21,46,15])))
        self.assertTrue(d.getCoords().isEqual(DataArrayDouble([0.,0.,0.,1.,0.,0.,0.,1.,0.,1.,1.,0.,0.,0.,1.,1.,0.,1.,0.,1.,1.,1.,1.,1.,0.3333333333333333,0.,0.3333333333333333,0.3333333333333333,0.3333333333333333,0.3333333333333333,0.,0.3333333333333333,0.3333333333333333,0.3333333333333333,0.3333333333333333,0.,0.6666666666666666,0.,0.6666666666666666,1.,0.3333333333333333,0.6666666666666666,0.6666666666666666,0.3333333333333333,1.,0.6666666666666666,0.3333333333333333,0.6666666666666666,1.,0.6666666666666666,0.3333333333333333,0.6666666666666666,0.6666666666666666,0.,0.6666666666666666,1.,0.3333333333333333,0.6666666666666666,0.6666666666666666,0.3333333333333333,0.3333333333333333,0.6666666666666666,1.,0.3333333333333333,0.6666666666666666,0.6666666666666666,0.3333333333333333,1.,0.6666666666666666,0.,0.6666666666666666,0.6666666666666666,0.5,0.,0.5,0.,0.,0.5,0.5,0.,0.,0.5,0.5,0.,0.,0.5,0.5,0.,0.5,0.,1.,0.,0.5,0.5,0.,1.,1.,0.5,0.5,1.,0.5,1.,0.5,0.5,1.,1.,0.5,0.,1.,1.,0.5,0.5,1.,0.,0.5,1.,0.5,0.5,1.,1.,0.,0.5,1.,0.,1.,0.5,0.25,0.25,0.25,0.75,0.25,0.75,0.75,0.75,0.25,0.25,0.75,0.75,0.5,0.5,0.5],47,3),1e-12))
        self.assertAlmostEqual(1.,d.getMeasureField(False).accumulate()[0],1e-13)
        pass

    def testDualMesh2D1(self):
        arr=DataArrayDouble(5) ; arr.iota()
        c=MEDCouplingCMesh() ; c.setCoords(arr,arr)
        m=c.buildUnstructured()
        m.simplexize(0)
        t=MEDCoupling1SGTUMesh(m)
        d=t.computeDualMesh()
        self.assertTrue(d.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,12,20,28,34,42,54,66,78,86,94,106,118,130,138,146,158,170,182,190,196,204,212,220,224])))
        self.assertTrue(d.getNodalConnectivity().isEqual(DataArrayInt([26,81,25,0,25,81,27,82,29,83,30,1,30,83,31,84,33,85,34,2,34,85,35,86,37,87,38,3,38,87,39,88,41,4,27,81,26,5,42,89,28,82,29,82,28,89,43,90,45,91,32,84,31,83,33,84,32,91,46,92,48,93,36,86,35,85,37,86,36,93,49,94,51,95,40,88,39,87,41,88,40,95,52,96,54,9,43,89,42,10,55,97,44,90,45,90,44,97,56,98,58,99,47,92,46,91,48,92,47,99,59,100,61,101,50,94,49,93,51,94,50,101,62,102,64,103,53,96,52,95,54,96,53,103,65,104,67,14,56,97,55,15,68,105,57,98,58,98,57,105,69,106,71,107,60,100,59,99,61,100,60,107,72,108,74,109,63,102,62,101,64,102,63,109,75,110,77,111,66,104,65,103,67,104,66,111,78,112,80,19,69,105,68,20,70,106,71,106,70,21,73,108,72,107,74,108,73,22,76,110,75,109,77,110,76,23,79,112,78,111,80,112,79,24])))
        self.assertTrue(d.getCoords().isEqual(DataArrayDouble([0.,0.,1.,0.,2.,0.,3.,0.,4.,0.,0.,1.,1.,1.,2.,1.,3.,1.,4.,1.,0.,2.,1.,2.,2.,2.,3.,2.,4.,2.,0.,3.,1.,3.,2.,3.,3.,3.,4.,3.,0.,4.,1.,4.,2.,4.,3.,4.,4.,4.,0.5,0.,0.,0.5,0.5,0.5,0.5,1.,1.,0.5,1.5,0.,1.5,0.5,1.5,1.,2.,0.5,2.5,0.,2.5,0.5,2.5,1.,3.,0.5,3.5,0.,3.5,0.5,3.5,1.,4.,0.5,0.,1.5,0.5,1.5,0.5,2.,1.,1.5,1.5,1.5,1.5,2.,2.,1.5,2.5,1.5,2.5,2.,3.,1.5,3.5,1.5,3.5,2.,4.,1.5,0.,2.5,0.5,2.5,0.5,3.,1.,2.5,1.5,2.5,1.5,3.,2.,2.5,2.5,2.5,2.5,3.,3.,2.5,3.5,2.5,3.5,3.,4.,2.5,0.,3.5,0.5,3.5,0.5,4.,1.,3.5,1.5,3.5,1.5,4.,2.,3.5,2.5,3.5,2.5,4.,3.,3.5,3.5,3.5,3.5,4.,4.,3.5,0.3333333333333333,0.3333333333333333,0.6666666666666666,0.6666666666666666,1.3333333333333333,0.3333333333333333,1.6666666666666665,0.6666666666666666,2.333333333333333,0.3333333333333333,2.6666666666666665,0.6666666666666666,3.333333333333333,0.3333333333333333,3.6666666666666665,0.6666666666666666,0.3333333333333333,1.3333333333333333,0.6666666666666666,1.6666666666666665,1.3333333333333333,1.3333333333333333,1.6666666666666665,1.6666666666666665,2.333333333333333,1.3333333333333333,2.6666666666666665,1.6666666666666665,3.333333333333333,1.3333333333333333,3.6666666666666665,1.6666666666666665,0.3333333333333333,2.333333333333333,0.6666666666666666,2.6666666666666665,1.3333333333333333,2.333333333333333,1.6666666666666665,2.6666666666666665,2.333333333333333,2.333333333333333,2.6666666666666665,2.6666666666666665,3.333333333333333,2.333333333333333,3.6666666666666665,2.6666666666666665,0.3333333333333333,3.333333333333333,0.6666666666666666,3.6666666666666665,1.3333333333333333,3.333333333333333,1.6666666666666665,3.6666666666666665,2.333333333333333,3.333333333333333,2.6666666666666665,3.6666666666666665,3.333333333333333,3.333333333333333,3.6666666666666665,3.6666666666666665],113,2),1e-12))
        self.assertAlmostEqual(16.,d.getMeasureField(False).accumulate()[0],1e-13)
        pass

    def testSwig2LoadBalanceBBox1(self):
        arr=DataArrayDouble(5) ; arr.iota()
        t=MEDCouplingCMesh() ; t.setCoords(arr,arr)
        arr=DataArrayDouble(16) ; arr.iota() ; arr*=2./15
        s=MEDCouplingCMesh() ; s.setCoords(arr,arr[:]) ; s.translate([2.,1.])
        #
        s1=s.build1SGTUnstructured()
        t1=t.build1SGTUnstructured()
        w=MEDCouplingPointSet.ComputeNbOfInteractionsWithSrcCells(s1,t1,1e-12)
        wExp=DataArrayInt([0,0,0,0,0,0,64,64,0,0,64,64,0,0,0,0])
        self.assertTrue(w.isEqual(wExp))
        slcs=w.splitInBalancedSlices(4)
        self.assertEqual(len(slcs),4)
        self.assertEqual(slcs,[slice(0,7,1),slice(7,8,1),slice(8,11,1),slice(11,16,1)])
        bbs=s1.getBoundingBoxForBBTree()
        bbt=t1.getBoundingBoxForBBTree()
        self.assertTrue(bbt.computeNbOfInteractionsWith(bbs,1e-12).isEqual(wExp))
        pass

    def testKrSpatialDiscretization2(self):
        srcPointCoordsXY=DataArrayDouble([0.8401877171547095,0.39438292681909304,0.7830992237586059,0.7984400334760733,0.9116473579367843,0.19755136929338396,0.335222755714889,0.768229594811904,0.2777747108031878,0.5539699557954305,0.47739705186216025,0.6288709247619244,0.36478447279184334,0.5134009101956155,0.9522297251747128,0.9161950680037007,0.6357117279599009,0.7172969294326831,0.14160255535580338,0.6069688762570586,0.01630057162432958,0.24288677062973696,0.13723157678601872,0.8041767542269904,0.15667908925408455,0.4009443942461835,0.12979044678145574,0.10880880202576929,0.998924518003559,0.21825690531090688,0.5129323944043984,0.8391122346926072,0.6126398325956612,0.29603161769734304,0.6375522677030192,0.5242871900667843,0.493582986990727,0.9727750238835695,0.29251678441302703,0.7713576977939148,0.5267449792133388,0.7699138362751873,0.4002286220901779,0.8915294520051822,0.2833147460051415,0.3524583472648907,0.8077245200088827,0.9190264739650424,0.06975527623191256,0.9493270753646861,0.5259953502221011,0.08605584785624214,0.19221384599442307,0.6632269270081198,0.8902326025488938,0.3488929352485076,0.06417132078864207,0.02002304886468828,0.4577017372742769,0.06309583832653977,0.23827995417559517,0.9706341316786754,0.9022080734848082,0.8509197867712563,0.2666657493760184,0.5397603407221662,0.3752069763723793,0.7602487363667454,0.5125353641400744,0.6677237607854063,0.5316064341606602,0.039280343353413204,0.4376375965949323,0.9318350562508382,0.9308097953585953,0.7209523430657351,0.28429340305006756,0.7385343149018168,0.6399788165651163,0.3540486797476414,0.687861390266503,0.16597416632155615,0.4401045276038835,0.880075236260926,0.829201093329676,0.3303371296871161,0.22896817104377232,0.8933724145839793,0.35036017855180435,0.6866699083180492,0.9564682529105192,0.5886401331930609,0.6573040395310633,0.8586763259296661,0.4395599194986559,0.9239697889070817,0.39843666665183225,0.8147668963366965,0.6842185252738271,0.9109720307919067,0.4824906566564416,0.21582495896882609,0.9502523741453198,0.9201282537170352,0.14766001475400292,0.8810621695039152,0.641080596317109,0.43195341826973177,0.6195964839400707,0.281059412416564,0.7860020980173732,0.3074578737409124,0.44703357920378145,0.22610662515559543,0.18753310953617705,0.27623467206779617,0.5564437553083728,0.4165012805799494,0.16960708618611428,0.9068039338601771,0.10317118843233734,0.1260753390966334,0.49544406658757667,0.7604752284290619,0.9847516650262995,0.9350039865518939,0.6844450168704823,0.3831883312124705,0.7497708824229291,0.36866354167864823,0.2941603620043771,0.2322615386137094,0.5844885006474743,0.24441273568403568,0.15238979186508328,0.7321485158671385,0.12547490472228962,0.7934703881821923,0.164101933671209,0.7450713891280216,0.07452980059875632,0.9501040316885822,0.05252926240327268,0.5215633798025378,0.1762106563785163,0.24006237240511102,0.797798051870334,0.732654411686889,0.6565636529850605,0.9674051385221095,0.6394583455470663,0.7597348418830591,0.09348047715308166,0.13490241166898162,0.5202100698464597,0.07823214171371988,0.06990639775521419,0.2046550862512808,0.4614204733918516,0.8196772801781433,0.5733186283955903,0.7555808353962288,0.05193881879185271,0.1578071285774033,0.9999935710802644,0.204328610656936,0.8899556444445419,0.12546847580255405,0.9977989993047895,0.054057577650089554,0.8705398649305757,0.07232879943788462,0.004161608873010431,0.9230691273338484,0.5938921792404224,0.180372265717188,0.16313149927329806,0.3916902306450951,0.9130266774040771,0.8196951527240198,0.35909536870154335,0.552485022485482,0.5794299941414176,0.452575845854625,0.687387434620125,0.09964006352221597,0.5308079880340062,0.7572938323753392,0.30429514977349675,0.9922284614258579,0.5769711125534824,0.877613778169087,0.7478092963564253,0.6289099313453351,0.03542090674649035,0.7478028669710285,0.8332385420022712,0.9253765511910322,0.8732713427735824,0.8310375408413995],100,2)
        srcFieldValsOnPoints=DataArrayDouble([0.7643742528498438,-0.023507696856211995,1.1082895131907775,0.6299357452572031,0.8892623544912389,0.72212114810697,0.9196401044320336,-0.759961711221917,0.40801932617748826,0.8441134300809151,0.982483804252809,0.6752368914020778,0.9924403977479798,1.1063334970204484,0.9403055261137516,0.3624481886322733,1.1344772505996308,0.7522965618948239,0.17077741651388564,0.6504551671311436,0.45843479588425423,0.41098905950326753,1.0681420394050904,-0.3483587903820091,0.5620151050607809,1.384969776596035,0.7948875141132845,0.7931192000237167,1.062498042490183,1.3709072529577366,0.44929346605311893,-0.4469683401788374,0.9035857424514101,0.6137249300593463,0.6355610879026966,1.4318174829507697,0.3097567072129551,-0.20515052260807165,0.6922559820922779,1.0341638749443423,1.3072652153341024,0.38511367353000436,0.9160514929274943,0.54513408530581,0.722252267913328,0.06684522818576251,0.10571899758067793,0.3193844999960903,0.5213532270828706,-0.04834998649603944,1.2408805068350615,-0.7632951295676795,0.5980054665011202,0.9064738717547436,1.1541070755096696,1.008234260272265,1.2225806960553827,1.0788560195121106,0.9818990282104452,0.5621951325841853,1.0796757508374188,0.5082872315589883,-0.9153702001062469,0.9560418838920791,0.9251098559152824,1.1603063610984021,1.2122303611181837,0.7379539363312343,0.6877611899207183,0.723966552446608,0.5596025827162566,0.8849725005989729,1.0908363665075547,0.08956512916455672,-0.10247645571248344,0.3236718069555875,1.069478546398975,1.3900071080692746,1.0322398863403262,0.45315515354558034,0.4249870238786733,1.030226761858634,0.974024629584669,1.2838885424020365,1.3451943506525155,1.4029933267831995,0.6025539675442462,1.2947650597767038,1.0006061239483002,-0.4017336259949164,0.8771165113201297,0.9158909024218246,1.403798605551443,0.4742904006425974,0.3671787905896653,0.20646491720419674,0.40739337434288925,0.7341932402033597,-0.4295893651836911,-0.3187777570661546],100,1)
        targetPointCoordsXY=DataArrayDouble([-0.5,-0.5,-0.5,-0.35,-0.5,-0.2,-0.5,-0.05,-0.5,0.1,-0.5,0.25,-0.5,0.4,-0.5,0.55,-0.5,0.7,-0.5,0.85,-0.5,1.0,-0.5,1.15,-0.5,1.3,-0.5,1.45,-0.35,-0.5,-0.35,-0.35,-0.35,-0.2,-0.35,-0.05,-0.35,0.1,-0.35,0.25,-0.35,0.4,-0.35,0.55,-0.35,0.7,-0.35,0.85,-0.35,1.0,-0.35,1.15,-0.35,1.3,-0.35,1.45,-0.2,-0.5,-0.2,-0.35,-0.2,-0.2,-0.2,-0.05,-0.2,0.1,-0.2,0.25,-0.2,0.4,-0.2,0.55,-0.2,0.7,-0.2,0.85,-0.2,1.0,-0.2,1.15,-0.2,1.3,-0.2,1.45,-0.05,-0.5,-0.05,-0.35,-0.05,-0.2,-0.05,-0.05,-0.05,0.1,-0.05,0.25,-0.05,0.4,-0.05,0.55,-0.05,0.7,-0.05,0.85,-0.05,1.0,-0.05,1.15,-0.05,1.3,-0.05,1.45,0.1,-0.5,0.1,-0.35,0.1,-0.2,0.1,-0.05,0.1,0.1,0.1,0.25,0.1,0.4,0.1,0.55,0.1,0.7,0.1,0.85,0.1,1.0,0.1,1.15,0.1,1.3,0.1,1.45,0.25,-0.5,0.25,-0.35,0.25,-0.2,0.25,-0.05,0.25,0.1,0.25,0.25,0.25,0.4,0.25,0.55,0.25,0.7,0.25,0.85,0.25,1.0,0.25,1.15,0.25,1.3,0.25,1.45,0.4,-0.5,0.4,-0.35,0.4,-0.2,0.4,-0.05,0.4,0.1,0.4,0.25,0.4,0.4,0.4,0.55,0.4,0.7,0.4,0.85,0.4,1.0,0.4,1.15,0.4,1.3,0.4,1.45,0.55,-0.5,0.55,-0.35,0.55,-0.2,0.55,-0.05,0.55,0.1,0.55,0.25,0.55,0.4,0.55,0.55,0.55,0.7,0.55,0.85,0.55,1.0,0.55,1.15,0.55,1.3,0.55,1.45,0.7,-0.5,0.7,-0.35,0.7,-0.2,0.7,-0.05,0.7,0.1,0.7,0.25,0.7,0.4,0.7,0.55,0.7,0.7,0.7,0.85,0.7,1.0,0.7,1.15,0.7,1.3,0.7,1.45,0.85,-0.5,0.85,-0.35,0.85,-0.2,0.85,-0.05,0.85,0.1,0.85,0.25,0.85,0.4,0.85,0.55,0.85,0.7,0.85,0.85,0.85,1.0,0.85,1.15,0.85,1.3,0.85,1.45,1.0,-0.5,1.0,-0.35,1.0,-0.2,1.0,-0.05,1.0,0.1,1.0,0.25,1.0,0.4,1.0,0.55,1.0,0.7,1.0,0.85,1.0,1.0,1.0,1.15,1.0,1.3,1.0,1.45,1.15,-0.5,1.15,-0.35,1.15,-0.2,1.15,-0.05,1.15,0.1,1.15,0.25,1.15,0.4,1.15,0.55,1.15,0.7,1.15,0.85,1.15,1.0,1.15,1.15,1.15,1.3,1.15,1.45,1.3,-0.5,1.3,-0.35,1.3,-0.2,1.3,-0.05,1.3,0.1,1.3,0.25,1.3,0.4,1.3,0.55,1.3,0.7,1.3,0.85,1.3,1.0,1.3,1.15,1.3,1.3,1.3,1.45,1.45,-0.5,1.45,-0.35,1.45,-0.2,1.45,-0.05,1.45,0.1,1.45,0.25,1.45,0.4,1.45,0.55,1.45,0.7,1.45,0.85,1.45,1.0,1.45,1.15,1.45,1.3,1.45,1.45],196,2)
        targetFieldValsExpected=DataArrayDouble([1.645976003316459, 1.454458180060204, 1.286087532859835, 1.147305389930914, 1.040143042030752, 0.9592075185603157, 0.8932542207607532, 0.8296417057622609, 0.7572539678257579, 0.6669048311361028, 0.551329882743212, 0.4064445075734602, 0.2323703965460786, 0.03253142054561309, 1.615321686989539, 1.414941300553572, 1.238383118538708, 1.096701655702075, 0.9955792747382535, 0.9271194507282707, 0.8741000712825546, 0.8201879508155141, 0.7537335933761495, 0.6656210809234322, 0.5470285414729397, 0.3927301586610237, 0.2044036897887453, -0.01181672742825013, 1.609602552867195, 1.400625195269133, 1.213287847440801, 1.065318574929208, 0.9717609562002842, 0.9182626517777217, 0.8760698972315855, 0.8258196104516153, 0.7586487405165288, 0.6686168424854784, 0.5434121624038266, 0.3741815029337978, 0.1661376046619205, -0.0704038088420833, 1.635421686625182, 1.422642113482769, 1.225977424080963, 1.066864693789366, 0.9864801043792362, 0.9486639217909161, 0.9075176697327381, 0.8471248730261529, 0.7660983406349626, 0.6675300501188994, 0.5320013361909732, 0.3404583135353376, 0.1074346390951333, -0.1520751802856468, 1.695346918429566, 1.489526279573347, 1.297678617961701, 1.139921240332637, 1.080508463804929, 1.036847769764088, 0.9687840669352359, 0.8790397822170175, 0.76938768351059, 0.6441978169925557, 0.4915328571013788, 0.2742929463574293, 0.0148214290833748, -0.2671755287427691, 1.782761788232491, 1.59423004798623, 1.422317125787222, 1.286999529473285, 1.20500638941831, 1.127058114031519, 1.022332539190471, 0.8945753999401338, 0.7469190939381181, 0.582396906110898, 0.4015920181411496, 0.1584700483835366, -0.1251860255418387, -0.4254052799545267, 1.881794862747652, 1.712890309994015, 1.557517508390291, 1.422727414977963, 1.308048056353061, 1.187569766723152, 1.03942150436647, 0.8677583087532357, 0.6766652050643343, 0.4703897480238999, 0.2497994532908829, -0.02005989176786582, -0.3224387891441491, -0.6331519303649853, 1.973114284621266, 1.820187301531605, 1.673403730111759, 1.528504440482262, 1.379693463484634, 1.207642134784147, 1.008217764780293, 0.7863328498822348, 0.5465383049529959, 0.2944879513187435, 0.03250657765404452, -0.2670900851421072, -0.5806516907976924, -0.8911331026431459, 2.038729888975378, 1.895652364645637, 1.751759791756183, 1.594035761810714, 1.403016809171641, 1.171403152610878, 0.913267035125007, 0.6343281031932027, 0.3434843176189371, 0.04195410032095204, -0.2645533663891493, -0.58577400250975, -0.8958218846257981, -1.192230697656513, 2.064018033720731, 1.922048791644444, 1.773847180028208, 1.600340336378483, 1.361620036333164, 1.060873411411508, 0.7373484802125152, 0.3868966266761109, 0.04316272760227413, -0.3009370030949727, -0.6505233805563486, -0.9669887470696283, -1.250005719852354, -1.519122595631787, 2.039938287785342, 1.887400820799651, 1.722008733683987, 1.523879290022419, 1.23834392230135, 0.8606985727866472, 0.4844892131548788, 0.08077959236877175, -0.3195742594962179, -0.726291368696764, -1.094357645641832, -1.359078900303776, -1.604725656501341, -1.845297168323687, 1.965762248218393, 1.791665198563286, 1.595056719739704, 1.353692777435502, 1.033006623003495, 0.6416349531117889, 0.2290046916364761, -0.1993180965088852, -0.6311618804827295, -1.051489875129883, -1.409404344854132, -1.681249363331096, -1.917859637689007, -2.145034400762945, 1.849053542205925, 1.648479366622312, 1.418493963148431, 1.141939527533839, 0.8042385795619003, 0.4127534639189761, -0.008572116677791453, -0.4428317297963555, -0.8745477268718713, -1.281769237471681, -1.635421857742795, -1.926210204560556, -2.175577364628722, -2.405762639746138, 1.701519686999922, 1.475879908746998, 1.219065416294153, 0.9203732349759972, 0.5740137315474942, 0.1856460506119944, -0.2298288912529738, -0.6558565521653752, -1.075391078040103, -1.469402631469075, -1.820558929095151, -2.123592211415966, -2.388177455227765, -2.628832075944413])
        coeffsExpected=DataArrayDouble([0.3953237723894342,-0.17220705170185724,0.620727139132215,-0.01938292763088709,-0.007524685306185282,0.0016277944443884584,-0.0005209587893117361,-1.8992696595839718,-0.13154330748345855,0.11248800965389728,-0.47310750305033406,0.03685741122098605,0.21362468750754374,0.8082608687799991,-0.6775548200221704,-0.027683208482275873,-0.007806877014495724,-0.013539239795959668,0.3478535665778018,0.005145793726360813,0.03708618549628136,-0.18235332489209385,-0.04517273339177797,-0.081755114492025,0.12791746560435255,0.09659355695676189,-0.024809653129318366,0.08327587452569823,-1.790380673650165,-0.10622983512164165,0.14989029282340274,0.05949513762355707,0.004548072841131278,0.011252095917834793,-0.004848057194721367,-0.2658537133108412,0.016651579133606154,-0.021640915366981317,0.008975511042160175,-0.021052213988815974,-0.09347841701844657,0.03533229488135717,-0.014556185287109863,-0.27228591670520086,0.002989987191209683,-0.5489428537951813,-0.02134456783001304,-0.22462281620064825,0.005230853443767429,-0.1894678262257301,0.0033140729457334884,5.295483062326795,-0.2724500716060311,0.026433905662192683,0.01368706308878908,-0.03014264855048227,0.053679001877659956,0.08109477254132096,-0.005004603067203444,0.016907143132293558,0.2105509502082437,0.003657404455024417,-4.904755847017426,0.01634808163992959,-0.008325515865305198,0.062188432751569676,-0.013114633511406406,0.11020519384963083,-0.008599402366091309,-0.012125149710784723,0.31723729052927313,-0.10298398036815914,-0.07250078775612204,0.39976713701763433,0.45897498107347223,0.01018626210400031,0.20163425809089347,0.19729093298588943,0.42863333455911523,0.015595097081693168,0.06060353651437489,-0.16379444813161725,-0.43290344196574165,-0.5931022701412187,1.1906610004748832,0.44418106894148945,0.06536220001548931,0.010261694323554562,-0.05943099382075491,-0.04939614579484797,0.002234505477641322,-0.011262130967449935,0.09644905007708474,-0.029518792883267808,0.41564004027396634,-0.18459770295961597,0.3100981306103734,-0.2509873737065425,0.5434321443668653,0.3009912967350914,1.9560655796099518,-0.7143435150084513,-1.5123449469879784])
        #
        nbOfInputPoints=100;
        f=MEDCouplingFieldDouble.New(ON_NODES_KR,ONE_TIME);
        mesh=MEDCoupling1SGTUMesh.New("aMesh",NORM_POINT1);
        mesh.setCoords(srcPointCoordsXY);
        f.setMesh(mesh);
        f.setArray(srcFieldValsOnPoints);
        f.checkCoherency();
        #
        res0=f.getValueOn([-0.5,-0.5]);
        self.assertAlmostEqual(targetFieldValsExpected.getIJ(0,0),res0[0],10)
        #
        valuesToTest=f.getValueOnMulti(targetPointCoordsXY);
        self.assertEqual(196,valuesToTest.getNumberOfTuples());
        self.assertEqual(1,valuesToTest.getNumberOfComponents());
        for i in xrange(40):
            self.assertAlmostEqual(targetFieldValsExpected[i],valuesToTest.getIJ(i,0),10)
            pass
        fd=f.getDiscretization()
        del f
        self.assertTrue(isinstance(fd,MEDCouplingFieldDiscretizationKriging))
        coeffs,isDrift=fd.computeVectorOfCoefficients(mesh,srcFieldValsOnPoints)
        self.assertEqual(3,isDrift)
        self.assertTrue(coeffsExpected.isEqual(coeffs,1e-8))
        # testing matrix
        pts3=[-0.5,-0.5,-0.5,-0.35,-0.35,-0.2]
        mesh.setCoords(srcPointCoordsXY[:4])
        m,nbCols=fd.computeEvaluationMatrixOnGivenPts(mesh,pts3)
        self.assertTrue(m.isEqual(DataArrayDouble([0.05768877688524917,-4.438982030395039,1.9495386255911573,3.431754627918642,0.11803848510231275,-4.138339658420563,1.6630742187104417,3.357226954607818,0.14630203028580618,-3.5156045565871734,1.414680070737206,2.954622455564169]),1e-12))
        if MEDCouplingHasNumPyBindings():
            import numpy as np
            m0=m.toNumPyArray() ; m0=m0.reshape(3,nbCols) ; m0=np.matrix(m0)
            srcFieldValsOnPoints2=DataArrayDouble(4,2) ; srcFieldValsOnPoints2[:,0]=srcFieldValsOnPoints[:4] ; srcFieldValsOnPoints2[:,1]=2*srcFieldValsOnPoints[:4]
            n0=srcFieldValsOnPoints2.toNumPyArray() ; n0=n0.reshape(4,2) ; n0=np.matrix(n0)
            #
            f=MEDCouplingFieldDouble.New(ON_NODES_KR,ONE_TIME) ;  f.setMesh(mesh) ; f.setArray(srcFieldValsOnPoints2) ; f.checkCoherency()
            self.assertTrue(DataArrayDouble(np.array((m0*n0))).isEqual(f.getValueOnMulti(pts3),1e-14))
            pass
        #
        pass
    
    # test the when input slice is all the same object is return by MEDCouplingMesh.buildPartRange
    def testSwig2MeshPartSlice1(self):
        a=DataArrayDouble(4) ; a.iota()
        c=MEDCouplingCMesh() ; c.setCoords(a,a) ; m=c.buildUnstructured()
        fc0=c.getMeasureField(False) ; fc1=fc0[:] ; fc2=fc0*fc1 ; fc2.setName(fc0.getName())
        self.assertEqual(fc0.getMesh().getHiddenCppPointer(),fc1.getMesh().getHiddenCppPointer())
        self.assertEqual(fc2.getMesh().getHiddenCppPointer(),fc1.getMesh().getHiddenCppPointer())
        self.assertTrue(fc2.isEqual(fc1,1e-12,1e-12))
        #
        fm0=m.getMeasureField(False) ; fm1=fm0[:] ; fm2=fm0*fm1 ; fm2.setName(fm0.getName())
        self.assertEqual(fm0.getMesh().getHiddenCppPointer(),fm1.getMesh().getHiddenCppPointer())
        self.assertEqual(fm2.getMesh().getHiddenCppPointer(),fm1.getMesh().getHiddenCppPointer())
        self.assertTrue(fm2.isEqual(fm1,1e-12,1e-12))
        pass

    # test the correct behaviour when attempting to aggregate two fields whose mesh is null
    def testSwig2MergeFieldsOnFieldsHavingNoMesh(self):
        a=DataArrayDouble(4) ; a.iota() ; a*=1.5
        c=MEDCouplingCMesh() ; c.setCoords(a,a) ; f1=c.getMeasureField(False)
        f1.setMesh(None) ; f2=f1.deepCpy() ; f2*=2
        f3=MEDCouplingFieldDouble.MergeFields(f1,f2)
        daExp=DataArrayDouble([2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5])
        self.assertTrue(f3.getArray().isEqual(daExp,1e-12))
        self.assertEqual(f3.getTypeOfField(),ON_CELLS)
        self.assertEqual(f3.getMesh(),None)
        f4=MEDCouplingFieldDouble.MergeFields([f1,f2])
        self.assertTrue(f4.getArray().isEqual(daExp,1e-12))
        self.assertEqual(f4.getTypeOfField(),ON_CELLS)
        self.assertEqual(f4.getMesh(),None)
        pass

    # test a simple node to cell convertion of a field
    def testSwig2NodeToCellDiscretization1(self):
        f=MEDCouplingFieldDouble(ON_NODES) ; f.setTime(1.1,2,3)
        a1=DataArrayDouble(4) ; a1.iota()
        a2=DataArrayDouble(3) ; a2.iota()
        m=MEDCouplingCMesh() ; m.setCoords(a1,a2)
        f.setMesh(m)
        arr=DataArrayDouble([21.,121.,20.,120.,19.,119.,18.,118.,17.,117.,16.,116.,15.,115.,14.,114.,13.,113.,12.,112.,11.,111.,10.,110.],12,2) ; arr.setInfoOnComponents(["aa [km]","bbb [kJ]"])
        f.setArray(arr) ; f.setName("toto")
        #
        f2=f.nodeToCellDiscretization()
        self.assertEqual(ON_CELLS,f2.getTypeOfField())
        self.assertEqual("toto",f2.getName())
        self.assertEqual([1.1,2,3],f2.getTime())
        self.assertEqual(["aa [km]","bbb [kJ]"],f2.getArray().getInfoOnComponents())
        self.assertEqual(6,f2.getArray().getNumberOfTuples())
        self.assertEqual(f.getMesh().getHiddenCppPointer(),f2.getMesh().getHiddenCppPointer())
        exp=DataArrayDouble([18.5,118.5,17.5,117.5,16.5,116.5,14.5,114.5,13.5,113.5,12.5,112.5],6,2) ; exp.setInfoOnComponents(["aa [km]","bbb [kJ]"])
        self.assertTrue(f2.getArray().isEqual(exp,1e-13))
        pass
    
    def testSwig2NonRegressionBugIntersectMeshes1(self):
        src=MEDCouplingUMesh("src",2)
        src.setCoords(DataArrayDouble([-2.5,-3,-2.5,3,2.5,3],3,2))
        src.allocateCells()
        src.insertNextCell(NORM_TRI3,[0,1,2])
        #
        trg=MEDCouplingUMesh("trg",2)
        trg.setCoords(DataArrayDouble([-2.5,-3.,0.,-3.,0.,-2.,-2.,0.,-2.25,0.,-2.5,0.,-2.5,-1.5,0.,-2.5,-1.25,-3.,-1.414213562373095,-1.414213562373095],10,2))
        trg.allocateCells()
        trg.insertNextCell(NORM_QPOLYG,[2,1,0,5,3,7,8,6,4,9])
        #
        a,b,c=MEDCouplingUMesh.Intersect2DMeshes(src,trg,1.0e-8)
        a.mergeNodes(1e-8)
        self.assertTrue(a.getCoords().isEqual(DataArrayDouble([-2.5,-3.,-2.5,3.,2.5,3.,0.,-3.,0.,-2.,-2.,0.,-2.25,0.,-2.5,0.,-2.5,-1.5,0.,-2.5,-1.25,-3.,-1.414213562373095,-1.414213562373095,-1.2803687993289596,-1.5364425591947515,-1.8901843996644798,-2.2682212795973755,-1.81117884244736,-0.8483107924994473,-2.5,1.5,0.,3.,0.6098156003355202,0.7317787204026243],18,2),1e-12))
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([32,12,0,7,5,13,8,6,14,32,7,1,2,12,5,15,16,17,14,6])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,9,20])))
        self.assertTrue(b.isEqual(DataArrayInt([0,0])))
        self.assertTrue(c.isEqual(DataArrayInt([0,-1])))
        pass

    def testSwig2MeshOrientCorrectly2DCells1(self):
        m=MEDCouplingUMesh("mesh",2)
        coo=DataArrayDouble([1.,0.,0.5,-0.1,0.,1.,0.,0.,0.07,0.5,0.59,0.5],6,2)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_TRI6,[3,0,2,1,5,4])
        m.insertNextCell(NORM_QPOLYG,[3,0,2,1,5,4])
        self.assertTrue(DataArrayDouble([-0.58093333350930543,-0.58093333350930543]).isEqual(m.getMeasureField(False).getArray(),1e-12))
        m.changeSpaceDimension(3)
        m.orientCorrectly2DCells([0.,0.,-1.],False)
        #
        m.checkCoherency()
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([6,3,2,0,4,5,1, 32,3,2,0,4,5,1])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,7,14])))
        m.changeSpaceDimension(2)
        self.assertTrue(DataArrayDouble([0.58093333350930543,0.58093333350930543]).isEqual(m.getMeasureField(False).getArray(),1e-12))
        pass

    def testSwig2Hexa8HavingFacesWarped1(self):
        """ This test is bases on a "error" of interpolation detected. After investigation cell #3 of src is warped that leads to the fact that when trg is 
        intersected with src the sum of intersection volume is greater than the volume of the trg cell.
        A test that can be done is to split the cell #3 of src into tetrohedrons and by summing all the volumes it does not fit the volume computed of cell#3 unsplitted (expect for
        GENERAL_24).
        """
        srcCoo=DataArrayDouble([0.15694071546650565,0.09383333333333337,6.920842121738133,0.15774332475430292,0.185486666666667,6.920682472824616,0.1585459340420992,0.27713999999999994,6.9205228239111,0.07427195882345167,0.05782666666666668,6.937285959830335,0.06343673343819695,0.11347333333333297,6.939441220162809,0.05260150805294228,0.16911999999999996,6.941596480495282,0.014076262238703396,0.04800666666666667,6.949259628344076,0.014076262238703396,0.07092000000000007,6.949259628344076,0.15407499632681992,0.09383333333333338,6.897607484780063,0.15489234394181514,0.18548666666666702,6.897567331066572,0.15570969155680933,0.27714,6.897527177353081,0.06988819198237989,0.05782666666666669,6.901743317269663,0.05885399917995321,0.11347333333333298,6.9022853924017955,0.047819806377526586,0.16912,6.902827467533927,0.0085871208577874,0.048006666666666684,6.9047548457815076,0.0085871208577874,0.07092000000000008,6.9047548457815076,0.153883333333333,0.09383333333333338,6.820902,0.154701666666667,0.18548666666666702,6.820902,0.15551999999999996,0.27714,6.820902,0.06959499999999999,0.05782666666666669,6.820902,0.058547499999999975,0.11347333333333298,6.820902,0.04749999999999999,0.16912,6.820902],22,3)
        src=MEDCouplingUMesh("TBmesh3D",3) ; src.setCoords(srcCoo)
        src.allocateCells()
        src.insertNextCell(NORM_HEXA8,[0,1,4,3,8,9,12,11])
        src.insertNextCell(NORM_HEXA8,[1,2,5,4,9,10,13,12])
        src.insertNextCell(NORM_HEXA8,[4,5,7,6,12,13,15,14])
        src.insertNextCell(NORM_HEXA8,[8,9,12,11,16,17,20,19])
        src.insertNextCell(NORM_HEXA8,[9,10,13,12,17,18,21,20])
        src.checkCoherency2()
        # trg is useless here but I keep it in case of MEDCouplingRemapper were expected to do something about warped NORM_HEXA8
        trgCoo=DataArrayDouble([0.0960891897852753,0.105088620541845,6.8598,0.0599574480546212,0.118434267436059,6.8598,0.113514510609589,0.14874473653263,6.8598,0.0831322609794463,0.167319109733883,6.8598,0.0960891897852753,0.105088620541845,6.92146666666667,0.0599574480546212,0.118434267436059,6.92146666666667,0.113514510609589,0.14874473653263,6.92146666666667,0.0831322609794463,0.167319109733883,6.92146666666667],8,3)
        trg=MEDCouplingUMesh("MESH",3) ; trg.setCoords(trgCoo)
        trg.allocateCells()
        trg.insertNextCell(NORM_HEXA8,[0,1,3,2,4,5,7,6])
        #
        srcFace=src.buildDescendingConnectivity()[0]
        conn=MEDCoupling1SGTUMesh(srcFace).getNodalConnectivity() ; conn.rearrange(4)
        eqFaces=srcFace.computePlaneEquationOf3DFaces()
        nodeIdInCell=3
        e=(srcFace.getCoords()[conn[:,nodeIdInCell]]*eqFaces[:,:-1]).sumPerTuple()+eqFaces[:,3]# e represent the error between the expected 'a*X+b*Y+c*Z+d' in eqFaces and 0. Closer e to 0. is closer the 4th point is to the plane built with the 3 first points
        lambd=-e/(eqFaces[:,:3]**2).sumPerTuple()
        pts=lambd*eqFaces[:,:-1]+srcFace.getCoords()[conn[:,nodeIdInCell]]#pts represent the projection of the last points of each NORM_QUAD4 to the plane defined by the 3 first points of the NORM_QUAD4 cell
        shouldBeZero=(pts*eqFaces[:,:-1]).sumPerTuple()+eqFaces[:,3]# this line is useless only to be sure that pts are on the plane.
        check=(pts-srcFace.getCoords()[conn[:,nodeIdInCell]]).magnitude() # check contains the distance of the last point to its plane
        idsToTest=check.getIdsNotInRange(0.,1e-10)
        self.assertTrue(idsToTest.isEqual(DataArrayInt([17,18,19,20,22,23,24])))
        idsToTest2=idsToTest.getIdsNotInRange(18,22)
        self.assertTrue(idsToTest2.isEqual(DataArrayInt([0,4,5,6])))
        idsToTest2.rearrange(2)
        self.assertTrue(idsToTest2.sumPerTuple().isEqual(DataArrayInt([4,11])))
        pass

    def testSwig2SortHexa8EachOther1(self):
        """
        testing MEDCoupling1SGTUMesh.sortHexa8EachOther method
        """
        coords1=DataArrayDouble([(-0.5,0.5,-0.5),(0.5,-0.5,-0.5),(-0.5,-0.5,0.5),(-0.5,-0.5,-0.5),(0.5,-0.5,0.5),(-0.5,0.5,0.5),(0.5,0.5,0.5),(0.5,0.5,-0.5)])
        m1=MEDCouplingUMesh("m1",3) ; m1.setCoords(coords1)
        m1.allocateCells() ; m1.insertNextCell(NORM_HEXA8,[7,1,3,0,6,4,2,5])
        m1.checkCoherency()
        #
        m2=m1.deepCpy() ; m2.setName("m2")
        #
        trs=[[0.,0.,-1.],[0.,0.,1.],[1.,0.,0.],[0.,-1.,0.],[-1.,0.,0.],[0.,1.,0.]]
        for i,t in enumerate(trs):
            for j in xrange(64):
                j2=(j//16) ; j1=((j%16)//4) ; j0=(j%4)
                m11=m1.deepCpy()
                m11.rotate([0.,0.,0.],[0.,0.,1.],float(j0)*pi/2)
                m11.rotate([0.,0.,0.],[0.,1.,0.],float(j1)*pi/2)
                m11.rotate([0.,0.,0.],[1.,0.,0.],float(j2)*pi/2)
                m11.translate(t)
                #
                m=MEDCouplingUMesh.MergeUMeshes(m2,m11)
                m.mergeNodes(1e-12)
                self.assertEqual(12,m.getNumberOfNodes())
                m=MEDCoupling1SGTUMesh(m)
                m.sortHexa8EachOther()
                tmp0=m.buildUnstructured().tetrahedrize(PLANAR_FACE_6)[0].buildUnstructured()
                self.assertEqual(20,tmp0.computeSkin().getNumberOfCells())
                pass
            pass
        pass

    def testSwig2normMinComputeAbs1(self):
        d=DataArrayDouble([4,-5,2,6.1,-7.33,1,-1,3e2,0.07,-0.009,-6,-1e30],4,3)
        d.setInfoOnComponents(["XX [m]","YYY [km]","ABSJJ [MW]"])
        d0=d.computeAbs()
        dExp=d.deepCpy() ; dExp.abs()
        self.assertTrue(dExp.isEqual(d0,1e-12))
        e=d0-DataArrayDouble([4,5,2,6.1,7.33,1,1,3e2,0.07,0.009,6,1e30],4,3)
        self.assertAlmostEqual(0.,e.normMin(),13)
        self.assertAlmostEqual(0.009,d.normMin(),13)
        #
        di=DataArrayInt([3,-12,5,6,14,16,-23,100,23,-1,0,-6],4,3)
        di.setInfoOnComponents(["XX [m]","YYY [km]","ABSJJ [MW]"])
        d0i=di.computeAbs()
        diExp=di.deepCpy() ; diExp.abs()
        self.assertTrue(diExp.isEqual(d0i))
        self.assertEqual([3,12,5,6,14,16,23,100,23,1,0,6],d0i.getValues())
        pass

    def testSwig2GetCellsContainingPointsForNonConvexPolygon1(self):
        coo=DataArrayDouble([-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,-0.5,0.,-0.5,0.,0.,0.5,0.,],7,2)
        m=MEDCouplingUMesh("Intersect2D",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_POLYGON,[6,3,4,5])
        m.insertNextCell(NORM_POLYGON,[4,0,1,2,6,5])
        m.checkCoherency2()
        #
        self.assertTrue(m.getCellsContainingPoint((0.4,-0.4),1e-12).isEqual(DataArrayInt([0])))
        self.assertTrue(m.getCellsContainingPoint((-0.4,-0.4),1e-12).isEqual(DataArrayInt([1])))
        self.assertTrue(m.getCellsContainingPoint((0.,-0.4),1e-12).isEqual(DataArrayInt([0,1])))
        pass
    
    def testSwig2GetCellsContainingPointsForNonConvexPolygon2(self):
        coo=DataArrayDouble([-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,-0.5,-2.0816681711721685e-17,-2.0816681711721685e-17,-0.17677669529663687,0.1767766952966369,0.,0.5,0.5,0.,0.17677669529663684,-0.17677669529663692,0.17677669529663692,0.17677669529663684,-0.17677669529663692,-0.17677669529663687,0.,-0.5,-0.5,0.,0.33838834764831843,-0.3383883476483185,-0.33838834764831843,0.33838834764831843,-0.21213203435596423,0.21213203435596426,0.2121320343559642,-0.2121320343559643,0.21213203435596426,0.2121320343559642,-0.21213203435596423,-0.21213203435596428,0.3560660171779821,-0.35606601717798214,-0.35606601717798214,0.35606601717798214,0.19445436482630052,-0.19445436482630063,-0.19445436482630055,0.19445436482630057,0.,0.27],24,2)
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_QPOLYG,[8,5,4,9])
        m.insertNextCell(NORM_QPOLYG,[5,8,4,10])
        m.insertNextCell(NORM_QPOLYG,[16,8,5,15,21,9,22,17])
        m.insertNextCell(NORM_QPOLYG,[15,1,2,3,16,20,6,7,19,17])
        m.insertNextCell(NORM_QPOLYG,[15,5,8,16,22,10,21,18])
        m.insertNextCell(NORM_QPOLYG,[16,3,0,1,15,19,11,12,20,18])
        m.checkCoherency2()
        self.assertTrue(m.getCellsContainingPoint([0.,0.27],1e-12).isEqual(DataArrayInt([2])))
        pass

    def testSwig2DAIGetIdsEqualTuple1(self):
        da=DataArrayInt([0,7,1,2,4,1,2,1,1,2,0,1,2,1,5,1,1,2],9,2)
        self.assertTrue(da.getIdsEqualTuple([1,2]).isEqual(DataArrayInt([1,4,8])))
        self.assertTrue(da.getIdsEqualTuple((1,2)).isEqual(DataArrayInt([1,4,8])))
        self.assertTrue(da.getIdsEqualTuple(DataArrayInt([1,2])).isEqual(DataArrayInt([1,4,8])))
        da.rearrange(3)
        self.assertRaises(InterpKernelException,da.getIdsEqualTuple,[1,2])# mismatch nb of compo (3) and nb of elts in input tuple (2)
        self.assertTrue(da.getIdsEqualTuple([2,0,1]).isEqual(DataArrayInt([3])))
        self.assertTrue(da.getIdsEqualTuple([2,0,7]).isEqual(DataArrayInt([])))
        da.rearrange(1)
        self.assertTrue(da.getIdsEqualTuple(2).isEqual(DataArrayInt([3,6,9,12,17])))
        self.assertTrue(da.getIdsEqualTuple(2).isEqual(da.getIdsEqual(2)))
        pass

    def testSwig2GaussNEStaticInfo1(self):
        self.assertTrue(DataArrayDouble(MEDCouplingFieldDiscretizationGaussNE.GetWeightArrayFromGeometricType(NORM_TRI3)).isEqual(DataArrayDouble([0.16666666666666666,0.16666666666666666,0.16666666666666666]),1e-12))
        self.assertTrue(DataArrayDouble(MEDCouplingFieldDiscretizationGaussNE.GetRefCoordsFromGeometricType(NORM_TRI3)).isEqual(DataArrayDouble([0.,0.,1.,0.,0.,1.]),1e-12))
        self.assertTrue(DataArrayDouble(MEDCouplingFieldDiscretizationGaussNE.GetLocsFromGeometricType(NORM_TRI3)).isEqual(DataArrayDouble([0.16666666666666666,0.16666666666666666,0.6666666666666667,0.16666666666666666,0.16666666666666666,0.6666666666666667]),1e-12))
        pass

    def testSwigReverseNodalConnOnStructuredMesh(self):
        # 1D - standard
        c=MEDCouplingCMesh() ; arr=DataArrayDouble(10) ; arr.iota()
        c.setCoordsAt(0,arr)
        rn,rni=c.getReverseNodalConnectivity()
        rn2,rni2=c.buildUnstructured().getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(DataArrayInt([0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8])))
        self.assertTrue(rni.isEqual(DataArrayInt([0,1,3,5,7,9,11,13,15,17,18])))
        self.assertTrue(rn.isEqual(rn2)) ; self.assertTrue(rni.isEqual(rni2))
        # 1D - limit
        c=MEDCouplingCMesh() ; arr=DataArrayDouble(1) ; arr.iota()
        c.setCoordsAt(0,arr)
        rn,rni=c.getReverseNodalConnectivity()
        rn2,rni2=c.buildUnstructured().getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(DataArrayInt([0])))
        self.assertTrue(rni.isEqual(DataArrayInt([0,1])))
        self.assertTrue(rn.isEqual(rn2)) ; self.assertTrue(rni.isEqual(rni2))
        # 1D - limit
        c=MEDCouplingCMesh() ; arr=DataArrayDouble(0) ; arr.iota()
        c.setCoordsAt(0,arr)
        rn,rni=c.getReverseNodalConnectivity()
        rn.isEqual(DataArrayInt([]))
        rni.isEqual(DataArrayInt([0]))
        # 2D - standard
        c=MEDCouplingCMesh() ; arr=DataArrayDouble(5) ; arr.iota() ; arr2=DataArrayDouble(4) ; arr.iota()
        c.setCoords(arr,arr2)
        rn,rni=c.getReverseNodalConnectivity()
        rn2,rni2=c.buildUnstructured().getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(DataArrayInt([0,0,1,1,2,2,3,3,0,4,0,1,4,5,1,2,5,6,2,3,6,7,3,7,4,8,4,5,8,9,5,6,9,10,6,7,10,11,7,11,8,8,9,9,10,10,11,11])))
        self.assertTrue(rni.isEqual(DataArrayInt([0,1,3,5,7,8,10,14,18,22,24,26,30,34,38,40,41,43,45,47,48])))
        self.assertTrue(rn.isEqual(rn2)) ; self.assertTrue(rni.isEqual(rni2))
        # 2D - limit
        c=MEDCouplingCMesh() ; arr=DataArrayDouble(10) ; arr.iota() ; arr2=DataArrayDouble(1) ; arr.iota()
        c.setCoords(arr,arr2)
        rn,rni=c.getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(DataArrayInt([0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8])))
        self.assertTrue(rni.isEqual(DataArrayInt([0,1,3,5,7,9,11,13,15,17,18])))
        # 2D - limit
        c=MEDCouplingCMesh() ; arr=DataArrayDouble(10) ; arr.iota() ; arr2=DataArrayDouble(1) ; arr.iota()
        c.setCoords(arr2,arr)
        rn,rni=c.getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(DataArrayInt([0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8])))
        self.assertTrue(rni.isEqual(DataArrayInt([0,1,3,5,7,9,11,13,15,17,18])))
        # 3D - standard
        c=MEDCouplingCMesh() ; arr0=DataArrayDouble(5) ; arr0.iota() ; arr1=DataArrayDouble(3) ; arr1.iota() ; arr2=DataArrayDouble(4) ; arr2.iota()
        c.setCoords(arr0,arr1,arr2)
        rn,rni=c.getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(DataArrayInt([0,0,1,1,2,2,3,3,0,4,0,1,4,5,1,2,5,6,2,3,6,7,3,7,4,4,5,5,6,6,7,7,0,8,0,1,8,9,1,2,9,10,2,3,10,11,3,11,0,4,8,12,0,1,4,5,8,9,12,13,1,2,5,6,9,10,13,14,2,3,6,7,10,11,14,15,3,7,11,15,4,12,4,5,12,13,5,6,13,14,6,7,14,15,7,15,8,16,8,9,16,17,9,10,17,18,10,11,18,19,11,19,8,12,16,20,8,9,12,13,16,17,20,21,9,10,13,14,17,18,21,22,10,11,14,15,18,19,22,23,11,15,19,23,12,20,12,13,20,21,13,14,21,22,14,15,22,23,15,23,16,16,17,17,18,18,19,19,16,20,16,17,20,21,17,18,21,22,18,19,22,23,19,23,20,20,21,21,22,22,23,23])))
        self.assertTrue(rni.isEqual(DataArrayInt([0,1,3,5,7,8,10,14,18,22,24,25,27,29,31,32,34,38,42,46,48,52,60,68,76,80,82,86,90,94,96,98,102,106,110,112,116,124,132,140,144,146,150,154,158,160,161,163,165,167,168,170,174,178,182,184,185,187,189,191,192])))
        rn2,rni2=c.buildUnstructured().getReverseNodalConnectivity()
        self.assertTrue(rn.isEqual(rn2)) ; self.assertTrue(rni.isEqual(rni2))
        pass

    def testSwig2CellToNodeDiscretization1(self):
        m=MEDCouplingCMesh() ; arr0=DataArrayDouble(5) ; arr0.iota() ; arr1=DataArrayDouble(4) ; arr1.iota() ; m.setCoords(arr0,arr1)
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(m) ; f.setTime(1.1,5,6)
        arr=DataArrayDouble(12) ; arr.iota()
        arr=DataArrayDouble.Meld(arr,arr+100.) ; arr.setInfoOnComponents(["aaa","bbb"])
        f.setArray(arr)
        f.checkCoherency()
        #
        ref=DataArrayDouble([0.,0.5,1.5,2.5,3.,2.,2.5,3.5,4.5,5.,6.,6.5,7.5,8.5,9.,8.,8.5,9.5,10.5,11.])
        ref=DataArrayDouble.Meld(ref,ref+100.) ; ref.setInfoOnComponents(["aaa","bbb"])
        f2=f.cellToNodeDiscretization()
        f2.checkCoherency()
        self.assertEqual(f2.getTime()[1:],[5,6])
        self.assertAlmostEqual(f2.getTime()[0],1.1,15)
        self.assertEqual(f2.getMesh().getHiddenCppPointer(),m.getHiddenCppPointer())
        self.assertTrue(f2.getArray().isEqual(ref,1e-12))
        rn,rni=m.getReverseNodalConnectivity()
        rni2=(rni.deltaShiftIndex()).convertToDblArr()
        arr2=(f.getArray()[rn]).accumulatePerChunck(rni)/rni2
        self.assertTrue(f2.getArray().isEqual(arr2,1e-12))
        del f2
        #
        u=m.buildUnstructured() ; f.setMesh(u) ; del m
        f3=f.cellToNodeDiscretization()
        f3.checkCoherency()
        self.assertEqual(f3.getTime()[1:],[5,6])
        self.assertAlmostEqual(f3.getTime()[0],1.1,15)
        self.assertEqual(f3.getMesh().getHiddenCppPointer(),u.getHiddenCppPointer())
        self.assertTrue(f3.getArray().isEqual(ref,1e-12))
        pass

    def testSwig2GetMeshSpaceDimensionCMesh1(self):
        c=MEDCouplingCMesh()
        arr0=DataArrayDouble([0,1,2])
        arr1=DataArrayDouble([0])
        c.setCoords(arr0,arr0,arr0)
        self.assertEqual(c.getMeshDimension(),3)
        self.assertEqual(c.getSpaceDimension(),3)
        #
        c.setCoords(arr0,arr0,arr1)
        self.assertEqual(c.getMeshDimension(),2)
        self.assertEqual(c.getSpaceDimension(),3)
        #
        c.setCoords(arr0,arr0)
        self.assertEqual(c.getMeshDimension(),2)
        self.assertEqual(c.getSpaceDimension(),2)
        #
        c.setCoords(arr0,arr1)
        self.assertEqual(c.getMeshDimension(),1)
        self.assertEqual(c.getSpaceDimension(),2)
        #
        c.setCoords(arr0)
        self.assertEqual(c.getMeshDimension(),1)
        self.assertEqual(c.getSpaceDimension(),1)
        #
        c.setCoords(arr1)
        self.assertEqual(c.getMeshDimension(),0)
        self.assertEqual(c.getSpaceDimension(),1)
        pass

    def testSwig2BuildSpreadZonesWithPolyOnQPolyg1(self):
        nx=6
        ny=6
        m=MEDCouplingCMesh()
        arr1=DataArrayDouble(nx) ; arr1.iota()
        arr2=DataArrayDouble(ny) ; arr2.iota()
        m.setCoords(arr1,arr2)
        m=m.buildUnstructured()
        da=DataArrayInt.Range(nx-1,(nx-1)*(ny-1),nx)
        m2=m[da] ; m2.simplexize(0)
        dan=da.buildComplement(m.getNumberOfCells())
        m1=m[dan]
        m=MEDCouplingUMesh.MergeUMeshesOnSameCoords(m1,m2)
        #
        m.convertLinearCellsToQuadratic()
        m1=m[::2] ; m2=m[1::2] ; m2.convertAllToPoly()
        m=MEDCouplingUMesh.MergeUMeshesOnSameCoords(m1,m2)
        p=m.buildSpreadZonesWithPoly()
        self.assertTrue(p.getNodalConnectivity().isEqual(DataArrayInt([32,1,0,6,12,18,24,30,31,32,33,34,35,29,23,17,11,5,4,3,2,36,37,94,62,72,83,84,86,89,99,92,93,82,71,60,51,49,46,43,40])))
        self.assertTrue(p.getNodalConnectivityIndex().isEqual(DataArrayInt([0,41])))
        self.assertTrue(p.getCoords().isEqual(DataArrayDouble([0.,0.,1.,0.,2.,0.,3.,0.,4.,0.,5.,0.,0.,1.,1.,1.,2.,1.,3.,1.,4.,1.,5.,1.,0.,2.,1.,2.,2.,2.,3.,2.,4.,2.,5.,2.,0.,3.,1.,3.,2.,3.,3.,3.,4.,3.,5.,3.,0.,4.,1.,4.,2.,4.,3.,4.,4.,4.,5.,4.,0.,5.,1.,5.,2.,5.,3.,5.,4.,5.,5.,5.,0.5,0.,0.,0.5,0.5,1.,1.,0.5,1.5,0.,1.5,1.,2.,0.5,2.5,0.,2.5,1.,3.,0.5,3.5,0.,3.5,1.,4.,0.5,4.5,0.,4.5,1.,5.,0.5,1.,1.5,1.5,2.,2.,1.5,2.5,2.,3.,1.5,3.5,2.,4.,1.5,4.5,2.,5.,1.5,0.5,2.,0.,2.5,0.5,3.,1.,2.5,2.,2.5,2.5,3.,3.,2.5,3.5,3.,4.,2.5,4.5,3.,5.,2.5,0.,3.5,0.5,4.,1.,3.5,1.5,3.,1.5,4.,2.,3.5,3.,3.5,3.5,4.,4.,3.5,4.5,4.,5.,3.5,0.,4.5,0.5,5.,1.,4.5,1.5,5.,2.,4.5,2.5,4.,2.5,5.,3.,4.5,4.,4.5,4.5,5.,5.,4.5,0.,1.5,0.5,1.5,1.5,2.5,2.5,3.5,3.5,4.5,3.5,5.0],100,2),1e-13))
        pass

    def testSwig2Conformize2D1(self):
        eps = 1.0e-8
        coo = [0.,-0.5,0.,0.,0.5,0.,0.5,-0.5,0.25,
               -0.1,0.25,0.,0.5,-0.1,0.,0.5,0.5,0.5,0.25,0.4,0.25,0.5,0.5,0.4]
        conn = [5,5,2,6,4,5,6,3,0,1,5,4,5,10,8,11,9,5,11,2,1,7,10,9]
        connI = [0,5,12,17,24]
        m = MEDCouplingUMesh("box",2)
        cooArr = DataArrayDouble(coo,len(coo)/2,2)
        m.setCoords(cooArr)
        m.setConnectivity(DataArrayInt(conn),DataArrayInt(connI))
        m.mergeNodes(eps)
        m.checkCoherency()
        self.assertTrue(m.conformize2D(eps).isEqual(DataArrayInt([3])))
        self.assertEqual(m.getCoords().getHiddenCppPointer(),cooArr.getHiddenCppPointer()) # check that coordinates remain the same here
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,5,2,6,4,5,6,3,0,1,5,4,5,10,8,11,9,5,11,2,5,1,7,10,9])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,12,17,25])))
        pass

    def testSwig2Conformize2D2(self):
        eps = 1.0e-8
        coo=DataArrayDouble([-10,-6,0,-6,0,0,7,0,-10,2,0,2,0,6,7,6,0,8,7,8,-10,12,-4,12,0,12,0,11,7,11,-4,16,0,16,7,16],18,2)
        conn=DataArrayInt([2,3,7,6, 13,16,17,14, 4,10,12,5, 9,14,13,8, 8,9,7,6, 5,4,0,1, 16,12,11,15])
        m=MEDCoupling1SGTUMesh("mesh",NORM_QUAD4)
        m.setCoords(coo)
        m.setNodalConnectivity(conn)
        m=m.buildUnstructured()
        self.assertTrue(m.conformize2D(eps).isEqual(DataArrayInt([0,1,2,5])))
        self.assertEqual(m.getCoords().getHiddenCppPointer(),coo.getHiddenCppPointer()) # check that coordinates remain the same here
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,2,3,7,6,5, 5,13,12,16,17,14, 5,4,10,11,12,13,8,6,5, 4,9,14,13,8, 4,8,9,7,6, 5,5,4,0,1,2, 4,16,12,11,15])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,6,12,21,26,31,37,42])))
        pass

    def testSwigSplit2DCells1(self):
        coo=DataArrayDouble([[0,0],[1,0],[1,1],[0,1],[0.5,0],[1,0.5],[0.5,1],[0.,0.5]])
        m=MEDCouplingUMesh("mesh",2)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_QUAD8,[0,1,2,3,4,5,6,7])
        _,d,di,_,_=m.buildDescendingConnectivity()
        subb=DataArrayInt([5])
        subbi=DataArrayInt([0,0,1,1,1])
        mid=DataArrayInt([-1,-1])
        midi=DataArrayInt([0,0,2,2,2])
        self.assertEqual(2,m.split2DCells(d,di,subb,subbi,mid,midi))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([32,0,1,5,2,3,4,8,9,6,7])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,11])))
        self.assertTrue(m.getCoords().isEqual(DataArrayDouble([[0,0],[1,0],[1,1],[0,1],[0.5,0],[1,0.5],[0.5,1],[0.,0.5],[1.,0.25],[1.,0.75]]),1e-12))
        pass

    def testSwig2Conformize2D3(self):
        eps = 1.0e-8
        coo=DataArrayDouble([-10,-6,0,-6,0,0,7,0,-10,2,0,2,0,6.5,7,6.5,0,8,7,8,-10,12,-4,12,0,12,0,11,7,11,-4,16,0,16,7,16],18,2)
        conn=DataArrayInt([2,3,7,6, 13,16,17,14, 4,10,12,5, 9,14,13,8, 8,9,7,6, 5,4,0,1, 16,12,11,15])
        m=MEDCoupling1SGTUMesh("mesh",NORM_QUAD4)
        m.setCoords(coo)
        m.setNodalConnectivity(conn)
        m=m.buildUnstructured()
        m.convertLinearCellsToQuadratic()
        self.assertTrue(m.conformize2D(eps).isEqual(DataArrayInt([0,1,2,5])))
        self.assertTrue(m.getCoords().getHiddenCppPointer()!=coo.getHiddenCppPointer()) # coordinates are not the same here contrary to testSwig2Conformize2D2 ...
        self.assertTrue(m.getCoords()[:18].isEqual(coo,1e-12)) # but the 18 first nodes are the same
        pass

    def testSwig2Conformize2D4(self):
        eps = 1.0e-8
        coo=DataArrayDouble([-10,-6,0,-6,0,0,7,0,-10,2,0,2,0,6.5,7,6.5,0,8,7,8,-10,12,-4,12,0,12,0,11,7,11,-4,16,0,16,7,16],18,2)
        conn=DataArrayInt([2,3,7,6, 13,16,17,14, 4,10,12,5, 9,14,13,8, 8,9,7,6, 5,4,0,1, 16,12,11,15])
        m=MEDCoupling1SGTUMesh("mesh",NORM_QUAD4)
        m.setCoords(coo)
        m.setNodalConnectivity(conn)
        m=m.buildUnstructured()
        m.convertLinearCellsToQuadratic()
        self.assertEqual(42,m.getNumberOfNodes())
        oldCoo=m.getCoords().deepCpy()
        m.conformize2D(eps)
        self.assertTrue(m.getCoords()[:42].isEqual(oldCoo,1e-12))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([32,2,3,7,6,5,18,19,20,42,43,32,13,12,16,17,14,44,38,23,24,25,32,4,10,11,12,13,8,6,5,26,45,39,44,31,34,42,29,8,9,14,13,8,30,25,31,32,8,8,9,7,6,32,33,20,34,32,5,4,0,1,2,29,35,36,46,43,8,16,12,11,15,38,39,40,41])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,11,22,39,48,57,68,77])))
        self.assertTrue(m.getCoords().isEqual(DataArrayDouble([[-10.,-6.0],[0.,-6.0],[0.,0.0],[7.,0.0],[-10.,2.0],[0.,2.0],[0.,6.5],[7.,6.5],[0.,8.0],[7.,8.0],[-10.,12.0],[-4.,12.0],[0.,12.0],[0.,11.0],[7.,11.0],[-4.,16.0],[0.,16.0],[7.,16.0],[3.5, 0.0],[7.,3.25],[3.5, 6.5],[0.,3.25],[0.,13.5],[3.5, 16.0],[7.,13.5],[3.5, 11.0],[-10.,7.0],[-5.,12.0],[0.,7.0],[-5.,2.0],[7.,9.5],[0.,9.5],[3.5, 8.0],[7.,7.25],[0.,7.25],[-10.,-2.0],[-5.,-6.0],[0.,-2.0],[0.,14.0],[-2.,12.0],[-4.,14.0],[-2.,16.0],[0.,4.25],[0.,1.0],[0.,11.5],[-7.,12.0],[0.,-3.]]),1e-12))
        pass

    def testSwig2Conformize2D5(self):
        eps=1e-8
        coo=DataArrayDouble([[2,2],[2,-6],[10,-2],[-2,-2],[6,0],[6,-4],[2,7],[2,4.5],[-1.4641016151377544,0],[-1.950753362380551,-1.3742621398390762],[-7,-3],[-0.8284271247461898,-4.82842712474619],[0.26794919243112281,3.5],[0,1.4641016151377548],[-4.4753766811902755,-2.1871310699195381],[-3.9142135623730949,-3.9142135623730949],[-1.8042260651806146,-3.23606797749979]])
        m=MEDCouplingUMesh("mesh",2)
        m.allocateCells()
        m.setCoords(coo)
        m.insertNextCell(NORM_TRI6,[1,2,0,5,4,3])
        m.insertNextCell(NORM_TRI6,[8,6,0,12,7,13])
        m.insertNextCell(NORM_TRI6,[11,9,10,16,14,15])
        self.assertTrue(m.conformize2D(eps).isEqual(DataArrayInt([0])))
        self.assertTrue(m.getCoords().isEqual(DataArrayDouble([2.,2.,2.,-6.,10.,-2.,-2.,-2.,6.,0.,6.,-4.,2.,7.,2.,4.5,-1.4641016151377544,0.,-1.950753362380551,-1.3742621398390762,-7.,-3.,-0.8284271247461898,-4.82842712474619,0.2679491924311228,3.5,8.881784197001252e-16,1.4641016151377548,-4.4753766811902755,-2.187131069919538,-3.914213562373095,-3.914213562373095,-1.8042260651806146,-3.236067977499789,-1.7705659643687133,-0.6647725630649153,0.46926627053963865,-5.695518130045146],19,2),1e-12))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([32,1,2,0,8,9,11,5,4,13,17,16,18,6,8,6,0,12,7,13,6,11,9,10,16,14,15])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,13,20,27])))
        pass

    def testSwigExtendedSlice1(self):
        d=DataArrayInt([5,6,7])
        self.assertTrue(d[2:].isEqual(DataArrayInt([7])))
        self.assertTrue(d[3:].isEqual(DataArrayInt([])))
        try:
            d[4:]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        d=DataArrayInt([5,6,7,8])
        self.assertEqual(d[-1],8)
        self.assertEqual(d[-4],5)
        try:
            d[-5]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        self.assertTrue(d[2::-1].isEqual(DataArrayInt([7,6,5])))
        self.assertTrue(d[0::-1].isEqual(DataArrayInt([5])))
        self.assertTrue(d[-1::-1].isEqual(DataArrayInt([8,7,6,5])))
        self.assertTrue(d[-3::-1].isEqual(DataArrayInt([6,5])))
        self.assertTrue(d[-5::-1].isEqual(DataArrayInt([])))
        try:
            d[-6::-1]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        d=DataArrayInt([])
        self.assertTrue(d[0:].isEqual(DataArrayInt([])))
        #
        d=DataArrayDouble([5,6,7])
        self.assertTrue(d[2:].isEqual(DataArrayDouble([7]),1e-12))
        self.assertTrue(d[3:].isEqual(DataArrayDouble([]),1e-12))
        try:
            d[4:]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        d=DataArrayDouble([5,6,7,8])
        self.assertAlmostEqual(d[-1],8.,12)
        self.assertAlmostEqual(d[-4],5.,12)
        try:
            d[-5]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        self.assertTrue(d[2::-1].isEqual(DataArrayDouble([7,6,5]),1e-12))
        self.assertTrue(d[0::-1].isEqual(DataArrayDouble([5]),1e-12))
        self.assertTrue(d[-1::-1].isEqual(DataArrayDouble([8,7,6,5]),1e-12))
        self.assertTrue(d[-3::-1].isEqual(DataArrayDouble([6,5]),1e-12))
        self.assertTrue(d[-5::-1].isEqual(DataArrayDouble([]),1e-12))
        try:
            d[-6::-1]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
        d=DataArrayDouble([])
        self.assertTrue(d[0:].isEqual(DataArrayDouble([]),1e-12))
        pass

    def testSwig2Hexa27GP1(self):
        """ This test focused on shape functions of hexa27.
        """
        coo=DataArrayDouble([[0.,2.,2.],[0.,0.,2.],[2.,0.,2.],[2.,2.,2.],[0.,2.,0.],[0.,0.,0.],[2.,0.,0.],[2.,2.,0.], [0.,1.,2.],[1.,0.,2.],[2.,1.,2.],[1.,2.,2.], [0.,1.,0.],[1.,0.,0.],[2.,1.,0.],[1.,2.,0.], [0.,2.,1.],[0.,0.,1.],[2.,0.,1.],[2.,2.,1.], [1.,1.,2.], [0.,1.,1.],[1.,0.,1.],[2.,1.,1.],[1.,2.,1.], [1.,1.,0.], [1.,1.,1.]])
        m=MEDCouplingUMesh("mesh",3) ; m.setCoords(coo)
        m.allocateCells()
        # the cell description is exactly those described in the description of HEXA27 in MED file 3.0.7 documentation
        m.insertNextCell(NORM_HEXA27,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26])
        refCoo=[-1.,-1.,-1.,-1.,1.,-1.,1.,1.,-1.,1.,-1.,-1.,-1.,-1.,1.,-1.,1.,1.,1.,1.,1.,1.,-1.,1.,-1.,0.,-1.,0.,1.,-1.,1.,0.,-1.,0.,-1.,-1.,-1.,0.,1.,0.,1.,1.,1.,0.,1.,0.,-1.,1.,-1.,-1.,0.,-1.,1.,0.,1.,1.,0.,1.,-1.,0.,0.,0.,-1.,-1.,0.,0.,0.,1.,0.,1.,0.,0.,0.,-1.,0.,0.,0.,1.,0.,0.,0.]
        weights=[0.1714677640603571,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.43895747599451346,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.43895747599451346,0.27434842249657115,0.43895747599451346,0.7023319615912209,0.43895747599451346,0.27434842249657115,0.43895747599451346,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.43895747599451346,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.1714677640603571]
        gCoords=[-0.774596669241483,-0.774596669241483,-0.774596669241483,-0.774596669241483,-0.774596669241483,0.0,-0.774596669241483,-0.774596669241483,0.774596669241483,-0.774596669241483,0.0,-0.774596669241483,-0.774596669241483,0.0,0.0,-0.774596669241483,0.0,0.774596669241483,-0.774596669241483,0.774596669241483,-0.774596669241483,-0.774596669241483,0.774596669241483,0.0,-0.774596669241483,0.774596669241483,0.774596669241483,0.0,-0.774596669241483,-0.774596669241483,0.0,-0.774596669241483,0.0,0.0,-0.774596669241483,0.774596669241483,0.0,0.0,-0.774596669241483,0.0,0.0,0.0,0.0,0.0,0.774596669241483,0.0,0.774596669241483,-0.774596669241483,0.0,0.774596669241483,0.0,0.0,0.774596669241483,0.774596669241483,0.774596669241483,-0.774596669241483,-0.774596669241483,0.774596669241483,-0.774596669241483,0.0,0.774596669241483,-0.774596669241483,0.774596669241483,0.774596669241483,0.0,-0.774596669241483,0.774596669241483,0.0,0.0,0.774596669241483,0.0,0.774596669241483,0.774596669241483,0.774596669241483,-0.774596669241483,0.774596669241483,0.774596669241483,0.0,0.774596669241483,0.774596669241483,0.774596669241483]
        fGauss=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fGauss.setName("fGauss")
        fGauss.setMesh(m)
        fGauss.setGaussLocalizationOnType(NORM_HEXA27,refCoo,gCoords,weights)
        arr=DataArrayDouble(fGauss.getNumberOfTuplesExpected()) ; arr.iota()
        fGauss.setArray(arr)
        arrOfDisc=fGauss.getLocalizationOfDiscr()
        # the test is here
        self.assertTrue(arrOfDisc.isEqual(DataArrayDouble([0.2254033307585172,1.7745966692414836,1.7745966692414834,0.22540333075851715,1.7745966692414834,1.,0.22540333075851715,1.7745966692414836,0.22540333075851715,0.22540333075851715,1.,1.7745966692414834,0.2254033307585171,1.,1.,0.22540333075851715,1.0000000000000002,0.2254033307585171,0.22540333075851715,0.22540333075851715,1.7745966692414838,0.22540333075851715,0.22540333075851715,1.,0.22540333075851715,0.22540333075851715,0.22540333075851715,1.,1.7745966692414832,1.7745966692414834,1.,1.774596669241483,1.,1.0000000000000002,1.7745966692414832,0.22540333075851712,1.,1.,1.774596669241483,1.,1.,1.,1.,1.,0.2254033307585171,1.,0.22540333075851715,1.7745966692414834,1.,0.2254033307585171,1.,1.0000000000000002,0.22540333075851715,0.2254033307585171,1.7745966692414834,1.7745966692414834,1.7745966692414836,1.7745966692414832,1.7745966692414834,1.0000000000000002,1.7745966692414834,1.7745966692414836,0.22540333075851712,1.7745966692414832,1.,1.7745966692414834,1.774596669241483,1.,1.,1.7745966692414832,1.0000000000000002,0.22540333075851712,1.7745966692414836,0.22540333075851715,1.7745966692414836,1.7745966692414832,0.22540333075851715,1.,1.7745966692414836,0.22540333075851715,0.22540333075851715],27,3),1e-12))
        #
        weights=27*[1]
        gCoords=refCoo
        fGauss.setGaussLocalizationOnType(NORM_HEXA27,refCoo,gCoords,weights)
        arrOfDisc2=fGauss.getLocalizationOfDiscr()
        self.assertTrue(arrOfDisc2.isEqual(coo,1e-12))
        pass

    def testSwig2Pyra13GP1(self):
        coo=DataArrayDouble([[0.,2.,0.],[2.,2.,0.],[2.,0.,0.],[0.,0.,0.],[1.,1.,2.],[1.,2.,0.],[2.,1.,0.],[1.,0.,0.],[0.,1.,0.],[0.5,1.5,1.],[1.5,1.5,1.],[1.5,0.5,1.],[0.5,0.5,1.]])
        m=MEDCouplingUMesh("mesh",3) ; m.setCoords(coo)
        m.allocateCells()
        # the cell description is exactly those described in the description of PYRA13 in MED file 3.0.7 documentation
        m.insertNextCell(NORM_PYRA13,[0,1,2,3,4,5,6,7,8,9,10,11,12])
        refCoords=[1.,0.,0.,0.,-1.,0.,-1.,0.,0.,0.,1.,0.,0.,0.,1.,0.5,-0.5,0.,-0.5,-0.5,0.,-0.5,0.5,0.,0.5,0.5,0.,0.5,0.,0.5,0.,-0.5,0.5,-0.5,0.,0.5,0.,0.5,0.5]
        gaussCoords=[0.,0.,0.5,0.21210450275,0.21210450275,0.5,-0.21210450275,0.21210450275,0.5,-0.21210450275,-0.21210450275,0.5,0.21210450275,-0.21210450275,0.5,0.,0.,0.07579099449999999,0.,0.,0.9242090055000001,0.5394929090572634,0.,0.17359176399999998,0.,0.5394929090572634,0.17359176399999998,-0.5394929090572634,0.,0.17359176399999998,0.,-0.5394929090572634,0.17359176399999998,0.1133235629427366,0.,0.826408236,0.,0.1133235629427366,0.826408236,-0.1133235629427366,0.,0.826408236,0.,-0.1133235629427366,0.826408236,0.5826406005183961,0.5826406005183961,-0.053206449499999975,-0.5826406005183961,0.5826406005183961,-0.053206449499999975,-0.5826406005183961,-0.5826406005183961,-0.053206449499999975,0.5826406005183961,-0.5826406005183961,-0.053206449499999975,0.5532064495,0.,0.5,0.,0.5532064495,0.5,-0.5532064495,0.,0.5,0.,-0.5532064495,0.5,-0.029434151018396033,-0.029434151018396033,1.0532064495,0.029434151018396033,-0.029434151018396033,1.0532064495,0.029434151018396033,0.029434151018396033,1.0532064495,-0.029434151018396033,0.029434151018396033,1.0532064495]
        weights=[0.0492545926875,0.031210562625,0.031210562625,0.031210562625,0.031210562625,0.10663554205740113,0.0007171281994273535,0.0816994048010844,0.0816994048010844,0.0816994048010844,0.0816994048010844,0.0036048554264914074,0.0036048554264914074,0.0036048554264914074,0.0036048554264914074,0.008958181586640837,0.008958181586640837,0.008958181586640837,0.008958181586640837,0.002018983875,0.002018983875,0.002018983875,0.002018983875,2.286237794882217e-05,2.286237794882217e-05,2.286237794882217e-05,2.286237794882217e-05]
        fGauss=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fGauss.setName("fGauss")
        fGauss.setMesh(m)
        fGauss.setGaussLocalizationOnType(NORM_PYRA13,refCoords,gaussCoords,weights)
        arr=DataArrayDouble(fGauss.getNumberOfTuplesExpected()) ; arr.iota()
        fGauss.setArray(arr)
        arrOfDisc=fGauss.getLocalizationOfDiscr()
        # the test is here
        self.assertTrue(arrOfDisc.isEqual(DataArrayDouble([1.,1.,1.,0.5757909945,1.,1.,1.,0.5757909945,1.,1.4242090055,1.,1.,1.,1.4242090055,1.,1.,1.,0.151581989,1.,1.,1.848418011,0.4605070909427367,1.5394929090572635,0.347183528,0.4605070909427367,0.4605070909427367,0.347183528,1.5394929090572638,0.4605070909427366,0.347183528,1.5394929090572635,1.5394929090572638,0.347183528,0.8866764370572636,1.1133235629427367,1.652816472,0.8866764370572636,0.8866764370572636,1.652816472,1.1133235629427367,0.8866764370572636,1.652816472,1.1133235629427365,1.1133235629427367,1.652816472,-0.16528120103679209,1.,-0.106412899,1.,-0.1652812010367921,-0.106412899,2.1652812010367914,1.,-0.106412899,1.,2.165281201036791,-0.106412899,0.4467935505,1.5532064495,1.,0.4467935505,0.4467935505,1.,1.5532064495,0.4467935505,1.,1.5532064495,1.5532064495,1.,1.0588683020367922,1.,2.106412899,1.,1.0588683020367922,2.106412899,0.9411316979632077,1.,2.106412899,1.,0.9411316979632078,2.106412899],27,3),1e-12))
        #
        weights=13*[1]
        gaussCoords=refCoords[:] ; gaussCoords[14]=0.9999999999999 # change z of point #4 0.999... instead of 1. because with shape function it leads to division by 0. !
        fGauss.setGaussLocalizationOnType(NORM_PYRA13,refCoords,gaussCoords,weights)
        arrOfDisc2=fGauss.getLocalizationOfDiscr()
        self.assertTrue(arrOfDisc2.isEqual(coo,1e-10)) # be less exigent 1e-10 instead of 1e-12 due to shape function sensitivity arount 0.,0.,1. !
        pass

    def testSwig2Tri7GP1(self):
        coo=DataArrayDouble([[0,0],[0,2],[2,0],[0,1],[1,1],[1,0],[0.6666666666666667,0.6666666666666667]])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo)
        m.allocateCells()
        # the cell description is exactly those described in the description of TRI7 in MED file 3.0.7 documentation
        m.insertNextCell(NORM_TRI7,range(7))
        refCoords=[0.,0.,1.,0.,0.,1.,0.5,0.,0.5,0.5,0.,0.5,0.3333333333333333,0.3333333333333333]
        gaussCoords=[0.3333333333333333,0.3333333333333333,0.470142064105115,0.470142064105115,0.05971587178977,0.470142064105115,0.470142064105115,0.05971587178977,0.101286507323456,0.101286507323456,0.797426985353088,0.101286507323456,0.101286507323456,0.797426985353088]
        weights=[0.062969590272413,0.062969590272413,0.062969590272413,0.066197076394253,0.066197076394253,0.066197076394253,0.1125]
        fGauss=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fGauss.setName("fGauss")
        fGauss.setMesh(m)
        fGauss.setGaussLocalizationOnType(NORM_TRI7,refCoords,gaussCoords,weights)
        arr=DataArrayDouble(fGauss.getNumberOfTuplesExpected()) ; arr.iota()
        fGauss.setArray(arr)
        arrOfDisc=fGauss.getLocalizationOfDiscr()
        self.assertTrue(arrOfDisc.isEqual(DataArrayDouble([0.666666666666667,0.666666666666667,0.9402841282102293,0.9402841282102293,0.9402841282102299,0.11943174357954002,0.11943174357953992,0.9402841282102299,0.20257301464691194,0.20257301464691196,0.20257301464691205,1.5948539707061757,1.5948539707061757,0.20257301464691202],7,2),1e-12))
        #
        weights=7*[1]
        gaussCoords=refCoords
        fGauss.setGaussLocalizationOnType(NORM_TRI7,refCoords,gaussCoords,weights)
        arrOfDisc2=fGauss.getLocalizationOfDiscr()
        self.assertTrue(arrOfDisc2.isEqual(coo,1e-12))
        pass

    def testSwig2StructuredDesc1(self):
        c=MEDCouplingCMesh()
        arr0=DataArrayDouble(3) ; arr0.iota()
        arr1=DataArrayDouble(4) ; arr1.iota()
        arr2=DataArrayDouble(5) ; arr2.iota()
        c.setCoords(arr0,arr1,arr2)
        #
        self.assertEqual(98,c.getNumberOfCellsOfSubLevelMesh())
        m=c.build1SGTSubLevelMesh()
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([0,12,15,3,12,24,27,15,24,36,39,27,36,48,51,39,3,15,18,6,15,27,30,18,27,39,42,30,39,51,54,42,6,18,21,9,18,30,33,21,30,42,45,33,42,54,57,45,1,13,16,4,13,25,28,16,25,37,40,28,37,49,52,40,4,16,19,7,16,28,31,19,28,40,43,31,40,52,55,43,7,19,22,10,19,31,34,22,31,43,46,34,43,55,58,46,2,14,17,5,14,26,29,17,26,38,41,29,38,50,53,41,5,17,20,8,17,29,32,20,29,41,44,32,41,53,56,44,8,20,23,11,20,32,35,23,32,44,47,35,44,56,59,47,0,12,13,1,12,24,25,13,24,36,37,25,36,48,49,37,1,13,14,2,13,25,26,14,25,37,38,26,37,49,50,38,3,15,16,4,15,27,28,16,27,39,40,28,39,51,52,40,4,16,17,5,16,28,29,17,28,40,41,29,40,52,53,41,6,18,19,7,18,30,31,19,30,42,43,31,42,54,55,43,7,19,20,8,19,31,32,20,31,43,44,32,43,55,56,44,9,21,22,10,21,33,34,22,33,45,46,34,45,57,58,46,10,22,23,11,22,34,35,23,34,46,47,35,46,58,59,47,0,1,4,3,3,4,7,6,6,7,10,9,1,2,5,4,4,5,8,7,7,8,11,10,12,13,16,15,15,16,19,18,18,19,22,21,13,14,17,16,16,17,20,19,19,20,23,22,24,25,28,27,27,28,31,30,30,31,34,33,25,26,29,28,28,29,32,31,31,32,35,34,36,37,40,39,39,40,43,42,42,43,46,45,37,38,41,40,40,41,44,43,43,44,47,46,48,49,52,51,51,52,55,54,54,55,58,57,49,50,53,52,52,53,56,55,55,56,59,58])))
        self.assertEqual(NORM_QUAD4,m.getCellModelEnum())
        #
        self.assertTrue(MEDCouplingStructuredMesh.Build1GTNodalConnectivityOfSubLevelMesh([3,7]).isEqual(DataArrayInt([0,3,3,6,6,9,9,12,12,15,15,18,1,4,4,7,7,10,10,13,13,16,16,19,2,5,5,8,8,11,11,14,14,17,17,20,0,1,1,2,3,4,4,5,6,7,7,8,9,10,10,11,12,13,13,14,15,16,16,17,18,19,19,20])))
        pass

    def testSwig2Colinearize2D1(self):
        coo=DataArrayDouble([-5.,0.,-1.,0.,4.,3.,7.,0.,1.,6.,1.,0.,-3.,0.,6.,1.,5.,0.,3.,0.],10,2)
        #
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_POLYGON,[5,9,8,3,7,2,4,0,6,1])
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([0])))
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,0,3,4])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4])))
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([])))
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,0,3,4])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4])))
        #
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_POLYGON,[8,3,7,2,4,0,6,1,5,9])
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([0])))
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,0,3,4])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4])))
        #
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_POLYGON,[3,7,2,4,0,6,1,5,9,8])
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([0])))
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,3,4,0])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4])))
        #
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo) ; m.allocateCells()
        m.insertNextCell(NORM_POLYGON,[4,0,6,1,5,9,8,3,7,2,])
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([0])))
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,4,0,3])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4])))
        ## false quadratic
        coo2=DataArrayDouble([(-5,0),(-1,0),(4,3),(7,0),(1,6),(1,0),(-3,0),(6,1),(5,0),(3,0),(2,0),(4,0),(6,0),(6.5,0.5),(5,2),(2.5,4.5),(-2,3),(-4,0),(-2,0),(0,0)])
        coo2.setInfoOnComponents(["aa","bbbb"])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo2) ; m.allocateCells()
        m.insertNextCell(NORM_QPOLYG,[5,9,8,3,7,2,4,0,6,1,10,11,12,13,14,15,16,17,18,19])
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([0])))
        self.assertNotEqual(refPtr,m.getCoords().getHiddenCppPointer())#not same coordinates here
        self.assertEqual(["aa","bbbb"],m.getCoords().getInfoOnComponents())
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(coo2.isEqual(m.getCoords()[:20],1e-12))
        self.assertTrue(m.getCoords()[20:].isEqualWithoutConsideringStr(DataArrayDouble([(1.,0.),(4.,3.)]),1e-12))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([32,0,3,4,20,21,16])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,7])))
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([])))
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([32,0,3,4,20,21,16])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,7])))
        # mix of quadratic and linear inside a QPOLYG cell
        coo2=DataArrayDouble([(-5,0),(-1,0),(7.,6.),(7,0),(1,6),(1,0),(-3,0),(8.2426406871192839,3),(5,0),(3,0),  (2,0),(4,0),(6,0),(7.9196888946291288,1.3764116995614091),(7.9196888946291288,4.6235883004385911),(4,7.2426406871192848),(-2,3),(-4,0),(-2,0),(0,0)])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo2) ; m.allocateCells()
        m.insertNextCell(NORM_QPOLYG,[5,9,8,3,7,2,4,0,6,1,10,11,12,13,14,15,16,17,18,19])
        refPtr=m.getCoords().getHiddenCppPointer()
        self.assertTrue(m.colinearize2D(1e-12).isEqual(DataArrayInt([0])))
        self.assertNotEqual(refPtr,m.getCoords().getHiddenCppPointer())#not same coordinates here
        self.assertTrue(coo2.isEqual(m.getCoords()[:20],1e-12))
        self.assertTrue(m.getCoords()[20:].isEqual(DataArrayDouble([(1.,0.),(7.,6.)]),1e-12))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([32,0,3,4,20,21,16])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,7])))
        pass

    def testSwig2BoundingBoxForBBTree1(self):
        """ This test appears simple but it checks that bounding box are correctly computed for quadratic polygons. It can help a lot to reduce the amount of intersections !
        """
        coo=DataArrayDouble([-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,-0.5,0.45,0.,0.3181980515339464,0.31819805153394637,0.,0.45,-0.31819805153394637,0.3181980515339464,-0.45,0.,-0.3181980515339465,-0.31819805153394637,0.,-0.45,0.3181980515339463,-0.3181980515339465,-0.5,0.0,0.0,0.5,0.5,0.0,0.0,-0.5,-0.4090990257669732,-0.4090990257669732,0.40909902576697316,-0.4090990257669732],18,2)
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_QPOLYG,[0,1,2,3,11,5,7,9,12,13,14,17,4,6,8,16])
        m.insertNextCell(NORM_QPOLYG,[3,0,9,11,15,16,10,17])
        self.assertTrue(m.getBoundingBoxForBBTree().isEqual(DataArrayDouble([-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,-0.31819805153394637],2,4),1e-12))
        pass

    def testSwig2CartBuildUnstructuredOnExoticCases1(self):
        """ Test focusing on traduction from cartesian to unstructured mesh when spaceDim greater than meshDim.
        """
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(4) ; arrY.iota()
        arrZ=DataArrayDouble(1) ; arrZ.iota()
        m.setCoords(arrX,arrY,arrZ)
        self.assertEqual(2,m.getMeshDimension())
        self.assertEqual(3,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([4,1,0,3,4,4,2,1,4,5,4,4,3,6,7,4,5,4,7,8,4,7,6,9,10,4,8,7,10,11])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25,30])))
        coo0=DataArrayDouble([(0,0,0),(1,0,0),(2,0,0),(0,1,0),(1,1,0),(2,1,0),(0,2,0),(1,2,0),(2,2,0),(0,3,0),(1,3,0),(2,3,0)])
        self.assertTrue(mu.getCoords().isEqual(coo0,1e-12))
        mu.writeVTK("tutu.vtu")
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(1) ; arrY.iota()
        arrZ=DataArrayDouble(4) ; arrZ.iota()
        m.setCoords(arrX,arrY,arrZ)
        self.assertEqual(2,m.getMeshDimension())
        self.assertEqual(3,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([4,1,0,3,4,4,2,1,4,5,4,4,3,6,7,4,5,4,7,8,4,7,6,9,10,4,8,7,10,11])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25,30])))
        coo1=DataArrayDouble([(0,0,0),(1,0,0),(2,0,0),(0,0,1),(1,0,1),(2,0,1),(0,0,2),(1,0,2),(2,0,2),(0,0,3),(1,0,3),(2,0,3)])
        self.assertTrue(mu.getCoords().isEqual(coo1,1e-12))
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(1) ; arrX.iota() ; arrX+=9
        arrY=DataArrayDouble(3) ; arrY.iota()
        arrZ=DataArrayDouble(4) ; arrZ.iota()
        m.setCoords(arrX,arrY,arrZ)
        self.assertEqual(2,m.getMeshDimension())
        self.assertEqual(3,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([4,1,0,3,4,4,2,1,4,5,4,4,3,6,7,4,5,4,7,8,4,7,6,9,10,4,8,7,10,11])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25,30])))
        coo2=DataArrayDouble([(9,0,0),(9,1,0),(9,2,0),(9,0,1),(9,1,1),(9,2,1),(9,0,2),(9,1,2),(9,2,2),(9,0,3),(9,1,3),(9,2,3)])
        self.assertTrue(mu.getCoords().isEqual(coo2,1e-12))
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(1) ; arrY.iota(7)
        arrZ=DataArrayDouble(1) ; arrZ.iota(8)
        m.setCoords(arrX,arrY,arrZ)
        self.assertEqual(1,m.getMeshDimension())
        self.assertEqual(3,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([1,0,1,1,1,2])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6])))
        coo3=DataArrayDouble([(0,7,8),(1,7,8),(2,7,8)])
        self.assertTrue(mu.getCoords().isEqual(coo3,1e-12))
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(1) ; arrX.iota(7)
        arrY=DataArrayDouble(1) ; arrY.iota(8)
        arrZ=DataArrayDouble(3) ; arrZ.iota()
        m.setCoords(arrX,arrY,arrZ)
        self.assertEqual(1,m.getMeshDimension())
        self.assertEqual(3,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([1,0,1,1,1,2])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6])))
        coo4=DataArrayDouble([(7,8,0),(7,8,1),(7,8,2)])
        self.assertTrue(mu.getCoords().isEqual(coo4,1e-12))
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(1) ; arrY.iota(7)
        m.setCoords(arrX,arrY)
        self.assertEqual(1,m.getMeshDimension())
        self.assertEqual(2,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([1,0,1,1,1,2])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6])))
        coo5=DataArrayDouble([(0,7),(1,7),(2,7)])
        self.assertTrue(mu.getCoords().isEqual(coo5,1e-12))
        #
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(1) ; arrX.iota(7)
        arrY=DataArrayDouble(3) ; arrY.iota()
        m.setCoords(arrX,arrY)
        self.assertEqual(1,m.getMeshDimension())
        self.assertEqual(2,m.getSpaceDimension())
        mu=m.buildUnstructured()
        self.assertTrue(mu.getNodalConnectivity().isEqual(DataArrayInt([1,0,1,1,1,2])))
        self.assertTrue(mu.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6])))
        coo6=DataArrayDouble([(7,0),(7,1),(7,2)])
        self.assertTrue(mu.getCoords().isEqual(coo6,1e-12))
        pass

    def testSwig2Colinearize2D2(self):
        """ simple non regression test but that has revealed a bug"""
        coo=DataArrayDouble([(0,0),(0,0.5),(0,1),(1,1),(1,0),(0.5,0)])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo)
        m.allocateCells() ; m.insertNextCell(NORM_POLYGON,[0,1,2,3,4,5])
        m.checkCoherency2()
        refPtr=m.getCoords().getHiddenCppPointer()
        #
        m.colinearize2D(1e-12)
        m.checkCoherency2()
        self.assertEqual(refPtr,m.getCoords().getHiddenCppPointer())
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,0,2,3,4])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5])))
        pass

    def testSwig2CheckAndPreparePermutation2(self):
        a=DataArrayInt([10003,9999999,5,67])
        self.assertTrue(DataArrayInt.CheckAndPreparePermutation(a).isEqual(DataArrayInt([2,3,0,1])))
        a=DataArrayInt([10003,-9999999,5,67])
        self.assertTrue(DataArrayInt.CheckAndPreparePermutation(a).isEqual(DataArrayInt([3,0,1,2])))
        a=DataArrayInt([])
        self.assertTrue(DataArrayInt.checkAndPreparePermutation(a).isEqual(DataArrayInt([])))
        pass

    def testSwig2ComputeNeighborsOfNodes1(self):
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(4) ; arrY.iota()
        arrZ=DataArrayDouble(5) ; arrZ.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY,arrZ) ; m=m.buildUnstructured()
        # 3D
        a,b=m.computeNeighborsOfNodes()
        self.assertTrue(a.isEqual(DataArrayInt([1,3,12,0,4,13,2,1,5,14,0,4,15,6,3,1,16,5,7,4,2,17,8,3,7,18,9,6,4,19,8,10,7,5,20,11,6,10,21,9,7,22,11,10,8,23,13,15,0,24,12,16,1,14,25,13,17,2,26,12,16,3,18,27,15,13,4,17,19,28,16,14,5,20,29,15,19,6,21,30,18,16,7,20,22,31,19,17,8,23,32,18,22,9,33,21,19,10,23,34,22,20,11,35,25,27,12,36,24,28,13,26,37,25,29,14,38,24,28,15,30,39,27,25,16,29,31,40,28,26,17,32,41,27,31,18,33,42,30,28,19,32,34,43,31,29,20,35,44,30,34,21,45,33,31,22,35,46,34,32,23,47,37,39,24,48,36,40,25,38,49,37,41,26,50,36,40,27,42,51,39,37,28,41,43,52,40,38,29,44,53,39,43,30,45,54,42,40,31,44,46,55,43,41,32,47,56,42,46,33,57,45,43,34,47,58,46,44,35,59,49,51,36,48,52,37,50,49,53,38,48,52,39,54,51,49,40,53,55,52,50,41,56,51,55,42,57,54,52,43,56,58,55,53,44,59,54,58,45,57,55,46,59,58,56,47])))
        self.assertTrue(b.isEqual(DataArrayInt([0,3,7,10,14,19,23,27,32,36,39,43,46,50,55,59,64,70,75,80,86,91,95,100,104,108,113,117,122,128,133,138,144,149,153,158,162,166,171,175,180,186,191,196,202,207,211,216,220,223,227,230,234,239,243,247,252,256,259,263,266])))
        # 2D
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY) ; m=m.buildUnstructured()
        a,b=m.computeNeighborsOfNodes()
        self.assertTrue(a.isEqual(DataArrayInt([1,3,0,4,2,1,5,0,4,6,3,1,5,7,4,2,8,3,7,9,6,4,8,10,7,5,11,6,10,9,7,11,10,8])))
        self.assertTrue(b.isEqual(DataArrayInt([0,2,5,7,10,14,17,20,24,27,29,32,34])))
        # 1D
        m=m.buildDescendingConnectivity()[0]
        a,b=m.computeNeighborsOfNodes()
        self.assertTrue(a.isEqual(DataArrayInt([1,3,0,4,2,1,5,0,4,6,3,1,5,7,4,2,8,3,7,9,6,4,8,10,7,5,11,6,10,9,7,11,10,8])))
        self.assertTrue(b.isEqual(DataArrayInt([0,2,5,7,10,14,17,20,24,27,29,32,34])))
        pass

    def setUp(self):
        pass
    pass

if __name__ == '__main__':
    unittest.main()
