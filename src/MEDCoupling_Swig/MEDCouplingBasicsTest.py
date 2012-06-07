#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

from MEDCoupling import *
import unittest
from math import pi,e,sqrt,cos,sin
from MEDCouplingDataForTest import MEDCouplingDataForTest

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
        self.assertTrue(mesh.getAllTypes()==[4])
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
        self.assertEqual(2,len(mesh.getAllTypes()));
        self.assertEqual(NORM_TRI3,mesh.getAllTypes()[0]);
        self.assertEqual(NORM_QUAD4,mesh.getAllTypes()[1]);
        self.assertEqual(1,len(subMesh.getAllTypes()));
        self.assertEqual(NORM_QUAD4,subMesh.getAllTypes()[0]);
        self.assertEqual(name,"PartOf_Toto");
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
        self.assertEqual(2,len(subMesh.getAllTypes()));
        self.assertEqual(NORM_TRI3,subMesh.getAllTypes()[0]);
        self.assertEqual(NORM_QUAD4,subMesh.getAllTypes()[1]);
        self.assertEqual(name,"PartOf_Toto");
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
        self.assertEqual(1,len(subMesh.getAllTypes()));
        self.assertEqual(NORM_QUAD4,subMesh.getAllTypes()[0]);
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
        self.assertEqual(2,len(subMesh.getAllTypes()));
        self.assertEqual(NORM_TRI3,subMesh.getAllTypes()[0]);
        self.assertEqual(NORM_QUAD4,subMesh.getAllTypes()[1]);
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
        self.assertEqual(2,len(subMesh.getAllTypes()));
        self.assertEqual(NORM_TRI3,subMesh.getAllTypes()[0]);
        self.assertEqual(NORM_QUAD4,subMesh.getAllTypes()[1]);
        self.assertEqual(3,subMesh.getNumberOfCells());
        pass
    
    def testZipCoords(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        self.assertEqual(2,len(mesh.getAllTypes()));
        self.assertEqual(2,mesh.getSpaceDimension());
        self.assertEqual(9,mesh.getNumberOfNodes());
        self.assertEqual(5,mesh.getNumberOfCells());
        oldConn=mesh.getNodalConnectivity().getValues()[0:mesh.getNodalConnectivity().getNbOfElems()];
        oldConnIndex=mesh.getNodalConnectivityIndex().getValues()[0:mesh.getNumberOfCells()+1]
        oldCoords=mesh.getCoords();
        mesh.zipCoords();
        self.assertEqual(2,len(mesh.getAllTypes()));
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
        self.assertEqual(NORM_QUAD4,subMesh.getAllTypes()[0]);
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
        self.assertEqual(NORM_QUAD4,subMesh.getAllTypes()[0]);
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
        self.assertEqual(f3.getTimeDiscretization(),NO_TIME);
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
        f1=m.fillFromAnalytic(ON_CELLS,1,"x+y");
        f1.checkCoherency();                    
        self.assertEqual(f1.getTypeOfField(),ON_CELLS);
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
        self.assertEqual(2,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        f1.applyFunc(1,"x+y");
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
        self.assertEqual(2,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        #
        f2=f1.clone(True);
        self.assertRaises(InterpKernelException, f2.applyFunc, 1, "a+b+c+d");
        self.assertRaises(InterpKernelException, f2.applyFunc, 1, "a/0");
        self.assertRaises(InterpKernelException, f2.applyFunc, "a/0");
        f2.applyFunc("abs(u)^2.4+2*u");
        self.assertEqual(f1.getTypeOfField(),ON_NODES);
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f3.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f3.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f4.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f4.getTimeDiscretization(),NO_TIME);
        tmp=f4.getArray().getValues();
        for i in xrange(len(tmp)):
            self.assertTrue(abs(tmp[i]-2.)<1.e-12)
            pass
        #
        f4=f2.buildNewTimeReprFromThis(ONE_TIME,False);
        f4.checkCoherency();
        self.assertEqual(f4.getTypeOfField(),ON_NODES);
        self.assertEqual(f4.getTimeDiscretization(),ONE_TIME);
        self.assertRaises(InterpKernelException,f1.__add__,f4);
        f5=f4.buildNewTimeReprFromThis(NO_TIME,False);
        self.assertEqual(f5.getTypeOfField(),ON_NODES);
        self.assertEqual(f5.getTimeDiscretization(),NO_TIME);
        f3=f1+f5;
        tmp=f3.getArray().getValues();
        values4=[-1.2,-0.2,0.8,-0.2,0.8,1.8,0.8,1.8,2.8]
        self.assertEqual(len(values3),len(tmp))
        for i in xrange(len(tmp)):
            self.assertTrue(abs(tmp[i]-values4[i])<1.e-12)
            pass
        #
        f4=f2.buildNewTimeReprFromThis(ONE_TIME,True);
        f4.checkCoherency();
        self.assertEqual(f4.getTypeOfField(),ON_NODES);
        self.assertEqual(f4.getTimeDiscretization(),ONE_TIME);
        self.assertRaises(InterpKernelException,f1.__add__,f4);
        f5=f4.buildNewTimeReprFromThis(NO_TIME,True);
        self.assertEqual(f5.getTypeOfField(),ON_NODES);
        self.assertEqual(f5.getTimeDiscretization(),NO_TIME);
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
        f1=m.fillFromAnalytic(ON_NODES,1,"x+y+z");
        f2=m.fillFromAnalytic(ON_NODES,1,"a*a+b+c*c");
        f3=f1/f2;
        f3.checkCoherency();
        self.assertEqual(f3.getTypeOfField(),ON_NODES);
        self.assertEqual(f3.getTimeDiscretization(),NO_TIME);
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
        self.assertEqual(f1.getTimeDiscretization(),NO_TIME);
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
        targetMesh.rotate(center,[],0.78539816339744830962);
        t1=None
        t2=None
        t1,t2=targetMesh.getCellsContainingPoints(pos,6,1e-12);
        self.assertEqual(6,t1.getNumberOfTuples());
        self.assertEqual(7,t2.getNumberOfTuples());
        self.assertEqual(list(t1.getValues()),expectedValues1);
        self.assertEqual(list(t2.getValues()),expectedValues2);
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
        self.assertEqual(3, len(dis));
        self.assertEqual(NORM_HEXA8, dis[0]);
        self.assertEqual(27, dis[1]);
        self.assertEqual(0, dis[2]);
        
        idsPerType=[]
        self.assertRaises(InterpKernelException, mesh1.checkTypeConsistencyAndContig, dis, idsPerType);
        dis[2]=-1;
        idsPerType=[]
        self.assertTrue(not mesh1.checkTypeConsistencyAndContig(dis, idsPerType));
        dis[0]=NORM_QUAD4;
        self.assertRaises(InterpKernelException, mesh1.checkTypeConsistencyAndContig, dis, idsPerType);
        
        dis[0]=NORM_HEXA8;
        dis[2]=0;
        ids=DataArrayInt.New();
        ids.alloc(10, 1);
        ids.fillWithValue(111);
        idsPerType=[ids];
        check=mesh1.checkTypeConsistencyAndContig(dis, idsPerType);
        self.assertTrue(check);
        self.assertTrue(check.isEqual(ids));
        
        code, idsInPflPerType, pfls=mesh1.splitProfilePerType(ids);
        self.assertEqual(3, len(code));
        self.assertEqual(NORM_HEXA8, code[0]);
        self.assertEqual(27, code[1]);
        self.assertEqual(0, code[2]);
        self.assertEqual(1, len(idsInPflPerType));
        self.assertEqual(1, len(pfls));
        self.assertTrue(idsInPflPerType[0].isEqual(ids));
        self.assertTrue(pfls[0].isEqual(ids));
        
        cells1=[0, 1, 25, 26]
        partMesh1=mesh1.buildPart(cells1)
        self.assertTrue(isinstance(partMesh1,MEDCouplingMesh))
        self.assertEqual(4, partMesh1.getNumberOfCellsWithType(NORM_HEXA8));
        self.assertEqual(64, mesh1.getNumberOfNodes());
        self.assertEqual(64, partMesh1.getNumberOfNodes());
        
        cells2=[25, 26]
        partMesh2, arr1=mesh1.buildPartAndReduceNodes(cells2)
        self.assertTrue(isinstance(partMesh2,MEDCouplingUMesh))
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
        types=mesh.getAllTypes();
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
        types2=mesh.getAllTypes();
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
        self.assertAlmostEqual(5.125,f.getWeightedAverageValue(),14);
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
        f2=f1.buildSubPart(part1);
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
        f1=m1.getEdgeRatioField();
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
        self.assertRaises(InterpKernelException, dbl2.reverse);
        self.assertRaises(InterpKernelException, dbl2.iota, 10.);
        
        dbl=DataArrayDouble.New();
        #DataArrayDouble not allocated yet
        self.assertRaises(InterpKernelException, dbl.iota, 10.);
        self.assertRaises(InterpKernelException, dbl.isUniform, 10., 1e-15);
        self.assertRaises(InterpKernelException, dbl.sort);
        self.assertRaises(InterpKernelException, dbl.reverse);
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
        c=mesh.getNodeIdsNearPoint(pts,1e-7);
        self.assertEqual([4,9,11],c.getValues());
        c,cI=mesh.getNodeIdsNearPoints(pts,3,1e-7);
        self.assertEqual([0,3,3,4],cI.getValues());
        self.assertEqual([4,9,11,6],c.getValues());
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
        d.alloc(5,1);
        d.iota();
        e=MEDCouplingCMesh.New();
        e.setCoordsAt(0,d);
        f=e.buildUnstructured();
        d2=f.getCoords().applyFunc("x*x/2");
        f.setCoords(d2);
        f.changeSpaceDimension(2);
        #
        center=[0.,0.]
        f.rotate(center,[],pi/3);
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
        types=m.getAllTypes();
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
        types=m.getAllTypes();
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
        s=da1.getDifferentValues(True);# API different from C++ because SWIG complains...
        expected1=[1,2,3,5]
        self.assertEqual(expected1,s);
        pass

    def testSwigErrorProtection3(self):
        da=DataArrayInt.New()
        da.setValues([1,2,3,4],4,3)
        self.assertEqual([1, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        da=DataArrayInt.New()
        da.setValues((1,2,3,4,4,3),4,3)
        self.assertEqual([1, 2, 3, 4, 4, 3, 0, 0, 0, 0, 0, 0],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        da.setValues(10*[1]+290*[2],4,3)
        self.assertEqual(10*[1]+[2,2],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        #
        da=DataArrayDouble.New()
        da.setValues([1,2,3.,4],4,3)
        self.assertEqual([1., 2., 3., 4., 0., 0., 0., 0., 0., 0., 0., 0.],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        da=DataArrayDouble.New()
        da.setValues((1,2,3,4.,4,3),4,3)
        self.assertEqual([1., 2., 3., 4., 4., 3., 0., 0., 0., 0., 0., 0.],da.getValues())
        self.assertEqual(3,da.getNumberOfComponents());
        self.assertEqual(4,da.getNumberOfTuples());
        da.setValues(10*[1]+290*[2],4,3)
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
        try:
            da2=da[5:8,-2]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
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
        try:
            da2=da[5:8,-2]
        except InterpKernelException as e:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
            pass
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
        types1=m.getAllTypes();
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

    def testUMeshFindCellsIdsOnBoundary1(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        da5=m.findCellsIdsOnBoundary();
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
        self.assertRaises(InterpKernelException,m1.fillFromAnalytic2,ON_NODES,1,"y+z");
        m1.getCoords().setInfoOnComponent(0,"x [m]");
        m1.getCoords().setInfoOnComponent(1,"y");
        m1.getCoords().setInfoOnComponent(2,"z");
        f1=m1.fillFromAnalytic2(ON_NODES,1,"y+z");
        self.assertEqual(1,f1.getNumberOfComponents());
        self.assertEqual(9,f1.getNumberOfTuples());
        expected1=[0.2, 0.7, 1.2, 0.7, 1.2, 1.7, 1.2, 1.7, 2.2]
        for i in xrange(9):
            self.assertAlmostEqual(expected1[i],f1.getArray().getIJ(0,i),12);
            pass
        pass

    def testFillFromAnalyticThree1(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        vs=3*[None];
        vs[0]="x"; vs[1]="Y"; vs[2]="z";
        self.assertRaises(InterpKernelException,m1.fillFromAnalytic3,ON_NODES,1,vs,"y+z");
        vs[1]="y";
        f1=m1.fillFromAnalytic3(ON_NODES,1,vs,"y+z");
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
        self.assertEqual(6,len(code));
        self.assertEqual(2,len(idsInPflPerType));
        expected1=[3,1,0, 4,2,1]
        for i in xrange(6):
            self.assertEqual(expected1[i],code[i]);
            pass
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
        self.assertEqual(6,len(code));
        self.assertEqual(2,len(idsInPflPerType));
        expected2=[3,1,0, 4,3,-1]
        for i in xrange(6):
            self.assertEqual(expected2[i],code[i]);
            pass
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
        self.assertEqual(6,len(code));
        self.assertEqual(2,len(idsInPflPerType));
        expected3=[3,2,0, 4,1,1]
        for i in xrange(6):
            self.assertEqual(expected3[i],code[i]);
            pass
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
        self.assertEqual(3,len(code));
        self.assertEqual(1,len(idsInPflPerType));
        expected4=[4,2,0]
        for i in xrange(3):
            self.assertEqual(expected4[i],code[i]);
            pass
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
        self.assertEqual(6,len(code));
        self.assertEqual(3,code[0]);
        self.assertEqual(2,code[1]);
        self.assertEqual(0,code[2]);
        self.assertEqual(4,code[3]);
        self.assertEqual(3,code[4]);
        self.assertEqual(0,code[5]);
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
        expected1=[14, 1, 2, 3, 4,
                   18, 5, 6, 7, 8, 9, 10, 11, 12,
                   14, 13, 14, 15, 16,
                   31, 17, 18, 19, -1, 20, 22, 21, -1, 17, 18, 21, 20, -1, 18, 19, 22, 21, -1, 19, 17, 20, 22,
                   16, 23, 24, 25, 26, 27, 28,
                   31, 29, 30, 31, 32, 33, -1, 34, 38, 37, 36, 35, -1, 29, 30, 35, 34, -1, 30, 31, 36, 35, -1, 31, 32, 37, 36, -1, 32, 33, 38, 37, -1, 33, 29, 34, 38,
                   18, 39, 40, 41, 42, 43, 44, 45, 46,
                   22, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
                   31, 59, 60, 61, 62, 63, 64, 65, -1, 66, 72, 71, 70, 69, 68, 67, -1, 59, 60, 67, 66, -1, 60, 61, 68, 67, -1, 61, 62, 69, 68, -1, 62, 63, 70, 69, -1, 63, 64, 71, 70, -1, 64, 65, 72, 71, -1, 65, 59, 66, 72];
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
        f2=f1.buildNewTimeReprFromThis(ONE_TIME,False)
        f2.setArray(bary)
        self.assertRaises(InterpKernelException,f2.copyTinyAttrFrom,f1)
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
        
        array12=[0.]*(6*4)
        da.setValues(array12,6,4) #bad NumberOfComponents
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
        b.setValues(arr2,5,1);
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
        expected1=[0,1,1,2,2]
        expected2=[0,0,1,1,2]
        self.assertEqual(5,d1.getNumberOfTuples());
        self.assertEqual(5,d2.getNumberOfTuples());
        self.assertEqual(5,m3.getNumberOfCells());
        self.assertEqual(22,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[5,17,1,16,12,5,18,1,17,13,5,19,2,18,13,5,20,2,19,14,5,21,3,20,14]
        expected4=[0,5,10,15,20,25]
        expected5=[-1.0,2.0,1.0,2.0,2.0,2.0,4.0,2.0,-1.0,4.0,1.0,4.0,2.0,4.0,4.0,4.0,-0.5,-1.5,1.5,-1.5,2.5,-1.5,4.5,-1.5,-0.5,2.5,1.5,2.5,2.5,2.5,4.5,2.5,-0.5,2.0,1.0,2.5,1.5,2.0,2.0,2.5,2.5,2.0,4.0,2.5]
        self.assertEqual(25,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(6,m3.getNodalConnectivityIndex().getNumberOfTuples());
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
        expected1=[0,1,1,2,3,3,4,5,5,6,7,7]
        expected2=[0,0,1,2,2,3,4,4,5,6,6,7]
        self.assertEqual(12,d1.getNumberOfTuples());
        self.assertEqual(12,d2.getNumberOfTuples());
        self.assertEqual(12,m3.getNumberOfCells());
        self.assertEqual(88,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[6,28,1,25,44,45,46,8,26,1,28,27,47,48,49,50,8,40,2,26,27,51,52,53,54,6,28,25,5,55,56,57,8,28,5,32,31,58,59,60,61,8,32,6,41,31,62,63,64,65,6,25,37,5,66,67,68,8,32,5,37,36,69,70,71,72,8,42,6,32,36,73,74,75,76,6,1,37,25,77,78,79,8,37,1,26,38,80,81,82,83,8,26,2,43,38,84,85,86,87]
        expected4=[0,7,16,25,32,41,50,57,66,75,82,91,100]
        expected5=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1.,-1.1,-1.,0.,-1.,1.1,-1.,1.7,-1.,1.118033988749895,1.,-1.118033988749895,1.,-1.118033988749895,-1.,1.118033988749895,-1.,0.7071067811865477,0.7071067811865476,0.5,0.,0.,0.5,1.05,0.,0.7071067811865475,0.7071067811865477,0.55,1.,1.1,0.5,1.4012585384440737,0.535233134659635,1.3,0.,1.1,0.5,1.1090169943749475,1.,0.,0.5,-0.5,0.,-0.7071067811865477,0.7071067811865476,-0.7071067811865475,0.7071067811865477,-1.05,0.,-1.1,0.5,-0.55,1.,-1.3,0.,-1.4012585384440737,0.5352331346596344,-1.1090169943749475,1.,-1.1,0.5,0.,-0.5,-0.7071067811865475,-0.7071067811865477,-0.5,0.,-1.05,0.,-0.7071067811865478,-0.7071067811865475,-0.55,-1.,-1.1,-0.5,-1.4012585384440732,-0.5352331346596354,-1.3,0.,-1.1,-0.5,-1.1090169943749475,-1.,0.7071067811865475,-0.7071067811865477,0.,-0.5,0.5,0.,0.7071067811865477,-0.7071067811865475,1.05,0.,1.1,-0.5,0.55,-1.,1.3,0.,1.4012585384440737,-0.535233134659635,1.1090169943749475,-1.,1.1,-0.5]
        self.assertEqual(100,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(13,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in xrange(176):
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
        expected1=[0,0,1,2,2,3,4,4,5,6,6,7]
        expected2=[0,1,1,2,3,3,4,5,5,6,7,7]
        self.assertEqual(12,d1.getNumberOfTuples());
        self.assertEqual(12,d2.getNumberOfTuples());
        self.assertEqual(12,m3.getNumberOfCells());
        self.assertEqual(88,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[6,16,15,18,44,45,46,8,18,2,1,16,47,48,49,50,8,17,1,2,40,51,52,53,54,6,18,15,20,55,56,57,8,20,7,6,18,58,59,60,61,8,41,6,7,21,62,63,64,65,6,20,15,22,66,67,68,8,22,11,7,20,69,70,71,72,8,21,7,11,42,73,74,75,76,6,22,15,16,77,78,79,8,16,1,13,22,80,81,82,83,8,43,13,1,17,84,85,86,87]
        expected4=[0,7,16,25,32,41,50,57,66,75,82,91,100]
        expected5=[0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1.,-1.1,-1.,0.,-1.,1.1,-1.,1.7,-1.,0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,1.1180339887498951,1.,-1.1180339887498951,1.,-1.1180339887498951,-1.,1.1180339887498951,-1.,0.5,0.,0.,0.5,0.7071067811865477,0.7071067811865476,0.55,1.,1.1,0.5,1.05,0.,0.7071067811865477,0.7071067811865477,1.3,0.,1.1,0.5,1.1090169943749475,1.,1.4012585384440737,0.535233134659635,0.,0.5,-0.5,0.,-0.7071067811865477,0.7071067811865476,-1.05,0.,-1.1,0.5,-0.55,1.,-0.7071067811865477,0.7071067811865477,-1.1090169943749475,1.,-1.1,0.5,-1.3,0.,-1.4012585384440737,0.5352331346596344,-0.5,0.,0.,-0.5,-0.7071067811865475,-0.7071067811865477,-0.55,-1.,-1.1,-0.5,-1.05,0.,-0.7071067811865479,-0.7071067811865476,-1.3,0.,-1.1,-0.5,-1.1090169943749475,-1.,-1.4012585384440734,-0.5352331346596354,0.,-0.5,0.5,0.,0.7071067811865475,-0.7071067811865477,1.05,0.,1.1,-0.5,0.55,-1.,0.7071067811865477,-0.7071067811865476,1.1090169943749475,-1.,1.1,-0.5,1.3,0.,1.4012585384440737,-0.535233134659635]
        self.assertEqual(100,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(13,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in xrange(176):
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
        expected2=[[0,1,6,5],[1,2,7,6],[2,3,8,7],[3,4,9,8],[4,5,11,10],[5,6,12,11],[6,7,13,12],[7,8,14,13],[8,9,16,15],[9,10,17,16],[10,11,18,17],[11,12,19,18]]
        for i in xrange(m1c.getNumberOfCells()):
            self.assertEqual(expected2[i],m1c.getNodeIdsOfCell(i))
            pass
        # test in 3D
        m1c.setCoordsAt(2,coordsZ);
        expected3=[[0,1,6,5,20,21,26,25],[1,2,7,6,21,22,27,26],[2,3,8,7,22,23,28,27],[3,4,9,8,23,24,29,28],[4,5,11,10,24,25,31,30],[5,6,12,11,25,26,32,31],[6,7,13,12,26,27,33,32],[7,8,14,13,27,28,34,33],[8,9,16,15,28,29,36,35],[9,10,17,16,29,30,37,36],[10,11,18,17,30,31,38,37],[11,12,19,18,31,32,39,38],[20,21,26,25,40,41,46,45],[21,22,27,26,41,42,47,46],[22,23,28,27,42,43,48,47],[23,24,29,28,43,44,49,48],[24,25,31,30,44,45,51,50],[25,26,32,31,45,46,52,51],[26,27,33,32,46,47,53,52],[27,28,34,33,47,48,54,53],[28,29,36,35,48,49,56,55],[29,30,37,36,49,50,57,56],[30,31,38,37,50,51,58,57],[31,32,39,38,51,52,59,58]]
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
            self.assertAlmostEqual(expected9[i],mesh3D_4.getCoords().getIJ(0,i),1e-12);
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
        self.assertRaises(InterpKernelException,da.__getitem__,-1)
        self.assertRaises(InterpKernelException,da.__getitem__,3)
        da=DataArrayInt.New([4,5,6,7,8,9],2,3)
        self.assertEqual(9,da[1,2])
        da=DataArrayDouble.New([4.1,5.2,6.3])
        self.assertAlmostEqual(5.2,da[1],12)
        self.assertRaises(InterpKernelException,da.__getitem__,-1)
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

    def setUp(self):
        pass
    pass

unittest.main()
