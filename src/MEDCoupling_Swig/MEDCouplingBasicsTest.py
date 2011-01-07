#  -*- coding: iso-8859-1 -*-
#  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

from MEDCoupling import *
import unittest
from math import pi,e,sqrt
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
        self.assertRaises(Exception,meshM1D.getMeshDimension);
        self.assertRaises(Exception,meshM1D.getNumberOfNodes);
        self.assertRaises(Exception,meshM1D.getNumberOfCells);
        self.assertRaises(Exception,meshM1D.setMeshDimension,-2)
        self.assertRaises(Exception,meshM1D.setMeshDimension,-10)
        meshM1D.setMeshDimension(-1);
        meshM1D.checkCoherency();
        self.assertEqual(meshM1D.getMeshDimension(),-1);
        self.assertEqual(meshM1D.getNumberOfCells(),1);
        self.assertRaises(Exception,meshM1D.getNumberOfNodes);
        self.assertRaises(Exception,meshM1D.getSpaceDimension);
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
        revNodal=DataArrayInt.New();
        revNodalIndx=DataArrayInt.New();
        mesh.getReverseNodalConnectivity(revNodal,revNodalIndx);
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
        self.assertEqual(expected1,list(boundaryNodes));
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
        subMesh=subMesh.buildPartOfMySelf(range(3),True);
        self.assertEqual("PartOf_Toto",subMesh.getName());
        pass
    def testBuildPartOfMySelfNode(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        tab1=[5,7]
        subMesh=mesh.buildPartOfMySelfNode(tab1[0:2],True);
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
        subMesh=mesh.buildPartOfMySelfNode(tab1[0:2],False);
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
        self.assertRaises(Exception,field.setNature,Integral);
        self.assertRaises(Exception,field.setNature,IntegralGlobConstraint);
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
        self.assertEqual(6,ret1.getNumberOfNodes());
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
        comm,commI=targetMesh.findCommonNodes(-1,1e-10);
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
        comm,commI=targetMesh.findCommonNodes(-1,1e-10);
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
        self.assertEqual(3,cells[0]);
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
        self.assertEqual(3,cells[0]);
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
        self.assertRaises(Exception,m.fillFromAnalytic,ON_NODES,1,"1./(x-0.2)");
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
        self.assertRaises(Exception,f1.__add__,f4);
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
        self.assertRaises(Exception,f1.__add__,f4);
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
        # self.assertRaises(Exception,f2.__imul__,f1);
        pass

    def testOperationsOnFields4(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        nbOfCells=m.getNumberOfCells();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,CONST_ON_TIME_INTERVAL);
        f1.setMesh(m);
        array=DataArrayDouble.New();
        f1.setArray(array);
        self.assertRaises(Exception,f1.setEndArray,array);
        self.assertRaises(Exception,f1.getEndArray);
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
        self.assertRaises(Exception,f1.getValueOn,pos,3.2)
        f2=MEDCouplingFieldDouble.New(ON_CELLS,LINEAR_TIME);
        f2.setMesh(m);
        f2.setArray(f1.getArray());
        f2.setStartTime(2.,3,0);
        f2.setEndTime(4.,13,0);
        self.assertRaises(Exception,f2.checkCoherency)
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
        self.assertRaises(Exception,f1.mergeNodes,1.e-10)
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
        self.assertTrue(expected1,arr1.getValues());
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
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
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
        self.assertEqual(6,len(t1));
        self.assertEqual(7,len(t2));
        expectedValues1=[0,4,3,0,1,2]
        expectedValues2=[0,1,2,3,4,5,6]
        self.assertEqual(list(t1),expectedValues1);
        self.assertEqual(list(t2),expectedValues2);
        #2D with no help of bounding box.
        center=[0.2,0.2]
        MEDCouplingPointSet.Rotate2DAlg(center,0.78539816339744830962,6,pos);
        targetMesh.rotate(center,[],0.78539816339744830962);
        t1=None
        t2=None
        t1,t2=targetMesh.getCellsContainingPoints(pos,6,1e-12);
        self.assertEqual(6,len(t1));
        self.assertEqual(7,len(t2));
        self.assertEqual(list(t1),expectedValues1);
        self.assertEqual(list(t2),expectedValues2);
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
        self.assertEqual(list(t1),expectedValues3);
        pos3=[0.2,0.2]
        t1=None
        t1=targetMesh.getCellsContainingPoint(pos3,1e-12);
        self.assertEqual(5,len(t1));
        expectedValues4=[0,1,2,3,4]
        self.assertEqual(list(t1),expectedValues4);
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
        self.assertEqual(list(t1),expectedValues5);
        pos6=[0., 50., 0.]
        t1=None
        t1=targetMesh.getCellsContainingPoint(pos6,1e-12);
        self.assertEqual(2,len(t1));
        expectedValues6=[0,2]
        self.assertEqual(list(t1),expectedValues6);
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
        #self.assertEqual(m1.getCoords()!=m2.getCoords());
        m1.tryToShareSameCoords(m2,1e-12);
        #self.assertEqual(m1.getCoords()==m2.getCoords());
        m1.tryToShareSameCoords(m2,1e-12);
        #self.assertEqual(m1.getCoords()==m2.getCoords());
        m2.tryToShareSameCoords(m1,1e-12);
        #self.assertEqual(m1.getCoords()==m2.getCoords());
        #
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_2();
        #self.assertEqual(m1.getCoords()!=m2.getCoords());
        m1.tryToShareSameCoords(m2,1e-12);
        #self.assertEqual(m1.getCoords()==m2.getCoords());
        m1.tryToShareSameCoords(m2,1e-12);
        #self.assertEqual(m1.getCoords()==m2.getCoords());
        m2.tryToShareSameCoords(m1,1e-12);
        #self.assertEqual(m1.getCoords()==m2.getCoords());
        #
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DSourceMesh_1();
        #self.assertEqual(m1.getCoords()!=m2.getCoords());
        self.assertRaises(Exception,m1.tryToShareSameCoords,m2,1e-12)
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
        self.assertEqual(0,f.getNbOfGaussLocalization());
        f.setGaussLocalizationOnType(NORM_TRI3,_refCoo1,_gsCoo1,_wg1);
        self.assertRaises(Exception,f.setGaussLocalizationOnType,NORM_QUAD4,_refCoo1,_gsCoo1,_wg1)
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
        self.assertRaises(Exception,f.checkCoherency);
        ids1=[0,1,3,4]
        self.assertRaises(Exception,f.setGaussLocalizationOnCells,ids1,_refCoo2,_gsCoo1,_wg1);
        self.assertEqual(0,f.getNbOfGaussLocalization());
        ids2=[0,4]
        f.setGaussLocalizationOnCells(ids2,_refCoo2,_gsCoo1,_wg1);
        self.assertEqual(1,f.getNbOfGaussLocalization());
        self.assertEqual(0,f.getGaussLocalizationIdOfOneCell(0));
        self.assertRaises(Exception,f.getGaussLocalizationIdOfOneCell,1);
        ids3=[1,2]
        f.setGaussLocalizationOnCells(ids3,_refCoo1,_gsCoo1,_wg1);
        self.assertEqual(2,f.getNbOfGaussLocalization());
        self.assertEqual(0,f.getGaussLocalizationIdOfOneCell(0));
        self.assertEqual(1,f.getGaussLocalizationIdOfOneCell(1));
        self.assertEqual(1,f.getGaussLocalizationIdOfOneCell(2));
        self.assertRaises(Exception,f.checkCoherency);#<- cell 3 has no localization
        ids4=[3]
        _gsCoo2=_gsCoo1;
        _wg2=_wg1;
        _gsCoo2[0]=0.8888777776666;
        _wg2[0]=0.1234567892377;
        f.setGaussLocalizationOnCells(ids4,_refCoo2,_gsCoo2,_wg2);
        self.assertEqual(3,f.getNbOfGaussLocalization());
        tmpIds=f.getCellIdsHavingGaussLocalization(0);
        self.assertEqual(ids2,list(tmpIds));
        self.assertRaises(Exception,f.checkCoherency);#<- it's always not ok because undelying array not with the good size.
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
        vec=[0.,0.,1.]
        self.assertRaises(Exception,m.are2DCellsNotCorrectlyOriented,vec,False);
        m.changeSpaceDimension(3);
        res1=m.are2DCellsNotCorrectlyOriented(vec,False);
        self.assertTrue(len(res1)==0);
        vec[2]=-1.;
        res1=m.are2DCellsNotCorrectlyOriented(vec,False);
        self.assertEqual(5,len(res1));
        #
        vec[2]=1.;
        # connectivity inversion
        conn=m.getNodalConnectivity().getValues();
        tmp=conn[11];
        conn[11]=conn[12];
        conn[12]=tmp;
        m.getNodalConnectivity().setValues(conn,len(conn),1)
        res1=m.are2DCellsNotCorrectlyOriented(vec,False);
        self.assertEqual(1,len(res1));
        self.assertEqual(2,res1[0]);
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
        vec=[0.,0.,-1.]
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
        expected5=[11.3068,27.3621,43.7881]
        for i in xrange(3):
            self.assertTrue(abs(expected5[i]-res[i])<1e-12);
            pass
        self.assertTrue(abs(expected5[0]-f1.normL1(0))<1e-12);
        self.assertTrue(abs(expected5[1]-f1.normL1(1))<1e-12);
        self.assertTrue(abs(expected5[2]-f1.normL1(2))<1e-12);
        #normL2
        res=f1.normL2();
        self.assertTrue(3,len(res))
        expected7=[9.0252562290496776, 21.545259176904789, 34.433193070059595]
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
            self.assertTrue(abs(sqrt(2.)*expected5[i]-res[i])<1e-12);
            pass
        res=f1.normL2();
        for i in xrange(3):
            self.assertTrue(abs(sqrt(sqrt(2.))*expected7[i]-res[i])<1e-12);
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
        self.assertTrue(expected1==types);
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
        self.assertTrue(expected2==types2);
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
        self.assertRaises(Exception,mesh1.checkGeoEquivalWith,mesh2,0,1e-12);#deepEqual fails
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
        self.assertRaises(Exception,mesh1.checkGeoEquivalWith,mesh2,0,1e-12);#deepEqual fails
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
        self.assertRaises(Exception,mesh1.checkGeoEquivalWith,mesh2,0,1e-12)
        self.assertTrue(cellCor==None);
        self.assertTrue(nodeCor==None);
        self.assertRaises(Exception,mesh1.checkGeoEquivalWith,mesh2,1,1e-12)
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
        m.renumberNodes(renum1,9);
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
        m2.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
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
        self.assertRaises(Exception,m1.tryToShareSameCoordsPermute,m2,1e-12);# <- here in this order the sharing is impossible.
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
        expected2=[7.,107.,9.,109.,17.,117.,10.,110.,11.,111.,12.,112.,13.,113.,15.,115.,14.,114.,16.,116.,8.,108.]
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
        renum=[0,2,1,3,4,5,6,8,7,9]
        mesh2.renumberCells(renum,False);
        #
        f2=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        f2.setMesh(mesh2);
        array=DataArrayDouble.New();
        arr2=[7.1,107.1,9.1,109.1,8.1,108.1,10.1,110.1,11.1,111.1,12.1,112.1,13.1,113.1,15.1,115.1,14.1,114.1,16.1,116.1]
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
        part2=[1,4,2,5]
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
        part3=[1,4,2,5,7]
        arrr=DataArrayInt.New();
        arrr.setValues(part3,5,1);
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
        part4=[1,4,2,5,7,8]
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
        self.assertRaises(Exception,f1.checkCoherency);#no end array specified !
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
        self.assertRaises(Exception,mesh1.checkGeoEquivalWith,mesh2,0,1e-12);
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
        self.assertRaises(Exception,f1.fillFromAnalytic,1,"y+x");
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
        self.assertRaises(Exception,f1.fillFromAnalytic,1,"1./(x-0.2)");
        pass

    def testFieldDoubleOpEqual1(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        self.assertRaises(Exception,f1.assign,0.07);
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
        mesh.insertNextCell(NORM_HEXA8,8,tmpConn[0:8])
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
        self.assertEqual(0,m.getCellContainingPoint([2.4],12));
        self.assertEqual(1,m.getCellContainingPoint([3.7],12));
        self.assertEqual(2,m.getCellContainingPoint([5.9],12));
        self.assertEqual(-1,m.getCellContainingPoint([10.3],12));
        self.assertEqual(-1,m.getCellContainingPoint([1.3],12));
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
        elems=m2.giveElemsInBoundingBox([3.5,6.,12.2,25.,0.,1.5],1e-7)
        self.assertEqual([1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17],elems)
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
        self.assertRaises(Exception,f1.checkCoherency);#no end array specified !
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
        self.assertRaises(Exception,a2.keepSelectedComponents,arr5V);
        self.assertRaises(Exception,a2.keepSelectedComponents,arr6V);
        self.assertRaises(Exception,a2.setSelectedComponents,a1,arr7V);
        arr7V=arr7V[0:3]
        self.assertRaises(Exception,a2.setSelectedComponents,a1,arr7V);
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
        self.assertEqual([4,9,11],c);
        c,cI=mesh.getNodeIdsNearPoints(pts,3,1e-7);
        self.assertEqual([0,3,3,4],cI);
        self.assertEqual([4,9,11,6],c);
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
        vec1=[10.,0.,0.]
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
        self.assertRaises(Exception,a.selectByTupleIdSafe,arr4);
        arr5=[4,2,0,6,7]
        self.assertRaises(Exception,a.selectByTupleIdSafe,arr5);
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
        self.assertRaises(Exception,c.selectByTupleIdSafe,arr4);
        self.assertRaises(Exception,c.selectByTupleIdSafe,arr5);
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
        self.assertRaises(Exception,m.rotate,[0.,0.,0.],(0.3,6,"1.2"),0.37)
        self.assertRaises(Exception,m.rotate,[0.,"0.",0.],[0.3,0.6,1.2],0.37)
        self.assertRaises(Exception,m.rotate,[0.,0.,0.],[0.3,'0.6',1.2],0.37)
        m2=m.buildPartOfMySelf([2,5],True)
        m3=m.buildPartOfMySelf((2,5),True)
        self.assertTrue(m2.isEqual(m3,1e-12))
        self.assertRaises(Exception,m.buildPartOfMySelf,[2,5.],True)
        da1=m.getCoords().keepSelectedComponents([1])
        da2=m.getCoords().keepSelectedComponents((1,))
        self.assertTrue(da1.isEqual(da2,1e-12))
        self.assertRaises(Exception,m.getCoords().keepSelectedComponents,["1"])
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
        self.assertRaises(Exception,targetMesh.insertNextCell,NORM_QUAD4,4,targetConn[0:4])
        targetMesh.setMeshDimension(2);
        targetMesh.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
        self.assertRaises(Exception,targetMesh.insertNextCell,NORM_TETRA4,4,targetConn[0:4])
        self.assertRaises(Exception,targetMesh.insertNextCell,NORM_SEG2,2,targetConn[0:2])
        self.assertRaises(Exception,targetMesh.insertNextCell,NORM_POINT1,1,targetConn[0:1])
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
        self.assertRaises(Exception,f2.__div__,f1)
        f3.checkCoherency();
        f1/=f2;
        #self.assertRaises(Exception,f2.__idiv__,f1) # mem leaks
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
        self.assertRaises(da1.rearrange(7),Exception);
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
        self.assertRaises(da2.rearrange(7),Exception);
        #
        da2.rearrange(1);
        self.assertTrue(ptr2==da2.getConstPointer());
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

    def testDARearrange1(self):
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
        self.assertRaises(Exception,a.buildPermutationArr,b)
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
        self.assertEqual(1,b.getArray().getNumberOfComponents());
        self.assertEqual(3,b.getArray().getNumberOfTuples());
        expected1=[0.125,0.25,0.25];
        for i in xrange(3):
            self.assertAlmostEqual(expected1[i],b.getArray().getIJ(0,i),14);
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
    
    def setUp(self):
        pass
    pass

unittest.main()
