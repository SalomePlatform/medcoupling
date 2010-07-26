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

from libMEDCoupling_Swig import *
import unittest
from math import pi
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
        self.failUnless(arr.isEqual(arr3,1e-14))
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
        self.failUnless(mesh.getName()=="mesh1")
        for i in range(nbOfCells):
            mesh.insertNextCell(NORM_QUAD4,4,tab4[4*i:4*(i+1)]);
            pass
        mesh.finishInsertingCells()
        self.failUnless(mesh.getNumberOfCells()==nbOfCells)
        self.failUnless(mesh.getNodalConnectivity().getNbOfElems()==30)
        self.failUnless(mesh.getNodalConnectivityIndex().getNbOfElems()==nbOfCells+1)
        myCoords=DataArrayDouble.New()
        myCoords.setValues(coords,nbOfNodes,3);
        self.failUnless(myCoords.getIJ(3,2)==-0.305)
        mesh.setCoords(myCoords);
        mesh.checkCoherency();
        self.failUnless(mesh.getAllTypes()==[4])
        myFalseConn=DataArrayInt.New()
        myFalseConn.setValues(tab4,6,4)
        self.failUnless(myFalseConn.getIJ(1,1)==3)
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
        self.failUnless(3==mesh.getSpaceDimension())
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
        targetMesh.insertNextCell(NORM_POINT0,0,[]);
        targetMesh.insertNextCell(NORM_POINT0,0,[]);
        targetMesh.insertNextCell(NORM_POINT0,0,[]);
        targetMesh.insertNextCell(NORM_POINT0,0,[]);
        targetMesh.insertNextCell(NORM_POINT0,0,[]);
        targetMesh.insertNextCell(NORM_POINT0,0,[]);
        targetMesh.insertNextCell(NORM_POINT0,0,[]);
        targetMesh.insertNextCell(NORM_POINT0,0,[]);
        targetMesh.finishInsertingCells();
        myCoords=DataArrayDouble.New();
        myCoords.setValues(targetCoords,9,3);
        targetMesh.setCoords(myCoords);
        self.failUnlessEqual(targetMesh.getSpaceDimension(),3)
        self.failUnlessEqual(targetMesh.getNumberOfCells(),8)
        self.failUnlessEqual(targetMesh.getNumberOfNodes(),9)
        self.failUnlessEqual(targetMesh.getMeshDimension(),0)
        pass
    def testMeshM1D(self):
        meshM1D=MEDCouplingUMesh.New();
        ## CPPUNIT_ASSERT_THROW(meshM1D->getMeshDimension(),INTERP_KERNEL::Exception);
        ## CPPUNIT_ASSERT_THROW(meshM1D->getNumberOfNodes(),INTERP_KERNEL::Exception);
        ## CPPUNIT_ASSERT_THROW(meshM1D->getNumberOfCells(),INTERP_KERNEL::Exception);
        ## CPPUNIT_ASSERT_THROW(meshM1D->setMeshDimension(-2),INTERP_KERNEL::Exception);
        ## CPPUNIT_ASSERT_THROW(meshM1D->setMeshDimension(-10),INTERP_KERNEL::Exception);
        ## CPPUNIT_ASSERT_THROW(meshM1D->checkCoherency(),INTERP_KERNEL::Exception);
        meshM1D.setMeshDimension(-1);
        meshM1D.checkCoherency();
        self.failUnlessEqual(meshM1D.getMeshDimension(),-1);
        self.failUnlessEqual(meshM1D.getNumberOfCells(),1);
        ## CPPUNIT_ASSERT_THROW(meshM1D.getNumberOfNodes(),INTERP_KERNEL::Exception);
        ## CPPUNIT_ASSERT_THROW(meshM1D.getSpaceDimension(),INTERP_KERNEL::Exception);
        cpy=meshM1D.clone(True);
        self.failUnless(cpy.isEqual(meshM1D,1e-12));
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
        self.failUnlessEqual(array.getIJ(3,2),7.);
        array2=array.deepCopy();
        self.failUnlessEqual(array2.getIJ(3,2),7.)
        #
        array3=DataArrayInt.New();
        array3.setValues(5*3*[17],5,3);
        self.failUnlessEqual(array3.getIJ(3,2),17);
        array4=array3.deepCopy();
        self.failUnlessEqual(array4.getIJ(3,2),17);
        pass
    def testRevNodal(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1()
        revNodal=DataArrayInt.New();
        revNodalIndx=DataArrayInt.New();
        mesh.getReverseNodalConnectivity(revNodal,revNodalIndx);
        revNodalExpected=[0,0,1,1,2,0,3,0,1,2,3,4,2,4,3,3,4,4];
        revNodalIndexExpected=[0,1,3,5,7,12,14,15,17,18];
        self.failUnlessEqual(revNodal.getNbOfElems(),18)
        self.failUnlessEqual(revNodalIndx.getNbOfElems(),10)
        self.failUnlessEqual(revNodal.getValues(),revNodalExpected)
        self.failUnlessEqual(revNodalIndx.getValues(),revNodalIndexExpected)
        pass
    def testConvertToPolyTypes(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        elts=[1,3];
        mesh.convertToPolyTypes(elts);
        mesh.checkCoherency();
        self.failUnlessEqual(5,mesh.getNumberOfCells());
        self.failUnlessEqual(23,mesh.getNodalConnectivity().getNumberOfTuples());
        expected1=[4, 0, 3, 4, 1, 5, 1, 4, 2, 3, 4, 5, 2, 5, 6, 7, 4, 3, 4, 7, 8, 5, 4]
        self.failUnlessEqual(expected1,mesh.getNodalConnectivity().getValues());
        #
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        mesh.convertToPolyTypes(elts);
        mesh.checkCoherency();
        self.failUnlessEqual(8,mesh.getNumberOfCells());
        self.failUnlessEqual(114,mesh.getNodalConnectivity().getNumberOfTuples());
        mesh.convertToPolyTypes(elts);
        mesh.checkCoherency();
        self.failUnlessEqual(8,mesh.getNumberOfCells());
        self.failUnlessEqual(114,mesh.getNodalConnectivity().getNumberOfTuples());
        pass
    def testDescConn2D(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        desc=DataArrayInt.New();
        descIndx=DataArrayInt.New();
        revDesc=DataArrayInt.New();
        revDescIndx=DataArrayInt.New();
        mesh2=mesh.buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
        mesh2.checkCoherency();
        self.failUnlessEqual(1,mesh2.getMeshDimension());
        self.failUnlessEqual(13,mesh2.getNumberOfCells());
        self.failUnlessEqual(14,revDescIndx.getNbOfElems()); self.failUnlessEqual(14,revDescIndx.getNumberOfTuples());
        self.failUnlessEqual(6,descIndx.getNbOfElems()); self.failUnlessEqual(6,descIndx.getNumberOfTuples());
        self.failUnlessEqual(18,desc.getNbOfElems()); self.failUnlessEqual(18,desc.getNumberOfTuples());
        self.failUnlessEqual(18,revDesc.getNbOfElems()); self.failUnlessEqual(18,revDesc.getNumberOfTuples());
        expected1=[0,1,2,3, 2,4,5, 6,7,4, 8,9,1,10, 11,12,6,9];
        self.failUnlessEqual(expected1,desc.getValues());
        expected2=[0,4,7,10,14,18];
        self.failUnlessEqual(expected2,descIndx.getValues());
        expected3=[0,1,3,5,6,8,9,11,12,13,15,16,17,18];
        self.failUnlessEqual(expected3,revDescIndx.getValues());
        expected4=[0, 0,3, 0,1, 0, 1,2, 1, 2,4, 2, 3, 3,4, 3, 4, 4];
        self.failUnlessEqual(expected4,revDesc.getValues());
        conn=mesh2.getNodalConnectivity();
        connIndex=mesh2.getNodalConnectivityIndex();
        expected5=[0,3,6,9,12,15,18,21,24,27,30,33,36,39];
        self.failUnlessEqual(expected5,connIndex.getValues());
        expected6=[1, 0, 3, 1, 3, 4, 1, 4, 1, 1, 1, 0, 1, 4, 2, 1, 2, 1, 1, 4, 5, 1, 5, 2, 1, 6, 7, 1, 7, 4, 1, 3, 6, 1, 7, 8, 1, 8, 5];
        self.failUnlessEqual(expected6,conn.getValues());
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
        self.failUnlessEqual(1,mesh2.getMeshDimension());
        self.failUnlessEqual(13,mesh2.getNumberOfCells());
        self.failUnlessEqual(14,revDescIndx.getNbOfElems()); self.failUnlessEqual(14,revDescIndx.getNumberOfTuples());
        self.failUnlessEqual(6,descIndx.getNbOfElems()); self.failUnlessEqual(6,descIndx.getNumberOfTuples());
        self.failUnlessEqual(18,desc.getNbOfElems()); self.failUnlessEqual(18,desc.getNumberOfTuples());
        self.failUnlessEqual(18,revDesc.getNbOfElems()); self.failUnlessEqual(18,revDesc.getNumberOfTuples());
        self.failUnlessEqual(expected1,desc.getValues());
        self.failUnlessEqual(expected2,descIndx.getValues());
        self.failUnlessEqual(expected3,revDescIndx.getValues());
        self.failUnlessEqual(expected4,revDesc.getValues());
        conn=mesh2.getNodalConnectivity();
        connIndex=mesh2.getNodalConnectivityIndex();
        self.failUnlessEqual(expected5,connIndex.getValues());
        self.failUnlessEqual(expected6,conn.getValues());
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
        self.failUnlessEqual(2,mesh2.getMeshDimension());
        self.failUnlessEqual(36,mesh2.getNumberOfCells());
        self.failUnlessEqual(37,revDescIndx.getNbOfElems()); self.failUnlessEqual(37,revDescIndx.getNumberOfTuples());
        self.failUnlessEqual(9,descIndx.getNbOfElems()); self.failUnlessEqual(9,descIndx.getNumberOfTuples());
        self.failUnlessEqual(48,desc.getNbOfElems()); self.failUnlessEqual(48,desc.getNumberOfTuples());
        self.failUnlessEqual(48,revDesc.getNbOfElems()); self.failUnlessEqual(48,revDesc.getNumberOfTuples());
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
        
        self.failUnlessEqual(expected1,descIndx.getValues());
        self.failUnlessEqual(expected2,desc.getValues());
        self.failUnlessEqual(expected3,revDescIndx.getValues());
        self.failUnlessEqual(expected4,revDesc.getValues());
        self.failUnlessEqual(expected5,mesh2.getNodalConnectivityIndex().getValues());
        self.failUnlessEqual(expected6,mesh2.getNodalConnectivity().getValues());
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
        self.failUnlessEqual(2,mesh2.getMeshDimension());
        self.failUnlessEqual(36,mesh2.getNumberOfCells());
        self.failUnlessEqual(37,revDescIndx.getNbOfElems()); self.failUnlessEqual(37,revDescIndx.getNumberOfTuples());
        self.failUnlessEqual(9,descIndx.getNbOfElems()); self.failUnlessEqual(9,descIndx.getNumberOfTuples());
        self.failUnlessEqual(48,desc.getNbOfElems()); self.failUnlessEqual(48,desc.getNumberOfTuples());
        self.failUnlessEqual(48,revDesc.getNbOfElems()); self.failUnlessEqual(48,revDesc.getNumberOfTuples());
        self.failUnlessEqual(expected1,descIndx.getValues());
        self.failUnlessEqual(expected2,desc.getValues());
        self.failUnlessEqual(expected3,revDescIndx.getValues());
        self.failUnlessEqual(expected4,revDesc.getValues());
        self.failUnlessEqual(expected5,mesh2.getNodalConnectivityIndex().getValues());
        self.failUnlessEqual(expected7,mesh2.getNodalConnectivity().getValues());
        pass
    def testFindBoundaryNodes(self):
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        boundaryNodes=mesh.findBoundaryNodes();
        expected1=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26];
        self.failUnlessEqual(expected1,boundaryNodes);
        pass
    def testBoundaryMesh(self):
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        mesh2=mesh.buildBoundaryMesh(False);
        self.failUnlessEqual(24,mesh2.getNumberOfCells());
        self.failUnlessEqual(26,mesh2.getNumberOfNodes());
        pass
    def testBuildPartOfMySelf(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        mesh.setName("Toto");
        tab1=[0,4]
        tab2=[0,2,3]
        #
        subMesh=mesh.buildPartOfMySelf(tab1,True);
        self.failUnless(isinstance(subMesh,MEDCouplingUMesh))
        name=subMesh.getName();
        self.failUnlessEqual(2,len(mesh.getAllTypes()));
        self.failUnlessEqual(NORM_TRI3,mesh.getAllTypes()[0]);
        self.failUnlessEqual(NORM_QUAD4,mesh.getAllTypes()[1]);
        self.failUnlessEqual(1,len(subMesh.getAllTypes()));
        self.failUnlessEqual(NORM_QUAD4,subMesh.getAllTypes()[0]);
        self.failUnlessEqual(name,"PartOf_Toto");
        self.failUnlessEqual(2,subMesh.getNumberOfCells());
        subConn=[4,0,3,4,1,4,7,8,5,4];
        subConnIndex=[0,5,10];
        self.failUnlessEqual(10,subMesh.getNodalConnectivity().getNbOfElems());
        self.failUnlessEqual(3,subMesh.getNodalConnectivityIndex().getNbOfElems());
        self.failUnlessEqual(subConn[0:10],subMesh.getNodalConnectivity().getValues());
        self.failUnlessEqual(subConnIndex[0:3],subMesh.getNodalConnectivityIndex().getValues());
        #
        subMesh=mesh.buildPartOfMySelf(tab2[0:3],True);
        self.failUnless(isinstance(subMesh,MEDCouplingUMesh));
        name=subMesh.getName();
        self.failUnlessEqual(2,len(subMesh.getAllTypes()));
        self.failUnlessEqual(NORM_TRI3,subMesh.getAllTypes()[0]);
        self.failUnlessEqual(NORM_QUAD4,subMesh.getAllTypes()[1]);
        self.failUnlessEqual(name,"PartOf_Toto");
        self.failUnlessEqual(3,subMesh.getNumberOfCells());
        subConn2=[4,0,3,4,1,3,4,5,2,4,6,7,4,3]
        subConnIndex2=[0,5,9,14]
        self.failUnlessEqual(14,subMesh.getNodalConnectivity().getNbOfElems());
        self.failUnlessEqual(4,subMesh.getNodalConnectivityIndex().getNbOfElems());
        self.failUnlessEqual(subConn2[0:14],subMesh.getNodalConnectivity().getValues());
        self.failUnlessEqual(subConnIndex2[0:4],subMesh.getNodalConnectivityIndex().getValues());
        pass
    def testBuildPartOfMySelfNode(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        tab1=[5,7]
        subMesh=mesh.buildPartOfMySelfNode(tab1[0:2],True);
        self.failUnless(isinstance(subMesh,MEDCouplingUMesh))
        self.failUnlessEqual(1,len(subMesh.getAllTypes()));
        self.failUnlessEqual(NORM_QUAD4,subMesh.getAllTypes()[0]);
        self.failUnlessEqual(1,subMesh.getNumberOfCells());
        self.failUnlessEqual(5,subMesh.getNodalConnectivity().getNbOfElems());
        self.failUnlessEqual(2,subMesh.getNodalConnectivityIndex().getNbOfElems());
        subConn=[4,7,8,5,4]
        subConnIndex=[0,5]
        self.failUnlessEqual(subConn[0:5],subMesh.getNodalConnectivity().getValues());
        self.failUnlessEqual(subConnIndex[0:2],subMesh.getNodalConnectivityIndex().getValues());
        #
        subMesh=mesh.buildPartOfMySelfNode(tab1[0:2],False);
        self.failUnless(isinstance(subMesh,MEDCouplingUMesh))
        self.failUnlessEqual(2,len(subMesh.getAllTypes()));
        self.failUnlessEqual(NORM_TRI3,subMesh.getAllTypes()[0]);
        self.failUnlessEqual(NORM_QUAD4,subMesh.getAllTypes()[1]);
        self.failUnlessEqual(3,subMesh.getNumberOfCells());
        self.failUnlessEqual(14,subMesh.getNodalConnectivity().getNbOfElems());
        self.failUnlessEqual(4,subMesh.getNodalConnectivityIndex().getNbOfElems());
        subConn2=[3,4,5,2,4,6,7,4,3,4,7,8,5,4]
        subConnIndex2=[0,4,9,14]
        self.failUnlessEqual(subConn2[0:14],subMesh.getNodalConnectivity().getValues());
        self.failUnlessEqual(subConnIndex2[0:4],subMesh.getNodalConnectivityIndex().getValues());
        #testing the case where length of tab2 is greater than max number of node per cell.
        tab2=[0,3,2,1,4,5,6]
        subMesh=mesh.buildPartOfMySelfNode(tab2[0:7],True);
        self.failUnless(isinstance(subMesh,MEDCouplingUMesh))
        self.failUnlessEqual(2,len(subMesh.getAllTypes()));
        self.failUnlessEqual(NORM_TRI3,subMesh.getAllTypes()[0]);
        self.failUnlessEqual(NORM_QUAD4,subMesh.getAllTypes()[1]);
        self.failUnlessEqual(3,subMesh.getNumberOfCells());
        pass
    def testZipCoords(self):
        mesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        self.failUnlessEqual(2,len(mesh.getAllTypes()));
        self.failUnlessEqual(2,mesh.getSpaceDimension());
        self.failUnlessEqual(9,mesh.getNumberOfNodes());
        self.failUnlessEqual(5,mesh.getNumberOfCells());
        oldConn=mesh.getNodalConnectivity().getValues()[0:mesh.getNodalConnectivity().getNbOfElems()];
        oldConnIndex=mesh.getNodalConnectivityIndex().getValues()[0:mesh.getNumberOfCells()+1]
        oldCoords=mesh.getCoords();
        mesh.zipCoords();
        self.failUnlessEqual(2,len(mesh.getAllTypes()));
        self.failUnlessEqual(2,mesh.getSpaceDimension());
        self.failUnlessEqual(9,mesh.getNumberOfNodes());
        self.failUnlessEqual(5,mesh.getNumberOfCells());
        self.failUnlessEqual(mesh.getCoords().getValues()[0:2*9],oldCoords.getValues());
        self.failUnlessEqual(oldConn,mesh.getNodalConnectivity().getValues());
        self.failUnlessEqual(oldConnIndex,mesh.getNodalConnectivityIndex().getValues());
        #
        tab1=[0,4]
        subMesh=mesh.buildPartOfMySelf(tab1,True);
        self.failUnless(isinstance(subMesh,MEDCouplingUMesh))
        traducer=subMesh.zipCoordsTraducer();
        expectedTraducer=[0,1,-1,2,3,4,-1,5,6]
        self.failUnlessEqual(expectedTraducer,traducer.getValues());
        self.failUnlessEqual(NORM_QUAD4,subMesh.getAllTypes()[0]);
        self.failUnlessEqual(2,subMesh.getNumberOfCells());
        subConn=[4,0,2,3,1,4,5,6,4,3]
        subConnIndex=[0,5,10]
        self.failUnlessEqual(7,subMesh.getNumberOfNodes());
        self.failUnlessEqual(10,subMesh.getNodalConnectivity().getNbOfElems());
        self.failUnlessEqual(3,subMesh.getNodalConnectivityIndex().getNbOfElems());
        self.failUnlessEqual(subConn,subMesh.getNodalConnectivity().getValues());
        self.failUnlessEqual(subConnIndex,subMesh.getNodalConnectivityIndex().getValues());
        #
        subMesh=mesh.buildPartOfMySelf(tab1,False);
        self.failUnless(isinstance(subMesh,MEDCouplingUMesh))
        self.failUnlessEqual(NORM_QUAD4,subMesh.getAllTypes()[0]);
        self.failUnlessEqual(2,subMesh.getNumberOfCells());
        self.failUnlessEqual(7,subMesh.getNumberOfNodes());
        self.failUnlessEqual(10,subMesh.getNodalConnectivity().getNbOfElems());
        self.failUnlessEqual(3,subMesh.getNodalConnectivityIndex().getNbOfElems());
        self.failUnlessEqual(subConn,subMesh.getNodalConnectivity().getValues());
        self.failUnlessEqual(subConnIndex,subMesh.getNodalConnectivityIndex().getValues());
        pass
    def testZipConnectivity(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells1=[2,3,4]
        m3=m2.buildPartOfMySelf(cells1,True);
        self.failUnless(isinstance(m3,MEDCouplingUMesh))
        m4=MEDCouplingDataForTest.build2DSourceMesh_1();
        m5=MEDCouplingUMesh.mergeUMeshes(m1,m3);
        m6=MEDCouplingUMesh.mergeUMeshes(m5,m4);
        #
        self.failUnlessEqual(10,m6.getNumberOfCells());
        self.failUnlessEqual(22,m6.getNumberOfNodes());
        (arr,areNodesMerged)=m6.mergeNodes(1e-13);
        self.failUnless(areNodesMerged);
        self.failUnlessEqual(10,m6.getNumberOfCells());
        self.failUnlessEqual(9,m6.getNumberOfNodes());
        #
        arr=m6.zipConnectivityTraducer(0);
        self.failUnlessEqual(7,m6.getNumberOfCells());
        m7=m6.clone(True);
        arr=m6.zipConnectivityTraducer(0);
        self.failUnless(m7.isEqual(m6,1e-12));
        self.failUnlessEqual(7,m6.getNumberOfCells());
        pass
    def testEqualMesh(self):
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1();
        mesh2=MEDCouplingDataForTest.build2DTargetMesh_1();
        #
        self.failUnless(mesh1.isEqual(mesh1,1e-12));
        #
        self.failUnless(mesh1.isEqual(mesh2,1e-12));
        self.failUnless(mesh2.isEqual(mesh1,1e-12));
        pt=mesh2.getCoords().getValues();
        tmp=pt[1]
        mesh2.getCoords().setIJ(0,1,5.999);
        self.failUnless(not mesh1.isEqual(mesh2,1e-12));
        self.failUnless(not mesh2.isEqual(mesh1,1e-12));
        mesh2.getCoords().setIJ(0,1,tmp);
        self.failUnless(mesh1.isEqual(mesh2,1e-12));
        self.failUnless(mesh2.isEqual(mesh1,1e-12));
        #
        pt2=mesh1.getNodalConnectivity().getValues();
        mesh1.getNodalConnectivity().setIJ(5,0,pt2[5]+1);
        self.failUnless(not mesh1.isEqual(mesh2,1e-12));
        self.failUnless(not mesh2.isEqual(mesh1,1e-12));
        mesh1.getNodalConnectivity().setIJ(5,0,pt2[5]);
        self.failUnless(mesh1.isEqual(mesh2,1e-12));
        self.failUnless(mesh2.isEqual(mesh1,1e-12));
        #
        pt2=mesh1.getNodalConnectivityIndex().getValues();
        mesh1.getNodalConnectivityIndex().setIJ(1,0,pt2[1]+1);
        self.failUnless(not mesh1.isEqual(mesh2,1e-12));
        self.failUnless(not mesh2.isEqual(mesh1,1e-12));
        mesh1.getNodalConnectivityIndex().setIJ(1,0,pt2[1]);
        self.failUnless(mesh1.isEqual(mesh2,1e-12));
        self.failUnless(mesh2.isEqual(mesh1,1e-12));
        #
        tmp3=mesh1.getName();
        mesh1.setName("lllll");
        self.failUnless(not mesh1.isEqual(mesh2,1e-12));
        self.failUnless(not mesh2.isEqual(mesh1,1e-12));
        mesh1.setName(tmp3);
        self.failUnless(mesh1.isEqual(mesh2,1e-12));
        self.failUnless(mesh2.isEqual(mesh1,1e-12));
        #
        tmp3=mesh2.getCoords().getInfoOnComponent(1);
        mesh2.getCoords().setInfoOnComponent(1,"kkkkkk");
        self.failUnless(not mesh1.isEqual(mesh2,1e-12));
        self.failUnless(not mesh2.isEqual(mesh1,1e-12));
        mesh2.getCoords().setInfoOnComponent(1,tmp3);
        self.failUnless(mesh1.isEqual(mesh2,1e-12));
        self.failUnless(mesh2.isEqual(mesh1,1e-12));
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
        self.failUnless(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        #
        fieldOnNodes1=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        self.failUnless(not fieldOnCells1.isEqual(fieldOnNodes1,1e-12,1e-15));
        self.failUnless(not fieldOnNodes1.isEqual(fieldOnCells1,1e-12,1e-15));
        #
        fieldOnCells2=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        self.failUnless(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        self.failUnless(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells1.setTime(4.,6,7);
        self.failUnless(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells2.setTime(4.,6,7);
        self.failUnless(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells1.setName("Power");
        self.failUnless(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells2.setName("Power");
        self.failUnless(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        #
        fieldOnCells1.setMesh(mesh1);
        self.failUnless(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells2.setMesh(mesh1);
        self.failUnless(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        arr=DataArrayDouble.New();
        arr.setName("popo");
        arr.setValues(mesh1.getNumberOfCells()*3*[6.],mesh1.getNumberOfCells(),3);
        fieldOnCells1.setArray(arr);
        self.failUnless(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        fieldOnCells2.setArray(arr);
        self.failUnless(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        #
        arr2=arr.deepCopy();
        fieldOnCells2.setArray(arr2);
        self.failUnless(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        arr.setIJ(1,2,6.1);
        self.failUnless(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        arr.setIJ(1,2,6.);
        self.failUnless(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        arr2.setName("popo2");
        self.failUnless(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        #
        arr2.setName("popo");
        self.failUnless(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        #
        arr2.setInfoOnComponent(2,"jjj");
        self.failUnless(not fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(not fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        arr.setInfoOnComponent(2,"jjj");
        self.failUnless(fieldOnCells1.isEqual(fieldOnCells2,1e-12,1e-15));
        self.failUnless(fieldOnCells2.isEqual(fieldOnCells1,1e-12,1e-15));
        pass

    def testNatureChecking(self):
        field=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        field.setNature(Integral);
        field.setNature(ConservativeVolumic);
        field.setNature(IntegralGlobConstraint);
        field=MEDCouplingFieldDouble.New(ON_NODES,NO_TIME);
        field.setNature(ConservativeVolumic);
        #self.failUnless_THROW(field.setNature(Integral),INTERP_KERNEL::Exception);
        #self.failUnless_THROW(field.setNature(IntegralGlobConstraint),INTERP_KERNEL::Exception);
        pass

    def testBuildSubMeshData(self):
        targetMesh=MEDCouplingDataForTest.build2DTargetMesh_1()
        #check buildSubMesh on field on cells
        fieldCells=MEDCouplingFieldDouble.New(ON_CELLS,NO_TIME);
        fieldCells.setMesh(targetMesh);
        elts=[1,2,4]
        ret1,di=fieldCells.buildSubMeshData(elts);
        self.failUnless(isinstance(ret1,MEDCouplingUMesh))
        self.failUnlessEqual(3,ret1.getNumberOfCells());
        self.failUnlessEqual(6,ret1.getNumberOfNodes());
        self.failUnlessEqual(3,di.getNumberOfTuples());
        self.failUnlessEqual(1,di.getNumberOfComponents());
        toCheck=di.getValues();
        self.failUnless(elts,toCheck);
        #check buildSubMesh on field on nodes
        fieldNodes=MEDCouplingFieldDouble.New(ON_NODES,NO_TIME);
        fieldNodes.setMesh(targetMesh);
        ret2,di=fieldNodes.buildSubMeshData(elts);
        self.failUnless(isinstance(ret2,MEDCouplingUMesh))
        self.failUnlessEqual(3,ret2.getNumberOfCells());
        self.failUnlessEqual(6,ret2.getNumberOfNodes());
        self.failUnlessEqual(6,di.getNumberOfTuples());
        self.failUnlessEqual(1,di.getNumberOfComponents());
        toCheck=di.getValues();
        expected=[1,2,4,5,7,8]
        self.failUnlessEqual(expected,toCheck);
        pass
    def testExtrudedMesh1(self):
        mesh3D,mesh2D=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        ext=MEDCouplingExtrudedMesh.New(mesh3D,mesh2D,1);
        self.failUnlessEqual(18,ext.getNumberOfCells());
        self.failUnlessEqual(60,ext.getNumberOfNodes());
        ids3D=ext.getMesh3DIds();
        ids3DExpected=[5,4,3,2,1,0, 11,10,9,8,7,6, 17,16,15,14,13,12]
        self.failUnlessEqual(18,ids3D.getNumberOfTuples());
        self.failUnlessEqual(1,ids3D.getNumberOfComponents());
        self.failUnlessEqual(ids3DExpected,ids3D.getValues());
        mesh1D=ext.getMesh1D();
        self.failUnlessEqual(4,mesh1D.getNumberOfNodes());
        self.failUnlessEqual(3,mesh1D.getNumberOfCells());
        mesh1DExpected=[0.66666666666666663, 1.4583333333333333, 0, 0.66666666666666663,
                        1.4583333333333333, 1, 0.66666666666666663, 1.4583333333333333,
                        2, 0.66666666666666663, 1.4583333333333333, 3]
        mesh1DCoords=mesh1D.getCoords();
        self.failUnlessEqual(4,mesh1DCoords.getNumberOfTuples());
        self.failUnlessEqual(3,mesh1DCoords.getNumberOfComponents());
        self.failUnlessEqual(mesh1DExpected,mesh1DCoords.getValues());
        conn1D=mesh1D.getNodalConnectivity();
        self.failUnlessEqual(9,conn1D.getNumberOfTuples());
        self.failUnlessEqual(1,conn1D.getNumberOfComponents());
        conn1DExpected=[1,0,1,1,1,2,1,2,3]
        self.failUnlessEqual(conn1DExpected,conn1D.getValues());
        pass

    def testExtrudedMesh3(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m1.changeSpaceDimension(3);
        m2=MEDCouplingDataForTest.buildCU1DMesh_U();
        m2.changeSpaceDimension(3);
        center=[0.,0.,0.]
        vector=[0.,1.,0.]
        m2.rotate(center,vector,-pi/2.);
        m3=m1.buildExtrudedMeshFromThis(m2,0);
        #
        m4=MEDCouplingExtrudedMesh.New(m3,m1,0);
        self.failUnlessEqual(15,m4.getNumberOfCells());
        self.failUnlessEqual(5,m4.getMesh2D().getNumberOfCells());
        self.failUnlessEqual(3,m4.getMesh1D().getNumberOfCells());
        m3DIds=m4.getMesh3DIds().getValues();
        self.failUnlessEqual(range(15),m3DIds);
        #some random in cells to check that extrusion alg find it correctly
        expected1=[1,3,2,0,6,5,7,10,11,8,12,9,14,13,4]
        m3.renumberCells(expected1,False);
        m4=MEDCouplingExtrudedMesh.New(m3,m1,0);
        self.failUnlessEqual(15,m4.getNumberOfCells());
        self.failUnlessEqual(5,m4.getMesh2D().getNumberOfCells());
        self.failUnlessEqual(3,m4.getMesh1D().getNumberOfCells());
        m3DIds=m4.getMesh3DIds().getValues();
        self.failUnlessEqual(expected1,m3DIds);
        #play with polygons and polyedrons
        cells=[2,3]
        m1.convertToPolyTypes(cells);
        m3=m1.buildExtrudedMeshFromThis(m2,0);
        self.failUnlessEqual(NORM_HEXA8,m3.getTypeOfCell(0));
        self.failUnlessEqual(NORM_PENTA6,m3.getTypeOfCell(1));
        self.failUnlessEqual(NORM_POLYHED,m3.getTypeOfCell(2));
        self.failUnlessEqual(NORM_POLYHED,m3.getTypeOfCell(3));
        self.failUnlessEqual(NORM_HEXA8,m3.getTypeOfCell(4));
        m3.renumberCells(expected1,False);
        m4=MEDCouplingExtrudedMesh.New(m3,m1,0);
        self.failUnlessEqual(15,m4.getNumberOfCells());
        self.failUnlessEqual(5,m4.getMesh2D().getNumberOfCells());
        self.failUnlessEqual(3,m4.getMesh1D().getNumberOfCells());
        m3DIds=m4.getMesh3DIds().getValues();
        self.failUnlessEqual(expected1,m3DIds);
        pass

    def testFindCommonNodes(self):
        targetMesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        comm,commI=targetMesh.findCommonNodes(1e-10);
        self.failUnlessEqual(1,commI.getNumberOfTuples());
        self.failUnlessEqual(0,comm.getNumberOfTuples());
        o2n,newNbOfNodes=targetMesh.buildNewNumberingFromCommNodesFrmt(comm,commI);
        self.failUnlessEqual(27,newNbOfNodes);
        self.failUnlessEqual(27,o2n.getNumberOfTuples());
        o2nExp1=range(27)
        self.failUnlessEqual(o2nExp1,o2n.getValues());
        #
        targetMesh=MEDCouplingDataForTest.build3DTargetMeshMergeNode_1();
        self.failUnlessEqual(31,targetMesh.getNumberOfNodes());
        comm,commI=targetMesh.findCommonNodes(1e-10);
        self.failUnlessEqual(3,commI.getNumberOfTuples());
        self.failUnlessEqual(6,comm.getNumberOfTuples());
        commExpected=[1,27,28,29,23,30]
        commIExpected=[0,4,6]
        self.failUnlessEqual(commExpected,comm.getValues());
        self.failUnlessEqual(commIExpected,commI.getValues());
        o2n,newNbOfNodes=targetMesh.buildNewNumberingFromCommNodesFrmt(comm,commI);
        self.failUnlessEqual(31,o2n.getNumberOfTuples());
        self.failUnlessEqual(27,newNbOfNodes);
        o2nExp2=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
                 21,22,23,24,25,26,1,1,1,23]
        self.failUnlessEqual(o2nExp2,o2n.getValues());
        #
        targetMesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        time=targetMesh.getTimeOfThis();
        o2n,areNodesMerged=targetMesh.mergeNodes(1e-10);
        targetMesh.updateTime();
        self.failUnlessEqual(time,targetMesh.getTimeOfThis());
        self.failUnless(not areNodesMerged);
        #
        targetMesh=MEDCouplingDataForTest.build3DTargetMeshMergeNode_1();
        time=targetMesh.getTimeOfThis();
        o2n,areNodesMerged=targetMesh.mergeNodes(1e-10);
        targetMesh.updateTime();
        self.failUnless(time!=targetMesh.getTimeOfThis());
        self.failUnless(areNodesMerged);
        connExp=[18,0,1,4,3,9,10,13,12, 18,1,2,5,4,10,11,14,13, 18,3,4,7,6,12,13,16,15,
                 18,4,5,8,7,13,14,17,16,
                 18,9,10,13,12,18,19,22,21, 18,10,11,14,13,19,20,23,22, 18,12,13,16,15,21,22,25,24,
                 18,13,14,17,16,22,23,26,25]
        self.failUnlessEqual(72,targetMesh.getNodalConnectivity().getNumberOfTuples());
        self.failUnlessEqual(connExp,targetMesh.getNodalConnectivity().getValues());
        self.failUnlessEqual(27,targetMesh.getCoords().getNumberOfTuples());
        coordsExp=[ 0., 0., 0., 50., 0., 0. , 200., 0., 0.  , 0., 50., 0., 50., 50., 0. ,
                    200., 50., 0.,   0., 200., 0., 50., 200., 0. , 200., 200., 0. ,
                    0., 0., 50., 50., 0., 50. , 200., 0., 50.  , 0., 50., 50., 50.,
                    50., 50. , 200., 50., 50.,   0., 200., 50., 50., 200., 50. ,
                    200., 200., 50. , 0., 0., 200., 50., 0., 200. , 200., 0., 200.  
                    , 0., 50., 200., 50., 50., 200. , 200., 50., 200., 
                    0., 200., 200., 50., 200., 200. , 200., 200., 200. ]
        self.failUnlessEqual(coordsExp,targetMesh.getCoords().getValues());
        # 2D
        targetMesh=MEDCouplingDataForTest.build2DTargetMeshMergeNode_1();
        self.failUnlessEqual(18,targetMesh.getNumberOfNodes());
        time=targetMesh.getTimeOfThis();
        o2n,areNodesMerged=targetMesh.mergeNodes(1e-10);
        self.failUnless(time!=targetMesh.getTimeOfThis());
        self.failUnless(areNodesMerged);
        self.failUnlessEqual(9,targetMesh.getNumberOfNodes());
        connExp2=[4,0,4,3,1, 3,1,3,2, 3,3,5,2, 4,4,6,7,3, 4,7,8,5,3]
        self.failUnlessEqual(23,targetMesh.getNodalConnectivity().getNumberOfTuples());
        self.failUnlessEqual(connExp2,targetMesh.getNodalConnectivity().getValues());
        coordsExp2=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, 0.2,0.2, -0.3,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7]
        self.failUnlessEqual(9,targetMesh.getCoords().getNumberOfTuples());
        self.failUnlessEqual(coordsExp2,targetMesh.getCoords().getValues());
        pass

    def testCheckButterflyCells(self):
        sourceMesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells=sourceMesh.checkButterflyCells();
        self.failUnlessEqual(0,len(cells));
        conn=sourceMesh.getNodalConnectivity()
        tmp=conn.getIJ(15,0)
        conn.setIJ(15,0,conn.getIJ(16,0))
        conn.setIJ(16,0,tmp)
        cells=sourceMesh.checkButterflyCells();
        self.failUnlessEqual(1,len(cells));
        self.failUnlessEqual(3,cells[0]);
        tmp=conn.getIJ(15,0)
        conn.setIJ(15,0,conn.getIJ(16,0))
        conn.setIJ(16,0,tmp)
        cells=sourceMesh.checkButterflyCells();
        self.failUnlessEqual(0,len(cells));
        # 3D surf
        sourceMesh=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        cells=sourceMesh.checkButterflyCells();
        self.failUnlessEqual(0,len(cells));
        conn=sourceMesh.getNodalConnectivity()
        tmp=conn.getIJ(15,0)
        conn.setIJ(15,0,conn.getIJ(16,0))
        conn.setIJ(16,0,tmp)
        cells=sourceMesh.checkButterflyCells();
        self.failUnlessEqual(1,len(cells));
        self.failUnlessEqual(3,cells[0]);
        tmp=conn.getIJ(15,0)
        conn.setIJ(15,0,conn.getIJ(16,0))
        conn.setIJ(16,0,tmp)
        cells=sourceMesh.checkButterflyCells();
        self.failUnlessEqual(0,len(cells));
        pass

    def testMergeMesh1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DSourceMesh_1();
        vec=[1.,0.]
        m2.translate(vec);
        m3=m1.mergeMyselfWith(m2);
        self.failUnless(isinstance(m3,MEDCouplingUMesh));
        m3.checkCoherency();
        m4=MEDCouplingDataForTest.build2DTargetMeshMerged_1();
        self.failUnless(m3.isEqual(m4,1.e-12));
        da,isMerged=m3.mergeNodes(1.e-12);
        self.failUnlessEqual(11,m3.getNumberOfNodes());
        self.failUnless(isMerged);
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
        m4=MEDCouplingUMesh.mergeUMeshesOnSameCoords(meshes);
        m4.checkCoherency();
        self.failUnlessEqual(15,m4.getNumberOfCells());
        cells1=[0,1,2,3,4]
        m1_1=m4.buildPartOfMySelf(cells1,True);
        m1_1.setName(m1.getName());
        self.failUnless(m1.isEqual(m1_1,1e-12));
        cells2=[5,6,7,8,9]
        m2_1=m4.buildPartOfMySelf(cells2,True);
        m2_1.setName(m2.getName());
        self.failUnless(m2.isEqual(m2_1,1e-12));
        cells3=[10,11,12,13,14]
        m3_1=m4.buildPartOfMySelf(cells3,True);
        m3_1.setName(m3.getName());
        self.failUnless(m3.isEqual(m3_1,1e-12));
        pass

    def testMergeField1(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DSourceMesh_1();
        vec=[1.,0.]
        m2.translate(vec);
        f1=m1.getMeasureField(True);
        f2=m2.getMeasureField(True);
        f3=MEDCouplingFieldDouble.mergeFields(f1,f2);
        f3.checkCoherency();
        m4=MEDCouplingDataForTest.build2DTargetMeshMerged_1();
        self.failUnless(f3.getMesh().isEqual(m4,1.e-12));
        name=f3.getName();
        self.failUnlessEqual(name,"MeasureOfMesh_");
        self.failUnlessEqual(f3.getTypeOfField(),ON_CELLS);
        self.failUnlessEqual(f3.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(1,f3.getNumberOfComponents());
        self.failUnlessEqual(7,f3.getNumberOfTuples());
        values=[0.25,0.125,0.125,0.25,0.25,0.5,0.5]
        tmp=f3.getArray().getValues();
        self.failUnlessEqual(len(values),len(tmp))
        for i in xrange(7):
            self.failUnless(abs(values[i]-tmp[i])<1e-12)
            pass
        pass

    def testFillFromAnalytic(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();             
        f1=m.fillFromAnalytic(ON_CELLS,1,"x+y");
        f1.checkCoherency();                    
        self.failUnlessEqual(f1.getTypeOfField(),ON_CELLS);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(1,f1.getNumberOfComponents());
        self.failUnlessEqual(5,f1.getNumberOfTuples());
        values1=[-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9]
        tmp=f1.getArray().getValues();
        self.failUnlessEqual(len(values1),len(tmp))
        for i in xrange(len(tmp)):
            self.failUnless(abs(tmp[i]-values1[i])<1.e-12)
            pass
        #
        f1=m.fillFromAnalytic(ON_NODES,1,"x+y");
        f1.checkCoherency();
        self.failUnlessEqual(f1.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(1,f1.getNumberOfComponents());
        self.failUnlessEqual(9,f1.getNumberOfTuples());
        values2=[-0.6,-0.1,0.4,-0.1,0.4,0.9,0.4,0.9,1.4]
        tmp=f1.getArray().getValues();
        self.failUnlessEqual(len(values2),len(tmp))
        for i in xrange(len(tmp)):
            self.failUnless(abs(tmp[i]-values2[i])<1.e-12)
            pass
        #
        f1=m.fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+(2*(x+y))*JVec");
        f1.checkCoherency();
        self.failUnlessEqual(f1.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(2,f1.getNumberOfComponents());
        self.failUnlessEqual(9,f1.getNumberOfTuples());
        values3=[-0.6,-1.2,-0.1,-0.2,0.4,0.8,-0.1,-0.2,0.4,0.8,0.9,1.8,0.4,0.8,0.9,1.8,1.4,2.8]
        tmp=f1.getArray().getValues();
        self.failUnlessEqual(len(values3),len(tmp))
        for i in xrange(len(tmp)):
            self.failUnless(abs(tmp[i]-values3[i])<1.e-12)
            pass
        values4=f1.accumulate();
        self.failUnless(abs(3.6-values4[0])<1.e-12);
        self.failUnless(abs(7.2-values4[1])<1.e-12);
        values4=f1.measureAccumulate(True);
        self.failUnless(abs(0.5-values4[0])<1.e-12);
        self.failUnless(abs(1.-values4[1])<1.e-12);
        #
        ## self.failUnlessEqual_THROW(f1=m.fillFromAnalytic(ON_NODES,1,func3),Exception);
        pass

    def testFillFromAnalytic2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=m.fillFromAnalytic(ON_CELLS,1,"y+x");
        f1.checkCoherency();
        self.failUnlessEqual(f1.getTypeOfField(),ON_CELLS);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(1,f1.getNumberOfComponents());
        self.failUnlessEqual(5,f1.getNumberOfTuples());
        values1=[-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9]
        tmp=f1.getArray().getValues();
        self.failUnlessEqual(len(values1),len(tmp))
        for i in xrange(len(values1)):
            self.failUnless(abs(values1[i]-tmp[i])<1.e-12);
            pass
        #
        f1=m.fillFromAnalytic(ON_NODES,1,"y+2*x");
        f1.checkCoherency();
        self.failUnlessEqual(f1.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(1,f1.getNumberOfComponents());
        self.failUnlessEqual(9,f1.getNumberOfTuples());
        values2=[-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1]
        tmp=f1.getArray().getValues();
        self.failUnlessEqual(len(values2),len(tmp))
        for i in xrange(len(values2)):
            self.failUnless(abs(values2[i]-tmp[i])<1.e-12);
            pass
        f1=m.fillFromAnalytic(ON_NODES,1,"2.*x+y");
        f1.checkCoherency();
        self.failUnlessEqual(f1.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(1,f1.getNumberOfComponents());
        self.failUnlessEqual(9,f1.getNumberOfTuples());
        tmp=f1.getArray().getValues();
        values2Bis=[-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1]
        self.failUnlessEqual(len(values2Bis),len(tmp))
        for i in xrange(len(values2Bis)):
            self.failUnless(abs(values2Bis[i]-tmp[i])<1.e-12);
            pass
        #
        f1=m.fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+2*(x+y)*JVec");
        f1.checkCoherency();
        self.failUnlessEqual(f1.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(2,f1.getNumberOfComponents());
        self.failUnlessEqual(9,f1.getNumberOfTuples());
        values3=[-0.6,-1.2,-0.1,-0.2,0.4,0.8,-0.1,-0.2,0.4,0.8,0.9,1.8,0.4,0.8,0.9,1.8,1.4,2.8]
        tmp=f1.getArray().getValues();
        self.failUnlessEqual(len(values3),len(tmp))
        for i in xrange(len(values3)):
            self.failUnless(abs(values3[i]-tmp[i])<1.e-12);
            pass
        values4=f1.accumulate();
        self.failUnless(abs(3.6-values4[0])<1.e-12);
        self.failUnless(abs(7.2-values4[1])<1.e-12);
        values4=f1.measureAccumulate(True);
        self.failUnless(abs(0.5-values4[0])<1.e-12);
        self.failUnless(abs(1.-values4[1])<1.e-12);
        pass

    def testApplyFunc(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=m.fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+(2*(x+y))*JVec");
        f1.checkCoherency();
        self.failUnlessEqual(f1.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(2,f1.getNumberOfComponents());
        self.failUnlessEqual(9,f1.getNumberOfTuples());
        f1.applyFunc(1,"x+y");
        self.failUnlessEqual(f1.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(1,f1.getNumberOfComponents());
        self.failUnlessEqual(9,f1.getNumberOfTuples());
        values1=[-1.8,-0.3,1.2,-0.3,1.2,2.7,1.2,2.7,4.2]
        tmp=f1.getArray().getValues();
        self.failUnlessEqual(len(values1),len(tmp))
        for i in xrange(len(tmp)):
            self.failUnless(abs(tmp[i]-values1[i])<1.e-12)
            pass
        pass

    def testApplyFunc2(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        f1=m.fillFromAnalytic(ON_NODES,2,"(x+y)*IVec+2*(x+y)*JVec");
        f1.checkCoherency();
        self.failUnlessEqual(f1.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(2,f1.getNumberOfComponents());
        self.failUnlessEqual(9,f1.getNumberOfTuples());
        #
        f2=f1.clone(True);
        f2.applyFunc("abs(u)^2.4+2*u");
        self.failUnlessEqual(f1.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(2,f1.getNumberOfComponents());
        self.failUnlessEqual(9,f1.getNumberOfTuples());
        values2=[-0.9065304805418678, -0.85105859001709905, -0.19601892829446504, -0.37898777756476987,
                 0.91090317490482353, 2.1853504664669781, -0.19601892829446504, -0.37898777756476987,
                 0.91090317490482353, 2.1853504664669781, 2.5765725275664879, 7.6987743736515295,
                 0.91090317490482353, 2.1853504664669781, 2.5765725275664879, 7.6987743736515295,
                 5.0423700574830965, 17.435300118916864]
        tmp=f2.getArray().getValues();
        self.failUnlessEqual(len(tmp),len(values2))
        for i in xrange(len(tmp)):
            self.failUnless(abs(tmp[i]-values2[i])<1.e-12)
            pass
        #
        f1.applyFunc(1,"x+y");
        self.failUnlessEqual(f1.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        self.failUnlessEqual(1,f1.getNumberOfComponents());
        self.failUnlessEqual(9,f1.getNumberOfTuples());
        values1=[-1.8,-0.3,1.2,-0.3,1.2,2.7,1.2,2.7,4.2]
        tmp=f1.getArray().getValues();
        self.failUnlessEqual(len(tmp),len(values1))
        for i in xrange(len(tmp)):
            self.failUnless(abs(tmp[i]-values1[i])<1.e-12)
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
        self.failUnlessEqual(f3.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f3.getTimeDiscretization(),NO_TIME);
        values1=[-1.2,-0.2,0.8,-0.2,0.8,1.8,0.8,1.8,2.8]
        tmp=f3.getArray().getValues();
        self.failUnlessEqual(len(values1),len(tmp))
        for i in xrange(len(tmp)):
            self.failUnless(abs(tmp[i]-values1[i])<1.e-12)
            pass
        #
        f3=f1*f2;
        f3.checkCoherency();
        self.failUnlessEqual(f3.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f3.getTimeDiscretization(),NO_TIME);
        values2=[0.36,0.01,0.16,0.01,0.16,0.81,0.16,0.81,1.96]
        tmp=f3.getArray().getValues();
        self.failUnlessEqual(len(values2),len(tmp))
        for i in xrange(len(tmp)):
            self.failUnless(abs(tmp[i]-values2[i])<1.e-12)
            pass
        #
        f3=f1+f2;
        f4=f1-f3;
        f4.checkCoherency();
        self.failUnlessEqual(f4.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f4.getTimeDiscretization(),NO_TIME);
        values3=[0.6,0.1,-0.4,0.1,-0.4,-0.9,-0.4,-0.9,-1.4]
        tmp=f4.getArray().getValues();
        self.failUnlessEqual(len(values3),len(tmp))
        for i in xrange(len(tmp)):
            self.failUnless(abs(tmp[i]-values3[i])<1.e-12)
            pass
        #
        f3=f1+f2;
        f4=f3/f2;
        f4.checkCoherency();
        self.failUnlessEqual(f4.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f4.getTimeDiscretization(),NO_TIME);
        tmp=f4.getArray().getValues();
        for i in xrange(len(tmp)):
            self.failUnless(abs(tmp[i]-2.)<1.e-12)
            pass
        #
        f4=f2.buildNewTimeReprFromThis(ONE_TIME,False);
        f4.checkCoherency();
        self.failUnlessEqual(f4.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f4.getTimeDiscretization(),ONE_TIME);
        ## self.failUnlessEqual_THROW(f3=f1+f4,Exception);
        f5=f4.buildNewTimeReprFromThis(NO_TIME,False);
        self.failUnlessEqual(f5.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f5.getTimeDiscretization(),NO_TIME);
        f3=f1+f5;
        tmp=f3.getArray().getValues();
        values4=[-1.2,-0.2,0.8,-0.2,0.8,1.8,0.8,1.8,2.8]
        self.failUnlessEqual(len(values3),len(tmp))
        for i in xrange(len(tmp)):
            self.failUnless(abs(tmp[i]-values4[i])<1.e-12)
            pass
        #
        f4=f2.buildNewTimeReprFromThis(ONE_TIME,True);
        f4.checkCoherency();
        self.failUnlessEqual(f4.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f4.getTimeDiscretization(),ONE_TIME);
        ## self.failUnlessEqual_THROW(f3=f1+f4,Exception);
        f5=f4.buildNewTimeReprFromThis(NO_TIME,True);
        self.failUnlessEqual(f5.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f5.getTimeDiscretization(),NO_TIME);
        f3=f1+f5;
        tmp=f3.getArray().getValues();
        values5=[-1.2,-0.2,0.8,-0.2,0.8,1.8,0.8,1.8,2.8]
        self.failUnlessEqual(len(values5),len(tmp))
        for i in xrange(len(tmp)):
            self.failUnless(abs(tmp[i]-values5[i])<1.e-12)
            pass
        pass

    def testOperationsOnFields2(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        f1=m.fillFromAnalytic(ON_NODES,1,"x+y+z");
        f2=m.fillFromAnalytic(ON_NODES,1,"a*a+b+c*c");
        f3=f1/f2;
        f3.checkCoherency();
        self.failUnlessEqual(f3.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f3.getTimeDiscretization(),NO_TIME);
        expected1=[-2.4999999999999991, 1.2162162162162162, 0.77868852459016391,
                   0.7407407407407407, 1.129032258064516, 0.81632653061224492,
                   0.86538461538461531, 1.0919540229885056, 0.84302325581395343]
        self.failUnlessEqual(1,f3.getNumberOfComponents());
        self.failUnlessEqual(9,f3.getNumberOfTuples());
        val=f3.getArray().getValues();
        for i in xrange(9):
            self.failUnless(abs(expected1[i]-val[i])<1.e-12);
        #
        f1=m.buildOrthogonalField();
        f2=m.fillFromAnalytic(ON_CELLS,1,"x");
        f3=f1*f2;
        expected2=[-0.035355339059327376,0.,0.035355339059327376, 0.2592724864350674,0.,-0.2592724864350674, 0.37712361663282529,0.,-0.37712361663282529, -0.035355339059327376,0.,0.035355339059327376, 0.31819805153394637,0.,-0.31819805153394637]
        val=f3.getArray().getValues();
        for i in xrange(15):
            self.failUnless(abs(expected2[i]-val[i])<1.e-12);
            pass
        #
        f3=f2*f1;
        val=f3.getArray().getValues();
        for i in xrange(15):
            self.failUnless(abs(expected2[i]-val[i])<1.e-12);
            pass
        pass

    def testOperationsOnFields3(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        f1=m.fillFromAnalytic(ON_NODES,1,"x+y+z");
        f2=m.fillFromAnalytic(ON_NODES,1,"a*a+b+c*c");
        f1/=f2
        f1.checkCoherency();
        self.failUnlessEqual(f1.getTypeOfField(),ON_NODES);
        self.failUnlessEqual(f1.getTimeDiscretization(),NO_TIME);
        expected1=[-2.4999999999999991, 1.2162162162162162, 0.77868852459016391,
                   0.7407407407407407, 1.129032258064516, 0.81632653061224492,
                   0.86538461538461531, 1.0919540229885056, 0.84302325581395343]
        self.failUnlessEqual(1,f1.getNumberOfComponents());
        self.failUnlessEqual(9,f1.getNumberOfTuples());
        val=f1.getArray().getValues();
        for i in xrange(9):
            self.failUnless(abs(expected1[i]-val[i])<1.e-12);
            pass
        #
        f1=m.buildOrthogonalField();
        f2=m.fillFromAnalytic(ON_CELLS,1,"x");
        f1*=f2
        expected2=[-0.035355339059327376,0.,0.035355339059327376, 0.2592724864350674,0.,-0.2592724864350674, 0.37712361663282529,0.,-0.37712361663282529, -0.035355339059327376,0.,0.035355339059327376, 0.31819805153394637,0.,-0.31819805153394637]
        val=f1.getArray().getValues();
        for i in xrange(15):
            self.failUnless(abs(expected2[i]-val[i])<1.e-12);
            pass
        #
        f1=m.buildOrthogonalField();
        ## self.failUnlessEqual_THROW(f2*=f1,INTERP_KERNEL::Exception);
        pass

    def testOperationsOnFields4(self):
        m=MEDCouplingDataForTest.build2DTargetMesh_1();
        nbOfCells=m.getNumberOfCells();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,CONST_ON_TIME_INTERVAL);
        f1.setMesh(m);
        array=DataArrayDouble.New();
        f1.setArray(array);
        ## self.failUnlessEqual_THROW(f1.setEndArray(array),Exception);
        ## self.failUnlessEqual_THROW(f1.getEndArray(),Exception);
        arr1=[0.,10.,20.,1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.]
        arr2=[5.,15.,25.,6.,16.,26.,7.,17.,27.,8.,18.,28.,9.,19.,29.]
        array.setValues(arr1,nbOfCells,3);
        f1.setStartTime(2.,0,0);
        f1.setEndTime(3.,0,0);
        f1.checkCoherency();
        pos=[0.3,-0.2]
        res=f1.getValueOn(pos);
        self.failUnless(abs(arr1[3]-res[0])<1.e-12);
        self.failUnless(abs(arr1[4]-res[1])<1.e-12);
        self.failUnless(abs(arr1[5]-res[2])<1.e-12);
        res=None
        res=f1.getValueOn(pos,2.2);
        self.failUnless(abs(arr1[3]-res[0])<1.e-12);
        self.failUnless(abs(arr1[4]-res[1])<1.e-12);
        self.failUnless(abs(arr1[5]-res[2])<1.e-12);
        res=None
        ## self.failUnlessEqual_THROW(f1.getValueOn(pos,3.2,res),Exception);
        f2=MEDCouplingFieldDouble.New(ON_CELLS,LINEAR_TIME);
        f2.setMesh(m);
        f2.setArray(f1.getArray());
        f2.setStartTime(2.,3,0);
        f2.setEndTime(4.,13,0);
        ## self.failUnlessEqual_THROW(f2.checkCoherency(),Exception);
        array2=DataArrayDouble.New();
        array2.setValues(arr2,nbOfCells,3);
        f2.setEndArray(array2);
        f2.checkCoherency();
        #
        res=None
        res=f2.getValueOn(pos,3.21);
        self.failUnless(abs(4.025-res[0])<1.e-12);
        self.failUnless(abs(14.025-res[1])<1.e-12);
        self.failUnless(abs(24.025-res[2])<1.e-12);
        f3=f2.clone(True);
        self.failUnless(f2.isEqual(f3,1e-12,1e-12));
        f3.getEndArray().setIJ(0,0,5.001);
        self.failUnless(not f2.isEqual(f3,1e-12,1e-12));
        self.failUnless(f2.isEqual(f3,1e-12,1e-2));
        f3.setStartTime(2.1,3,0);
        self.failUnless(not f2.isEqual(f3,1e-12,1e-2));
        f3.setStartTime(2.,3,0);
        self.failUnless(f2.isEqual(f3,1e-12,1e-2));
        f3.setStartTime(2.,4,0);
        self.failUnless(not f2.isEqual(f3,1e-12,1e-2));
        f3.setStartTime(2.,3,1);
        self.failUnless(not f2.isEqual(f3,1e-12,1e-2));
        f3.setStartTime(2.,3,0);
        self.failUnless(f2.isEqual(f3,1e-12,1e-2));
        f3.setEndTime(4.1,13,0);
        self.failUnless(not f2.isEqual(f3,1e-12,1e-2));
        f3.setEndTime(4.,13,0);
        self.failUnless(f2.isEqual(f3,1e-12,1e-2));
        f3.setEndTime(4.,14,0);
        self.failUnless(not f2.isEqual(f3,1e-12,1e-2));
        f3.setEndTime(4.,13,1);
        self.failUnless(not f2.isEqual(f3,1e-12,1e-2));
        f3.setEndTime(4.,13,0);
        self.failUnless(f2.isEqual(f3,1e-12,1e-2));
        f4=f2+f2
        res=None
        res=f4.getValueOn(pos,3.21);
        self.failUnless(abs(8.05-res[0])<1.e-12);
        self.failUnless(abs(28.05-res[1])<1.e-12);
        self.failUnless(abs(48.05-res[2])<1.e-12);
        f4+=f2;
        res=None
        res=f4.getValueOn(pos,3.21);
        self.failUnless(abs(12.075-res[0])<1.e-12);
        self.failUnless(abs(42.075-res[1])<1.e-12);
        self.failUnless(abs(72.075-res[2])<1.e-12);
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
        ## self.failUnlessEqual_THROW(f1.mergeNodes(1e-10),Exception);
        pass

    def testCheckConsecutiveCellTypes(self):
        sourceMesh=MEDCouplingDataForTest.build2DSourceMesh_1();
        targetMesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        self.failUnless(sourceMesh.checkConsecutiveCellTypes());
        self.failUnless(not targetMesh.checkConsecutiveCellTypes());
        pass

    def testRearrange2ConsecutiveCellTypes(self):
        m1_1=MEDCouplingDataForTest.build2DSourceMesh_1();
        m2_1=MEDCouplingDataForTest.build2DTargetMesh_1();
        arr1=m1_1.rearrange2ConsecutiveCellTypes();
        m1_2=MEDCouplingDataForTest.build2DSourceMesh_1();
        self.failUnless(m1_2.isEqual(m1_1,1e-12));
        expected1=[0,1]
        self.failUnlessEqual(2,arr1.getNumberOfTuples());
        self.failUnlessEqual(1,arr1.getNumberOfComponents());
        self.failUnless(expected1,arr1.getValues());
        expected2=[0,3,4,1,2]
        arr1=m2_1.rearrange2ConsecutiveCellTypes();
        self.failUnlessEqual(5,arr1.getNumberOfTuples());
        self.failUnlessEqual(1,arr1.getNumberOfComponents());
        self.failUnlessEqual(expected2,arr1.getValues());
        m2_2=MEDCouplingDataForTest.build2DTargetMesh_1();
        self.failUnlessEqual(5,arr1.getNumberOfTuples());
        self.failUnlessEqual(1,arr1.getNumberOfComponents());
        self.failUnlessEqual(expected2,arr1.getValues());
        self.failUnless(not m2_2.isEqual(m2_1,1e-12));
        m2_2.renumberCells(expected2,False);
        self.failUnless(m2_2.isEqual(m2_1,1e-12));
        pass

    def testSplitByType(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        v=m1.splitByType();
        self.failUnlessEqual(3,len(v));
        m2=MEDCouplingUMesh.mergeUMeshesOnSameCoords(v);
        m2.setName(m1.getName());
        self.failUnless(m1.isEqual(m2,1.e-12));
        pass

    def testFuseUMeshesOnSameCoords(self):
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        cells1=[2,3,4]
        m3=m2.buildPartOfMySelf(cells1,True);
        self.failUnless(isinstance(m3,MEDCouplingUMesh))
        cells2=[1,2,4]
        m4=m2.buildPartOfMySelf(cells2,True);
        self.failUnless(isinstance(m4,MEDCouplingUMesh))
        cells3=[1,2]
        m5=m2.buildPartOfMySelf(cells3,True);
        self.failUnless(isinstance(m5,MEDCouplingUMesh))
        meshes=[m3,m4,m5]
        #
        m7,corr=MEDCouplingUMesh.fuseUMeshesOnSameCoords(meshes,0);
        self.failUnlessEqual(4,m7.getNumberOfCells());
        self.failUnlessEqual(3,len(corr));
        expectedVals1=[3,3,2]
        expectedVals2=[[0,1,2],[3,0,2],[3,0]]
        for i in xrange(3):
            arr=corr[i];
            self.failUnlessEqual(1,arr.getNumberOfComponents());
            nbOfVals=expectedVals1[i];
            self.failUnlessEqual(nbOfVals,arr.getNumberOfTuples());
            vals=arr.getValues();
            self.failUnlessEqual(expectedVals2[i],vals);
            pass
        arr2,fidsOfGroups=DataArrayInt.makePartition(corr,m7.getNumberOfCells());
        fidExp=[5,1,3,4]
        fidsGrp=[[1,3,5],[3,4,5],[4,5]]
        self.failUnlessEqual(3,len(fidsOfGroups));
        self.failUnlessEqual(1,arr2.getNumberOfComponents());
        self.failUnlessEqual(4,arr2.getNumberOfTuples());
        self.failUnlessEqual(fidExp,arr2.getValues());
        for i in xrange(3):
            nbOfVals=expectedVals1[i];
            self.failUnlessEqual(fidsOfGroups[i],fidsGrp[i]);
            pass
        pass

    def testFuseUMeshesOnSameCoords2(self):
        m1,m2=MEDCouplingDataForTest.build3DExtrudedUMesh_1();
        part1=[2,3,6,4,10]
        m3=m1.buildPartOfMySelf(part1,True);
        part2=[5,6,4,7]
        m4=m1.buildPartOfMySelf(part2,True);
        meshes=[m1,m3,m3,m4]
        m5,corr=MEDCouplingUMesh.fuseUMeshesOnSameCoords(meshes,0);
        self.failUnlessEqual(18,m5.getNumberOfCells());
        exp2=[
            [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],
            [2,3,6,4,10],
            [2,3,6,4,10],
            [5,6,4,7]]
        i=0;
        for it in corr:
            self.failUnlessEqual(exp2[i],it.getValues());
            i+=1
            pass
        pass

    def testBuildOrthogonalField(self):
        targetMesh=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        field=targetMesh.buildOrthogonalField();
        expected=[0.70710678118654746,0.,-0.70710678118654746]
        self.failUnlessEqual(5,field.getNumberOfTuples());
        self.failUnlessEqual(3,field.getNumberOfComponents());
        vals=field.getArray().getValues();
        for i in xrange(15):
            self.failUnless(abs(expected[i%3]-vals[i])<1e-12);
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
        self.failUnlessEqual(1,field.getNumberOfTuples());
        self.failUnlessEqual(3,field.getNumberOfComponents());
        vals=field.getArray().getValues();
        self.failUnless(abs(-0.70710678118654746-vals[0])<1e-12);
        self.failUnless(abs(0.-vals[1])<1e-12);
        self.failUnless(abs(0.70710678118654746-vals[2])<1e-12);
        pass

    def testGetCellsContainingPoint(self):
        targetMesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        pos=[0.,0.,0.4,0.4,0.,0.4,0.1,0.1,0.25,0.,0.65,0.]
        #2D basic
        t1,t2=targetMesh.getCellsContainingPoints(pos,6,1e-12);
        self.failUnlessEqual(6,len(t1));
        self.failUnlessEqual(7,len(t2));
        expectedValues1=[0,4,3,0,1,2]
        expectedValues2=[0,1,2,3,4,5,6]
        self.failUnlessEqual(t1,expectedValues1);
        self.failUnlessEqual(t2,expectedValues2);
        #2D with no help of bounding box.
        center=[0.2,0.2]
        MEDCouplingPointSet.rotate2DAlg(center,0.78539816339744830962,6,pos);
        targetMesh.rotate(center,[],0.78539816339744830962);
        t1=None
        t2=None
        t1,t2=targetMesh.getCellsContainingPoints(pos,6,1e-12);
        self.failUnlessEqual(6,len(t1));
        self.failUnlessEqual(7,len(t2));
        self.failUnlessEqual(t1,expectedValues1);
        self.failUnlessEqual(t2,expectedValues2);
        #2D outside
        pos1bis=[-0.3303300858899107,-0.11819805153394641]
        self.failUnlessEqual(-1,targetMesh.getCellContainingPoint(pos1bis,1e-12));
        #test limits 2D
        targetMesh=MEDCouplingDataForTest.build2DTargetMesh_1();
        pos2=[0.2,-0.05]
        t1=None
        t1=targetMesh.getCellsContainingPoint(pos2,1e-12)
        self.failUnlessEqual(2,len(t1));
        expectedValues3=[0,1]
        self.failUnlessEqual(t1,expectedValues3);
        pos3=[0.2,0.2]
        t1=None
        t1=targetMesh.getCellsContainingPoint(pos3,1e-12);
        self.failUnlessEqual(5,len(t1));
        expectedValues4=[0,1,2,3,4]
        self.failUnlessEqual(t1,expectedValues4);
        self.failUnlessEqual(0,targetMesh.getCellContainingPoint(pos3,1e-12));
        #3D
        targetMesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        pos4=[25.,25.,25.]
        self.failUnlessEqual(0,targetMesh.getCellContainingPoint(pos4,1e-12));
        pos5=[50.,50.,50.]
        t1=None
        t1=targetMesh.getCellsContainingPoint(pos5,1e-12);
        self.failUnlessEqual(8,len(t1));
        expectedValues5=[0,1,2,3,4,5,6,7]
        self.failUnlessEqual(t1,expectedValues5);
        pos6=[0., 50., 0.]
        t1=None
        t1=targetMesh.getCellsContainingPoint(pos6,1e-12);
        self.failUnlessEqual(2,len(t1));
        expectedValues6=[0,2]
        self.failUnlessEqual(t1,expectedValues6);
        #3D outside
        pos7=[-1.0,-1.0,0.]
        self.failUnlessEqual(-1,targetMesh.getCellContainingPoint(pos7,1e-12));
        #3D outside 2
        center2=[0.,0.,0.]
        vec2=[0.,-1.,0.]
        targetMesh.rotate(center2,vec2,0.78539816339744830962);
        pos8=[-25.,25.,12.]
        self.failUnlessEqual(-1,targetMesh.getCellContainingPoint(pos8,1e-12));
        pass

    def testGetValueOn1(self):
        # not implemented yet
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
        self.failUnlessEqual(2,len(res))
        self.failUnless(abs(8.-res[0])<1e-12);
        self.failUnless(abs(18.-res[1])<1e-12);
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
        self.failUnlessEqual(2,len(res))
        self.failUnless(abs(17.5-res[0])<1e-12);
        self.failUnless(abs(27.5-res[1])<1e-12);
        pos3=[0.033333333333333326,0.36666666666666664]
        res=None
        res=fieldOnNodes.getValueOn(pos3);
        self.failUnlessEqual(2,len(res))
        self.failUnless(abs(18.666666666666667-res[0])<1e-12);
        self.failUnless(abs(28.666666666666667-res[1])<1e-12);
        pass

    def testCMesh0(self):
        # not implemented yet
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
        self.failUnlessEqual(18,len(val))
        for i in xrange(18):
            self.failUnless(abs(expected1[i]-val[i])<1e-12);
            pass
        pass

    def testTryToShareSameCoords(self):
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        #self.failUnlessEqual(m1.getCoords()!=m2.getCoords());
        m1.tryToShareSameCoords(m2,1e-12);
        #self.failUnlessEqual(m1.getCoords()==m2.getCoords());
        m1.tryToShareSameCoords(m2,1e-12);
        #self.failUnlessEqual(m1.getCoords()==m2.getCoords());
        m2.tryToShareSameCoords(m1,1e-12);
        #self.failUnlessEqual(m1.getCoords()==m2.getCoords());
        #
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_2();
        #self.failUnlessEqual(m1.getCoords()!=m2.getCoords());
        m1.tryToShareSameCoords(m2,1e-12);
        #self.failUnlessEqual(m1.getCoords()==m2.getCoords());
        m1.tryToShareSameCoords(m2,1e-12);
        #self.failUnlessEqual(m1.getCoords()==m2.getCoords());
        m2.tryToShareSameCoords(m1,1e-12);
        #self.failUnlessEqual(m1.getCoords()==m2.getCoords());
        #
        m1=MEDCouplingDataForTest.build2DTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DSourceMesh_1();
        #self.failUnlessEqual(m1.getCoords()!=m2.getCoords());
        ## self.failUnlessEqual_THROW(m1.tryToShareSameCoords(m2,1e-12),Exception);
        pass

    def testFindNodeOnPlane(self):
        mesh=MEDCouplingDataForTest.build3DTargetMesh_1();
        pt=[300.,300.,0.]
        v=[0.,0.,2.]
        n=mesh.findNodesOnPlane(pt,v,1e-12);
        self.failUnlessEqual(9,len(n));
        m3dSurf=mesh.buildFacePartOfMySelfNode(n,True);
        self.failUnless(isinstance(m3dSurf,MEDCouplingUMesh))
        me=MEDCouplingExtrudedMesh.New(mesh,m3dSurf,0);
        da=me.getMesh3DIds();
        self.failUnlessEqual(8,me.getNumberOfCells());
        expected=[0,1,2,3,4,5,6,7]
        val=da.getValues();
        self.failUnlessEqual(expected,val);
        pass

    def testRenumberCells(self):
        m=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        self.failUnless(m.isEqual(m2,0));
        arr=[12,3,25,2,26]
        m.renumberCells(arr,True);
        self.failUnless(not m.isEqual(m2,0));
        self.failUnlessEqual(NORM_QUAD4,m.getTypeOfCell(0));
        self.failUnlessEqual(NORM_TRI3,m.getTypeOfCell(1));
        self.failUnlessEqual(NORM_QUAD4,m.getTypeOfCell(2));
        self.failUnlessEqual(NORM_TRI3,m.getTypeOfCell(3));
        self.failUnlessEqual(NORM_QUAD4,m.getTypeOfCell(4));
        arr2=[5,-1,-5,4,8]
        m.renumberCells(arr2,True);
        self.failUnless(m.isEqual(m2,0));
        pass

    def testChangeSpaceDimension(self):
        m1=MEDCouplingDataForTest.build3DSurfTargetMesh_1();
        m2=MEDCouplingDataForTest.build2DTargetMesh_1();
        #
        self.failUnlessEqual(3,m1.getSpaceDimension());
        m1.changeSpaceDimension(2);
        self.failUnlessEqual(2,m1.getSpaceDimension());
        m1.setName(m2.getName());
        self.failUnless(m1.isEqual(m2,1e-12));
        m1.changeSpaceDimension(3);
        self.failUnlessEqual(3,m1.getSpaceDimension());
        expected=[-0.3,-0.3,0., 0.2,-0.3,0., 0.7,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.7,0.2,0., -0.3,0.7,0., 0.2,0.7,0., 0.7,0.7,0.]
        val=m1.getCoords().getValues();
        for i in xrange(27):
            self.failUnless(abs(expected[i]-val[i])<1e-14);
            pass
        pass

    def setUp(self):
        pass
    pass

unittest.main()
