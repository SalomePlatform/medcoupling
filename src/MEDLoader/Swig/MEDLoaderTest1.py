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
# Author : Anthony Geay (CEA/DEN)

import MEDLoader
import unittest
from math import pi,e,sqrt
from MEDLoaderDataForTest import MEDLoaderDataForTest

class MEDLoaderTest1(unittest.TestCase):
    def testMesh1DRW(self):
        mesh=MEDLoaderDataForTest.build1DMesh_1();
        mesh.checkConsistencyLight();
        MEDLoader.WriteUMesh("Pyfile1.med",mesh,True);
        mesh_rw=MEDLoader.ReadUMeshFromFile("Pyfile1.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMesh2DCurveRW(self):
        mesh=MEDLoaderDataForTest.build2DCurveMesh_1();
        mesh.checkConsistencyLight();
        MEDLoader.WriteUMesh("Pyfile2.med",mesh,True);
        mesh_rw=MEDLoader.ReadUMeshFromFile("Pyfile2.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMesh2DRW(self):
        mesh=MEDLoaderDataForTest.build2DMesh_1();
        mesh.checkConsistencyLight();
        MEDLoader.WriteUMesh("Pyfile3.med",mesh,True);
        mesh_rw=MEDLoader.ReadUMeshFromFile("Pyfile3.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMesh3DSurfRW(self):
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        mesh.checkConsistencyLight();
        MEDLoader.WriteUMesh("Pyfile4.med",mesh,True);
        mesh_rw=MEDLoader.ReadUMeshFromFile("Pyfile4.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMesh3DRW(self):
        mesh=MEDLoaderDataForTest.build3DMesh_1();
        mesh.checkConsistencyLight();
        MEDLoader.WriteUMesh("Pyfile5.med",mesh,True);
        mesh_rw=MEDLoader.ReadUMeshFromFile("Pyfile5.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testFieldRW1(self):
        f1=MEDLoaderDataForTest.buildVecFieldOnCells_1();
        MEDLoader.WriteField("Pyfile6.med",f1,True);
        f2=MEDLoader.ReadFieldCell("Pyfile6.med",f1.getMesh().getName(),0,f1.getName(),0,1);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        #
        f1=MEDLoaderDataForTest.buildVecFieldOnNodes_1();
        MEDLoader.WriteField("Pyfile7.med",f1,True);
        f2=MEDLoader.ReadFieldNode("Pyfile7.med",f1.getMesh().getName(),0,f1.getName(),2,3);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        self.assertRaises(Exception,MEDLoader.ReadFieldCell,"Pyfile7.med",f1.getMesh().getName(),0,f1.getName(),2,3);
        pass

    def testFieldRW2(self):
        fileName="Pyfile8.med";
        VAL1=12345.67890314;
        VAL2=-1111111111111.;
        f1=MEDLoaderDataForTest.buildVecFieldOnCells_1();
        MEDLoader.WriteField(fileName,f1,True);
        f1.setTime(10.,8,9);
        f1.getArray().setIJ(0,0,VAL1);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.setTime(10.14,18,19);
        f1.getArray().setIJ(0,0,VAL2);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        #retrieving time steps...
        f2=MEDLoader.ReadFieldCell(fileName,f1.getMesh().getName(),0,f1.getName(),8,9);
        f1.setTime(10.,8,9);
        f1.getArray().setIJ(0,0,VAL1);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        f2=MEDLoader.ReadFieldCell(fileName,f1.getMesh().getName(),0,f1.getName(),0,1);
        f3=MEDLoaderDataForTest.buildVecFieldOnCells_1();
        self.assertTrue(f3.isEqual(f2,1e-12,1e-12));
        f2=MEDLoader.ReadFieldCell(fileName,f1.getMesh().getName(),0,f1.getName(),18,19);
        f1.setTime(10.14,18,19);
        f1.getArray().setIJ(0,0,VAL2);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        #test of throw on invalid (dt,it)
        self.assertRaises(Exception,MEDLoader.ReadFieldCell,fileName,f1.getMesh().getName(),0,f1.getName(),28,19);
        #ON NODES
        f1=MEDLoaderDataForTest.buildVecFieldOnNodes_1();
        fileName2="Pyfile9.med";
        MEDLoader.WriteField(fileName2,f1,True);
        f1.setTime(110.,108,109);
        tmp=f1.getArray().getPointer();
        f1.getArray().setIJ(0,3,VAL1);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName2,f1);
        f1.setTime(210.,208,209);
        f1.getArray().setIJ(0,3,VAL2);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName2,f1);
        f2=MEDLoader.ReadFieldNode(fileName2,f1.getMesh().getName(),0,f1.getName(),108,109);
        f1.setTime(110.,108,109);
        f1.getArray().setIJ(0,3,VAL1);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        f2=MEDLoader.ReadFieldNode(fileName2,f1.getMesh().getName(),0,f1.getName(),2,3);
        f3=MEDLoaderDataForTest.buildVecFieldOnNodes_1();
        self.assertTrue(f3.isEqual(f2,1e-12,1e-12));
        f2=MEDLoader.ReadFieldNode(fileName2,f1.getMesh().getName(),0,f1.getName(),208,209);
        f1.setTime(210.,208,209);
        f1.getArray().setIJ(0,3,VAL2);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        pass

    #
    # Multi field in a same file, but this field has several
    #
    def testFieldRW3(self):
        fileName="Pyfile11.med";
        VAL1=12345.67890314;
        VAL2=-1111111111111.;
        name1="AField";
        name3="AMesh1";
        f1=MEDLoaderDataForTest.buildVecFieldOnCells_1();
        f1.getMesh().setName(name3);
        f1.setName(name1);
        f1.setTime(10.,8,9);
        tmp=f1.getArray().getPointer();
        f1.getArray().setIJ(0,0,VAL1);
        MEDLoader.WriteField(fileName,f1,True);
        f1.setTime(10.14,18,19);
        f1.getArray().setIJ(0,0,VAL2);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.getMesh().setName(name3);
        f1.setTime(10.55,28,29);
        f1.getArray().setIJ(0,0,3*VAL1);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        vec=MEDLoader.GetMeshNamesOnField(fileName,name1);
        f1.setTime(10.66,38,39);
        f1.getArray().setIJ(0,0,3*VAL2);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.setTime(10.77,48,49);
        f1.getArray().setIJ(0,0,4*VAL2);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        #ON NODES
        f1=MEDLoaderDataForTest.buildVecFieldOnNodes_1();
        f1.setName(name1);
        f1.getMesh().setName(name3);
        f1.setTime(110.,8,9);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.setTime(110.,108,109);
        tmp=f1.getArray().getPointer();
        f1.getArray().setIJ(0,3,VAL1);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.setTime(210.,208,209);
        f1.getArray().setIJ(0,3,VAL2);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        #
        it1=MEDLoader.GetCellFieldIterations(fileName,name3,name1);
        self.assertEqual(5,len(it1));
        self.assertEqual(8,it1[0][0]); self.assertEqual(9,it1[0][1]);
        self.assertEqual(18,it1[1][0]); self.assertEqual(19,it1[1][1]);
        self.assertEqual(28,it1[2][0]); self.assertEqual(29,it1[2][1]);
        self.assertEqual(38,it1[3][0]); self.assertEqual(39,it1[3][1]);
        self.assertEqual(48,it1[4][0]); self.assertEqual(49,it1[4][1]);
        it3=MEDLoader.GetNodeFieldIterations(fileName,name3,name1);
        self.assertEqual(3,len(it3));
        self.assertEqual(8,it3[0][0]); self.assertEqual(9,it3[0][1]);
        self.assertEqual(108,it3[1][0]); self.assertEqual(109,it3[1][1]);
        self.assertEqual(208,it3[2][0]); self.assertEqual(209,it3[2][1]);
        #
        #
        f1=MEDLoader.ReadFieldCell(fileName,name3,0,name1,8,9);
        self.assertAlmostEqual(VAL1,f1.getArray().getIJ(0,0),13);
        f1=MEDLoader.ReadFieldCell(fileName,name3,0,name1,18,19);
        self.assertAlmostEqual(VAL2,f1.getArray().getIJ(0,0),13);
        f1=MEDLoader.ReadFieldCell(fileName,name3,0,name1,28,29);
        self.assertAlmostEqual(3*VAL1,f1.getArray().getIJ(0,0),13);
        f1=MEDLoader.ReadFieldCell(fileName,name3,0,name1,38,39);
        self.assertAlmostEqual(3*VAL2,f1.getArray().getIJ(0,0),13);
        f1=MEDLoader.ReadFieldCell(fileName,name3,0,name1,48,49);
        self.assertAlmostEqual(4*VAL2,f1.getArray().getIJ(0,0),13);
        #
        f1=MEDLoader.ReadFieldNode(fileName,name3,0,name1,8,9);
        self.assertAlmostEqual(71.,f1.getArray().getIJ(0,3),13);
        f1=MEDLoader.ReadFieldNode(fileName,name3,0,name1,108,109);
        self.assertAlmostEqual(VAL1,f1.getArray().getIJ(0,3),13);
        f1=MEDLoader.ReadFieldNode(fileName,name3,0,name1,208,209);
        self.assertAlmostEqual(VAL2,f1.getArray().getIJ(0,3),13);
        pass

    def testMultiMeshRW1(self):
        fileName="Pyfile10.med";
        mesh1=MEDLoaderDataForTest.build3DMesh_1();
        part1=[1,2,4,13,15]
        mesh2=mesh1.buildPartOfMySelf(part1,True);
        mesh2.setName("mesh2");
        part2=[3,4,13,14]
        mesh3=mesh1.buildPartOfMySelf(part2,True);
        mesh3.setName("mesh3");
        mesh4=MEDLoader.MEDCouplingUMesh.New();
        mesh4.setName("mesh4");
        mesh4.setMeshDimension(3);
        mesh4.allocateCells(1);
        conn=[0,11,1,3]
        mesh4.insertNextCell(MEDLoader.NORM_TETRA4,4,conn[0:4])
        mesh4.finishInsertingCells();
        mesh4.setCoords(mesh1.getCoords());
        meshes=[mesh1,mesh2,mesh3,mesh4]
        mnane="3DToto";
        MEDLoader.WriteUMeshesPartition(fileName,mnane,meshes,True);
        #
        mesh5=MEDLoader.ReadUMeshFromFile(fileName,mnane);
        mesh1.setName(mnane);
        part3=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
        mesh6=mesh5.buildPartOfMySelf(part3,True);
        mesh6.setName(mnane);
        self.assertTrue(mesh6.isEqual(mesh1,1e-12));
        grps=MEDLoader.GetMeshGroupsNames(fileName,mnane);
        self.assertEqual(4,len(grps));
        grps.index("mesh2");
        grps.index("mesh3");
        grps.index("mesh4");
        grps.index("3DMesh_1");
        #
        vec=("mesh2",);
        mesh2_2=MEDLoader.ReadUMeshFromGroups(fileName,mnane,0,vec);
        self.assertTrue(mesh2_2.isEqual(mesh2,1e-12));
        vec=["mesh3"];
        mesh3_2=MEDLoader.ReadUMeshFromGroups(fileName,mnane,0,vec);
        self.assertTrue(mesh3_2.isEqual(mesh3,1e-12));
        vec=["mesh4"];
        mesh4_2=MEDLoader.ReadUMeshFromGroups(fileName,mnane,0,vec);
        self.assertTrue(mesh4_2.isEqual(mesh4,1e-12));
        vec="3DMesh_1";
        mesh1_2=MEDLoader.ReadUMeshFromGroups(fileName,mnane,0,vec);
        mesh1.setName("3DMesh_1");
        self.assertTrue(mesh1_2.isEqual(mesh1,1e-12));
        #
        vec=["Family_-3","Family_-5"];
        mesh2_2=MEDLoader.ReadUMeshFromFamilies(fileName,mnane,0,vec);
        mesh2_2.setName("mesh2");
        self.assertTrue(mesh2_2.isEqual(mesh2,1e-12));
        #
        ret=MEDLoader.GetMeshFamiliesNamesOnGroup(fileName,"3DToto","3DMesh_1");
        self.assertEqual(4,len(ret));
        ref=['Family_-3','Family_-4','Family_-2','Family_-5']
        self.assertIn(ref[0],ret);
        self.assertIn(ref[1],ret);
        self.assertIn(ref[2],ret);
        self.assertIn(ref[3],ret);
        #
        ret1=MEDLoader.GetMeshGroupsNamesOnFamily(fileName,"3DToto","Family_-3");
        self.assertEqual(2,len(ret1));
        self.assertEqual(ret1[0],"3DMesh_1");
        self.assertEqual(ret1[1],"mesh2");
        pass

    def testFieldProfilRW1(self):
        fileName="Pyfile12.med";
        mesh1=MEDLoaderDataForTest.build3DMesh_1();
        da,b,newNbOfNodes=mesh1.mergeNodes(1e-12);
        MEDLoader.WriteUMesh(fileName,mesh1,True);
        part1=[1,2,4,13,15]
        mesh2=mesh1.buildPartOfMySelf(part1,True);
        mesh2.setName(mesh1.getName());#<- important for the test
        #
        nbOfCells=mesh2.getNumberOfCells();
        self.assertEqual(5,nbOfCells);
        f1=MEDLoader.MEDCouplingFieldDouble.New(MEDLoader.ON_CELLS,MEDLoader.ONE_TIME);
        f1.setName("VectorFieldOnCells");
        f1.setMesh(mesh2);
        array=MEDLoader.DataArrayDouble.New();
        array.alloc(nbOfCells,2);
        f1.setArray(array);
        arr1=[71.,171.,10.,110.,20.,120.,30.,130.,40.,140.]
        array.setValues(arr1,nbOfCells,2);
        f1.setTime(3.14,2,7);
        f1.checkConsistencyLight();
        #
        MEDLoader.WriteField(fileName,f1,False);#<- False important for the test
        #
        f2=MEDLoader.ReadFieldCell(fileName,f1.getMesh().getName(),0,f1.getName(),2,7);
        tt=MEDLoader.GetTypesOfField(fileName,f1.getMesh().getName(),f1.getName());
        self.assertEqual(tt,[MEDLoader.ON_CELLS]);
        f2.checkConsistencyLight();
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        #
        pass

    def testFieldGaussRW1(self):
        fileName="Pyfile13.med";
        f1=MEDLoaderDataForTest.buildVecFieldOnGauss_1();
        MEDLoader.WriteField(fileName,f1,True);
        f2=MEDLoader.ReadField(MEDLoader.ON_GAUSS_PT,fileName,f1.getMesh().getName(),0,f1.getName(),1,5);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        pass

    def testFieldGaussNERW1(self):
        fileName="Pyfile14.med";
        f1=MEDLoaderDataForTest.buildVecFieldOnGaussNE_1();
        MEDLoader.WriteField(fileName,f1,True);
        self.assertEqual([MEDLoader.ON_GAUSS_NE],MEDLoader.GetTypesOfField(fileName,'2DMesh_2','MyFieldOnGaussNE')) #Bug 22/5/2014
        f2=MEDLoader.ReadField(MEDLoader.ON_GAUSS_NE,fileName,f1.getMesh().getName(),0,f1.getName(),1,5);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        pass

    def testMesh3DSurfShuffleRW(self):
        fileName="Pyfile15.med";
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        renumber1=[2,5,1,0,3,4]
        mesh.renumberCells(renumber1,False);
        mesh.checkConsistencyLight();
        MEDLoader.WriteUMesh(fileName,mesh,True);
        mesh_rw=MEDLoader.ReadUMeshFromFile(fileName,mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMultiFieldShuffleRW1(self):
        fileName="Pyfile17.med";
        m=MEDLoaderDataForTest.build3DMesh_2();
        self.assertEqual(20,m.getNumberOfCells());
        self.assertEqual(45,m.getNumberOfNodes());
        polys=[1,4,6]
        m.convertToPolyTypes(polys);
        renum=[1,3,2,8,9,12,13,16,19,0,4,7,5,15,14,17,10,18,6,11]
        m.renumberCells(renum,False);
        m.orientCorrectlyPolyhedrons();
        # Writing
        MEDLoader.WriteUMesh(fileName,m,True);
        f1Tmp=m.getMeasureField(False);
        f1=f1Tmp.buildNewTimeReprFromThis(MEDLoader.ONE_TIME,False);
        f1.setTime(0.,1,2);
        f_1=f1.cloneWithMesh(True);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.applyFunc("2*x");
        f1.setTime(0.01,3,4);
        f_2=f1.cloneWithMesh(True);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.applyFunc("2*x/3");
        f1.setTime(0.02,5,6);
        f_3=f1.cloneWithMesh(True);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        # Reading
        its=[(1,2),(3,4),(5,6)];
        fs=MEDLoader.ReadFieldsOnSameMesh(MEDLoader.ON_CELLS,fileName,f_1.getMesh().getName(),0,f_1.getName(),its);
        self.assertEqual(3,len(fs));
        self.assertTrue(fs[0].isEqual(f_1,1e-12,1e-12));
        self.assertTrue(fs[1].isEqual(f_2,1e-12,1e-12));
        self.assertTrue(fs[2].isEqual(f_3,1e-12,1e-12));
        pass

    def testWriteUMeshesRW1(self):
        fileName="Pyfile18.med";
        m3d=MEDLoaderDataForTest.build3DMesh_2();
        pt=[0.,0.,-0.3]
        vec=[0.,0.,1.]
        nodes=m3d.findNodesOnPlane(pt,vec,1e-12);
        m2d=m3d.buildFacePartOfMySelfNode(nodes,True);
        renumber=[1,2,0,4,3]
        m2d.renumberCells(renumber,False);
        m2d.setName("ExampleOfMultiDimW");
        meshes=[m2d,m3d]
        MEDLoader.WriteUMeshes(fileName,meshes,True);
        m3d_bis=MEDLoader.ReadUMeshFromFile(fileName,m2d.getName(),0);
        self.assertTrue(not m3d_bis.isEqual(m3d,1e-12));
        m3d_bis.setName(m3d.getName());
        self.assertTrue(m3d_bis.isEqual(m3d,1e-12));
        m2d_bis=MEDLoader.ReadUMeshFromFile(fileName,m2d.getName(),-1);#-1 for faces
        self.assertTrue(m2d_bis.isEqual(m2d,1e-12));
        # Creation of a field on faces.
        f1=MEDLoader.MEDCouplingFieldDouble.New(MEDLoader.ON_CELLS,MEDLoader.ONE_TIME);
        f1.setName("FieldOnFacesShuffle");
        f1.setMesh(m2d);
        array=MEDLoader.DataArrayDouble.New();
        arr1=[71.,171.,10.,110.,20.,120.,30.,130.,40.,140.]
        array.setValues(arr1,m2d.getNumberOfCells(),2);
        array.setInfoOnComponent(0,"plkj [mm]");
        array.setInfoOnComponent(1,"pqqqss [mm]");
        f1.setArray(array);
        tmp=array.setValues(arr1,m2d.getNumberOfCells(),2);
        f1.setTime(3.14,2,7);
        f1.checkConsistencyLight();
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f2=MEDLoader.ReadFieldCell(fileName,f1.getMesh().getName(),-1,f1.getName(),2,7);
        self.assertTrue(f2.isEqual(f1,1e-12,1e-12));
        pass

    def testFieldNodeProfilRW1(self):
        fileName="Pyfile19.med";
        fileName2="Pyfile20.med";
        m=MEDLoaderDataForTest.build2DMesh_1();
        nbOfNodes=m.getNumberOfNodes();
        MEDLoader.WriteUMesh(fileName,m,True);
        f1=MEDLoader.MEDCouplingFieldDouble.New(MEDLoader.ON_NODES,MEDLoader.ONE_TIME);
        f1.setName("VFieldOnNodes");
        f1.setMesh(m);
        array=MEDLoader.DataArrayDouble.New();
        arr1=[1.,101.,2.,102.,3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.]
        array.setValues(arr1,nbOfNodes,2);
        f1.setArray(array);
        array.setInfoOnComponent(0,"tyty [mm]");
        array.setInfoOnComponent(1,"uiop [MW]");
        f1.setTime(3.14,2,7);
        f1.checkConsistencyLight();
        arr2=[1,4]
        f2=f1.buildSubPart(arr2);
        f2.getMesh().setName(f1.getMesh().getName());
        MEDLoader.WriteField(fileName,f2,False);#<- False important for the test
        #
        f3=MEDLoader.ReadFieldNode(fileName,f2.getMesh().getName(),0,f2.getName(),2,7);
        f3.checkConsistencyLight();
        self.assertTrue(f3.isEqual(f2,1e-12,1e-12));
        #
        arr3=[1,3,0,5,2,4]
        f2.renumberNodes(arr3);
        MEDLoader.WriteUMesh(fileName2,m,True);
        MEDLoader.WriteField(fileName2,f2,False);#<- False important for the test
        f3=MEDLoader.ReadFieldNode(fileName2,f2.getMesh().getName(),0,f2.getName(),2,7);
        f3.checkConsistencyLight();
        self.assertTrue(f3.isEqual(f2,1e-12,1e-12));
        #
        pass

    def testFieldNodeProfilRW2(self):
        fileName="Pyfile23.med";
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        MEDLoader.WriteUMesh(fileName,mesh,True);
        #
        f1=MEDLoader.MEDCouplingFieldDouble.New(MEDLoader.ON_NODES,MEDLoader.ONE_TIME);
        f1.setName("FieldMix");
        f1.setMesh(mesh);
        arr2=[1071.,1171.,1010.,1110.,1020.,1120.,1030.,1130.,1040.,1140.,1050.,1150.,
              1060.,1160.,1070.,1170.,1080.,1180.,1090.,1190.,1091.,1191.,1092.,1192.];
        array=MEDLoader.DataArrayDouble.New();
        array.setValues(arr2,12,2);
        f1.setArray(array);
        array.setInfoOnComponent(0,"plkj [mm]");
        array.setInfoOnComponent(1,"pqqqss [mm]");
        tmp=array.getPointer();
        f1.setTime(3.17,2,7);
        #
        renumArr=[3,7,2,1,5,11,10,0,9,6,8,4]
        f1.renumberNodes(renumArr);
        f1.checkConsistencyLight();
        MEDLoader.WriteField(fileName,f1,False);#<- False important for the test
        f2=MEDLoader.ReadFieldNode(fileName,f1.getMesh().getName(),0,f1.getName(),2,7);
        self.assertTrue(f2.isEqual(f1,1e-12,1e-12));
        #
        pass

    def testMixCellAndNodesFieldRW1(self):
        fileName="Pyfile21.med";
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        f1=MEDLoader.MEDCouplingFieldDouble.New(MEDLoader.ON_CELLS,MEDLoader.ONE_TIME);
        f1.setName("FieldMix");
        f1.setMesh(mesh);
        array=MEDLoader.DataArrayDouble.New();
        f1.setArray(array);
        arr1=[71.,171.,10.,110.,20.,120.,30.,130.,40.,140.,50.,150.]
        array.setValues(arr1,6,2);
        array.setInfoOnComponent(0,"plkj [mm]");
        array.setInfoOnComponent(1,"pqqqss [mm]");
        f1.setTime(3.14,2,7);
        f1.checkConsistencyLight();
        #
        f2=MEDLoader.MEDCouplingFieldDouble.New(MEDLoader.ON_NODES,MEDLoader.ONE_TIME);
        f2.setName("FieldMix");
        f2.setMesh(mesh);
        array=MEDLoader.DataArrayDouble.New();
        f2.setArray(array);
        arr2=[1071.,1171.,1010.,1110.,1020.,1120.,1030.,1130.,1040.,1140.,1050.,1150.,
              1060.,1160.,1070.,1170.,1080.,1180.,1090.,1190.,1091.,1191.,1092.,1192.]
        array.setValues(arr2,12,2)
        array.setInfoOnComponent(0,"plkj [mm]");
        array.setInfoOnComponent(1,"pqqqss [mm]");
        f2.setTime(3.14,2,7);
        f2.checkConsistencyLight();
        #
        MEDLoader.WriteField(fileName,f1,True);
        ts=MEDLoader.GetTypesOfField(fileName,f1.getMesh().getName(),f1.getName());
        self.assertEqual(1,len(ts));
        self.assertEqual(MEDLoader.ON_CELLS,ts[0]);
        fs=MEDLoader.GetAllFieldNamesOnMesh(fileName,f1.getMesh().getName());
        self.assertEqual(1,len(fs));
        self.assertTrue(fs[0]=="FieldMix");
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f2);
        fs=MEDLoader.GetAllFieldNamesOnMesh(fileName,f1.getMesh().getName());
        self.assertEqual(1,len(fs));
        self.assertTrue(fs[0]=="FieldMix");
        #
        ts=MEDLoader.GetTypesOfField(fileName,f1.getMesh().getName(),f1.getName());
        self.assertEqual(2,len(ts));
        self.assertEqual(MEDLoader.ON_NODES,ts[0]);
        self.assertEqual(MEDLoader.ON_CELLS,ts[1]);
        #
        f3=MEDLoader.ReadFieldNode(fileName,f1.getMesh().getName(),0,f1.getName(),2,7);
        self.assertTrue(f3.isEqual(f2,1e-12,1e-12));
        f3=MEDLoader.ReadFieldCell(fileName,f1.getMesh().getName(),0,f1.getName(),2,7);
        self.assertTrue(f3.isEqual(f1,1e-12,1e-12));
        #
        pass

    def testGetAllFieldNamesRW1(self):
        fileName="Pyfile22.med";
        mesh=MEDLoaderDataForTest.build2DMesh_2();
        f1=MEDLoader.MEDCouplingFieldDouble.New(MEDLoader.ON_NODES,MEDLoader.ONE_TIME);
        f1.setName("Field1");
        f1.setTime(3.44,5,6);
        f1.setMesh(mesh);
        f1.fillFromAnalytic(2,"x+y");
        MEDLoader.WriteField(fileName,f1,True);
        f1.setTime(1002.3,7,8);
        f1.fillFromAnalytic(2,"x+77.*y");
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.setName("Field2");
        MEDLoader.WriteField(fileName,f1,False);
        f1.setName("Field3");
        mesh.setName("2DMesh_2Bis");
        MEDLoader.WriteField(fileName,f1,False);
        f1=MEDLoader.MEDCouplingFieldDouble.New(MEDLoader.ON_CELLS,MEDLoader.ONE_TIME);
        f1.setName("Field8");
        f1.setTime(8.99,7,9);
        f1.setMesh(mesh);
        f1.fillFromAnalytic(3,"3*x+y");
        MEDLoader.WriteField(fileName,f1,False);
        fs=MEDLoader.GetAllFieldNames(fileName);
        self.assertEqual(4,len(fs));
        self.assertTrue(fs[0]=="Field1");
        self.assertTrue(fs[1]=="Field2");
        self.assertTrue(fs[2]=="Field3");
        self.assertTrue(fs[3]=="Field8");
        pass

    def testBigNbOfCompoNonReg(self):
        fileName="Pyfile57.med"
        m=MEDLoader.MEDCouplingCMesh() ; m.setCoords(MEDLoader.DataArrayDouble([0,1,2,3]),MEDLoader.DataArrayDouble([0,1]),MEDLoader.DataArrayDouble([0,1]))
        m=m.buildUnstructured() ; m.setName("TinyMesh")
        f=MEDLoader.MEDCouplingFieldDouble(MEDLoader.ON_CELLS) ; f.setMesh(m)
        nbOfCompo=4100
        arr=MEDLoader.DataArrayDouble(nbOfCompo*3) ; arr.iota()
        arr.rearrange(nbOfCompo)
        arr.setInfoOnComponents(["c%i" % (i) for i in range(nbOfCompo)])
        f.setArray(arr)
        f.setName("FieldBigCompo")
        MEDLoader.WriteField(fileName,f,True)
        f2=MEDLoader.ReadFieldCell(fileName,m.getName(),0,f.getName(),-1,-1)
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        pass

    def testMultiMeshTypeWrite0(self):
        fname="Pyfile73.med"
        m=MEDLoader.MEDCoupling1SGTUMesh("mesh",MEDLoader.NORM_QUAD4) ; m.allocateCells()
        m.insertNextCell([0,2,1,3])
        m.setCoords(MEDLoader.DataArrayDouble([0.,0.,1.,1.,1.,0.,0.,1.],4,2))
        #
        ms = [m.deepCopy() for i in range(4)]
        for i,elt in enumerate(ms):
            elt.translate([float(i)*1.5,0.])
            pass
        #
        m0=MEDLoader.MEDCoupling1SGTUMesh.Merge1SGTUMeshes(ms)
        f=m0.getMeasureField(False) ; f.getArray().setInfoOnComponents(["ABC [defg]"])
        MEDLoader.WriteField(fname,f,True)
        #
        fRead=MEDLoader.ReadFieldCell(fname,"merge",0,f.getName(),-1,-1)
        fRead.setMesh(MEDLoader.MEDCoupling1SGTUMesh(fRead.getMesh()))
        self.assertTrue(f.isEqual(fRead,1e-12,1e-12))
        #
        m0=m0.buildUnstructured() ; m0.convertAllToPoly()
        m0=MEDLoader.MEDCoupling1DGTUMesh(m0)
        f=m0.getMeasureField(False) ; f.getArray().setInfoOnComponents(["ABC [defg]"])
        MEDLoader.WriteField(fname,f,True)
        #
        fRead=MEDLoader.ReadFieldCell(fname,"merge",0,f.getName(),-1,-1)
        fRead.setMesh(MEDLoader.MEDCoupling1DGTUMesh(fRead.getMesh()))
        self.assertTrue(f.isEqual(fRead,1e-12,1e-12))
        #
        m0=MEDLoader.MEDCouplingCMesh()
        arr=MEDLoader.DataArrayDouble(4) ; arr.iota()
        m0.setCoords(arr,arr)
        m0.setName("mesh")
        f=m0.getMeasureField(False) ; f.getArray().setInfoOnComponents(["ABC [defg]"])
        MEDLoader.WriteField(fname,f,True)
        #
        fRead=MEDLoader.ReadFieldCell(fname,"mesh",0,f.getName(),-1,-1)
        self.assertTrue(f.isEqual(fRead,1e-12,1e-12))
        #
        c=m0.buildUnstructured().getCoords()
        m0=MEDLoader.MEDCouplingCurveLinearMesh("mesh")
        m0.setNodeGridStructure([4,4])
        m0.setCoords(c)
        f=m0.getMeasureField(False) ; f.getArray().setInfoOnComponents(["ABC [defg]"])
        MEDLoader.WriteField(fname,f,True)
        #
        fRead=MEDLoader.ReadFieldCell(fname,"mesh",0,f.getName(),-1,-1)
        self.assertTrue(f.isEqual(fRead,1e-12,1e-12))
        pass
    
    def testMultiMeshTypeWrite1(self):
        fname="Pyfile74.med"
        m=MEDLoader.MEDCoupling1SGTUMesh("mesh",MEDLoader.NORM_QUAD4) ; m.allocateCells()
        m.insertNextCell([0,2,1,3])
        m.setCoords(MEDLoader.DataArrayDouble([0.,0.,1.,1.,1.,0.,0.,1.],4,2))
        #
        ms = [m.deepCopy() for i in range(4)]
        for i,elt in enumerate(ms):
            elt.translate([float(i)*1.5,0.])
            pass
        m0=MEDLoader.MEDCoupling1SGTUMesh.Merge1SGTUMeshes(ms)
        MEDLoader.WriteMesh(fname,m0,True)
        #
        mRead=MEDLoader.ReadMeshFromFile(fname,"merge",0)
        self.assertTrue(isinstance(mRead,MEDLoader.MEDCouplingUMesh))
        mRead=MEDLoader.MEDCoupling1SGTUMesh(mRead)
        self.assertTrue(m0.isEqual(mRead,1e-12))
        #
        m0=m0.buildUnstructured() ; m0.convertAllToPoly()
        m0=MEDLoader.MEDCoupling1DGTUMesh(m0)
        MEDLoader.WriteMesh(fname,m0,True)
        #
        mRead=MEDLoader.ReadMeshFromFile(fname,"merge",0)
        mRead=MEDLoader.MEDCoupling1DGTUMesh(mRead)
        self.assertTrue(m0.isEqual(mRead,1e-12))
        #
        m0=MEDLoader.MEDCouplingCMesh()
        arr=MEDLoader.DataArrayDouble(4) ; arr.iota()
        m0.setCoords(arr,arr)
        m0.setName("mesh")
        MEDLoader.WriteMesh(fname,m0,True)
        #
        mRead=MEDLoader.ReadMeshFromFile(fname,0)
        self.assertTrue(isinstance(mRead,MEDLoader.MEDCouplingCMesh))
        self.assertTrue(m0.isEqual(mRead,1e-12))
        #
        c=m0.buildUnstructured().getCoords()
        m0=MEDLoader.MEDCouplingCurveLinearMesh("mesh")
        m0.setNodeGridStructure([4,4])
        m0.setCoords(c)
        MEDLoader.WriteMesh(fname,m0,True)
        #
        mRead=MEDLoader.ReadMeshFromFile(fname,0)
        self.assertTrue(isinstance(mRead,MEDLoader.MEDCouplingCurveLinearMesh))
        self.assertTrue(m0.isEqual(mRead,1e-12))
        pass

    def testChangeGroupName(self):
        """ This test is a non regression test on MEDFileUMesh.changeGroupName thanks to Alliance.
        """
        mfd=MEDLoaderDataForTest.buildAMEDFileDataWithGroupOnOneFamilyForSauv()
        mesh = mfd.getMeshes().getMeshAtPos(0)
        mesh.changeGroupName("grp0_LM1", "xonall1")
        self.assertTrue("xonall1" in mesh.getGroupsNames())
        pass

    def testFieldWithTooLongName(self):
        """ This test is a non regression test, to check that in basic API the policies are taken into account.
        """
        fname="Pyfile75.med"
        # Coordinates
        coords = [0.,0., 0.,1., 1.,1., 1.,0.]
        # lvl 0 connectivity
        conn2D   = [1,2,3,4]
        # lvl 0 mesh
        m=MEDLoader.MEDCouplingUMesh.New("mesh",2)
        m.allocateCells(1)
        m.insertNextCell(MEDLoader.NORM_QUAD4,4,conn2D)
        m.finishInsertingCells()
        # assigning coordinates
        meshCoords=MEDLoader.DataArrayDouble.New()
        meshCoords.setValues(coords, 4, 2)
        m.setCoords(meshCoords)
        #
        f=MEDLoader.MEDCouplingFieldDouble.New(MEDLoader.ON_CELLS,MEDLoader.ONE_TIME)
        f.setMesh(m)
        d=MEDLoader.DataArrayDouble.New()
        d.alloc(1,1)
        d.iota(1.)
        # seting a long name
        d.setInfoOnComponent(0,"CONCENTRATION of I129")
        f.setArray(d)
        f.setName("field")
        #
        mm=MEDLoader.MEDFileUMesh()
        MEDLoader.SetTooLongStrPolicy(2)
        MEDLoader.AssignStaticWritePropertiesTo(mm)
        self.assertEqual(2,mm.getTooLongStrPolicy())
        MEDLoader.SetTooLongStrPolicy(0)
        MEDLoader.AssignStaticWritePropertiesTo(mm)
        self.assertEqual(0,mm.getTooLongStrPolicy())
        del mm
        #
        MEDLoader.SetTooLongStrPolicy(2)
        self.assertRaises(MEDLoader.InterpKernelException,MEDLoader.WriteField,fname,f,True)# the component name is too long + policy 2 -> throw
        f.getArray().setInfoOnComponent(0,'I129')
        MEDLoader.WriteField(fname,f,True)
        pass

    def testUsingAlreadyWrittenMesh2(self):
        """ This test focuses on MEDLoader.WriteFieldUsingAlreadyWrittenMesh with mesh different from UMesh.
        """
        fname="Pyfile76.med"
        mesh=MEDLoader.MEDCouplingCMesh("mesh")
        arrX=MEDLoader.DataArrayDouble([0,1,2,3])
        arrY=MEDLoader.DataArrayDouble([0,2,3,5,7])
        arrZ=MEDLoader.DataArrayDouble([7])
        mesh.setCoords(arrX,arrY,arrZ)
        #
        f1=MEDLoader.MEDCouplingFieldDouble(MEDLoader.ON_NODES) ; f1.setName("f1")
        f1.setMesh(mesh)
        arr=MEDLoader.DataArrayDouble(20) ; arr.iota()
        f1.setArray(arr)
        f1.checkConsistencyLight()
        #
        f2=MEDLoader.MEDCouplingFieldDouble(MEDLoader.ON_NODES) ; f2.setName("f2")
        f2.setMesh(mesh)
        arr=MEDLoader.DataArrayDouble(20) ; arr.iota() ; arr*=3
        f2.setArray(arr)
        f2.checkConsistencyLight()
        #
        f11=f1.deepCopy() ; (f11.getArray())[:]*=4 ; f11.setTime(1.1,5,6)
        #
        MEDLoader.WriteMesh(fname,f1.getMesh(),True)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f1)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f2)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f11)
        ##
        f1r=MEDLoader.ReadFieldNode(fname,"mesh",0,"f1",-1,-1);
        self.assertTrue(f1.isEqual(f1r,1e-12,1e-12))
        self.assertTrue(f1r.getArray().isEqual(MEDLoader.DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.]),1e-12))
        f2r=MEDLoader.ReadFieldNode(fname,"mesh",0,"f2",-1,-1);
        self.assertTrue(f2.isEqual(f2r,1e-12,1e-12))
        self.assertTrue(f2r.getArray().isEqual(MEDLoader.DataArrayDouble([0.,3.,6.,9.,12.,15.,18.,21.,24.,27.,30.,33.,36.,39.,42.,45.,48.,51.,54.,57.]),1e-12))
        f3r=MEDLoader.ReadFieldNode(fname,"mesh",0,"f1",5,6);
        self.assertTrue(f11.isEqual(f3r,1e-12,1e-12))
        self.assertTrue(f3r.getArray().isEqual(MEDLoader.DataArrayDouble([0.,4.,8.,12.,16.,20.,24.,28.,32.,36.,40.,44.,48.,52.,56.,60.,64.,68.,72.,76.]),1e-12))
        pass

    def testEasyFieldRead1(self):
        fname="Pyfile111.med"
        arr=MEDLoader.DataArrayDouble(4) ; arr.iota()
        m=MEDLoader.MEDCouplingCMesh() ; m.setCoords(arr,arr)
        m=m.buildUnstructured()
        m.setName("mesh")
        f=MEDLoader.MEDCouplingFieldDouble(MEDLoader.ON_CELLS)
        f.setName("field")
        f.setTime(3.,1,2)
        da=MEDLoader.DataArrayDouble([2,3,4,5,6,7,8,9,10])
        f.setArray(da) ; f.setMesh(m)
        MEDLoader.WriteField(fname,f,True)
        #
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname,"field"),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname,"field",1,2),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field",1,2),1e-12,1e-12))
        #
        f3=f.deepCopy()
        f3.setArray(f.getArray()+30)
        f3.setName("field2")
        f3.setTime(5.,4,5)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f3)
        #
        self.assertRaises(Exception,MEDLoader.ReadField,fname) # because several fields in fname now
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname,"field"),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname,"field",1,2),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field",1,2),1e-12,1e-12))
        self.assertTrue(f3.isEqual(MEDLoader.ReadField(fname,"field2"),1e-12,1e-12))
        self.assertTrue(f3.isEqual(MEDLoader.ReadField(fname,"field2",4,5),1e-12,1e-12))
        self.assertTrue(f3.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field2",4,5),1e-12,1e-12))
        #
        f2=f.deepCopy()
        f2.setTime(4.,3,4)
        f2.setArray(f2.getArray()+10)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f2)
        #
        self.assertRaises(Exception,MEDLoader.ReadField,fname) # because unique field "field" has more than one time step
        self.assertRaises(Exception,MEDLoader.ReadField,fname,"field") # because unique field "field" has more than one time step
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname,"field",1,2),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field",1,2),1e-12,1e-12))
        self.assertTrue(f2.isEqual(MEDLoader.ReadField(fname,"field",3,4),1e-12,1e-12))
        self.assertTrue(f2.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field",3,4),1e-12,1e-12))
        #
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname,"field",1,2),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field",1,2),1e-12,1e-12))
        self.assertTrue(f2.isEqual(MEDLoader.ReadField(fname,"field",3,4),1e-12,1e-12))
        self.assertTrue(f2.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field",3,4),1e-12,1e-12))
        #
        self.assertTrue(f3.isEqual(MEDLoader.ReadField(fname,"field2"),1e-12,1e-12))
        self.assertTrue(f3.isEqual(MEDLoader.ReadField(fname,"field2",4,5),1e-12,1e-12))
        self.assertRaises(Exception,MEDLoader.ReadField,fname,"field2",5,5) # invalid time step
        self.assertTrue(f3.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field2",4,5),1e-12,1e-12))
        self.assertRaises(Exception,MEDLoader.ReadField,MEDLoader.ON_CELLS,fname,"mesh",0,"field2",5,5) # invalid time step
        # Test on profile - restart from scratch
        mm=MEDLoader.MEDFileUMesh()
        mm[0]=m
        mm.write(fname,2)
        #
        pfl = MEDLoader.DataArrayInt(list(range(8)))
        pfl.setName("PFL")
        #
        f=MEDLoader.MEDCouplingFieldDouble(MEDLoader.ON_CELLS)
        f.setName("field")
        f.setTime(3.,1,2)
        da=MEDLoader.DataArrayDouble([2,3,4,5,6,7,8,9])
        f.setArray(da) ; f.setMesh(m[pfl])
        f.checkConsistencyLight()
        #
        f1ts=MEDLoader.MEDFileField1TS()
        f1ts.setFieldProfile(f,mm,0,pfl)
        f1ts.write(fname,0)
        #
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname,"field"),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname,"field",1,2),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field",1,2),1e-12,1e-12))
        #
        f3=f.deepCopy()
        f3.setArray(f.getArray()+30)
        f3.setName("field2")
        f3.setTime(5.,4,5)
        f1ts=MEDLoader.MEDFileField1TS()
        f1ts.setFieldProfile(f3,mm,0,pfl)
        f1ts.write(fname,0)
        #
        self.assertRaises(Exception,MEDLoader.ReadField,fname) # because several fields in fname now
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname,"field"),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname,"field",1,2),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field",1,2),1e-12,1e-12))
        self.assertTrue(f3.isEqual(MEDLoader.ReadField(fname,"field2"),1e-12,1e-12))
        self.assertTrue(f3.isEqual(MEDLoader.ReadField(fname,"field2",4,5),1e-12,1e-12))
        self.assertTrue(f3.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field2",4,5),1e-12,1e-12))
        #
        f2=f.deepCopy()
        f2.setTime(4.,3,4)
        f2.setArray(f2.getArray()+10)
        f1ts=MEDLoader.MEDFileField1TS()
        f1ts.setFieldProfile(f2,mm,0,pfl)
        f1ts.write(fname,0)
        #
        self.assertRaises(Exception,MEDLoader.ReadField,fname) # because unique field "field" has more than one time step
        self.assertRaises(Exception,MEDLoader.ReadField,fname,"field") # because unique field "field" has more than one time step
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname,"field",1,2),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field",1,2),1e-12,1e-12))
        self.assertTrue(f2.isEqual(MEDLoader.ReadField(fname,"field",3,4),1e-12,1e-12))
        self.assertTrue(f2.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field",3,4),1e-12,1e-12))
        #
        self.assertTrue(f.isEqual(MEDLoader.ReadField(fname,"field",1,2),1e-12,1e-12))
        self.assertTrue(f.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field",1,2),1e-12,1e-12))
        self.assertTrue(f2.isEqual(MEDLoader.ReadField(fname,"field",3,4),1e-12,1e-12))
        self.assertTrue(f2.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field",3,4),1e-12,1e-12))
        #
        self.assertTrue(f3.isEqual(MEDLoader.ReadField(fname,"field2"),1e-12,1e-12))
        self.assertTrue(f3.isEqual(MEDLoader.ReadField(fname,"field2",4,5),1e-12,1e-12))
        self.assertRaises(Exception,MEDLoader.ReadField,fname,"field2",5,5) # invalid time step
        self.assertTrue(f3.isEqual(MEDLoader.ReadField(MEDLoader.ON_CELLS,fname,"mesh",0,"field2",4,5),1e-12,1e-12))
        self.assertRaises(Exception,MEDLoader.ReadField,MEDLoader.ON_CELLS,fname,"mesh",0,"field2",5,5) # invalid time step
        pass
    
    pass

if __name__ == "__main__":
  unittest.main()
