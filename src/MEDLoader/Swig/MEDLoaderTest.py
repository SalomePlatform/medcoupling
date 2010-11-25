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

from libMEDLoader_Swig import *
import unittest
from math import pi,e,sqrt
from MEDLoaderDataForTest import MEDLoaderDataForTest

class MEDLoaderTest(unittest.TestCase):
    def testMesh1DRW(self):
        mesh=MEDLoaderDataForTest.build1DMesh_1();
        mesh.checkCoherency();
        MEDLoader.WriteUMesh("Pyfile1.med",mesh,True);
        mesh_rw=MEDLoader.ReadUMeshFromFile("Pyfile1.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMesh2DCurveRW(self):
        mesh=MEDLoaderDataForTest.build2DCurveMesh_1();
        mesh.checkCoherency();
        MEDLoader.WriteUMesh("Pyfile2.med",mesh,True);
        mesh_rw=MEDLoader.ReadUMeshFromFile("Pyfile2.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMesh2DRW(self):
        mesh=MEDLoaderDataForTest.build2DMesh_1();
        mesh.checkCoherency();
        MEDLoader.WriteUMesh("Pyfile3.med",mesh,True);
        mesh_rw=MEDLoader.ReadUMeshFromFile("Pyfile3.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMesh3DSurfRW(self):
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        mesh.checkCoherency();
        MEDLoader.WriteUMesh("Pyfile4.med",mesh,True);
        mesh_rw=MEDLoader.ReadUMeshFromFile("Pyfile4.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMesh3DRW(self):
        mesh=MEDLoaderDataForTest.build3DMesh_1();
        mesh.checkCoherency();
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
        name2="AMesh2";
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
        f1.getMesh().setName(name2);
        f1.setTime(10.55,28,29);
        f1.getArray().setIJ(0,0,3*VAL1);
        MEDLoader.WriteField(fileName,f1,False);
        vec=MEDLoader.GetMeshNamesOnField(fileName,name1);
        self.assertEqual(2,len(vec));
        self.assertTrue(vec[0]==name3);
        self.assertTrue(vec[1]==name2);
        f1.setTime(10.66,38,39);
        f1.getArray().setIJ(0,0,3*VAL2);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.setTime(10.77,48,49);
        f1.getArray().setIJ(0,0,4*VAL2);
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        #ON NODES
        f1=MEDLoaderDataForTest.buildVecFieldOnNodes_1();
        f1.setName(name1);
        f1.getMesh().setName(name2);
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
        self.assertEqual(2,len(it1));
        self.assertEqual(8,it1[0][0]); self.assertEqual(9,it1[0][1]);
        self.assertEqual(18,it1[1][0]); self.assertEqual(19,it1[1][1]);
        it2=MEDLoader.GetCellFieldIterations(fileName,name2,name1);
        self.assertEqual(3,len(it2));
        self.assertEqual(28,it2[0][0]); self.assertEqual(29,it2[0][1]);
        self.assertEqual(38,it2[1][0]); self.assertEqual(39,it2[1][1]);
        self.assertEqual(48,it2[2][0]); self.assertEqual(49,it2[2][1]);
        it3=MEDLoader.GetNodeFieldIterations(fileName,name2,name1);
        self.assertEqual(3,len(it3));
        self.assertEqual(8,it3[0][0]); self.assertEqual(9,it3[0][1]);
        self.assertEqual(108,it3[1][0]); self.assertEqual(109,it3[1][1]);
        self.assertEqual(208,it3[2][0]); self.assertEqual(209,it3[2][1]);
        it4=MEDLoader.GetNodeFieldIterations(fileName,name3,name1);
        self.assertTrue(len(it4)==0);
        #
        #
        f1=MEDLoader.ReadFieldCell(fileName,name3,0,name1,8,9);
        self.assertAlmostEqual(VAL1,f1.getArray().getIJ(0,0),13);
        f1=MEDLoader.ReadFieldCell(fileName,name3,0,name1,18,19);
        self.assertAlmostEqual(VAL2,f1.getArray().getIJ(0,0),13);
        f1=MEDLoader.ReadFieldCell(fileName,name2,0,name1,28,29);
        self.assertAlmostEqual(3*VAL1,f1.getArray().getIJ(0,0),13);
        f1=MEDLoader.ReadFieldCell(fileName,name2,0,name1,38,39);
        self.assertAlmostEqual(3*VAL2,f1.getArray().getIJ(0,0),13);
        f1=MEDLoader.ReadFieldCell(fileName,name2,0,name1,48,49);
        self.assertAlmostEqual(4*VAL2,f1.getArray().getIJ(0,0),13);
        #
        f1=MEDLoader.ReadFieldNode(fileName,name2,0,name1,8,9);
        self.assertAlmostEqual(71.,f1.getArray().getIJ(0,3),13);
        f1=MEDLoader.ReadFieldNode(fileName,name2,0,name1,108,109);
        self.assertAlmostEqual(VAL1,f1.getArray().getIJ(0,3),13);
        f1=MEDLoader.ReadFieldNode(fileName,name2,0,name1,208,209);
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
        mesh4=MEDCouplingUMesh.New();
        mesh4.setName("mesh4");
        mesh4.setMeshDimension(3);
        mesh4.allocateCells(1);
        conn=[0,11,1,3]
        mesh4.insertNextCell(NORM_TETRA4,4,conn[0:4])
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
        vec=["mesh2"];
        mesh2_2=MEDLoader.ReadUMeshFromGroups(fileName,mnane,0,vec);
        self.assertTrue(mesh2_2.isEqual(mesh2,1e-12));
        vec=["mesh3"];
        mesh3_2=MEDLoader.ReadUMeshFromGroups(fileName,mnane,0,vec);
        self.assertTrue(mesh3_2.isEqual(mesh3,1e-12));
        vec=["mesh4"];
        mesh4_2=MEDLoader.ReadUMeshFromGroups(fileName,mnane,0,vec);
        self.assertTrue(mesh4_2.isEqual(mesh4,1e-12));
        vec=["3DMesh_1"];
        mesh1_2=MEDLoader.ReadUMeshFromGroups(fileName,mnane,0,vec);
        mesh1.setName("3DMesh_1");
        self.assertTrue(mesh1_2.isEqual(mesh1,1e-12));
        #
        vec=["Family_4","Family_2"];
        mesh2_2=MEDLoader.ReadUMeshFromFamilies(fileName,mnane,0,vec);
        mesh2_2.setName("mesh2");
        self.assertTrue(mesh2_2.isEqual(mesh2,1e-12));
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
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setName("VectorFieldOnCells");
        f1.setMesh(mesh2);
        array=DataArrayDouble.New();
        array.alloc(nbOfCells,2);
        f1.setArray(array);
        arr1=[71.,171.,10.,110.,20.,120.,30.,130.,40.,140.]
        array.setValues(arr1,nbOfCells,2);
        f1.setTime(3.14,2,7);
        f1.checkCoherency();
        #
        MEDLoader.WriteField(fileName,f1,False);#<- False important for the test
        #
        f2=MEDLoader.ReadFieldCell(fileName,f1.getMesh().getName(),0,f1.getName(),2,7);
        f2.checkCoherency();
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        #
        pass

    def testFieldGaussRW1(self):
        fileName="Pyfile13.med";
        f1=MEDLoaderDataForTest.buildVecFieldOnGauss_1();
        MEDLoader.WriteField(fileName,f1,True);
        f2=MEDLoader.ReadField(ON_GAUSS_PT,fileName,f1.getMesh().getName(),0,f1.getName(),1,5);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        pass

    def testFieldGaussNERW1(self):
        fileName="Pyfile14.med";
        f1=MEDLoaderDataForTest.buildVecFieldOnGaussNE_1();
        MEDLoader.WriteField(fileName,f1,True);
        f2=MEDLoader.ReadField(ON_GAUSS_NE,fileName,f1.getMesh().getName(),0,f1.getName(),1,5);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        pass

    def testMesh3DSurfShuffleRW(self):
        fileName="Pyfile15.med";
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        renumber1=[2,5,1,0,3,4]
        mesh.renumberCells(renumber1,False);
        mesh.checkCoherency();
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
        f1=f1Tmp.buildNewTimeReprFromThis(ONE_TIME,False);
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
        fs=MEDLoader.ReadFieldsOnSameMesh(ON_CELLS,fileName,f_1.getMesh().getName(),0,f_1.getName(),its);
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
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setName("FieldOnFacesShuffle");
        f1.setMesh(m2d);
        array=DataArrayDouble.New();
        arr1=[71.,171.,10.,110.,20.,120.,30.,130.,40.,140.]
        array.setValues(arr1,m2d.getNumberOfCells(),2);
        array.setInfoOnComponent(0,"plkj (mm)");
        array.setInfoOnComponent(1,"pqqqss (mm)");
        f1.setArray(array);
        tmp=array.setValues(arr1,m2d.getNumberOfCells(),2);
        f1.setTime(3.14,2,7);
        f1.checkCoherency();
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
        f1=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
        f1.setName("VFieldOnNodes");
        f1.setMesh(m);
        array=DataArrayDouble.New();
        arr1=[1.,101.,2.,102.,3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.]
        array.setValues(arr1,nbOfNodes,2);
        f1.setArray(array);
        array.setInfoOnComponent(0,"tyty (mm)");
        array.setInfoOnComponent(1,"uiop (MW)");
        f1.setTime(3.14,2,7);
        f1.checkCoherency();
        arr2=[2,4,5,3,6,7]
        f2=f1.buildSubPart(arr2);
        f2.getMesh().setName(f1.getMesh().getName());
        MEDLoader.WriteField(fileName,f2,False);#<- False important for the test
        #
        f3=MEDLoader.ReadFieldNode(fileName,f2.getMesh().getName(),0,f2.getName(),2,7);
        f3.checkCoherency();
        self.assertTrue(f3.isEqual(f2,1e-12,1e-12));
        #
        arr3=[1,3,0,5,2,4]
        f2.renumberNodes(arr3);
        MEDLoader.WriteUMesh(fileName2,m,True);
        MEDLoader.WriteField(fileName2,f2,False);#<- False important for the test
        f3=MEDLoader.ReadFieldNode(fileName2,f2.getMesh().getName(),0,f2.getName(),2,7);
        f3.checkCoherency();
        self.assertTrue(f3.isEqual(f2,1e-12,1e-12));
        #
        pass

    def testFieldNodeProfilRW2(self):
        fileName="Pyfile23.med";
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        MEDLoader.WriteUMesh(fileName,mesh,True);
        #
        f1=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
        f1.setName("FieldMix");
        f1.setMesh(mesh);
        arr2=[1071.,1171.,1010.,1110.,1020.,1120.,1030.,1130.,1040.,1140.,1050.,1150.,
              1060.,1160.,1070.,1170.,1080.,1180.,1090.,1190.,1091.,1191.,1092.,1192.];
        array=DataArrayDouble.New();
        array.setValues(arr2,12,2);
        f1.setArray(array);
        array.setInfoOnComponent(0,"plkj (mm)");
        array.setInfoOnComponent(1,"pqqqss (mm)");
        tmp=array.getPointer();
        f1.setTime(3.17,2,7);
        #
        renumArr=[3,7,2,1,5,11,10,0,9,6,8,4]
        f1.renumberNodes(renumArr);
        f1.checkCoherency();
        MEDLoader.WriteField(fileName,f1,False);#<- False important for the test
        f2=MEDLoader.ReadFieldNode(fileName,f1.getMesh().getName(),0,f1.getName(),2,7);
        self.assertTrue(f2.isEqual(f1,1e-12,1e-12));
        #
        pass

    def testMixCellAndNodesFieldRW1(self):
        fileName="Pyfile21.med";
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setName("FieldMix");
        f1.setMesh(mesh);
        array=DataArrayDouble.New();
        f1.setArray(array);
        arr1=[71.,171.,10.,110.,20.,120.,30.,130.,40.,140.,50.,150.]
        array.setValues(arr1,6,2);
        array.setInfoOnComponent(0,"plkj (mm)");
        array.setInfoOnComponent(1,"pqqqss (mm)");
        f1.setTime(3.14,2,7);
        f1.checkCoherency();
        #
        f2=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
        f2.setName("FieldMix");
        f2.setMesh(mesh);
        array=DataArrayDouble.New();
        f2.setArray(array);
        arr2=[1071.,1171.,1010.,1110.,1020.,1120.,1030.,1130.,1040.,1140.,1050.,1150.,
              1060.,1160.,1070.,1170.,1080.,1180.,1090.,1190.,1091.,1191.,1092.,1192.]
        array.setValues(arr2,12,2)
        array.setInfoOnComponent(0,"plkj (mm)");
        array.setInfoOnComponent(1,"pqqqss (mm)");
        f2.setTime(3.17,2,7);
        f2.checkCoherency();
        #
        MEDLoader.WriteField(fileName,f1,True);
        ts=MEDLoader.GetTypesOfField(fileName,f1.getName(),f1.getMesh().getName());
        self.assertEqual(1,len(ts));
        self.assertEqual(ON_CELLS,ts[0]);
        fs=MEDLoader.GetAllFieldNamesOnMesh(fileName,f1.getMesh().getName());
        self.assertEqual(1,len(fs));
        self.assertTrue(fs[0]=="FieldMix");
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f2);
        fs=MEDLoader.GetAllFieldNamesOnMesh(fileName,f1.getMesh().getName());
        self.assertEqual(1,len(fs));
        self.assertTrue(fs[0]=="FieldMix");
        #
        ts=MEDLoader.GetTypesOfField(fileName,f1.getName(),f1.getMesh().getName());
        self.assertEqual(2,len(ts));
        self.assertEqual(ON_NODES,ts[0]);
        self.assertEqual(ON_CELLS,ts[1]);
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
        f1=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME);
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
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
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
    pass

unittest.main()
