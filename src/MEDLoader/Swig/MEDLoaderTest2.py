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

from MEDLoader import *
import unittest
from math import pi,e,sqrt
from MEDLoaderDataForTest import MEDLoaderDataForTest

class MEDLoaderTest2(unittest.TestCase):
    def testMesh1DRW(self):
        mesh=MEDLoaderDataForTest.build1DMesh_1();
        mesh.checkConsistencyLight();
        WriteUMeshDep("Pyfile1.med",mesh,False);
        mesh_rw=ReadUMeshFromFile("Pyfile1.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMesh2DCurveRW(self):
        mesh=MEDLoaderDataForTest.build2DCurveMesh_1();
        mesh.checkConsistencyLight();
        WriteUMeshDep("Pyfile2.med",mesh,False);
        mesh_rw=ReadUMeshFromFile("Pyfile2.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMesh2DRW(self):
        mesh=MEDLoaderDataForTest.build2DMesh_1();
        mesh.checkConsistencyLight();
        WriteUMeshDep("Pyfile3.med",mesh,False);
        mesh_rw=ReadUMeshFromFile("Pyfile3.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMesh3DSurfRW(self):
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        mesh.checkConsistencyLight();
        WriteUMeshDep("Pyfile4.med",mesh,False);
        mesh_rw=ReadUMeshFromFile("Pyfile4.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testMesh3DRW(self):
        mesh=MEDLoaderDataForTest.build3DMesh_1();
        mesh.checkConsistencyLight();
        WriteUMeshDep("Pyfile5.med",mesh,False);
        mesh_rw=ReadUMeshFromFile("Pyfile5.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    def testFieldRW1(self):
        f1=MEDLoaderDataForTest.buildVecFieldOnCells_1();
        WriteFieldDep("Pyfile6.med",f1,False);
        f2=ReadFieldCell("Pyfile6.med",f1.getMesh().getName(),0,f1.getName(),0,1);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        #
        f1=MEDLoaderDataForTest.buildVecFieldOnNodes_1();
        WriteFieldDep("Pyfile7.med",f1,False);
        f2=ReadFieldNode("Pyfile7.med",f1.getMesh().getName(),0,f1.getName(),2,3);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        pass

    def testFieldRW2(self):
        fileName="Pyfile8.med";
        VAL1=12345.67890314;
        VAL2=-1111111111111.;
        f1=MEDLoaderDataForTest.buildVecFieldOnCells_1();
        WriteFieldDep(fileName,f1,False);
        f1.setTime(10.,8,9);
        f1.getArray().setIJ(0,0,VAL1);
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.setTime(10.14,18,19);
        f1.getArray().setIJ(0,0,VAL2);
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        #retrieving time steps...
        f2=ReadFieldCell(fileName,f1.getMesh().getName(),0,f1.getName(),8,9);
        f1.setTime(10.,8,9);
        f1.getArray().setIJ(0,0,VAL1);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        f2=ReadFieldCell(fileName,f1.getMesh().getName(),0,f1.getName(),0,1);
        f3=MEDLoaderDataForTest.buildVecFieldOnCells_1();
        self.assertTrue(f3.isEqual(f2,1e-12,1e-12));
        f2=ReadFieldCell(fileName,f1.getMesh().getName(),0,f1.getName(),18,19);
        f1.setTime(10.14,18,19);
        f1.getArray().setIJ(0,0,VAL2);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        #ON NODES
        f1=MEDLoaderDataForTest.buildVecFieldOnNodes_1();
        fileName2="Pyfile9.med";
        WriteFieldDep(fileName2,f1,False);
        f1.setTime(110.,108,109);
        tmp=f1.getArray().getPointer();
        f1.getArray().setIJ(0,3,VAL1);
        WriteFieldUsingAlreadyWrittenMesh(fileName2,f1);
        f1.setTime(210.,208,209);
        f1.getArray().setIJ(0,3,VAL2);
        WriteFieldUsingAlreadyWrittenMesh(fileName2,f1);
        f2=ReadFieldNode(fileName2,f1.getMesh().getName(),0,f1.getName(),108,109);
        f1.setTime(110.,108,109);
        f1.getArray().setIJ(0,3,VAL1);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12));
        f2=ReadFieldNode(fileName2,f1.getMesh().getName(),0,f1.getName(),2,3);
        f3=MEDLoaderDataForTest.buildVecFieldOnNodes_1();
        self.assertTrue(f3.isEqual(f2,1e-12,1e-12));
        f2=ReadFieldNode(fileName2,f1.getMesh().getName(),0,f1.getName(),208,209);
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
        WriteFieldDep(fileName,f1,False);
        f1.setTime(10.14,18,19);
        f1.getArray().setIJ(0,0,VAL2);
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.getMesh().setName(name3);
        f1.setTime(10.55,28,29);
        f1.getArray().setIJ(0,0,3*VAL1);
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.setTime(10.66,38,39);
        f1.getArray().setIJ(0,0,3*VAL2);
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.setTime(10.77,48,49);
        f1.getArray().setIJ(0,0,4*VAL2);
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        #ON NODES
        f1=MEDLoaderDataForTest.buildVecFieldOnNodes_1();
        f1.setName(name1);
        f1.getMesh().setName(name3);
        f1.setTime(110.,8,9);
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.setTime(110.,108,109);
        tmp=f1.getArray().getPointer();
        f1.getArray().setIJ(0,3,VAL1);
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.setTime(210.,208,209);
        f1.getArray().setIJ(0,3,VAL2);
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        #
        it1=GetCellFieldIterations(fileName,name3,name1);
        self.assertEqual(5,len(it1));
        self.assertEqual(8,it1[0][0]); self.assertEqual(9,it1[0][1]);
        self.assertEqual(18,it1[1][0]); self.assertEqual(19,it1[1][1]);
        self.assertEqual(28,it1[2][0]); self.assertEqual(29,it1[2][1]);
        self.assertEqual(38,it1[3][0]); self.assertEqual(39,it1[3][1]);
        self.assertEqual(48,it1[4][0]); self.assertEqual(49,it1[4][1]);
        it3=GetNodeFieldIterations(fileName,name3,name1);
        self.assertEqual(3,len(it3));
        self.assertEqual(8,it3[0][0]); self.assertEqual(9,it3[0][1]);
        self.assertEqual(108,it3[1][0]); self.assertEqual(109,it3[1][1]);
        self.assertEqual(208,it3[2][0]); self.assertEqual(209,it3[2][1]);
        #
        #
        f1=ReadFieldCell(fileName,name3,0,name1,8,9);
        self.assertAlmostEqual(VAL1,f1.getArray().getIJ(0,0),13);
        f1=ReadFieldCell(fileName,name3,0,name1,18,19);
        self.assertAlmostEqual(VAL2,f1.getArray().getIJ(0,0),13);
        f1=ReadFieldCell(fileName,name3,0,name1,28,29);
        self.assertAlmostEqual(3*VAL1,f1.getArray().getIJ(0,0),13);
        f1=ReadFieldCell(fileName,name3,0,name1,38,39);
        self.assertAlmostEqual(3*VAL2,f1.getArray().getIJ(0,0),13);
        f1=ReadFieldCell(fileName,name3,0,name1,48,49);
        self.assertAlmostEqual(4*VAL2,f1.getArray().getIJ(0,0),13);
        #
        f1=ReadFieldNode(fileName,name3,0,name1,8,9);
        self.assertAlmostEqual(71.,f1.getArray().getIJ(0,3),13);
        f1=ReadFieldNode(fileName,name3,0,name1,108,109);
        self.assertAlmostEqual(VAL1,f1.getArray().getIJ(0,3),13);
        f1=ReadFieldNode(fileName,name3,0,name1,208,209);
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
        WriteUMeshesPartitionDep(fileName,mnane,meshes,False);
        #
        mesh5=ReadUMeshFromFile(fileName,mnane);
        mesh1.setName(mnane);
        part3=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]
        mesh6=mesh5.buildPartOfMySelf(part3,True);
        mesh6.setName(mnane);
        self.assertTrue(mesh6.isEqual(mesh1,1e-12));
        grps=GetMeshGroupsNames(fileName,mnane);
        self.assertEqual(4,len(grps));
        grps.index("mesh2");
        grps.index("mesh3");
        grps.index("mesh4");
        grps.index("3DMesh_1");
        #
        vec=["mesh2"];
        mesh2_2=ReadUMeshFromGroups(fileName,mnane,0,vec);
        self.assertTrue(mesh2_2.isEqual(mesh2,1e-12));
        vec=["mesh3"];
        mesh3_2=ReadUMeshFromGroups(fileName,mnane,0,vec);
        self.assertTrue(mesh3_2.isEqual(mesh3,1e-12));
        vec=["mesh4"];
        mesh4_2=ReadUMeshFromGroups(fileName,mnane,0,vec);
        self.assertTrue(mesh4_2.isEqual(mesh4,1e-12));
        vec=["3DMesh_1"];
        mesh1_2=ReadUMeshFromGroups(fileName,mnane,0,vec);
        mesh1.setName("3DMesh_1");
        self.assertTrue(mesh1_2.isEqual(mesh1,1e-12));
        #
        vec=["Family_-5","Family_-3"];
        mesh2_2=ReadUMeshFromFamilies(fileName,mnane,0,vec);
        mesh2_2.setName("mesh2");
        self.assertTrue(mesh2_2.isEqual(mesh2,1e-12));
        pass

    def testMesh3DSurfShuffleRW(self):
        fileName="Pyfile15.med";
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        renumber1=[2,5,1,0,3,4]
        mesh.renumberCells(renumber1,False);
        mesh.checkConsistencyLight();
        WriteUMeshDep(fileName,mesh,False);
        mesh_rw=ReadUMeshFromFile(fileName,mesh.getName(),0);
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
        WriteUMeshDep(fileName,m,False);
        f1Tmp=m.getMeasureField(False);
        f1=f1Tmp.buildNewTimeReprFromThis(ONE_TIME,False);
        f1.setTime(0.,1,2);
        f_1=f1.cloneWithMesh(True);
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.applyFunc("2*x");
        f1.setTime(0.01,3,4);
        f_2=f1.cloneWithMesh(True);
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f1.applyFunc("2*x/3");
        f1.setTime(0.02,5,6);
        f_3=f1.cloneWithMesh(True);
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        # Reading
        its=[(1,2),(3,4),(5,6)];
        fs=ReadFieldsOnSameMesh(ON_CELLS,fileName,f_1.getMesh().getName(),0,f_1.getName(),its);
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
        WriteUMeshes(fileName,meshes,False);
        m3d_bis=ReadUMeshFromFile(fileName,m2d.getName(),0);
        self.assertTrue(not m3d_bis.isEqual(m3d,1e-12));
        m3d_bis.setName(m3d.getName());
        self.assertTrue(m3d_bis.isEqual(m3d,1e-12));
        m2d_bis=ReadUMeshFromFile(fileName,m2d.getName(),-1);#-1 for faces
        self.assertTrue(m2d_bis.isEqual(m2d,1e-12));
        # Creation of a field on faces.
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME);
        f1.setName("FieldOnFacesShuffle");
        f1.setMesh(m2d);
        array=DataArrayDouble.New();
        arr1=[71.,171.,10.,110.,20.,120.,30.,130.,40.,140.]
        array.setValues(arr1,m2d.getNumberOfCells(),2);
        array.setInfoOnComponent(0,"plkj [mm]");
        array.setInfoOnComponent(1,"pqqqss [mm]");
        f1.setArray(array);
        tmp=array.setValues(arr1,m2d.getNumberOfCells(),2);
        f1.setTime(3.14,2,7);
        f1.checkConsistencyLight();
        WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
        f2=ReadFieldCell(fileName,f1.getMesh().getName(),-1,f1.getName(),2,7);
        self.assertTrue(f2.isEqual(f1,1e-12,1e-12));
        pass
    pass

if __name__ == "__main__":
  unittest.main()
