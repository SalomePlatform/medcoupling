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
# Author : Anthony Geay (CEA/DEN)

from MEDLoader import *
import unittest
from math import pi,e,sqrt
from MEDLoaderDataForTest import MEDLoaderDataForTest,WriteInTmpDir
from MEDLoaderDataForTest import TestWriteUMeshesRW1,TestMultiFieldShuffleRW1

class MEDLoaderTest2(unittest.TestCase):
    @WriteInTmpDir
    def testMesh1DRW(self):
        mesh=MEDLoaderDataForTest.build1DMesh_1();
        mesh.checkConsistencyLight();
        WriteUMeshDep("Pyfile1.med",mesh,False);
        mesh_rw=ReadUMeshFromFile("Pyfile1.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    @WriteInTmpDir
    def testMesh2DCurveRW(self):
        mesh=MEDLoaderDataForTest.build2DCurveMesh_1();
        mesh.checkConsistencyLight();
        WriteUMeshDep("Pyfile2.med",mesh,False);
        mesh_rw=ReadUMeshFromFile("Pyfile2.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    @WriteInTmpDir
    def testMesh2DRW(self):
        mesh=MEDLoaderDataForTest.build2DMesh_1();
        mesh.checkConsistencyLight();
        WriteUMeshDep("Pyfile3.med",mesh,False);
        mesh_rw=ReadUMeshFromFile("Pyfile3.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    @WriteInTmpDir
    def testMesh3DSurfRW(self):
        mesh=MEDLoaderDataForTest.build3DSurfMesh_1();
        mesh.checkConsistencyLight();
        WriteUMeshDep("Pyfile4.med",mesh,False);
        mesh_rw=ReadUMeshFromFile("Pyfile4.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    @WriteInTmpDir
    def testMesh3DRW(self):
        mesh=MEDLoaderDataForTest.build3DMesh_1();
        mesh.checkConsistencyLight();
        WriteUMeshDep("Pyfile5.med",mesh,False);
        mesh_rw=ReadUMeshFromFile("Pyfile5.med",mesh.getName(),0);
        self.assertTrue(mesh.isEqual(mesh_rw,1e-12));
        pass

    @WriteInTmpDir
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

    @WriteInTmpDir
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
    @WriteInTmpDir
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

    @WriteInTmpDir
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

    @WriteInTmpDir
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

    @WriteInTmpDir
    def testMultiFieldShuffleRW1(self):
        TestMultiFieldShuffleRW1(self)
        pass

    @WriteInTmpDir
    def testWriteUMeshesRW1(self):
        TestWriteUMeshesRW1(self)
        pass
    pass

if __name__ == "__main__":
  unittest.main()
