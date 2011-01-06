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

from MEDLoader import *
import unittest
from math import pi,e,sqrt
from MEDLoaderDataForTest import MEDLoaderDataForTest

class MEDLoaderTest(unittest.TestCase):
    def testMEDMesh1(self):
        fileName="Pyfile18.med"
        mname="ExampleOfMultiDimW"
        medmesh=MEDFileUMesh.New(fileName,mname)
        self.assertEqual((0,-1),medmesh.getNonEmptyLevels())
        m1_0=medmesh.getLevel0Mesh()
        m1_1=MEDLoader.ReadUMeshFromFile(fileName,mname,0)
        self.assertTrue(m1_0.isEqual(m1_1,1e-12));
        m2_0=medmesh.getLevelM1Mesh()
        m2_1=MEDLoader.ReadUMeshFromFile(fileName,mname,-1)
        self.assertTrue(m2_0.isEqual(m2_1,1e-12));
        pass
    def testMEDMesh2(self):
        fileName="Pyfile10.med"
        mname="3DToto"
        outFileName="MEDFileMesh1.med"
        medmesh=MEDFileUMesh.New(fileName,mname)
        self.assertEqual((0,),medmesh.getNonEmptyLevels())
        m1_0=medmesh.getLevel0Mesh()
        m1_1=MEDLoader.ReadUMeshFromFile(fileName,mname,0)
        self.assertTrue(m1_0.isEqual(m1_1,1e-12));
        g1_0=medmesh.getGroup(0,"mesh2")
        g1_1=MEDLoader.ReadUMeshFromGroups(fileName,mname,0,["mesh2"]);
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getGroup(0,"mesh3")
        g1_1=MEDLoader.ReadUMeshFromGroups(fileName,mname,0,["mesh3"]);
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getGroups(0,["mesh3","mesh2"])
        g1_1=MEDLoader.ReadUMeshFromGroups(fileName,mname,0,["mesh3","mesh2"]);
        g1_1.setName(g1_0.getName())
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getFamily(0,"Family_2")
        g1_1=MEDLoader.ReadUMeshFromFamilies(fileName,mname,0,["Family_2"]);
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getFamilies(0,["Family_2","Family_4"])
        g1_1=MEDLoader.ReadUMeshFromFamilies(fileName,mname,0,["Family_2","Family_4"]);
        g1_1.setName(g1_0.getName())
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        medmesh.write(outFileName,2);
        self.assertEqual([2,3,5,14,16],medmesh.getGroupArr(0,"mesh2").getValues());
        self.assertEqual([2,3,16],medmesh.getFamilyArr(0,"Family_2").getValues());
        self.assertEqual([2,3,5,14,16],medmesh.getFamiliesArr(0,["Family_4","Family_2"]).getValues());
        self.assertEqual([19,2,3,4,5,14,15,16],medmesh.getGroupsArr(0,["mesh2","mesh4","mesh3"]).getValues());
        famn=medmesh.getFamilyNameGivenId(0)
        self.assertEqual(range(60),medmesh.getNodeFamilyArr(famn).getValues());
        #without renum
        self.assertEqual([2,3,5,14,16],medmesh.getGroupArr(0,"mesh2",False).getValues());
        self.assertEqual([2,3,16],medmesh.getFamilyArr(0,"Family_2",False).getValues());
        self.assertEqual([2,3,5,14,16],medmesh.getFamiliesArr(0,["Family_4","Family_2"],False).getValues());
        self.assertEqual([0,2,3,4,5,14,15,16],medmesh.getGroupsArr(0,["mesh2","mesh3","mesh4"],False).getValues());
        self.assertEqual(range(60),medmesh.getNodeFamilyArr(famn,False).getValues());
        pass

    # this tests emulates MEDMEM ( Except that it works ! ) The permutation are NOT taken into account
    def testMEDMesh3(self):
        outFileName="MEDFileMesh3.med"
        c=DataArrayDouble.New()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ];
        targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        c.setValues(coords,9,2)
        m=MEDCouplingUMesh.New();
        m.setMeshDimension(2);
        m.allocateCells(5);
        m.insertNextCell(NORM_TRI3,3,targetConn[4:7])
        m.insertNextCell(NORM_TRI3,3,targetConn[7:10])
        m.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
        m.insertNextCell(NORM_QUAD4,4,targetConn[10:14])
        m.insertNextCell(NORM_QUAD4,4,targetConn[14:18])
        m.finishInsertingCells();
        m.setCoords(c)
        m.checkCoherency()
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(1);
        m1.allocateCells(3);
        m1.insertNextCell(NORM_SEG2,2,[1,4])
        m1.insertNextCell(NORM_SEG2,2,[3,6])
        m1.insertNextCell(NORM_SEG3,3,[2,8,5])
        m1.finishInsertingCells();
        m1.setCoords(c)
        m1.checkCoherency()
        m2=MEDCouplingUMesh.New();
        m2.setMeshDimension(0);
        m2.allocateCells(4);
        m2.insertNextCell(NORM_POINT1,1,[1])
        m2.insertNextCell(NORM_POINT1,1,[3])
        m2.insertNextCell(NORM_POINT1,1,[2])
        m2.insertNextCell(NORM_POINT1,1,[6])
        m2.finishInsertingCells();
        m2.setCoords(c)
        m2.checkCoherency()
        #
        mm=MEDFileUMesh.New()
        mm.setName("MyFirstMEDCouplingMEDmesh")
        mm.setDescription("IHopeToConvinceLastMEDMEMUsers")
        mm.setCoords(c)
        mm.setMeshAtLevelOld(-1,m1);
        mm.setMeshAtLevelOld(0,m);
        mm.setMeshAtLevelOld(-2,m2);
        # playing with groups
        g1_2=DataArrayInt.New()
        g1_2.setValues([1,3],2,1)
        g1_2.setName("G1")
        g2_2=DataArrayInt.New()
        g2_2.setValues([1,2,3],3,1)
        g2_2.setName("G2")
        mm.setGroupsAtLevel(0,[g1_2,g2_2],False)
        g1_1=DataArrayInt.New()
        g1_1.setValues([0,1,2],3,1)
        g1_1.setName("G1")
        g2_1=DataArrayInt.New()
        g2_1.setValues([0,2],2,1)
        g2_1.setName("G2")
        mm.setGroupsAtLevel(-1,[g1_1,g2_1],False)
        g1_N=DataArrayInt.New()
        g1_N.setValues(range(8),8,1)
        g1_N.setName("G1")
        g2_N=DataArrayInt.New()
        g2_N.setValues(range(9),9,1)
        g2_N.setName("G2")
        mm.setGroupsAtLevel(1,[g1_N,g2_N],False)
        # check content of mm
        t=mm.getGroupArr(0,"G1",False)
        self.assertTrue(g1_2.isEqual(t));
        t=mm.getGroupArr(0,"G2",False)
        self.assertTrue(g2_2.isEqual(t));
        t=mm.getGroupArr(-1,"G1",False)
        self.assertTrue(g1_1.isEqual(t));
        t=mm.getGroupArr(-1,"G2",False)
        self.assertTrue(g2_1.isEqual(t));
        t=mm.getGroupArr(1,"G1",False)
        self.assertTrue(g1_N.isEqual(t));
        t=mm.getGroupArr(1,"G2",False)
        self.assertTrue(g2_N.isEqual(t));
        #
        mm.write(outFileName,2);
        pass

    # this test is the testMEDMesh3 except that permutation is dealed here
    def testMEDMesh4(self):
        outFileName="MEDFileMesh4.med"
        c=DataArrayDouble.New()
        coords=[-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ];
        targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        c.setValues(coords,9,2)
        m=MEDCouplingUMesh.New();
        m.setMeshDimension(2);
        m.allocateCells(5);
        m.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
        m.insertNextCell(NORM_TRI3,3,targetConn[4:7])
        m.insertNextCell(NORM_TRI3,3,targetConn[7:10])
        m.insertNextCell(NORM_QUAD4,4,targetConn[10:14])
        m.insertNextCell(NORM_QUAD4,4,targetConn[14:18])
        m.finishInsertingCells();
        m.setCoords(c)
        m.checkCoherency()
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(1);
        m1.allocateCells(3);
        m1.insertNextCell(NORM_SEG2,2,[1,4])
        m1.insertNextCell(NORM_SEG3,3,[2,8,5])
        m1.insertNextCell(NORM_SEG2,2,[3,6])
        m1.finishInsertingCells();
        m1.setCoords(c)
        m1.checkCoherency()
        m2=MEDCouplingUMesh.New();
        m2.setMeshDimension(0);
        m2.allocateCells(4);
        m2.insertNextCell(NORM_POINT1,1,[1])
        m2.insertNextCell(NORM_POINT1,1,[3])
        m2.insertNextCell(NORM_POINT1,1,[2])
        m2.insertNextCell(NORM_POINT1,1,[6])
        m2.finishInsertingCells();
        m2.setCoords(c)
        m2.checkCoherency()
        #
        mm=MEDFileUMesh.New()
        mm.setName("My2ndMEDCouplingMEDmesh")
        mm.setDescription("ThisIsImpossibleToDoWithMEDMEM")
        mm.setCoords(c)
        renumNode=DataArrayInt.New()
        renumNode.setValues([10,11,12,13,14,15,16,17,18],9,1)
        mm.setRenumArr(1,renumNode)
        mm.setMeshAtLevel(-1,m1);
        mm.setMeshAtLevel(0,m);
        mm.setMeshAtLevel(-2,m2);
        # playing with groups
        g1_2=DataArrayInt.New()
        g1_2.setValues([2,3],2,1)
        g1_2.setName("G1")
        g2_2=DataArrayInt.New()
        g2_2.setValues([2,0,3],3,1)
        g2_2.setName("G2")
        mm.setGroupsAtLevel(0,[g1_2,g2_2],True)
        g1_1=DataArrayInt.New()
        g1_1.setValues([0,2,1],3,1)
        g1_1.setName("G1")
        g2_1=DataArrayInt.New()
        g2_1.setValues([0,2],2,1)
        g2_1.setName("G2")
        mm.setGroupsAtLevel(-1,[g1_1,g2_1],True)
        g1_N=DataArrayInt.New()
        g1_N.setValues([10,11,12,13,14,15,16,17],8,1)
        g1_N.setName("G1")
        g2_N=DataArrayInt.New()
        g2_N.setValues([10,11,12,13,14,15,16,17,18],9,1)
        g2_N.setName("G2")
        mm.setGroupsAtLevel(1,[g1_N,g2_N],True)
        # check content of mm
        t=mm.getGroupArr(0,"G1",True)
        self.assertTrue(g1_2.isEqual(t));
        t=mm.getGroupArr(0,"G2",True)
        self.assertTrue(g2_2.isEqual(t));
        t=mm.getGroupArr(-1,"G1",True)
        self.assertTrue(g1_1.isEqual(t));
        t=mm.getGroupArr(-1,"G2",True)
        self.assertTrue(g2_1.isEqual(t));
        #
        mm.write(outFileName,2);
        pass
    pass

unittest.main()
