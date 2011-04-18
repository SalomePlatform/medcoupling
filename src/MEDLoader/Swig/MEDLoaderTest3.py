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
        medmesh=MEDFileMesh.New(fileName,mname)
        self.assertEqual((0,-1),medmesh.getNonEmptyLevels())
        m1_0=medmesh.getLevel0Mesh(True)
        m1_1=MEDLoader.ReadUMeshFromFile(fileName,mname,0)
        self.assertTrue(m1_0.isEqual(m1_1,1e-12));
        m2_0=medmesh.getLevelM1Mesh(True)
        m2_1=MEDLoader.ReadUMeshFromFile(fileName,mname,-1)
        self.assertTrue(m2_0.isEqual(m2_1,1e-12));
        pass
    def testMEDMesh2(self):
        fileName="Pyfile10.med"
        mname="3DToto"
        outFileName="MEDFileMesh1.med"
        medmesh=MEDFileUMesh.New(fileName,mname)
        self.assertEqual((0,),medmesh.getNonEmptyLevels())
        m1_0=medmesh.getLevel0Mesh(True)
        m1_1=MEDLoader.ReadUMeshFromFile(fileName,mname,0)
        self.assertTrue(m1_0.isEqual(m1_1,1e-12));
        g1_0=medmesh.getGroup(0,"mesh2",True)
        g1_1=MEDLoader.ReadUMeshFromGroups(fileName,mname,0,["mesh2"]);
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getGroup(0,"mesh3",True)
        g1_1=MEDLoader.ReadUMeshFromGroups(fileName,mname,0,["mesh3"]);
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getGroups(0,["mesh3","mesh2"])
        g1_1=MEDLoader.ReadUMeshFromGroups(fileName,mname,0,["mesh3","mesh2"]);
        g1_1.setName(g1_0.getName())
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getFamily(0,"Family_2",True)
        g1_1=MEDLoader.ReadUMeshFromFamilies(fileName,mname,0,["Family_2"]);
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getFamilies(0,["Family_2","Family_4"],True)
        g1_1=MEDLoader.ReadUMeshFromFamilies(fileName,mname,0,["Family_2","Family_4"]);
        g1_1.setName(g1_0.getName())
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        medmesh.write(outFileName,2);
        self.assertEqual([2,3,5,14,16],medmesh.getGroupArr(0,"mesh2",True).getValues());
        self.assertEqual([2,3,16],medmesh.getFamilyArr(0,"Family_2",True).getValues());
        self.assertEqual([2,3,5,14,16],medmesh.getFamiliesArr(0,["Family_4","Family_2"],True).getValues());
        self.assertEqual([19,2,3,4,5,14,15,16],medmesh.getGroupsArr(0,["mesh2","mesh4","mesh3"],True).getValues());
        famn=medmesh.getFamilyNameGivenId(0)
        self.assertRaises(InterpKernelException,medmesh.getNodeFamilyArr,famn,True);
        #without renum
        self.assertEqual([2,3,5,14,16],medmesh.getGroupArr(0,"mesh2",False).getValues());
        self.assertEqual([2,3,16],medmesh.getFamilyArr(0,"Family_2",False).getValues());
        self.assertEqual([2,3,5,14,16],medmesh.getFamiliesArr(0,["Family_4","Family_2"],False).getValues());
        self.assertEqual([0,2,3,4,5,14,15,16],medmesh.getGroupsArr(0,["mesh2","mesh3","mesh4"],False).getValues());
        self.assertRaises(InterpKernelException,medmesh.getNodeFamilyArr,famn,False);
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
        c.setInfoOnComponent(0,"abcdef [km]")
        c.setInfoOnComponent(1,"ghij [MW]")
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
        mm.setRenumFieldArr(1,renumNode)
        mm.setMeshAtLevel(-1,m1,True);
        mm.setMeshAtLevel(0,m,True);
        mm.setMeshAtLevel(-2,m2,True);
        mm.removeMeshAtLevel(-2)
        mm.setMeshAtLevel(-2,m2,True);
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
        mm2=MEDFileMesh.New(outFileName)
        res=mm.isEqual(mm2,1e-12)
        self.assertTrue(res[0])
        pass

    #testing persistence of retrieved arrays
    def testMEDMesh5(self):
        fileName="Pyfile18.med"
        mname="ExampleOfMultiDimW"
        medmesh=MEDFileUMesh.New(fileName,mname)
        m1_0=medmesh.getLevel0Mesh(True)
        da1=medmesh.getFamilyFieldAtLevel(0)
        del medmesh
        self.assertEqual(20,m1_0.getNumberOfCells())
        self.assertEqual(20,da1.getNumberOfTuples())
        pass

    def testMEDMesh6(self):
        outFileName="MEDFileMesh5.med"
        m=MEDFileCMesh.New()
        m.setTime(2.3,-1,-1)
        m1=MEDCouplingCMesh.New();
        da=DataArrayDouble.New()
        da.setValues([0.,1.,2.],3,1)
        da.setInfoOnComponent(0,"XX [mm]")
        m1.setCoordsAt(0,da)
        da=DataArrayDouble.New()
        da.setValues([0.,1.2],2,1)
        da.setInfoOnComponent(0,"YY [km]")
        m1.setCoordsAt(1,da)
        da=DataArrayDouble.New()
        da.setValues([0.,1.3],2,1)
        da.setInfoOnComponent(0,"ZZ [um]")
        m1.setCoordsAt(2,da)
        m.setMesh(m1)
        m.setName("myFirstCartMesh")
        m.setDescription("mmmmpppppppp")
        m.setTimeValue(2.3)
        m.setTimeUnit("ms")
        da=DataArrayInt.New()
        da.setValues([0,0,1,0,1,2,4,3,0,1,2,2],12,1)
        m.setFamilyFieldArr(1,da)
        m.setFamilyId("family1",1)
        da=m.getFamilyArr(1,"family1")
        expected1=[2,4,9]
        self.assertEqual(expected1,da.getValues())
        m.write(outFileName,2);
        mm=MEDFileMesh.New(outFileName)
        self.assertTrue(m.isEqual(mm,1e-12)[0])
        self.assertEqual(expected1,mm.getFamilyArr(1,"family1").getValues())
        m2=mm.getMesh()
        tt=m.getTime()
        m1.setTime(tt[0],tt[1],tt[2])
        m1.setName(m.getName())
        m1.setTimeUnit(m.getTimeUnit())
        m1.setDescription(m.getDescription())
        self.assertTrue(m2.isEqual(m1,1e-12));
        pass

    def testMEDMesh7(self):
       fileName="Pyfile24.med"
       m2,m1,m0,f2,f1,f0,p,n2,n1,n0,fns,fids,grpns,famIdsPerGrp=MEDLoaderDataForTest.buildMultiLevelMesh_1()
       m=MEDFileUMesh.New()
       m.setCoords(m2.getCoords())
       m.setMeshAtLevel(0,m2)
       m.setMeshAtLevel(-1,m1)
       m.setMeshAtLevel(-2,m0)
       m.setFamilyFieldArr(0,f2)
       m.setFamilyFieldArr(-1,f1)
       m.setFamilyFieldArr(-2,f0)
       m.setFamilyFieldArr(1,p)
       m.setRenumFieldArr(0,n2)
       m.setRenumFieldArr(-1,n1)
       m.setRenumFieldArr(-2,n0)
       nbOfFams=len(fns)
       for i in xrange(nbOfFams):
           m.addFamily(fns[i],fids[i])
           pass
       nbOfGrps=len(grpns)
       for i in xrange(nbOfGrps):
           m.setFamiliesIdsOnGroup(grpns[i],famIdsPerGrp[i])
           pass
       m.setName(m2.getName())
       m.setDescription(m2.getDescription())
       m.write(fileName,2)
       pass
    pass

unittest.main()
