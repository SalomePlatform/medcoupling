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
        m.insertNextCell(NORM_POLYGON,4,targetConn[10:14])
        m.insertNextCell(NORM_POLYGON,4,targetConn[14:18])
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
        mm.setMeshAtLevel(-1,m1);
        mm.setMeshAtLevel(0,m);
        mm.setMeshAtLevel(-2,m2);
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
        mm.createGroupOnAll(0,"GrpOnAllCell")
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
        self.assertTrue(mm.existsGroup("GrpOnAllCell"));
        t=mm.getGroupArr(0,"GrpOnAllCell")
        self.assertTrue(t.getValues()==range(5))
        #
        mm.write(outFileName,2);
        #
        mm=MEDFileMesh.New(outFileName)
        mbis=mm.getMeshAtLevel(0)
        m.setName(mm.getName()) ; m.setDescription(mm.getDescription())
        self.assertTrue(m.isEqual(mbis,1e-12));
        #
        self.assertEqual(([[(3, 2), (4, 1), (5, 8)], [(1, 2), (2, 1)], [(0, 4)]], 2, 2, 9),MEDLoader.GetUMeshGlobalInfo(outFileName,"MyFirstMEDCouplingMEDmesh"))
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
        self.assertTrue(not mm.existsGroup("GrpOnAllCell"));
        #
        mm.write(outFileName,2);
        mm2=MEDFileMesh.New(outFileName)
        res=mm.isEqual(mm2,1e-12)
        self.assertTrue(res[0])
        l=list(mm2.getFamiliesOnGroup("G2")) ; l.sort()
        self.assertEqual(['Family_10','Family_11','Family_3','Family_4','Family_7'],l)
        mm2.keepFamIdsOnlyOnLevs([3],[-1])
        for lev in mm.getGrpNonEmptyLevelsExt("G2"):
            self.assertEqual(mm.getGroupArr(lev,"G2").getValues(),mm2.getGroupArr(lev,"G2").getValues())
            pass
        l=list(mm2.getFamiliesOnGroup("G2")) ; l.sort()
        self.assertEqual(['Family_10','Family_11','Family_12','Family_3','Family_4','Family_7'],l)
        #
        self.assertEqual([7,7,6],mm2.getFamilyFieldAtLevel(-1).getValues())
        mm2.getFamilyFieldAtLevel(-1).setIJ(1,0,8)
        self.assertEqual([7,8,6],mm2.getFamilyFieldAtLevel(-1).getValues())
        self.assertTrue(not mm2.existsFamily("Family_8"))
        mm2.createGroupOnAll(-1,"GrpOnAllFace")
        self.assertTrue(mm2.existsFamily("Family_8"))
        self.assertEqual(range(3),mm2.getGroupArr(-1,"GrpOnAllFace").getValues())
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
        m.setTime(-1,-1,2.3)
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
        m1.setTime(tt[2],tt[0],tt[1])
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
        #
        self.assertEqual((-1,),m.getGrpNonEmptyLevels("A2A4"))
        self.assertEqual((),m.getGrpNonEmptyLevels("A1"))
        self.assertEqual((-2,),m.getGrpNonEmptyLevels("AP2"))
        self.assertEqual((-1,-2),m.getGrpsNonEmptyLevels(["A2A4","AP2"]))
        self.assertEqual((-1,),m.getFamNonEmptyLevels('A4A3____________________________'))
        self.assertEqual((0,),m.getFamNonEmptyLevels('MESH____DALT3___DALLE___________'))
        self.assertEqual((0,-1,),m.getFamsNonEmptyLevels(['MESH____DALT3___DALLE___________','A4A3____________________________']))
        self.assertEqual(('A1A2','A2A4','A3A1','A3C5','A4A3','B1C1','B2B4','B3B1','B4C3','C1C4','C2B2','C3C2','C4B3','C5A4'),m.getGroupsOnSpecifiedLev(-1))
        self.assertEqual(('DALLE','DALQ1','DALQ2','DALT3','MESH'),m.getGroupsOnSpecifiedLev(0))
        #
        m.write(fileName,2)
        pass

    def funcToTestDelItem(self,ff):
        del ff[[0.02,(3,4)]]
        pass

    #emulation of pointe.med file.
    def testMEDField1(self):
        mm=MEDFileMesh.New("Pyfile17.med")
        mm.write("Pyfile17_bis.med",2)
        ff=MEDFileFieldMultiTS("Pyfile17.med")
        tsExpected=[[1,2],[3,4],[5,6]]
        for pos,f1ts in enumerate(ff):
            self.assertEqual(tsExpected[pos],f1ts.getTime()[:2])
            self.assertEqual(type(f1ts),MEDFileField1TS)
            pass
        self.assertEqual("MeasureOfMesh_Extruded",ff.getName())
        self.assertEqual([3,4],ff[1].getTime()[:-1])
        self.assertEqual([3,4],ff[3,4].getTime()[:-1])
        self.assertEqual([3,4],ff[0.01].getTime()[:-1])
        ff.write("Pyfile17_bis.med",0)
        #
        ts=ff.getTimeSteps() ; ts=[elt[:-1] for elt in ts]
        self.assertEqual([(1,2),(3,4),(5,6)],ts)
        self.funcToTestDelItem(ff)
        ts=ff.getTimeSteps() ; ts=[elt[:-1] for elt in ts]
        self.assertEqual([(1,2)],ts)
        pass

    #profiles
    def testMEDField2(self):
        mm=MEDFileMesh.New("Pyfile19.med")
        mm.write("Pyfile19_bis.med",2)
        ff=MEDFileFieldMultiTS.New("Pyfile19.med","VFieldOnNodes")
        ff.write("Pyfile19_bis.med",0)
        self.assertEqual([('tyty','mm'),('uiop','MW')],MEDLoader.GetComponentsNamesOfField("Pyfile19_bis.med","VFieldOnNodes"))
        pass

    #gauss points
    def testMEDField3(self):
        mm=MEDFileMesh.New("Pyfile13.med")
        mm.write("Pyfile13_bis.med",2)
        ff=MEDFileFieldMultiTS.New("Pyfile13.med","MyFirstFieldOnGaussPoint")
        ff.write("Pyfile13_bis.med",0)
        ff=MEDFileField1TS.New("Pyfile13.med","MyFirstFieldOnGaussPoint",1,5)
        f=ff.getFieldAtLevel(ON_GAUSS_PT,0)
        f2=MEDLoader.ReadFieldGauss("Pyfile13.med",'2DMesh_2',0,'MyFirstFieldOnGaussPoint',1,5)
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        pass

    #gauss NE
    def testMEDField4(self):
        mm=MEDFileMesh.New("Pyfile14.med")
        mm.write("Pyfile14_bis.med",2)
        ff=MEDFileFieldMultiTS.New("Pyfile14.med","MyFieldOnGaussNE")
        ff.write("Pyfile14_bis.med",0)
        ff=MEDFileField1TS.New("Pyfile14.med","MyFieldOnGaussNE",1,5)
        f=ff.getFieldAtLevel(ON_GAUSS_NE,0)
        f2=MEDLoader.ReadFieldGaussNE("Pyfile14.med",'2DMesh_2',0,"MyFieldOnGaussNE",1,5)
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        pass

    # MEDField get/set on pointe.med
    def testMEDField5(self):
        ff=MEDFileField1TS.New("Pyfile17.med","MeasureOfMesh_Extruded",1,2)
        f=ff.getFieldAtLevel(ON_CELLS,0)
        f2=MEDLoader.ReadFieldCell("Pyfile17.med","Extruded",0,"MeasureOfMesh_Extruded",1,2)
        self.assertTrue(f.getMesh().getCoords().isEqual(f2.getMesh().getCoords(),1e-12))
        f.getMesh().tryToShareSameCoords(f2.getMesh(),1e-12)
        f.changeUnderlyingMesh(f2.getMesh(),22,1e-12)
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        # no with renumbering
        f=ff.getFieldAtLevel(ON_CELLS,0,1)
        f2=MEDLoader.ReadFieldCell("Pyfile17.med","Extruded",0,"MeasureOfMesh_Extruded",1,2)
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        f=ff.getFieldAtLevel(ON_CELLS,0,3)
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        f=ff.getFieldAtLevel(ON_CELLS,0,2)
        self.assertTrue(not f.isEqual(f2,1e-12,1e-12))
        f.changeUnderlyingMesh(f2.getMesh(),12,1e-12)
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        pass

    # MEDField get/set on profiles nodes
    def testMEDField6(self):
        ff=MEDFileFieldMultiTS.New("Pyfile7.med","VectorFieldOnNodes")
        its=ff.getIterations()
        self.assertRaises(InterpKernelException,ff.getFieldAtLevel,ON_CELLS,its[0][0],its[0][1],0)# request on cell and it is not on cells
        f=ff.getFieldAtLevel(ON_NODES,its[0][0],its[0][1],0)
        f2=MEDLoader.ReadFieldNode("Pyfile7.med",'3DSurfMesh_1',0,"VectorFieldOnNodes",its[0][0],its[0][1])
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        ff=MEDFileFieldMultiTS.New("Pyfile19.med","VFieldOnNodes")
        its=ff.getIterations()
        f=ff.getFieldAtLevel(ON_NODES,its[0][0],its[0][1],0)
        f2=MEDLoader.ReadFieldNode("Pyfile19.med",'2DMesh_1',0,"VFieldOnNodes",its[0][0],its[0][1])
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        self.assertRaises(InterpKernelException,ff.getFieldAtLevel,ON_CELLS,its[0][0],its[0][1],0)# request on cell and it is not on cells
        self.assertRaises(InterpKernelException,ff.getFieldAtLevel,ON_NODES,its[0][0],its[0][1],0,1)#request renumber following mesh : it is on profile !
        pass

    # MEDField get/set on profiles cells
    def testMEDField7(self):
        ff=MEDFileFieldMultiTS.New("Pyfile12.med","VectorFieldOnCells")
        its=ff.getIterations()
        f=ff.getFieldAtLevel(ON_CELLS,its[0][0],its[0][1],0)
        f2=MEDLoader.ReadFieldCell("Pyfile12.med",'3DMesh_1',0,"VectorFieldOnCells",its[0][0],its[0][1])
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        pass

    #first test of assignation. No profile and types sorted by type.
    def testMEDField8(self):
        fname="Pyfile25.med"
        f1=MEDLoaderDataForTest.buildVecFieldOnCells_1();
        m1=f1.getMesh()
        mm1=MEDFileUMesh.New()
        mm1.setCoords(m1.getCoords())
        mm1.setMeshAtLevel(0,m1)
        mm1.setName(m1.getName())
        mm1.write(fname,2)
        ff1=MEDFileField1TS.New()
        ff1.setFieldNoProfileSBT(f1)
        ff1.write(fname,0)
        f2=MEDLoader.ReadFieldCell(fname,f1.getMesh().getName(),0,f1.getName(),f1.getTime()[1],f1.getTime()[2]);
        itt,orr,ti=ff1.getTime()
        self.assertEqual(0,itt); self.assertEqual(1,orr); self.assertAlmostEqual(2.,ti,14);
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12))
        ff1.setTime(3,4,2.3)
        itt,orr,ti=ff1.getTime()
        self.assertEqual(3,itt); self.assertEqual(4,orr); self.assertAlmostEqual(2.3,ti,14);
        da,infos=ff1.getUndergroundDataArrayExt()
        f2.getArray().setName(da.getName())#da has the same name than f2
        self.assertTrue(da.isEqual(f2.getArray(),1e-12))
        self.assertEqual([((3, 0), (0, 2)), ((4, 0), (2, 4)), ((6, 0), (4, 5)), ((5, 0), (5, 6))],infos)
        #
        fname="Pyfile26.med"
        f1=MEDLoaderDataForTest.buildVecFieldOnNodes_1();
        m1=f1.getMesh()
        mm1=MEDFileUMesh.New()
        mm1.setCoords(m1.getCoords())
        mm1.setMeshAtLevel(0,m1)
        mm1.setName(m1.getName())
        mm1.write(fname,2)
        ff1=MEDFileField1TS.New()
        ff1.setFieldNoProfileSBT(f1)
        nv=1456.
        da=ff1.getUndergroundDataArray().setIJ(0,0,nv)
        ff1.write(fname,0)
        f2=MEDLoader.ReadFieldNode(fname,f1.getMesh().getName(),0,f1.getName(),f1.getTime()[1],f1.getTime()[2])
        self.assertTrue(not f1.isEqual(f2,1e-12,1e-12))
        f1.getArray().setIJ(0,0,nv)
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12))
        #
        fname="Pyfile27.med"
        f1=MEDLoaderDataForTest.buildVecFieldOnGaussNE_1();
        m1=f1.getMesh()
        mm1=MEDFileUMesh.New()
        mm1.setCoords(m1.getCoords())
        mm1.setMeshAtLevel(0,m1)
        mm1.setName(m1.getName())
        mm1.write(fname,2)
        ff1=MEDFileField1TS.New()
        ff1.setFieldNoProfileSBT(f1)
        ff1.write(fname,0)
        f2=MEDLoader.ReadFieldGaussNE(fname,f1.getMesh().getName(),0,f1.getName(),f1.getTime()[1],f1.getTime()[2])
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12))
        da,infos=ff1.getUndergroundDataArrayExt()
        f2.getArray().setName(da.getName())#da has the same name than f2
        self.assertTrue(da.isEqual(f2.getArray(),1e-12))
        self.assertEqual([((3, 0), (0, 6)), ((4, 0), (6, 14)), ((6, 0), (14, 20))],infos)
        #
        fname="Pyfile28.med"
        f1=MEDLoaderDataForTest.buildVecFieldOnGauss_2_Simpler();
        m1=f1.getMesh()
        mm1=MEDFileUMesh.New()
        mm1.setCoords(m1.getCoords())
        mm1.setMeshAtLevel(0,m1)
        mm1.setName(m1.getName())
        mm1.write(fname,2)
        ff1=MEDFileField1TS.New()
        ff1.setFieldNoProfileSBT(f1)
        ff1.write(fname,0)
        ff2=MEDFileField1TS.New(fname,f1.getName(),f1.getTime()[1],f1.getTime()[2])
        f2=ff2.getFieldAtLevel(ON_GAUSS_PT,0)
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12))
        sbt=ff2.getFieldSplitedByType2()
        loc1=ff2.getLocalization("Loc_MyFirstFieldOnGaussPoint_NORM_TRI6_5")
        self.assertEqual("Loc_MyFirstFieldOnGaussPoint_NORM_TRI6_5",loc1.getName())
        self.assertEqual((-1, 1,-1,-1,1,-1,-1,0,0,-1,0,0),loc1.getRefCoords())
        self.assertEqual(6,loc1.getNumberOfPointsInCells())
        self.assertEqual(3,loc1.getNumberOfGaussPoints())
        self.assertEqual(2,loc1.getDimension())
        da,infos=ff2.getUndergroundDataArrayExt()
        f2.getArray().setName(da.getName())#da has the same name than f2
        self.assertTrue(da.isEqual(f2.getArray(),1e-12))
        self.assertEqual(53,da.getNumberOfTuples())
        self.assertEqual([((3, 0), (0, 18)), ((3, 1), (18, 30)), ((3, 2), (30, 36)), ((4, 0), (36, 42)), ((4, 1), (42, 44)), ((6, 0), (44, 53))],infos)
        #
        pass
    
    def testMEDFileData1(self):
        fname="Pyfile29.med"
        d=MEDFileData.New()
        #
        m1=MEDLoaderDataForTest.build1DMesh_1()
        mm1=MEDFileUMesh.New() ; mm1.setCoords(m1.getCoords()) ; mm1.setMeshAtLevel(0,m1) ; mm1.setName(m1.getName())
        mmm1=MEDFileMeshMultiTS.New() ;
        mmm1.setOneTimeStep(mm1)
        m2=MEDLoaderDataForTest.build2DCurveMesh_1()
        mm2=MEDFileUMesh.New() ; mm2.setCoords(m2.getCoords()) ; mm2.setMeshAtLevel(0,m2) ; mm2.setName(m2.getName())
        mmm2=MEDFileMeshMultiTS.New() ; mmm2.setOneTimeStep(mm2)
        ms=MEDFileMeshes.New(); ms.setMeshAtPos(0,mm1) ; ms.setMeshAtPos(1,mm2)
        d.setMeshes(ms)
        for name,mmm in zip(["1DMesh_1","2DCurveMesh_1"],ms):
            self.assertEqual(name,mmm.getName())
            self.assertEqual(type(mmm),MEDFileUMesh)
            pass
        self.assertEqual(('1DMesh_1', '2DCurveMesh_1'),d.getMeshes().getMeshesNames())
        #
        ff1=MEDFileFieldMultiTS.New()
        ff21=MEDFileFieldMultiTS.New()
        ff22=MEDFileFieldMultiTS.New()
        f1=m1.getMeasureField(True) ; f1.setName("f1") ; f1=f1.buildNewTimeReprFromThis(ONE_TIME,False)
        f1.getArray().setInfoOnComponent(0,"power [kW]")
        ff1.appendFieldNoProfileSBT(f1)
        f21=m2.getMeasureField(True) ; f21.setName("f21") ; f21=f21.buildNewTimeReprFromThis(ONE_TIME,False)
        f21.getArray().setInfoOnComponent(0,"sta [mm]") ;
        ff21.appendFieldNoProfileSBT(f21)
        f22=f21.deepCpy() ; f22.setName("f22") ; f22=f22.buildNewTimeReprFromThis(ONE_TIME,False) ;
        f22.applyFunc(2,"3*x*IVec+2*x*JVec")
        f22.getArray().setInfoOnComponent(0,"distance [km]") ; f22.getArray().setInfoOnComponent(1,"displacement [cm]")
        ff22.appendFieldNoProfileSBT(f22)
        fs=MEDFileFields.New()
        fs.pushField(ff1) ; fs.pushField(ff21) ; fs.pushField(ff22)
        for name,fmts in zip(["f1","f21","f22"],fs):
            self.assertEqual(name,fmts.getName())
            pass
        d.setFields(fs)
        #
        fname2="Pyfile29_2.med"
        d.write(fname2,2)
        #
        d2=MEDFileData.New(fname2)
        self.assertEqual(2,d2.getNumberOfMeshes())
        self.assertEqual(3,d2.getNumberOfFields())
        self.assertTrue(isinstance(d2.getMeshes().getMeshAtPos(0),MEDFileUMesh))
        self.assertTrue(isinstance(d2.getMeshes()[0],MEDFileUMesh))
        self.assertTrue(isinstance(d2.getMeshes()['2DCurveMesh_1'],MEDFileUMesh))
        m1bis=d2.getMeshes().getMeshAtPos(0).getMeshAtLevel(0)
        self.assertTrue(m1.isEqual(m1bis,1e-12))
        self.assertEqual(('f1', 'f21', 'f22'),d2.getFields().getFieldsNames())
        self.assertEqual([(-1,-1,0.0)],d2.getFields().getFieldAtPos(2).getTimeSteps())
        self.assertEqual([(-1,-1,0.0)],d2.getFields()[2].getTimeSteps())
        self.assertEqual([(-1,-1,0.0)],d2.getFields().getFieldWithName("f21").getTimeSteps())
        self.assertEqual([(-1,-1,0.0)],d2.getFields()["f21"].getTimeSteps())
        pass
    
    def testMEDField9(self):
        # first test field profile WR. Full type but with some type missing
        fname="Pyfile30.med"
        m1=MEDLoaderDataForTest.build2DMesh_3()
        mm1=MEDFileUMesh.New() ; mm1.setCoords(m1.getCoords()) ; mm1.setMeshAtLevel(0,m1) ;
        mm1.write(fname,2)
        ff1=MEDFileField1TS.New()
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f1.setName("F1")
        d=DataArrayDouble.New() ; d.alloc(2*9,1) ; d.iota(7.); d.rearrange(2); d.setInfoOnComponent(0,"sigX [MPa]") ; d.setInfoOnComponent(1,"sigY [GPa]")
        f1.setArray(d) # note that f1 is NOT defined fully (no mesh !). It is not a bug of test it is too test that MEDFileField1TS.appendFieldProfile is NOT sensible of that.
        da=DataArrayInt.New(); da.alloc(9,1) ; da.iota(0) ; da.setName("sup1")
        #
        ff1.setFieldProfile(f1,mm1,0,da)
        ff1.changePflsNames([(["sup1_NORM_QUAD4"],"ForV650")])
        ff1.write(fname,0)
        #
        vals,pfl=ff1.getFieldWithProfile(ON_CELLS,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))# profiles names cannot be contracted in pfl array name
        self.assertTrue(vals.isEqual(d,1e-14))
        #
        ff2=MEDFileField1TS.New(fname,f1.getName(),-1,-1)
        sbt=ff2.getFieldSplitedByType2()
        self.assertEqual(3,sbt[0][0])#TRI3
        self.assertEqual(0,sbt[0][1][0][0])#CELL For TRI3
        self.assertEqual("",sbt[0][1][0][2])#no profile For TRI3
        self.assertEqual([7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],sbt[0][1][0][1].getValues())# values for TRI3
        self.assertEqual(4,sbt[1][0])#QUAD4
        self.assertEqual(0,sbt[1][1][0][0])#CELL For QUAD4
        self.assertEqual("ForV650",sbt[1][1][0][2])# profile For QUAD4
        self.assertEqual([19, 20, 21, 22, 23, 24],sbt[1][1][0][1].getValues())# values for QUAD4
        self.assertEqual([0],ff2.getTypesOfFieldAvailable())
        vals,pfl=ff2.getFieldWithProfile(ON_CELLS,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(d,1e-14))
        pass
    
    def testMEDField10(self):
        fname="Pyfile31.med"
        m1=MEDLoaderDataForTest.build2DMesh_1()
        m1.renumberCells([0,1,4,2,3,5],False)
        mm1=MEDFileUMesh.New() ; mm1.setCoords(m1.getCoords()) ; mm1.setMeshAtLevel(0,m1) ; mm1.setName(m1.getName())
        mm1.write(fname,2)
        ff1=MEDFileFieldMultiTS.New()
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f1.setName("F2")
        d=DataArrayDouble.New() ; d.alloc(2*4,1) ; d.iota(7.); d.rearrange(2); d.setInfoOnComponent(0,"sigX [MPa]") ; d.setInfoOnComponent(1,"sigY [GPa]")
        f1.setArray(d) # note that f1 is NOT defined fully (no mesh !). It is not a bug of test it is too test that MEDFileField1TS.appendFieldProfile is NOT sensible of that.
        da=DataArrayInt.New(); da.setValues([0,1,2,4],4,1) ; da.setName("sup2")
        #
        ff1.appendFieldProfile(f1,mm1,0,da)
        f1.setTime(1.2,1,2) ; e=d.applyFunc("2*x") ; e.copyStringInfoFrom(d) ; f1.setArray(e) ;
        ff1.appendFieldProfile(f1,mm1,0,da)
        ff1.write(fname,0)
        #
        vals,pfl=ff1.getFieldWithProfile(ON_CELLS,1,2,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(e,1e-14))
        vals,pfl=ff1.getFieldWithProfile(ON_CELLS,-1,-1,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(d,1e-14))
        #
        ff2=MEDFileFieldMultiTS.New(fname,f1.getName())
        self.assertEqual([(-1,-1,0.0), (1,2,1.2)],ff2.getTimeSteps())
        vals,pfl=ff2.getFieldWithProfile(ON_CELLS,1,2,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(e,1e-14))
        vals,pfl=ff2.getFieldWithProfile(ON_CELLS,-1,-1,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(d,1e-14))
        pass
    
    # idem testMEDField9 method except that here testing profile on nodes and not on cells.
    def testMEDField11(self):
        fname="Pyfile32.med"
        m1=MEDLoaderDataForTest.build2DMesh_1()
        m1.renumberCells([0,1,4,2,3,5],False)
        mm1=MEDFileUMesh.New() ; mm1.setCoords(m1.getCoords()) ; mm1.setMeshAtLevel(0,m1) ;
        mm1.write(fname,2)
        ff1=MEDFileField1TS.New()
        f1=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME) ; f1.setName("F1Node")
        d=DataArrayDouble.New() ; d.alloc(2*6,1) ; d.iota(7.); d.rearrange(2); d.setInfoOnComponent(0,"sigX [MPa]") ; d.setInfoOnComponent(1,"sigY [GPa]")
        f1.setArray(d) # note that f1 is NOT defined fully (no mesh !). It is not a bug of test it is too test that MEDFileField1TS.appendFieldProfile is NOT sensible of that.
        da=DataArrayInt.New(); da.setValues([1,2,4,5,7,8],6,1) ; da.setName("sup1Node")
        #
        ff1.setFieldProfile(f1,mm1,0,da)
        self.assertEqual(ff1.getNonEmptyLevels(),(-1, []))
        ff1.write(fname,0)
        #
        vals,pfl=ff1.getFieldWithProfile(ON_NODES,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(d,1e-14))
        ## #
        ff2=MEDFileField1TS.New(fname,f1.getName(),-1,-1)
        vals,pfl=ff2.getFieldWithProfile(ON_NODES,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(d,1e-14))
        pass

    def testMEDField12(self):
        fname="Pyfile33.med"
        m1=MEDLoaderDataForTest.build2DMesh_1()
        m1.renumberCells([0,1,4,2,3,5],False)
        mm1=MEDFileUMesh.New() ; mm1.setCoords(m1.getCoords()) ; mm1.setMeshAtLevel(0,m1) ;
        mm1.write(fname,2)
        ff1=MEDFileFieldMultiTS.New()
        f1=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME) ; f1.setName("F1Node")
        d=DataArrayDouble.New() ; d.alloc(2*6,1) ; d.iota(7.); d.rearrange(2); d.setInfoOnComponent(0,"sigX [MPa]") ; d.setInfoOnComponent(1,"sigY [GPa]")
        f1.setArray(d) # note that f1 is NOT defined fully (no mesh !). It is not a bug of test it is too test that MEDFileField1TS.appendFieldProfile is NOT sensible of that.
        da=DataArrayInt.New(); da.setValues([1,2,4,5,7,8],6,1) ; da.setName("sup1Node")
        #
        ff1.appendFieldProfile(f1,mm1,0,da)
        f1.setTime(1.2,1,2) ; e=d.applyFunc("2*x") ; e.copyStringInfoFrom(d) ; f1.setArray(e) ;
        ff1.appendFieldProfile(f1,mm1,0,da)
        ff1.write(fname,0)
        #
        vals,pfl=ff1.getFieldWithProfile(ON_NODES,1,2,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(e,1e-14))
        vals,pfl=ff1.getFieldWithProfile(ON_NODES,-1,-1,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(d,1e-14))
        #
        ff2=MEDFileFieldMultiTS.New(fname,f1.getName())
        vals,pfl=ff2.getFieldWithProfile(ON_NODES,1,2,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(e,1e-14))
        vals,pfl=ff2.getFieldWithProfile(ON_NODES,-1,-1,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(d,1e-14))
        pass

    def testMEDField13(self):
        fname="Pyfile34.med"
        m1=MEDLoaderDataForTest.build2DMesh_1()
        m1.renumberCells([0,1,4,2,3,5],False)
        tmp=m1.getName();
        m1=m1.buildPartOfMySelf(range(5),True) ; m1.setName(tmp) # suppression of last cell that is a polygon
        mm1=MEDFileUMesh.New() ; mm1.setCoords(m1.getCoords()) ; mm1.setMeshAtLevel(0,m1) ;
        mm1.write(fname,2)
        ff1=MEDFileField1TS.New()
        f1=MEDCouplingFieldDouble.New(ON_GAUSS_NE,ONE_TIME) ; f1.setName("F3Node")
        d=DataArrayDouble.New() ; d.alloc(2*11,1) ; d.iota(7.); d.rearrange(2); d.setInfoOnComponent(0,"sigX [MPa]") ; d.setInfoOnComponent(1,"sigY [GPa]")
        f1.setArray(d) # note that f1 is NOT defined fully (no mesh !). It is not a bug of test it is too test that MEDFileField1TS.appendFieldProfile is NOT sensible of that.
        da=DataArrayInt.New(); da.setValues([0,2,3],3,1) ; da.setName("sup1NodeElt")
        #
        ff1.setFieldProfile(f1,mm1,0,da)
        ff1.write(fname,0)
        #
        vals,pfl=ff1.getFieldWithProfile(ON_GAUSS_NE,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(d,1e-14))
        #
        ff2=MEDFileField1TS.New(fname,f1.getName(),-1,-1)
        vals,pfl=ff2.getFieldWithProfile(ON_GAUSS_NE,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(d,1e-14))
        pass

    def testMEDField14(self):
        fname="Pyfile35.med"
        m1=MEDLoaderDataForTest.build2DMesh_1()
        m1.renumberCells([0,1,4,2,3,5],False)
        tmp=m1.getName();
        m1=m1.buildPartOfMySelf(range(5),True) ; m1.setName(tmp) # suppression of last cell that is a polygon
        mm1=MEDFileUMesh.New() ; mm1.setCoords(m1.getCoords()) ; mm1.setMeshAtLevel(0,m1) ;
        mm1.write(fname,2)
        ff1=MEDFileFieldMultiTS.New()
        f1=MEDCouplingFieldDouble.New(ON_GAUSS_NE,ONE_TIME) ; f1.setName("F4Node")
        d=DataArrayDouble.New() ; d.alloc(2*11,1) ; d.iota(7.); d.rearrange(2); d.setInfoOnComponent(0,"sigX [MPa]") ; d.setInfoOnComponent(1,"sigY [GPa]")
        f1.setArray(d) # note that f1 is NOT defined fully (no mesh !). It is not a bug of test it is too test that MEDFileField1TS.appendFieldProfile is NOT sensible of that.
        da=DataArrayInt.New(); da.setValues([0,2,3],3,1) ; da.setName("sup1NodeElt")
        #
        ff1.appendFieldProfile(f1,mm1,0,da)
        f1.setTime(1.2,1,2) ; e=d.applyFunc("2*x") ; e.copyStringInfoFrom(d) ; f1.setArray(e) ;
        ff1.appendFieldProfile(f1,mm1,0,da)
        ff1.write(fname,0)
        #
        vals,pfl=ff1.getFieldWithProfile(ON_GAUSS_NE,-1,-1,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(d,1e-14))
        vals,pfl=ff1.getFieldWithProfile(ON_GAUSS_NE,1,2,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(e,1e-14))
        self.assertEqual([[3],[3]],ff1.getTypesOfFieldAvailable())
        #
        ff2=MEDFileFieldMultiTS.New(fname,f1.getName())
        vals,pfl=ff1.getFieldWithProfile(ON_GAUSS_NE,-1,-1,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(d,1e-14))
        vals,pfl=ff1.getFieldWithProfile(ON_GAUSS_NE,1,2,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))
        self.assertTrue(vals.isEqual(e,1e-14))
        pass
    # Tricky test of the case of in a MED file containing a Field on GAUSS_NE is lying on a profile that is reality represents all the geom entities of a level.
    # By default when using setFieldProfile method such profile is not created because it is not useful ! So here a trick is used to force MEDLoader to do that
    # for the necessity of the test ! The idea is too create artificially a mesh having one more fictious cell per type and to roll back right after !
    def testMEDField15(self):
        fname="Pyfile36.med"
        m0=MEDLoaderDataForTest.build2DMesh_1()
        m0.renumberCells([0,1,4,2,3,5],False)
        tmp=m0.getName();
        m1=m0.buildPartOfMySelf([0,1,1,2,3,3,4,4],True) ; m1.setName(tmp) # suppression of last cell that is a polygon and creation of one more cell per type
        mm1=MEDFileUMesh.New() ; mm1.setCoords(m1.getCoords()) ; mm1.setMeshAtLevel(0,m1) ;
        ff1=MEDFileField1TS.New()
        f1=MEDCouplingFieldDouble.New(ON_GAUSS_NE,ONE_TIME) ; f1.setName("F4Node")
        d=DataArrayDouble.New() ; d.alloc(2*20,1) ; d.iota(7.); d.rearrange(2); d.setInfoOnComponent(0,"sigX [MPa]") ; d.setInfoOnComponent(1,"sigY [GPa]")
        f1.setArray(d) # note that f1 is NOT defined fully (no mesh !). It is not a bug of test it is too test that MEDFileField1TS.appendFieldProfile is NOT sensible of that.
        da=DataArrayInt.New(); da.setValues([0,1,3,4,6],5,1) ; da.setName("sup1NodeElt")
        #
        ff1.setFieldProfile(f1,mm1,0,da)
        m1=m0.buildPartOfMySelf(range(5),True) ; m1.setName(tmp) ; mm1.setMeshAtLevel(0,m1) ;
        mm1.write(fname,2)
        ff1.write(fname,0)
        f1=ff1.getFieldOnMeshAtLevel(ON_GAUSS_NE,m1,0)
        f2,p1=ff1.getFieldWithProfile(ON_GAUSS_NE,0,mm1) ; f2.setName("")
        self.assertTrue(p1.isIdentity())
        self.assertEqual(5,p1.getNumberOfTuples())
        self.assertTrue(f1.getArray().isEqual(f2,1e-12))
        pass
    # Test for getFieldAtTopLevel method
    def testMEDField16(self):
        fname="Pyfile37.med"
        f1=MEDLoaderDataForTest.buildVecFieldOnCells_1();
        m1=f1.getMesh()
        mm1=MEDFileUMesh.New()
        mm1.setCoords(m1.getCoords())
        mm1.setMeshAtLevel(0,m1)
        mm1.setName(m1.getName())
        ff1=MEDFileField1TS.New()
        ff1.setFieldNoProfileSBT(f1)
        m2=m1.buildDescendingConnectivity()[0]
        m2.sortCellsInMEDFileFrmt()
        m2.setName(m1.getName())
        mm1.setMeshAtLevel(-1,m2)
        mm1.write(fname,2)
        f2=m2.getMeasureField(True)
        dd=DataArrayDouble.New()
        dd.alloc(f2.getArray().getNumberOfTuples(),3)
        dd[:,0]=f2.getArray()
        dd[:,1]=2*f2.getArray()
        dd[:,2]=3*f2.getArray()
        f2=f2.buildNewTimeReprFromThis(ONE_TIME,False)
        f2.setArray(dd)
        f2.copyTinyStringsFrom(f1)
        f2.copyTinyAttrFrom(f1)
        ff1.setFieldNoProfileSBT(f2)
        ff1.write(fname,0)
        # Reading Pyfile37.med
        ff2=MEDFileField1TS.New(fname,f2.getName(),0,1)
        f1bis=ff2.getFieldAtLevel(ON_CELLS,0)
        self.assertTrue(f1.isEqual(f1bis,1e-12,1e-12))
        f1bis=ff2.getFieldAtLevel(ON_CELLS,-1)
        self.assertTrue(f2.isEqual(f1bis,1e-12,1e-12))
        f1bis=ff2.getFieldAtTopLevel(ON_CELLS)
        self.assertTrue(f1.isEqual(f1bis,1e-12,1e-12))
        # More complex
        fname="Pyfile38.med"
        mm1.write(fname,2)
        ff1=MEDFileField1TS.New()
        ff1.setFieldNoProfileSBT(f2)
        ff1.write(fname,0)
        ff2=MEDFileField1TS.New(fname,f2.getName(),0,1)
        f1bis=ff2.getFieldAtTopLevel(ON_CELLS)
        self.assertTrue(f2.isEqual(f1bis,1e-12,1e-12))
        pass

    # Non regression test to check that globals are correctly appended on MEDFileFields::setFieldAtPos
    def testMEDField17(self):
        fname="Pyfile39.med"
        m1=MEDLoaderDataForTest.build2DMesh_1()
        m1.renumberCells([0,1,4,2,3,5],False)
        mm1=MEDFileUMesh.New() ; mm1.setCoords(m1.getCoords()) ; mm1.setMeshAtLevel(0,m1) ; mm1.setName(m1.getName())
        mm1.write(fname,2)
        ffs=MEDFileFields.New()
        ff1=MEDFileFieldMultiTS.New()
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f1.setName("F2")
        d=DataArrayDouble.New() ; d.alloc(2*4,1) ; d.iota(7.); d.rearrange(2); d.setInfoOnComponent(0,"sigX [MPa]") ; d.setInfoOnComponent(1,"sigY [GPa]")
        f1.setArray(d) # note that f1 is NOT defined fully (no mesh !). It is not a bug of test it is too test that MEDFileField1TS.appendFieldProfile is NOT sensible of that.
        da=DataArrayInt.New(); da.setValues([0,1,2,4],4,1) ; da.setName("sup2")
        #
        ff1.appendFieldProfile(f1,mm1,0,da)
        f1.setTime(1.2,1,2) ; e=d.applyFunc("2*x") ; e.copyStringInfoFrom(d) ; f1.setArray(e) ;
        ff1.appendFieldProfile(f1,mm1,0,da)
        ffs.resize(1)
        ffs.setFieldAtPos(0,ff1)
        ffs.write(fname,0)
        #
        ffsr=MEDFileFields.New(fname)
        ff3=ffsr.getFieldAtPos(0)
        f4=ff3.getFieldAtTopLevel(ON_CELLS,1,2)
        self.assertTrue(f4.getArray().isEqual(f1.getArray(),1e-12))
        pass

    # Non regression test to check that globals are correctly appended on MEDFileFields::setFieldAtPos
    def testMEDField18(self):
        fname="Pyfile40.med"
        m1=MEDLoaderDataForTest.build2DMesh_1()
        m1.renumberCells([0,1,4,2,3,5],False)
        mm1=MEDFileUMesh.New() ; mm1.setCoords(m1.getCoords()) ; mm1.setMeshAtLevel(0,m1) ; mm1.setName(m1.getName())
        mm1.write(fname,2)
        ffs=MEDFileFields.New()
        ff1=MEDFileFieldMultiTS.New()
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f1.setName("F2")
        d=DataArrayDouble.New() ; d.alloc(2*4,1) ; d.iota(7.); d.rearrange(2); d.setInfoOnComponent(0,"sigX [MPa]") ; d.setInfoOnComponent(1,"sigY [GPa]")
        f1.setArray(d) # note that f1 is NOT defined fully (no mesh !). It is not a bug of test it is too test that MEDFileField1TS.appendFieldProfile is NOT sensible of that.
        da=DataArrayInt.New(); da.setValues([0,1,2,4],4,1) ; da.setName("sup2")
        #
        ff1.appendFieldProfile(f1,mm1,0,da)
        f1.setTime(1.2,1,2) ; e=d.applyFunc("2*x") ; e.copyStringInfoFrom(d) ; f1.setArray(e) ;
        ff1.appendFieldProfile(f1,mm1,0,da)
        ffs.pushField(ff1)
        ffs.write(fname,0)
        #
        ffsr=MEDFileFields.New(fname)
        ff3=ffsr.getFieldAtPos(0)
        f4=ff3.getFieldAtTopLevel(ON_CELLS,1,2)
        self.assertTrue(f4.getArray().isEqual(f1.getArray(),1e-12))
        pass

    def testMEDFieldBug1(self):
        fname="Pyfile13.med"
        d=MEDFileData.New(fname)
        self.assertEqual(('GP_MyFirstFieldOnGaussPoint0', 'GP_MyFirstFieldOnGaussPoint1', 'GP_MyFirstFieldOnGaussPoint2'),d.getFields().getFieldAtPos(0).getLocs())
        pass

    def testMEDMesh8(self):
        m=MEDLoaderDataForTest.build1DMesh_1()
        m.convertQuadraticCellsToLinear()
        mm=MEDFileUMesh.New()
        mm.setMeshAtLevel(0,m)
        g1=DataArrayInt.New() ; g1.setValues([0,2],2,1) ; g1.setName("g1")
        g2=DataArrayInt.New() ; g2.setValues([1,3],2,1) ; g2.setName("g2")
        g3=DataArrayInt.New() ; g3.setValues([1,2,3],3,1) ; g3.setName("g3")
        mm.setGroupsAtLevel(0,[g1,g2],False)
        self.assertEqual(('g1','g2'),mm.getGroupsNames())
        self.assertEqual(('Family_2','Family_3'),mm.getFamiliesNames())
        self.assertEqual(('Family_2',),mm.getFamiliesOnGroup('g1'))
        self.assertEqual(('Family_3',),mm.getFamiliesOnGroup('g2'))
        mm.assignFamilyNameWithGroupName()
        self.assertEqual(('g1','g2'),mm.getGroupsNames())
        self.assertEqual(('g1','g2'),mm.getFamiliesNames())
        self.assertEqual(('g1',),mm.getFamiliesOnGroup('g1'))
        self.assertEqual(('g2',),mm.getFamiliesOnGroup('g2'))
        #
        mm=MEDFileUMesh.New()
        mm.setMeshAtLevel(0,m)
        mm.setGroupsAtLevel(0,[g1,g2,g3],False)
        self.assertEqual(('g1','g2','g3'),mm.getGroupsNames())
        self.assertEqual(('Family_2', 'Family_4', 'Family_5'),mm.getFamiliesNames())
        self.assertEqual(('Family_2', 'Family_4'),mm.getFamiliesOnGroup('g1'))
        self.assertEqual(('Family_5',),mm.getFamiliesOnGroup('g2'))
        self.assertEqual(('Family_4','Family_5',),mm.getFamiliesOnGroup('g3'))
        mm.assignFamilyNameWithGroupName() # here it does nothing because no such group-family bijection found
        self.assertEqual(('g1','g2','g3'),mm.getGroupsNames())
        self.assertEqual(('Family_2', 'Family_4', 'Family_5'),mm.getFamiliesNames())
        self.assertEqual(('Family_2', 'Family_4'),mm.getFamiliesOnGroup('g1'))
        self.assertEqual(('Family_5',),mm.getFamiliesOnGroup('g2'))
        self.assertEqual(('Family_4','Family_5',),mm.getFamiliesOnGroup('g3'))
        mm.changeFamilyId(5,6)
        g=mm.getGroupArr(0,"g3")
        self.assertTrue(g.isEqual(g3));
        g=mm.getGroupArr(0,"g2")
        self.assertTrue(g.isEqual(g2));
        g=mm.getGroupArr(0,"g1")
        self.assertTrue(g.isEqual(g1));
        pass
    
    # bug detected by gauthier
    def testMEDLoaderMEDLoaderNSReadFieldDoubleDataInMedFile(self):
        fname="Pyfile41.med"
        f1=MEDLoaderDataForTest.buildVecFieldOnCells_1();
        m1=f1.getMesh()
        mm1=MEDFileUMesh.New()
        mm1.setCoords(m1.getCoords())
        mm1.setMeshAtLevel(0,m1)
        mm1.write(fname,2)
        ff1=MEDFileField1TS.New()
        ff1.setFieldNoProfileSBT(f1)
        ff1.write(fname,0)
        # writing mesh1 and field1, now creation of mesh2 and field2
        f2=f1.deepCpy()
        m2=f2.getMesh()
        m2.translate([0.5,0.6,0.7])
        m2.setName("3DSurfMesh_2")
        f2.getArray()[:]*=2.
        f2.setName("VectorFieldOnCells2")
        mm2=MEDFileUMesh.New()
        mm2.setCoords(m2.getCoords())
        mm2.setMeshAtLevel(0,m2)
        mm2.write(fname,0)
        ff2=MEDFileField1TS.New()
        ff2.setFieldNoProfileSBT(f2)
        ff2.write(fname,0)
        #
        f3=MEDLoader.ReadFieldCell(fname,"3DSurfMesh_1",0,"VectorFieldOnCells",0,1)
        self.assertTrue(f3.isEqual(f1,1e-12,1e-12))
        f4=MEDLoader.ReadFieldCell(fname,"3DSurfMesh_2",0,"VectorFieldOnCells2",0,1)
        self.assertTrue(f4.isEqual(f2,1e-12,1e-12))
        pass

    def testMEDLoaderMultiLevelCellField1(self):
        fname="Pyfile42.med"
        m2,m1,m0,f2,f1,f0,p,n2,n1,n0,fns,fids,grpns,famIdsPerGrp=MEDLoaderDataForTest.buildMultiLevelMesh_1()
        m=MEDFileUMesh.New()
        m.setCoords(m2.getCoords())
        m.setMeshAtLevel(0,m2)
        m.setMeshAtLevel(-1,m1)
        m.setMeshAtLevel(-2,m0)
        m.write(fname,2)
        #
        FieldName1="Field1"
        compNames1=["comp1","comp2","comp3"]
        ff1=MEDFileField1TS.New()
        da2=DataArrayDouble.New()
        da2.alloc(m2.getNumberOfCells()*len(compNames1),1)
        da2.iota(7.)
        da2.rearrange(len(compNames1))
        da2.setInfoOnComponents(compNames1)
        f2=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f2.setName(FieldName1) ; f2.setArray(da2) ; f2.setMesh(m2) ; f2.checkCoherency()
        ff1.setFieldNoProfileSBT(f2)
        self.assertEqual(ff1.getNonEmptyLevels(),(2, [0]))
        da0=DataArrayDouble.New()
        da0.alloc(m0.getNumberOfCells()*len(compNames1),1)
        da0.iota(190.)
        da0.rearrange(len(compNames1))
        da0.setInfoOnComponents(compNames1)
        f0=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f0.setName(FieldName1) ; f0.setArray(da0) ; f0.setMesh(m0) ; f0.checkCoherency()
        ff1.setFieldNoProfileSBT(f0)
        self.assertEqual(ff1.getNonEmptyLevels(),(2, [0,-2]))
        da1=DataArrayDouble.New()
        da1.alloc(m1.getNumberOfCells()*len(compNames1),1)
        da1.iota(90.)
        da1.rearrange(len(compNames1))
        da1.setInfoOnComponents(compNames1)
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f1.setName(FieldName1) ; f1.setArray(da1) ; f1.setMesh(m1) ; f1.checkCoherency()
        ff1.setFieldNoProfileSBT(f1)
        self.assertEqual(ff1.getNonEmptyLevels(),(2, [0,-1,-2]))
        #
        ff1.write(fname,0)
        #
        FieldName2="Field2"
        compNames2=["comp11","comp22"]
        ff2=MEDFileField1TS.New()
        da0=DataArrayDouble.New()
        da0.alloc(m0.getNumberOfCells()*2,1)
        da0.iota(-190.)
        da0.rearrange(2)
        da0.setInfoOnComponents(compNames2)
        f0=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f0.setName(FieldName2) ; f0.setArray(da0) ; f0.setMesh(m0) ; f0.checkCoherency()
        ff2.setFieldNoProfileSBT(f0)
        self.assertEqual(ff2.getNonEmptyLevels(),(0, [0]))
        da1=DataArrayDouble.New()
        da1.alloc(m1.getNumberOfCells()*len(compNames2),1)
        da1.iota(-90.)
        da1.rearrange(len(compNames2))
        da1.setInfoOnComponents(compNames2)
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f1.setName(FieldName2) ; f1.setArray(da1) ; f1.setMesh(m1) ; f1.checkCoherency()
        ff2.setFieldNoProfileSBT(f1)
        self.assertEqual(ff2.getNonEmptyLevels(),(1, [0,-1]))
        #
        ff2.write(fname,0)
        #
        ff1=MEDFileField1TS.New(fname,FieldName1,-1,-1)
        self.assertEqual(ff1.getNonEmptyLevels(),(2, [0,-1,-2]))
        self.assertEqual(ff1.getFieldSplitedByType(),[(0, [(0, (0, 4), '', '')]), (1, [(0, (4, 84), '', '')]), (3, [(0, (84, 148), '', '')]), (4, [(0, (148, 212), '', '')])])
        ff2=MEDFileField1TS.New(fname,FieldName2,-1,-1)
        self.assertEqual(ff2.getNonEmptyLevels(),(1, [0,-1]))
        self.assertEqual(ff2.getFieldSplitedByType(),[(0, [(0, (0, 4), '', '')]), (1, [(0, (4, 84), '', '')])])
        pass

    def testFieldOnPflRetrieveOnMdimRelMax1(self):
        fname="Pyfile43.med"
        m2,m1,m0,f2,f1,f0,p,n2,n1,n0,fns,fids,grpns,famIdsPerGrp=MEDLoaderDataForTest.buildMultiLevelMesh_1()
        m=MEDFileUMesh.New()
        m.setMeshAtLevel(0,m2)
        m.setMeshAtLevel(-1,m1)
        m.setMeshAtLevel(-2,m0)
        f=MEDFileField1TS.New()
        ff=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME)
        ff.setName("NodeFieldPfl")
        arr=DataArrayDouble.New() ; arr.setValues([1.,10.,100.,2.,20.,200.],2,3)
        ff.setArray(arr)
        pfl=DataArrayInt.New() ; pfl.setValues([2,3],2,1) ; pfl.setName("PflNode")
        f.setFieldProfile(ff,m,-2,pfl)
        tes0=f.getFieldOnMeshAtLevel(ON_NODES,-1,m)
        self.assertEqual(ON_NODES,tes0.getTypeOfField())
        self.assertEqual(1,tes0.getMesh().getMeshDimension())
        self.assertEqual(1,tes0.getMesh().getNumberOfCells())
        self.assertEqual(2,tes0.getMesh().getNumberOfNodes())
        self.assertEqual([1,0,1],tes0.getMesh().getNodalConnectivity().getValues())
        self.assertEqual([0,3],tes0.getMesh().getNodalConnectivityIndex().getValues())
        self.assertEqual(2,tes0.getArray().getNumberOfTuples())
        self.assertEqual(3,tes0.getArray().getNumberOfComponents())
        expected1=[1.,10.,100.,2.,20.,200.]
        nodeCoordsWithValue1=[10.,2.5,0.]
        nodeCoordsWithValue2=[10.,3.75,0.]
        for i in xrange(3):
            self.assertAlmostEqual(nodeCoordsWithValue1[i],tes0.getMesh().getCoordinatesOfNode(0)[i],13);
            self.assertAlmostEqual(nodeCoordsWithValue2[i],tes0.getMesh().getCoordinatesOfNode(1)[i],13);
            pass
        for i in xrange(6):
            self.assertAlmostEqual(expected1[i],tes0.getArray().getIJ(0,i),13);
            pass
        del tes0
        #
        tes1=f.getFieldOnMeshAtLevel(ON_NODES,1,m)
        self.assertEqual(ON_CELLS,tes1.getTypeOfField())# it is not a bug even if ON_NODES has been sepecified
        self.assertEqual(0,tes1.getMesh().getMeshDimension())
        self.assertEqual(2,tes1.getMesh().getNumberOfCells())
        self.assertEqual(135,tes1.getMesh().getNumberOfNodes())
        self.assertEqual([0,2,0,3],tes1.getMesh().getNodalConnectivity().getValues())
        self.assertEqual([0,2,4],tes1.getMesh().getNodalConnectivityIndex().getValues())
        self.assertEqual(2,tes1.getArray().getNumberOfTuples())
        self.assertEqual(3,tes1.getArray().getNumberOfComponents())
        for i in xrange(6):
            self.assertAlmostEqual(expected1[i],tes1.getArray().getIJ(0,i),13);
            pass
        m.write(fname,2)
        f.write(fname,0)
        #
        pfl=DataArrayInt.New() ; pfl.setValues([3,2],2,1) ; pfl.setName("PflNode")
        f=MEDFileField1TS.New()
        f.setFieldProfile(ff,m,-2,pfl)
        tes2=f.getFieldOnMeshAtLevel(ON_NODES,-1,m)
        self.assertEqual(ON_NODES,tes2.getTypeOfField())
        self.assertEqual(1,tes2.getMesh().getMeshDimension())
        self.assertEqual(1,tes2.getMesh().getNumberOfCells())
        self.assertEqual(2,tes2.getMesh().getNumberOfNodes())
        self.assertEqual([1,0,1],tes2.getMesh().getNodalConnectivity().getValues())
        self.assertEqual([0,3],tes2.getMesh().getNodalConnectivityIndex().getValues())
        self.assertEqual(2,tes2.getArray().getNumberOfTuples())
        self.assertEqual(3,tes2.getArray().getNumberOfComponents())
        expected2=[2.,20.,200.,1.,10.,100.]
        for i in xrange(3):
            self.assertAlmostEqual(nodeCoordsWithValue1[i],tes2.getMesh().getCoordinatesOfNode(0)[i],13);
            self.assertAlmostEqual(nodeCoordsWithValue2[i],tes2.getMesh().getCoordinatesOfNode(1)[i],13);
            pass
        for i in xrange(6):
            self.assertAlmostEqual(expected2[i],tes2.getArray().getIJ(0,i),13);#compare tes2 and tes3
            pass
        #
        tes3=f.getFieldOnMeshAtLevel(ON_NODES,1,m)
        self.assertEqual(ON_CELLS,tes3.getTypeOfField())# it is not a bug even if ON_NODES has been sepecified
        self.assertEqual(0,tes3.getMesh().getMeshDimension())
        self.assertEqual(2,tes3.getMesh().getNumberOfCells())
        self.assertEqual(135,tes3.getMesh().getNumberOfNodes())
        self.assertEqual([0,3,0,2],tes3.getMesh().getNodalConnectivity().getValues())
        self.assertEqual([0,2,4],tes3.getMesh().getNodalConnectivityIndex().getValues())
        self.assertEqual(2,tes3.getArray().getNumberOfTuples())
        self.assertEqual(3,tes3.getArray().getNumberOfComponents())
        for i in xrange(6):
            self.assertAlmostEqual(expected1[i],tes3.getArray().getIJ(0,i),13);
            pass
        pass

    def testDuplicateNodesOnM1Group1(self):
        fname="Pyfile44.med"
        m=MEDCouplingCMesh.New()
        m.setCoordsAt(0,DataArrayDouble.New([0.,1.1,2.3,3.6,5.,6.5]))
        m.setCoordsAt(1,DataArrayDouble.New([0.,1.1,2.3,3.6,5.]))
        m=m.buildUnstructured() ; m.setName("AnthonyDuplicate")
        m.getCoords().setInfoOnComponents(["X [km]","Z [mm]"])
        m2=m.buildDescendingConnectivity()[0][[8,11,14,20,21,22,23,24,25,26,31,32,33,34,35,36,37]]
        m2.setName(m.getName())
        grp=DataArrayInt.New([4,6,8]) ; grp.setName("Grp")
        grp2=DataArrayInt.New([9,16]) ; grp2.setName("Grp2")
        mm=MEDFileUMesh.New()
        mm.setMeshAtLevel(0,m)
        mm.setMeshAtLevel(-1,m2)
        mm.setGroupsAtLevel(-1,[grp,grp2])
        grpNode=DataArrayInt.New([4,21,23]) ; grpNode.setName("GrpNode")
        mm.setGroupsAtLevel(1,[grpNode])
        ref0=[4,15,14,20,21,4,16,15,21,22,4,17,16,22,23]
        ref1=[4,9,8,14,15,4,10,9,15,16,4,11,10,16,17]
        ref2=[4,9,8,14,30,4,10,9,30,31,4,11,10,31,32]
        #
        self.assertEqual(30,mm.getNumberOfNodes())
        self.assertEqual(ref0,mm.getMeshAtLevel(0)[[12,13,14]].getNodalConnectivity().getValues())
        self.assertEqual(ref1,mm.getMeshAtLevel(0)[[7,8,9]].getNodalConnectivity().getValues())
        #
        nodes,cells,cells2=mm.duplicateNodesOnM1Group("Grp")
        self.assertEqual([15,16,17],nodes.getValues());
        self.assertEqual([7,8,9],cells.getValues());
        self.assertEqual([12,13,14],cells2.getValues());
        self.assertEqual(33,mm.getNumberOfNodes())
        self.assertEqual([4,6,8],mm.getGroupArr(-1,"Grp").getValues())
        self.assertEqual([9,16],mm.getGroupArr(-1,"Grp2").getValues())
        self.assertEqual([4,21,23],mm.getGroupArr(1,"GrpNode").getValues())
        self.assertEqual([17,18,19],mm.getGroupArr(-1,"Grp_dup").getValues())
        self.assertEqual(ref0,mm.getMeshAtLevel(0)[[12,13,14]].getNodalConnectivity().getValues())#cells 7,8,9 and 12,13,14 are lying on "Grp" but only 7,8 and 9 are renumbered
        self.assertEqual(ref2,mm.getMeshAtLevel(0)[[7,8,9]].getNodalConnectivity().getValues())#
        self.assertRaises(InterpKernelException,mm.getGroup(-1,"Grp_dup").checkGeoEquivalWith,mm.getGroup(-1,"Grp"),2,1e-12);# Grp_dup and Grp are not equal considering connectivity only
        mm.getGroup(-1,"Grp_dup").checkGeoEquivalWith(mm.getGroup(-1,"Grp"),12,1e-12)# Grp_dup and Grp are equal considering connectivity and coordinates
        refValues=DataArrayDouble.New([1.21,1.32,1.43,1.54,1.65,1.32,1.44,1.56,1.68,1.8,1.43,1.56,1.69,1.82,1.95,1.54,1.68,1.82,1.96,2.1])
        valsToTest=mm.getMeshAtLevel(0).getMeasureField(True).getArray() ; delta=(valsToTest-refValues) ; delta.abs()
        self.assertTrue(delta.getMaxValue()[0]<1e-12)
        #
        mm.getCoords()[-len(nodes):]+=[0.,-0.3]
        self.assertRaises(InterpKernelException,mm.getGroup(-1,"Grp_dup").checkGeoEquivalWith,mm.getGroup(-1,"Grp"),12,1e-12);
        refValues2=refValues[:] ; refValues2[7:10]=[1.365,1.26,1.35]
        valsToTest=mm.getMeshAtLevel(0).getMeasureField(True).getArray() ; delta=(valsToTest-refValues2) ; delta.abs()
        self.assertTrue(delta.getMaxValue()[0]<1e-12)
        mm.write(fname,2)
        pass

    def testDuplicateNodesOnM1Group2(self):
        fname="Pyfile45.med"
        m=MEDCouplingCMesh.New()
        m.setCoordsAt(0,DataArrayDouble.New([0.,1.1,2.3,3.6,5.,6.5]))
        m.setCoordsAt(1,DataArrayDouble.New([0.,1.1,2.3,3.6,5.]))
        m=m.buildUnstructured() ; m.setName("AnthonyDuplicate")
        m.getCoords().setInfoOnComponents(["X [km]","Z [mm]"])
        m2=m.buildDescendingConnectivity()[0][[8,11,14,20,21,22,23,24,25,26,31,32,33,34,35,36,37]]
        m2.setName(m.getName())
        grp=DataArrayInt.New([4,6]) ; grp.setName("Grp")
        grp2=DataArrayInt.New([9,16]) ; grp2.setName("Grp2")
        mm=MEDFileUMesh.New()
        mm.setMeshAtLevel(0,m)
        mm.setMeshAtLevel(-1,m2)
        mm.setGroupsAtLevel(-1,[grp,grp2])
        grpNode=DataArrayInt.New([4,21,23]) ; grpNode.setName("GrpNode")
        mm.setGroupsAtLevel(1,[grpNode])
        ref0=[4,15,14,20,21,4,16,15,21,22,4,17,16,22,23]
        ref1=[4,9,8,14,15,4,10,9,15,16]
        ref2=[4,9,8,14,30,4,10,9,30,16]
        #
        self.assertEqual(30,mm.getNumberOfNodes())
        self.assertEqual(ref0,mm.getMeshAtLevel(0)[[12,13,14]].getNodalConnectivity().getValues())
        self.assertEqual(ref1,mm.getMeshAtLevel(0)[[7,8]].getNodalConnectivity().getValues())
        #
        nodes,cells,cells2=mm.duplicateNodesOnM1Group("Grp")
        self.assertEqual([15],nodes.getValues());
        self.assertEqual([7,8],cells.getValues());
        self.assertEqual([12,13],cells2.getValues());
        self.assertEqual(31,mm.getNumberOfNodes())
        self.assertEqual([4,6],mm.getGroupArr(-1,"Grp").getValues())
        self.assertEqual([9,16],mm.getGroupArr(-1,"Grp2").getValues())
        self.assertEqual([4,21,23],mm.getGroupArr(1,"GrpNode").getValues())
        self.assertEqual([17,18],mm.getGroupArr(-1,"Grp_dup").getValues())
        self.assertEqual(ref0,mm.getMeshAtLevel(0)[[12,13,14]].getNodalConnectivity().getValues())#cells 7,8,9 and 12,13,14 are lying on "Grp" but only 7,8 and 9 are renumbered
        self.assertEqual(ref2,mm.getMeshAtLevel(0)[[7,8]].getNodalConnectivity().getValues())#
        self.assertRaises(InterpKernelException,mm.getGroup(-1,"Grp_dup").checkGeoEquivalWith,mm.getGroup(-1,"Grp"),2,1e-12);# Grp_dup and Grp are not equal considering connectivity only
        mm.getGroup(-1,"Grp_dup").checkGeoEquivalWith(mm.getGroup(-1,"Grp"),12,1e-12)# Grp_dup and Grp are equal considering connectivity and coordinates
        refValues=DataArrayDouble.New([1.21,1.32,1.43,1.54,1.65,1.32,1.44,1.56,1.68,1.8,1.43,1.56,1.69,1.82,1.95,1.54,1.68,1.82,1.96,2.1])
        valsToTest=mm.getMeshAtLevel(0).getMeasureField(True).getArray() ; delta=(valsToTest-refValues) ; delta.abs()
        self.assertTrue(delta.getMaxValue()[0]<1e-12)
        #
        mm.getCoords()[-len(nodes):]+=[0.,-0.3]
        self.assertRaises(InterpKernelException,mm.getGroup(-1,"Grp_dup").checkGeoEquivalWith,mm.getGroup(-1,"Grp"),12,1e-12);
        refValues2=refValues[:] ; refValues2[7:9]=[1.365,1.47]
        valsToTest=mm.getMeshAtLevel(0).getMeasureField(True).getArray() ; delta=(valsToTest-refValues2) ; delta.abs()
        self.assertTrue(delta.getMaxValue()[0]<1e-12)
        mm.write(fname,2)       
        pass

    def testBasicConstructors(self):
        fname="Pyfile18.med"
        m=MEDFileMesh(fname)
        m=MEDFileMesh(fname,"ExampleOfMultiDimW",-1,-1)
        m=MEDFileMesh(fname)
        m=MEDFileUMesh(fname,"ExampleOfMultiDimW",-1,-1)
        m=MEDFileUMesh(fname)
        m=MEDFileUMesh()
        self.testMEDMesh6()
        m=MEDFileCMesh("MEDFileMesh5.med")
        m=MEDFileCMesh("MEDFileMesh5.med","myFirstCartMesh",-1,-1)
        m=MEDFileCMesh()
        m=MEDFileMeshMultiTS()
        m=MEDFileMeshMultiTS(fname)
        m=MEDFileMeshMultiTS(fname,"ExampleOfMultiDimW")
        m=MEDFileMeshes()
        m=MEDFileMeshes(fname)
        m=MEDFileField1TS()
        m=MEDFileField1TS(fname,"FieldOnFacesShuffle",2,7)
        m=MEDFileFieldMultiTS()
        m=MEDFileFieldMultiTS(fname,"FieldOnFacesShuffle")
        m=MEDFileFields()
        m=MEDFileFields(fname)
        m=MEDFileData()
        m=MEDFileData(fname)
        #
        m=DataArrayInt() ; m=DataArrayInt(5,2) ; m=DataArrayInt([6,5,4,3,2,1],3,2)
        m=DataArrayDouble() ; m=DataArrayDouble(5,2) ; m=DataArrayDouble([6,5,4,3,2,1],3,2)
        m=MEDCouplingUMesh("jjj",2) ; m=MEDCouplingUMesh()
        m=MEDCouplingCMesh()
        m=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
        m=MEDCouplingFieldTemplate(ON_NODES)
        m=MEDCouplingMultiFields([])
        m=MEDCouplingFieldOverTime([])
        pass

    # This is a non regression test. When a field lies partially on a mesh but fully on one of its geometric type.
    def testBugSemiPartialField(self):
        fname="Pyfile46.med"
        m=MEDLoaderDataForTest.build2DMesh_3()
        m=m[:10] ; m.setName("mesh")
        f=m.getMeasureField(ON_CELLS)
        f=f.buildNewTimeReprFromThis(ONE_TIME,False)
        f.setTime(5.5,3,4)
        f.setName("SemiPartialField")
        #
        f1=f[:6] ; f1.getMesh().setName(m.getName())
        f2=f[6:] ; f2.getMesh().setName(m.getName())
        #
        mm=MEDFileUMesh.New()
        mm.setMeshAtLevel(0,m)
        ff=MEDFileField1TS.New()
        ff.setFieldProfile(f1,mm,0,DataArrayInt.Range(0,6,1)) # no name on profile -> normally it is an error but in this special case
        mm.write(fname,2)
        ff.write(fname,0)
        #
        ff2=MEDFileField1TS.New(fname,f.getName(),f.getTime()[1],f.getTime()[2])
        fread=ff2.getFieldOnMeshAtLevel(ON_CELLS,0,mm)
        fread2=ff2.getFieldAtLevel(ON_CELLS,0)
        #
        fread.checkCoherency()
        fread2.checkCoherency()
        self.assertTrue(fread.isEqual(f1,1e-12,1e-12))
        self.assertTrue(fread2.isEqual(f1,1e-12,1e-12))
        pass

    def testUnPolyze1(self):
        fname="Pyfile47.med"
        mm=MEDLoaderDataForTest.buildMLMeshUnPolyze(self)
        ref=[13,14,14,12,12,12,12,12,12,12,12,13,12,14,14,13,15,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12]
        self.assertEqual(ref,mm.getFamilyFieldAtLevel(1).getValues())
        self.assertEqual(mm.unPolyze()[:3],(True,[(3,2,0),(4,3,2),(5,4,5),(14,2,9),(16,3,11),(31,2,14)],[(3,3,0),(4,3,3),(5,3,6),(14,3,9),(16,3,12),(18,1,15)]))
        mm.write(fname,2)
        self.assertEqual(mm.getGroupArr(0,"grp0_L0").getValues(),[0,1,2,6])
        self.assertEqual(mm.getGroupArr(0,"grp1_L0").getValues(),[1,3,4,5,6])
        self.assertEqual(mm.getGroupArr(-1,"grp0_LM1").getValues(),[1,2,3,4,5])
        self.assertEqual(mm.getGroupArr(-1,"grp1_LM1").getValues(),[3,4,5,6])
        self.assertEqual(mm.getGroupArr(-1,"grp2_LM1").getValues(),[2,6,7,8])
        self.assertEqual(mm.getGroupArr(1,"grp0_Node").getValues(),[0,11,15,16])
        self.assertEqual(mm.getGroupArr(1,"grp1_Node").getValues(),[1,2,13,14,16])
        self.assertEqual(mm.getFamilyFieldAtLevel(1).getValues(),ref)
        # to test
        mm.setRenumFieldArr(0,None)
        mm.setFamilyFieldArr(-1,None)
        pass

    def testUnPolyze2(self):
        fname="Pyfile48.med"
        mfd=MEDFileData.New()
        mm=MEDLoaderDataForTest.buildMLMeshUnPolyze(self)
        meshes=MEDFileMeshes.New()
        meshes.pushMesh(mm)
        mfd.setMeshes(meshes)
        fields=MEDFileFields.New()
        mfd.setFields(fields)
        ff=MEDFileFieldMultiTS.New()
        fields.pushField(ff)
        #
        f0_0=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME) ; f0_0.setName("f0")
        f0_0.setTime(9.5,3,4)
        da=DataArrayDouble.New(38*2) ; da.iota(6.) ; da.rearrange(2) ; da.setInfoOnComponents(["Power [MW]","Density [kg/m^3]"])
        f0_0.setArray(da)
        f0_0.setMesh(mm.getMeshAtLevel(0))
        ff.appendFieldNoProfileSBT(f0_0)
        ff0=ff.getTimeStepAtPos(0)
        f0_1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f0_1.setName("f0")
        f0_1.setTime(9.5,3,4)
        pfl=DataArrayInt.New([1,4,5,6]) ; pfl.setName("pfltest")
        f0_1.setMesh(mm.getMeshAtLevel(0)[pfl])
        da=DataArrayDouble.New([1401.,101401.,1602.,101602.,3100.,103100.,3101.,103101.],4,2) ; da.setInfoOnComponents(["Power [MW]","Density [kg/m^3]"])
        f0_1.setArray(da)
        ff0.setFieldProfile(f0_1,mm,0,pfl)
        f0_2=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f0_2.setName("f0")#provoquer error
        f0_2.setTime(9.5,3,4)
        pfl2=DataArrayInt.New([0,1,2,3,4,5,6,8]) ; pfl2.setName("pfltestM1")
        da=DataArrayDouble.New([300.,100300.,301.,100301.,400.,100400.,401.,100401.,402.,100402.,3200.,103200.,3201.,103201.,3203.,103203.],8,2) ; da.setInfoOnComponents(["Power [MW]","Density [kg/m^3]"])#provoquer error
        f0_2.setMesh(mm.getMeshAtLevel(-1)[pfl2])
        f0_2.setArray(da)
        ff0.setFieldProfile(f0_2,mm,-1,pfl2)
        mfd.getFields().shallowCpyGlobs(ff0)
        #
        mfd.unPolyzeMeshes()
        #
        fmts=mfd.getFields()[0]
        self.assertEqual(fmts.getNumberOfTS(),1)
        self.assertEqual(fmts.getTimeSteps(),[(3,4,9.5)])
        arr,entry=fmts.getUndergroundDataArrayExt(3,4)
        self.assertEqual(entry,[((3,0),(38,40)),((4,0),(40,43)),((5,0),(43,46)),((14,0),(46,48)),((16,0),(48,49)),((18,0),(49,50)),((40,0),(0,38))])
        self.assertTrue(arr[38:40].isEqualWithoutConsideringStr(DataArrayDouble([300.0,100300.0,301.0,100301.0],2,2),1e-8))
        self.assertTrue(arr[40:43].isEqualWithoutConsideringStr(DataArrayDouble([400.0,100400.0,401.0,100401.0,402.0,100402.0],3,2),1e-8))
        self.assertTrue(arr[43:46].isEqualWithoutConsideringStr(DataArrayDouble([3200.0,103200.0,3201.0,103201.0,3203.0,103203.0],3,2),1e-8))
        self.assertTrue(arr[46:48].isEqualWithoutConsideringStr(DataArrayDouble([1401.0,101401.0,3100.0,103100.0],2,2),1e-8))
        self.assertTrue(arr[48:49].isEqualWithoutConsideringStr(DataArrayDouble([1602.0,101602.0],1,2),1e-8))
        self.assertTrue(arr[49:50].isEqualWithoutConsideringStr(DataArrayDouble([3101.0,103101.0],1,2),1e-8))
        self.assertEqual(('NewPfl_0','NewPfl_1','NewPfl_2'),fmts.getPflsReallyUsed())
        self.assertEqual([(3,[(0,(38,40),'NewPfl_0','')]),(4,[(0,(40,43),'','')]),(5,[(0,(43,46),'','')]),(14,[(0,(46,48),'NewPfl_1','')]),(16,[(0,(48,49),'NewPfl_2','')]),(18,[(0,(49,50),'','')]),(40,[(1,(0,38),'','')])],fmts.getFieldSplitedByType(3,4))
        self.assertEqual(fmts.getProfile("NewPfl_0").getValues(),[0,1])
        self.assertEqual(fmts.getProfile("NewPfl_1").getValues(),[1,2])
        self.assertEqual(fmts.getProfile("NewPfl_2").getValues(),[2])
        ftest0=fmts.getFieldOnMeshAtLevel(ON_CELLS,3,4,0,mfd.getMeshes()[0])
        self.assertTrue(ftest0.getArray().isEqualWithoutConsideringStr(DataArrayDouble([1401.,101401.,3100.,103100.,1602.,101602.,3101.,103101.],4,2),1e-8))
        self.assertEqual(ftest0.getMesh().getNodalConnectivity().getValues(),[14,4,5,6,7,14,26,27,28,29,16,20,21,22,23,24,25,18,30,31,32,33,34,35,36,37])
        self.assertEqual(ftest0.getMesh().getNodalConnectivityIndex().getValues(),[0,5,10,17,26])
        ftest1=fmts.getFieldOnMeshAtLevel(ON_CELLS,3,4,-1,mfd.getMeshes()[0])
        self.assertTrue(ftest1.getArray().isEqualWithoutConsideringStr(DataArrayDouble([300.,100300.,301.,100301.,400.,100400.,401.,100401.,402.,100402.,3200.,103200.,3201.,103201.,3203.,103203.]),1e-8))
        self.assertEqual(ftest1.getMesh().getNodalConnectivity().getValues(),[3,0,1,2,3,3,4,5,4,6,7,8,9,4,10,11,12,13,4,14,15,16,17,5,18,19,20,21,22,5,23,24,25,26,27,5,31,32,33,34,35,36,37])
        self.assertEqual(ftest1.getMesh().getNodalConnectivityIndex().getValues(),[0,4,8,13,18,23,29,35,43])
        #
        mfd.write(fname,2)
        pass
    pass

unittest.main()
