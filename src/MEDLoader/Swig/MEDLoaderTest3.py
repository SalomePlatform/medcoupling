#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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

class MEDLoaderTest(unittest.TestCase):
    def testMEDMesh1(self):
        fileName="Pyfile18.med"
        mname="ExampleOfMultiDimW"
        medmesh=MEDFileMesh.New(fileName,mname)
        self.assertRaises(InterpKernelException,MEDFileMesh.New,fileName,"")
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
        g1_0=medmesh.getFamily(0,"Family_-3",True)
        g1_1=MEDLoader.ReadUMeshFromFamilies(fileName,mname,0,["Family_-3"]);
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getFamilies(0,["Family_-3","Family_-5"],True)
        g1_1=MEDLoader.ReadUMeshFromFamilies(fileName,mname,0,["Family_-3","Family_-5"]);
        g1_1.setName(g1_0.getName())
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        medmesh.write(outFileName,2);
        self.assertEqual([1,2,4,13,15],medmesh.getGroupArr(0,"mesh2",True).getValues());
        self.assertEqual([1,2,15],medmesh.getFamilyArr(0,"Family_-3",True).getValues());
        self.assertEqual([1,2,4,13,15],medmesh.getFamiliesArr(0,["Family_-5","Family_-3"],True).getValues());
        self.assertEqual([18,1,2,3,4,13,14,15],medmesh.getGroupsArr(0,["mesh2","mesh4","mesh3"],True).getValues());
        famn=medmesh.getFamilyNameGivenId(0)
        self.assertRaises(InterpKernelException,medmesh.getNodeFamilyArr,famn,True);
        #without renum
        self.assertEqual([2,3,5,14,16],medmesh.getGroupArr(0,"mesh2").getValues());
        self.assertEqual([2,3,16],medmesh.getFamilyArr(0,"Family_-3").getValues());
        self.assertEqual([2,3,5,14,16],medmesh.getFamiliesArr(0,["Family_-5","Family_-3"]).getValues());
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
        self.assertTrue(mm.getUnivNameWrStatus())
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
        mmCpy=mm.deepCpy()
        self.assertTrue(mm.isEqual(mmCpy,1e-12)[0]) ; del mm
        mmCpy.write(outFileName,2);
        #
        mm=MEDFileMesh.New(outFileName)
        #
        self.assertEqual([NORM_TRI3,NORM_QUAD4,NORM_POLYGON],mm.getGeoTypesAtLevel(0))
        self.assertEqual([NORM_SEG2,NORM_SEG3],mm.getGeoTypesAtLevel(-1))
        self.assertEqual([NORM_POINT1],mm.getGeoTypesAtLevel(-2))
        mm0=mm.getDirectUndergroundSingleGeoTypeMesh(NORM_POLYGON)
        self.assertTrue(isinstance(mm0,MEDCoupling1DGTUMesh))
        self.assertTrue(mm0.getNodalConnectivity().isEqual(DataArrayInt([6,7,4,3,7,8,5,4])))
        self.assertTrue(mm0.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8])))
        lmm=mm.getDirectUndergroundSingleGeoTypeMeshes(0)
        self.assertEqual(3,len(lmm))
        self.assertTrue(isinstance(lmm[0],MEDCoupling1SGTUMesh))
        self.assertTrue(isinstance(lmm[1],MEDCoupling1SGTUMesh))
        self.assertTrue(isinstance(lmm[2],MEDCoupling1DGTUMesh))
        #
        self.assertTrue(mm.getUnivNameWrStatus())
        self.assertTrue(isinstance(mm.getUnivName(),str))
        self.assertTrue(len(mm.getUnivName())!=0)
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
        self.assertEqual(['Family_-3','Family_-4','Family_-7','Family_10','Family_11'],l)
        mm2.keepFamIdsOnlyOnLevs([3],[-1])
        for lev in mm.getGrpNonEmptyLevelsExt("G2"):
            self.assertEqual(mm.getGroupArr(lev,"G2").getValues(),mm2.getGroupArr(lev,"G2").getValues())
            pass
        l=list(mm2.getFamiliesOnGroup("G2")) ; l.sort()
        self.assertEqual(['Family_-3','Family_-4','Family_-7','Family_10','Family_11'],l)
        #
        self.assertEqual([-7,-7,-6],mm2.getFamilyFieldAtLevel(-1).getValues())
        mm2.getFamilyFieldAtLevel(-1).setIJ(1,0,-8)
        self.assertEqual([-7,-8,-6],mm2.getFamilyFieldAtLevel(-1).getValues())
        self.assertTrue(not mm2.existsFamily("Family_-8"))
        mm2.createGroupOnAll(-1,"GrpOnAllFace")
        self.assertTrue(mm2.existsFamily("Family_-8"))
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
        self.assertTrue(m.getUnivNameWrStatus())
        m.write(outFileName,2);
        mm=MEDFileMesh.New(outFileName)
        self.assertEqual([NORM_HEXA8],mm.getGeoTypesAtLevel(0))
        self.assertTrue(isinstance(mm,MEDFileCMesh))
        self.assertTrue(isinstance(mm.getUnivName(),str))
        self.assertTrue(len(mm.getUnivName())!=0)
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
        self.assertRaises(InterpKernelException,MEDFileField1TS,fileName)#throw because no field in file fileName
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
        self.assertEqual(3,len(ff))
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
        ff=MEDFileFieldMultiTS.New("Pyfile19.med")
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
        ff3=MEDFileField1TS.New("Pyfile13.med","MyFirstFieldOnGaussPoint")
        f3=ff3.getFieldAtLevel(ON_GAUSS_PT,0)
        self.assertTrue(f.isEqual(f3,1e-12,1e-12))
        ff4=MEDFileField1TS.New("Pyfile13.med")
        f4=ff4.getFieldAtLevel(ON_GAUSS_PT,0)
        self.assertTrue(f.isEqual(f4,1e-12,1e-12))
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
        f1InvalidCpy=f1.deepCpy()
        f1InvalidCpy.setDiscretization(MEDCouplingFieldDiscretizationGauss())
        f1InvalidCpy2=f1.deepCpy()
        f1InvalidCpy2.setDiscretization(MEDCouplingFieldDiscretizationGauss())
        m1=f1.getMesh()
        mm1=MEDFileUMesh.New()
        mm1.setCoords(m1.getCoords())
        mm1.setMeshAtLevel(0,m1)
        mm1.setName(m1.getName())
        mm1.write(fname,2)
        ff1=MEDFileField1TS.New()
        self.assertRaises(InterpKernelException,ff1.setFieldNoProfileSBT,f1InvalidCpy) # fails because no Gauss localization per cell set !*
        f1InvalidCpy2.getDiscretization().setArrayOfDiscIds(f1.getDiscretization().getArrayOfDiscIds()) # fails because no Gauss localization set whereas gauss locid per cell given !
        self.assertRaises(InterpKernelException,ff1.setFieldNoProfileSBT,f1InvalidCpy2)
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
        ff1=ff1.deepCpy()
        ff1.write(fname,0)
        #
        vals,pfl=ff1.getFieldWithProfile(ON_CELLS,0,mm1) ; vals.setName("")
        self.assertTrue(pfl.isEqualWithoutConsideringStr(da))# profiles names cannot be contracted in pfl array name
        self.assertTrue(vals.isEqual(d,1e-14))
        #
        ff2=MEDFileField1TS.New(fname,f1.getName(),-1,-1)
        ff3=MEDFileField1TS.New(fname,f1.getName(),-1,-1)
        ff2.deepCpyGlobs(ff3)
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
        ff1=ff1.deepCpy()
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
        ffs=ffs.deepCpy()
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
        self.assertEqual(('Loc_MyFirstFieldOnGaussPoint_NORM_QUAD4_1','Loc_MyFirstFieldOnGaussPoint_NORM_TRI3_0','Loc_MyFirstFieldOnGaussPoint_NORM_TRI6_2'),d.getFields().getFieldAtPos(0).getLocs())
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
        self.assertEqual(('Family_-2','Family_-3'),mm.getFamiliesNames())
        self.assertEqual(('Family_-2',),mm.getFamiliesOnGroup('g1'))
        self.assertEqual(('Family_-3',),mm.getFamiliesOnGroup('g2'))
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
        self.assertEqual(('Family_-2', 'Family_-4', 'Family_-5'),mm.getFamiliesNames())
        self.assertEqual(('Family_-2', 'Family_-4'),mm.getFamiliesOnGroup('g1'))
        self.assertEqual(('Family_-5',),mm.getFamiliesOnGroup('g2'))
        self.assertEqual(('Family_-4','Family_-5',),mm.getFamiliesOnGroup('g3'))
        mm.assignFamilyNameWithGroupName() # here it does nothing because no such group-family bijection found
        self.assertEqual(('g1','g2','g3'),mm.getGroupsNames())
        self.assertEqual(('Family_-2', 'Family_-4', 'Family_-5'),mm.getFamiliesNames())
        self.assertEqual(('Family_-2', 'Family_-4'),mm.getFamiliesOnGroup('g1'))
        self.assertEqual(('Family_-5',),mm.getFamiliesOnGroup('g2'))
        self.assertEqual(('Family_-4','Family_-5',),mm.getFamiliesOnGroup('g3'))
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
        m=MEDFileMesh.New(fname)
        m=MEDFileMesh.New(fname,"ExampleOfMultiDimW",-1,-1)
        m=MEDFileMesh.New(fname)
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
        self.assertEqual(mm.unPolyze()[:3],(True,[[3,2,0],[4,3,2],[5,4,5],[14,2,9],[16,3,11],[31,2,14]],[[3,3,0],[4,3,3],[5,3,6],[14,3,9],[16,3,12],[18,1,15]]))
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

    def testGaussWriteOnPfl1(self):
        fname="Pyfile49.med"
        fname2="Pyfile50.med"
        coords=DataArrayDouble([0.,0.,0.,1.,1.,1.,1.,0.,0.,0.5,0.5,1.,1.,0.5,0.5,0.],8,2)
        mQ8=MEDCouplingUMesh("",2) ; mQ8.setCoords(coords)
        mQ8.allocateCells(1)
        mQ8.insertNextCell(NORM_QUAD8,range(8))
        mQ8.finishInsertingCells()
        mQ4=MEDCouplingUMesh("",2) ; mQ4.setCoords(coords)
        mQ4.allocateCells(1)
        mQ4.insertNextCell(NORM_QUAD4,range(4))
        mQ4.finishInsertingCells()
        mT3=MEDCouplingUMesh("",2) ; mT3.setCoords(coords)
        mT3.allocateCells(1)
        mT3.insertNextCell(NORM_TRI3,range(3))
        mT3.finishInsertingCells()
        
        tr=[[0.,4.],[2.,4.],[4.,4.],[6.,4.],[8.,4.],[10.,4.],[12.,4.],[14.,4.],[16.,4.],[18.,4.],[20.,4.],[0.,0.],[2.,0.], [0.,2.],[2.,2.],[4.,2.],[6.,2.],[8.,2.],[10.,2.],[12.,2.]]
        ms=11*[mT3]+2*[mQ4]+7*[mQ8]
        ms[:]=(elt.deepCpy() for elt in ms)
        for m,t in zip(ms,tr):
            d=m.getCoords() ; d+= t
            pass
        m=MEDCouplingUMesh.MergeUMeshes(ms)
        m.setName("mesh")
        m2=m[:13] ; m2.setName(m.getName())
        ### Use case 1 : Pfl on all tri3 and on all quad4. If we were on CELLS or GAUSS_NE no pfl were needed. But here 2 discs in tri3.
        ### So here 2 pfls will be created (pfl_TRI3_loc_0 and pfl_TRI3_loc_1)
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME)
        f.setMesh(m2)
        f.setTime(4.5,1,2)
        da=DataArrayDouble(34) ; da.iota(3.)
        f.setArray(da)
        f.setName("fieldCellOnPflWithoutPfl")
        fInvalid=f.deepCpy()
        f.setGaussLocalizationOnCells([0,1,2,3,4,5,6,7,8],[0.,0.,1.,0.,1.,1.],[0.3,0.3,0.7,0.7],[0.8,0.2])
        f.setGaussLocalizationOnCells([9,10],[0.,0.,1.,0.,1.,1.],[0.3,0.3,0.7,0.7,0.8,0.8],[0.8,0.07,0.13])
        f.setGaussLocalizationOnCells([11,12],[0.,0.,1.,0.,1.,1.,0.,1.],[0.3,0.3,0.7,0.7,0.8,0.8,0.8,0.8,0.8,0.8],[0.8,0.07,0.1,0.01,0.02])
        f.checkCoherency()
        fInvalid2=fInvalid.deepCpy()
        fInvalid2.getDiscretization().setArrayOfDiscIds(f.getDiscretization().getArrayOfDiscIds())
        #
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        mm.write(fname,2)
        #
        f1ts=MEDFileField1TS.New()
        pfl=DataArrayInt(range(13)) ; pfl.setName("pfl")
        self.assertRaises(InterpKernelException,f1ts.setFieldProfile,fInvalid,mm,0,pfl) # fails because no Gauss localization per cell set !
        self.assertRaises(InterpKernelException,f1ts.setFieldProfile,fInvalid2,mm,0,pfl) # fails because no Gauss localization set whereas gauss locid per cell given !
        f1ts.setFieldProfile(f,mm,0,pfl)
        f1ts.write(fname,0)
        #
        self.assertEqual(f1ts.getPfls(),('pfl_NORM_TRI3_loc_0', 'pfl_NORM_TRI3_loc_1'))
        self.assertEqual(f1ts.getPflsReallyUsed(),('pfl_NORM_TRI3_loc_0', 'pfl_NORM_TRI3_loc_1'))
        da1=DataArrayInt([0,1,2,3,4,5,6,7,8]) ; da1.setName("pfl_NORM_TRI3_loc_0")
        self.assertTrue(f1ts.getProfile("pfl_NORM_TRI3_loc_0").isEqual(da1))
        da1=DataArrayInt([9,10]) ; da1.setName("pfl_NORM_TRI3_loc_1")
        self.assertTrue(f1ts.getProfile("pfl_NORM_TRI3_loc_1").isEqual(da1))
        self.assertEqual(f1ts.getLocs(),('Loc_fieldCellOnPflWithoutPfl_NORM_TRI3_0', 'Loc_fieldCellOnPflWithoutPfl_NORM_TRI3_1', 'Loc_fieldCellOnPflWithoutPfl_NORM_QUAD4_2'))
        self.assertEqual(f1ts.getLocsReallyUsed(),('Loc_fieldCellOnPflWithoutPfl_NORM_TRI3_0', 'Loc_fieldCellOnPflWithoutPfl_NORM_TRI3_1', 'Loc_fieldCellOnPflWithoutPfl_NORM_QUAD4_2'))
        #
        dataRead=MEDFileData.New(fname)
        mRead=dataRead.getMeshes()[0]
        f1tsRead=dataRead.getFields()[0][0]
        f1tsRead.getFieldOnMeshAtLevel(ON_GAUSS_PT,0,mRead)
        f2=f1tsRead.getFieldOnMeshAtLevel(ON_GAUSS_PT,0,mRead)
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        f2_bis=MEDLoader.ReadFieldGauss(fname,m.getName(),0,f.getName(),f.getTime()[1],f.getTime()[2])
        f2_bis.checkCoherency()
        self.assertTrue(f.isEqual(f2_bis,1e-12,1e-12))
        #
        MEDLoader.WriteField(fname2,f,True)
        f2_ter=MEDLoader.ReadFieldGauss(fname2,m.getName(),0,f.getName(),f.getTime()[1],f.getTime()[2])
        self.assertTrue(f.isEqual(f2_ter,1e-12,1e-12))
        ## Use case 2 : Pfl on part tri3 with 2 disc and on part quad8 with 1 disc
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME)
        pfl=DataArrayInt([1,2,5,6,8,9,15,16,17,18]) ; pfl.setName("pfl2")
        m2=m[pfl] ; m2.setName(m.getName())
        f.setMesh(m2)
        f.setTime(4.5,1,2)
        da=DataArrayDouble(35) ; da.iota(3.)
        f.setArray(da)
        f.setName("fieldCellOnPflWithoutPfl2")
        f.setGaussLocalizationOnCells([0,1,3],[0.,0.,1.,0.,1.,1.],[0.3,0.3,0.7,0.7],[0.8,0.2])
        f.setGaussLocalizationOnCells([2,4,5],[0.,0.,1.,0.,1.,1.],[0.3,0.3,0.7,0.7,0.8,0.8],[0.8,0.07,0.13])
        f.setGaussLocalizationOnCells([6,7,8,9],[0.,0.,1.,0.,1.,1.,0.,1.,0.5,0.,1.,0.5,0.5,1.,0.,0.5],[0.3,0.3,0.7,0.7,0.8,0.8,0.8,0.8,0.8,0.8],[0.8,0.07,0.1,0.01,0.02])
        f.checkCoherency()
        #
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        mm.write(fname,2)
        f1ts=MEDFileField1TS.New()
        f1ts.setFieldProfile(f,mm,0,pfl)
        self.assertEqual(f1ts.getPfls(),('pfl2_NORM_TRI3_loc_0','pfl2_NORM_TRI3_loc_1','pfl2_NORM_QUAD8_loc_2'))
        self.assertEqual(f1ts.getProfile("pfl2_NORM_TRI3_loc_0").getValues(),[1,2,6])
        self.assertEqual(f1ts.getProfile("pfl2_NORM_TRI3_loc_1").getValues(),[5,8,9])
        self.assertEqual(f1ts.getProfile("pfl2_NORM_QUAD8_loc_2").getValues(),[2,3,4,5])
        f1ts.write(fname,0)
        dataRead=MEDFileData.New(fname)
        mRead=dataRead.getMeshes()[0]
        f1tsRead=dataRead.getFields()[0][0]
        f1tsRead.getFieldOnMeshAtLevel(ON_GAUSS_PT,0,mRead)
        f3=f1tsRead.getFieldOnMeshAtLevel(ON_GAUSS_PT,0,mRead)
        f3.renumberCells([0,1,3,2,4,5,6,7,8,9])
        self.assertTrue(f.isEqual(f3,1e-12,1e-12))
        f3_bis=MEDLoader.ReadFieldGauss(fname,m.getName(),0,f.getName(),f.getTime()[1],f.getTime()[2])
        f3_bis.renumberCells([0,1,3,2,4,5,6,7,8,9])
        self.assertTrue(f.isEqual(f3_bis,1e-12,1e-12))
        #
        MEDLoader.WriteField(fname2,f,True)
        f3_ter=MEDLoader.ReadFieldGauss(fname2,m.getName(),0,f.getName(),f.getTime()[1],f.getTime()[2])
        f3_ter.renumberCells([0,1,3,2,4,5,6,7,8,9])
        self.assertTrue(f.isEqual(f3_ter,1e-12,1e-12))
        ## Use case 3 : no pfl but creation of pfls due to gauss pts
        f=MEDCouplingFieldDouble.New(ON_GAUSS_PT,ONE_TIME)
        f.setMesh(m)
        f.setTime(4.5,1,2)
        da=DataArrayDouble(60) ; da.iota(3.)
        f.setArray(da)
        f.setName("fieldCellWithoutPfl")
        f.setGaussLocalizationOnCells([0,1,2,3,4,5,6,7,8],[0.,0.,1.,0.,1.,1.],[0.3,0.3,0.7,0.7],[0.8,0.2])
        f.setGaussLocalizationOnCells([9,10],[0.,0.,1.,0.,1.,1.],[0.3,0.3,0.7,0.7,0.8,0.8],[0.8,0.07,0.13])
        f.setGaussLocalizationOnCells([11,12],[0.,0.,1.,0.,1.,1.,0.,1.],[0.3,0.3,0.7,0.7,0.8,0.8,0.8,0.8,0.8,0.8],[0.8,0.07,0.1,0.01,0.02])
        f.setGaussLocalizationOnCells([13,14,15,17,18],[0.,0.,1.,0.,1.,1.,0.,1.,0.5,0.,1.,0.5,0.5,1.,0.,0.5],[0.3,0.3,0.7,0.7,0.8,0.8,0.8,0.8],[0.8,0.1,0.03,0.07])
        f.setGaussLocalizationOnCells([16,19],[0.,0.,1.,0.,1.,1.,0.,1.,0.5,0.,1.,0.5,0.5,1.,0.,0.5],[0.3,0.3,0.7,0.7,0.8,0.8],[0.8,0.1,0.1])
        f.checkCoherency()
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m) 
        f1ts=MEDFileField1TS.New()
        f1ts.setFieldNoProfileSBT(f)
        self.assertEqual(f1ts.getPfls(),('Pfl_fieldCellWithoutPfl_NORM_TRI3_0','Pfl_fieldCellWithoutPfl_NORM_TRI3_1','Pfl_fieldCellWithoutPfl_NORM_QUAD8_3','Pfl_fieldCellWithoutPfl_NORM_QUAD8_4'))
        self.assertEqual(f1ts.getProfile("Pfl_fieldCellWithoutPfl_NORM_TRI3_0").getValues(),[0,1,2,3,4,5,6,7,8])
        self.assertEqual(f1ts.getProfile("Pfl_fieldCellWithoutPfl_NORM_TRI3_1").getValues(),[9,10])
        self.assertEqual(f1ts.getProfile("Pfl_fieldCellWithoutPfl_NORM_QUAD8_3").getValues(),[0,1,2,4,5])
        self.assertEqual(f1ts.getProfile("Pfl_fieldCellWithoutPfl_NORM_QUAD8_4").getValues(),[3,6])
        mm.write(fname,2)
        f1ts.write(fname,0)
        #
        dataRead=MEDFileData.New(fname)
        mRead=dataRead.getMeshes()[0]
        f1tsRead=dataRead.getFields()[0][0]
        f3=f1tsRead.getFieldOnMeshAtLevel(ON_GAUSS_PT,0,mRead)
        f3.renumberCells([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,16,19])
        self.assertTrue(f.isEqual(f3,1e-12,1e-12))
        f3_bis=MEDLoader.ReadFieldGauss(fname,m.getName(),0,f.getName(),f.getTime()[1],f.getTime()[2])
        f3_bis.renumberCells([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,16,19])
        self.assertTrue(f.isEqual(f3_bis,1e-12,1e-12))
        #
        MEDLoader.WriteField(fname2,f,True)
        f3_ter=MEDLoader.ReadFieldGauss(fname2,m.getName(),0,f.getName(),f.getTime()[1],f.getTime()[2])
        f3_ter.renumberCells([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,16,19])
        self.assertTrue(f.isEqual(f3_ter,1e-12,1e-12))
        pass

    # Testing profile on nodes when the profile is identity but not on all nodes.
    def testMEDFieldPflOnNode1(self):
        fname="Pyfile51.med"
        coo=DataArrayDouble([0.,0.,0.5,0.,1.,0.,0.,0.5,0.5,0.5,1.,0.5,0.,1.,0.5,1.,1.,1.],9,2)
        m0=MEDCouplingUMesh("Mesh",2)
        m0.allocateCells(5)
        m0.insertNextCell(NORM_TRI3,[1,4,2])
        m0.insertNextCell(NORM_TRI3,[4,5,2])
        m0.insertNextCell(NORM_QUAD4,[0,3,4,1])
        m0.insertNextCell(NORM_QUAD4,[3,6,7,4])
        m0.insertNextCell(NORM_QUAD4,[4,7,8,5])
        m0.finishInsertingCells()
        m0.setCoords(coo)
        m1=MEDCouplingUMesh(m0.getName(),1)
        m1.allocateCells(9)
        conn1=[0,1,0,3,3,4,4,1,5,4,2,4,1,2,3,6,5,8]
        for i in xrange(9):
            m1.insertNextCell(NORM_SEG2,conn1[2*i:2*i+2])
            pass
        m1.finishInsertingCells()
        m1.setCoords(coo)
        #
        m=MEDFileUMesh()
        m.setMeshAtLevel(0,m0)
        m.setMeshAtLevel(-1,m1)
        #
        dt=3 ; it=2 ; tim=4.5
        fieldNode0=MEDCouplingFieldDouble(ON_NODES,ONE_TIME)
        fieldNode0.setName("fieldNode0")
        fieldNode0.setTime(tim,dt,it)
        pfl0=DataArrayInt([0,1,2,3,4]) ; pfl0.setName("PflIdentity0") # important to keep like that
        arr=DataArrayDouble([10,11,12,13,14])
        fieldNode0.setArray(arr)
        f0=MEDFileField1TS()
        f0.setFieldProfile(fieldNode0,m,0,pfl0)
        m.write(fname,2) ; f0.write(fname,0)
        fieldNode1=MEDCouplingFieldDouble(ON_NODES,ONE_TIME)
        fieldNode1.setName("fieldNode1")
        fieldNode1.setTime(tim,dt,it)
        pfl1=DataArrayInt([0,1,2,3,4,5,6]) ; pfl1.setName("PflIdentity1")
        arr1=DataArrayDouble([20,21,22,23,24,25,26])
        fieldNode1.setArray(arr1)
        f1=MEDFileField1TS()
        f1.setFieldProfile(fieldNode1,m,-1,pfl1)
        f1.write(fname,0)
        del m,f0,m0,m1,f1
        ## Reading from file
        m=MEDFileMesh.New(fname)
        m0=m.getMeshAtLevel(0)
        m00=m0.deepCpy() ; m00=m00[[0,2]] ; m00.setName(m.getName()) ; m00.zipCoords()
        fieldNode0.setMesh(m00)
        f0=MEDFileField1TS.New(fname,fieldNode0.getName(),dt,it)
        ff0_1=f0.getFieldOnMeshAtLevel(ON_NODES,m0)
        ff0_1.checkCoherency()
        self.assertTrue(ff0_1.isEqual(fieldNode0,1e-12,1e-12))
        ff0_2=f0.getFieldAtLevel(ON_NODES,0)
        ff0_2.checkCoherency()
        self.assertTrue(ff0_2.isEqual(fieldNode0,1e-12,1e-12))
        ff0_3=f0.getFieldOnMeshAtLevel(ON_NODES,0,m)
        ff0_3.checkCoherency()
        self.assertTrue(ff0_3.isEqual(fieldNode0,1e-12,1e-12))
        ff0_4=MEDLoader.ReadFieldNode(fname,m.getName(),0,fieldNode0.getName(),dt,it)
        ff0_4.checkCoherency()
        self.assertTrue(ff0_4.isEqual(fieldNode0,1e-12,1e-12))
        f1=MEDFileField1TS.New(fname,fieldNode1.getName(),dt,it)
        m1=m.getMeshAtLevel(-1)
        m10=m1.deepCpy() ; m10=m10[[0,1,2,3,4,5,6,7]] ; m10.setName(m.getName()) ; m10.zipCoords()
        fieldNode1.setMesh(m10)
        ff1_1=f1.getFieldOnMeshAtLevel(ON_NODES,m1)
        ff1_1.checkCoherency()
        self.assertTrue(ff1_1.isEqual(fieldNode1,1e-12,1e-12))
        ff1_2=f1.getFieldAtLevel(ON_NODES,-1)
        ff1_2.checkCoherency()
        self.assertTrue(ff1_2.isEqual(fieldNode1,1e-12,1e-12))
        ff1_3=f1.getFieldOnMeshAtLevel(ON_NODES,-1,m)
        ff1_3.checkCoherency()
        self.assertTrue(ff1_3.isEqual(fieldNode1,1e-12,1e-12))
        ff1_4=MEDLoader.ReadFieldNode(fname,m.getName(),-1,fieldNode1.getName(),dt,it)
        ff1_4.checkCoherency()
        self.assertTrue(ff1_4.getMesh().isEqual(m10,1e-12))
        self.assertRaises(InterpKernelException,f1.getFieldOnMeshAtLevel,ON_NODES,m0) # error because impossible to build a sub mesh at level 0 lying on nodes [0,1,2,3,4,5,6]
        self.assertRaises(InterpKernelException,f1.getFieldAtLevel,ON_NODES,0) # error because impossible to build a sub mesh at level 0 lying on nodes [0,1,2,3,4,5,6]
        self.assertRaises(InterpKernelException,f1.getFieldOnMeshAtLevel,ON_NODES,0,m) # error because impossible to build a sub mesh at level 0 lying on nodes [0,1,2,3,4,5,6]
        arr_r,pfl1_r=f1.getFieldWithProfile(ON_NODES,-1,m)
        arr_r.setName(fieldNode1.getArray().getName())
        self.assertTrue(arr_r.isEqual(fieldNode1.getArray(),1e-12))
        pfl1_r.setName(pfl1.getName())
        self.assertTrue(pfl1_r.isEqual(pfl1))
        pass
    
        # Testing profile on nodes when the profile is identity but not on all nodes.
    def testMEDFieldPflOnCell1(self):
        fname="Pyfile52.med"
        coo=DataArrayDouble([0.,0.,0.5,0.,1.,0.,0.,0.5,0.5,0.5,1.,0.5,0.,1.,0.5,1.,1.,1.],9,2)
        m0=MEDCouplingUMesh("Mesh",2)
        m0.allocateCells(5)
        m0.insertNextCell(NORM_TRI3,[1,4,2])
        m0.insertNextCell(NORM_TRI3,[4,5,2])
        m0.insertNextCell(NORM_QUAD4,[0,3,4,1])
        m0.insertNextCell(NORM_QUAD4,[3,6,7,4])
        m0.insertNextCell(NORM_QUAD4,[4,7,8,5])
        m0.finishInsertingCells()
        m0.setCoords(coo)
        m1=MEDCouplingUMesh(m0.getName(),1)
        m1.allocateCells(9)
        conn1=[0,1,0,3,3,4,4,1,5,4,2,4,1,2,3,6,5,8]
        for i in xrange(9):
            m1.insertNextCell(NORM_SEG2,conn1[2*i:2*i+2])
            pass
        m1.finishInsertingCells()
        m1.setCoords(coo)
        #
        m=MEDFileUMesh()
        m.setMeshAtLevel(0,m0)
        m.setMeshAtLevel(-1,m1)
        #
        dt=3 ; it=2 ; tim=4.5
        fieldCell0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
        fieldCell0.setName("fieldCell0")
        fieldCell0.setTime(tim,dt,it)
        pfl0=DataArrayInt([0,1,2]) ; pfl0.setName("PflIdentity0") # important to keep like that
        arr=DataArrayDouble([10,11,12])
        fieldCell0.setArray(arr)
        f0=MEDFileField1TS()
        f0.setFieldProfile(fieldCell0,m,0,pfl0)
        m.write(fname,2) ; f0.write(fname,0)
        fieldCell1=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME)
        fieldCell1.setName("fieldCell1")
        fieldCell1.setTime(tim,dt,it)
        pfl1=DataArrayInt([0,1,2,3,4,5,6]) ; pfl1.setName("PflIdentity1")
        arr1=DataArrayDouble([20,21,22,23,24,25,26])
        fieldCell1.setArray(arr1)
        f1=MEDFileField1TS()
        f1.setFieldProfile(fieldCell1,m,-1,pfl1)
        f1.write(fname,0)
        del m,f0,m0,m1,f1
        ## Reading from file
        m=MEDFileMesh.New(fname)
        m0=m.getMeshAtLevel(0)
        m00=m0.deepCpy() ; m00=m00[pfl0] ; m00.setName(m.getName())
        fieldCell0.setMesh(m00)
        f0=MEDFileField1TS.New(fname,fieldCell0.getName(),dt,it)
        ff0_1=f0.getFieldOnMeshAtLevel(ON_CELLS,m0)
        ff0_1.checkCoherency()
        self.assertTrue(ff0_1.isEqual(fieldCell0,1e-12,1e-12))
        ff0_2=f0.getFieldAtLevel(ON_CELLS,0)
        ff0_2.checkCoherency()
        self.assertTrue(ff0_2.isEqual(fieldCell0,1e-12,1e-12))
        ff0_3=f0.getFieldOnMeshAtLevel(ON_CELLS,0,m)
        ff0_3.checkCoherency()
        self.assertTrue(ff0_3.isEqual(fieldCell0,1e-12,1e-12))
        ff0_4=MEDLoader.ReadFieldCell(fname,m.getName(),0,fieldCell0.getName(),dt,it)
        ff0_4.checkCoherency()
        self.assertTrue(ff0_4.isEqual(fieldCell0,1e-12,1e-12))
        f1=MEDFileField1TS.New(fname,fieldCell1.getName(),dt,it)
        m1=m.getMeshAtLevel(-1)
        m10=m1.deepCpy() ; m10=m10[pfl1] ; m10.setName(m.getName())
        fieldCell1.setMesh(m10)
        ff1_1=f1.getFieldOnMeshAtLevel(ON_CELLS,m1)
        ff1_1.checkCoherency()
        self.assertTrue(ff1_1.isEqual(fieldCell1,1e-12,1e-12))
        ff1_2=f1.getFieldAtLevel(ON_CELLS,-1)
        ff1_2.checkCoherency()
        self.assertTrue(ff1_2.isEqual(fieldCell1,1e-12,1e-12))
        ff1_3=f1.getFieldOnMeshAtLevel(ON_CELLS,-1,m)
        ff1_3.checkCoherency()
        self.assertTrue(ff1_3.isEqual(fieldCell1,1e-12,1e-12))
        ff1_4=MEDLoader.ReadFieldCell(fname,m.getName(),-1,fieldCell1.getName(),dt,it)
        ff1_4.checkCoherency()
        self.assertTrue(ff1_4.getMesh().isEqual(m10,1e-12))
        self.assertRaises(InterpKernelException,f1.getFieldOnMeshAtLevel,ON_CELLS,m0) # error because impossible to build a sub mesh at level 0 lying on cells [0,1,2,3,4,5,6]
        self.assertRaises(InterpKernelException,f1.getFieldAtLevel,ON_CELLS,0) # error because impossible to build a sub mesh at level 0 lying on cells [0,1,2,3,4,5,6]
        self.assertRaises(InterpKernelException,f1.getFieldOnMeshAtLevel,ON_CELLS,0,m) # error because impossible to build a sub mesh at level 0 lying on cells [0,1,2,3,4,5,6]
        arr_r,pfl1_r=f1.getFieldWithProfile(ON_CELLS,-1,m)
        arr_r.setName(fieldCell1.getArray().getName())
        self.assertTrue(arr_r.isEqual(fieldCell1.getArray(),1e-12))
        pfl1_r.setName(pfl1.getName())
        self.assertTrue(pfl1_r.isEqual(pfl1))
        pass

    def testMEDFileUMeshZipCoords1(self):
        m=MEDFileUMesh()
        coo=DataArrayDouble(30) ; coo.iota(1.) ; coo.rearrange(3) ; coo.setInfoOnComponents(["aaa [b]","cc [dd]", "e [fff]"])
        m0=MEDCouplingUMesh("toto",2) ; m0.allocateCells(0) ; m0.insertNextCell(NORM_TRI3,[1,2,3]) ; m0.insertNextCell(NORM_QUAD4,[2,4,3,4]) ; m0.insertNextCell(NORM_POLYGON,[1,6,6,6,2])
        m1=MEDCouplingUMesh("toto",1) ; m1.allocateCells(0) ; m1.insertNextCell(NORM_SEG2,[1,6]) ; m1.insertNextCell(NORM_SEG2,[7,3])
        m2=MEDCouplingUMesh("toto",0) ; m2.allocateCells(0) ; m2.insertNextCell(NORM_POINT1,[2]) ; m2.insertNextCell(NORM_POINT1,[6]) ; m2.insertNextCell(NORM_POINT1,[8])
        m0.setCoords(coo) ; m.setMeshAtLevel(0,m0)
        m1.setCoords(coo) ; m.setMeshAtLevel(-1,m1)
        m2.setCoords(coo) ; m.setMeshAtLevel(-2,m2)
        numCoo=DataArrayInt(10) ; numCoo.iota(3) ; m.setRenumFieldArr(1,numCoo)
        famCoo=DataArrayInt(10) ; famCoo.iota(4) ; m.setFamilyFieldArr(1,famCoo)
        da=DataArrayInt([20,30,40]) ; m.setRenumFieldArr(0,da) ; da=DataArrayInt([200,300,400]) ; m.setFamilyFieldArr(0,da)
        da=DataArrayInt([50,60]) ; m.setRenumFieldArr(-1,da) ; da=DataArrayInt([500,600]) ; m.setFamilyFieldArr(-1,da)
        da=DataArrayInt([70,80,90]) ; m.setRenumFieldArr(-2,da) ; da=DataArrayInt([700,800,900]) ; m.setFamilyFieldArr(-2,da)
        o2n=m.zipCoords()
        self.assertTrue(o2n.isEqual(DataArrayInt([-1,0,1,2,3,-1,4,5,6,-1])))
        self.assertTrue(m.getNumberFieldAtLevel(1).isEqual(DataArrayInt([4,5,6,7,9,10,11])))
        self.assertTrue(m.getFamilyFieldAtLevel(1).isEqual(DataArrayInt([5,6,7,8,10,11,12])))
        self.assertTrue(m.getMeshAtLevel(0).getNodalConnectivity().isEqual(DataArrayInt([3,0,1,2,4,1,3,2,3,5,0,4,4,4,1])))
        self.assertTrue(m.getMeshAtLevel(0).getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,9,15])))
        self.assertTrue(m.getMeshAtLevel(-1).getNodalConnectivity().isEqual(DataArrayInt([1,0,4,1,5,2])))
        self.assertTrue(m.getMeshAtLevel(-1).getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6])))
        self.assertTrue(m.getMeshAtLevel(-2).getNodalConnectivity().isEqual(DataArrayInt([0,1,0,4,0,6])))
        self.assertTrue(m.getMeshAtLevel(-2).getNodalConnectivityIndex().isEqual(DataArrayInt([0,2,4,6])))
        pass

    def testMEDUMeshAddNodeGroup1(self):
        fname="Pyfile53.med"
        m=MEDFileUMesh()
        coo=DataArrayDouble(39) ; coo.iota(1.) ; coo.rearrange(3) ; coo.setInfoOnComponents(["aaa [b]","cc [dd]", "e [fff]"])
        m0=MEDCouplingUMesh("toto",2) ; m0.allocateCells(0) ; m0.insertNextCell(NORM_TRI3,[1,2,3]) ; m0.insertNextCell(NORM_QUAD4,[2,4,3,4]) ; m0.insertNextCell(NORM_POLYGON,[1,6,6,6,2])
        m1=MEDCouplingUMesh("toto",1) ; m1.allocateCells(0) ; m1.insertNextCell(NORM_SEG2,[1,6]) ; m1.insertNextCell(NORM_SEG2,[7,3])
        m2=MEDCouplingUMesh("toto",0) ; m2.allocateCells(0) ; m2.insertNextCell(NORM_POINT1,[2]) ; m2.insertNextCell(NORM_POINT1,[6]) ; m2.insertNextCell(NORM_POINT1,[8])
        m0.setCoords(coo) ; m.setMeshAtLevel(0,m0)
        m1.setCoords(coo) ; m.setMeshAtLevel(-1,m1)
        m2.setCoords(coo) ; m.setMeshAtLevel(-2,m2)
        #
        mm=m.deepCpy()
        famCoo=DataArrayInt([0,2,0,3,2,0,-1,0,0,0,0,-1,3]) ; mm.setFamilyFieldArr(1,famCoo)
        da0=DataArrayInt([0,0,0]) ; mm.setFamilyFieldArr(0,da0)
        da1=DataArrayInt([0,3]) ; mm.setFamilyFieldArr(-1,da1)
        da2=DataArrayInt([0,0,0]) ; mm.setFamilyFieldArr(-2,da2)
        mm.setFamilyId("MyFam",2)
        mm.setFamilyId("MyOtherFam",3)
        mm.setFamilyId("MyOther-1",-1)
        mm.setFamiliesOnGroup("grp0",["MyOtherFam"])
        mm.setFamiliesOnGroup("grpA",["MyOther-1"])
        #
        daTest=DataArrayInt([1,3,4,6,9,10,12]) ; daTest.setName("grp1")
        mm.addNodeGroup(daTest)
        self.assertTrue(mm.getGroupArr(1,daTest.getName()).isEqual(daTest))
        self.assertTrue(mm.getFamilyFieldAtLevel(1).isEqual(DataArrayInt([6,2,6,8,2,6,5,6,6,7,7,4,8])))
        for lev,arr in [(0,da0),(-1,da1),(-2,da2)]:
            self.assertTrue(mm.getFamilyFieldAtLevel(lev).isEqual(arr))
            pass
        self.assertEqual(mm.getFamiliesNames(),('Family_4','Family_5','Family_7','Family_8','MyFam','MyOther-1','MyOtherFam'))
        self.assertEqual(mm.getGroupsNames(),('grp0','grp1','grpA'))
        self.assertEqual(mm.getFamilyNameGivenId(3),'MyOtherFam')
        self.assertEqual(mm.getFamilyNameGivenId(2),'MyFam')
        for famName,famId in [('Family_4',4),('Family_5',5),('Family_7',7),('Family_8',8)]:
            self.assertEqual(mm.getFamilyNameGivenId(famId),famName)
            pass
        self.assertEqual(mm.getFamiliesOnGroup("grp0"),('MyOtherFam','Family_8'))
        da=DataArrayInt([3,12]) ; da.setName("grp0")
        self.assertTrue(mm.getGroupArr(1,"grp0").isEqual(da))
        da.setValues([1])
        self.assertTrue(mm.getGroupArr(-1,"grp0").isEqual(da))
        mm.write(fname,2)
        mm=MEDFileMesh.New(fname)
        self.assertTrue(mm.getGroupArr(1,daTest.getName()).isEqual(daTest))
        self.assertTrue(mm.getFamilyFieldAtLevel(1).isEqual(DataArrayInt([6,2,6,8,2,6,5,6,6,7,7,4,8])))
        for lev,arr in [(0,da0),(-1,da1),(-2,da2)]:
            self.assertTrue(mm.getFamilyFieldAtLevel(lev).isEqual(arr))
            pass
        self.assertEqual(mm.getFamiliesNames(),('FAMILLE_ZERO','Family_4','Family_5','Family_7','Family_8','MyFam','MyOther-1','MyOtherFam'))
        self.assertEqual(mm.getGroupsNames(),('grp0','grp1','grpA'))
        self.assertEqual(mm.getFamilyNameGivenId(3),'MyOtherFam')
        self.assertEqual(mm.getFamilyNameGivenId(2),'MyFam')
        for famName,famId in [('Family_4',4),('Family_5',5),('Family_7',7),('Family_8',8)]:
            self.assertEqual(mm.getFamilyNameGivenId(famId),famName)
            pass
        self.assertEqual(mm.getFamiliesOnGroup("grp0"),('Family_8','MyOtherFam'))
        da=DataArrayInt([3,12]) ; da.setName("grp0")
        self.assertTrue(mm.getGroupArr(1,"grp0").isEqual(da))
        da.setValues([1])
        self.assertTrue(mm.getGroupArr(-1,"grp0").isEqual(da))
        pass

    def testMEDUMeshAddGroup1(self):
        fname="Pyfile54.med"
        m=MEDFileUMesh()
        coo=DataArrayDouble(9) ; coo.iota(1.) ; coo.rearrange(3) ; coo.setInfoOnComponents(["aaa [b]","cc [dd]", "e [fff]"])
        m0=MEDCouplingUMesh("toto",2) ; m0.allocateCells(0)
        for i in xrange(7):
            m0.insertNextCell(NORM_TRI3,[1,2,1])
            pass
        for i in xrange(4):
            m0.insertNextCell(NORM_QUAD4,[1,1,2,0])
            pass
        for i in xrange(2):
            m0.insertNextCell(NORM_POLYGON,[0,0,1,1,2,2])
            pass
        m1=MEDCouplingUMesh("toto",1) ; m1.allocateCells(0) ; m1.insertNextCell(NORM_SEG2,[1,6]) ; m1.insertNextCell(NORM_SEG2,[7,3])
        m2=MEDCouplingUMesh("toto",0) ; m2.allocateCells(0) ; m2.insertNextCell(NORM_POINT1,[2]) ; m2.insertNextCell(NORM_POINT1,[6]) ; m2.insertNextCell(NORM_POINT1,[8])
        m0.setCoords(coo) ; m.setMeshAtLevel(0,m0)
        m1.setCoords(coo) ; m.setMeshAtLevel(-1,m1)
        m2.setCoords(coo) ; m.setMeshAtLevel(-2,m2)
        #
        mm=m.deepCpy()
        famCoo=DataArrayInt([0,2,0,3,2,0,-1,0,0,0,0,-1,3]) ; mm.setFamilyFieldArr(0,famCoo)
        da0=DataArrayInt([0,0,0]) ; mm.setFamilyFieldArr(1,da0)
        da1=DataArrayInt([0,3]) ; mm.setFamilyFieldArr(-1,da1)
        da2=DataArrayInt([0,0,0]) ; mm.setFamilyFieldArr(-2,da2)
        mm.setFamilyId("MyFam",2)
        mm.setFamilyId("MyOtherFam",3)
        mm.setFamilyId("MyOther-1",-1)
        mm.setFamiliesOnGroup("grp0",["MyOtherFam"])
        mm.setFamiliesOnGroup("grpA",["MyOther-1"])
        #
        daTest=DataArrayInt([1,3,4,6,9,10,12]) ; daTest.setName("grp1")
        mm.addGroup(0,daTest)
        self.assertTrue(mm.getGroupArr(0,daTest.getName()).isEqual(daTest))
        self.assertTrue(mm.getFamilyFieldAtLevel(0).isEqual(DataArrayInt([-6,2,-6,-8,2,-6,-5,-6,-6,-7,-7,-4,-8])))
        for lev,arr in [(1,da0),(-1,da1),(-2,da2)]:
            self.assertTrue(mm.getFamilyFieldAtLevel(lev).isEqual(arr))
            pass
        self.assertEqual(mm.getFamiliesNames(),('Family_-4','Family_-5','Family_-7','Family_-8','MyFam','MyOther-1','MyOtherFam'))
        self.assertEqual(mm.getGroupsNames(),('grp0','grp1','grpA'))
        self.assertEqual(mm.getFamilyNameGivenId(3),'MyOtherFam')
        self.assertEqual(mm.getFamilyNameGivenId(2),'MyFam')
        for famName,famId in [('Family_-4',-4),('Family_-5',-5),('Family_-7',-7),('Family_-8',-8)]:
            self.assertEqual(mm.getFamilyNameGivenId(famId),famName)
            pass
        self.assertEqual(mm.getFamiliesOnGroup("grp0"),('MyOtherFam','Family_-8'))
        da=DataArrayInt([3,12]) ; da.setName("grp0")
        self.assertTrue(mm.getGroupArr(0,"grp0").isEqual(da))
        da.setValues([1])
        self.assertTrue(mm.getGroupArr(-1,"grp0").isEqual(da))
        mm.write(fname,2)
        mm=MEDFileMesh.New(fname)
        self.assertTrue(mm.getGroupArr(0,daTest.getName()).isEqual(daTest))
        self.assertTrue(mm.getFamilyFieldAtLevel(0).isEqual(DataArrayInt([-6,2,-6,-8,2,-6,-5,-6,-6,-7,-7,-4,-8])))
        for lev,arr in [(1,da0),(-1,da1),(-2,da2)]:
            self.assertTrue(mm.getFamilyFieldAtLevel(lev).isEqual(arr))
            pass
        self.assertEqual(mm.getFamiliesNames(),('FAMILLE_ZERO','Family_-4','Family_-5','Family_-7','Family_-8','MyFam','MyOther-1','MyOtherFam'))
        self.assertEqual(mm.getGroupsNames(),('grp0','grp1','grpA'))
        self.assertEqual(mm.getFamilyNameGivenId(3),'MyOtherFam')
        self.assertEqual(mm.getFamilyNameGivenId(2),'MyFam')
        for famName,famId in [('Family_-4',-4),('Family_-5',-5),('Family_-7',-7),('Family_-8',-8)]:
            self.assertEqual(mm.getFamilyNameGivenId(famId),famName)
            pass
        self.assertEqual(mm.getFamiliesOnGroup("grp0"),('Family_-8','MyOtherFam'))
        da=DataArrayInt([3,12]) ; da.setName("grp0")
        self.assertTrue(mm.getGroupArr(0,"grp0").isEqual(da))
        da.setValues([1])
        self.assertTrue(mm.getGroupArr(-1,"grp0").isEqual(da))
        pass

    def testHeapMem1(self):
        a=DataArrayInt() ; aa=a.getHeapMemorySize()
        a.alloc(0,1)
        strMulFac=a.getHeapMemorySize()-aa ; del a ; del aa
        #
        m=MEDCouplingCMesh()
        arr=DataArrayDouble(10,1) ; arr.iota(0)
        m.setCoords(arr,arr)
        m=m.buildUnstructured()
        m.setName("mm")
        f=m.getMeasureField(ON_CELLS)
        self.assertIn(m.getHeapMemorySize(),xrange(3552-100,3552+100+4*strMulFac))
        self.assertIn(f.getHeapMemorySize(),xrange(4215-100,4215+100+8*strMulFac))
        #
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        self.assertIn(mm.getHeapMemorySize(),xrange(3889-100,3889+100+10*strMulFac))
        ff=MEDFileField1TS()
        ff.setFieldNoProfileSBT(f)
        self.assertIn(ff.getHeapMemorySize(),xrange(771-40,771+21+(4+1)*strMulFac))
        #
        fff=MEDFileFieldMultiTS()
        fff.appendFieldNoProfileSBT(f)
        self.assertIn(fff.getHeapMemorySize(),xrange(815-50,815+30+(6+2)*strMulFac))
        f.setTime(1.,0,-1)
        fff.appendFieldNoProfileSBT(f)
        self.assertIn(fff.getHeapMemorySize(),xrange(1594-90,1594+50+(10+1)*strMulFac))
        self.assertIn(fff[0,-1].getHeapMemorySize(),xrange(771-40,771+20+(4+1)*strMulFac))
        f2=f[:50]
        f2.setTime(2.,1,-1)
        pfl=DataArrayInt.Range(0,50,1) ; pfl.setName("pfl")
        fff.appendFieldProfile(f2,mm,0,pfl)
        self.assertIn(fff.getHeapMemorySize(),xrange(2348-130,2348+100+(10+2)*strMulFac))
        self.assertIn(fff.getProfile("pfl").getHeapMemorySize(),xrange(204-10,204+10+2*strMulFac))
        self.assertIn(fff[1,-1].getHeapMemorySize(),xrange(738-50,738+30+4*strMulFac))
        pass

    def testCurveLinearMesh1(self):
        fname="Pyfile55.med"
        mesh=MEDCouplingCurveLinearMesh();
        mesh.setTime(2.3,4,5);
        mesh.setTimeUnit("us");
        mesh.setName("Example of Cuve linear mesh");
        mesh.setDescription("buildCLMesh");
        a1=DataArrayDouble(3*20,1);
        a1.iota(7.) ; a1.rearrange(3);
        mesh.setCoords(a1);
        mesh.setNodeGridStructure([4,5]);
        mesh.checkCoherency();
        #
        m=MEDFileCurveLinearMesh()
        m.setMesh(mesh)
        d=DataArrayInt(20) ; d.iota(4)
        m.setFamilyFieldArr(1,d)
        d3=DataArrayInt(20) ; d3.iota(400)
        m.setRenumFieldArr(1,d3)
        d2=DataArrayInt(12) ; d2.iota(40)
        m.setFamilyFieldArr(0,d2)
        d4=DataArrayInt(21) ; d4.iota(4000)
        self.assertRaises(InterpKernelException,m.setRenumFieldArr,1,d4)
        d4.popBackSilent()
        m.setRenumFieldArr(1,d4)
        m.write(fname,2)
        #
        m1=MEDFileCurveLinearMesh(fname)
        mm=m1.getMesh()
        self.assertTrue(mm.isEqual(mesh,1e-12))
        self.assertEqual(mm.getSpaceDimension(),3)
        self.assertEqual(mm.getSpaceDimensionOnNodeStruct(),2)
        #
        m1=MEDFileMesh.New(fname)
        self.assertTrue(isinstance(m1,MEDFileCurveLinearMesh))
        self.assertTrue(isinstance(m1.getUnivName(),str))
        self.assertTrue(len(m1.getUnivName())!=0)
        self.assertTrue(m1.getMesh().isEqual(mesh,1e-12))
        pass

    def testParameters1(self):
        fname="Pyfile56.med"
        m=MEDCouplingCMesh() ; arr=DataArrayDouble([0.,1.2,3.5]) ; m.setCoords(arr,arr) ; m.setName("mesh")
        mm=MEDFileCMesh() ; mm.setMesh(m)
        ms=MEDFileMeshes() ; ms.pushMesh(mm)
        data=MEDFileData()
        p=MEDFileParameters()
        data.setParams(p) ; data.setMeshes(ms)
        pts=MEDFileParameterMultiTS()
        pts.setName("A") ; pts.setDescription("An example of parameter") ; pts.setTimeUnit("ms")
        pts.appendValue(1,2,3.4,567.89)
        pts.appendValue(2,3,5.6,999.123)
        pts2=pts.deepCpy() ; pts2.setName("B") ; pts2.setDescription("A second example")
        p.pushParam(pts) ; p.pushParam(pts2)
        data.write(fname,2)
        p2=MEDFileParameters(fname)
        self.assertTrue(p.isEqual(p2,1e-14)[0])
        self.assertAlmostEqual(p[1][1,2].getValue(),567.89,13)
        p3=p.deepCpy()
        pts4=pts2.deepCpy()
        pts3=pts2.deepCpy()
        self.assertTrue(pts3.isEqual(pts2,1e-14)[0])
        pts2.eraseTimeStepIds([0])
        self.assertTrue(not pts3.isEqual(pts2,1e-14)[0])
        del pts3[[3.4]]
        self.assertTrue(pts3.isEqual(pts2,1e-14)[0])
        self.assertRaises(InterpKernelException,p[1].__getitem__,(1,2))
        self.assertRaises(InterpKernelException,p["B"].__getitem__,(1,2))
        self.assertAlmostEqual(p[0][1,2].getValue(),567.89,13)
        self.assertAlmostEqual(p["A"][1,2].getValue(),567.89,13)
        p=p3
        self.assertTrue(p.isEqual(p2,1e-14)[0])
        self.assertTrue(p2["B"].isEqual(pts,1e-14)[0])
        self.assertTrue(not p2["B"].isEqual(pts2,1e-14)[0])
        self.assertAlmostEqual(p2[0][1,2].getValue(),567.89,13)
        self.assertEqual(p.getParamsNames(),('A','B'))
        ptsr=MEDFileParameterMultiTS(fname,"B")
        self.assertTrue(ptsr.isEqual(pts4,1e-14)[0])
        ptsr=MEDFileParameterMultiTS(fname)
        self.assertTrue(ptsr.isEqual(pts,1e-14)[0])
        p1tsr=MEDFileParameterDouble1TS(fname)
        self.assertEqual(p1tsr.getName(),"A")
        self.assertAlmostEqual(p1tsr.getValue(),567.89,13)
        p1tsr=MEDFileParameterDouble1TS(fname,"B")
        self.assertEqual(p1tsr.getName(),"B")
        self.assertAlmostEqual(p1tsr.getValue(),567.89,13)
        p1tsr=MEDFileParameterDouble1TS(fname,"B",2,3)
        self.assertEqual(p1tsr.getName(),"B")
        self.assertAlmostEqual(p1tsr.getValue(),999.123,13)
        data2=MEDFileData(fname)
        self.assertEqual(2,data2.getNumberOfParams())
        self.assertAlmostEqual(data2.getParams()["B"][1,2].getValue(),567.89,13)
        pass

    def testNamesOnCellAndNodesInMeshes1(self):
        fname="Pyfile58.med"
        fname2="Pyfile59.med"
        m=MEDLoaderDataForTest.build3DSurfMesh_1()
        m1=m.buildDescendingConnectivity()[0]
        m1.sortCellsInMEDFileFrmt()
        #
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        mm.setMeshAtLevel(-1,m1)
        namesCellL0=DataArrayAsciiChar(6,16)
        namesCellL0[:]=["CellL0#%.3d      "%(i) for i in xrange(6)]
        mm.setNameFieldAtLevel(0,namesCellL0)
        namesCellL1=DataArrayAsciiChar.Aggregate([namesCellL0,namesCellL0,namesCellL0.substr(2)])
        namesCellL1[:]=["CellLM1#%.3d     "%(i) for i in xrange(16)]
        mm.setNameFieldAtLevel(-1,namesCellL1)
        namesNodes=namesCellL1.substr(4,16)
        namesNodes[:]=["Node#%.3d        "%(i) for i in xrange(12)]
        mm.setNameFieldAtLevel(1,namesNodes)
        mm.write(fname,2)
        #
        mmr=MEDFileMesh.New(fname)
        self.assertTrue(mm.getNameFieldAtLevel(0).isEqual(DataArrayAsciiChar(["CellL0#%.3d      "%(i) for i in xrange(6)])))
        self.assertTrue(mm.getNameFieldAtLevel(-1).isEqual(DataArrayAsciiChar(["CellLM1#%.3d     "%(i) for i in xrange(16)])))
        self.assertTrue(mm.getNameFieldAtLevel(1).isEqual(DataArrayAsciiChar(["Node#%.3d        "%(i) for i in xrange(12)])))
        self.assertTrue(mm.isEqual(mmr,1e-12)[0])
        mmr.getNameFieldAtLevel(1).setIJ(0,0,'M')
        self.assertTrue(not mm.isEqual(mmr,1e-12)[0])
        mmr.getNameFieldAtLevel(1).setIJ(0,0,'N')
        self.assertTrue(mm.isEqual(mmr,1e-12)[0])
        mmCpy=mm.deepCpy()
        self.assertTrue(mm.isEqual(mmCpy,1e-12)[0])
        # remove names on nodes
        mmCpy.setNameFieldAtLevel(1,None)
        self.assertTrue(not mm.isEqual(mmCpy,1e-12)[0])
        mm.setNameFieldAtLevel(1,None)
        self.assertTrue(mm.isEqual(mmCpy,1e-12)[0])
        mm.setNameFieldAtLevel(-1,None)
        mm.write(fname,2)
        mmr=MEDFileMesh.New(fname)
        self.assertEqual(mmr.getNameFieldAtLevel(1),None)
        self.assertTrue(mmr.getNameFieldAtLevel(0).isEqual(DataArrayAsciiChar(["CellL0#%.3d      "%(i) for i in xrange(6)])))
        self.assertEqual(mmr.getNameFieldAtLevel(-1),None)
        #
        c=MEDCouplingCMesh()
        arr=DataArrayDouble([0.,1.1,2.3])
        c.setCoords(arr,arr)
        c.setName("cmesh")
        cc=MEDFileCMesh()
        cc.setMesh(c)
        cc.setNameFieldAtLevel(0,DataArrayAsciiChar(["Cell#%.3d        "%(i) for i in xrange(4)]))
        cc.setNameFieldAtLevel(1,DataArrayAsciiChar(["Node#%.3d        "%(i) for i in xrange(9)]))
        cc.write(fname2,2)
        ccr=MEDFileMesh.New(fname2)
        self.assertTrue(ccr.getNameFieldAtLevel(0).isEqual(DataArrayAsciiChar(["Cell#%.3d        "%(i) for i in xrange(4)])))
        self.assertTrue(ccr.getNameFieldAtLevel(1).isEqual(DataArrayAsciiChar(["Node#%.3d        "%(i) for i in xrange(9)])))
        self.assertTrue(cc.isEqual(ccr,1e-12)[0])
        ccr.getNameFieldAtLevel(1).setIJ(0,0,'M')
        self.assertTrue(not cc.isEqual(ccr,1e-12)[0])
        ccr.getNameFieldAtLevel(1).setIJ(0,0,'N')
        self.assertTrue(cc.isEqual(ccr,1e-12)[0])
        ccCpy=cc.deepCpy()
        self.assertTrue(cc.isEqual(ccCpy,1e-12)[0])
        pass

    def testToExportInExamples1(self):
        m=MEDCouplingCMesh()
        arr=DataArrayDouble([0.,1.,2.,3.,4.])
        m.setCoords(arr,arr)
        m=m.buildUnstructured() ; m.setName("mesh")
        grp1=DataArrayInt([0,1,2,4,5,6,8,9,10,12,13,14]) ; grp1.setName("grp1")
        grp2=DataArrayInt([3,7,11,15]) ; grp2.setName("grp2")
        m2=m.computeSkin()
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        mm.setMeshAtLevel(-1,m2)
        mm.setGroupsAtLevel(0,[grp1,grp2])
        mm.write("example.med",2)
        #
        m0=mm.getMeshAtLevel(0)
        m1=mm.getMeshAtLevel(-1)
        grp1=mm.getGroupArr(0,"grp1")
        grp2=mm.getGroupArr(0,"grp2")
        grps=[grp1,grp2]
        whichGrp=DataArrayInt(m0.getNumberOfCells())
        whichGrp.fillWithValue(-1)
        for grpId,grp in enumerate(grps):
            whichGrp[grp]=grpId
            pass
        a,b,bI,c,cI=m0.buildDescendingConnectivity()
        e,f=a.areCellsIncludedIn(m1,2)
        self.assertTrue(e)
        c2,c2I=MEDCouplingUMesh.ExtractFromIndexedArrays(f,c,cI)
        self.assertTrue(c2I.deltaShiftIndex().isUniform(1))
        c2.transformWithIndArr(whichGrp)
        splitOfM1=len(grps)*[None]
        for grpId,grp in enumerate(grps):
            tmp=c2.getIdsEqual(grpId)
            splitOfM1[grpId]=tmp
            pass
        splitOfM1[0].isEqual(DataArrayInt([0,1,2,3,6,8,10,11,12,13]))
        splitOfM1[1].isEqual(DataArrayInt([4,5,7,9,14,15]))
        pass

    def testBugCorrection1(self):
        fs=MEDFileFields()
        fs.resize(3)
        self.assertEqual(fs[0],None)
        self.assertEqual(3,len(fs))
        pass

    def testCompareMEDFilesContainingOnlyFieldsOnCell1(self):
        f1Name="Pyfile60.med"
        f2Name="Pyfile61.med"
        d1=MEDLoaderDataForTest.buildACompleteMEDDataStructureWithFieldsOnCells_1()
        d1.write(f1Name,2)
        d2=MEDLoaderDataForTest.buildACompleteMEDDataStructureWithFieldsOnCells_1()
        d2.write(f2Name,2)
        # reading and compare
        d1=MEDFileData(f1Name) ; d2=MEDFileData(f2Name)
        for mn in d1.getMeshes().getMeshesNames():
            m1=d1.getMeshes()[mn]
            m2=d2.getMeshes()[mn]
            for lev in m1.getNonEmptyLevels():
                grpsNames=m1.getGroupsOnSpecifiedLev(lev)
                for grpName in grpsNames:
                    self.assertTrue(m1.getGroupArr(lev,grpName).isEqual(m2.getGroupArr(lev,grpName))) # compare groups
                    pass
                pass
            pass
        for fieldn in d1.getFields().getFieldsNames():
            f1=d1.getFields()[fieldn]
            f2=d2.getFields()[fieldn]
            for it,order,tim in f1.getTimeSteps():
                f1t=f1[it,order]
                f2t=f2[it,order]
                if len(f1t.getPflsReallyUsed())!=0:
                    # profile case
                    for lev in f1t.getNonEmptyLevels()[1]:
                        arr1,pfl1=f1t.getFieldWithProfile(ON_CELLS,lev,m1)
                        arr2,pfl2=f2t.getFieldWithProfile(ON_CELLS,lev,m2)
                        self.assertTrue(pfl1.isEqual(pfl2))
                        self.assertTrue(arr1.isEqual(arr2,1e-10))
                        pass
                    pass
                else:
                    # no profile case
                    for lev in f1t.getNonEmptyLevels()[1]:
                        f1mc=f1t.getFieldOnMeshAtLevel(ON_CELLS,lev,m1)
                        f2mc=f2t.getFieldOnMeshAtLevel(ON_CELLS,lev,m2)
                        self.assertTrue(f1mc.isEqual(f2mc,1e-10,1e-10))
                        pass
                    pass
                pass
            pass
        pass

    def testNonRegBugNormalizeFamIdsMEDFile1(self):
        m=MEDCouplingCMesh()
        arr=DataArrayDouble([0.,1.,2.,3.,4.])
        m.setCoords(arr,arr,arr)
        m=m.buildUnstructured()
        m2=m.buildDescendingConnectivity()[0]
        m.setName("mesh")
        g1=DataArrayInt([0,1,2,3]) ; g1.setName("g1")
        g2=DataArrayInt([2,3,5,6]) ; g2.setName("g2")
        g1Face=DataArrayInt([20,21,22,23]) ; g1Face.setName("g1Face")
        g2Face=DataArrayInt([22,23,25,26]) ; g2Face.setName("g2Face")
        g1Node=DataArrayInt([10,11,12,13]) ; g1Node.setName("g1Node")
        g2Node=DataArrayInt([12,13,15,16]) ; g2Node.setName("g2Node")
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        mm.setGroupsAtLevel(0,[g1,g2])
        s1=set(mm.getFamiliesOnGroup("g1")) ; s2=set(mm.getFamiliesOnGroup("g2"))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g1"),(0,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g2"),(0,))
        mm.normalizeFamIdsMEDFile()
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g1"),(0,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g2"),(0,))
        self.assertTrue(mm.getGroupArr(0,"g1").isEqual(g1))
        self.assertTrue(mm.getGroupArr(0,"g2").isEqual(g2))
        self.assertEqual(s1,set(mm.getFamiliesOnGroup("g1")))
        self.assertEqual(s2,set(mm.getFamiliesOnGroup("g2")))
        for g in mm.getGroupsOnSpecifiedLev(0):
            for f in mm.getFamiliesIdsOnGroup(g):
                self.assertTrue(f<0)
                pass
            pass
        #
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        mm.setMeshAtLevel(-1,m2)
        mm.setGroupsAtLevel(0,[g1,g2])
        mm.setGroupsAtLevel(-1,[g1Face,g2Face])
        s1=set(mm.getFamiliesOnGroup("g1")) ; s2=set(mm.getFamiliesOnGroup("g2"))
        s3=set(mm.getFamiliesOnGroup("g1Face")) ; s4=set(mm.getFamiliesOnGroup("g2Face"))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g1"),(0,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g2"),(0,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g1Face"),(-1,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g2Face"),(-1,))
        mm.normalizeFamIdsMEDFile()
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g1"),(0,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g2"),(0,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g1Face"),(-1,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g2Face"),(-1,))
        self.assertTrue(mm.getGroupArr(0,"g1").isEqual(g1))
        self.assertTrue(mm.getGroupArr(0,"g2").isEqual(g2))
        self.assertTrue(mm.getGroupArr(-1,"g1Face").isEqual(g1Face))
        self.assertTrue(mm.getGroupArr(-1,"g2Face").isEqual(g2Face))
        self.assertEqual(s1,set(mm.getFamiliesOnGroup("g1")))
        self.assertEqual(s2,set(mm.getFamiliesOnGroup("g2")))
        self.assertEqual(s3,set(mm.getFamiliesOnGroup("g1Face")))
        self.assertEqual(s4,set(mm.getFamiliesOnGroup("g2Face")))
        for lev in [0,-1]:
            for g in mm.getGroupsOnSpecifiedLev(lev):
                for f in mm.getFamiliesIdsOnGroup(g):
                    self.assertTrue(f<0)
                    pass
                pass
            pass
         #
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        mm.setMeshAtLevel(-1,m2)
        mm.setGroupsAtLevel(0,[g1,g2])
        mm.setGroupsAtLevel(-1,[g1Face,g2Face])
        mm.setGroupsAtLevel(1,[g1Node,g2Node])
        s1=set(mm.getFamiliesOnGroup("g1")) ; s2=set(mm.getFamiliesOnGroup("g2"))
        s3=set(mm.getFamiliesOnGroup("g1Face")) ; s4=set(mm.getFamiliesOnGroup("g2Face"))
        s5=set(mm.getFamiliesOnGroup("g1Node")) ; s6=set(mm.getFamiliesOnGroup("g2Node"))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g1"),(0,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g2"),(0,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g1Face"),(-1,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g2Face"),(-1,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g1Node"),(1,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g2Node"),(1,))
        mm.normalizeFamIdsMEDFile()
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g1"),(0,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g2"),(0,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g1Face"),(-1,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g2Face"),(-1,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g1Node"),(1,))
        self.assertEqual(mm.getGrpNonEmptyLevelsExt("g2Node"),(1,))
        self.assertTrue(mm.getGroupArr(0,"g1").isEqual(g1))
        self.assertTrue(mm.getGroupArr(0,"g2").isEqual(g2))
        self.assertTrue(mm.getGroupArr(-1,"g1Face").isEqual(g1Face))
        self.assertTrue(mm.getGroupArr(-1,"g2Face").isEqual(g2Face))
        self.assertTrue(mm.getGroupArr(1,"g1Node").isEqual(g1Node))
        self.assertTrue(mm.getGroupArr(1,"g2Node").isEqual(g2Node))
        self.assertEqual(s1,set(mm.getFamiliesOnGroup("g1")))
        self.assertEqual(s2,set(mm.getFamiliesOnGroup("g2")))
        self.assertEqual(s3,set(mm.getFamiliesOnGroup("g1Face")))
        self.assertEqual(s4,set(mm.getFamiliesOnGroup("g2Face")))
        self.assertEqual(s5,set(mm.getFamiliesOnGroup("g1Node")))
        self.assertEqual(s6,set(mm.getFamiliesOnGroup("g2Node")))
        for lev in [0,-1]:
            for g in mm.getGroupsOnSpecifiedLev(lev):
                for f in mm.getFamiliesIdsOnGroup(g):
                    self.assertTrue(f<0)
                    pass
                pass
            pass
        for g in mm.getGroupsOnSpecifiedLev(1):
            for f in mm.getFamiliesIdsOnGroup(g):
                self.assertTrue(f>0)
                pass
            pass
        pass
    
    def testNonRegressionMantis22212ChangeGrpName(self):
        fileName="Pyfile62.med"
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
        #
        mm0=MEDFileMesh.New(fileName)
        mm1=MEDFileMesh.New(fileName)
        groupNamesIni=MEDLoader.GetMeshGroupsNames(fileName,"ma")
        for name in groupNamesIni:
            mm1.changeGroupName(name,name+'N')
            pass
        mm1.write(fileName,2)
        del mm1
        #
        mm2=MEDFileMesh.New(fileName)
        for name in groupNamesIni:
            for lev in mm0.getGrpNonEmptyLevelsExt(name):
                arr0=mm0.getGroupArr(lev,name)
                arr2=mm2.getGroupArr(lev,name+'N')
                arr0.setName(name+'N')
                self.assertTrue(arr0.isEqual(arr2))
                pass
            pass
        pass

    def testInt32InMEDFileFieldStar1(self):
        fname="Pyfile63.med"
        f1=MEDLoaderDataForTest.buildVecFieldOnCells_1();
        arr=f1.getArray().convertToIntArr()
        f1.setArray(None)
        m1=f1.getMesh()
        mm1=MEDFileUMesh.New()
        mm1.setCoords(m1.getCoords())
        mm1.setMeshAtLevel(0,m1)
        mm1.setName(m1.getName())
        mm1.write(fname,2)
        ff1=MEDFileIntField1TS()
        ff1.setFieldNoProfileSBT(f1,arr)
        a,b=ff1.getFieldOnMeshAtLevel(0,ON_CELLS,mm1)
        self.assertEqual(b.getInfoOnComponents(),['power [MW/m^3]','density [g/cm^3]','temperature [K]'])
        self.assertTrue(b.isEqual(arr))
        self.assertTrue(a.isEqual(f1,1e-12,1e-12))
        ff1.write(fname,0)
        ff2=MEDFileAnyTypeField1TS.New(fname)
        self.assertEqual(ff2.getName(),"VectorFieldOnCells")
        self.assertEqual(ff2.getTime(),[0,1,2.0])
        self.assertTrue(isinstance(ff2,MEDFileIntField1TS))
        a,b=ff1.getFieldOnMeshAtLevel(0,ON_CELLS,mm1)
        self.assertEqual(b.getInfoOnComponents(),['power [MW/m^3]','density [g/cm^3]','temperature [K]'])
        self.assertTrue(b.isEqual(arr))
        self.assertTrue(a.isEqual(f1,1e-12,1e-12))
        ff2.setTime(1,2,3.)
        c=ff2.getUndergroundDataArray() ; c*=2
        ff2.write(fname,0) # 2 time steps in 
        ffs1=MEDFileAnyTypeFieldMultiTS.New(fname,"VectorFieldOnCells")
        self.assertEqual(ffs1.getTimeSteps(),[(0, 1, 2.0), (1, 2, 3.0)])
        self.assertEqual(len(ffs1),2)
        self.assertTrue(isinstance(ffs1,MEDFileIntFieldMultiTS))
        a,b=ffs1[2.].getFieldOnMeshAtLevel(0,ON_CELLS,mm1)
        self.assertTrue(b.isEqual(arr))
        self.assertTrue(a.isEqual(f1,1e-12,1e-12))
        a,b=ffs1[2.].getFieldOnMeshAtLevel(0,ON_CELLS,mm1)
        self.assertTrue(b.isEqual(arr))
        self.assertTrue(a.isEqual(f1,1e-12,1e-12))
        it=ffs1.__iter__() ; it.next() ; ff2bis=it.next()
        a,b=ff2bis.getFieldOnMeshAtLevel(0,ON_CELLS,mm1)
        self.assertTrue(b.isEqual(2*arr))
        f1.setTime(3.,1,2)
        self.assertTrue(a.isEqual(f1,1e-12,1e-12))
        bc=DataArrayInt(6,3) ; bc[:]=0 ; bc.setInfoOnComponents(['power [MW/m^3]','density [g/cm^3]','temperature [K]'])
        for it in ffs1:
            a,b=it.getFieldOnMeshAtLevel(0,ON_CELLS,mm1)
            bc+=b
            pass
        self.assertTrue(bc.isEqual(3*arr))
        nf1=MEDCouplingFieldDouble(ON_NODES)
        nf1.setTime(9.,10,-1)
        nf1.setMesh(f1.getMesh())
        narr=DataArrayInt(12,2) ; narr.setInfoOnComponents(["aa [u1]","bbbvv [ppp]"]) ; narr[:,0]=range(12) ; narr[:,1]=2*narr[:,0]
        nf1.setName("VectorFieldOnNodes")
        nff1=MEDFileIntField1TS.New()
        nff1.setFieldNoProfileSBT(nf1,narr)
        self.assertEqual(nff1.getInfo(),('aa [u1]','bbbvv [ppp]'))
        self.assertEqual(nff1.getTime(),[10,-1,9.0])
        nff1.write(fname,0)
        #
        nf2=MEDCouplingFieldDouble(ON_NODES)
        nf2.setTime(19.,20,-11)
        nf2.setMesh(f1.getMesh())
        narr2=DataArrayInt(8,2) ; narr.setInfoOnComponents(["aapfl [u1]","bbbvvpfl [ppp]"]) ; narr2[:,0]=range(8) ; narr2[:,0]+=10  ; narr2[:,1]=3*narr2[:,0]
        nf2.setName("VectorFieldOnNodesPfl") ; narr2.setName(nf2.getName())
        nff2=MEDFileIntField1TS.New()
        npfl=DataArrayInt([1,2,4,5,6,7,10,11]) ; npfl.setName("npfl")
        nff2.setFieldProfile(nf2,narr2,mm1,0,npfl)
        nff2.getFieldWithProfile(ON_NODES,0,mm1)
        a,b=nff2.getFieldWithProfile(ON_NODES,0,mm1) ; b.setName(npfl.getName())
        self.assertTrue(b.isEqual(npfl))
        self.assertTrue(a.isEqual(narr2))
        nff2.write(fname,0)
        nff2bis=MEDFileIntField1TS(fname,"VectorFieldOnNodesPfl")
        a,b=nff2bis.getFieldWithProfile(ON_NODES,0,mm1) ; b.setName(npfl.getName())
        self.assertTrue(b.isEqual(npfl))
        self.assertTrue(a.isEqual(narr2))
        #
        nf3=MEDCouplingFieldDouble(ON_NODES)
        nf3.setName("VectorFieldOnNodesDouble")
        nf3.setTime(29.,30,-21)
        nf3.setMesh(f1.getMesh())
        nf3.setArray(f1.getMesh().getCoords())
        nff3=MEDFileField1TS.New()
        nff3.setFieldNoProfileSBT(nf3)
        nff3.write(fname,0)
        fs=MEDFileFields(fname)
        self.assertEqual(len(fs),4)
        ffs=[it for it in fs]
        self.assertTrue(isinstance(ffs[0],MEDFileIntFieldMultiTS))
        self.assertTrue(isinstance(ffs[1],MEDFileIntFieldMultiTS))
        self.assertTrue(isinstance(ffs[2],MEDFileFieldMultiTS))
        self.assertTrue(isinstance(ffs[3],MEDFileIntFieldMultiTS))
        #
        self.assertTrue(fs["VectorFieldOnCells"][0].getUndergroundDataArray().isEqualWithoutConsideringStr(arr))
        self.assertTrue(fs["VectorFieldOnCells"][1,2].getUndergroundDataArray().isEqualWithoutConsideringStr(2*arr))
        self.assertTrue(fs["VectorFieldOnNodesPfl"][0].getUndergroundDataArray().isEqualWithoutConsideringStr(narr2))
        self.assertTrue(fs["VectorFieldOnNodes"][9.].getUndergroundDataArray().isEqualWithoutConsideringStr(narr))
        self.assertTrue(fs["VectorFieldOnNodesDouble"][29.].getUndergroundDataArray().isEqualWithoutConsideringStr(f1.getMesh().getCoords(),1e-12))
        #
        nf3_read=MEDFileFieldMultiTS(fname,"VectorFieldOnNodesDouble")
        self.assertTrue(nf3_read[29.].getUndergroundDataArray().isEqualWithoutConsideringStr(f1.getMesh().getCoords(),1e-12))
        self.assertRaises(InterpKernelException,MEDFileIntFieldMultiTS.New,fname,"VectorFieldOnNodesDouble")# exception because trying to read a double field with int instance
        self.assertRaises(InterpKernelException,MEDFileFieldMultiTS.New,fname,"VectorFieldOnNodes")# exception because trying to read a int field with double instance
        MEDFileField1TS.New(fname,"VectorFieldOnNodesDouble",30,-21)
        self.assertRaises(InterpKernelException,MEDFileIntField1TS.New,fname,"VectorFieldOnNodesDouble",30,-21)# exception because trying to read a double field with int instance
        MEDFileIntField1TS.New(fname,"VectorFieldOnNodes",10,-1)
        self.assertRaises(InterpKernelException,MEDFileField1TS.New,fname,"VectorFieldOnNodes",10,-1)# exception because trying to read a double field with int instance
        #
        self.assertEqual(fs.getMeshesNames(),('3DSurfMesh_1','3DSurfMesh_1','3DSurfMesh_1','3DSurfMesh_1'))
        self.assertTrue(fs.changeMeshNames([('3DSurfMesh_1','3DSurfMesh')]))
        self.assertEqual(fs.getMeshesNames(),('3DSurfMesh','3DSurfMesh','3DSurfMesh','3DSurfMesh'))
        self.assertTrue(not fs.changeMeshNames([('3DSurfMesh_1','3DSurfMesh')]))
        pass

    def testMEDFileFields1(self):
        fname="Pyfile64.med"
        f1=MEDCouplingFieldDouble(ON_NODES)
        f1.setTime(0.001,0,-1) ; f1.setTimeUnit("us")
        c=DataArrayDouble(12) ; c.iota(); m=MEDCouplingCMesh() ; m.setCoordsAt(0,c) ; m.setName("mesh")
        mm=MEDFileCMesh() ; mm.setMesh(m) ; mm.write(fname,2)
        f1.setMesh(m)
        arr=DataArrayDouble(12,2) ; arr.setInfoOnComponents(["aa [u1]","bbbvv [ppp]"]) ; arr[:,0]=range(12) ; arr[:,1]=2*arr[:,0]
        f1.setArray(arr)
        f1.setName("Field1")
        ff1=MEDFileField1TS.New()
        ff1.setFieldNoProfileSBT(f1)
        self.assertEqual(ff1.getDtUnit(),"us")
        ff1.write(fname,0)
        f1.setTime(1.001,1,-1) ; ff1=MEDFileField1TS.New() ; ff1.setFieldNoProfileSBT(f1) ; ff1.write(fname,0)
        f1.setTime(2.001,2,-1) ; ff1=MEDFileField1TS.New() ; ff1.setFieldNoProfileSBT(f1) ; ff1.write(fname,0)
        #
        self.assertEqual(MEDFileFields(fname).getCommonIterations(),([(0,-1),(1,-1),(2,-1)],False))
        ff1s=MEDFileFieldMultiTS(fname,"Field1")
        ff1s.setName("Field2")
        ff1s.write(fname,0)
        self.assertEqual(MEDFileFields(fname).getCommonIterations(),([(0,-1),(1,-1),(2,-1)],False))
        f1.setTime(3.001,3,-1) ; ff1=MEDFileField1TS.New() ; ff1.setFieldNoProfileSBT(f1) ; ff1.write(fname,0)
        self.assertEqual(MEDFileFields(fname).getCommonIterations(),([(0,-1),(1,-1),(2,-1)],True))
        self.assertEqual(MEDFileFields(fname).partOfThisLyingOnSpecifiedTimeSteps([(1,-1)]).getCommonIterations(),([(1,-1)],False))
        self.assertEqual(MEDFileFields(fname).partOfThisNotLyingOnSpecifiedTimeSteps([(1,-1)]).getCommonIterations(),([(0,-1),(2,-1)],True))
        f1.setName("Field2") ; f1.setTime(3.001,3,-1) ; ff1=MEDFileField1TS.New() ; ff1.setFieldNoProfileSBT(f1) ; ff1.write(fname,0)
        self.assertEqual(MEDFileFields(fname).getCommonIterations(),([(0,-1),(1,-1),(2,-1),(3,-1)],False))
        self.assertEqual(MEDFileFields(fname)[1].getDtUnit(),"us")
        pass

    # Multi time steps and multi fields management without Globals (profiles, locs) aspects
    def testMEDFileFields2(self):
        fname="Pyfile65.med"
        # to check that all is initialize 
        MEDFileField1TS().__str__()
        MEDFileFieldMultiTS().__str__()
        # building a mesh containing 4 tri3 + 5 quad4
        tri=MEDCouplingUMesh("tri",2)
        tri.allocateCells() ; tri.insertNextCell(NORM_TRI3,[0,1,2])
        tri.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,0.)]))
        tris=[tri.deepCpy() for i in xrange(4)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads=[quad.deepCpy() for i in xrange(5)]
        for i,elt in enumerate(quads): elt.translate([5+i,0])
        quads=MEDCouplingUMesh.MergeUMeshes(quads)
        m=MEDCouplingUMesh.MergeUMeshes(tris,quads)
        m.setName("mesh") ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        #
        fmts0_0=MEDFileFieldMultiTS()
        fmts0_1=MEDFileFieldMultiTS()
        # time steps
        for i in xrange(10):
            infos1=["aa [bb]","ccc [ddd]"] ; name1="1stField"
            d=DataArrayDouble(18) ; d.iota(i*10) ; d.rearrange(2) ; d.setInfoOnComponents(infos1)
            f=MEDCouplingFieldDouble(ON_CELLS) ; f.setName(name1) ; f.setArray(d) ; f.setMesh(m)
            f.setTime(float(i+1)+0.1,i+1,-i-1)
            fmts0_0.appendFieldNoProfileSBT(f)
            f1ts=MEDFileField1TS() ; f1ts.setFieldNoProfileSBT(f) ; fmts0_1.pushBackTimeStep(f1ts)
            self.assertEqual(fmts0_1.getName(),name1)
            self.assertEqual(fmts0_0.getInfo(),('aa [bb]','ccc [ddd]'))
            self.assertEqual(fmts0_1.getInfo(),('aa [bb]','ccc [ddd]'))
            if i>1:
                # components names have been modified to generate errors
                d.setInfoOnComponents(['aa [bb]','eee [dd]'])
                self.assertRaises(InterpKernelException,fmts0_0.appendFieldNoProfileSBT,f)
                self.assertRaises(InterpKernelException,f1ts.setInfo,['aa [bb]'])#throw because mismatch of number of components
                f1ts.setInfo(['aa [bb]','eee [dd]'])
                self.assertRaises(InterpKernelException,fmts0_1.pushBackTimeStep,f1ts)
                pass
            # add a mismatch of nb of compos
            pass
        fmts0_2=fmts0_0.deepCpy()
        fmts0_3=fmts0_0.deepCpy()
        fmts0_4=fmts0_0.deepCpy()
        fmts0_5=fmts0_0.shallowCpy()
        self.assertTrue(len(fmts0_0)==10 and len(fmts0_1)==10 and len(fmts0_2)==10 and len(fmts0_3)==10 and len(fmts0_4)==10 and len(fmts0_5)==10)
        del fmts0_2[::2]
        self.assertTrue(len(fmts0_2)==5 and fmts0_2.getIterations()==[(2,-2),(4,-4),(6,-6),(8,-8),(10,-10)])
        del fmts0_3[[1.1,(6,-6),9]]
        self.assertTrue(len(fmts0_3)==7 and fmts0_3.getIterations()==[(2,-2),(3,-3),(4,-4),(5,-5),(7,-7),(8,-8),(9,-9)])
        fmts0_6=fmts0_4[[1.1,(6,-6),8]]
        self.assertTrue(isinstance(fmts0_6,MEDFileFieldMultiTS))
        self.assertTrue(len(fmts0_6)==3 and fmts0_6.getIterations()==[(1,-1),(6,-6),(9,-9)])
        fmts0_7=fmts0_4[::-3]
        self.assertTrue(isinstance(fmts0_7,MEDFileFieldMultiTS))
        self.assertTrue(len(fmts0_7)==4 and fmts0_7.getIterations()==[(10,-10),(7,-7),(4,-4),(1,-1)])
        #
        fs0=MEDFileFields()
        fs0.pushField(fmts0_0)
        fmts0_2.setName("2ndField") ; fs0.pushField(fmts0_2)
        fmts0_3.setName("3rdField") ; fs0.pushField(fmts0_3)
        fmts0_4.setName("4thField") ; fs0.pushField(fmts0_4)
        self.assertTrue(len(fs0)==4 and fs0.getFieldsNames()==('1stField','2ndField','3rdField','4thField'))
        fs0.write(fname,2)
        fs0=MEDFileFields(fname)
        self.assertEqual(fs0.getCommonIterations(),([(2,-2),(4,-4),(8,-8)],True))
        fs1=fs0.partOfThisLyingOnSpecifiedTimeSteps(fs0.getCommonIterations()[0])
        self.assertTrue(fs1.getFieldsNames()==('1stField','2ndField','3rdField','4thField') and fs1.getCommonIterations()==([(2,-2),(4,-4),(8,-8)],False))
        del fs1[["2ndField",3]]
        self.assertTrue(fs1.getFieldsNames()==('1stField','3rdField') and fs1.getCommonIterations()==([(2,-2),(4,-4),(8,-8)],False))
        fs2=fs0[[0,"4thField"]]
        self.assertTrue(isinstance(fs2,MEDFileFields))
        self.assertEqual(fs2.getFieldsNames(),('1stField','4thField'))
        #
        mm=MEDFileUMesh() ; mm.setMeshAtLevel(0,m) ; mm.write(fname,0)
        pass

    # Multi time steps and multi fields management with Globals (profiles, locs) aspects
    def testMEDFileFields3(self):
        fname="Pyfile66.med"
        # building a mesh containing 4 tri3 + 5 quad4
        tri=MEDCouplingUMesh("tri",2)
        tri.allocateCells() ; tri.insertNextCell(NORM_TRI3,[0,1,2])
        tri.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,0.)]))
        tris=[tri.deepCpy() for i in xrange(4)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads=[quad.deepCpy() for i in xrange(5)]
        for i,elt in enumerate(quads): elt.translate([5+i,0])
        quads=MEDCouplingUMesh.MergeUMeshes(quads)
        m=MEDCouplingUMesh.MergeUMeshes(tris,quads)
        m.setName("mesh") ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        #
        mm=MEDFileUMesh() ; mm.setMeshAtLevel(0,m) ; mm.write(fname,2)
        #
        pfl=DataArrayInt([0,1,2,3,4,5,6]) ; pfl.setName("pfl")
        pfl2=DataArrayInt([0,1,2,3,4,5,6,8]) ; pfl2.setName("pfl2")
        fmts0_0=MEDFileFieldMultiTS()
        fmts0_1=MEDFileFieldMultiTS()
        # time steps
        for i in xrange(10):
            infos1=["aa [bb]","ccc [ddd]"] ; name1="1stField"
            d=DataArrayDouble(14) ; d.iota(i*10) ; d.rearrange(2) ; d.setInfoOnComponents(infos1)
            f=MEDCouplingFieldDouble(ON_CELLS) ; f.setName(name1) ; f.setArray(d) ; f.setMesh(m)
            f.setTime(float(i+1)+0.1,i+1,-i-1)
            fmts0_0.appendFieldProfile(f,mm,0,pfl)
            f1ts=MEDFileField1TS() ; f1ts.setFieldProfile(f,mm,0,pfl) ; fmts0_1.pushBackTimeStep(f1ts)
            self.assertEqual(fmts0_0.getInfo(),('aa [bb]','ccc [ddd]'))
            self.assertEqual(fmts0_1.getInfo(),('aa [bb]','ccc [ddd]'))
            pass
        #
        self.assertEqual(fmts0_0.getPfls(),10*('pfl_NORM_QUAD4',))
        self.assertEqual(fmts0_1.getPfls(),('pfl_NORM_QUAD4',))
        fmts0_0.zipPflsNames()
        self.assertEqual(fmts0_0.getPfls(),('pfl_NORM_QUAD4',))
        self.assertTrue(fmts0_1.getProfile("pfl_NORM_QUAD4").isEqual(fmts0_0.getProfile("pfl_NORM_QUAD4")))
        fmts0_2=fmts0_0.deepCpy()
        fmts0_3=fmts0_0.deepCpy()
        fmts0_4=fmts0_0.deepCpy()
        fs0=MEDFileFields()
        fs0.pushField(fmts0_0)
        fmts0_2.setName("2ndField") ; fs0.pushField(fmts0_2)
        fmts0_3.setName("3rdField") ; fs0.pushField(fmts0_3)
        fmts0_4.setName("4thField") ; fs0.pushField(fmts0_4)
        self.assertEqual(fs0.getPfls(),('pfl_NORM_QUAD4',))
        #
        fmts0_5=MEDFileFieldMultiTS()
        for i in xrange(7):
            infos1=["aa [bb]","ccc [ddd]"] ; name1="1stField"
            d=DataArrayDouble(16) ; d.iota(i*10) ; d.rearrange(2) ; d.setInfoOnComponents(infos1)
            f=MEDCouplingFieldDouble(ON_CELLS) ; f.setName(name1) ; f.setArray(d) ; f.setMesh(m)
            f.setTime(float(i+1)+0.1,i+1,-i-1)
            f1ts=MEDFileField1TS() ; f1ts.setFieldProfile(f,mm,0,pfl2) ; fmts0_5.pushBackTimeStep(f1ts)
            pass
        fmts0_5.setName("5thField") ; fs0.pushField(fmts0_5)
        self.assertEqual(fs0.getPfls(),('pfl_NORM_QUAD4','pfl2_NORM_QUAD4'))
        fs0.checkGlobsCoherency()
        fs0.write(fname,0)
        pass
    
    def testSplitComponents1(self):
        fname="Pyfile67.med"
        # building a mesh containing 4 tri3 + 5 quad4
        tri=MEDCouplingUMesh("tri",2)
        tri.allocateCells() ; tri.insertNextCell(NORM_TRI3,[0,1,2])
        tri.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,0.)]))
        tris=[tri.deepCpy() for i in xrange(4)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads=[quad.deepCpy() for i in xrange(5)]
        for i,elt in enumerate(quads): elt.translate([5+i,0])
        quads=MEDCouplingUMesh.MergeUMeshes(quads)
        m=MEDCouplingUMesh.MergeUMeshes(tris,quads)
        m.setName("mesh") ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        #
        mm=MEDFileUMesh() ; mm.setMeshAtLevel(0,m) ; mm.write(fname,2)
        #
        pfl=DataArrayInt([0,1,2,3,4,5,6]) ; pfl.setName("pfl")
        pfl2=DataArrayInt([0,1,2,3,4,5,6,8]) ; pfl2.setName("pfl2")
        fs=MEDFileFields()
        fmts0_1=MEDFileFieldMultiTS()
        # time steps
        infos1=['aa [bb]','ccc [ddd]',"ZZZZ [MW*s]"]
        for i in xrange(10):
            name1="1stField"
            d=DataArrayDouble(21) ; d.iota(i*10) ; d.rearrange(3) ; d.setInfoOnComponents(infos1)
            f=MEDCouplingFieldDouble(ON_CELLS) ; f.setName(name1) ; f.setArray(d) ; f.setMesh(m)
            f.setTime(float(i+1)+0.1,i+1,-i-1)
            f1ts=MEDFileField1TS() ; f1ts.setFieldProfile(f,mm,0,pfl) ; fmts0_1.pushBackTimeStep(f1ts)
            self.assertEqual(fmts0_1.getInfo(),tuple(infos1))
            pass
        fs.pushField(fmts0_1)
        self.assertEqual(1,len(fs))
        l=fmts0_1.splitComponents()
        self.assertEqual(3,len(l))
        for elt in l: self.assertEqual(10,len(elt))
        for elt in l: self.assertTrue(isinstance(elt,MEDFileFieldMultiTS))
        for elt in l:
            elt.setName("%s_%s"%(elt.getName(),DataArray.GetVarNameFromInfo(elt.getInfo()[0])))
            pass
        fs.pushFields(l)
        self.assertEqual(4,len(fs))
        for elt in fs: self.assertEqual(10,len(elt))
        self.assertEqual(fs.getPfls(),('pfl_NORM_QUAD4',))
        self.assertEqual(fs.getPflsReallyUsed(),('pfl_NORM_QUAD4',))
        #
        fs.write(fname,0) ; del fs
        #
        fs1=MEDFileFields(fname)
        self.assertEqual(fs1.getPfls(),('pfl_NORM_QUAD4',))
        self.assertEqual(fs1.getPflsReallyUsed(),('pfl_NORM_QUAD4',))
        self.assertEqual(4,len(fs1))
        for i in xrange(10):
            for j,fieldName in enumerate(['1stField_aa','1stField_ccc','1stField_ZZZZ']):
                f1ts=fs1[fieldName][i]
                f=f1ts.getFieldOnMeshAtLevel(ON_CELLS,0,mm)
                d=DataArrayDouble(21) ; d.iota(i*10) ; d.rearrange(3) ; d=d[:,j] ; d.setInfoOnComponent(0,infos1[j])
                self.assertTrue(d.isEqual(f.getArray(),1e-13))
                pass
            f1ts=fs1["1stField"][i]
            f=f1ts.getFieldOnMeshAtLevel(ON_CELLS,0,mm)
            d=DataArrayDouble(21) ; d.iota(i*10) ; d.rearrange(3) ; d.setInfoOnComponents(infos1)
            self.assertTrue(d.isEqual(f.getArray(),1e-13))
            pass
        pass

    def testMEDFileFieldConvertTo1(self):
        fname="Pyfile68.med"
        # building a mesh containing 4 tri3 + 5 quad4
        tri=MEDCouplingUMesh("tri",2)
        tri.allocateCells() ; tri.insertNextCell(NORM_TRI3,[0,1,2])
        tri.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,0.)]))
        tris=[tri.deepCpy() for i in xrange(4)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads=[quad.deepCpy() for i in xrange(5)]
        for i,elt in enumerate(quads): elt.translate([5+i,0])
        quads=MEDCouplingUMesh.MergeUMeshes(quads)
        m=MEDCouplingUMesh.MergeUMeshes(tris,quads)
        m.setName("mesh") ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        mm=MEDFileUMesh() ; mm.setMeshAtLevel(0,m)
        #
        ff0=MEDFileField1TS()
        f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m) ; arr=DataArrayDouble(m.getNumberOfCells()*2) ; arr.iota() ; arr.rearrange(2) ; arr.setInfoOnComponents(["X [km]","YY [mm]"]) ; f0.setArray(arr) ; f0.setName("FieldCell")
        f0.checkCoherency()
        ff0.setFieldNoProfileSBT(f0)
        #
        fspExp=[(3,[(0,(0,4),'','')]),(4,[(0,(4,9),'','')])]
        self.assertEqual(ff0.getFieldSplitedByType(),fspExp)
        #
        ff0i=ff0.convertToInt()
        self.assertEqual(ff0i.getFieldSplitedByType(),fspExp)
        self.assertTrue(arr.convertToIntArr().isEqual(ff0i.getUndergroundDataArray()))
        #
        ff1=ff0i.convertToDouble()
        self.assertTrue(ff1.getUndergroundDataArray().isEqual(ff0.getUndergroundDataArray(),1e-13))
        self.assertEqual(ff1.getFieldSplitedByType(),fspExp)
        # With profiles
        del arr,f0,ff0,ff1,ff0i,fspExp
        ff0=MEDFileField1TS()
        f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m[:7]) ; arr=DataArrayDouble(7*2) ; arr.iota() ; arr.rearrange(2) ; arr.setInfoOnComponents(["XX [pm]","YYY [hm]"]) ; f0.setArray(arr) ; f0.setName("FieldCellPfl")
        f0.checkCoherency()
        pfl=DataArrayInt.Range(0,7,1) ; pfl.setName("pfl")
        ff0.setFieldProfile(f0,mm,0,pfl)
        fspExp=[(3,[(0,(0,4),'','')]),(4,[(0,(4,7),'pfl_NORM_QUAD4','')])]
        self.assertEqual(ff0.getFieldSplitedByType(),fspExp)
        #
        ff0i=ff0.convertToInt()
        self.assertTrue(isinstance(ff0i,MEDFileIntField1TS))
        self.assertEqual(ff0i.getFieldSplitedByType(),fspExp)
        self.assertTrue(arr.convertToIntArr().isEqual(ff0i.getUndergroundDataArray()))
        #
        ff1=ff0i.convertToDouble()
        self.assertTrue(isinstance(ff1,MEDFileField1TS))
        self.assertTrue(ff1.getUndergroundDataArray().isEqual(ff0.getUndergroundDataArray(),1e-13))
        self.assertEqual(ff1.getFieldSplitedByType(),fspExp)
        ## MultiTimeSteps
        ff0=MEDFileFieldMultiTS()
        f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m[:7]) ; arr=DataArrayDouble(7*2) ; arr.iota() ; arr.rearrange(2) ; arr.setInfoOnComponents(["X [km]","YY [mm]"]) ; f0.setArray(arr) ; f0.setName("FieldCellMTime") ; f0.setTime(0.1,0,10)
        f0.checkCoherency()
        ff0.appendFieldProfile(f0,mm,0,pfl)
        f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m[:7]) ; arr=DataArrayDouble(7*2) ; arr.iota(100) ; arr.rearrange(2) ; arr.setInfoOnComponents(["X [km]","YY [mm]"]) ; f0.setArray(arr) ; f0.setName("FieldCellMTime") ; f0.setTime(1.1,1,11)
        f0.checkCoherency()
        ff0.appendFieldProfile(f0,mm,0,pfl)
        f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m[:7]) ; arr=DataArrayDouble(7*2) ; arr.iota(200) ; arr.rearrange(2) ; arr.setInfoOnComponents(["X [km]","YY [mm]"]) ; f0.setArray(arr) ; f0.setName("FieldCellMTime") ; f0.setTime(2.1,2,12)
        f0.checkCoherency()
        ff0.appendFieldProfile(f0,mm,0,pfl)
        ff1=ff0.convertToInt()
        self.assertTrue(isinstance(ff1,MEDFileIntFieldMultiTS))
        self.assertEqual(ff1.getTimeSteps(),[(0,10,0.1),(1,11,1.1),(2,12,2.1)])
        for delt,(dt,it,t) in  zip([0,100,200],ff1.getTimeSteps()):
            self.assertEqual(ff1.getFieldSplitedByType(dt,it),fspExp)
            arr=ff1.getUndergroundDataArray(dt,it)
            arr.isEqualWithoutConsideringStr(DataArrayInt.Range(delt,delt+7,1))
            pass
        self.assertEqual(ff1.getPfls(),('pfl_NORM_QUAD4', 'pfl_NORM_QUAD4', 'pfl_NORM_QUAD4'))
        #
        mm.write(fname,2)
        ff1.write(fname,0)
        #
        ff1=ff1.convertToDouble()
        self.assertTrue(isinstance(ff1,MEDFileFieldMultiTS))
        self.assertEqual(ff1.getTimeSteps(),[(0,10,0.1),(1,11,1.1),(2,12,2.1)])
        for delt,(dt,it,t) in  zip([0,100,200],ff1.getTimeSteps()):
            self.assertEqual(ff1.getFieldSplitedByType(dt,it),fspExp)
            arr=ff1.getUndergroundDataArray(dt,it)
            arr.isEqualWithoutConsideringStr(DataArrayInt.Range(delt,delt+7,1).convertToDblArr(),1e-14)
            pass
        self.assertEqual(ff1.getPfls(),('pfl_NORM_QUAD4', 'pfl_NORM_QUAD4', 'pfl_NORM_QUAD4'))
        #
        ff1=MEDFileAnyTypeFieldMultiTS.New(fname,"FieldCellMTime")
        self.assertTrue(isinstance(ff1,MEDFileIntFieldMultiTS))
        self.assertEqual(ff1.getTimeSteps(),[(0,10,0.1),(1,11,1.1),(2,12,2.1)])
        for delt,(dt,it,t) in  zip([0,100,200],ff1.getTimeSteps()):
            self.assertTrue(ff1.getFieldSplitedByType(dt,it),fspExp)
            arr=ff1.getUndergroundDataArray(dt,it)
            arr.isEqualWithoutConsideringStr(DataArrayInt.Range(delt,delt+7,1))
            pass
        self.assertEqual(ff1.getPfls(),('pfl_NORM_QUAD4',))
        pass

    def testMEDFileFieldPartialLoading(self):
        fname="Pyfile69.med"
        #
        a=DataArrayInt() ; aa=a.getHeapMemorySize()
        a.alloc(0,1)
        strMulFac=a.getHeapMemorySize()-aa ; del a ; del aa
        # building a mesh containing 30 tri3 + 40 quad4
        tri=MEDCouplingUMesh("tri",2)
        tri.allocateCells() ; tri.insertNextCell(NORM_TRI3,[0,1,2])
        tri.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,0.)]))
        tris=[tri.deepCpy() for i in xrange(30)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads=[quad.deepCpy() for i in xrange(40)]
        for i,elt in enumerate(quads): elt.translate([40+i,0])
        quads=MEDCouplingUMesh.MergeUMeshes(quads)
        m=MEDCouplingUMesh.MergeUMeshes(tris,quads)
        m.setName("mesh") ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        mm=MEDFileUMesh() ; mm.setMeshAtLevel(0,m) ; mm.write(fname,2)
        #
        ff0=MEDFileField1TS()
        f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m) ; arr=DataArrayDouble(m.getNumberOfCells()*2) ; arr.iota() ; arr.rearrange(2) ; arr.setInfoOnComponents(["X [km]","YY [mm]"]) ; f0.setArray(arr) ; f0.setName("FieldCell")
        f0.checkCoherency()
        ff0.setFieldNoProfileSBT(f0)
        ff0.write(fname,0)
        #
        fspExp=[(3,[(0,(0,30),'','')]),(4,[(0,(30,70),'','')])]
        self.assertEqual(ff0.getFieldSplitedByType(),fspExp)
        # With profiles
        ff0=MEDFileField1TS()
        f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m[:50]) ; arr=DataArrayDouble(50*2) ; arr.iota() ; arr.rearrange(2) ; arr.setInfoOnComponents(["XX [pm]","YYY [hm]"]) ; f0.setArray(arr) ; f0.setName("FieldCellPfl")
        f0.checkCoherency()
        pfl=DataArrayInt.Range(0,50,1) ; pfl.setName("pfl")
        ff0.setFieldProfile(f0,mm,0,pfl)
        fspExp=[(3,[(0,(0,30),'','')]),(4,[(0,(30,50),'pfl_NORM_QUAD4','')])]
        self.assertEqual(ff0.getFieldSplitedByType(),fspExp)
        ff0.write(fname,0)
        #
        ff0=MEDFileField1TS(fname,False)
        self.assertEqual(ff0.getName(),"FieldCell")
        self.assertTrue(not ff0.getUndergroundDataArray().isAllocated())
        self.assertEqual(ff0.getUndergroundDataArray().getInfoOnComponents(),['X [km]','YY [mm]'])
        heap_memory_ref=ff0.getHeapMemorySize()
        self.assertIn(heap_memory_ref,xrange(182,298+2*strMulFac))
        ff0.loadArrays() ##
        arr=DataArrayDouble(140) ; arr.iota() ; arr.rearrange(2)
        self.assertTrue(ff0.getUndergroundDataArray().isEqualWithoutConsideringStr(arr,1e-14))
        self.assertEqual(ff0.getHeapMemorySize()-heap_memory_ref,70*8*2)
        #
        ff0=MEDFileField1TS(fname,"FieldCellPfl",False)
        self.assertEqual(ff0.getUndergroundDataArray().getInfoOnComponents(),["XX [pm]","YYY [hm]"])
        heap_memory_ref=ff0.getHeapMemorySize()
        self.assertIn(heap_memory_ref,xrange(350,415+6*strMulFac))
        ff0.loadArrays() ##
        arr=DataArrayDouble(100) ; arr.iota() ; arr.rearrange(2)
        self.assertTrue(ff0.getUndergroundDataArray().isEqualWithoutConsideringStr(arr,1e-14))
        self.assertEqual(ff0.getHeapMemorySize()-heap_memory_ref,50*8*2)
        ff0.loadArrays() ##
        self.assertTrue(ff0.getUndergroundDataArray().isEqualWithoutConsideringStr(arr,1e-14))
        self.assertEqual(ff0.getHeapMemorySize()-heap_memory_ref,50*8*2)
        ff0.getUndergroundDataArray().setIJ(30,1,5.5)
        self.assertTrue(not ff0.getUndergroundDataArray().isEqualWithoutConsideringStr(arr,1e-14))
        ff0.loadArrays() ##
        self.assertTrue(ff0.getUndergroundDataArray().isEqualWithoutConsideringStr(arr,1e-14))
        ff0.getUndergroundDataArray().setIJ(30,1,5.5)
        self.assertTrue(not ff0.getUndergroundDataArray().isEqualWithoutConsideringStr(arr,1e-14))
        ff0.loadArraysIfNecessary() ##
        self.assertEqual(ff0.getUndergroundDataArray().getIJ(30,1),5.5)
        self.assertTrue(not ff0.getUndergroundDataArray().isEqualWithoutConsideringStr(arr,1e-14))
        heap_memory_ref=ff0.getHeapMemorySize()
        self.assertIn(heap_memory_ref,xrange(1100,1215+2*strMulFac))
        ff0.unloadArrays()
        hmd=ff0.getHeapMemorySize()-heap_memory_ref
        self.assertEqual(hmd,-800) # -50*8*2
        ff0.loadArrays() ##
        self.assertEqual(ff0.getHeapMemorySize()-heap_memory_ref,0)
        #
        ff0=MEDFileField1TS(fname,"FieldCellPfl",-1,-1,False)
        heap_memory_ref=ff0.getHeapMemorySize()
        self.assertIn(heap_memory_ref,xrange(299,415+6*strMulFac))
        ff0.loadArrays() ##
        self.assertTrue(ff0.getUndergroundDataArray().isEqualWithoutConsideringStr(arr,1e-14))
        self.assertEqual(ff0.getHeapMemorySize()-heap_memory_ref,50*8*2)
        #
        fieldName="FieldCellMultiTS"
        ff0=MEDFileFieldMultiTS()
        for t in xrange(20):
            f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m) ; arr=DataArrayDouble(m.getNumberOfCells()*2) ; arr.iota(float(t+1000)) ; arr.rearrange(2) ; arr.setInfoOnComponents(["X [km]","YY [mm]"]) ; f0.setArray(arr) ; f0.setName(fieldName)
            f0.setTime(float(t)+0.1,t,100+t)
            f0.checkCoherency()
            ff0.appendFieldNoProfileSBT(f0)
            pass
        ff0.write(fname,0)
        #
        ff0=MEDFileAnyTypeFieldMultiTS.New(fname,fieldName,False)
        heap_memory_ref=ff0.getHeapMemorySize()
        self.assertIn(heap_memory_ref,xrange(5536,5956+(80+26)*strMulFac))
        ff0.loadArrays()
        self.assertEqual(ff0.getHeapMemorySize()-heap_memory_ref,20*70*8*2)
        del ff0
        #
        ffs=MEDFileFields(fname,False)
        heap_memory_ref=ffs.getHeapMemorySize()
        self.assertIn(heap_memory_ref,xrange(5335,6687+(80+50)*strMulFac))
        ffs.loadArrays()
        self.assertEqual(ffs.getHeapMemorySize()-heap_memory_ref,20*70*8*2+70*8*2+50*8*2)
        pass

    def testMEDFileMeshReadSelector1(self):
        mrs=MEDFileMeshReadSelector()
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs.__str__() ; mrs.__repr__()
        #
        mrs=MEDFileMeshReadSelector(0)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(1)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(2)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(3)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(4)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(5)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(6)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(7)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(8)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(9)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(10)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(11)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(12)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(13)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(14)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(15)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(16)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(17)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(18)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(19)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(20)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(21)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(22)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(23)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(24)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(25)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(26)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(27)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(28)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(29)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(30)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(31)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and not mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(32)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(33)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(34)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(35)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(36)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(37)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(38)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(39)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(40)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(41)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(42)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(43)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(44)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(45)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(46)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(47)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and not mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(48)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(49)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(50)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(51)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(52)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(53)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(54)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(55)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and not mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(56)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(57)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(58)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(59)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and not mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(60)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(61)
        self.assertTrue(mrs.isCellFamilyFieldReading() and not mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(62)
        self.assertTrue(not mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        mrs=MEDFileMeshReadSelector(63)
        self.assertTrue(mrs.isCellFamilyFieldReading() and mrs.isNodeFamilyFieldReading() and mrs.isCellNameFieldReading() and mrs.isNodeNameFieldReading() and mrs.isCellNumFieldReading() and mrs.isNodeNumFieldReading())
        #
        mrs=MEDFileMeshReadSelector(63)
        mrs.setCellFamilyFieldReading(False)
        self.assertEqual(mrs.getCode(),62)
        mrs.setCellFamilyFieldReading(True)
        self.assertEqual(mrs.getCode(),63)
        mrs.setNodeFamilyFieldReading(False)
        self.assertEqual(mrs.getCode(),61)
        mrs.setNodeFamilyFieldReading(True)
        self.assertEqual(mrs.getCode(),63)
        mrs.setCellNameFieldReading(False)
        self.assertEqual(mrs.getCode(),59)
        mrs.setCellNameFieldReading(True)
        self.assertEqual(mrs.getCode(),63)
        mrs.setNodeNameFieldReading(False)
        self.assertEqual(mrs.getCode(),55)
        mrs.setNodeNameFieldReading(True)
        self.assertEqual(mrs.getCode(),63)
        mrs.setCellNumFieldReading(False)
        self.assertEqual(mrs.getCode(),47)
        mrs.setCellNumFieldReading(True)
        self.assertEqual(mrs.getCode(),63)
        mrs.setNodeNumFieldReading(False)
        self.assertEqual(mrs.getCode(),31)
        mrs.setNodeNumFieldReading(True)
        self.assertEqual(mrs.getCode(),63)
        pass
    
    def testPartialReadOfMeshes(self):
        fname="Pyfile70.med"
        # building a mesh containing 4 tri3 + 5 quad4
        tri=MEDCouplingUMesh("tri",2)
        tri.allocateCells() ; tri.insertNextCell(NORM_TRI3,[0,1,2])
        tri.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,0.)]))
        tris=[tri.deepCpy() for i in xrange(4)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads=[quad.deepCpy() for i in xrange(5)]
        for i,elt in enumerate(quads): elt.translate([5+i,0])
        quads=MEDCouplingUMesh.MergeUMeshes(quads)
        m=MEDCouplingUMesh.MergeUMeshes(tris,quads)
        m.setName("mesh") ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        m1=m.buildDescendingConnectivity()[0]
        mm=MEDFileUMesh() ; mm.setMeshes([m,m1])
        #
        grp0=DataArrayInt([1,2,3,5,6]) ; grp0.setName("grp0")
        grp1=DataArrayInt([1,2,3,5,7,8]) ; grp1.setName("grp1")
        mm.setGroupsAtLevel(0,[grp0,grp1])
        grp2=DataArrayInt.Range(0,32,2) ; grp2.setName("grp2")
        grp3=DataArrayInt.Range(1,32,7) ; grp3.setName("grp3")
        mm.setGroupsAtLevel(-1,[grp2,grp3])
        grp4=DataArrayInt.Range(0,32,2) ; grp4.setName("grp4")
        grp5=DataArrayInt.Range(1,32,7) ; grp5.setName("grp5")
        mm.setGroupsAtLevel(1,[grp4,grp5])
        mm.setRenumFieldArr(0,DataArrayInt.Range(2,11,1))
        mm.setRenumFieldArr(-1,DataArrayInt.Range(3,35,1))
        mm.setRenumFieldArr(1,DataArrayInt.Range(4,36,1))
        #
        mm.write(fname,2)
        ##
        mm=MEDFileMesh.New(fname,"mesh",-1,-1,MEDFileMeshReadSelector())
        b4_ref_heap_mem=mm.getHeapMemorySize()
        mm.getMeshAtLevel(0)## please let this line : force to move 1GTUMesh -> UMesh
        mm.getMeshAtLevel(-1)## please let this line : force to move 1GTUMesh -> UMesh
        ref_heap_mem=mm.getHeapMemorySize()
        # check the gain of memory using 1GTUMesh instead of UMesh
        self.assertTrue(ref_heap_mem-b4_ref_heap_mem>=(32+9)*4*2-32)# 32+9=nbCells 4=sizeof(int) 2=the types+index -32=loss linked to vector
        #
        mm=MEDFileMesh.New(fname,MEDFileMeshReadSelector(0))
        self.assertEqual(len(mm.getGroupsNames()),0)
        self.assertTrue(mm.getMeshAtLevel(0).isEqual(m,1e-13))
        self.assertTrue(mm.getMeshAtLevel(-1).isEqual(m1,1e-13))
        self.assertTrue(mm.getFamilyFieldAtLevel(0) is None)
        self.assertTrue(mm.getFamilyFieldAtLevel(-1) is None)
        self.assertTrue(mm.getFamilyFieldAtLevel(1) is None)
        self.assertTrue(mm.getNumberFieldAtLevel(0) is None)
        self.assertTrue(mm.getNumberFieldAtLevel(-1) is None)
        self.assertTrue(mm.getNumberFieldAtLevel(1) is None)
        delta1=ref_heap_mem-mm.getHeapMemorySize()
        self.assertTrue(delta1>=4*(32+9)*3+32*4*3)
        #
        mm=MEDFileMesh.New(fname,MEDFileMeshReadSelector(1))
        self.assertEqual(len(mm.getGroupsNames()),6)
        self.assertTrue(mm.getMeshAtLevel(0).isEqual(m,1e-13))
        self.assertTrue(mm.getMeshAtLevel(-1).isEqual(m1,1e-13))
        self.assertTrue(mm.getFamilyFieldAtLevel(0)!=None)
        self.assertTrue(mm.getFamilyFieldAtLevel(-1)!=None)
        self.assertTrue(mm.getFamilyFieldAtLevel(1) is None)
        self.assertTrue(mm.getNumberFieldAtLevel(0) is None)
        self.assertTrue(mm.getNumberFieldAtLevel(-1) is None)
        self.assertTrue(mm.getNumberFieldAtLevel(1) is None)
        delta2=ref_heap_mem-mm.getHeapMemorySize()
        self.assertTrue(delta2<delta1)
        self.assertTrue(delta2>=4*(32+9)*1+32*4*3)
        #
        mm=MEDFileUMesh(fname,MEDFileMeshReadSelector(3))
        self.assertEqual(len(mm.getGroupsNames()),6)
        self.assertTrue(mm.getMeshAtLevel(0).isEqual(m,1e-13))
        self.assertTrue(mm.getMeshAtLevel(-1).isEqual(m1,1e-13))
        self.assertTrue(mm.getFamilyFieldAtLevel(0)!=None)
        self.assertTrue(mm.getFamilyFieldAtLevel(-1)!=None)
        self.assertTrue(mm.getFamilyFieldAtLevel(1)!=None)
        self.assertTrue(mm.getNumberFieldAtLevel(0) is None)
        self.assertTrue(mm.getNumberFieldAtLevel(-1) is None)
        self.assertTrue(mm.getNumberFieldAtLevel(1) is None)
        delta3=ref_heap_mem-mm.getHeapMemorySize()
        self.assertTrue(delta3<delta2)
        self.assertTrue(delta3>=4*(32+9)*1+32*4*1)
        #
        mm=MEDFileUMesh(fname,"mesh",-1,-1,MEDFileMeshReadSelector(19))
        self.assertEqual(len(mm.getGroupsNames()),6)
        self.assertTrue(mm.getMeshAtLevel(0).isEqual(m,1e-13))
        self.assertTrue(mm.getMeshAtLevel(-1).isEqual(m1,1e-13))
        self.assertTrue(mm.getFamilyFieldAtLevel(0)!=None)
        self.assertTrue(mm.getFamilyFieldAtLevel(-1)!=None)
        self.assertTrue(mm.getFamilyFieldAtLevel(1)!=None)
        self.assertTrue(mm.getNumberFieldAtLevel(0)!=None)
        self.assertTrue(mm.getNumberFieldAtLevel(-1)!=None)
        self.assertTrue(mm.getNumberFieldAtLevel(1) is None)
        delta4=ref_heap_mem-mm.getHeapMemorySize()
        self.assertTrue(delta4<delta3)
        self.assertTrue(delta4>=32*4*2)
        #
        mm=MEDFileUMesh.New(fname,"mesh",-1,-1,MEDFileMeshReadSelector(51))
        self.assertEqual(len(mm.getGroupsNames()),6)
        self.assertTrue(mm.getMeshAtLevel(0).isEqual(m,1e-13))
        self.assertTrue(mm.getMeshAtLevel(-1).isEqual(m1,1e-13))
        self.assertTrue(mm.getFamilyFieldAtLevel(0)!=None)
        self.assertTrue(mm.getFamilyFieldAtLevel(-1)!=None)
        self.assertTrue(mm.getFamilyFieldAtLevel(1)!=None)
        self.assertTrue(mm.getNumberFieldAtLevel(0)!=None)
        self.assertTrue(mm.getNumberFieldAtLevel(-1)!=None)
        self.assertTrue(mm.getNumberFieldAtLevel(1)!=None)
        delta5=ref_heap_mem-mm.getHeapMemorySize()
        self.assertTrue(delta5<delta4)
        self.assertEqual(delta5,0)
        pass

    # this test checks that setFieldProfile perform a check of the array length
    # compared to the profile length. This test also checks that mesh attribute of field
    # is not used by setFieldProfile (because across this test mesh is equal to None)
    def testCheckCompatibilityPfl1(self):
        # building a mesh containing 4 tri3 + 5 quad4
        tri=MEDCouplingUMesh("tri",2)
        tri.allocateCells() ; tri.insertNextCell(NORM_TRI3,[0,1,2])
        tri.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,0.)]))
        tris=[tri.deepCpy() for i in xrange(4)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads=[quad.deepCpy() for i in xrange(5)]
        for i,elt in enumerate(quads): elt.translate([5+i,0])
        quads=MEDCouplingUMesh.MergeUMeshes(quads)
        m=MEDCouplingUMesh.MergeUMeshes(tris,quads)
        m.setName("mesh") ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        m1=m.buildDescendingConnectivity()[0]
        mm=MEDFileUMesh() ; mm.setMeshes([m,m1])
        #
        f1ts=MEDFileField1TS()
        f=MEDCouplingFieldDouble(ON_NODES)
        vals=DataArrayDouble(7) ; vals.iota(1000)
        f.setArray(vals)
        f.setName("anonymous") # f has no mesh it is not a bug
        pfl=DataArrayInt([0,1,2,3,4,5,6]) ; pfl.setName("pfl")
        f1ts.setFieldProfile(f,mm,0,pfl)
        #
        f1ts=MEDFileField1TS()
        f=MEDCouplingFieldDouble(ON_NODES)
        vals=DataArrayDouble(8) ; vals.iota(1000)
        f.setArray(vals)
        f.setName("anonymous") # f has no mesh it is not a bug
        pfl=DataArrayInt([0,1,2,3,4,5,6]) ; pfl.setName("pfl")
        self.assertRaises(InterpKernelException,f1ts.setFieldProfile,f,mm,0,pfl)
        #
        f1ts=MEDFileField1TS()
        f=MEDCouplingFieldDouble(ON_CELLS)
        vals=DataArrayDouble(7) ; vals.iota(1000)
        f.setArray(vals)
        f.setName("anonymous") # f has no mesh it is not a bug
        pfl=DataArrayInt([1,2,3,5,6,7,8]) ; pfl.setName("pfl")
        f1ts.setFieldProfile(f,mm,0,pfl)
        self.assertTrue(f1ts.getUndergroundDataArray().isEqual(vals,1e-10))
        #
        f1ts=MEDFileField1TS()
        f=MEDCouplingFieldDouble(ON_GAUSS_PT)
        vals=DataArrayDouble(27) ; vals.iota(1000)
        f.setArray(vals)
        f.setName("anonymous") # f has no mesh it is not a bug
        pfl=DataArrayInt([1,2,3,5,6,7,8]) ; pfl.setName("pfl")
        f.setMesh(m[pfl])
        f.setGaussLocalizationOnCells([0,1],[0.,0.,1.,0.,1.,1.],[0.3,0.3,0.7,0.7,0.1,0.1],[0.3,0.6,0.1])
        f.setGaussLocalizationOnCells([2],[0.,0.,1.,0.,1.,1.],[0.3,0.3],[1.])
        f.setGaussLocalizationOnCells([3,4,5,6],[0.,0.,1.,0.,1.,1.,0.,1.],[0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.5,0.5],[0.2,0.3,0.4,0.07,0.03])
        f.setMesh(None)
        f1ts.setFieldProfile(f,mm,0,pfl)
        self.assertTrue(f1ts.getUndergroundDataArray().isEqual(vals,1e-10))
        vals=DataArrayDouble(26) ; vals.iota(1040) ; f.setArray(vals)
        self.assertRaises(InterpKernelException,f1ts.setFieldProfile,f,mm,0,pfl)
        vals=DataArrayDouble(27) ; vals.iota(1000)
        self.assertTrue(f1ts.getUndergroundDataArray().isEqual(vals,1e-10))
        #
        f1ts=MEDFileField1TS()
        f=MEDCouplingFieldDouble(ON_GAUSS_NE)
        vals=DataArrayDouble(25) ; vals.iota(1000)
        f.setArray(vals)
        f.setName("anonymous") # f has no mesh it is not a bug
        pfl=DataArrayInt([1,2,3,5,6,7,8]) ; pfl.setName("pfl")
        f1ts.setFieldProfile(f,mm,0,pfl)
        self.assertTrue(f1ts.getUndergroundDataArray().isEqual(vals,1e-10))
        vals2=DataArrayDouble(26) ; vals2.iota(1050)
        f.setArray(vals2)
        self.assertRaises(InterpKernelException,f1ts.setFieldProfile,f,mm,0,pfl)
        self.assertTrue(f1ts.getUndergroundDataArray().isEqual(vals,1e-10))
        #
        f1ts=MEDFileField1TS()
        self.assertRaises(InterpKernelException,f1ts.setFieldProfile,f,mm,0,pfl)
        self.assertRaises(InterpKernelException,f1ts.setFieldProfile,f,mm,0,pfl)
        f.setArray(vals)
        f1ts.setFieldProfile(f,mm,0,pfl)
        pass
    
    def testWRMeshWithNoCells(self):
        fname="Pyfile71.med"
        a=DataArrayDouble(4) ; a.iota()
        c=MEDCouplingCMesh() ; c.setCoords(a,a) ; m0=c.buildUnstructured()
        m00=MEDCouplingUMesh("mesh",1) ; m00.setCoords(m0.getCoords()) ; m00.allocateCells(0)
        m=MEDFileUMesh()
        m.setMeshAtLevel(0,m00)
        m.setRenumFieldArr(1,DataArrayInt(range(10,26)))
        m.setFamilyFieldArr(1,DataArrayInt([-1,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,0,-1,-3,-3,-3]))
        m.write(fname,2)
        del m,a,c,m0,m00
        #
        m=MEDFileMesh.New(fname)
        self.assertEqual((),m.getNonEmptyLevels())
        self.assertTrue(m.getCoords().isEqual(DataArrayDouble([(0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(2,1),(3,1),(0,2),(1,2),(2,2),(3,2),(0,3),(1,3),(2,3),(3,3)]),1e-12))
        self.assertTrue(m.getNumberFieldAtLevel(1).isEqual(DataArrayInt(range(10,26))))
        self.assertTrue(m.getFamilyFieldAtLevel(1).isEqual(DataArrayInt([-1,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,0,-1,-3,-3,-3])))
        pass

    #@unittest.skipUnless(False,"requires Vadim's green light")
    def testWRQPolyg1(self):
        fname="Pyfile72.med"
        m=MEDCoupling1SGTUMesh("mesh",NORM_QUAD4) ; m.allocateCells()
        m.insertNextCell([0,2,1,3])
        m.setCoords(DataArrayDouble([0.,0.,1.,1.,1.,0.,0.,1.],4,2))
        #
        ms=[m.deepCpy() for i in xrange(4)]
        for i,elt in enumerate(ms):
            elt.translate([float(i)*1.5,0.])
            pass
        m0=MEDCoupling1SGTUMesh.Merge1SGTUMeshes(ms).buildUnstructured()
        m0.convertAllToPoly()
        #
        ms=[m.deepCpy() for i in xrange(5)]
        for i,elt in enumerate(ms):
            elt.translate([float(i)*1.5,1.5])
            pass
        m1=MEDCoupling1SGTUMesh.Merge1SGTUMeshes(ms).buildUnstructured()
        m1.convertAllToPoly()
        m1.convertLinearCellsToQuadratic()
        #
        m=MEDCouplingUMesh.MergeUMeshes(m0,m1)
        ##
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        grp0=DataArrayInt([0,2,3]) ; grp0.setName("grp0")
        grp1=DataArrayInt([4,6,7]) ; grp1.setName("grp1")
        grp2=DataArrayInt([0,1,2,4,5,6]) ; grp2.setName("grp2")
        mm.setGroupsAtLevel(0,[grp0,grp1,grp2])
        ##
        mm.write(fname,2)
        del mm
        #
        mm_read=MEDFileUMesh(fname)
        self.assertTrue(mm_read.getGroupArr(0,"grp0").isEqual(grp0))
        self.assertTrue(mm_read.getGroupArr(0,"grp1").isEqual(grp1))
        self.assertTrue(mm_read.getGroupArr(0,"grp2").isEqual(grp2))
        self.assertTrue(mm_read.getMeshAtLevel(0).isEqual(m,1e-12))
        ##
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setName("MyFirstField")
        f.setMesh(m)
        arr0=DataArrayDouble(9) ; arr0.iota()
        arr1=DataArrayDouble(9) ; arr1.iota(100)
        arr=DataArrayDouble.Meld(arr0,arr1) ; arr.setInfoOnComponents(["mm [kg]","sds [m]"])
        f.setArray(arr) ; f.checkCoherency()
        f.setTime(5.6,1,2)
        ff=MEDFileField1TS()
        ff.setFieldNoProfileSBT(f)
        ff.write(fname,0)
        ##
        ff_read=MEDFileField1TS(fname)
        f_read=ff_read.getFieldOnMeshAtLevel(ON_CELLS,0,mm_read)
        self.assertTrue(f_read.isEqual(f,1e-12,1e-12))
        pass

    def testLoadIfNecessaryOnFromScratchFields0(self):
        """
        This test checks that a call to loadArraysIfNecessary works (does nothing) on field data structure whatever its level 1TS, MTS, Fields.
        """
        fname="Pyfile77.med"
        coords=DataArrayDouble([(0,0,0),(2,1,0),(1,0,0),(1,1,0),(2,0,0),(0,1,0)])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coords)
        m.allocateCells()
        m.insertNextCell(NORM_QUAD4,[0,5,3,2])
        m.insertNextCell(NORM_QUAD4,[4,2,3,1])
        m.finishInsertingCells()
        #
        mm=MEDFileUMesh() ; mm.setMeshAtLevel(0,m)
        ms=MEDFileMeshes() ; ms.pushMesh(mm)
        fs=MEDFileFields()
        arrs=4*[None]
        #
        ff0=MEDFileFieldMultiTS() ; fs.pushField(ff0)
        f0=MEDCouplingFieldDouble(ON_GAUSS_NE) ; f0.setMesh(m) ; f0.setTimeUnit("ms")
        f0.setTime(1.1,1,1)
        f0.setName("myELNOField")
        arrs[0]=DataArrayDouble([7,5,3,1,5,3,1,7]) ; arrs[0].setInfoOnComponent(0,"Comp0")
        f0.setArray(arrs[0])
        ff0.appendFieldNoProfileSBT(f0)
        #
        f0.setTime(2.2,2,1)
        arrs[1]=DataArrayDouble([1,7,5,3,7,5,3,1]) ; arrs[1].setInfoOnComponent(0,"Comp0")
        f0.setArray(arrs[1])
        ff0.appendFieldNoProfileSBT(f0)
        #
        f0.setTime(3.3,3,1)
        arrs[2]=DataArrayDouble([3,1,7,5,1,7,5,3]) ; arrs[2].setInfoOnComponent(0,"Comp0")
        f0.setArray(arrs[2])
        ff0.appendFieldNoProfileSBT(f0)
        #
        f0.setTime(4.4,4,1)
        arrs[3]=DataArrayDouble([5,3,1,7,3,1,7,5]) ; arrs[3].setInfoOnComponent(0,"Comp0")
        f0.setArray(arrs[3])
        ff0.appendFieldNoProfileSBT(f0)
        #
        for i,arr in enumerate(arrs):
            self.assertTrue(fs[0][i].getUndergroundDataArray().isEqual(arr,1e-12))
            fs[0][i].loadArraysIfNecessary()
            self.assertTrue(fs[0][i].getUndergroundDataArray().isEqual(arr,1e-12))
            pass
        fs.loadArraysIfNecessary()
        for i,arr in enumerate(arrs):
            self.assertTrue(fs[0][i].getUndergroundDataArray().isEqual(arr,1e-12))
            pass
        fs[0].loadArraysIfNecessary()
        for i,arr in enumerate(arrs):
            self.assertTrue(fs[0][i].getUndergroundDataArray().isEqual(arr,1e-12))
            pass
        pass
    
    def testField1TSSetFieldNoProfileSBTPerGeoTypes(self):
        """ This test is very important, because the same mechanism is used by the MEDReader to generate a field on all the mesh without any processing and memory.
        """
        fname="Pyfile78.med"
        coords=DataArrayDouble([-0.3,-0.3,0., 0.2,-0.3,0., 0.7,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.7,0.2,0., -0.3,0.7,0., 0.2,0.7,0., 0.7,0.7,0. ],9,3)
        targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4];
        m0=MEDCouplingUMesh("mesh",3) ; m0.setCoords(coords)
        m0.allocateCells()
        for elt in [[0,1,2,3],[1,2,3,4],[2,3,4,5],[3,4,5,6],[4,5,6,7],[5,6,7,8]]:#6
            m0.insertNextCell(NORM_TETRA4,elt)
            pass
        for elt in [[0,1,2,3,4],[1,2,3,4,5],[2,3,4,5,6],[3,4,5,6,7],[4,5,6,7,8]]:#5
            m0.insertNextCell(NORM_PYRA5,elt)
            pass
        for elt in [[0,1,2,3,4,5],[1,2,3,4,5,6],[2,3,4,5,6,7],[3,4,5,6,7,8]]:#4
            m0.insertNextCell(NORM_PENTA6,elt)
            pass
        m0.checkCoherency2()
        m1=MEDCouplingUMesh(); m1.setName("mesh")
        m1.setMeshDimension(2);
        m1.allocateCells(5);
        m1.insertNextCell(NORM_TRI3,3,targetConn[4:7]);
        m1.insertNextCell(NORM_TRI3,3,targetConn[7:10]);
        m1.insertNextCell(NORM_QUAD4,4,targetConn[0:4]);
        m1.insertNextCell(NORM_QUAD4,4,targetConn[10:14]);
        m1.insertNextCell(NORM_QUAD4,4,targetConn[14:18]);
        m1.setCoords(coords);
        m3=MEDCouplingUMesh("mesh",0) ; m3.setCoords(coords)
        m3.allocateCells()
        m3.insertNextCell(NORM_POINT1,[2])
        m3.insertNextCell(NORM_POINT1,[3])
        m3.insertNextCell(NORM_POINT1,[4])
        m3.insertNextCell(NORM_POINT1,[5])
        #
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m0)
        mm.setMeshAtLevel(-1,m1)
        mm.setMeshAtLevel(-3,m3)
        mm.write(fname,2)
        #### The file is written only with one mesh and no fields. Let's put a field on it geo types per geo types.
        mm=MEDFileMesh.New(fname)
        fs=MEDFileFields()
        fmts=MEDFileFieldMultiTS()
        f1ts=MEDFileField1TS()
        for lev in mm.getNonEmptyLevels():
            for gt in mm.getGeoTypesAtLevel(lev):
                p0=mm.getDirectUndergroundSingleGeoTypeMesh(gt)
                f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(p0)
                arr=DataArrayDouble(f.getNumberOfTuplesExpected()) ; arr.iota()
                f.setArray(arr) ; f.setName("f0")
                f1ts.setFieldNoProfileSBT(f)
                pass
            pass
        self.assertEqual(mm.getNonEmptyLevels(),(0,-1,-3))
        for lev in [0,-1,-3]:
            mm.getDirectUndergroundSingleGeoTypeMeshes(lev) # please let this line, it is for the test to emulate that
            pass
        fmts.pushBackTimeStep(f1ts)
        fs.pushField(fmts)
        fs.write(fname,0)
        del fs,fmts,f1ts
        #### The file contains now one mesh and one cell field with all cells wathever their level ang type fetched.
        fs=MEDFileFields(fname)
        self.assertEqual(len(fs),1)
        self.assertEqual(len(fs[0]),1)
        f1ts=fs[0][0]
        self.assertEqual(f1ts.getFieldSplitedByType(),[(0,[(0,(0,4),'','')]),(3,[(0,(4,6),'','')]),(4,[(0,(6,9),'','')]),(14,[(0,(9,15),'','')]),(15,[(0,(15,20),'','')]),(16,[(0,(20,24),'','')])])
        self.assertTrue(f1ts.getUndergroundDataArray().isEqual(DataArrayDouble([0,1,2,3,0,1,0,1,2,0,1,2,3,4,5,0,1,2,3,4,0,1,2,3]),1e-12))
        pass

    def testMEDFileUMeshSetName(self):
        """ This test is a small but important one for MEDReader in sauv mode. When .sauv file is loaded the convertion is performed in memory and a preparation is done then.
        This preparation makes access to internal MEDCouplingMesh pointers whose name must be updated.
        """
        fname="Pyfile79.med"
        targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4];
        mm=MEDFileUMesh()
        m0=MEDCouplingUMesh() ; m0.setMeshDimension(2) # important no name here.
        coords=DataArrayDouble([-0.3,-0.3,0., 0.2,-0.3,0., 0.7,-0.3,0., -0.3,0.2,0., 0.2,0.2,0., 0.7,0.2,0., -0.3,0.7,0., 0.2,0.7,0., 0.7,0.7,0. ],9,3)
        m0.allocateCells(5);
        m0.insertNextCell(NORM_TRI3,3,targetConn[4:7]);
        m0.insertNextCell(NORM_TRI3,3,targetConn[7:10]);
        m0.insertNextCell(NORM_QUAD4,4,targetConn[0:4]);
        m0.insertNextCell(NORM_QUAD4,4,targetConn[10:14]);
        m0.insertNextCell(NORM_QUAD4,4,targetConn[14:18]);
        m0.setCoords(coords);
        mm.setMeshAtLevel(0,m0)
        m2=MEDCouplingUMesh() ; m2.setMeshDimension(0) ; m2.setCoords(coords) # important no name here.
        m2.allocateCells()
        m2.insertNextCell(NORM_POINT1,[2])
        m2.insertNextCell(NORM_POINT1,[3])
        m2.insertNextCell(NORM_POINT1,[4])
        m2.insertNextCell(NORM_POINT1,[5])
        mm.setMeshAtLevel(-2,m2)
        self.assertEqual(mm.getName(),"")
        self.assertEqual(mm.getMeshAtLevel(0).getName(),"")
        mm.forceComputationOfParts()
        self.assertEqual(mm.getDirectUndergroundSingleGeoTypeMesh(NORM_TRI3).getName(),"")
        mm.setName("abc")
        self.assertEqual(mm.getName(),"abc")
        self.assertEqual(mm.getDirectUndergroundSingleGeoTypeMesh(NORM_TRI3).getName(),"abc")
        self.assertEqual(mm.getDirectUndergroundSingleGeoTypeMesh(NORM_QUAD4).getName(),"abc")
        self.assertEqual(mm.getDirectUndergroundSingleGeoTypeMesh(NORM_POINT1).getName(),"abc")
        self.assertEqual(mm.getMeshAtLevel(0).getName(),"abc")
        pass

    def testMEDFileFieldsUnloadArraysWithoutDataLoss1(self):
        fileName="Pyfile80.med"
        m=MEDCouplingCMesh() ; m.setName("cmesh")
        arr=DataArrayDouble(6) ; arr.iota()
        m.setCoords(arr,arr)
        nbCells=m.getNumberOfCells()
        self.assertEqual(25,nbCells)
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setName("FieldOnCell") ; f.setMesh(m)
        arr=DataArrayDouble(nbCells) ; arr.iota()
        mm=MEDFileCMesh()
        mm.setMesh(m)
        #
        fmts=MEDFileFieldMultiTS()
        #
        for i in xrange(nbCells):
            t=(float(i)+0.1,i+1,-i-2)
            f.setTime(*t)
            arr2=DataArrayDouble(nbCells)
            perm=DataArrayInt(nbCells) ; perm.iota(i) ; perm%=nbCells
            arr2[perm]=arr
            f.setArray(arr2)
            f1ts=MEDFileField1TS()
            f1ts.setFieldNoProfileSBT(f)
            fmts.pushBackTimeStep(f1ts)
            pass
        fmts.unloadArraysWithoutDataLoss()
        self.assertTrue(fmts[0].getUndergroundDataArray().isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.]),1e-12))
        fs=MEDFileFields() ; fs.pushField(fmts)
        self.assertTrue(fs[0][0].getUndergroundDataArray().isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.]),1e-12))
        fs.unloadArraysWithoutDataLoss()
        self.assertTrue(fs[0][0].getUndergroundDataArray().isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.]),1e-12))
        f1ts=fs[0][0]
        self.assertTrue(f1ts.getUndergroundDataArray().isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.]),1e-12))
        f1ts.unloadArraysWithoutDataLoss()
        self.assertTrue(f1ts.getUndergroundDataArray().isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.]),1e-12))
        mm.write(fileName,2)
        fs.write(fileName,0)
        del m,fmts,mm,f,f1ts
        ##
        mm=MEDFileMesh.New(fileName)
        fmts=MEDFileFieldMultiTS(fileName)
        self.assertTrue(fmts[0].getUndergroundDataArray().isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.]),1e-12))
        fmts.unloadArraysWithoutDataLoss()
        self.assertTrue(not fmts[0].getUndergroundDataArray().isAllocated())
        fmts.loadArraysIfNecessary()
        self.assertTrue(fmts[0].getUndergroundDataArray().isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.]),1e-12))
        del mm,fmts
        fs=MEDFileFields(fileName)
        self.assertTrue(fs[0][0].getUndergroundDataArray().isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.]),1e-12))
        fs.unloadArraysWithoutDataLoss()
        self.assertTrue(not fs[0][0].getUndergroundDataArray().isAllocated())
        fs.loadArraysIfNecessary()
        self.assertTrue(fs[0][0].getUndergroundDataArray().isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.]),1e-12))
        del fs
        f1ts=MEDFileField1TS(fileName)
        self.assertTrue(f1ts.getUndergroundDataArray().isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.]),1e-12))
        f1ts.unloadArraysWithoutDataLoss()
        self.assertTrue(not f1ts.getUndergroundDataArray().isAllocated())
        f1ts.loadArraysIfNecessary()
        self.assertTrue(f1ts.getUndergroundDataArray().isEqual(DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.]),1e-12))
        pass

    def testMEDFileUMeshLoadPart1(self):
        """ This method tests MEDFileUMesh.LoadPart that loads only a part of a specified mesh in a MED file. The part is specfied using a slice of cell ids. Only nodes on which cells lies are loaded to reduce at most the amount of
        memory of the returned instance.
        """
        fileName="Pyfile81.med"
        arr=DataArrayDouble(6) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
        m=m.buildUnstructured()
        m.setName("Mesh")
        m.changeSpaceDimension(3,0.)
        infos=["aa [b]","cc [de]","gg [klm]"]
        m.getCoords().setInfoOnComponents(infos)
        m.checkCoherency2()
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        m1=MEDCouplingCMesh() ; m1.setCoords(arr) ; m1.setName("Mesh") 
        m1=m1.buildUnstructured() ; m1.setCoords(m.getCoords())
        mm.setMeshAtLevel(-1,m1)
        renum0=DataArrayInt([3,6,7,10,11,0,2,1,9,8,5,4,12,13,14,24,23,22,21,20,19,18,17,16,15])
        famField0=DataArrayInt([-3,-6,-7,-10,-11,0,-2,-1,-9,-8,-5,-4,-12,-13,-14,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15])
        namesCellL0=DataArrayAsciiChar(25,16)
        namesCellL0[:]=["Cell#%.3d        "%(i) for i in xrange(25)]
        renumM1=DataArrayInt([3,4,0,2,1])
        famFieldM1=DataArrayInt([-3,-4,0,-2,-1])
        mm.setRenumFieldArr(0,renum0)
        mm.setFamilyFieldArr(0,famField0)
        mm.setNameFieldAtLevel(0,namesCellL0)
        mm.setRenumFieldArr(-1,renumM1)
        mm.setFamilyFieldArr(-1,famFieldM1)
        renum1=DataArrayInt([13,16,17,20,21,10,12,11,19,18,15,14,22,23,24,34,33,32,31,30,29,28,27,26,25,45,44,43,42,41,40,39,38,37,36,35])
        famField1=DataArrayInt([-13,-16,-17,-20,-21,-10,-12,-11,-19,-18,-15,-14,-22,-23,-24,-34,-33,-32,-31,-30,-29,-28,-27,-26,-25,-45,-44,-43,-42,-41,-40,-39,-38,-37,-36,-35])
        namesNodes=DataArrayAsciiChar(36,16)
        namesNodes[:]=["Node#%.3d        "%(i) for i in xrange(36)]
        mm.setRenumFieldArr(1,renum1)
        mm.setFamilyFieldArr(1,famField1)
        mm.setNameFieldAtLevel(1,namesNodes)
        mm.setFamilyId("Fam7",77)
        mm.setFamilyId("Fam8",88)
        mm.setGroupsOnFamily("Fam7",["Grp0","Grp1"])
        mm.setGroupsOnFamily("Fam8",["Grp1","Grp2"])
        mm.write(fileName,2)
        #
        mm0=MEDFileUMesh.LoadPartOf(fileName,"Mesh",[NORM_QUAD4],[0,10,1])
        self.assertEqual(mm0.getAllGeoTypes(),[NORM_QUAD4])
        self.assertTrue(mm0.getDirectUndergroundSingleGeoTypeMesh(NORM_QUAD4).getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,2,1,7,8,3,2,8,9,4,3,9,10,5,4,10,11,7,6,12,13,8,7,13,14,9,8,14,15,10,9,15,16,11,10,16,17])))
        coo=DataArrayDouble([(0,0,0),(1,0,0),(2,0,0),(3,0,0),(4,0,0),(5,0,0),(0,1,0),(1,1,0),(2,1,0),(3,1,0),(4,1,0),(5,1,0),(0,2,0),(1,2,0),(2,2,0),(3,2,0),(4,2,0),(5,2,0)]) ; coo.setInfoOnComponents(infos)
        self.assertTrue(mm0.getCoords().isEqual(coo,1e-12))
        self.assertTrue(mm0.getFamilyFieldAtLevel(0).isEqual(famField0[:10]))
        self.assertTrue(mm0.getNumberFieldAtLevel(0).isEqual(renum0[:10]))
        self.assertTrue(mm0.getNameFieldAtLevel(0).isEqual(namesCellL0[:10]))
        self.assertTrue(mm0.getFamilyFieldAtLevel(1).isEqual(famField1[:18]))
        self.assertTrue(mm0.getNumberFieldAtLevel(1).isEqual(renum1[:18]))
        self.assertTrue(mm0.getNameFieldAtLevel(1).isEqual(namesNodes[:18]))
        #
        mm1=MEDFileUMesh.LoadPartOf(fileName,"Mesh",[NORM_QUAD4],[11,25,1])
        self.assertEqual(mm1.getAllGeoTypes(),[NORM_QUAD4])
        self.assertTrue(mm1.getDirectUndergroundSingleGeoTypeMesh(NORM_QUAD4).getNodalConnectivity().isEqual(DataArrayInt([1,0,6,7,2,1,7,8,3,2,8,9,4,3,9,10,6,5,11,12,7,6,12,13,8,7,13,14,9,8,14,15,10,9,15,16,12,11,17,18,13,12,18,19,14,13,19,20,15,14,20,21,16,15,21,22])))
        coo=DataArrayDouble([(1,2,0),(2,2,0),(3,2,0),(4,2,0),(5,2,0),(0,3,0),(1,3,0),(2,3,0),(3,3,0),(4,3,0),(5,3,0),(0,4,0),(1,4,0),(2,4,0),(3,4,0),(4,4,0),(5,4,0),(0,5,0),(1,5,0),(2,5,0),(3,5,0),(4,5,0),(5,5,0)]) ; coo.setInfoOnComponents(infos)
        self.assertTrue(mm1.getCoords().isEqual(coo,1e-12))
        self.assertTrue(mm1.getFamilyFieldAtLevel(0).isEqual(famField0[11:]))
        self.assertTrue(mm1.getNumberFieldAtLevel(0).isEqual(renum0[11:]))
        self.assertTrue(mm1.getNameFieldAtLevel(0).isEqual(namesCellL0[11:]))
        self.assertTrue(mm1.getFamilyFieldAtLevel(1).isEqual(famField1[13:]))
        self.assertTrue(mm1.getNumberFieldAtLevel(1).isEqual(renum1[13:]))
        self.assertTrue(mm1.getNameFieldAtLevel(1).isEqual(namesNodes[13:]))
        #
        mm2=MEDFileUMesh.LoadPartOf(fileName,"Mesh",[NORM_SEG2,NORM_QUAD4],[0,5,1,1,10,1])
        self.assertEqual(mm2.getAllGeoTypes(),[NORM_QUAD4,NORM_SEG2])
        self.assertTrue(mm2.getFamilyFieldAtLevel(0).isEqual(famField0[1:10]))
        self.assertTrue(mm2.getNumberFieldAtLevel(0).isEqual(renum0[1:10]))
        self.assertTrue(mm2.getNameFieldAtLevel(0).isEqual(namesCellL0[1:10]))
        self.assertTrue(mm2.getFamilyFieldAtLevel(-1).isEqual(famFieldM1))
        self.assertTrue(mm2.getNumberFieldAtLevel(-1).isEqual(renumM1))
        self.assertTrue(mm2.getNameFieldAtLevel(-1) is None)
        self.assertTrue(mm2.getDirectUndergroundSingleGeoTypeMesh(NORM_QUAD4).getNodalConnectivity().isEqual(DataArrayInt([2,1,7,8,3,2,8,9,4,3,9,10,5,4,10,11,7,6,12,13,8,7,13,14,9,8,14,15,10,9,15,16,11,10,16,17])))
        self.assertTrue(mm2.getDirectUndergroundSingleGeoTypeMesh(NORM_SEG2).getNodalConnectivity().isEqual(DataArrayInt([0,1,1,2,2,3,3,4,4,5])))
        coo=DataArrayDouble([(0,0,0),(1,0,0),(2,0,0),(3,0,0),(4,0,0),(5,0,0),(0,1,0),(1,1,0),(2,1,0),(3,1,0),(4,1,0),(5,1,0),(0,2,0),(1,2,0),(2,2,0),(3,2,0),(4,2,0),(5,2,0)]) ; coo.setInfoOnComponents(infos)
        self.assertTrue(mm2.getCoords().isEqual(coo,1e-12))
        self.assertTrue(mm2.getFamilyFieldAtLevel(1).isEqual(famField1[:18]))
        self.assertTrue(mm2.getNumberFieldAtLevel(1).isEqual(renum1[:18]))
        self.assertTrue(mm2.getNameFieldAtLevel(1).isEqual(namesNodes[:18]))
        pass

    def testMEDFileFieldsLoadPart1(self):
        """This method tests partial loading on fields on CELL. It is the same principle than those in testMEDFileUMeshLoadPart1.
        """
        fileName="Pyfile82.med"
        meshName="Mesh"
        compos=["aa [kg]","bbb [m/s]"]
        arr=DataArrayDouble(6) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
        m=m.buildUnstructured()
        m.setName(meshName)
        m.changeSpaceDimension(3,0.)
        infos=["aa [b]","cc [de]","gg [klm]"]
        m.getCoords().setInfoOnComponents(infos)
        m.checkCoherency2()
        f=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f.setMesh(m)
        f.setName("Field")
        arr=DataArrayDouble(25,2) ; arr.setInfoOnComponents(compos)
        arr[:,0]=range(25)
        arr[:,1]=range(100,125)
        f.setArray(arr)
        MEDLoader.WriteField(fileName,f,2)
        f=MEDCouplingFieldDouble(ON_NODES,ONE_TIME) ; f.setMesh(m)
        f.setName("FieldNode")
        arr=DataArrayDouble(36,2) ; arr.setInfoOnComponents(compos)
        arr[:,0]=range(200,236)
        arr[:,1]=range(300,336)
        f.setArray(arr)
        f.checkCoherency()
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f)
        #
        ms=MEDFileMeshes()
        mm=MEDFileUMesh.LoadPartOf(fileName,meshName,[NORM_QUAD4],[0,6,1])
        ms.pushMesh(mm)
        fs=MEDFileFields.LoadPartOf(fileName,False,ms)
        self.assertEqual(fs[1][0].getFieldSplitedByType(),[(40,[(1,(0,14),'','')])])
        #
        ms=MEDFileMeshes()
        mm=MEDFileUMesh.LoadPartOf(fileName,meshName,[NORM_QUAD4],[3,15,1])
        ms.pushMesh(mm)
        fs=MEDFileFields.LoadPartOf(fileName,False,ms)
        fs=fs.deepCpy()
        fs[0][0].loadArrays()
        arr=DataArrayDouble(12,2) ; arr[:,0]=range(3,15) ; arr[:,1]=range(103,115)
        arr.setInfoOnComponents(compos)
        self.assertTrue(fs[0][0].getUndergroundDataArray().isEqual(arr,1e-12))
        fs[1][0].loadArrays()
        arr=DataArrayDouble(21,2) ; arr[:,0]=range(203,224) ; arr[:,1]=range(303,324)
        arr.setInfoOnComponents(compos)
        self.assertTrue(fs[1][0].getUndergroundDataArray().isEqual(arr,1e-12))
        pass

    def testMEDFileWithoutCells1(self):
        fileName="Pyfile83.med"
        coo=DataArrayDouble([(0,0,0),(1,0,0),(2,0,0)])
        coo.setInfoOnComponents(["aa [m]","bbb [s]","cccc [m/s]"])
        mm=MEDFileUMesh()
        mm.setCoords(coo)
        mm.setName("mesh")
        mm.write(fileName,2)
        #
        mm=MEDFileMesh.New(fileName)
        self.assertEqual(mm.getName(),"mesh")
        self.assertTrue(mm.getCoords().isEqual(coo,1e-12))
        pass

    def testZipCoordsWithLoadPart1(self):
        """ Test close to Pyfile82.med except that here zipCoords on MEDFileUMesh is invoked here to see if the PartDef is correctly updated.
        """
        fileName="Pyfile84.med"
        meshName="Mesh"
        compos=["aa [kg]","bbb [m/s]"]
        arr=DataArrayDouble(6) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
        m=m.buildUnstructured()
        m.setName(meshName)
        m.changeSpaceDimension(3,0.)
        infos=["aa [b]","cc [de]","gg [klm]"]
        m.getCoords().setInfoOnComponents(infos)
        m.checkCoherency2()
        f=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f.setMesh(m)
        f.setName("Field")
        arr=DataArrayDouble(25,2) ; arr.setInfoOnComponents(compos)
        arr[:,0]=range(25)
        arr[:,1]=range(100,125)
        f.setArray(arr)
        MEDLoader.WriteField(fileName,f,2)
        f=MEDCouplingFieldDouble(ON_NODES,ONE_TIME) ; f.setMesh(m)
        f.setName("FieldNode")
        arr=DataArrayDouble(36,2) ; arr.setInfoOnComponents(compos)
        arr[:,0]=range(200,236)
        arr[:,1]=range(300,336)
        f.setArray(arr)
        f.checkCoherency()
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fileName,f)
        #
        ms=MEDFileMeshes()
        mm=MEDFileUMesh.LoadPartOf(fileName,meshName,[NORM_QUAD4],[4,6,1])
        ms.pushMesh(mm)
        spd=mm.getPartDefAtLevel(0,NORM_QUAD4)
        self.assertEqual(spd.getSlice(),slice(4,6,1))
        spd=mm.getPartDefAtLevel(1)
        self.assertEqual(spd.getSlice(),slice(4,14,1))
        self.assertTrue(spd.getNumberOfElems()==10 and spd.getNumberOfElems()==mm.getNumberOfNodes())
        mm.zipCoords() # <- The important line is here !
        spd=mm.getPartDefAtLevel(0,NORM_QUAD4)
        self.assertEqual(spd.getSlice(),slice(4,6,1))
        spd=mm.getPartDefAtLevel(1)
        self.assertTrue(spd.getNumberOfElems()==8 and spd.getNumberOfElems()==mm.getNumberOfNodes())
        self.assertTrue(spd.toDAI().isEqual(DataArrayInt([4,5,6,7,10,11,12,13])))
        fs=MEDFileFields.LoadPartOf(fileName,False,ms)
        fs[0][0].loadArrays()
        arr=DataArrayDouble([(4,104),(5,105)])
        arr.setInfoOnComponents(compos)
        self.assertTrue(fs[0][0].getUndergroundDataArray().isEqual(arr,1e-12))
        fs[1][0].loadArrays()
        arr=DataArrayDouble([(204,304),(205,305),(206,306),(207,307),(210,310),(211,311),(212,312),(213,313)])
        arr.setInfoOnComponents(compos)
        self.assertTrue(fs[1][0].getUndergroundDataArray().isEqual(arr,1e-12))
        pass
    
    def testMEDFileCMeshSetGroupsAtLevel(self):
        """ Non regression test to check that setGroupsAtLevel is available with MEDFileCMesh.
        """
        m=MEDCouplingCMesh() ; m.setCoords(DataArrayDouble([0,1,2,3,4]),DataArrayDouble([0,1,2,3,4]))
        m.setName("Mesh")
        mm=MEDFileCMesh() ; mm.setMesh(m)
        grp=DataArrayInt([1,3,4,5,7]) ; grp.setName("MyAssembly")
        mm.setGroupsAtLevel(0,[grp])
        self.assertTrue(mm.getFamilyFieldAtLevel(0).isEqual(DataArrayInt([-1,-2,-1,-2,-2,-2,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1])))
        pass

    def testMEDFileUMeshBuildExtrudedMesh1(self):
        """ New functionality of MEDFileUMesh.buildExtrudedMesh."""
        fileName="PyFile85.med"
        meshName2D="Mesh"
        meshName1D="Mesh1D"
        meshName3DOut="Mesh3D"
        #
        d1=DataArrayInt([0,4,20,24])
        d2=DataArrayInt([0,1,2,3,7,8,12,13,17,18,19,20])
        #
        a=DataArrayDouble(6) ; a.iota()
        m=MEDCouplingCMesh() ; m.setCoords(a,a)
        m=m.buildUnstructured()
        d1c=d1.buildComplement(m.getNumberOfCells())
        m=m[d1c] ; m.zipCoords()
        m0=m[d2] ; m1=m[d2.buildComplement(m.getNumberOfCells())]
        m0.simplexize(0)
        m=MEDCouplingUMesh.MergeUMeshesOnSameCoords([m0,m1])
        m.setName(meshName2D)
        mMinus1,a,b,c,d=m.buildDescendingConnectivity()
        e=d.deltaShiftIndex().getIdsEqual(1)
        #
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m) ; mm.setMeshAtLevel(-1,mMinus1)
        grp0=DataArrayInt([0,1,2,3,4,5,24,25,26]) ; grp0.setName("grp0")
        mm.setGroupsAtLevel(0,[grp0])
        grp1=e ; grp1.setName("grp1")
        mm.setGroupsAtLevel(-1,[grp1])
        mm.write(fileName,2)
        #
        a=DataArrayDouble(10) ; a.iota()
        tmp=MEDCouplingCMesh() ; tmp.setCoords(a) ; tmp=tmp.buildUnstructured()
        tmp.setName(meshName1D)
        tmp.changeSpaceDimension(3)
        tmp.setCoords(tmp.getCoords()[:,[1,2,0]])
        mm1D=MEDFileUMesh()
        mm1D.setMeshAtLevel(0,tmp)
        mm1D.write(fileName,0)
        # test is here !
        mm2D=MEDFileMesh.New(fileName,meshName2D)
        mm1D=MEDFileMesh.New(fileName,meshName1D)
        m1D=mm1D.getMeshAtLevel(0)
        mm3D=mm2D.buildExtrudedMesh(m1D,0)
        #
        self.assertEqual(mm3D.getName(),mm2D.getName())
        self.assertEqual(mm3D.getNumberOfCellsAtLevel(0),297)
        self.assertEqual(mm3D.getNumberOfCellsAtLevel(-1),576)
        self.assertEqual(mm3D.getGroupsNames(),('grp0','grp1'))
        self.assertEqual(mm3D.getGrpNonEmptyLevels("grp0"),(0,))
        self.assertEqual(mm3D.getGrpNonEmptyLevels("grp1"),(-1,))
        d=DataArrayDouble([(1.,0.,0.),(2.,0.,0.),(3.,0.,0.),(4.,0.,0.),(0.,1.,0.),(1.,1.,0.),(2.,1.,0.),(3.,1.,0.),(4.,1.,0.),(5.,1.,0.),(0.,2.,0.),(1.,2.,0.),(2.,2.,0.),(3.,2.,0.),(4.,2.,0.),(5.,2.,0.),(0.,3.,0.),(1.,3.,0.),(2.,3.,0.),(3.,3.,0.),(4.,3.,0.),(5.,3.,0.),(0.,4.,0.),(1.,4.,0.),(2.,4.,0.),(3.,4.,0.),(4.,4.,0.),(5.,4.,0.),(1.,5.,0.),(2.,5.,0.),(3.,5.,0.),(4.,5.,0.),(1.,0.,1.),(2.,0.,1.),(3.,0.,1.),(4.,0.,1.),(0.,1.,1.),(1.,1.,1.),(2.,1.,1.),(3.,1.,1.),(4.,1.,1.),(5.,1.,1.),(0.,2.,1.),(1.,2.,1.),(2.,2.,1.),(3.,2.,1.),(4.,2.,1.),(5.,2.,1.),(0.,3.,1.),(1.,3.,1.),(2.,3.,1.),(3.,3.,1.),(4.,3.,1.),(5.,3.,1.),(0.,4.,1.),(1.,4.,1.),(2.,4.,1.),(3.,4.,1.),(4.,4.,1.),(5.,4.,1.),(1.,5.,1.),(2.,5.,1.),(3.,5.,1.),(4.,5.,1.),(1.,0.,2.),(2.,0.,2.),(3.,0.,2.),(4.,0.,2.),(0.,1.,2.),(1.,1.,2.),(2.,1.,2.),(3.,1.,2.),(4.,1.,2.),(5.,1.,2.),(0.,2.,2.),(1.,2.,2.),(2.,2.,2.),(3.,2.,2.),(4.,2.,2.),(5.,2.,2.),(0.,3.,2.),(1.,3.,2.),(2.,3.,2.),(3.,3.,2.),(4.,3.,2.),(5.,3.,2.),(0.,4.,2.),(1.,4.,2.),(2.,4.,2.),(3.,4.,2.),(4.,4.,2.),(5.,4.,2.),(1.,5.,2.),(2.,5.,2.),(3.,5.,2.),(4.,5.,2.),(1.,0.,3.),(2.,0.,3.),(3.,0.,3.),(4.,0.,3.),(0.,1.,3.),(1.,1.,3.),(2.,1.,3.),(3.,1.,3.),(4.,1.,3.),(5.,1.,3.),(0.,2.,3.),(1.,2.,3.),(2.,2.,3.),(3.,2.,3.),(4.,2.,3.),(5.,2.,3.),(0.,3.,3.),(1.,3.,3.),(2.,3.,3.),(3.,3.,3.),(4.,3.,3.),(5.,3.,3.),(0.,4.,3.),(1.,4.,3.),(2.,4.,3.),(3.,4.,3.),(4.,4.,3.),(5.,4.,3.),(1.,5.,3.),(2.,5.,3.),(3.,5.,3.),(4.,5.,3.),(1.,0.,4.),(2.,0.,4.),(3.,0.,4.),(4.,0.,4.),(0.,1.,4.),(1.,1.,4.),(2.,1.,4.),(3.,1.,4.),(4.,1.,4.),(5.,1.,4.),(0.,2.,4.),(1.,2.,4.),(2.,2.,4.),(3.,2.,4.),(4.,2.,4.),(5.,2.,4.),(0.,3.,4.),(1.,3.,4.),(2.,3.,4.),(3.,3.,4.),(4.,3.,4.),(5.,3.,4.),(0.,4.,4.),(1.,4.,4.),(2.,4.,4.),(3.,4.,4.),(4.,4.,4.),(5.,4.,4.),(1.,5.,4.),(2.,5.,4.),(3.,5.,4.),(4.,5.,4.),(1.,0.,5.),(2.,0.,5.),(3.,0.,5.),(4.,0.,5.),(0.,1.,5.),(1.,1.,5.),(2.,1.,5.),(3.,1.,5.),(4.,1.,5.),(5.,1.,5.),(0.,2.,5.),(1.,2.,5.),(2.,2.,5.),(3.,2.,5.),(4.,2.,5.),(5.,2.,5.),(0.,3.,5.),(1.,3.,5.),(2.,3.,5.),(3.,3.,5.),(4.,3.,5.),(5.,3.,5.),(0.,4.,5.),(1.,4.,5.),(2.,4.,5.),(3.,4.,5.),(4.,4.,5.),(5.,4.,5.),(1.,5.,5.),(2.,5.,5.),(3.,5.,5.),(4.,5.,5.),(1.,0.,6.),(2.,0.,6.),(3.,0.,6.),(4.,0.,6.),(0.,1.,6.),(1.,1.,6.),(2.,1.,6.),(3.,1.,6.),(4.,1.,6.),(5.,1.,6.),(0.,2.,6.),(1.,2.,6.),(2.,2.,6.),(3.,2.,6.),(4.,2.,6.),(5.,2.,6.),(0.,3.,6.),(1.,3.,6.),(2.,3.,6.),(3.,3.,6.),(4.,3.,6.),(5.,3.,6.),(0.,4.,6.),(1.,4.,6.),(2.,4.,6.),(3.,4.,6.),(4.,4.,6.),(5.,4.,6.),(1.,5.,6.),(2.,5.,6.),(3.,5.,6.),(4.,5.,6.),(1.,0.,7.),(2.,0.,7.),(3.,0.,7.),(4.,0.,7.),(0.,1.,7.),(1.,1.,7.),(2.,1.,7.),(3.,1.,7.),(4.,1.,7.),(5.,1.,7.),(0.,2.,7.),(1.,2.,7.),(2.,2.,7.),(3.,2.,7.),(4.,2.,7.),(5.,2.,7.),(0.,3.,7.),(1.,3.,7.),(2.,3.,7.),(3.,3.,7.),(4.,3.,7.),(5.,3.,7.),(0.,4.,7.),(1.,4.,7.),(2.,4.,7.),(3.,4.,7.),(4.,4.,7.),(5.,4.,7.),(1.,5.,7.),(2.,5.,7.),(3.,5.,7.),(4.,5.,7.),(1.,0.,8.),(2.,0.,8.),(3.,0.,8.),(4.,0.,8.),(0.,1.,8.),(1.,1.,8.),(2.,1.,8.),(3.,1.,8.),(4.,1.,8.),(5.,1.,8.),(0.,2.,8.),(1.,2.,8.),(2.,2.,8.),(3.,2.,8.),(4.,2.,8.),(5.,2.,8.),(0.,3.,8.),(1.,3.,8.),(2.,3.,8.),(3.,3.,8.),(4.,3.,8.),(5.,3.,8.),(0.,4.,8.),(1.,4.,8.),(2.,4.,8.),(3.,4.,8.),(4.,4.,8.),(5.,4.,8.),(1.,5.,8.),(2.,5.,8.),(3.,5.,8.),(4.,5.,8.),(1.,0.,9.),(2.,0.,9.),(3.,0.,9.),(4.,0.,9.),(0.,1.,9.),(1.,1.,9.),(2.,1.,9.),(3.,1.,9.),(4.,1.,9.),(5.,1.,9.),(0.,2.,9.),(1.,2.,9.),(2.,2.,9.),(3.,2.,9.),(4.,2.,9.),(5.,2.,9.),(0.,3.,9.),(1.,3.,9.),(2.,3.,9.),(3.,3.,9.),(4.,3.,9.),(5.,3.,9.),(0.,4.,9.),(1.,4.,9.),(2.,4.,9.),(3.,4.,9.),(4.,4.,9.),(5.,4.,9.),(1.,5.,9.),(2.,5.,9.),(3.,5.,9.),(4.,5.,9.)])
        self.assertTrue(mm3D.getCoords().isEqual(d,1e-12))
        d=DataArrayInt([16,1,0,5,33,32,37,16,1,5,6,33,37,38,16,2,1,6,34,33,38,16,2,6,7,34,38,39,16,3,2,7,35,34,39,16,3,7,8,35,39,40,16,5,4,10,37,36,42,16,5,10,11,37,42,43,16,9,8,14,41,40,46,16,9,14,15,41,46,47,16,11,10,16,43,42,48,16,11,16,17,43,48,49,16,15,14,20,47,46,52,16,15,20,21,47,52,53,16,17,16,22,49,48,54,16,17,22,23,49,54,55,16,21,20,26,53,52,58,16,21,26,27,53,58,59,16,24,23,28,56,55,60,16,24,28,29,56,60,61,16,25,24,29,57,56,61,16,25,29,30,57,61,62,16,26,25,30,58,57,62,16,26,30,31,58,62,63,16,33,32,37,65,64,69,16,33,37,38,65,69,70,16,34,33,38,66,65,70,16,34,38,39,66,70,71,16,35,34,39,67,66,71,16,35,39,40,67,71,72,16,37,36,42,69,68,74,16,37,42,43,69,74,75,16,41,40,46,73,72,78,16,41,46,47,73,78,79,16,43,42,48,75,74,80,16,43,48,49,75,80,81,16,47,46,52,79,78,84,16,47,52,53,79,84,85,16,49,48,54,81,80,86,16,49,54,55,81,86,87,16,53,52,58,85,84,90,16,53,58,59,85,90,91,16,56,55,60,88,87,92,16,56,60,61,88,92,93,16,57,56,61,89,88,93,16,57,61,62,89,93,94,16,58,57,62,90,89,94,16,58,62,63,90,94,95,16,65,64,69,97,96,101,16,65,69,70,97,101,102,16,66,65,70,98,97,102,16,66,70,71,98,102,103,16,67,66,71,99,98,103,16,67,71,72,99,103,104,16,69,68,74,101,100,106,16,69,74,75,101,106,107,16,73,72,78,105,104,110,16,73,78,79,105,110,111,16,75,74,80,107,106,112,16,75,80,81,107,112,113,16,79,78,84,111,110,116,16,79,84,85,111,116,117,16,81,80,86,113,112,118,16,81,86,87,113,118,119,16,85,84,90,117,116,122,16,85,90,91,117,122,123,16,88,87,92,120,119,124,16,88,92,93,120,124,125,16,89,88,93,121,120,125,16,89,93,94,121,125,126,16,90,89,94,122,121,126,16,90,94,95,122,126,127,16,97,96,101,129,128,133,16,97,101,102,129,133,134,16,98,97,102,130,129,134,16,98,102,103,130,134,135,16,99,98,103,131,130,135,16,99,103,104,131,135,136,16,101,100,106,133,132,138,16,101,106,107,133,138,139,16,105,104,110,137,136,142,16,105,110,111,137,142,143,16,107,106,112,139,138,144,16,107,112,113,139,144,145,16,111,110,116,143,142,148,16,111,116,117,143,148,149,16,113,112,118,145,144,150,16,113,118,119,145,150,151,16,117,116,122,149,148,154,16,117,122,123,149,154,155,16,120,119,124,152,151,156,16,120,124,125,152,156,157,16,121,120,125,153,152,157,16,121,125,126,153,157,158,16,122,121,126,154,153,158,16,122,126,127,154,158,159,16,129,128,133,161,160,165,16,129,133,134,161,165,166,16,130,129,134,162,161,166,16,130,134,135,162,166,167,16,131,130,135,163,162,167,16,131,135,136,163,167,168,16,133,132,138,165,164,170,16,133,138,139,165,170,171,16,137,136,142,169,168,174,16,137,142,143,169,174,175,16,139,138,144,171,170,176,16,139,144,145,171,176,177,16,143,142,148,175,174,180,16,143,148,149,175,180,181,16,145,144,150,177,176,182,16,145,150,151,177,182,183,16,149,148,154,181,180,186,16,149,154,155,181,186,187,16,152,151,156,184,183,188,16,152,156,157,184,188,189,16,153,152,157,185,184,189,16,153,157,158,185,189,190,16,154,153,158,186,185,190,16,154,158,159,186,190,191,16,161,160,165,193,192,197,16,161,165,166,193,197,198,16,162,161,166,194,193,198,16,162,166,167,194,198,199,16,163,162,167,195,194,199,16,163,167,168,195,199,200,16,165,164,170,197,196,202,16,165,170,171,197,202,203,16,169,168,174,201,200,206,16,169,174,175,201,206,207,16,171,170,176,203,202,208,16,171,176,177,203,208,209,16,175,174,180,207,206,212,16,175,180,181,207,212,213,16,177,176,182,209,208,214,16,177,182,183,209,214,215,16,181,180,186,213,212,218,16,181,186,187,213,218,219,16,184,183,188,216,215,220,16,184,188,189,216,220,221,16,185,184,189,217,216,221,16,185,189,190,217,221,222,16,186,185,190,218,217,222,16,186,190,191,218,222,223,16,193,192,197,225,224,229,16,193,197,198,225,229,230,16,194,193,198,226,225,230,16,194,198,199,226,230,231,16,195,194,199,227,226,231,16,195,199,200,227,231,232,16,197,196,202,229,228,234,16,197,202,203,229,234,235,16,201,200,206,233,232,238,16,201,206,207,233,238,239,16,203,202,208,235,234,240,16,203,208,209,235,240,241,16,207,206,212,239,238,244,16,207,212,213,239,244,245,16,209,208,214,241,240,246,16,209,214,215,241,246,247,16,213,212,218,245,244,250,16,213,218,219,245,250,251,16,216,215,220,248,247,252,16,216,220,221,248,252,253,16,217,216,221,249,248,253,16,217,221,222,249,253,254,16,218,217,222,250,249,254,16,218,222,223,250,254,255,16,225,224,229,257,256,261,16,225,229,230,257,261,262,16,226,225,230,258,257,262,16,226,230,231,258,262,263,16,227,226,231,259,258,263,16,227,231,232,259,263,264,16,229,228,234,261,260,266,16,229,234,235,261,266,267,16,233,232,238,265,264,270,16,233,238,239,265,270,271,16,235,234,240,267,266,272,16,235,240,241,267,272,273,16,239,238,244,271,270,276,16,239,244,245,271,276,277,16,241,240,246,273,272,278,16,241,246,247,273,278,279,16,245,244,250,277,276,282,16,245,250,251,277,282,283,16,248,247,252,280,279,284,16,248,252,253,280,284,285,16,249,248,253,281,280,285,16,249,253,254,281,285,286,16,250,249,254,282,281,286,16,250,254,255,282,286,287,16,257,256,261,289,288,293,16,257,261,262,289,293,294,16,258,257,262,290,289,294,16,258,262,263,290,294,295,16,259,258,263,291,290,295,16,259,263,264,291,295,296,16,261,260,266,293,292,298,16,261,266,267,293,298,299,16,265,264,270,297,296,302,16,265,270,271,297,302,303,16,267,266,272,299,298,304,16,267,272,273,299,304,305,16,271,270,276,303,302,308,16,271,276,277,303,308,309,16,273,272,278,305,304,310,16,273,278,279,305,310,311,16,277,276,282,309,308,314,16,277,282,283,309,314,315,16,280,279,284,312,311,316,16,280,284,285,312,316,317,16,281,280,285,313,312,317,16,281,285,286,313,317,318,16,282,281,286,314,313,318,16,282,286,287,314,318,319,18,6,5,11,12,38,37,43,44,18,7,6,12,13,39,38,44,45,18,8,7,13,14,40,39,45,46,18,12,11,17,18,44,43,49,50,18,13,12,18,19,45,44,50,51,18,14,13,19,20,46,45,51,52,18,18,17,23,24,50,49,55,56,18,19,18,24,25,51,50,56,57,18,20,19,25,26,52,51,57,58,18,38,37,43,44,70,69,75,76,18,39,38,44,45,71,70,76,77,18,40,39,45,46,72,71,77,78,18,44,43,49,50,76,75,81,82,18,45,44,50,51,77,76,82,83,18,46,45,51,52,78,77,83,84,18,50,49,55,56,82,81,87,88,18,51,50,56,57,83,82,88,89,18,52,51,57,58,84,83,89,90,18,70,69,75,76,102,101,107,108,18,71,70,76,77,103,102,108,109,18,72,71,77,78,104,103,109,110,18,76,75,81,82,108,107,113,114,18,77,76,82,83,109,108,114,115,18,78,77,83,84,110,109,115,116,18,82,81,87,88,114,113,119,120,18,83,82,88,89,115,114,120,121,18,84,83,89,90,116,115,121,122,18,102,101,107,108,134,133,139,140,18,103,102,108,109,135,134,140,141,18,104,103,109,110,136,135,141,142,18,108,107,113,114,140,139,145,146,18,109,108,114,115,141,140,146,147,18,110,109,115,116,142,141,147,148,18,114,113,119,120,146,145,151,152,18,115,114,120,121,147,146,152,153,18,116,115,121,122,148,147,153,154,18,134,133,139,140,166,165,171,172,18,135,134,140,141,167,166,172,173,18,136,135,141,142,168,167,173,174,18,140,139,145,146,172,171,177,178,18,141,140,146,147,173,172,178,179,18,142,141,147,148,174,173,179,180,18,146,145,151,152,178,177,183,184,18,147,146,152,153,179,178,184,185,18,148,147,153,154,180,179,185,186,18,166,165,171,172,198,197,203,204,18,167,166,172,173,199,198,204,205,18,168,167,173,174,200,199,205,206,18,172,171,177,178,204,203,209,210,18,173,172,178,179,205,204,210,211,18,174,173,179,180,206,205,211,212,18,178,177,183,184,210,209,215,216,18,179,178,184,185,211,210,216,217,18,180,179,185,186,212,211,217,218,18,198,197,203,204,230,229,235,236,18,199,198,204,205,231,230,236,237,18,200,199,205,206,232,231,237,238,18,204,203,209,210,236,235,241,242,18,205,204,210,211,237,236,242,243,18,206,205,211,212,238,237,243,244,18,210,209,215,216,242,241,247,248,18,211,210,216,217,243,242,248,249,18,212,211,217,218,244,243,249,250,18,230,229,235,236,262,261,267,268,18,231,230,236,237,263,262,268,269,18,232,231,237,238,264,263,269,270,18,236,235,241,242,268,267,273,274,18,237,236,242,243,269,268,274,275,18,238,237,243,244,270,269,275,276,18,242,241,247,248,274,273,279,280,18,243,242,248,249,275,274,280,281,18,244,243,249,250,276,275,281,282,18,262,261,267,268,294,293,299,300,18,263,262,268,269,295,294,300,301,18,264,263,269,270,296,295,301,302,18,268,267,273,274,300,299,305,306,18,269,268,274,275,301,300,306,307,18,270,269,275,276,302,301,307,308,18,274,273,279,280,306,305,311,312,18,275,274,280,281,307,306,312,313,18,276,275,281,282,308,307,313,314])
        self.assertTrue(mm3D.getMeshAtLevel(0).getNodalConnectivity().isEqual(d))
        d=DataArrayInt([0,7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140,147,154,161,168,175,182,189,196,203,210,217,224,231,238,245,252,259,266,273,280,287,294,301,308,315,322,329,336,343,350,357,364,371,378,385,392,399,406,413,420,427,434,441,448,455,462,469,476,483,490,497,504,511,518,525,532,539,546,553,560,567,574,581,588,595,602,609,616,623,630,637,644,651,658,665,672,679,686,693,700,707,714,721,728,735,742,749,756,763,770,777,784,791,798,805,812,819,826,833,840,847,854,861,868,875,882,889,896,903,910,917,924,931,938,945,952,959,966,973,980,987,994,1001,1008,1015,1022,1029,1036,1043,1050,1057,1064,1071,1078,1085,1092,1099,1106,1113,1120,1127,1134,1141,1148,1155,1162,1169,1176,1183,1190,1197,1204,1211,1218,1225,1232,1239,1246,1253,1260,1267,1274,1281,1288,1295,1302,1309,1316,1323,1330,1337,1344,1351,1358,1365,1372,1379,1386,1393,1400,1407,1414,1421,1428,1435,1442,1449,1456,1463,1470,1477,1484,1491,1498,1505,1512,1521,1530,1539,1548,1557,1566,1575,1584,1593,1602,1611,1620,1629,1638,1647,1656,1665,1674,1683,1692,1701,1710,1719,1728,1737,1746,1755,1764,1773,1782,1791,1800,1809,1818,1827,1836,1845,1854,1863,1872,1881,1890,1899,1908,1917,1926,1935,1944,1953,1962,1971,1980,1989,1998,2007,2016,2025,2034,2043,2052,2061,2070,2079,2088,2097,2106,2115,2124,2133,2142,2151,2160,2169,2178,2187,2196,2205,2214,2223,2232,2241])
        self.assertTrue(mm3D.getMeshAtLevel(0).getNodalConnectivityIndex().isEqual(d))
        d=DataArrayInt([4,1,0,32,33,4,0,5,37,32,4,5,1,33,37,4,5,6,38,37,4,6,1,33,38,4,2,1,33,34,4,6,2,34,38,4,6,7,39,38,4,7,2,34,39,4,3,2,34,35,4,7,3,35,39,4,7,8,40,39,4,8,3,35,40,4,5,4,36,37,4,4,10,42,36,4,10,5,37,42,4,10,11,43,42,4,11,5,37,43,4,9,8,40,41,4,8,14,46,40,4,14,9,41,46,4,14,15,47,46,4,15,9,41,47,4,10,16,48,42,4,16,11,43,48,4,16,17,49,48,4,17,11,43,49,4,14,20,52,46,4,20,15,47,52,4,20,21,53,52,4,21,15,47,53,4,16,22,54,48,4,22,17,49,54,4,22,23,55,54,4,23,17,49,55,4,20,26,58,52,4,26,21,53,58,4,26,27,59,58,4,27,21,53,59,4,24,23,55,56,4,23,28,60,55,4,28,24,56,60,4,28,29,61,60,4,29,24,56,61,4,25,24,56,57,4,29,25,57,61,4,29,30,62,61,4,30,25,57,62,4,26,25,57,58,4,30,26,58,62,4,30,31,63,62,4,31,26,58,63,4,11,12,44,43,4,12,6,38,44,4,12,13,45,44,4,13,7,39,45,4,13,14,46,45,4,17,18,50,49,4,18,12,44,50,4,18,19,51,50,4,19,13,45,51,4,19,20,52,51,4,24,18,50,56,4,25,19,51,57,4,33,32,64,65,4,32,37,69,64,4,37,33,65,69,4,37,38,70,69,4,38,33,65,70,4,34,33,65,66,4,38,34,66,70,4,38,39,71,70,4,39,34,66,71,4,35,34,66,67,4,39,35,67,71,4,39,40,72,71,4,40,35,67,72,4,37,36,68,69,4,36,42,74,68,4,42,37,69,74,4,42,43,75,74,4,43,37,69,75,4,41,40,72,73,4,40,46,78,72,4,46,41,73,78,4,46,47,79,78,4,47,41,73,79,4,42,48,80,74,4,48,43,75,80,4,48,49,81,80,4,49,43,75,81,4,46,52,84,78,4,52,47,79,84,4,52,53,85,84,4,53,47,79,85,4,48,54,86,80,4,54,49,81,86,4,54,55,87,86,4,55,49,81,87,4,52,58,90,84,4,58,53,85,90,4,58,59,91,90,4,59,53,85,91,4,56,55,87,88,4,55,60,92,87,4,60,56,88,92,4,60,61,93,92,4,61,56,88,93,4,57,56,88,89,4,61,57,89,93,4,61,62,94,93,4,62,57,89,94,4,58,57,89,90,4,62,58,90,94,4,62,63,95,94,4,63,58,90,95,4,43,44,76,75,4,44,38,70,76,4,44,45,77,76,4,45,39,71,77,4,45,46,78,77,4,49,50,82,81,4,50,44,76,82,4,50,51,83,82,4,51,45,77,83,4,51,52,84,83,4,56,50,82,88,4,57,51,83,89,4,65,64,96,97,4,64,69,101,96,4,69,65,97,101,4,69,70,102,101,4,70,65,97,102,4,66,65,97,98,4,70,66,98,102,4,70,71,103,102,4,71,66,98,103,4,67,66,98,99,4,71,67,99,103,4,71,72,104,103,4,72,67,99,104,4,69,68,100,101,4,68,74,106,100,4,74,69,101,106,4,74,75,107,106,4,75,69,101,107,4,73,72,104,105,4,72,78,110,104,4,78,73,105,110,4,78,79,111,110,4,79,73,105,111,4,74,80,112,106,4,80,75,107,112,4,80,81,113,112,4,81,75,107,113,4,78,84,116,110,4,84,79,111,116,4,84,85,117,116,4,85,79,111,117,4,80,86,118,112,4,86,81,113,118,4,86,87,119,118,4,87,81,113,119,4,84,90,122,116,4,90,85,117,122,4,90,91,123,122,4,91,85,117,123,4,88,87,119,120,4,87,92,124,119,4,92,88,120,124,4,92,93,125,124,4,93,88,120,125,4,89,88,120,121,4,93,89,121,125,4,93,94,126,125,4,94,89,121,126,4,90,89,121,122,4,94,90,122,126,4,94,95,127,126,4,95,90,122,127,4,75,76,108,107,4,76,70,102,108,4,76,77,109,108,4,77,71,103,109,4,77,78,110,109,4,81,82,114,113,4,82,76,108,114,4,82,83,115,114,4,83,77,109,115,4,83,84,116,115,4,88,82,114,120,4,89,83,115,121,4,97,96,128,129,4,96,101,133,128,4,101,97,129,133,4,101,102,134,133,4,102,97,129,134,4,98,97,129,130,4,102,98,130,134,4,102,103,135,134,4,103,98,130,135,4,99,98,130,131,4,103,99,131,135,4,103,104,136,135,4,104,99,131,136,4,101,100,132,133,4,100,106,138,132,4,106,101,133,138,4,106,107,139,138,4,107,101,133,139,4,105,104,136,137,4,104,110,142,136,4,110,105,137,142,4,110,111,143,142,4,111,105,137,143,4,106,112,144,138,4,112,107,139,144,4,112,113,145,144,4,113,107,139,145,4,110,116,148,142,4,116,111,143,148,4,116,117,149,148,4,117,111,143,149,4,112,118,150,144,4,118,113,145,150,4,118,119,151,150,4,119,113,145,151,4,116,122,154,148,4,122,117,149,154,4,122,123,155,154,4,123,117,149,155,4,120,119,151,152,4,119,124,156,151,4,124,120,152,156,4,124,125,157,156,4,125,120,152,157,4,121,120,152,153,4,125,121,153,157,4,125,126,158,157,4,126,121,153,158,4,122,121,153,154,4,126,122,154,158,4,126,127,159,158,4,127,122,154,159,4,107,108,140,139,4,108,102,134,140,4,108,109,141,140,4,109,103,135,141,4,109,110,142,141,4,113,114,146,145,4,114,108,140,146,4,114,115,147,146,4,115,109,141,147,4,115,116,148,147,4,120,114,146,152,4,121,115,147,153,4,129,128,160,161,4,128,133,165,160,4,133,129,161,165,4,133,134,166,165,4,134,129,161,166,4,130,129,161,162,4,134,130,162,166,4,134,135,167,166,4,135,130,162,167,4,131,130,162,163,4,135,131,163,167,4,135,136,168,167,4,136,131,163,168,4,133,132,164,165,4,132,138,170,164,4,138,133,165,170,4,138,139,171,170,4,139,133,165,171,4,137,136,168,169,4,136,142,174,168,4,142,137,169,174,4,142,143,175,174,4,143,137,169,175,4,138,144,176,170,4,144,139,171,176,4,144,145,177,176,4,145,139,171,177,4,142,148,180,174,4,148,143,175,180,4,148,149,181,180,4,149,143,175,181,4,144,150,182,176,4,150,145,177,182,4,150,151,183,182,4,151,145,177,183,4,148,154,186,180,4,154,149,181,186,4,154,155,187,186,4,155,149,181,187,4,152,151,183,184,4,151,156,188,183,4,156,152,184,188,4,156,157,189,188,4,157,152,184,189,4,153,152,184,185,4,157,153,185,189,4,157,158,190,189,4,158,153,185,190,4,154,153,185,186,4,158,154,186,190,4,158,159,191,190,4,159,154,186,191,4,139,140,172,171,4,140,134,166,172,4,140,141,173,172,4,141,135,167,173,4,141,142,174,173,4,145,146,178,177,4,146,140,172,178,4,146,147,179,178,4,147,141,173,179,4,147,148,180,179,4,152,146,178,184,4,153,147,179,185,4,161,160,192,193,4,160,165,197,192,4,165,161,193,197,4,165,166,198,197,4,166,161,193,198,4,162,161,193,194,4,166,162,194,198,4,166,167,199,198,4,167,162,194,199,4,163,162,194,195,4,167,163,195,199,4,167,168,200,199,4,168,163,195,200,4,165,164,196,197,4,164,170,202,196,4,170,165,197,202,4,170,171,203,202,4,171,165,197,203,4,169,168,200,201,4,168,174,206,200,4,174,169,201,206,4,174,175,207,206,4,175,169,201,207,4,170,176,208,202,4,176,171,203,208,4,176,177,209,208,4,177,171,203,209,4,174,180,212,206,4,180,175,207,212,4,180,181,213,212,4,181,175,207,213,4,176,182,214,208,4,182,177,209,214,4,182,183,215,214,4,183,177,209,215,4,180,186,218,212,4,186,181,213,218,4,186,187,219,218,4,187,181,213,219,4,184,183,215,216,4,183,188,220,215,4,188,184,216,220,4,188,189,221,220,4,189,184,216,221,4,185,184,216,217,4,189,185,217,221,4,189,190,222,221,4,190,185,217,222,4,186,185,217,218,4,190,186,218,222,4,190,191,223,222,4,191,186,218,223,4,171,172,204,203,4,172,166,198,204,4,172,173,205,204,4,173,167,199,205,4,173,174,206,205,4,177,178,210,209,4,178,172,204,210,4,178,179,211,210,4,179,173,205,211,4,179,180,212,211,4,184,178,210,216,4,185,179,211,217,4,193,192,224,225,4,192,197,229,224,4,197,193,225,229,4,197,198,230,229,4,198,193,225,230,4,194,193,225,226,4,198,194,226,230,4,198,199,231,230,4,199,194,226,231,4,195,194,226,227,4,199,195,227,231,4,199,200,232,231,4,200,195,227,232,4,197,196,228,229,4,196,202,234,228,4,202,197,229,234,4,202,203,235,234,4,203,197,229,235,4,201,200,232,233,4,200,206,238,232,4,206,201,233,238,4,206,207,239,238,4,207,201,233,239,4,202,208,240,234,4,208,203,235,240,4,208,209,241,240,4,209,203,235,241,4,206,212,244,238,4,212,207,239,244,4,212,213,245,244,4,213,207,239,245,4,208,214,246,240,4,214,209,241,246,4,214,215,247,246,4,215,209,241,247,4,212,218,250,244,4,218,213,245,250,4,218,219,251,250,4,219,213,245,251,4,216,215,247,248,4,215,220,252,247,4,220,216,248,252,4,220,221,253,252,4,221,216,248,253,4,217,216,248,249,4,221,217,249,253,4,221,222,254,253,4,222,217,249,254,4,218,217,249,250,4,222,218,250,254,4,222,223,255,254,4,223,218,250,255,4,203,204,236,235,4,204,198,230,236,4,204,205,237,236,4,205,199,231,237,4,205,206,238,237,4,209,210,242,241,4,210,204,236,242,4,210,211,243,242,4,211,205,237,243,4,211,212,244,243,4,216,210,242,248,4,217,211,243,249,4,225,224,256,257,4,224,229,261,256,4,229,225,257,261,4,229,230,262,261,4,230,225,257,262,4,226,225,257,258,4,230,226,258,262,4,230,231,263,262,4,231,226,258,263,4,227,226,258,259,4,231,227,259,263,4,231,232,264,263,4,232,227,259,264,4,229,228,260,261,4,228,234,266,260,4,234,229,261,266,4,234,235,267,266,4,235,229,261,267,4,233,232,264,265,4,232,238,270,264,4,238,233,265,270,4,238,239,271,270,4,239,233,265,271,4,234,240,272,266,4,240,235,267,272,4,240,241,273,272,4,241,235,267,273,4,238,244,276,270,4,244,239,271,276,4,244,245,277,276,4,245,239,271,277,4,240,246,278,272,4,246,241,273,278,4,246,247,279,278,4,247,241,273,279,4,244,250,282,276,4,250,245,277,282,4,250,251,283,282,4,251,245,277,283,4,248,247,279,280,4,247,252,284,279,4,252,248,280,284,4,252,253,285,284,4,253,248,280,285,4,249,248,280,281,4,253,249,281,285,4,253,254,286,285,4,254,249,281,286,4,250,249,281,282,4,254,250,282,286,4,254,255,287,286,4,255,250,282,287,4,235,236,268,267,4,236,230,262,268,4,236,237,269,268,4,237,231,263,269,4,237,238,270,269,4,241,242,274,273,4,242,236,268,274,4,242,243,275,274,4,243,237,269,275,4,243,244,276,275,4,248,242,274,280,4,249,243,275,281,4,257,256,288,289,4,256,261,293,288,4,261,257,289,293,4,261,262,294,293,4,262,257,289,294,4,258,257,289,290,4,262,258,290,294,4,262,263,295,294,4,263,258,290,295,4,259,258,290,291,4,263,259,291,295,4,263,264,296,295,4,264,259,291,296,4,261,260,292,293,4,260,266,298,292,4,266,261,293,298,4,266,267,299,298,4,267,261,293,299,4,265,264,296,297,4,264,270,302,296,4,270,265,297,302,4,270,271,303,302,4,271,265,297,303,4,266,272,304,298,4,272,267,299,304,4,272,273,305,304,4,273,267,299,305,4,270,276,308,302,4,276,271,303,308,4,276,277,309,308,4,277,271,303,309,4,272,278,310,304,4,278,273,305,310,4,278,279,311,310,4,279,273,305,311,4,276,282,314,308,4,282,277,309,314,4,282,283,315,314,4,283,277,309,315,4,280,279,311,312,4,279,284,316,311,4,284,280,312,316,4,284,285,317,316,4,285,280,312,317,4,281,280,312,313,4,285,281,313,317,4,285,286,318,317,4,286,281,313,318,4,282,281,313,314,4,286,282,314,318,4,286,287,319,318,4,287,282,314,319,4,267,268,300,299,4,268,262,294,300,4,268,269,301,300,4,269,263,295,301,4,269,270,302,301,4,273,274,306,305,4,274,268,300,306,4,274,275,307,306,4,275,269,301,307,4,275,276,308,307,4,280,274,306,312,4,281,275,307,313])
        self.assertTrue(mm3D.getMeshAtLevel(-1).getNodalConnectivity().isEqual(d))
        d=DataArrayInt(577) ; d.iota(0) ; d*=5
        self.assertTrue(mm3D.getMeshAtLevel(-1).getNodalConnectivityIndex().isEqual(d))
        #
        self.assertEqual(mm3D.getGroupArr(0,"grp0").getName(),"grp0")
        self.assertEqual(mm3D.getGroupArr(-1,"grp1").getName(),"grp1")
        self.assertTrue(mm3D.getGroupArr(0,"grp0").isEqualWithoutConsideringStr(DataArrayInt([0,1,2,3,4,5,24,25,26,27,28,29,48,49,50,51,52,53,72,73,74,75,76,77,96,97,98,99,100,101,120,121,122,123,124,125,144,145,146,147,148,149,168,169,170,171,172,173,192,193,194,195,196,197,216,217,218,225,226,227,234,235,236,243,244,245,252,253,254,261,262,263,270,271,272,279,280,281,288,289,290])))
        self.assertTrue(mm3D.getGroupArr(-1,"grp1").isEqualWithoutConsideringStr(DataArrayInt([0,1,5,9,12,13,14,18,22,23,30,31,33,37,38,40,42,46,50,51,64,65,69,73,76,77,78,82,86,87,94,95,97,101,102,104,106,110,114,115,128,129,133,137,140,141,142,146,150,151,158,159,161,165,166,168,170,174,178,179,192,193,197,201,204,205,206,210,214,215,222,223,225,229,230,232,234,238,242,243,256,257,261,265,268,269,270,274,278,279,286,287,289,293,294,296,298,302,306,307,320,321,325,329,332,333,334,338,342,343,350,351,353,357,358,360,362,366,370,371,384,385,389,393,396,397,398,402,406,407,414,415,417,421,422,424,426,430,434,435,448,449,453,457,460,461,462,466,470,471,478,479,481,485,486,488,490,494,498,499,512,513,517,521,524,525,526,530,534,535,542,543,545,549,550,552,554,558,562,563])))
        pass
    pass

unittest.main()
