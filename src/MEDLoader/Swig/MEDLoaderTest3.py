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
# Author : Anthony Geay (EDF R&D)

from MEDLoader import *
import unittest
import platform
from math import pi,e,sqrt
from MEDLoaderDataForTest import MEDLoaderDataForTest
from distutils.version import LooseVersion

import sys
if sys.version_info.major < 3:
    import cPickle as pickle
else:
    import pickle

class MEDLoaderTest3(unittest.TestCase):
    def testMEDMesh1(self):
        fileName="Pyfile18.med"
        mname="ExampleOfMultiDimW"
        medmesh=MEDFileMesh.New(fileName,mname)
        self.assertRaises(InterpKernelException,MEDFileMesh.New,fileName,"")
        self.assertEqual((0,-1),medmesh.getNonEmptyLevels())
        m1_0=medmesh.getLevel0Mesh(True)
        m1_1=ReadUMeshFromFile(fileName,mname,0)
        self.assertTrue(m1_0.isEqual(m1_1,1e-12));
        m2_0=medmesh.getLevelM1Mesh(True)
        m2_1=ReadUMeshFromFile(fileName,mname,-1)
        self.assertTrue(m2_0.isEqual(m2_1,1e-12));
        pass

    def testMEDMesh2(self):
        fileName="Pyfile10.med"
        mname="3DToto"
        outFileName="MEDFileMesh1.med"
        medmesh=MEDFileUMesh.New(fileName,mname)
        self.assertEqual((0,),medmesh.getNonEmptyLevels())
        m1_0=medmesh.getLevel0Mesh(True)
        m1_1=ReadUMeshFromFile(fileName,mname,0)
        self.assertTrue(m1_0.isEqual(m1_1,1e-12));
        g1_0=medmesh.getGroup(0,"mesh2",True)
        g1_1=ReadUMeshFromGroups(fileName,mname,0,["mesh2"]);
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getGroup(0,"mesh3",True)
        g1_1=ReadUMeshFromGroups(fileName,mname,0,["mesh3"]);
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getGroups(0,["mesh3","mesh2"])
        g1_1=ReadUMeshFromGroups(fileName,mname,0,["mesh3","mesh2"]);
        g1_1.setName(g1_0.getName())
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getFamily(0,"Family_-3",True)
        g1_1=ReadUMeshFromFamilies(fileName,mname,0,["Family_-3"]);
        self.assertTrue(g1_0.isEqual(g1_1,1e-12));
        g1_0=medmesh.getFamilies(0,["Family_-3","Family_-5"],True)
        g1_1=ReadUMeshFromFamilies(fileName,mname,0,["Family_-3","Family_-5"]);
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
        m.checkConsistencyLight()
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(1);
        m1.allocateCells(3);
        m1.insertNextCell(NORM_SEG2,2,[1,4])
        m1.insertNextCell(NORM_SEG2,2,[3,6])
        m1.insertNextCell(NORM_SEG3,3,[2,8,5])
        m1.finishInsertingCells();
        m1.setCoords(c)
        m1.checkConsistencyLight()
        m2=MEDCouplingUMesh.New();
        m2.setMeshDimension(0);
        m2.allocateCells(4);
        m2.insertNextCell(NORM_POINT1,1,[1])
        m2.insertNextCell(NORM_POINT1,1,[3])
        m2.insertNextCell(NORM_POINT1,1,[2])
        m2.insertNextCell(NORM_POINT1,1,[6])
        m2.finishInsertingCells();
        m2.setCoords(c)
        m2.checkConsistencyLight()
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
        g1_N.setValues(list(range(8)),8,1)
        g1_N.setName("G1")
        g2_N=DataArrayInt.New()
        g2_N.setValues(list(range(9)),9,1)
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
        self.assertTrue(t.getValues()==list(range(5)))
        #
        mmCpy=mm.deepCopy()
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
        self.assertEqual(([[(3, 2), (4, 1), (5, 8)], [(1, 2), (2, 1)], [(0, 4)]], 2, 2, 9),GetUMeshGlobalInfo(outFileName,"MyFirstMEDCouplingMEDmesh"))
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
        m.checkConsistencyLight()
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(1);
        m1.allocateCells(3);
        m1.insertNextCell(NORM_SEG2,2,[1,4])
        m1.insertNextCell(NORM_SEG3,3,[2,8,5])
        m1.insertNextCell(NORM_SEG2,2,[3,6])
        m1.finishInsertingCells();
        m1.setCoords(c)
        m1.checkConsistencyLight()
        m2=MEDCouplingUMesh.New();
        m2.setMeshDimension(0);
        m2.allocateCells(4);
        m2.insertNextCell(NORM_POINT1,1,[1])
        m2.insertNextCell(NORM_POINT1,1,[3])
        m2.insertNextCell(NORM_POINT1,1,[2])
        m2.insertNextCell(NORM_POINT1,1,[6])
        m2.finishInsertingCells();
        m2.setCoords(c)
        m2.checkConsistencyLight()
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
        self.assertEqual(list(range(3)),mm2.getGroupArr(-1,"GrpOnAllFace").getValues())
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
        self.assertTrue(m[0].isEqual(m1,1e-12))
        self.assertTrue(isinstance(m[0],MEDCouplingCMesh))
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
        for i in range(nbOfFams):
            m.addFamily(fns[i],fids[i])
            pass
        nbOfGrps=len(grpns)
        for i in range(nbOfGrps):
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
        self.assertEqual([('tyty','mm'),('uiop','MW')],GetComponentsNamesOfField("Pyfile19_bis.med","VFieldOnNodes"))
        pass

    #gauss points
    def testMEDField3(self):
        mm=MEDFileMesh.New("Pyfile13.med")
        mm.write("Pyfile13_bis.med",2)
        ff=MEDFileFieldMultiTS.New("Pyfile13.med","MyFirstFieldOnGaussPoint")
        ff.write("Pyfile13_bis.med",0)
        ff=MEDFileField1TS.New("Pyfile13.med","MyFirstFieldOnGaussPoint",1,5)
        f=ff.getFieldAtLevel(ON_GAUSS_PT,0)
        f2=ReadFieldGauss("Pyfile13.med",'2DMesh_2',0,'MyFirstFieldOnGaussPoint',1,5)
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
        f2=ReadFieldGaussNE("Pyfile14.med",'2DMesh_2',0,"MyFieldOnGaussNE",1,5)
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        pass

    # MEDField get/set on pointe.med
    def testMEDField5(self):
        ff=MEDFileField1TS.New("Pyfile17.med","MeasureOfMesh_Extruded",1,2)
        f=ff.getFieldAtLevel(ON_CELLS,0)
        f2=ReadFieldCell("Pyfile17.med","Extruded",0,"MeasureOfMesh_Extruded",1,2)
        self.assertTrue(f.getMesh().getCoords().isEqual(f2.getMesh().getCoords(),1e-12))
        f.getMesh().tryToShareSameCoords(f2.getMesh(),1e-12)
        f.changeUnderlyingMesh(f2.getMesh(),22,1e-12)
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        # no with renumbering
        f=ff.getFieldAtLevel(ON_CELLS,0,1)
        f2=ReadFieldCell("Pyfile17.med","Extruded",0,"MeasureOfMesh_Extruded",1,2)
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
        f2=ReadFieldNode("Pyfile7.med",'3DSurfMesh_1',0,"VectorFieldOnNodes",its[0][0],its[0][1])
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        ff=MEDFileFieldMultiTS.New("Pyfile19.med","VFieldOnNodes")
        its=ff.getIterations()
        f=ff.getFieldAtLevel(ON_NODES,its[0][0],its[0][1],0)
        f2=ReadFieldNode("Pyfile19.med",'2DMesh_1',0,"VFieldOnNodes",its[0][0],its[0][1])
        self.assertTrue(f.isEqual(f2,1e-12,1e-12))
        self.assertRaises(InterpKernelException,ff.getFieldAtLevel,ON_CELLS,its[0][0],its[0][1],0)# request on cell and it is not on cells
        self.assertRaises(InterpKernelException,ff.getFieldAtLevel,ON_NODES,its[0][0],its[0][1],0,1)#request renumber following mesh : it is on profile !
        pass

    # MEDField get/set on profiles cells
    def testMEDField7(self):
        ff=MEDFileFieldMultiTS.New("Pyfile12.med","VectorFieldOnCells")
        its=ff.getIterations()
        f=ff.getFieldAtLevel(ON_CELLS,its[0][0],its[0][1],0)
        f2=ReadFieldCell("Pyfile12.med",'3DMesh_1',0,"VectorFieldOnCells",its[0][0],its[0][1])
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
        f2=ReadFieldCell(fname,f1.getMesh().getName(),0,f1.getName(),f1.getTime()[1],f1.getTime()[2]);
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
        f2=ReadFieldNode(fname,f1.getMesh().getName(),0,f1.getName(),f1.getTime()[1],f1.getTime()[2])
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
        f2=ReadFieldGaussNE(fname,f1.getMesh().getName(),0,f1.getName(),f1.getTime()[1],f1.getTime()[2])
        self.assertTrue(f1.isEqual(f2,1e-12,1e-12))
        da,infos=ff1.getUndergroundDataArrayExt()
        f2.getArray().setName(da.getName())#da has the same name than f2
        self.assertTrue(da.isEqual(f2.getArray(),1e-12))
        self.assertEqual([((3, 0), (0, 6)), ((4, 0), (6, 14)), ((6, 0), (14, 20))],infos)
        #
        fname="Pyfile28.med"
        f1=MEDLoaderDataForTest.buildVecFieldOnGauss_2_Simpler();
        f1InvalidCpy=f1.deepCopy()
        f1InvalidCpy.setDiscretization(MEDCouplingFieldDiscretizationGauss())
        f1InvalidCpy2=f1.deepCopy()
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
        f22=f21.deepCopy() ; f22.setName("f22") ; f22=f22.buildNewTimeReprFromThis(ONE_TIME,False) ;
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
        ff1=ff1.deepCopy()
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
        ff1=ff1.deepCopy()
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
        m1=m1.buildPartOfMySelf(list(range(5)),True) ; m1.setName(tmp) # suppression of last cell that is a polygon
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
        m1=m1.buildPartOfMySelf(list(range(5)),True) ; m1.setName(tmp) # suppression of last cell that is a polygon
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
        m1=m0.buildPartOfMySelf(list(range(5)),True) ; m1.setName(tmp) ; mm1.setMeshAtLevel(0,m1) ;
        mm1.write(fname,2)
        ff1.write(fname,0)
        f1=ff1.getFieldOnMeshAtLevel(ON_GAUSS_NE,m1,0)
        f2,p1=ff1.getFieldWithProfile(ON_GAUSS_NE,0,mm1) ; f2.setName("")
        self.assertTrue(p1.isIota(5))
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
        ffs=ffs.deepCopy()
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
        f2=f1.deepCopy()
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
        f3=ReadFieldCell(fname,"3DSurfMesh_1",0,"VectorFieldOnCells",0,1)
        self.assertTrue(f3.isEqual(f1,1e-12,1e-12))
        f4=ReadFieldCell(fname,"3DSurfMesh_2",0,"VectorFieldOnCells2",0,1)
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
        f2=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f2.setName(FieldName1) ; f2.setArray(da2) ; f2.setMesh(m2) ; f2.checkConsistencyLight()
        ff1.setFieldNoProfileSBT(f2)
        self.assertEqual(ff1.getNonEmptyLevels(),(2, [0]))
        da0=DataArrayDouble.New()
        da0.alloc(m0.getNumberOfCells()*len(compNames1),1)
        da0.iota(190.)
        da0.rearrange(len(compNames1))
        da0.setInfoOnComponents(compNames1)
        f0=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f0.setName(FieldName1) ; f0.setArray(da0) ; f0.setMesh(m0) ; f0.checkConsistencyLight()
        ff1.setFieldNoProfileSBT(f0)
        self.assertEqual(ff1.getNonEmptyLevels(),(2, [0,-2]))
        da1=DataArrayDouble.New()
        da1.alloc(m1.getNumberOfCells()*len(compNames1),1)
        da1.iota(90.)
        da1.rearrange(len(compNames1))
        da1.setInfoOnComponents(compNames1)
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f1.setName(FieldName1) ; f1.setArray(da1) ; f1.setMesh(m1) ; f1.checkConsistencyLight()
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
        f0=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f0.setName(FieldName2) ; f0.setArray(da0) ; f0.setMesh(m0) ; f0.checkConsistencyLight()
        ff2.setFieldNoProfileSBT(f0)
        self.assertEqual(ff2.getNonEmptyLevels(),(0, [0]))
        da1=DataArrayDouble.New()
        da1.alloc(m1.getNumberOfCells()*len(compNames2),1)
        da1.iota(-90.)
        da1.rearrange(len(compNames2))
        da1.setInfoOnComponents(compNames2)
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME) ; f1.setName(FieldName2) ; f1.setArray(da1) ; f1.setMesh(m1) ; f1.checkConsistencyLight()
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
        for i in range(3):
            self.assertAlmostEqual(nodeCoordsWithValue1[i],tes0.getMesh().getCoordinatesOfNode(0)[i],13);
            self.assertAlmostEqual(nodeCoordsWithValue2[i],tes0.getMesh().getCoordinatesOfNode(1)[i],13);
            pass
        for i in range(6):
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
        for i in range(6):
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
        for i in range(3):
            self.assertAlmostEqual(nodeCoordsWithValue1[i],tes2.getMesh().getCoordinatesOfNode(0)[i],13);
            self.assertAlmostEqual(nodeCoordsWithValue2[i],tes2.getMesh().getCoordinatesOfNode(1)[i],13);
            pass
        for i in range(6):
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
        for i in range(6):
            self.assertAlmostEqual(expected1[i],tes3.getArray().getIJ(0,i),13);
            pass
        pass

    def testBuildInnerBoundaryAlongM1Group1(self):
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
        nodes,cells,cells2=mm.buildInnerBoundaryAlongM1Group("Grp")
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

    def testBuildInnerBoundaryAlongM1Group2(self):
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
        nodes,cells,cells2=mm.buildInnerBoundaryAlongM1Group("Grp")
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

    def testBuildInnerBoundaryAlongM1Group3(self):
        """ Test buildInnerBoundaryAlongM1Group() with *non-connex* cracks """
        fname = "Pyfile73.med"
        m = MEDCouplingCMesh.New()
        m.setCoordsAt(0, DataArrayDouble([0.0,1.1,2.3,3.6,5.0]))
        m.setCoordsAt(1, DataArrayDouble([0.,1.,2.]))
        m = m.buildUnstructured(); m.setName("simple")
        m2 = m.buildDescendingConnectivity()[0]
        m2.setName(m.getName())
            
        # A crack in two non connected parts of the mesh:
        grpSeg = DataArrayInt([3,19]) ; grpSeg.setName("Grp") 

        mm = MEDFileUMesh.New()
        mm.setMeshAtLevel(0,m)
        mm.setMeshAtLevel(-1,m2)
        mm.setGroupsAtLevel(-1,[grpSeg])
        nodes, cellsMod, cellsNotMod = mm.buildInnerBoundaryAlongM1Group("Grp")
        self.assertEqual([1,13],nodes.getValues());
        self.assertEqual([0,6],cellsMod.getValues());
        self.assertEqual([1,7],cellsNotMod.getValues());
        self.assertEqual(17,mm.getNumberOfNodes())
        self.assertEqual([3,19],mm.getGroupArr(-1,"Grp").getValues())
        self.assertEqual([22,23],mm.getGroupArr(-1,"Grp_dup").getValues())
        ref0=[4, 15, 0, 5, 6, 4, 8, 7, 12, 16]
        ref1=[4, 2, 1, 6, 7, 4, 9, 8, 13, 14]
        self.assertEqual(ref0,mm.getMeshAtLevel(0)[[0,6]].getNodalConnectivity().getValues())
        self.assertEqual(ref1,mm.getMeshAtLevel(0)[[1,7]].getNodalConnectivity().getValues())
        self.assertRaises(InterpKernelException,mm.getGroup(-1,"Grp_dup").checkGeoEquivalWith,mm.getGroup(-1,"Grp"),2,1e-12);# Grp_dup and Grp are not equal considering connectivity only
        mm.getGroup(-1,"Grp_dup").checkGeoEquivalWith(mm.getGroup(-1,"Grp"),12,1e-12)# Grp_dup and Grp are equal considering connectivity and coordinates

        refValues=DataArrayDouble([1.1, 1.2, 1.3, 1.4, 1.1, 1.2, 1.3, 1.4])
        valsToTest=mm.getMeshAtLevel(0).getMeasureField(True).getArray() ; delta=(valsToTest-refValues) ; delta.abs()
        self.assertTrue(delta.getMaxValue()[0]<1e-10)
        #
        mm.getCoords()[-len(nodes):]+=[0.,-0.3]
        self.assertRaises(InterpKernelException,mm.getGroup(-1,"Grp_dup").checkGeoEquivalWith,mm.getGroup(-1,"Grp"),12,1e-12);
        refValues2=refValues[:] ; refValues2[0] = 1.265; refValues2[6] = 1.105
        valsToTest=mm.getMeshAtLevel(0).getMeasureField(True).getArray() ;     delta=(valsToTest-refValues2) ; delta.abs()
        self.assertTrue(delta.getMaxValue()[0]<1e-12)
        mm.write(fname,2)   

    def testBuildInnerBoundaryAlongM1Group4(self):
        """ Test case where cells touch the M1 group on some nodes only and not on full egdes (triangle mesh for ex)
        """
        coo = DataArrayDouble([0.,0., 1.,0., 2.,0., 3.,0.,
                               0.,1., 1.,1., 2.,1., 3.,1.,
                               0.,2., 1.,2., 2.,2., 3.,2.], 12, 2)
        conn = [3,0,4,1,  3,1,4,5,
                3,5,9,10, 3,5,10,6,
                3,2,6,7,  3,2,7,3,
                3,4,8,9,  3,4,9,5,
                3,1,5,6,  3,1,6,2,
                3,6,10,11,3,6,11,7]
        # Only TRI3:
        connI = DataArrayInt()
        connI.alloc(13, 1); connI.iota(); connI *= 4
        m2 = MEDCouplingUMesh("2D", 2)
        m2.setCoords(coo)
        m2.setConnectivity(DataArrayInt(conn), connI)
        m2.checkConsistency()
        m1, _, _, _, _ = m2.buildDescendingConnectivity()
        grpIds = DataArrayInt([9,11]); grpIds.setName("group")
        grpIds2 = DataArrayInt([0,1]); grpIds2.setName("group2")
        mfu = MEDFileUMesh()
        mfu.setMeshAtLevel(0, m2)
        mfu.setMeshAtLevel(-1, m1)
        mfu.setGroupsAtLevel(-1, [grpIds, grpIds2])
        nNod = m2.getNumberOfNodes()
        nodesDup, cells1, cells2 = mfu.buildInnerBoundaryAlongM1Group("group")
        m2_bis = mfu.getMeshAtLevel(0)
        m2_bis.checkConsistency()
        m1_bis = mfu.getMeshAtLevel(-1)
        m1_bis.checkConsistency()
        self.assertEqual(nNod+2, mfu.getNumberOfNodes())
        self.assertEqual(nNod+2, m2_bis.getNumberOfNodes())
        self.assertEqual(nNod+2, m1_bis.getNumberOfNodes())
        self.assertEqual([6,7], nodesDup.getValues())
        self.assertEqual([2.,1., 3.,1.], m2_bis.getCoords()[nNod:].getValues())
        self.assertEqual(set([3,10,11]), set(cells1.getValues()))
        self.assertEqual(set([8,9,4,5]), set(cells2.getValues()))
        self.assertEqual([9,11],mfu.getGroupArr(-1,"group").getValues())
        self.assertEqual([23,24],mfu.getGroupArr(-1,"group_dup").getValues())
        self.assertEqual([0,1],mfu.getGroupArr(-1,"group2").getValues())
#         mfu.getMeshAtLevel(0).writeVTK("/tmp/mfu_M0.vtu")
        ref0 =[3, 5, 10, 12, 3, 12, 10, 11, 3, 12, 11, 13]
        ref1 =[3, 2, 6, 7, 3, 2, 7, 3, 3, 1, 5, 6, 3, 1, 6, 2]
        self.assertEqual(ref0,mfu.getMeshAtLevel(0)[[3,10,11]].getNodalConnectivity().getValues())
        self.assertEqual(ref1,mfu.getMeshAtLevel(0)[[4,5,8,9]].getNodalConnectivity().getValues())
        self.assertRaises(InterpKernelException,mfu.getGroup(-1,"group_dup").checkGeoEquivalWith,mfu.getGroup(-1,"group"),2,1e-12) # Grp_dup and Grp are not equal considering connectivity only
        mfu.getGroup(-1,"group_dup").checkGeoEquivalWith(mfu.getGroup(-1,"group"),12,1e-12)# Grp_dup and Grp are equal considering connectivity and coordinates
        m_bis0 = mfu.getMeshAtLevel(-1)
        m_desc, _, _, _, _ = m_bis0.buildDescendingConnectivity()
        m_bis0.checkDeepEquivalOnSameNodesWith(mfu.getMeshAtLevel(-1), 2, 9.9999999)

    def testBuildInnerBoundary5(self):
        """ Full 3D test with tetras only. In this case a tri from the group is not duplicated because it is made only
        of non duplicated nodes. The tri in question is hence not part of the final new "dup" group. """
        coo = DataArrayDouble([200.0, 200.0, 0.0, 200.0, 200.0, 200.0, 200.0, 0.0, 200.0, 200.0, 0.0, 0.0, 0.0, 200.0, 0.0, 0.0, 200.0, 200.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        200.0, 400.0, 200.0, 0.0, 400.0, 200.0, 200.0, 400.0, 0.0, 0.0, 400.0, 0.0, 200.0, 0.0, 100.00000000000016, 200.0, 63.15203310314546, 200.0, 200.0, 134.45205700643342,
         200.0, 200.0, 200.0, 100.00000000000016, 200.0, 63.15203310314546, 0.0, 200.0, 134.45205700643342, 0.0, 200.0, 0.0, 100.00000000000016, 0.0, 63.15203310314546, 
         200.0, 0.0, 134.45205700643342, 200.0, 0.0, 200.0, 100.00000000000016, 0.0, 63.15203310314546, 0.0, 0.0, 134.45205700643342, 0.0, 0.0, 200.0, 200.0, 100.02130053568538, 
         0.0, 200.0, 100.00938163175135, 200.0, 0.0, 100.02130053568538, 0.0, 0.0, 100.00938163175135, 299.3058739933347, 200.0, 200.0, 400.0, 98.68100542924483, 
         200.0, 302.8923433403344, 0.0, 200.0, 302.8923433403344, 200.0, 0.0, 400.0, 100.00000000000016, 0.0, 302.8923433403344, 0.0, 0.0, 400.0, 200.0, 98.55126825835082, 
         400.0, 0.0, 100.02162286181577, 99.31624553977466, 99.99999998882231, 200.0, 99.31624576683302, 100.00000010178034, 0.0, 99.31624560596512, 200.0, 100.0050761312483,
         99.31624560612883, 0.0, 100.00507613125338, 200.0, 99.99999995813045, 100.00950673487786, 0.0, 99.99999989928207, 100.0041870621175, 301.29063354383015, 
         100.0000000093269, 0.0, 301.29063360689975, 0.0, 100.00957769061164, 140.52853868782435, 99.99999963972768, 100.00509135751312, 297.87779091770784, 
         97.16750463405486, 97.18018457127863], 46, 3)
        c0 = [14, 45, 31, 21, 42, 14, 37, 38, 20, 44, 14, 39, 36, 41, 44, 14, 5, 25, 12, 13, 14, 38, 36, 44, 41, 14, 21, 20, 24, 44, 14, 38, 25, 41, 19, 14, 37, 38, 44, 41, 14, 16, 27,
         39, 41, 14, 21, 45, 26, 40, 14, 39, 37, 44, 41, 14, 14, 15, 24, 44, 14, 25, 38, 41, 13, 14, 27, 18, 6, 22, 14, 38, 36, 41, 13, 14, 44, 14, 15, 36, 14, 44, 23, 39, 26, 14,
         21,26, 23, 44, 14, 38, 44, 14, 24, 14, 39, 37, 41, 22, 14, 21, 33, 45, 42, 14, 27, 22, 39, 41, 14, 23, 26, 21, 3, 14, 27, 18, 22, 41, 14, 39, 36, 44, 17, 14, 21, 26, 44, 40,
         14, 39, 37, 22, 23, 14, 37, 38, 41, 19, 14, 25, 12, 13, 41, 14, 30, 26, 43, 45, 14, 38, 36, 13, 14, 14, 12, 36, 13, 41, 14, 20, 44, 21, 37, 14, 16, 36, 12, 41, 14, 39, 36,
         17, 16, 14, 44, 20, 24, 38, 14, 27, 16, 12, 41, 14, 26, 15, 17, 44, 14, 19, 18, 41, 37, 14, 40, 45, 26, 15, 14, 37, 38, 19, 20, 14, 17, 15, 26, 2, 14, 39, 36, 16, 41, 14,
         24, 21, 44, 40, 14, 16, 7, 27, 12, 14, 22, 18, 37, 41, 14, 21, 31, 45, 24, 14, 44, 40, 15, 24, 14, 24, 45, 15, 28, 14, 44, 40, 26, 15, 14, 24, 20, 21, 0, 14, 38, 36, 14,
         44, 14, 39, 37, 23, 44, 14, 45, 31, 42, 32, 14, 25, 18, 19, 4, 14, 36, 44, 17, 15, 14, 25, 19, 18, 41, 14, 24, 15, 14, 1, 14, 45, 24, 34, 28, 14, 35, 45, 30, 43, 14, 17,
         44, 39, 26, 14, 44, 23, 21, 37, 14, 30, 45, 29, 15, 14, 45, 35, 33, 43, 14, 30, 15, 26, 45, 14, 31, 21, 0, 24, 14, 33, 35, 32, 10, 14, 29, 45, 34, 28, 14, 32, 45, 34,
         29, 14, 45, 31, 32, 34, 14, 33, 26, 45, 43, 14, 45, 31, 34, 24, 14, 33, 26, 21, 45, 14, 11, 30, 35, 29, 14, 33, 35, 45, 32, 14, 33, 45, 42, 32, 14, 32, 8, 34, 31, 14,
         21, 26, 33, 3, 14, 35, 45, 32, 29, 14, 29, 34, 9, 28, 14, 15, 45, 24, 40, 14, 29, 45, 28, 15, 14, 21, 24, 45, 40, 14, 24, 15, 1, 28, 14, 35, 45, 29, 30, 14, 26, 15,
         30, 2]
        cI0 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185,
         190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 
         360, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410, 415, 420, 425, 430]
        m3 = MEDCouplingUMesh("3D", 3)
        m3.setCoords(coo)
        m3.setConnectivity(DataArrayInt(c0), DataArrayInt(cI0))
        m3.checkConsistency()
        m2, _, _, _, _ = m3.buildDescendingConnectivity()
        grpIds = DataArrayInt([36,74]); grpIds.setName("group")
        mfu = MEDFileUMesh()
        mfu.setMeshAtLevel(0, m3)
        mfu.setMeshAtLevel(-1, m2)
        grpIds3D = DataArrayInt([0,1]); grpIds3D.setName("group_3d")
        mfu.setGroupsAtLevel(0, [grpIds3D])  # just to check preservation of 3D group
        mfu.setGroupsAtLevel(-1, [grpIds])
        nNod = m3.getNumberOfNodes()
        nodesDup, cells1, cells2 = mfu.buildInnerBoundaryAlongM1Group("group")
        m3_bis = mfu.getMeshAtLevel(0)
        m3_bis.checkConsistency()
        m2_bis = mfu.getMeshAtLevel(-1)
        m2_bis.checkConsistency()
        self.assertEqual(nNod+1, mfu.getNumberOfNodes())
        self.assertEqual(nNod+1, m3_bis.getNumberOfNodes())
        self.assertEqual(nNod+1, m2_bis.getNumberOfNodes())
        self.assertEqual([3], nodesDup.getValues())
        self.assertEqual(m3_bis.getCoords()[3].getValues(), m3_bis.getCoords()[nNod:].getValues())
        self.assertEqual(set([22]), set(cells1.getValues()))
        self.assertEqual(set([77]), set(cells2.getValues()))
        self.assertEqual([36,74],mfu.getGroupArr(-1,"group").getValues())
        self.assertEqual([0,1],mfu.getGroupArr(0,"group_3d").getValues())
        self.assertEqual([213],mfu.getGroupArr(-1,"group_dup").getValues())  # here only one cell has been duplicated
        m_bis0 = mfu.getMeshAtLevel(-1)
        m_desc, _, _, _, _ = m_bis0.buildDescendingConnectivity()
        m_bis0.checkDeepEquivalOnSameNodesWith(mfu.getMeshAtLevel(-1), 2, 9.9999999)
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
        f=m.getMeasureField(False)
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
        fread.checkConsistencyLight()
        fread2.checkConsistencyLight()
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
        mQ8.insertNextCell(NORM_QUAD8,list(range(8)))
        mQ8.finishInsertingCells()
        mQ4=MEDCouplingUMesh("",2) ; mQ4.setCoords(coords)
        mQ4.allocateCells(1)
        mQ4.insertNextCell(NORM_QUAD4,list(range(4)))
        mQ4.finishInsertingCells()
        mT3=MEDCouplingUMesh("",2) ; mT3.setCoords(coords)
        mT3.allocateCells(1)
        mT3.insertNextCell(NORM_TRI3,list(range(3)))
        mT3.finishInsertingCells()
        
        tr=[[0.,4.],[2.,4.],[4.,4.],[6.,4.],[8.,4.],[10.,4.],[12.,4.],[14.,4.],[16.,4.],[18.,4.],[20.,4.],[0.,0.],[2.,0.], [0.,2.],[2.,2.],[4.,2.],[6.,2.],[8.,2.],[10.,2.],[12.,2.]]
        ms=11*[mT3]+2*[mQ4]+7*[mQ8]
        ms[:]=(elt.deepCopy() for elt in ms)
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
        fInvalid=f.deepCopy()
        f.setGaussLocalizationOnCells([0,1,2,3,4,5,6,7,8],[0.,0.,1.,0.,1.,1.],[0.3,0.3,0.7,0.7],[0.8,0.2])
        f.setGaussLocalizationOnCells([9,10],[0.,0.,1.,0.,1.,1.],[0.3,0.3,0.7,0.7,0.8,0.8],[0.8,0.07,0.13])
        f.setGaussLocalizationOnCells([11,12],[0.,0.,1.,0.,1.,1.,0.,1.],[0.3,0.3,0.7,0.7,0.8,0.8,0.8,0.8,0.8,0.8],[0.8,0.07,0.1,0.01,0.02])
        f.checkConsistencyLight()
        fInvalid2=fInvalid.deepCopy()
        fInvalid2.getDiscretization().setArrayOfDiscIds(f.getDiscretization().getArrayOfDiscIds())
        #
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        mm.write(fname,2)
        #
        f1ts=MEDFileField1TS.New()
        pfl=DataArrayInt(list(range(13))) ; pfl.setName("pfl")
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
        f2_bis=ReadFieldGauss(fname,m.getName(),0,f.getName(),f.getTime()[1],f.getTime()[2])
        f2_bis.checkConsistencyLight()
        self.assertTrue(f.isEqual(f2_bis,1e-12,1e-12))
        #
        WriteField(fname2,f,True)
        f2_ter=ReadFieldGauss(fname2,m.getName(),0,f.getName(),f.getTime()[1],f.getTime()[2])
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
        f.checkConsistencyLight()
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
        f3_bis=ReadFieldGauss(fname,m.getName(),0,f.getName(),f.getTime()[1],f.getTime()[2])
        f3_bis.renumberCells([0,1,3,2,4,5,6,7,8,9])
        self.assertTrue(f.isEqual(f3_bis,1e-12,1e-12))
        #
        WriteField(fname2,f,True)
        f3_ter=ReadFieldGauss(fname2,m.getName(),0,f.getName(),f.getTime()[1],f.getTime()[2])
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
        f.checkConsistencyLight()
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
        f3_bis=ReadFieldGauss(fname,m.getName(),0,f.getName(),f.getTime()[1],f.getTime()[2])
        f3_bis.renumberCells([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,16,19])
        self.assertTrue(f.isEqual(f3_bis,1e-12,1e-12))
        #
        WriteField(fname2,f,True)
        f3_ter=ReadFieldGauss(fname2,m.getName(),0,f.getName(),f.getTime()[1],f.getTime()[2])
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
        for i in range(9):
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
        m00=m0.deepCopy() ; m00=m00[[0,2]] ; m00.setName(m.getName()) ; m00.zipCoords()
        fieldNode0.setMesh(m00)
        f0=MEDFileField1TS.New(fname,fieldNode0.getName(),dt,it)
        ff0_1=f0.getFieldOnMeshAtLevel(ON_NODES,m0)
        ff0_1.checkConsistencyLight()
        self.assertTrue(ff0_1.isEqual(fieldNode0,1e-12,1e-12))
        ff0_2=f0.getFieldAtLevel(ON_NODES,0)
        ff0_2.checkConsistencyLight()
        self.assertTrue(ff0_2.isEqual(fieldNode0,1e-12,1e-12))
        ff0_3=f0.getFieldOnMeshAtLevel(ON_NODES,0,m)
        ff0_3.checkConsistencyLight()
        self.assertTrue(ff0_3.isEqual(fieldNode0,1e-12,1e-12))
        ff0_4=ReadFieldNode(fname,m.getName(),0,fieldNode0.getName(),dt,it)
        ff0_4.checkConsistencyLight()
        self.assertTrue(ff0_4.isEqual(fieldNode0,1e-12,1e-12))
        f1=MEDFileField1TS.New(fname,fieldNode1.getName(),dt,it)
        m1=m.getMeshAtLevel(-1)
        m10=m1.deepCopy() ; m10=m10[[0,1,2,3,4,5,6,7]] ; m10.setName(m.getName()) ; m10.zipCoords()
        fieldNode1.setMesh(m10)
        ff1_1=f1.getFieldOnMeshAtLevel(ON_NODES,m1)
        ff1_1.checkConsistencyLight()
        self.assertTrue(ff1_1.isEqual(fieldNode1,1e-12,1e-12))
        ff1_2=f1.getFieldAtLevel(ON_NODES,-1)
        ff1_2.checkConsistencyLight()
        self.assertTrue(ff1_2.isEqual(fieldNode1,1e-12,1e-12))
        ff1_3=f1.getFieldOnMeshAtLevel(ON_NODES,-1,m)
        ff1_3.checkConsistencyLight()
        self.assertTrue(ff1_3.isEqual(fieldNode1,1e-12,1e-12))
        ff1_4=ReadFieldNode(fname,m.getName(),-1,fieldNode1.getName(),dt,it)
        ff1_4.checkConsistencyLight()
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
        for i in range(9):
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
        m00=m0.deepCopy() ; m00=m00[pfl0] ; m00.setName(m.getName())
        fieldCell0.setMesh(m00)
        f0=MEDFileField1TS.New(fname,fieldCell0.getName(),dt,it)
        ff0_1=f0.getFieldOnMeshAtLevel(ON_CELLS,m0)
        ff0_1.checkConsistencyLight()
        self.assertTrue(ff0_1.isEqual(fieldCell0,1e-12,1e-12))
        ff0_2=f0.getFieldAtLevel(ON_CELLS,0)
        ff0_2.checkConsistencyLight()
        self.assertTrue(ff0_2.isEqual(fieldCell0,1e-12,1e-12))
        ff0_3=f0.getFieldOnMeshAtLevel(ON_CELLS,0,m)
        ff0_3.checkConsistencyLight()
        self.assertTrue(ff0_3.isEqual(fieldCell0,1e-12,1e-12))
        ff0_4=ReadFieldCell(fname,m.getName(),0,fieldCell0.getName(),dt,it)
        ff0_4.checkConsistencyLight()
        self.assertTrue(ff0_4.isEqual(fieldCell0,1e-12,1e-12))
        f1=MEDFileField1TS.New(fname,fieldCell1.getName(),dt,it)
        m1=m.getMeshAtLevel(-1)
        m10=m1.deepCopy() ; m10=m10[pfl1] ; m10.setName(m.getName())
        fieldCell1.setMesh(m10)
        ff1_1=f1.getFieldOnMeshAtLevel(ON_CELLS,m1)
        ff1_1.checkConsistencyLight()
        self.assertTrue(ff1_1.isEqual(fieldCell1,1e-12,1e-12))
        ff1_2=f1.getFieldAtLevel(ON_CELLS,-1)
        ff1_2.checkConsistencyLight()
        self.assertTrue(ff1_2.isEqual(fieldCell1,1e-12,1e-12))
        ff1_3=f1.getFieldOnMeshAtLevel(ON_CELLS,-1,m)
        ff1_3.checkConsistencyLight()
        self.assertTrue(ff1_3.isEqual(fieldCell1,1e-12,1e-12))
        ff1_4=ReadFieldCell(fname,m.getName(),-1,fieldCell1.getName(),dt,it)
        ff1_4.checkConsistencyLight()
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
        mm=m.deepCopy()
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
        for i in range(7):
            m0.insertNextCell(NORM_TRI3,[1,2,1])
            pass
        for i in range(4):
            m0.insertNextCell(NORM_QUAD4,[1,1,2,0])
            pass
        for i in range(2):
            m0.insertNextCell(NORM_POLYGON,[0,0,1,1,2,2])
            pass
        m1=MEDCouplingUMesh("toto",1) ; m1.allocateCells(0) ; m1.insertNextCell(NORM_SEG2,[1,6]) ; m1.insertNextCell(NORM_SEG2,[7,3])
        m2=MEDCouplingUMesh("toto",0) ; m2.allocateCells(0) ; m2.insertNextCell(NORM_POINT1,[2]) ; m2.insertNextCell(NORM_POINT1,[6]) ; m2.insertNextCell(NORM_POINT1,[8])
        m0.setCoords(coo) ; m.setMeshAtLevel(0,m0)
        m1.setCoords(coo) ; m.setMeshAtLevel(-1,m1)
        m2.setCoords(coo) ; m.setMeshAtLevel(-2,m2)
        #
        mm=m.deepCopy()
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
        f=m.getMeasureField(False)
        self.assertIn(m.getHeapMemorySize(), list(range(3552 - 100, 3552 + 100 + 4 * strMulFac)))
        self.assertIn(f.getHeapMemorySize(), list(range(4215 - 100, 4215 + 100 + 8 * strMulFac)))
        #
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        self.assertIn(mm.getHeapMemorySize(), list(range(3889 - 100, 4225 + 100 + 10 * strMulFac)))
        ff=MEDFileField1TS()
        ff.setFieldNoProfileSBT(f)
        self.assertIn(ff.getHeapMemorySize(), list(range(771 - 40, 871 + 21 + (4 + 1) * strMulFac)))
        #
        fff=MEDFileFieldMultiTS()
        fff.appendFieldNoProfileSBT(f)
        self.assertIn(fff.getHeapMemorySize(), list(range(815 - 50, 915 + 30 + (6 + 2) * strMulFac)))
        f.setTime(1.,0,-1)
        fff.appendFieldNoProfileSBT(f)
        self.assertIn(fff.getHeapMemorySize(), list(range(1594 - 90, 1794 + 50 + (10 + 1) * strMulFac)))
        self.assertIn(fff[0, -1].getHeapMemorySize(), list(range(771 - 40, 871 + 20 + (4 + 1) * strMulFac)))
        f2=f[:50]
        f2.setTime(2.,1,-1)
        pfl=DataArrayInt.Range(0,50,1) ; pfl.setName("pfl")
        fff.appendFieldProfile(f2,mm,0,pfl)
        self.assertIn(fff.getHeapMemorySize(), list(range(2348 - 130, 2608 + 100 + (10 + 2) * strMulFac)))
        self.assertIn(fff.getProfile("pfl").getHeapMemorySize(), list(range(204 - 10, 204 + 10 + 2 * strMulFac)))
        self.assertIn(fff[1, -1].getHeapMemorySize(), list(range(738 - 50, 838 + 30 + 4 * strMulFac)))
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
        mesh.checkConsistencyLight();
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
        pts2=pts.deepCopy() ; pts2.setName("B") ; pts2.setDescription("A second example")
        p.pushParam(pts) ; p.pushParam(pts2)
        data.write(fname,2)
        p2=MEDFileParameters(fname)
        self.assertTrue(p.isEqual(p2,1e-14)[0])
        self.assertAlmostEqual(p[1][1,2].getValue(),567.89,13)
        p3=p.deepCopy()
        pts4=pts2.deepCopy()
        pts3=pts2.deepCopy()
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
        namesCellL0[:] = ["CellL0#%.3d      " % (i) for i in range(6)]
        mm.setNameFieldAtLevel(0,namesCellL0)
        namesCellL1=DataArrayAsciiChar.Aggregate([namesCellL0,namesCellL0,namesCellL0.subArray(2)])
        namesCellL1[:] = ["CellLM1#%.3d     " % (i) for i in range(16)]
        mm.setNameFieldAtLevel(-1,namesCellL1)
        namesNodes=namesCellL1.subArray(4,16)
        namesNodes[:] = ["Node#%.3d        " % (i) for i in range(12)]
        mm.setNameFieldAtLevel(1,namesNodes)
        mm.write(fname,2)
        #
        mmr=MEDFileMesh.New(fname)
        self.assertTrue(mm.getNameFieldAtLevel(0).isEqual(DataArrayAsciiChar(["CellL0#%.3d      " % (i) for i in range(6)])))
        self.assertTrue(mm.getNameFieldAtLevel(-1).isEqual(DataArrayAsciiChar(["CellLM1#%.3d     " % (i) for i in range(16)])))
        self.assertTrue(mm.getNameFieldAtLevel(1).isEqual(DataArrayAsciiChar(["Node#%.3d        " % (i) for i in range(12)])))
        self.assertTrue(mm.isEqual(mmr,1e-12)[0])
        mmr.getNameFieldAtLevel(1).setIJ(0,0,'M')
        self.assertTrue(not mm.isEqual(mmr,1e-12)[0])
        mmr.getNameFieldAtLevel(1).setIJ(0,0,'N')
        self.assertTrue(mm.isEqual(mmr,1e-12)[0])
        mmCpy=mm.deepCopy()
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
        self.assertTrue(mmr.getNameFieldAtLevel(0).isEqual(DataArrayAsciiChar(["CellL0#%.3d      " % (i) for i in range(6)])))
        self.assertEqual(mmr.getNameFieldAtLevel(-1),None)
        #
        c=MEDCouplingCMesh()
        arr=DataArrayDouble([0.,1.1,2.3])
        c.setCoords(arr,arr)
        c.setName("cmesh")
        cc=MEDFileCMesh()
        cc.setMesh(c)
        cc.setNameFieldAtLevel(0, DataArrayAsciiChar(["Cell#%.3d        " % (i) for i in range(4)]))
        cc.setNameFieldAtLevel(1, DataArrayAsciiChar(["Node#%.3d        " % (i) for i in range(9)]))
        cc.write(fname2,2)
        ccr=MEDFileMesh.New(fname2)
        self.assertTrue(ccr.getNameFieldAtLevel(0).isEqual(DataArrayAsciiChar(["Cell#%.3d        " % (i) for i in range(4)])))
        self.assertTrue(ccr.getNameFieldAtLevel(1).isEqual(DataArrayAsciiChar(["Node#%.3d        " % (i) for i in range(9)])))
        self.assertTrue(cc.isEqual(ccr,1e-12)[0])
        ccr.getNameFieldAtLevel(1).setIJ(0,0,'M')
        self.assertTrue(not cc.isEqual(ccr,1e-12)[0])
        ccr.getNameFieldAtLevel(1).setIJ(0,0,'N')
        self.assertTrue(cc.isEqual(ccr,1e-12)[0])
        ccCpy=cc.deepCopy()
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
            tmp=c2.findIdsEqual(grpId)
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
        for i in range(nbOfFams):
            m.addFamily(fns[i],fids[i])
            pass
        nbOfGrps=len(grpns)
        for i in range(nbOfGrps):
            m.setFamiliesIdsOnGroup(grpns[i],famIdsPerGrp[i])
            pass
        m.setName(m2.getName())
        m.setDescription(m2.getDescription())
        m.write(fileName,2)
        #
        mm0=MEDFileMesh.New(fileName)
        mm1=MEDFileMesh.New(fileName)
        groupNamesIni=GetMeshGroupsNames(fileName,"ma")
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
        f1=f1.convertToIntField()
        m1=f1.getMesh()
        mm1=MEDFileUMesh.New()
        mm1.setCoords(m1.getCoords())
        mm1.setMeshAtLevel(0,m1)
        mm1.setName(m1.getName())
        mm1.write(fname,2)
        ff1=MEDFileIntField1TS()
        ff1.setFieldNoProfileSBT(f1)
        a=ff1.getFieldOnMeshAtLevel(0,ON_CELLS,mm1)
        self.assertEqual(a.getArray().getInfoOnComponents(),['power [MW/m^3]','density [g/cm^3]','temperature [K]'])
        self.assertTrue(a.isEqual(f1,1e-12,0))
        ff1.write(fname,0)
        ff2=MEDFileAnyTypeField1TS.New(fname)
        self.assertEqual(ff2.getName(),"VectorFieldOnCells")
        self.assertEqual(ff2.getTime(),[0,1,2.0])
        self.assertTrue(isinstance(ff2,MEDFileIntField1TS))
        a=ff1.getFieldOnMeshAtLevel(ON_CELLS,0,mm1)
        self.assertEqual(a.getArray().getInfoOnComponents(),['power [MW/m^3]','density [g/cm^3]','temperature [K]'])
        self.assertTrue(a.isEqual(f1,1e-12,0))
        ff2.setTime(1,2,3.)
        c=ff2.getUndergroundDataArray() ; c*=2
        ff2.write(fname,0) # 2 time steps in 
        ffs1=MEDFileAnyTypeFieldMultiTS.New(fname,"VectorFieldOnCells")
        self.assertEqual(ffs1.getTimeSteps(),[(0, 1, 2.0), (1, 2, 3.0)])
        self.assertEqual(len(ffs1),2)
        self.assertTrue(isinstance(ffs1,MEDFileIntFieldMultiTS))
        a=ffs1[2.].getFieldOnMeshAtLevel(ON_CELLS,0,mm1)
        self.assertTrue(a.isEqual(f1,1e-12,0))
        a=ffs1.getFieldOnMeshAtLevel(ON_CELLS,0,1,0,mm1)
        self.assertTrue(a.isEqual(f1,1e-12,0))
        it=ffs1.__iter__() ; it.next() ; ff2bis=it.next()
        a=ff2bis.getFieldOnMeshAtLevel(0,ON_CELLS,mm1)
        self.assertTrue(a.getArray().isEqual(2*f1.getArray()))
        f1.setTime(3.,1,2) ; f1.getArray()[:]*=2
        self.assertTrue(a.isEqual(f1,1e-12,0)) ; f1.getArray()[:]/=2
        bc=DataArrayInt(6,3) ; bc[:]=0 ; bc.setInfoOnComponents(['power [MW/m^3]','density [g/cm^3]','temperature [K]'])
        for it in ffs1:
            a=it.getFieldOnMeshAtLevel(ON_CELLS,0,mm1)
            bc+=a.getArray()
            pass
        self.assertTrue(bc.isEqual(3*f1.getArray()))
        nf1=MEDCouplingFieldInt(ON_NODES)
        nf1.setTime(9.,10,-1)
        nf1.setMesh(f1.getMesh())
        narr=DataArrayInt(12,2) ; narr.setInfoOnComponents(["aa [u1]","bbbvv [ppp]"]) ; narr[:,0]=list(range(12)) ; narr[:,1]=2*narr[:,0]
        nf1.setName("VectorFieldOnNodes") ; nf1.setArray(narr)
        nff1=MEDFileIntField1TS.New()
        nff1.setFieldNoProfileSBT(nf1)
        self.assertEqual(nff1.getInfo(),('aa [u1]','bbbvv [ppp]'))
        self.assertEqual(nff1.getTime(),[10,-1,9.0])
        nff1.write(fname,0)
        #
        nf2=MEDCouplingFieldInt(ON_NODES)
        nf2.setTime(19.,20,-11)
        nf2.setMesh(f1.getMesh())
        narr2=DataArrayInt(8,2) ; narr.setInfoOnComponents(["aapfl [u1]","bbbvvpfl [ppp]"]) ; narr2[:,0]=list(range(8)) ; narr2[:,0]+=10  ; narr2[:,1]=3*narr2[:,0]
        nf2.setName("VectorFieldOnNodesPfl") ; narr2.setName(nf2.getName()) ; nf2.setArray(narr2)
        nff2=MEDFileIntField1TS.New()
        npfl=DataArrayInt([1,2,4,5,6,7,10,11]) ; npfl.setName("npfl")
        nff2.setFieldProfile(nf2,mm1,0,npfl)
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
        self.assertTrue(fs["VectorFieldOnCells"][0].getUndergroundDataArray().isEqualWithoutConsideringStr(f1.getArray()))
        self.assertTrue(fs["VectorFieldOnCells"][1,2].getUndergroundDataArray().isEqualWithoutConsideringStr(2*f1.getArray()))
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
        arr=DataArrayDouble(12,2) ; arr.setInfoOnComponents(["aa [u1]","bbbvv [ppp]"]) ; arr[:,0]=list(range(12)) ; arr[:,1]=2*arr[:,0]
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
        tris = [tri.deepCopy() for i in range(4)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads = [quad.deepCopy() for i in range(5)]
        for i,elt in enumerate(quads): elt.translate([5+i,0])
        quads=MEDCouplingUMesh.MergeUMeshes(quads)
        m=MEDCouplingUMesh.MergeUMeshes(tris,quads)
        m.setName("mesh") ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        #
        fmts0_0=MEDFileFieldMultiTS()
        fmts0_1=MEDFileFieldMultiTS()
        # time steps
        for i in range(10):
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
        fmts0_2=fmts0_0.deepCopy()
        fmts0_3=fmts0_0.deepCopy()
        fmts0_4=fmts0_0.deepCopy()
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
        tris = [tri.deepCopy() for i in range(4)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads = [quad.deepCopy() for i in range(5)]
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
        for i in range(10):
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
        fmts0_2=fmts0_0.deepCopy()
        fmts0_3=fmts0_0.deepCopy()
        fmts0_4=fmts0_0.deepCopy()
        fs0=MEDFileFields()
        fs0.pushField(fmts0_0)
        fmts0_2.setName("2ndField") ; fs0.pushField(fmts0_2)
        fmts0_3.setName("3rdField") ; fs0.pushField(fmts0_3)
        fmts0_4.setName("4thField") ; fs0.pushField(fmts0_4)
        self.assertEqual(fs0.getPfls(),('pfl_NORM_QUAD4',))
        #
        fmts0_5=MEDFileFieldMultiTS()
        for i in range(7):
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
        tris = [tri.deepCopy() for i in range(4)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads = [quad.deepCopy() for i in range(5)]
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
        for i in range(10):
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
        for i in range(10):
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
        tris = [tri.deepCopy() for i in range(4)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads = [quad.deepCopy() for i in range(5)]
        for i,elt in enumerate(quads): elt.translate([5+i,0])
        quads=MEDCouplingUMesh.MergeUMeshes(quads)
        m=MEDCouplingUMesh.MergeUMeshes(tris,quads)
        m.setName("mesh") ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        mm=MEDFileUMesh() ; mm.setMeshAtLevel(0,m)
        #
        ff0=MEDFileField1TS()
        f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m) ; arr=DataArrayDouble(m.getNumberOfCells()*2) ; arr.iota() ; arr.rearrange(2) ; arr.setInfoOnComponents(["X [km]","YY [mm]"]) ; f0.setArray(arr) ; f0.setName("FieldCell")
        f0.checkConsistencyLight()
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
        f0.checkConsistencyLight()
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
        f0.checkConsistencyLight()
        ff0.appendFieldProfile(f0,mm,0,pfl)
        f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m[:7]) ; arr=DataArrayDouble(7*2) ; arr.iota(100) ; arr.rearrange(2) ; arr.setInfoOnComponents(["X [km]","YY [mm]"]) ; f0.setArray(arr) ; f0.setName("FieldCellMTime") ; f0.setTime(1.1,1,11)
        f0.checkConsistencyLight()
        ff0.appendFieldProfile(f0,mm,0,pfl)
        f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m[:7]) ; arr=DataArrayDouble(7*2) ; arr.iota(200) ; arr.rearrange(2) ; arr.setInfoOnComponents(["X [km]","YY [mm]"]) ; f0.setArray(arr) ; f0.setName("FieldCellMTime") ; f0.setTime(2.1,2,12)
        f0.checkConsistencyLight()
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
        tris = [tri.deepCopy() for i in range(30)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads = [quad.deepCopy() for i in range(40)]
        for i,elt in enumerate(quads): elt.translate([40+i,0])
        quads=MEDCouplingUMesh.MergeUMeshes(quads)
        m=MEDCouplingUMesh.MergeUMeshes(tris,quads)
        m.setName("mesh") ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        mm=MEDFileUMesh() ; mm.setMeshAtLevel(0,m) ; mm.write(fname,2)
        #
        ff0=MEDFileField1TS()
        f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m) ; arr=DataArrayDouble(m.getNumberOfCells()*2) ; arr.iota() ; arr.rearrange(2) ; arr.setInfoOnComponents(["X [km]","YY [mm]"]) ; f0.setArray(arr) ; f0.setName("FieldCell")
        f0.checkConsistencyLight()
        ff0.setFieldNoProfileSBT(f0)
        ff0.write(fname,0)
        #
        fspExp=[(3,[(0,(0,30),'','')]),(4,[(0,(30,70),'','')])]
        self.assertEqual(ff0.getFieldSplitedByType(),fspExp)
        # With profiles
        ff0=MEDFileField1TS()
        f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m[:50]) ; arr=DataArrayDouble(50*2) ; arr.iota() ; arr.rearrange(2) ; arr.setInfoOnComponents(["XX [pm]","YYY [hm]"]) ; f0.setArray(arr) ; f0.setName("FieldCellPfl")
        f0.checkConsistencyLight()
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
        self.assertIn(heap_memory_ref, list(range(182, 540 + 2 * strMulFac)))
        ff0.loadArrays() ##
        arr=DataArrayDouble(140) ; arr.iota() ; arr.rearrange(2)
        self.assertTrue(ff0.getUndergroundDataArray().isEqualWithoutConsideringStr(arr,1e-14))
        self.assertEqual(ff0.getHeapMemorySize()-heap_memory_ref,70*8*2)
        #
        ff0=MEDFileField1TS(fname,"FieldCellPfl",False)
        self.assertEqual(ff0.getUndergroundDataArray().getInfoOnComponents(),["XX [pm]","YYY [hm]"])
        heap_memory_ref=ff0.getHeapMemorySize()
        self.assertIn(heap_memory_ref, list(range(350, 600 + 6 * strMulFac)))
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
        self.assertIn(heap_memory_ref, list(range(1100, 1400 + 2 * strMulFac)))
        ff0.unloadArrays()
        hmd=ff0.getHeapMemorySize()-heap_memory_ref
        self.assertEqual(hmd,-800) # -50*8*2
        ff0.loadArrays() ##
        self.assertEqual(ff0.getHeapMemorySize()-heap_memory_ref,0)
        #
        ff0=MEDFileField1TS(fname,"FieldCellPfl",-1,-1,False)
        heap_memory_ref=ff0.getHeapMemorySize()
        self.assertIn(heap_memory_ref, list(range(299, 620 + 6 * strMulFac)))
        ff0.loadArrays() ##
        self.assertTrue(ff0.getUndergroundDataArray().isEqualWithoutConsideringStr(arr,1e-14))
        self.assertEqual(ff0.getHeapMemorySize()-heap_memory_ref,50*8*2)
        #
        fieldName="FieldCellMultiTS"
        ff0=MEDFileFieldMultiTS()
        for t in range(20):
            f0=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f0.setMesh(m) ; arr=DataArrayDouble(m.getNumberOfCells()*2) ; arr.iota(float(t+1000)) ; arr.rearrange(2) ; arr.setInfoOnComponents(["X [km]","YY [mm]"]) ; f0.setArray(arr) ; f0.setName(fieldName)
            f0.setTime(float(t)+0.1,t,100+t)
            f0.checkConsistencyLight()
            ff0.appendFieldNoProfileSBT(f0)
            pass
        ff0.write(fname,0)
        #
        ff0=MEDFileAnyTypeFieldMultiTS.New(fname,fieldName,False)
        heap_memory_ref=ff0.getHeapMemorySize()
        self.assertIn(heap_memory_ref, list(range(5536, 9212 + (80 + 26 + 1) * strMulFac)))
        ff0.loadArrays()
        self.assertEqual(ff0.getHeapMemorySize()-heap_memory_ref,20*70*8*2)
        del ff0
        #
        ffs=MEDFileFields(fname,False)
        heap_memory_ref=ffs.getHeapMemorySize()
        self.assertIn(heap_memory_ref, list(range(5335, 10031 + (80 + 50 + len(ffs)) * strMulFac)))
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
        tris = [tri.deepCopy() for i in range(4)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads = [quad.deepCopy() for i in range(5)]
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
        tris = [tri.deepCopy() for i in range(4)]
        for i,elt in enumerate(tris): elt.translate([i,0])
        tris=MEDCouplingUMesh.MergeUMeshes(tris)
        quad=MEDCouplingUMesh("quad",2)
        quad.allocateCells() ; quad.insertNextCell(NORM_QUAD4,[0,1,2,3])
        quad.setCoords(DataArrayDouble([(0.,0.),(0.,1.),(1.,1.),(1.,0.)]))
        quads = [quad.deepCopy() for i in range(5)]
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
        m.setRenumFieldArr(1,DataArrayInt(list(range(10,26))))
        m.setFamilyFieldArr(1,DataArrayInt([-1,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,0,-1,-3,-3,-3]))
        m.write(fname,2)
        del m,a,c,m0,m00
        #
        m=MEDFileMesh.New(fname)
        self.assertEqual((),m.getNonEmptyLevels())
        self.assertTrue(m.getCoords().isEqual(DataArrayDouble([(0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(2,1),(3,1),(0,2),(1,2),(2,2),(3,2),(0,3),(1,3),(2,3),(3,3)]),1e-12))
        self.assertTrue(m.getNumberFieldAtLevel(1).isEqual(DataArrayInt(list(range(10,26)))))
        self.assertTrue(m.getFamilyFieldAtLevel(1).isEqual(DataArrayInt([-1,-1,-1,-1,-1,-2,-2,-2,-2,-2,-2,0,-1,-3,-3,-3])))
        pass

    def testWRQPolyg1(self):
        fname="Pyfile72.med"
        m=MEDCoupling1SGTUMesh("mesh",NORM_QUAD4) ; m.allocateCells()
        m.insertNextCell([0,2,1,3])
        m.setCoords(DataArrayDouble([0.,0.,1.,1.,1.,0.,0.,1.],4,2))
        #
        ms = [m.deepCopy() for i in range(4)]
        for i,elt in enumerate(ms):
            elt.translate([float(i)*1.5,0.])
            pass
        m0=MEDCoupling1SGTUMesh.Merge1SGTUMeshes(ms).buildUnstructured()
        m0.convertAllToPoly()
        #
        ms = [m.deepCopy() for i in range(5)]
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
        f.setArray(arr) ; f.checkConsistencyLight()
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
        m0.checkConsistency()
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
        for i in range(nbCells):
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
        m.checkConsistency()
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        m1=MEDCouplingCMesh() ; m1.setCoords(arr) ; m1.setName("Mesh") 
        m1=m1.buildUnstructured() ; m1.setCoords(m.getCoords())
        mm.setMeshAtLevel(-1,m1)
        renum0=DataArrayInt([3,6,7,10,11,0,2,1,9,8,5,4,12,13,14,24,23,22,21,20,19,18,17,16,15])
        famField0=DataArrayInt([-3,-6,-7,-10,-11,0,-2,-1,-9,-8,-5,-4,-12,-13,-14,-24,-23,-22,-21,-20,-19,-18,-17,-16,-15])
        namesCellL0=DataArrayAsciiChar(25,16)
        namesCellL0[:] = ["Cell#%.3d        " % (i) for i in range(25)]
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
        namesNodes[:] = ["Node#%.3d        " % (i) for i in range(36)]
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
        m.checkConsistency()
        f=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f.setMesh(m)
        f.setName("Field")
        arr=DataArrayDouble(25,2) ; arr.setInfoOnComponents(compos)
        arr[:,0]=list(range(25))
        arr[:,1]=list(range(100,125))
        f.setArray(arr)
        WriteField(fileName,f,True)
        f=MEDCouplingFieldDouble(ON_NODES,ONE_TIME) ; f.setMesh(m)
        f.setName("FieldNode")
        arr=DataArrayDouble(36,2) ; arr.setInfoOnComponents(compos)
        arr[:,0]=list(range(200,236))
        arr[:,1]=list(range(300,336))
        f.setArray(arr)
        f.checkConsistencyLight()
        WriteFieldUsingAlreadyWrittenMesh(fileName,f)
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
        fs=fs.deepCopy()
        fs[0][0].loadArrays()
        arr = DataArrayDouble(12, 2) ; arr[:, 0] = list(range(3, 15)) ; arr[:, 1] = list(range(103, 115))
        arr.setInfoOnComponents(compos)
        self.assertTrue(fs[0][0].getUndergroundDataArray().isEqual(arr,1e-12))
        fs[1][0].loadArrays()
        arr = DataArrayDouble(21, 2) ; arr[:, 0] = list(range(203, 224)) ; arr[:, 1] = list(range(303, 324))
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
        m.checkConsistency()
        f=MEDCouplingFieldDouble(ON_CELLS,ONE_TIME) ; f.setMesh(m)
        f.setName("Field")
        arr=DataArrayDouble(25,2) ; arr.setInfoOnComponents(compos)
        arr[:,0]=list(range(25))
        arr[:,1]=list(range(100,125))
        f.setArray(arr)
        WriteField(fileName,f,True)
        f=MEDCouplingFieldDouble(ON_NODES,ONE_TIME) ; f.setMesh(m)
        f.setName("FieldNode")
        arr=DataArrayDouble(36,2) ; arr.setInfoOnComponents(compos)
        arr[:,0]=list(range(200,236))
        arr[:,1]=list(range(300,336))
        f.setArray(arr)
        f.checkConsistencyLight()
        WriteFieldUsingAlreadyWrittenMesh(fileName,f)
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
        fileName="Pyfile85.med"
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
        e=d.deltaShiftIndex().findIdsEqual(1)
        #
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m) ; mm.setMeshAtLevel(-1,mMinus1)
        grp0=DataArrayInt([0,1,2,3,4,5,24,25,26]) ; grp0.setName("grp0")
        mm.setGroupsAtLevel(0,[grp0])
        grp1=e ; grp1.setName("grp1")
        mm.setGroupsAtLevel(-1,[grp1])
        mm.write(fileName,2)
        #
        a=DataArrayDouble(3) ; a.iota()
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
        self.assertEqual(mm3D.getNumberOfCellsAtLevel(0),66)
        self.assertEqual(mm3D.getNumberOfCellsAtLevel(-1),194)
        self.assertEqual(mm3D.getGroupsNames(),('grp0','grp0_extruded','grp0_top','grp1','grp1_extruded','grp1_top'))
        self.assertEqual(mm3D.getGrpNonEmptyLevels("grp0"),(-1,))
        self.assertEqual(mm3D.getGrpNonEmptyLevels("grp0_top"),(-1,))
        self.assertEqual(mm3D.getGrpNonEmptyLevels("grp0_extruded"),(0,))
        self.assertEqual(mm3D.getGrpNonEmptyLevels("grp1"),(-2,))
        self.assertEqual(mm3D.getGrpNonEmptyLevels("grp1_top"),(-2,))
        self.assertEqual(mm3D.getGrpNonEmptyLevels("grp1_extruded"),(-1,))
        d=DataArrayDouble([(1.,0.,0.),(2.,0.,0.),(3.,0.,0.),(4.,0.,0.),(0.,1.,0.),(1.,1.,0.),(2.,1.,0.),(3.,1.,0.),(4.,1.,0.),(5.,1.,0.),(0.,2.,0.),(1.,2.,0.),(2.,2.,0.),(3.,2.,0.),(4.,2.,0.),(5.,2.,0.),(0.,3.,0.),(1.,3.,0.),(2.,3.,0.),(3.,3.,0.),(4.,3.,0.),(5.,3.,0.),(0.,4.,0.),(1.,4.,0.),(2.,4.,0.),(3.,4.,0.),(4.,4.,0.),(5.,4.,0.),(1.,5.,0.),(2.,5.,0.),(3.,5.,0.),(4.,5.,0.),(1.,0.,1.),(2.,0.,1.),(3.,0.,1.),(4.,0.,1.),(0.,1.,1.),(1.,1.,1.),(2.,1.,1.),(3.,1.,1.),(4.,1.,1.),(5.,1.,1.),(0.,2.,1.),(1.,2.,1.),(2.,2.,1.),(3.,2.,1.),(4.,2.,1.),(5.,2.,1.),(0.,3.,1.),(1.,3.,1.),(2.,3.,1.),(3.,3.,1.),(4.,3.,1.),(5.,3.,1.),(0.,4.,1.),(1.,4.,1.),(2.,4.,1.),(3.,4.,1.),(4.,4.,1.),(5.,4.,1.),(1.,5.,1.),(2.,5.,1.),(3.,5.,1.),(4.,5.,1.),(1.,0.,2.),(2.,0.,2.),(3.,0.,2.),(4.,0.,2.),(0.,1.,2.),(1.,1.,2.),(2.,1.,2.),(3.,1.,2.),(4.,1.,2.),(5.,1.,2.),(0.,2.,2.),(1.,2.,2.),(2.,2.,2.),(3.,2.,2.),(4.,2.,2.),(5.,2.,2.),(0.,3.,2.),(1.,3.,2.),(2.,3.,2.),(3.,3.,2.),(4.,3.,2.),(5.,3.,2.),(0.,4.,2.),(1.,4.,2.),(2.,4.,2.),(3.,4.,2.),(4.,4.,2.),(5.,4.,2.),(1.,5.,2.),(2.,5.,2.),(3.,5.,2.),(4.,5.,2.)])
        self.assertTrue(mm3D.getCoords().isEqual(d,1e-12))
        d=DataArrayInt([16,1,0,5,33,32,37,16,1,5,6,33,37,38,16,2,1,6,34,33,38,16,2,6,7,34,38,39,16,3,2,7,35,34,39,16,3,7,8,35,39,40,16,5,4,10,37,36,42,16,5,10,11,37,42,43,16,9,8,14,41,40,46,16,9,14,15,41,46,47,16,11,10,16,43,42,48,16,11,16,17,43,48,49,16,15,14,20,47,46,52,16,15,20,21,47,52,53,16,17,16,22,49,48,54,16,17,22,23,49,54,55,16,21,20,26,53,52,58,16,21,26,27,53,58,59,16,24,23,28,56,55,60,16,24,28,29,56,60,61,16,25,24,29,57,56,61,16,25,29,30,57,61,62,16,26,25,30,58,57,62,16,26,30,31,58,62,63,16,33,32,37,65,64,69,16,33,37,38,65,69,70,16,34,33,38,66,65,70,16,34,38,39,66,70,71,16,35,34,39,67,66,71,16,35,39,40,67,71,72,16,37,36,42,69,68,74,16,37,42,43,69,74,75,16,41,40,46,73,72,78,16,41,46,47,73,78,79,16,43,42,48,75,74,80,16,43,48,49,75,80,81,16,47,46,52,79,78,84,16,47,52,53,79,84,85,16,49,48,54,81,80,86,16,49,54,55,81,86,87,16,53,52,58,85,84,90,16,53,58,59,85,90,91,16,56,55,60,88,87,92,16,56,60,61,88,92,93,16,57,56,61,89,88,93,16,57,61,62,89,93,94,16,58,57,62,90,89,94,16,58,62,63,90,94,95,18,6,5,11,12,38,37,43,44,18,7,6,12,13,39,38,44,45,18,8,7,13,14,40,39,45,46,18,12,11,17,18,44,43,49,50,18,13,12,18,19,45,44,50,51,18,14,13,19,20,46,45,51,52,18,18,17,23,24,50,49,55,56,18,19,18,24,25,51,50,56,57,18,20,19,25,26,52,51,57,58,18,38,37,43,44,70,69,75,76,18,39,38,44,45,71,70,76,77,18,40,39,45,46,72,71,77,78,18,44,43,49,50,76,75,81,82,18,45,44,50,51,77,76,82,83,18,46,45,51,52,78,77,83,84,18,50,49,55,56,82,81,87,88,18,51,50,56,57,83,82,88,89,18,52,51,57,58,84,83,89,90])
        self.assertTrue(mm3D[0].getNodalConnectivity().isEqual(d))
        d=DataArrayInt([0,7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140,147,154,161,168,175,182,189,196,203,210,217,224,231,238,245,252,259,266,273,280,287,294,301,308,315,322,329,336,345,354,363,372,381,390,399,408,417,426,435,444,453,462,471,480,489,498])
        self.assertTrue(mm3D[0].getNodalConnectivityIndex().isEqual(d))
        d=DataArrayInt([3,1,0,5,3,1,5,6,3,2,1,6,3,2,6,7,3,3,2,7,3,3,7,8,3,5,4,10,3,5,10,11,3,9,8,14,3,9,14,15,3,11,10,16,3,11,16,17,3,15,14,20,3,15,20,21,3,17,16,22,3,17,22,23,3,21,20,26,3,21,26,27,3,24,23,28,3,24,28,29,3,25,24,29,3,25,29,30,3,26,25,30,3,26,30,31,3,65,64,69,3,65,69,70,3,66,65,70,3,66,70,71,3,67,66,71,3,67,71,72,3,69,68,74,3,69,74,75,3,73,72,78,3,73,78,79,3,75,74,80,3,75,80,81,3,79,78,84,3,79,84,85,3,81,80,86,3,81,86,87,3,85,84,90,3,85,90,91,3,88,87,92,3,88,92,93,3,89,88,93,3,89,93,94,3,90,89,94,3,90,94,95,4,1,0,32,33,4,0,5,37,32,4,5,1,33,37,4,5,6,38,37,4,6,1,33,38,4,2,1,33,34,4,6,2,34,38,4,6,7,39,38,4,7,2,34,39,4,3,2,34,35,4,7,3,35,39,4,7,8,40,39,4,8,3,35,40,4,5,4,36,37,4,4,10,42,36,4,10,5,37,42,4,10,11,43,42,4,11,5,37,43,4,9,8,40,41,4,8,14,46,40,4,14,9,41,46,4,14,15,47,46,4,15,9,41,47,4,10,16,48,42,4,16,11,43,48,4,16,17,49,48,4,17,11,43,49,4,14,20,52,46,4,20,15,47,52,4,20,21,53,52,4,21,15,47,53,4,16,22,54,48,4,22,17,49,54,4,22,23,55,54,4,23,17,49,55,4,20,26,58,52,4,26,21,53,58,4,26,27,59,58,4,27,21,53,59,4,24,23,55,56,4,23,28,60,55,4,28,24,56,60,4,28,29,61,60,4,29,24,56,61,4,25,24,56,57,4,29,25,57,61,4,29,30,62,61,4,30,25,57,62,4,26,25,57,58,4,30,26,58,62,4,30,31,63,62,4,31,26,58,63,4,11,12,44,43,4,12,6,38,44,4,12,13,45,44,4,13,7,39,45,4,13,14,46,45,4,17,18,50,49,4,18,12,44,50,4,18,19,51,50,4,19,13,45,51,4,19,20,52,51,4,24,18,50,56,4,25,19,51,57,4,33,32,64,65,4,32,37,69,64,4,37,33,65,69,4,37,38,70,69,4,38,33,65,70,4,34,33,65,66,4,38,34,66,70,4,38,39,71,70,4,39,34,66,71,4,35,34,66,67,4,39,35,67,71,4,39,40,72,71,4,40,35,67,72,4,37,36,68,69,4,36,42,74,68,4,42,37,69,74,4,42,43,75,74,4,43,37,69,75,4,41,40,72,73,4,40,46,78,72,4,46,41,73,78,4,46,47,79,78,4,47,41,73,79,4,42,48,80,74,4,48,43,75,80,4,48,49,81,80,4,49,43,75,81,4,46,52,84,78,4,52,47,79,84,4,52,53,85,84,4,53,47,79,85,4,48,54,86,80,4,54,49,81,86,4,54,55,87,86,4,55,49,81,87,4,52,58,90,84,4,58,53,85,90,4,58,59,91,90,4,59,53,85,91,4,56,55,87,88,4,55,60,92,87,4,60,56,88,92,4,60,61,93,92,4,61,56,88,93,4,57,56,88,89,4,61,57,89,93,4,61,62,94,93,4,62,57,89,94,4,58,57,89,90,4,62,58,90,94,4,62,63,95,94,4,63,58,90,95,4,43,44,76,75,4,44,38,70,76,4,44,45,77,76,4,45,39,71,77,4,45,46,78,77,4,49,50,82,81,4,50,44,76,82,4,50,51,83,82,4,51,45,77,83,4,51,52,84,83,4,56,50,82,88,4,57,51,83,89,4,6,5,11,12,4,7,6,12,13,4,8,7,13,14,4,12,11,17,18,4,13,12,18,19,4,14,13,19,20,4,18,17,23,24,4,19,18,24,25,4,20,19,25,26,4,70,69,75,76,4,71,70,76,77,4,72,71,77,78,4,76,75,81,82,4,77,76,82,83,4,78,77,83,84,4,82,81,87,88,4,83,82,88,89,4,84,83,89,90])
        self.assertTrue(mm3D[-1].getNodalConnectivity().isEqual(d))
        d=DataArrayInt([0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,192,197,202,207,212,217,222,227,232,237,242,247,252,257,262,267,272,277,282,287,292,297,302,307,312,317,322,327,332,337,342,347,352,357,362,367,372,377,382,387,392,397,402,407,412,417,422,427,432,437,442,447,452,457,462,467,472,477,482,487,492,497,502,507,512,517,522,527,532,537,542,547,552,557,562,567,572,577,582,587,592,597,602,607,612,617,622,627,632,637,642,647,652,657,662,667,672,677,682,687,692,697,702,707,712,717,722,727,732,737,742,747,752,757,762,767,772,777,782,787,792,797,802,807,812,817,822,827,832,837,842,847,852,857,862,867,872,877,882,887,892,897,902,907,912,917,922])
        self.assertTrue(mm3D[-1].getNodalConnectivityIndex().isEqual(d))
        d=DataArrayInt([1,1,0,1,0,5,1,5,1,1,5,6,1,6,1,1,2,1,1,6,2,1,6,7,1,7,2,1,3,2,1,7,3,1,7,8,1,8,3,1,5,4,1,4,10,1,10,5,1,10,11,1,11,5,1,9,8,1,8,14,1,14,9,1,14,15,1,15,9,1,10,16,1,16,11,1,16,17,1,17,11,1,14,20,1,20,15,1,20,21,1,21,15,1,16,22,1,22,17,1,22,23,1,23,17,1,20,26,1,26,21,1,26,27,1,27,21,1,24,23,1,23,28,1,28,24,1,28,29,1,29,24,1,25,24,1,29,25,1,29,30,1,30,25,1,26,25,1,30,26,1,30,31,1,31,26,1,11,12,1,12,6,1,12,13,1,13,7,1,13,14,1,17,18,1,18,12,1,18,19,1,19,13,1,19,20,1,24,18,1,25,19,1,65,64,1,64,69,1,69,65,1,69,70,1,70,65,1,66,65,1,70,66,1,70,71,1,71,66,1,67,66,1,71,67,1,71,72,1,72,67,1,69,68,1,68,74,1,74,69,1,74,75,1,75,69,1,73,72,1,72,78,1,78,73,1,78,79,1,79,73,1,74,80,1,80,75,1,80,81,1,81,75,1,78,84,1,84,79,1,84,85,1,85,79,1,80,86,1,86,81,1,86,87,1,87,81,1,84,90,1,90,85,1,90,91,1,91,85,1,88,87,1,87,92,1,92,88,1,92,93,1,93,88,1,89,88,1,93,89,1,93,94,1,94,89,1,90,89,1,94,90,1,94,95,1,95,90,1,75,76,1,76,70,1,76,77,1,77,71,1,77,78,1,81,82,1,82,76,1,82,83,1,83,77,1,83,84,1,88,82,1,89,83])
        self.assertTrue(mm3D[-2].getNodalConnectivity().isEqual(d))
        d=DataArrayInt(129) ; d.iota() ; d*=3
        self.assertTrue(mm3D[-2].getNodalConnectivityIndex().isEqual(d))
        #
        self.assertEqual(mm3D.getGroupArr(-1,"grp0").getName(),"grp0")
        self.assertEqual(mm3D.getGroupArr(-2,"grp1").getName(),"grp1")
        self.assertTrue(mm3D.getGroupArr(-1,"grp0").isEqualWithoutConsideringStr(DataArrayInt([0,1,2,3,4,5,176,177,178])))
        self.assertTrue(mm3D.getGroupArr(-1,"grp0_top").isEqualWithoutConsideringStr(DataArrayInt([24,25,26,27,28,29,185,186,187])))
        self.assertTrue(mm3D.getGroupArr(-2,"grp1").isEqualWithoutConsideringStr(DataArrayInt([0,1,5,9,12,13,14,18,22,23,30,31,33,37,38,40,42,46,50,51])))
        self.assertTrue(mm3D.getGroupArr(-2,"grp1_top").isEqualWithoutConsideringStr(DataArrayInt([64,65,69,73,76,77,78,82,86,87,94,95,97,101,102,104,106,110,114,115])))
        self.assertTrue(mm3D.getGroupArr(0,"grp0_extruded").isEqualWithoutConsideringStr(DataArrayInt([0,1,2,3,4,5,24,25,26,27,28,29,48,49,50,57,58,59])))
        self.assertTrue(mm3D.getGroupArr(-1,"grp1_extruded").isEqualWithoutConsideringStr(DataArrayInt([48,49,53,57,60,61,62,66,70,71,78,79,81,85,86,88,90,94,98,99,112,113,117,121,124,125,126,130,134,135,142,143,145,149,150,152,154,158,162,163])))
        mm3D.setName("MeshExtruded")
        mm3D.write(fileName,0)
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def testMEDFileUMeshPickeling1(self):
        outFileName="Pyfile86.med"
        c=DataArrayDouble([-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 ],9,2)
        c.setInfoOnComponents(["aa","bbb"])
        targetConn=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
        m=MEDCouplingUMesh();
        m.setMeshDimension(2);
        m.allocateCells(5);
        m.insertNextCell(NORM_TRI3,3,targetConn[4:7])
        m.insertNextCell(NORM_TRI3,3,targetConn[7:10])
        m.insertNextCell(NORM_QUAD4,4,targetConn[0:4])
        m.insertNextCell(NORM_POLYGON,4,targetConn[10:14])
        m.insertNextCell(NORM_POLYGON,4,targetConn[14:18])
        m.finishInsertingCells();
        m.setCoords(c)
        m.checkConsistencyLight()
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(1);
        m1.allocateCells(3);
        m1.insertNextCell(NORM_SEG2,2,[1,4])
        m1.insertNextCell(NORM_SEG2,2,[3,6])
        m1.insertNextCell(NORM_SEG3,3,[2,8,5])
        m1.finishInsertingCells();
        m1.setCoords(c)
        m1.checkConsistencyLight()
        m2=MEDCouplingUMesh.New();
        m2.setMeshDimension(0);
        m2.allocateCells(4);
        m2.insertNextCell(NORM_POINT1,1,[1])
        m2.insertNextCell(NORM_POINT1,1,[3])
        m2.insertNextCell(NORM_POINT1,1,[2])
        m2.insertNextCell(NORM_POINT1,1,[6])
        m2.finishInsertingCells();
        m2.setCoords(c)
        m2.checkConsistencyLight()
        #
        mm=MEDFileUMesh.New()
        self.assertTrue(mm.getUnivNameWrStatus())
        mm.setName("MyFirstMEDCouplingMEDmesh")
        mm.setDescription("IHopeToConvinceLastMEDMEMUsers")
        mm.setCoords(c)
        mm[-1]=m1;
        mm[0]=m;
        mm.setRenumFieldArr(0,DataArrayInt([32,41,50,56,7]))
        mm[-2]=m2;
        mm.setRenumFieldArr(-2,DataArrayInt([102,52,45,63]))
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
        g1_N.setValues(list(range(8)),8,1)
        g1_N.setName("G1")
        g2_N=DataArrayInt.New()
        g2_N.setValues(list(range(9)),9,1)
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
        #
        st=pickle.dumps(mm,pickle.HIGHEST_PROTOCOL)
        mm2=pickle.loads(st)
        self.assertTrue(mm.isEqual(mm2,1e-12)[0])
        self.assertEqual(mm.getAxisType(),AX_CART)
        #
        mm.setAxisType(AX_CYL)
        st=pickle.dumps(mm,pickle.HIGHEST_PROTOCOL)
        mm2=pickle.loads(st)
        self.assertTrue(mm.isEqual(mm2,1e-12)[0])
        self.assertEqual(mm2.getAxisType(),AX_CYL)
        pass

    def testMEDFileFieldsLoadSpecificEntities1(self):
        nbNodes=11
        fieldName="myField"
        fileName="Pyfile87.med"
        nbPdt=10
        meshName="Mesh"
        #
        m=MEDCouplingCMesh()
        arr=DataArrayDouble(nbNodes) ; arr.iota()
        m.setCoords(arr)
        m=m.buildUnstructured()
        m.setName(meshName)
        #
        fmts=MEDFileFieldMultiTS()
        for i in range(nbPdt):
            f=MEDCouplingFieldDouble(ON_NODES)
            f.setMesh(m)
            arr=DataArrayDouble(nbNodes) ; arr.iota() ; arr*=i
            f.setArray(arr)
            f.setName(fieldName)
            f.setTime(float(i),i,0)
            fmts.appendFieldNoProfileSBT(f)
            pass
        #
        mm=MEDFileUMesh() ; mm[0]=m
        fmts.write(fileName,2)
        mm.write(fileName,0)
        #
        fs=MEDFileFields(fileName,False)
        fs2=MEDFileFields.LoadSpecificEntities(fileName,[(ON_NODES,NORM_ERROR)],False)
        fs.loadArraysIfNecessary()
        fs2.loadArraysIfNecessary()
        for i in range(nbPdt):
            self.assertTrue(fs[fieldName][i].getUndergroundDataArray().isEqual(fs2[fieldName][i].getUndergroundDataArray(),1e-12))
            pass
        m1=MEDCouplingCMesh() ; m1.setCoords(DataArrayDouble([0,1,2,3]),DataArrayDouble([0,1])) ; m1=m1.buildUnstructured() ; m1.simplexize(0)
        m2=MEDCouplingCMesh() ; m2.setCoords(DataArrayDouble([3,4,5]),DataArrayDouble([0,1])) ; m2=m2.buildUnstructured()
        m3=MEDCouplingUMesh.MergeUMeshes(m1,m2) ; m3.setName(meshName)
        fmts=MEDFileFieldMultiTS()
        for i in range(nbPdt):
            f=MEDCouplingFieldDouble(ON_CELLS)
            f.setMesh(m3)
            arr=DataArrayDouble(8) ; arr.iota() ; arr*=i
            f.setArray(arr)
            f.setName(fieldName)
            f.setTime(float(i),i,0)
            fmts.appendFieldNoProfileSBT(f)
            pass
        mm=MEDFileUMesh() ; mm[0]=m3
        del mm[0]
        self.assertEqual(mm.getNonEmptyLevels(),())
        mm[0]=m3
        self.assertEqual(mm.getNonEmptyLevels(),(0,))
        fmts.write(fileName,2)
        fs=MEDFileFields(fileName,False)
        fs2=MEDFileFields.LoadSpecificEntities(fileName,[(ON_CELLS,NORM_TRI3)],False)
        fs3=MEDFileFieldMultiTS.LoadSpecificEntities(fileName,fieldName,[(ON_CELLS,NORM_QUAD4)],False)
        fs4=MEDFileFields.LoadSpecificEntities(fileName,[(ON_CELLS,NORM_TRI3),(ON_CELLS,NORM_QUAD4)],False)
        fs.loadArraysIfNecessary()
        fs2.loadArraysIfNecessary()
        fs3.loadArraysIfNecessary()
        fs4.loadArraysIfNecessary()
        for i in range(nbPdt):
            self.assertTrue(fs[fieldName][i].getUndergroundDataArray()[:6].isEqual(fs2[fieldName][i].getUndergroundDataArray(),1e-12))
            self.assertTrue(fs[fieldName][i].getUndergroundDataArray()[6:8].isEqual(fs3[i].getUndergroundDataArray(),1e-12))
            self.assertTrue(fs[fieldName][i].getUndergroundDataArray().isEqual(fs4[fieldName][i].getUndergroundDataArray(),1e-12))
            pass
        pass

    def testMEDFileLotsOfTSRW1(self):
        nbNodes=11
        fieldName="myField"
        fileName="Pyfile88.med"
        nbPdt=300 # <- perftest = 30000
        meshName="Mesh"
        #
        maxPdt=100 # <- optimum = 500
        m=MEDCouplingCMesh()
        arr=DataArrayDouble(nbNodes) ; arr.iota()
        m.setCoords(arr)
        m=m.buildUnstructured()
        m.setName(meshName)
        #
        nbOfField=nbPdt//maxPdt
        fs=MEDFileFields()
        for j in range(nbOfField):
            fmts=MEDFileFieldMultiTS()
            s=DataArray.GetSlice(slice(0,nbPdt,1),j,nbOfField)
            for i in range(s.start, s.stop, s.step):
                f=MEDCouplingFieldDouble(ON_NODES)
                f.setMesh(m)
                arr=DataArrayDouble(nbNodes) ; arr.iota() ; arr*=i
                f.setArray(arr)
                f.setName("%s_%d"%(fieldName,j))
                f.setTime(float(i),i,0)
                fmts.appendFieldNoProfileSBT(f)
                pass
            fs.pushField(fmts)
            pass
        #
        mm=MEDFileUMesh() ; mm[0]=m
        fs.write(fileName,2)
        mm.write(fileName,0)
        ############
        def appendInDict(d,key,val):
            if key in d:
                d[key].append(val)
            else:
                d[key]=[val]
            pass
        import re
        allFields=GetAllFieldNames(fileName)
        allFieldsDict={}
        pat=re.compile("([\d]+)([\s\S]+)$")
        for st in allFields:
            stRev=st[::-1]
            m=pat.match(stRev)
            if m:
                appendInDict(allFieldsDict,m.group(2)[::-1],m.group(1)[::-1])
                pass
            else:
                appendInDict(allFieldsDict,st,'')
                pass
            pass
        fs2=MEDFileFields()
        for k in allFieldsDict:
            if allFieldsDict[k]!=['']:
                allFieldsDict[k]=sorted(allFieldsDict[k],key=lambda x: int(x))
                pass
            fmts2=[]
            for it in allFieldsDict[k]:
                fmts2.append(MEDFileFieldMultiTS.LoadSpecificEntities(fileName,k+it,[(ON_NODES,NORM_ERROR)]))
                pass
            fmts2.reverse()
            zeResu=fmts2.pop()
            nbIter=len(fmts2)
            for ii in range(nbIter):
                zeResu.pushBackTimeSteps(fmts2.pop())
                pass
            zeResu.setName(k)
            fs2.pushField(zeResu)
            pass
        self.assertEqual(fs2[0].getTimeSteps(), [(i, 0, float(i)) for i in range(nbPdt)])
        pass
    
    def testMEDFileMeshRearrangeFamIds1(self):
        """ Test for bug EDF10720. The aim of this test is the call of MEDFileMesh.rearrangeFamilies."""
        fileName="Pyfile89.med"
        meshName='Maillage_2'
        mm=MEDFileUMesh()
        coords=DataArrayDouble([(0.,0.,0.),(0.,0.,200.),(0.,200.,200.),(0.,200.,0.),(200.,0.,0.),(200.,0.,200.),(200.,200.,200.),(200.,200.,0.),(0.,0.,100.),(0.,100.,200.),(0.,200.,100.),(0.,100.,0.),(200.,0.,100.),(200.,100.,200.),(200.,200.,100.),(200.,100.,0.),(100.,0.,0.),(100.,0.,200.),(100.,200.,0.),(100.,200.,200.),(0.,116.87743909766768,83.12256090233232),(200.,116.87743909766768,83.12256090233232),(116.87743909766769,0.,116.87743909766769),(116.87743909766769,200.,116.87743909766769),(116.87743909766769,116.87743909766769,0.),(116.87743909766769,116.87743909766769,200.),(63.3851584383713,56.1391811199829,119.728314479261),(138.008709441123,116.039297556044,119.903790959468)])
        #
        c0=DataArrayInt([14,1,26,9,8,14,17,26,1,8,14,27,26,17,22,14,26,16,20,8,14,8,0,16,11,14,16,20,11,24,14,25,20,26,27,14,22,26,24,27,14,26,16,22,24,14,8,26,22,17,14,20,9,25,26,14,19,20,25,23,14,23,6,27,25,14,19,23,10,20,14,27,22,21,24,14,27,21,14,18,14,26,9,25,17,14,13,27,25,17,14,27,18,24,21,14,22,21,15,12,14,27,20,24,18,14,23,25,27,20,14,13,27,6,25,14,23,27,6,14,14,15,16,22,12,14,27,17,13,22,14,22,27,21,13,14,24,16,22,15,14,24,18,7,21,14,12,4,15,16,14,22,12,5,13,14,8,26,16,22,14,13,27,21,14,14,20,18,10,3,14,14,27,18,23,14,14,27,6,13,14,21,22,13,12,14,25,26,17,27,14,19,9,25,20,14,26,24,20,16,14,22,24,15,21,14,9,26,1,17,14,23,27,18,20,14,20,11,18,3,14,14,18,21,7,14,19,2,9,10,14,19,23,25,6,14,18,23,20,10,14,20,26,8,9,14,22,13,5,17,14,24,11,18,20,14,21,15,7,24,14,19,20,10,9,14,20,26,27,24,14,16,8,11,20])
        c0i=DataArrayInt([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275])
        m0=MEDCouplingUMesh(meshName,3) ; m0.setCoords(coords)
        m0.setConnectivity(c0,c0i)
        mm[0]=m0
        #
        c1=DataArrayInt([3,8,20,11,3,8,9,20,3,9,2,10,3,20,9,10,3,0,8,11,3,9,8,1,3,20,10,3,3,11,20,3,3,15,21,12,3,5,12,13,3,21,13,12,3,15,12,4,3,14,6,13,3,14,13,21,3,7,14,21,3,7,21,15,3,5,22,12,3,4,12,16,3,17,1,8,3,16,8,0,3,5,17,22,3,12,22,16,3,22,17,8,3,16,22,8,3,10,2,19,3,7,18,14,3,14,23,6,3,3,10,18,3,23,19,6,3,18,23,14,3,10,19,23,3,10,23,18,3,3,18,11,3,7,24,18,3,15,4,16,3,11,16,0,3,7,15,24,3,18,24,11,3,24,15,16,3,11,24,16,3,9,19,2,3,19,25,6,3,17,5,13,3,1,17,9,3,25,13,6,3,9,25,19,3,17,13,25,3,17,25,9])
        c1i=DataArrayInt([0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,192])
        m1=MEDCouplingUMesh(meshName,2) ; m1.setCoords(coords)
        m1.setConnectivity(c1,c1i)
        mm[-1]=m1
        #
        c2=DataArrayInt([0,8,8,1,1,9,9,2,3,10,10,2,0,11,11,3,4,12,12,5,5,13,13,6,7,14,14,6,4,15,15,7,0,16,16,4,1,17,17,5,3,18,18,7,2,19,19,6])
        m2=MEDCoupling1SGTUMesh(meshName,NORM_SEG2)
        m2.setNodalConnectivity(c2) ; m2.setCoords(coords)
        mm[-2]=m2.buildUnstructured()
        #
        ref0=DataArrayInt(55) ; ref0[:]=0
        mm.setFamilyFieldArr(0,ref0)
        mm.setFamilyFieldArr(1,DataArrayInt([0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]))
        ref1=DataArrayInt([0,0,0,0,0,0,0,0,-6,-6,-6,-6,-6,-6,-6,-6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        mm.setFamilyFieldArr(-1,ref1)
        ref2=DataArrayInt([0,0,-7,-7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        mm.setFamilyFieldArr(-2,ref2)
        #
        for f,fid in (('FAMILLE_ZERO',0),('FAM_-6_Groupe_1',-6),('FAM_-7_Groupe_2',-7),('FAM_2_Groupe_3',2)):
            mm.setFamilyId(f,fid)
        for grp,fams in [('Groupe_1',('FAM_-6_Groupe_1',)),('Groupe_2',('FAM_-7_Groupe_2',)),('Groupe_3',('FAM_2_Groupe_3',))]:
            mm.setFamiliesOnGroup(grp,fams)
        mm.write(fileName,2)
        #
        mm=MEDFileMesh.New(fileName)
        grp=mm.getGroup(-1,"Groupe_1")
        dai=grp.computeFetchedNodeIds()
        dai.setName("TOTO")
        mm.addGroup(1,dai)
        mm.rearrangeFamilies() # <- the aim of the test
        self.assertTrue(dai.isEqual(mm.getGroupArr(1,"TOTO")))
        self.assertTrue(mm.getFamilyFieldAtLevel(0).isEqual(ref0))
        self.assertTrue(mm.getFamilyFieldAtLevel(-1).isEqual(ref1))
        self.assertTrue(mm.getFamilyFieldAtLevel(-2).isEqual(ref2))
        self.assertTrue(mm.getFamilyFieldAtLevel(1).isEqual(DataArrayInt([0,0,2,0,9,9,9,9,0,0,0,0,9,9,9,9,0,0,0,0,0,9,0,0,0,0,0,0])))
        allGrps=[('Groupe_1',('FAM_-6_Groupe_1',)),('Groupe_2',('FAM_-7_Groupe_2',)),('Groupe_3',('FAM_2_Groupe_3',)),('TOTO',('Family_9',))]
        allFams=[('FAMILLE_ZERO',0),('FAM_-6_Groupe_1',-6),('FAM_-7_Groupe_2',-7),('FAM_2_Groupe_3',2),('Family_9',9)]
        self.assertEqual(list(mm.getGroupsNames()),[elt[0] for elt in allGrps])
        for elt,fams in allGrps:
            self.assertEqual(mm.getFamiliesOnGroup(elt),fams)
        self.assertEqual(list(mm.getFamiliesNames()),[elt[0] for elt in allFams])
        for elt,eltId in allFams:
            self.assertEqual(mm.getFamilyId(elt),eltId)
        pass

    def testNonRegrCMeshSetFieldPfl1(self):
        """ Non regression test. For structured mesh, push a false partial field in MEDFileField1TS using setFieldProfile."""
        ff=MEDFileField1TS()
        meshName="mesh"
        mm=MEDFileCMesh()
        m=MEDCouplingCMesh() ; arr=DataArrayDouble(5) ; arr.iota()
        m.setCoords(arr)
        m.setName(meshName)
        mm.setMesh(m)
        field=MEDCouplingFieldDouble(ON_CELLS)
        field.setMesh(m)
        field.setArray(DataArrayDouble([1.2,2.3,3.4,4.5]))
        field.setName("Field")
        field.checkConsistencyLight()
        pfl=DataArrayInt([0,1,2,3]) ; pfl.setName("TUTU") #<- false profile because defined on all cells !
        ff.setFieldProfile(field,mm,0,pfl) # <- bug was revealed here !
        self.assertEqual(ff.getPfls(),())
        field2=ff.getFieldOnMeshAtLevel(ON_CELLS,0,mm)
        self.assertTrue(field.isEqual(field2,1e-12,1e-12))
        del ff,mm,field,field2,pfl
        # same with unstructured mesh
        ff=MEDFileField1TS()
        meshName="mesh"
        mm=MEDFileUMesh()
        m=MEDCouplingCMesh() ; arr=DataArrayDouble(5) ; arr.iota()
        m.setCoords(arr)
        m.setName(meshName)
        m=m.buildUnstructured()
        mm[0]=m
        field=MEDCouplingFieldDouble(ON_CELLS)
        field.setMesh(m)
        field.setArray(DataArrayDouble([1.2,2.3,3.4,4.5]))
        field.setName("Field")
        field.checkConsistencyLight()
        pfl=DataArrayInt([0,1,2,3]) ; pfl.setName("TUTU")
        ff.setFieldProfile(field,mm,0,pfl)
        self.assertEqual(ff.getPfls(),())
        field2=ff.getFieldOnMeshAtLevel(ON_CELLS,0,mm)
        self.assertTrue(field.isEqual(field2,1e-12,1e-12))
        pass

    def testMEDFileUMeshLinearToQuadraticAndRev1(self):
        meshName="mesh"
        fileName="Pyfile90.med"
        fileName2="Pyfile91.med"
        arr=DataArrayDouble(5) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
        m=m.buildUnstructured()
        d=DataArrayInt([3,7,11,15])
        m1=m[d]
        m1.simplexize(0)
        m2=m[d.buildComplement(m.getNumberOfCells())]
        m=MEDCouplingUMesh.MergeUMeshesOnSameCoords(m1,m2)
        m.changeSpaceDimension(3,0.)
        arr=DataArrayDouble(3) ; arr.iota()
        m1D=MEDCouplingCMesh() ; m1D.setCoords(arr) ; m1D=m1D.buildUnstructured() ; m1D.changeSpaceDimension(3,0.)
        m1D.setCoords(m1D.getCoords()[:,[1,2,0]])
        delta=m.getNumberOfNodes()*(m1D.getNumberOfNodes()-1)
        m3D=m.buildExtrudedMesh(m1D,0)
        m3D.sortCellsInMEDFileFrmt()
        m3D.setName(meshName)
        m2D=m ; m2D.setCoords(m3D.getCoords()) ; m2D.shiftNodeNumbersInConn(delta) ; m2D.setName(meshName) ; m2D.checkConsistency()
        m1D=m2D.computeSkin() ; m1D.setName(meshName)
        m0D=MEDCouplingUMesh.Build0DMeshFromCoords(m3D.getCoords()) ; m0D.setName(meshName) ; m0D=m0D[[2,4,10]]
        #
        mm=MEDFileUMesh()
        mm[0]=m3D ; mm[-1]=m2D ; mm[-2]=m1D ; mm[-3]=m0D
        grpEdge0=DataArrayInt([1,2,3,5]) ; grpEdge0.setName("East")
        grpEdge1=DataArrayInt([0,1]) ; grpEdge1.setName("Corner1")
        grpFaceSouth=DataArrayInt([0,1,8,9,10]) ; grpFaceSouth.setName("SouthFace")
        grpFaceNorth=DataArrayInt([6,7,17,18,19]) ; grpFaceNorth.setName("NorthFace")
        diagFace=DataArrayInt([0,1,13,15,17]) ; diagFace.setName("DiagFace")
        vol1=DataArrayInt([20,21,23,24]) ; vol1.setName("vol1")
        vol2=DataArrayInt([2,3,4,5,21,24]) ; vol2.setName("vol2")
        mm.setGroupsAtLevel(0,[vol1,vol2])
        mm.setGroupsAtLevel(-1,[grpFaceSouth,grpFaceNorth,diagFace])
        mm.setGroupsAtLevel(-2,[grpEdge0,grpEdge1])
        #
        mmOut1=mm.linearToQuadratic(0,0.)
        mmOut1.write(fileName2,2)
        mmOut2=mmOut1.quadraticToLinear(0.)
        self.assertTrue(mm.isEqual(mmOut2,1e-12)[0])
        pass

    def testMEDFileMeshAddGroup1(self):
        m=MEDCouplingCMesh()
        arrX=DataArrayDouble(9) ; arrX.iota()
        arrY=DataArrayDouble(4) ; arrY.iota()
        m.setCoords(arrX,arrY)
        m.setName("mesh")
        mm=MEDFileCMesh()
        mm.setMesh(m)
        grp0=DataArrayInt([3,5,6,21,22]) ; grp0.setName("grp0")
        mm.addGroup(0,grp0)
        grp1=DataArrayInt([3,4,5,8,18,19,22]) ; grp1.setName("grp1")
        mm.addGroup(0,grp1)
        grp2=DataArrayInt([0,1,2,10,11]) ; grp2.setName("grp2")
        mm.addGroup(0,grp2)
        grp3=DataArrayInt([23]) ; grp3.setName("grp3")
        mm.addGroup(0,grp3)
        for grp in [grp0,grp1,grp2,grp3]:
            self.assertTrue(mm.getGroupArr(0,grp.getName()).isEqual(grp))
        self.assertEqual(mm.getGroupsNames(),('grp0','grp1','grp2','grp3'))
        delta=12
        for grp in [grp0,grp1,grp2,grp3]:
            grpNode=grp.deepCopy() ; grpNode+=delta ; grpNode.setName("%s_node"%grp.getName())
            mm.addGroup(1,grpNode)
        self.assertEqual(mm.getGroupsNames(),('grp0','grp0_node','grp1','grp1_node','grp2','grp2_node','grp3','grp3_node'))
        for grp in [grp0,grp1,grp2,grp3]:
            self.assertTrue(mm.getGroupArr(0,grp.getName()).isEqual(grp))
        for grp in [grp0,grp1,grp2,grp3]:
            grpExp=grp+delta ; grpExp.setName("%s_node"%grp.getName())
            self.assertTrue(mm.getGroupArr(1,"%s_node"%grp.getName()).isEqual(grpExp))
        mm.normalizeFamIdsMEDFile()
        for grp in [grp0,grp1,grp2,grp3]:
            self.assertTrue(mm.getGroupArr(0,grp.getName()).isEqual(grp))
        for grp in [grp0,grp1,grp2,grp3]:
            grpExp=grp+delta ; grpExp.setName("%s_node"%grp.getName())
            self.assertTrue(mm.getGroupArr(1,"%s_node"%grp.getName()).isEqual(grpExp))
        pass
    
    def testMEDFileJoint1(self):
        fileName="Pyfile92.med"
        coo=DataArrayDouble([(0,0,0),(1,0,0),(2,0,0)])
        coo.setInfoOnComponents(["x [cm]","y [cm]","z [cm]"])
        mm=MEDFileUMesh()
        mm.setCoords(coo)
        mm.setName("maa1")
        mm.setDescription("un maillage")
        mm.write(fileName,2)
        node_correspond=MEDFileJointCorrespondence(DataArrayInt([1,2,3,4,5,6,7,8]))
        cell_correspond=MEDFileJointCorrespondence(DataArrayInt([9,10,11,12]),NORM_TRI3,NORM_TRI3)
        one_step_joint=MEDFileJointOneStep()
        one_step_joint.pushCorrespondence(cell_correspond)
        one_step_joint.pushCorrespondence(node_correspond)
        one_joint=MEDFileJoint()
        one_joint.pushStep(one_step_joint)
        one_joint.setLocalMeshName("maa1")
        one_joint.setRemoteMeshName("maa1")
        one_joint.setDescription("joint_description")
        one_joint.setJointName("joint_1")
        one_joint.setDomainNumber(1)
        self.assertEqual( one_joint.getLocalMeshName(), "maa1")
        self.assertEqual( one_joint.getRemoteMeshName(), "maa1")
        self.assertEqual( one_joint.getDescription(), "joint_description")
        self.assertEqual( one_joint.getJointName(), "joint_1")
        self.assertEqual( one_joint.getDomainNumber(), 1)
        joints=MEDFileJoints()
        joints.pushJoint(one_joint);
        joints.write(fileName,0)
        # read back
        jointsR=MEDFileJoints(fileName,mm.getName())
        self.assertEqual( jointsR.getNumberOfJoints(), 1 )
        jR = jointsR.getJointAtPos(0)
        self.assertTrue( jR.isEqual( one_joint ))
        self.assertRaises( InterpKernelException, jointsR.getJointAtPos,1)
        self.assertRaises( InterpKernelException, jointsR.destroyJointAtPos,1)
        jointsR.destroyJointAtPos(0)
        pass
    
    def testMEDFileJoint2(self):
        fileNameWr="Pyfile93.med"
        coo=DataArrayDouble([(0,0,0),(1,0,0),(2,0,0)])
        coo.setInfoOnComponents(["x [cm]","y [cm]","z [cm]"])
        mm=MEDFileUMesh()
        mm.setCoords(coo)
        mm.setName("maa1")
        mm.setDescription("un maillage")
        node_correspond=MEDFileJointCorrespondence(DataArrayInt([13,14,15,16]))
        cell_correspond=MEDFileJointCorrespondence(DataArrayInt([17,18]),NORM_TETRA4,NORM_PENTA6)
        one_step_joint=MEDFileJointOneStep()
        two_step_joint=MEDFileJointOneStep()
        one_joint=MEDFileJoint()
        two_joint=MEDFileJoint()
        one_step_joint.pushCorrespondence(node_correspond)
        one_joint.pushStep(one_step_joint)
        two_step_joint.pushCorrespondence(cell_correspond)
        two_step_joint.pushCorrespondence(node_correspond)
        two_joint.pushStep(two_step_joint)
        one_joint.setLocalMeshName("maa1")
        one_joint.setRemoteMeshName("maa1")
        one_joint.setDescription("joint_description_1")
        one_joint.setJointName("joint_1")
        one_joint.setDomainNumber(1)
        two_joint.setLocalMeshName("maa1")
        two_joint.setRemoteMeshName("maa1")
        two_joint.setDescription("joint_description_2")
        two_joint.setJointName("joint_2")
        two_joint.setDomainNumber(2)
        joints=MEDFileJoints()
        joints.pushJoint(one_joint)
        joints.pushJoint(two_joint)
        mm.setJoints( joints )
        mm.write(fileNameWr,2)
        #
        mm=MEDFileMesh.New(fileNameWr)
        self.assertEqual( mm.getNumberOfJoints(), 2)
        jointsR = mm.getJoints();
        self.assertEqual( jointsR.getMeshName(), mm.getName() )
        self.assertEqual( len( jointsR ), 2 )
        jointR1 = jointsR[0]
        jointR2 = jointsR[1]
        self.assertFalse( jointR1 is None )
        self.assertFalse( jointR2 is None )
        self.assertTrue( jointR1.isEqual( one_joint ))
        self.assertTrue( jointR2.isEqual( two_joint ))
        pass

    def testMEDFileJoint1(self):
        node_correspond=MEDFileJointCorrespondence(DataArrayInt([1,2,3,4,5,6,7,8]))
        cell_correspond=MEDFileJointCorrespondence(DataArrayInt([9,10,11,12]),NORM_TRI3,NORM_TRI3)
        cell_correspon2=MEDFileJointCorrespondence(DataArrayInt([9,10,11]),NORM_TRI3,NORM_TRI3)
        cell_correspon3=MEDFileJointCorrespondence(DataArrayInt([9,10,11,12]),NORM_TRI3,NORM_QUAD4)
        joint1st_1=MEDFileJointOneStep()
        joint1st_1.pushCorrespondence(cell_correspond)
        joint1st_1.pushCorrespondence(node_correspond)
        joint1st_2=MEDFileJointOneStep()
        joint1st_2.pushCorrespondence(cell_correspond)
        joint1st_2.pushCorrespondence(node_correspond)
        joint1st_3=MEDFileJointOneStep()
        joint1st_3.pushCorrespondence(node_correspond)
        joint1st_3.pushCorrespondence(cell_correspond)
        joint1st_4=MEDFileJointOneStep()
        joint1st_4.pushCorrespondence(cell_correspond)
        joint1st_5=MEDFileJointOneStep()
        joint1st_5.pushCorrespondence(cell_correspon2)
        joint1st_6=MEDFileJointOneStep()
        joint1st_6.pushCorrespondence(cell_correspon3)
        self.assertTrue( joint1st_1.isEqual( joint1st_2 ))
        self.assertTrue( joint1st_1.isEqual( joint1st_3 ))
        self.assertFalse( joint1st_1.isEqual( joint1st_4 ))
        self.assertFalse( joint1st_4.isEqual( joint1st_5 ))
        self.assertFalse( joint1st_4.isEqual( joint1st_6 ))
        one_joint=MEDFileJoint()
        one_joint.pushStep(joint1st_1)
        one_joint.setLocalMeshName("maa1")
        one_joint.setRemoteMeshName("maa2")
        one_joint.setDescription("joint_description")
        one_joint.setJointName("joint_1")
        one_joint.setDomainNumber(1)
        self.assertEqual( "maa1", one_joint.getLocalMeshName())
        self.assertEqual( "maa2", one_joint.getRemoteMeshName())
        self.assertEqual( "joint_description", one_joint.getDescription())
        self.assertEqual( 1, one_joint.getDomainNumber())
        self.assertEqual( "joint_1", one_joint.getJointName())
        pass

    @unittest.skipUnless('linux'==platform.system().lower(),"stderr redirection not ported on Windows ?")
    def testMEDFileSafeCall0(self):
        """ EDF11242 : check status of MED file calls to detect problems immediately. Sorry this test generates awful messages !"""
        fname="Pyfile94.med"
        errfname="Pyfile94.err"
        class StdOutRedirect(object):
            def __init__(self,fileName):
                import os,sys
                sys.stderr.flush()
                self.stdoutOld=os.dup(2)
                self.fdOfSinkFile=os.open(fileName,os.O_CREAT | os.O_RDWR)
                fd2=os.dup2(self.fdOfSinkFile,2)
                self.origPyVal=sys.stderr
                class FlushFile(object):
                    def __init__(self,f):
                        self.f=f
                    def write(self,st):
                        self.f.write(st)
                        self.f.flush()
                    def flush(self):
                        return self.f.flush()
                    def isatty(self):
                        return self.f.isatty()
                sys.stderr=FlushFile(os.fdopen(self.fdOfSinkFile,"w"))
            def __del__(self):
                import os,sys
                sys.stderr=self.origPyVal
                if sys.version_info.major >= 3:
                    self.fdOfSinkFile.close()
                    pass
                #os.fsync(self.fdOfSinkFile)
                os.fsync(2)
                os.dup2(self.stdoutOld,2)
                os.close(self.stdoutOld)
        import os
        # first clean file if needed
        if os.path.exists(fname):
            os.remove(fname)
            pass
        # second : build a file from scratch
        m=MEDCouplingCMesh()
        arr=DataArrayDouble(11) ; arr.iota()
        m.setCoords(arr,arr)
        mm=MEDFileCMesh()
        mm.setMesh(m)
        mm.setName("mesh")
        mm.write(fname,2)
        # third : change permissions to remove write access on created file
        os.chmod(fname, 0o444)
        # four : try to append data on file -> check that it raises Exception
        f=MEDCouplingFieldDouble(ON_CELLS)
        f.setName("field")
        f.setMesh(m)
        f.setArray(DataArrayDouble(100))
        f.getArray()[:]=100.
        f.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldNoProfileSBT(f)
        # redirect stderr
        tmp=StdOutRedirect(errfname)
        self.assertRaises(InterpKernelException,f1ts.write,fname,0) # it should raise !
        del tmp
        #
        if os.path.exists(errfname):
            os.remove(errfname)
        #
        pass

    def testUnivStatus1(self):
        """ Non regression test to check the effectiveness of univ write status."""
        fname="Pyfile95.med"
        arr=DataArrayDouble(10) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr) ; m.setName("mesh")
        mm=MEDFileCMesh() ; mm.setMesh(m)
        mm.setUnivNameWrStatus(False) # test is here
        mm.write(fname,2)
        mm=MEDFileCMesh(fname)
        self.assertEqual(mm.getUnivName(),"")
        mm.setUnivNameWrStatus(True)
        mm.write(fname,2)
        mm=MEDFileCMesh(fname)
        self.assertTrue(mm.getUnivName()!="")
        pass
    
    def testEmptyMesh(self):
      """ MEDLoader should be able to consistently write and read an empty mesh (coords array
      with 0 tuples """
      fname = "Pyfile96.med" 
      m = MEDCouplingUMesh('toto', 2)
      m.setCoords(DataArrayDouble([], 0, 2))
      m.setConnectivity(DataArrayInt([]), DataArrayInt([0]))
      mfu = MEDFileUMesh()
      mfu.setMeshAtLevel(0, m)
      mfu.write(fname, 2)
      mfu2 = MEDFileUMesh(fname)
      self.assertEqual('toto', mfu2.getName())
      lvl = mfu2.getNonEmptyLevels()
      self.assertEqual((), lvl)

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def testMEDFileUMeshPickeling2(self):
      """ Check that pickalization can be performed on a unpickalized instance. Non regression test."""
      name="Mesh_1"
      grpName1="HAUT"
      grpName2="BASE"
      hauteur=1.
      nbOfNodesPerAxis=3
      arr=DataArrayDouble(nbOfNodesPerAxis) ; arr.iota() ; arr/=(nbOfNodesPerAxis-1) ; arr*=hauteur
      m=MEDCouplingCMesh() ; m.setCoords(arr,arr,arr) ; m=m.buildUnstructured() ; m.setName(name)
      mesh=MEDFileUMesh() ; mesh[0]=m
      m1=m.computeSkin() ; mesh[-1]=m1
      #
      bary1=m1.computeCellCenterOfMass()[:,2]
      grp1=bary1.findIdsInRange(hauteur-1e-12,hauteur+1e-12) ; grp1.setName(grpName1)
      grp2=bary1.findIdsInRange(0.-1e-12,0.+1e-12) ; grp2.setName(grpName2)
      mesh.setGroupsAtLevel(-1,[grp1,grp2])
      
      st=pickle.dumps(mesh,2)
      mm=pickle.loads(st)
      st2=pickle.dumps(mm,2)
      mm2=pickle.loads(st2)
      self.assertTrue(mesh.isEqual(mm2,1e-12)[0])
      pass

    def testMEDFileEquivalence1(self):
      """ First check of equivalence implementation in MEDFileMesh"""
      fileName="Pyfile97.med"
      meshName="M_01"
      mm=MEDFileUMesh()
      coo=DataArrayDouble([(0,0,0),(6,0,0),(19,0,0),(36,0,0),(0,4,0),(6,4,0),(19,4,0),(36,4,0),(0,13,0),(6,13,0),(19,13,0),(36,13,0),(0,24,0),(6,24,0),(19,24,0),(36,24,0),(0,0,6),(6,0,6),(19,0,6),(36,0,6),(0,4,6),(6,4,6),(19,4,6),(36,4,6),(0,13,6),(6,13,6),(19,13,6),(36,13,6),(0,24,6),(6,24,6),(19,24,6),(36,24,6),(6,0,3),(6,2,0),(12.5,0,0),(19,0,3),(19,2,0),(6,4,3),(12.5,4,0),(19,4,3),(6,2,6),(12.5,0,6),(19,2,6),(12.5,4,6),(6,2,3),(12.5,0,3),(12.5,2,0),(19,2,3),(12.5,4,3),(12.5,2,6),(12.5,2,3)])
      coo.setInfoOnComponents(["X [Sans_unite]","Y [Sans_unite]","Z [Sans_unite]"])
      connQ4=DataArrayInt([1,17,21,5,2,18,22,6,21,5,6,22,1,32,44,33,17,40,44,32,21,37,44,40,5,33,44,37,2,35,47,36,18,42,47,35,22,39,47,42,6,36,47,39,21,37,48,43,5,38,48,37,6,39,48,38,22,43,48,39])
      m1=MEDCoupling1SGTUMesh(meshName,NORM_QUAD4) ; m1.setCoords(coo) ; m1.setNodalConnectivity(connQ4) ; mm[-1]=m1
      connH8=DataArrayInt([20,16,17,21,4,0,1,5,22,18,19,23,6,2,3,7,24,20,21,25,8,4,5,9,25,21,22,26,9,5,6,10,26,22,23,27,10,6,7,11,28,24,25,29,12,8,9,13,29,25,26,30,13,9,10,14,30,26,27,31,14,10,11,15,21,40,49,43,37,44,50,48,40,17,41,49,44,32,45,50,49,41,18,42,50,45,35,47,43,49,42,22,48,50,47,39,44,32,45,50,33,1,34,46,37,44,50,48,5,33,46,38,48,50,47,39,38,46,36,6,50,45,35,47,46,34,2,36])
      m0=MEDCoupling1SGTUMesh(meshName,NORM_HEXA8) ; m0.setCoords(coo) ; m0.setNodalConnectivity(connH8) ; mm[0]=m0
      mm.getFamilyFieldAtLevel(-1)[:]=-2
      mm.getFamilyFieldAtLevel(0)[:]=0
      mm.addFamily("HOMARD________-1",-1)
      mm.addFamily("HOMARD________-2",-2)
      mm.addFamily("HOMARD________-3",-3)
      mm.setFamiliesIdsOnGroup("HOMARD",[-1,-2,-3])
      
      eqName="MAILLES_A_RECOLLER_APRES_HOMARD"
      descEq="Cette equivalence decrit les mailles a recoller. Dans chaque correspondance, le premier numero est celui de la maille coupee ; le second numero est celui d'une des petites mailles en regard."
      mm.initializeEquivalences()
      eqs=mm.getEquivalences()
      eq0=eqs.appendEmptyEquivalenceWithName(eqName)
      eq0.setDescription(descEq)
      corr=DataArrayInt([(0,3),(0,4),(0,5),(0,6),(1,7),(1,8),(1,9),(1,10),(2,11),(2,12),(2,13),(2,14)])
      eq0.setArray(-1,corr)
      self.assertEqual(eq0.getCell().size(),1)
      self.assertTrue(eq0.getCell().getArray(NORM_QUAD4).isEqual(corr))
      eq0.getCell().clear()
      self.assertEqual(eq0.getCell().size(),0)
      eq0.getCell().setArrayForType(NORM_QUAD4,corr)
      self.assertEqual(eq0.getCell().size(),1)
      self.assertTrue(eq0.getCell().getArray(NORM_QUAD4).isEqual(corr))
      mm.killEquivalences()
      mm.initializeEquivalences()
      eqs=mm.getEquivalences()
      eq0=eqs.appendEmptyEquivalenceWithName(eqName)
      eq0.setDescription(descEq)
      c=eq0.initCell()
      c.setArrayForType(NORM_QUAD4,corr)
      self.assertEqual(eq0.getCell().size(),1)
      self.assertTrue(eq0.getCell().getArray(NORM_QUAD4).isEqual(corr))
      mm2=mm.deepCopy()
      self.assertTrue(mm.isEqual(mm2,1e-12)[0])
      self.assertEqual(mm2.getEquivalences().size(),1)
      self.assertTrue(mm2.getEquivalences().getEquivalence(0).getCell().getArray(NORM_QUAD4).isEqual(corr))
      mm2.getEquivalences().getEquivalence(0).getCell().getArray(NORM_QUAD4)[0,0]=2
      self.assertTrue(not mm.isEqual(mm2,1e-12)[0])
      mm2.getEquivalences().getEquivalence(0).getCell().getArray(NORM_QUAD4)[0,0]=0
      self.assertTrue(mm.isEqual(mm2,1e-12)[0])
      mm.write(fileName,2)
      #
      mm3=MEDFileMesh.New(fileName)
      self.assertTrue(mm.isEqual(mm3,1e-12)[0])
      pass

    def testMEDFileForFamiliesPlayer1(self):
      """Non regression bug EDF11911. For serial killers using same family name to store both cells and nodes ! Only sky is the limit."""
      fileName="Pyfile98.med"
      meshName="mesh"
      magicSt="%s%%04i"%(MEDFileMesh.GetMagicFamilyStr())
      arr=DataArrayDouble(4) ; arr.iota()
      m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
      m=m.buildUnstructured()
      mm=MEDFileUMesh()
      mm[0]=m
      mm.setName(meshName)
      mm.setFamilyId("FAMILLE_ZERO",0)
      mm.getFamilyFieldAtLevel(0)[-3:]=-4
      mm.setFamilyId("RIDF%s"%(magicSt%0),-4)
      mm.setGroupsOnFamily("RIDF%s"%(magicSt%0),["RID"])
      d=DataArrayInt(16) ; d[:]=0 ; d[[1,2,4,5]]=3
      mm.setFamilyFieldArr(1,d)
      mm.setFamilyId("RIDF%s"%(magicSt%1),3)
      mm.setGroupsOnFamily("RIDF%s"%(magicSt%1),["RID"])
      self.assertEqual(mm.getFamiliesNames(),("FAMILLE_ZERO",'RIDF!/__\\!0000','RIDF!/__\\!0001'))
      self.assertEqual(mm.getFamiliesNamesWithFilePointOfView(),("FAMILLE_ZERO","RIDF","RIDF")) # <- the aim of test is here !
      self.assertEqual(mm.getFamiliesIdsOnGroup("RID"),(-4,3))
      mm.write(fileName,2)
      # now read such funny file !
      mm2=MEDFileMesh.New(fileName) # <- normaly mdump of Pyfile98.med must contain only RID and FAMILLE_ZERO families.
      self.assertTrue(mm.isEqual(mm2,1e-16))
      self.assertEqual(mm2.getFamiliesNames(),("FAMILLE_ZERO",'RIDF!/__\\!0000','RIDF!/__\\!0001'))
      self.assertEqual(mm2.getFamiliesNamesWithFilePointOfView(),("FAMILLE_ZERO","RIDF","RIDF"))
      self.assertEqual(mm2.getFamiliesIdsOnGroup("RID"),(-4,3))# <- very important too !
      pass

    def testCartesianizer1(self):
      """ This test is advanced to be sure that no unnecessary copies had been made during cartesianization process. """
      # UMesh non cart
      arr=DataArrayDouble(4) ; arr.iota() ; m=MEDCouplingCMesh() ; m.setCoords(arr,arr) ; m=m.buildUnstructured()
      mm=MEDFileUMesh() ; mm[0]=m ; mm.forceComputationOfParts()
      d0=DataArrayInt(16) ; d0[:]=0
      d1=DataArrayInt(9)  ; d1[:]=0
      mm.setFamilyFieldArr(0,d1) ; mm.setFamilyFieldArr(1,d0)
      mm.setName("a") ; mm.setDescription("b") ; mm.setTime(3,4,5.) ; mm.addFamily("c",-4) ; mm.setFamiliesOnGroup("d",["c"]) ; mm.setTimeUnit("ms")
      ref0=mm.getCoords().getHiddenCppPointer()
      ref1=mm[0].getNodalConnectivity().getHiddenCppPointer()
      self.assertEqual(ref0,mm[0].getCoords().getHiddenCppPointer())
      ref2=mm[0].getNodalConnectivityIndex().getHiddenCppPointer()
      ref3=mm.getDirectUndergroundSingleGeoTypeMesh(NORM_QUAD4).getNodalConnectivity().getHiddenCppPointer()
      self.assertEqual(ref0,mm.getDirectUndergroundSingleGeoTypeMesh(NORM_QUAD4).getCoords().getHiddenCppPointer())
      mm.setAxisType(AX_CYL) #<- important
      mm2=mm.cartesianize() # the trigger
      self.assertEqual(mm2.getAxisType(),AX_CART)
      mm.setAxisType(AX_CART) # this is here only to avoid complaints
      self.assertTrue(isinstance(mm2,MEDFileUMesh))
      self.assertTrue(mm.getHiddenCppPointer()!=mm2.getHiddenCppPointer())
      self.assertTrue(ref0==mm.getCoords().getHiddenCppPointer()) # <- here important
      self.assertTrue(ref0!=mm2.getCoords().getHiddenCppPointer()) # <- here important
      self.assertEqual(mm2.getCoords().getHiddenCppPointer(),mm2[0].getCoords().getHiddenCppPointer())
      self.assertEqual(mm2.getCoords().getHiddenCppPointer(),mm2.getDirectUndergroundSingleGeoTypeMesh(NORM_QUAD4).getCoords().getHiddenCppPointer())
      self.assertEqual(mm2[0].getNodalConnectivity().getHiddenCppPointer(),ref1) # <- here very important
      self.assertEqual(mm2[0].getNodalConnectivityIndex().getHiddenCppPointer(),ref2) # <- here very important
      self.assertEqual(mm2.getDirectUndergroundSingleGeoTypeMesh(NORM_QUAD4).getNodalConnectivity().getHiddenCppPointer(),ref3) # <- here very important
      self.assertEqual(mm2.getName(),mm.getName())
      self.assertEqual(mm2.getDescription(),mm.getDescription())
      self.assertEqual(mm2.getTime(),mm.getTime())
      self.assertEqual(mm2.getTime(),mm.getTime())
      self.assertEqual(mm2.getTimeUnit(),mm.getTimeUnit())
      self.assertEqual(mm2.getGroupsNames(),mm.getGroupsNames())
      self.assertEqual(mm2.getFamiliesNames(),mm.getFamiliesNames())
      self.assertEqual([mm2.getFamilyId(elt) for elt in mm2.getFamiliesNames()],[mm.getFamilyId(elt2) for elt2 in mm.getFamiliesNames()])
      self.assertEqual(mm.getFamilyFieldAtLevel(0).getHiddenCppPointer(),d1.getHiddenCppPointer())
      self.assertEqual(mm2.getFamilyFieldAtLevel(0).getHiddenCppPointer(),d1.getHiddenCppPointer()) # <- here very important
      self.assertEqual(mm.getFamilyFieldAtLevel(1).getHiddenCppPointer(),d0.getHiddenCppPointer())
      self.assertEqual(mm2.getFamilyFieldAtLevel(1).getHiddenCppPointer(),d0.getHiddenCppPointer()) # <- here very important
      # UMesh cart
      mm.setAxisType(AX_CART)
      mm2=mm.cartesianize() # the trigger
      self.assertEqual(mm2.getAxisType(),AX_CART)
      self.assertTrue(isinstance(mm2,MEDFileUMesh))
      self.assertTrue(mm.getHiddenCppPointer()==mm2.getHiddenCppPointer()) # optimization
      # CurveLinearMesh non cart
      arr=DataArrayDouble(4) ; arr.iota() ; m=MEDCouplingCMesh() ; m.setCoords(arr,arr) ; m=m.buildCurveLinear()
      mm=MEDFileCurveLinearMesh() ; mm.setMesh(m) ; mm.setAxisType(AX_CYL) #<- important
      mm.setFamilyFieldArr(0,d1) ; mm.setFamilyFieldArr(1,d0)
      mm.setName("a") ; mm.setDescription("b") ; mm.setTime(3,4,5.) ; mm.addFamily("c",-4) ; mm.setFamiliesOnGroup("d",["c"]) ; mm.setTimeUnit("ms")
      ref0=mm.getMesh().getCoords().getHiddenCppPointer()
      mm2=mm.cartesianize() # the trigger
      self.assertEqual(mm2.getAxisType(),AX_CART)
      self.assertTrue(isinstance(mm2,MEDFileCurveLinearMesh))
      self.assertTrue(mm.getHiddenCppPointer()!=mm2.getHiddenCppPointer())
      self.assertTrue(ref0==mm.getMesh().getCoords().getHiddenCppPointer()) # <- here important
      self.assertTrue(ref0!=mm2.getMesh().getCoords().getHiddenCppPointer()) # <- here important
      self.assertEqual(mm2.getMesh().getNodeGridStructure(),mm.getMesh().getNodeGridStructure())
      self.assertEqual(mm2.getName(),mm.getName())
      self.assertEqual(mm2.getDescription(),mm.getDescription())
      self.assertEqual(mm2.getTime(),mm.getTime())
      self.assertEqual(mm2.getTime(),mm.getTime())
      self.assertEqual(mm2.getTimeUnit(),mm.getTimeUnit())
      self.assertEqual(mm2.getGroupsNames(),mm.getGroupsNames())
      self.assertEqual(mm2.getFamiliesNames(),mm.getFamiliesNames())
      self.assertEqual([mm2.getFamilyId(elt) for elt in mm2.getFamiliesNames()],[mm.getFamilyId(elt2) for elt2 in mm.getFamiliesNames()])
      self.assertEqual(mm.getFamilyFieldAtLevel(0).getHiddenCppPointer(),d1.getHiddenCppPointer())
      self.assertEqual(mm2.getFamilyFieldAtLevel(0).getHiddenCppPointer(),d1.getHiddenCppPointer()) # <- here very important
      self.assertEqual(mm.getFamilyFieldAtLevel(1).getHiddenCppPointer(),d0.getHiddenCppPointer())
      self.assertEqual(mm2.getFamilyFieldAtLevel(1).getHiddenCppPointer(),d0.getHiddenCppPointer()) # <- here very important
      # CurveLinearMesh cart
      mm.setAxisType(AX_CART)
      mm2=mm.cartesianize() # the trigger
      self.assertEqual(mm2.getAxisType(),AX_CART)
      self.assertTrue(isinstance(mm2,MEDFileCurveLinearMesh))
      self.assertTrue(mm.getHiddenCppPointer()==mm2.getHiddenCppPointer()) # optimization
      # CMesh non cart
      arr=DataArrayDouble(4) ; arr.iota() ; m=MEDCouplingCMesh() ; m.setCoords(arr,arr)
      mm=MEDFileCMesh() ; mm.setMesh(m) ; mm.setAxisType(AX_CYL) #<- important
      mm.setFamilyFieldArr(0,d1) ; mm.setFamilyFieldArr(1,d0)
      mm.setName("a") ; mm.setDescription("b") ; mm.setTime(3,4,5.) ; mm.addFamily("c",-4) ; mm.setFamiliesOnGroup("d",["c"]) ; mm.setTimeUnit("ms")
      mm2=mm.cartesianize() # the trigger
      self.assertEqual(mm2.getAxisType(),AX_CART)
      self.assertTrue(isinstance(mm2,MEDFileCurveLinearMesh))
      self.assertEqual(mm2.getMesh().getNodeGridStructure(),mm.getMesh().getNodeGridStructure())
      self.assertEqual(mm2.getName(),mm.getName())
      self.assertEqual(mm2.getDescription(),mm.getDescription())
      self.assertEqual(mm2.getTime(),mm.getTime())
      self.assertEqual(mm2.getTime(),mm.getTime())
      self.assertEqual(mm2.getTimeUnit(),mm.getTimeUnit())
      self.assertEqual(mm2.getGroupsNames(),mm.getGroupsNames())
      self.assertEqual(mm2.getFamiliesNames(),mm.getFamiliesNames())
      self.assertEqual([mm2.getFamilyId(elt) for elt in mm2.getFamiliesNames()],[mm.getFamilyId(elt2) for elt2 in mm.getFamiliesNames()])
      self.assertEqual(mm.getFamilyFieldAtLevel(0).getHiddenCppPointer(),d1.getHiddenCppPointer())
      self.assertEqual(mm2.getFamilyFieldAtLevel(0).getHiddenCppPointer(),d1.getHiddenCppPointer()) # <- here very important
      self.assertEqual(mm.getFamilyFieldAtLevel(1).getHiddenCppPointer(),d0.getHiddenCppPointer())
      self.assertEqual(mm2.getFamilyFieldAtLevel(1).getHiddenCppPointer(),d0.getHiddenCppPointer()) # <- here very important
      # CMesh cart
      mm.setAxisType(AX_CART)
      mm2=mm.cartesianize() # the trigger
      self.assertEqual(mm2.getAxisType(),AX_CART)
      self.assertTrue(isinstance(mm2,MEDFileCMesh))
      self.assertTrue(mm.getHiddenCppPointer()==mm2.getHiddenCppPointer()) # optimization
      pass

    def testCheckCoherency(self):
      m2 = MEDCouplingUMesh("2d", 2)
      m2.setCoords(DataArrayDouble([(0.0, 1.0)] * 4, 4,2))  # whatever
      m2.setConnectivity(DataArrayInt([NORM_TRI3, 0,1,2,NORM_TRI3, 1,2,3]), DataArrayInt(([0,4,8])))
      m1 , _, _ , _, _ = m2.buildDescendingConnectivity()
      mum = MEDFileUMesh()
      mum.setMeshAtLevel(0, m2)
      mum.setMeshAtLevel(-1, m1)
      mum.checkConsistency()
      mum2 = mum.deepCopy()

      # Nodes
      arr = DataArrayInt([2]*4)
      mum.setFamilyFieldArr(1, arr); arr.reAlloc(35);
      self.assertRaises(InterpKernelException, mum.checkConsistency)
      mum=mum2; mum2=mum.deepCopy();
      arr = DataArrayInt([2]*4)
      mum.setRenumFieldArr(1, arr); arr.reAlloc(35);
      self.assertRaises(InterpKernelException, mum.checkConsistency)
      mum=mum2; mum2=mum.deepCopy();
      mum.setRenumFieldArr(1, DataArrayInt([2]*4))
      self.assertRaises(InterpKernelException, mum.checkConsistency)
      mum=mum2; mum2=mum.deepCopy();
      arr = DataArrayAsciiChar(['tutu           x']*4)
      mum.setNameFieldAtLevel(1, arr); arr.reAlloc(35);
      self.assertRaises(InterpKernelException, mum.checkConsistency)

      # 2D
      mum=mum2; mum2=mum.deepCopy();
      arr = DataArrayInt([2]*2)
      mum.setFamilyFieldArr(0, arr); arr.reAlloc(35);
      self.assertRaises(InterpKernelException, mum.checkConsistency)
      mum=mum2; mum2=mum.deepCopy();
      arr = DataArrayInt([2]*2)
      mum.setRenumFieldArr(0, arr); arr.reAlloc(35);
      self.assertRaises(InterpKernelException, mum.checkConsistency)
      mum=mum2; mum2=mum.deepCopy();
      mum.setRenumFieldArr(0, DataArrayInt([2]*2))
      self.assertRaises(InterpKernelException, mum.checkConsistency)
      mum=mum2; mum2=mum.deepCopy();
      arr = DataArrayAsciiChar(['tutu           x']*2)
      mum.setNameFieldAtLevel(0, arr); arr.reAlloc(35);
      self.assertRaises(InterpKernelException, mum.checkConsistency)

      # 1D
      mum=mum2; mum2=mum.deepCopy();
      arr = DataArrayInt([2]*5)
      mum.setFamilyFieldArr(-1, arr); arr.reAlloc(35);
      self.assertRaises(InterpKernelException, mum.checkConsistency)
      mum=mum2; mum2=mum.deepCopy();
      arr = DataArrayInt([2]*5)
      mum.setRenumFieldArr(-1, arr); arr.reAlloc(35);
      self.assertRaises(InterpKernelException, mum.checkConsistency)
      mum=mum2; mum2=mum.deepCopy();
      mum.setRenumFieldArr(-1, DataArrayInt([2]*5))
      self.assertRaises(InterpKernelException, mum.checkConsistency)
      mum=mum2; mum2=mum.deepCopy();
      arr = DataArrayAsciiChar(['tutu           x']*5)
      mum.setNameFieldAtLevel(-1, arr); arr.reAlloc(35);
      self.assertRaises(InterpKernelException, mum.checkConsistency)

    def testCheckSMESHConsistency(self):
      m2 = MEDCouplingUMesh("2d", 2)
      m2.setCoords(DataArrayDouble([(0.0, 1.0)] * 4, 4,2))  # whatever
      m2.setConnectivity(DataArrayInt([NORM_TRI3, 0,1,2,NORM_TRI3, 1,2,3]), DataArrayInt(([0,4,8])))
      m1 , _, _ , _, _ = m2.buildDescendingConnectivity()
      mum = MEDFileUMesh()
      mum.setMeshAtLevel(0, m2)
      mum.setMeshAtLevel(-1, m1)
      mum.checkConsistency()
      mum.checkSMESHConsistency()
      n2 = DataArrayInt(m2.getNumberOfCells(), 1); n2.iota(1)
      n1 = DataArrayInt(m1.getNumberOfCells(), 1); n1.iota(1)
      mum.setRenumFieldArr(0, n2)
      mum.setRenumFieldArr(-1, n1)
      self.assertRaises(InterpKernelException, mum.checkSMESHConsistency)
      mum.setRenumFieldArr(-1, n1+100)
      mum.checkSMESHConsistency()
      pass

    def testClearNodeAndCellNumbers(self):
      m2 = MEDCouplingUMesh("2d", 2)
      m2.setCoords(DataArrayDouble([(0.0, 1.0)] * 4, 4,2))  # whatever
      m2.setConnectivity(DataArrayInt([NORM_TRI3, 0,1,2,NORM_TRI3, 1,2,3]), DataArrayInt(([0,4,8])))
      m1 , _, _ , _, _ = m2.buildDescendingConnectivity()
      mum = MEDFileUMesh()
      mum.setMeshAtLevel(0, m2)
      mum.setMeshAtLevel(-1, m1)
      mum.checkConsistency()
      n2 = DataArrayInt(m2.getNumberOfCells(), 1); n2.iota(1)
      n1 = DataArrayInt(m1.getNumberOfCells(), 1); n1.iota(1)
      mum.setRenumFieldArr(0, n2)
      mum.setRenumFieldArr(-1, n1)
      mum.clearNodeAndCellNumbers()
      mum.checkSMESHConsistency()
      pass

    def testCMeshSetFamilyFieldArrNull(self):
      meshName="mesh"
      fname="Pyfile99.med"
      arrX=DataArrayDouble([0,1,2,3])
      arrY=DataArrayDouble([0,1,2])
      m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY) ; m.setName(meshName)
      mm=MEDFileCMesh() ; mm.setMesh(m)
      famCellIds=DataArrayInt([0,-2,-2,-1,-2,0])
      famNodeIds=DataArrayInt([0,0,0,3,4,1,2,7,2,1,0,0])
      mm.setFamilyFieldArr(0,famCellIds)
      mm.setFamilyFieldArr(1,famNodeIds)
      mm.write(fname,2)
      mm=MEDFileMesh.New(fname)
      self.assertTrue(mm.getFamilyFieldAtLevel(0) is not None)
      self.assertTrue(mm.getFamilyFieldAtLevel(1) is not None)
      mm.setFamilyFieldArr(0,None)#<- bug was here
      mm.setFamilyFieldArr(1,None)#<- bug was here
      self.assertTrue(mm.getFamilyFieldAtLevel(0) is None)
      self.assertTrue(mm.getFamilyFieldAtLevel(1) is None)
      mm3=mm.deepCopy()
      self.assertTrue(mm3.getFamilyFieldAtLevel(0) is None)
      self.assertTrue(mm3.getFamilyFieldAtLevel(1) is None)
      mm.write(fname,2)
      mm2=MEDFileMesh.New(fname)
      self.assertTrue(mm2.getFamilyFieldAtLevel(0) is None)
      self.assertTrue(mm2.getFamilyFieldAtLevel(1) is None)
      pass

    def testAppendFieldProfileOnIntField(self):
      fname="Pyfile100.med"
      arrX=DataArrayDouble([0,1,2,3])
      arrY=DataArrayDouble([0,1,2])
      mesh=MEDCouplingCMesh() ; mesh.setCoords(arrX,arrY) ; mesh.setName("Mesh")
      mm=MEDFileCMesh()
      mm.setMesh(mesh)
      #
      fmts=MEDFileIntFieldMultiTS()
      pflName="PFL"
      pfl=DataArrayInt([1,3,5]) ; pfl.setName(pflName)
      f=MEDCouplingFieldInt(ON_CELLS) ; f.setMesh(mesh)
      fieldName="FieldOnCell"
      f.setTime(1.2,1,1) ; f.setName(fieldName)
      arr=DataArrayInt([101,102,103]) ; f.setArray(arr)
      fmts.appendFieldProfile(f,mm,0,pfl)
      #
      mm.write(fname,2)
      fmts.write(fname,0)
      #
      mm=MEDFileMesh.New(fname)
      fmts=MEDFileAnyTypeFieldMultiTS.New(fname)
      self.assertTrue(isinstance(fmts,MEDFileIntFieldMultiTS))
      self.assertEqual(fmts.getName(),fieldName)
      self.assertEqual(len(fmts),1)
      f1ts=fmts[0]
      ftest,pfltest=f1ts.getFieldWithProfile(ON_CELLS,0,mm)
      self.assertEqual(pfltest.getName(),pflName)
      self.assertEqual(ftest.getName(),fieldName)
      self.assertTrue(ftest.isEqualWithoutConsideringStr(arr))
      ftest2=f1ts.getFieldOnMeshAtLevel(ON_CELLS,0,mm)
      self.assertTrue(ftest2.getArray().isEqualWithoutConsideringStr(arr))
      self.assertEqual(ftest2.getTime(),f.getTime())
      self.assertEqual(ftest2.getMesh().getNumberOfCells(),len(arr))
      pass

    def testMEDFileFieldEasyField1(self):
      """Check for all spatial discretization of field (cells,nodes,elno,gauss) for double field that all is OK. Here no profile and only top level is considered."""
      ## Basic test on cells on top level
      fname="Pyfile101.med"
      fieldName="field1"
      mm=MEDFileUMesh()
      coo=DataArrayDouble([(3,2,1),(8,7,6),(5,9,10)])
      m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo)
      m.allocateCells()
      m.insertNextCell(NORM_TRI3,[0,1,2])
      m.insertNextCell(NORM_TRI3,[3,4,5])
      m.insertNextCell(NORM_TRI3,[6,7,8])
      m.insertNextCell(NORM_TRI3,[9,10,11])
      m.insertNextCell(NORM_QUAD4,[100,101,102,103])
      m.insertNextCell(NORM_QUAD4,[104,105,106,107])
      mm[0]=m
      mm.write(fname,2)
      arr0=DataArrayDouble([10,11,12,13,100,101])
      f=MEDCouplingFieldDouble(ON_CELLS) ; f.setArray(arr0) ; f.setMesh(m)
      f.setName(fieldName) ; f.setTime(2.,6,7)
      f0=f.deepCopy()
      ff=MEDFileFieldMultiTS() ; ff.appendFieldNoProfileSBT(f)
      ff.write(fname,0)
      arr2=arr0+1000 ; f.setArray(arr2)
      f.setTime(3.,8,9) ; ff=MEDFileField1TS() ; ff.setFieldNoProfileSBT(f)
      ff.write(fname,0)
      f1=f.deepCopy()
      ##
      mm=MEDFileMesh.New(fname)
      f1ts=MEDFileField1TS(fname,fieldName,6,7)
      ftst0=f1ts.field(mm)
      self.assertTrue(f0.isEqual(ftst0,1e-12,1e-12))
      f1ts=MEDFileField1TS(fname,fieldName,8,9)
      ftst1=f1ts.field(mm)
      self.assertTrue(f1.isEqual(ftst1,1e-12,1e-12))
      fmts=MEDFileFieldMultiTS(fname,fieldName)
      self.assertTrue(f1.isEqual(fmts.field(8,9,mm),1e-12,1e-12))
      ## Basic test on nodes on top level
      f2=MEDCouplingFieldDouble(ON_NODES) ; arr2=DataArrayDouble([200,201,202]) ; arr2.setInfoOnComponent(0,"tutu") ; f2.setArray(arr2) ; f2.setMesh(m) ; f2.setTime(22.,23,24)
      f2.setName(fieldName)
      mm.write(fname,2)
      ff=MEDFileField1TS() ; ff.setFieldNoProfileSBT(f2) ; ff.write(fname,0)
      #
      mm=MEDFileMesh.New(fname)
      f1ts=MEDFileField1TS(fname,fieldName,23,24)
      self.assertTrue(f2.isEqual(f1ts.field(mm),1e-12,1e-12))
      fmts=MEDFileFieldMultiTS(fname,fieldName)
      self.assertTrue(f2.isEqual(fmts.field(23,24,mm),1e-12,1e-12))
      ## Node on elements
      f3=MEDCouplingFieldDouble(ON_GAUSS_NE) ; f3.setMesh(m) ; arr3=DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]) ; f3.setArray(arr3) ; f3.setTime(0.5,2,3)
      f3.setName(fieldName) ; f3.checkConsistencyLight()
      mm.write(fname,2) ; ff=MEDFileField1TS() ; ff.setFieldNoProfileSBT(f3) ; ff.write(fname,0)
      #
      mm=MEDFileMesh.New(fname)
      f1ts=MEDFileField1TS(fname,fieldName,2,3)
      self.assertTrue(f3.isEqual(f1ts.field(mm),1e-12,1e-12))
      ## Gauss
      f4=MEDCouplingFieldDouble(ON_GAUSS_PT) ; f4.setMesh(m) ; f4.setName(fieldName)
      f4.setGaussLocalizationOnType(NORM_TRI3,[0.,0.,1.,0.,1.,1.],[0.1,0.1, 0.2,0.2, 0.3,0.3, 0.4,0.4, 0.5,0.5],[0.2,0.3,0.1,0.05,0.35])
      f4.setGaussLocalizationOnType(NORM_QUAD4,[0.,0.,1.,0.,1.,1.,0.,1.],[0.3,0.4, 0.6,0.7],[0.7,0.3]) ; f4.setTime(0.25,4,5)
      arr4=DataArrayDouble([0,1,2,3,4 ,10,11,12,13,14, 20,21,22,23,24, 30,31,32,33,34, 45,46, 55,56]) ; arr4.setInfoOnComponent(0,"abc") ; f4.setArray(arr4)
      f4.checkConsistencyLight()
      mm.write(fname,2) ; ff=MEDFileField1TS() ; ff.setFieldNoProfileSBT(f4) ; ff.write(fname,0)
      #
      mm=MEDFileMesh.New(fname)
      f1ts=MEDFileField1TS(fname,fieldName,4,5)
      self.assertTrue(f4.isEqual(f1ts.field(mm),1e-12,1e-12))
      pass
  
    def testMEDFileFieldEasyField2(self):
        """Same thantestMEDFileFieldEasyField1 except that here intfields are considered.
        Check for all spatial discretization of field (cells,nodes,elno,gauss) for int field that all is OK. Here no profile and only top level is considered."""
        ## Basic test on cells on top level
        fname="Pyfile102.med"
        fieldName="field1"
        mm=MEDFileUMesh()
        coo=DataArrayDouble([(3,2,1),(8,7,6),(5,9,10)])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_TRI3,[0,1,2])
        m.insertNextCell(NORM_TRI3,[3,4,5])
        m.insertNextCell(NORM_TRI3,[6,7,8])
        m.insertNextCell(NORM_TRI3,[9,10,11])
        m.insertNextCell(NORM_QUAD4,[100,101,102,103])
        m.insertNextCell(NORM_QUAD4,[104,105,106,107])
        mm[0]=m
        mm.write(fname,2)
        arr0=DataArrayInt([10,11,12,13,100,101])
        f=MEDCouplingFieldInt(ON_CELLS) ; f.setArray(arr0) ; f.setMesh(m)
        f.setName(fieldName) ; f.setTime(2.,6,7)
        f0=f.deepCopy()
        ff=MEDFileIntFieldMultiTS() ; ff.appendFieldNoProfileSBT(f)
        ff.write(fname,0)
        arr2=arr0+1000 ; f.setArray(arr2)
        f.setTime(3.,8,9) ; ff=MEDFileIntField1TS() ; ff.setFieldNoProfileSBT(f)
        ff.write(fname,0)
        f1=f.deepCopy()
        ##
        mm=MEDFileMesh.New(fname)
        f1ts=MEDFileIntField1TS(fname,fieldName,6,7)
        ftst0=f1ts.field(mm)
        self.assertTrue(f0.isEqual(ftst0,1e-12,0))
        f1ts=MEDFileIntField1TS(fname,fieldName,8,9)
        ftst1=f1ts.field(mm)
        self.assertTrue(f1.isEqual(ftst1,1e-12,0))
        fmts=MEDFileIntFieldMultiTS(fname,fieldName)
        self.assertTrue(f1.isEqual(fmts.field(8,9,mm),1e-12,0))
        ## Basic test on nodes on top level
        f2=MEDCouplingFieldInt(ON_NODES) ; arr2=DataArrayInt([200,201,202]) ; arr2.setInfoOnComponent(0,"tutu") ; f2.setArray(arr2) ; f2.setMesh(m) ; f2.setTime(22.,23,24)
        f2.setName(fieldName)
        mm.write(fname,2)
        ff=MEDFileIntField1TS() ; ff.setFieldNoProfileSBT(f2) ; ff.write(fname,0)
        #
        mm=MEDFileMesh.New(fname)
        f1ts=MEDFileIntField1TS(fname,fieldName,23,24)
        self.assertTrue(f2.isEqual(f1ts.field(mm),1e-12,0))
        fmts=MEDFileIntFieldMultiTS(fname,fieldName)
        self.assertTrue(f2.isEqual(fmts.field(23,24,mm),1e-12,0))
        ## Node on elements
        f3=MEDCouplingFieldInt(ON_GAUSS_NE) ; f3.setMesh(m) ; arr3=DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]) ; f3.setArray(arr3) ; f3.setTime(0.5,2,3)
        f3.setName(fieldName) ; f3.checkConsistencyLight()
        mm.write(fname,2) ; ff=MEDFileIntField1TS() ; ff.setFieldNoProfileSBT(f3) ; ff.write(fname,0)
        #
        mm=MEDFileMesh.New(fname)
        f1ts=MEDFileIntField1TS(fname,fieldName,2,3)
        self.assertTrue(f3.isEqual(f1ts.field(mm),1e-12,0))
        ## Gauss
        f4=MEDCouplingFieldInt(ON_GAUSS_PT) ; f4.setMesh(m) ; f4.setName(fieldName)
        f4.setGaussLocalizationOnType(NORM_TRI3,[0.,0.,1.,0.,1.,1.],[0.1,0.1, 0.2,0.2, 0.3,0.3, 0.4,0.4, 0.5,0.5],[0.2,0.3,0.1,0.05,0.35])
        f4.setGaussLocalizationOnType(NORM_QUAD4,[0.,0.,1.,0.,1.,1.,0.,1.],[0.3,0.4, 0.6,0.7],[0.7,0.3]) ; f4.setTime(0.25,4,5)
        arr4=DataArrayInt([0,1,2,3,4 ,10,11,12,13,14, 20,21,22,23,24, 30,31,32,33,34, 45,46, 55,56]) ; arr4.setInfoOnComponent(0,"abc") ; f4.setArray(arr4)
        f4.checkConsistencyLight()
        mm.write(fname,2) ; ff=MEDFileIntField1TS() ; ff.setFieldNoProfileSBT(f4) ; ff.write(fname,0)
        #
        mm=MEDFileMesh.New(fname)
        f1ts=MEDFileIntField1TS(fname,fieldName,4,5)
        self.assertTrue(f4.isEqual(f1ts.field(mm),1e-12,0))
        pass

    def testMEDFileFieldEasyField3(self):
        """Here a multi level mesh. And field on cells lying on different level of this mesh. Show how "field" method deal with that. Here on field double are considered."""
        fname="Pyfile103.med"
        fieldName="field1"
        mm=MEDFileUMesh()
        coo=DataArrayDouble([(3,2,1),(8,7,6),(5,9,10)])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_TRI3,[0,1,2])
        m.insertNextCell(NORM_TRI3,[3,4,5])
        m.insertNextCell(NORM_TRI3,[6,7,8])
        m.insertNextCell(NORM_TRI3,[9,10,11])
        m.insertNextCell(NORM_QUAD4,[100,101,102,103])
        m.insertNextCell(NORM_QUAD4,[104,105,106,107])
        mm[-1]=m
        m0=MEDCouplingUMesh("mesh",3) ; m0.setCoords(coo)
        m0.allocateCells()
        m0.insertNextCell(NORM_TETRA4,[3,2,5,0])
        m0.insertNextCell(NORM_TETRA4,[7,6,3,2])
        mm[0]=m0
        mm.write(fname,2)
        # start slowly
        f1=MEDCouplingFieldDouble(ON_CELLS) ; f1.setName(fieldName) ; f1.setArray(DataArrayDouble([(0,100),(1,101)])) ; f1.setMesh(mm[0]) ; f1.setTime(4.,1,2)
        f1ts=MEDFileField1TS() ; f1ts.setFieldNoProfileSBT(f1) ; f1ts.write(fname,0)
        #
        mm=MEDFileMesh.New(fname) ; f1ts=MEDFileField1TS(fname,fieldName,1,2)
        self.assertTrue(f1.isEqual(f1ts.field(mm),1e-12,1e-12))
        # here f1 lying on level -1 not 0 check if "field" method detect it !
        f1=MEDCouplingFieldDouble(ON_CELLS) ; f1.setName(fieldName) ; f1.setArray(DataArrayDouble([(0,100),(1,101),(0,100),(1,101),(0,100),(1,101)]))
        f1.setMesh(mm[-1]) # -1 is very important
        f1.setTime(16.,3,4)
        f1.checkConsistencyLight()
        mm.write(fname,2)
        f1ts=MEDFileField1TS() ; f1ts.setFieldNoProfileSBT(f1) ; f1ts.write(fname,0)
        #
        mm=MEDFileMesh.New(fname) ; f1ts=MEDFileField1TS(fname,fieldName,3,4)
        self.assertTrue(f1.isEqual(f1ts.field(mm),1e-12,1e-12))
        # nodes on elements
        f3=MEDCouplingFieldDouble(ON_GAUSS_NE)
        f3.setMesh(mm[-1]) # this line is important
        arr3=DataArrayDouble([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]) ; f3.setArray(arr3) ; f3.setTime(0.5,2,3)
        f3.setName(fieldName) ; f3.checkConsistencyLight()
        mm.write(fname,2) ; ff=MEDFileField1TS() ; ff.setFieldNoProfileSBT(f3) ; ff.write(fname,0)
        #
        mm=MEDFileMesh.New(fname) ; f1ts=MEDFileField1TS(fname,fieldName,2,3)
        self.assertTrue(f3.isEqual(f1ts.field(mm),1e-12,1e-12))
        # gauss
        f4=MEDCouplingFieldDouble(ON_GAUSS_PT)
        f4.setMesh(mm[-1]) # this line is important
        f4.setName(fieldName)
        f4.setGaussLocalizationOnType(NORM_TRI3,[0.,0.,1.,0.,1.,1.],[0.1,0.1, 0.2,0.2, 0.3,0.3, 0.4,0.4, 0.5,0.5],[0.2,0.3,0.1,0.05,0.35])
        f4.setGaussLocalizationOnType(NORM_QUAD4,[0.,0.,1.,0.,1.,1.,0.,1.],[0.3,0.4, 0.6,0.7],[0.7,0.3]) ; f4.setTime(0.25,4,5)
        arr4=DataArrayDouble([0,1,2,3,4 ,10,11,12,13,14, 20,21,22,23,24, 30,31,32,33,34, 45,46, 55,56]) ; arr4.setInfoOnComponent(0,"abc") ; f4.setArray(arr4)
        f4.checkConsistencyLight()
        mm.write(fname,2) ; ff=MEDFileField1TS() ; ff.setFieldNoProfileSBT(f4) ; ff.write(fname,0)
        mm=MEDFileMesh.New(fname) ; f1ts=MEDFileField1TS(fname,fieldName,4,5)
        self.assertTrue(f4.isEqual(f1ts.field(mm),1e-12,1e-12))
        pass

    def testMEDFileFieldEasyField4(self):
        """ Same than testMEDFileFieldEasyField3 but with integers"""
        fname="Pyfile104.med"
        fieldName="field1"
        mm=MEDFileUMesh()
        coo=DataArrayDouble([(3,2,1),(8,7,6),(5,9,10)])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_TRI3,[0,1,2])
        m.insertNextCell(NORM_TRI3,[3,4,5])
        m.insertNextCell(NORM_TRI3,[6,7,8])
        m.insertNextCell(NORM_TRI3,[9,10,11])
        m.insertNextCell(NORM_QUAD4,[100,101,102,103])
        m.insertNextCell(NORM_QUAD4,[104,105,106,107])
        mm[-1]=m
        m0=MEDCouplingUMesh("mesh",3) ; m0.setCoords(coo)
        m0.allocateCells()
        m0.insertNextCell(NORM_TETRA4,[3,2,5,0])
        m0.insertNextCell(NORM_TETRA4,[7,6,3,2])
        mm[0]=m0
        mm.write(fname,2)
        # start slowly
        f1=MEDCouplingFieldInt(ON_CELLS) ; f1.setName(fieldName) ; f1.setArray(DataArrayInt([(0,100),(1,101)])) ; f1.setMesh(mm[0]) ; f1.setTime(4.,1,2)
        f1ts=MEDFileIntField1TS() ; f1ts.setFieldNoProfileSBT(f1) ; f1ts.write(fname,0)
        #
        mm=MEDFileMesh.New(fname) ; f1ts=MEDFileIntField1TS(fname,fieldName,1,2)
        self.assertTrue(f1.isEqual(f1ts.field(mm),1e-12,0))
        # here f1 lying on level -1 not 0 check if "field" method detect it !
        f1=MEDCouplingFieldInt(ON_CELLS) ; f1.setName(fieldName) ; f1.setArray(DataArrayInt([(0,100),(1,101),(0,100),(1,101),(0,100),(1,101)]))
        f1.setMesh(mm[-1]) # -1 is very important
        f1.setTime(16.,3,4)
        f1.checkConsistencyLight()
        mm.write(fname,2)
        f1ts=MEDFileIntField1TS() ; f1ts.setFieldNoProfileSBT(f1) ; f1ts.write(fname,0)
        #
        mm=MEDFileMesh.New(fname) ; f1ts=MEDFileIntField1TS(fname,fieldName,3,4)
        self.assertTrue(f1.isEqual(f1ts.field(mm),1e-12,0))
        # nodes on elements
        f3=MEDCouplingFieldInt(ON_GAUSS_NE)
        f3.setMesh(mm[-1]) # this line is important
        arr3=DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]) ; f3.setArray(arr3) ; f3.setTime(0.5,2,3)
        f3.setName(fieldName) ; f3.checkConsistencyLight()
        mm.write(fname,2) ; ff=MEDFileIntField1TS() ; ff.setFieldNoProfileSBT(f3) ; ff.write(fname,0)
        #
        mm=MEDFileMesh.New(fname) ; f1ts=MEDFileIntField1TS(fname,fieldName,2,3)
        self.assertTrue(f3.isEqual(f1ts.field(mm),1e-12,0))
        # gauss
        f4=MEDCouplingFieldInt(ON_GAUSS_PT)
        f4.setMesh(mm[-1]) # this line is important
        f4.setName(fieldName)
        f4.setGaussLocalizationOnType(NORM_TRI3,[0.,0.,1.,0.,1.,1.],[0.1,0.1, 0.2,0.2, 0.3,0.3, 0.4,0.4, 0.5,0.5],[0.2,0.3,0.1,0.05,0.35])
        f4.setGaussLocalizationOnType(NORM_QUAD4,[0.,0.,1.,0.,1.,1.,0.,1.],[0.3,0.4, 0.6,0.7],[0.7,0.3]) ; f4.setTime(0.25,4,5)
        arr4=DataArrayInt([0,1,2,3,4 ,10,11,12,13,14, 20,21,22,23,24, 30,31,32,33,34, 45,46, 55,56]) ; arr4.setInfoOnComponent(0,"abc") ; f4.setArray(arr4)
        f4.checkConsistencyLight()
        mm.write(fname,2) ; ff=MEDFileIntField1TS() ; ff.setFieldNoProfileSBT(f4) ; ff.write(fname,0)
        mm=MEDFileMesh.New(fname) ; f1ts=MEDFileIntField1TS(fname,fieldName,4,5)
        self.assertTrue(f4.isEqual(f1ts.field(mm),1e-12,0))
        pass

    def testMEDFileFieldEasyField5(self):
        """More and more difficult now look at how profiles are managed by "field" method."""
        fname="Pyfile105.med"
        fieldName="field1"
        mm=MEDFileUMesh()
        coo=DataArrayDouble([(3,2,1),(8,7,6),(5,9,10)])
        m=MEDCouplingUMesh("mesh",2) ; m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_TRI3,[0,1,2])
        m.insertNextCell(NORM_TRI3,[3,4,5])
        m.insertNextCell(NORM_TRI3,[6,7,8])
        m.insertNextCell(NORM_TRI3,[9,10,11])
        m.insertNextCell(NORM_QUAD4,[100,101,102,103])
        m.insertNextCell(NORM_QUAD4,[104,105,106,107])
        mm[0]=m
        mm.write(fname,2)
        pfl=DataArrayInt([0,2,3,5]) ; pfl.setName("pfl")
        m2=m.deepCopy()[pfl] ; m2.setName(m.getName())
        #
        arr0=DataArrayDouble([10,11,12,13])
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setArray(arr0) ; f.setMesh(m2)
        f.setName(fieldName) ; f.setTime(2.,6,7) ; f.checkConsistencyLight()
        ff=MEDFileFieldMultiTS() ; ff.appendFieldProfile(f,mm,0,pfl) # ff is a field on profile
        ff.write(fname,0)
        #
        mm=MEDFileMesh.New(fname) ; f1ts=MEDFileField1TS(fname,fieldName,6,7)
        self.assertTrue(f.isEqual(f1ts.field(mm),1e-12,1e-12))
        # more complicated -> multi level
        m0=MEDCouplingUMesh("mesh",3) ; m0.setCoords(coo)
        m0.allocateCells()
        m0.insertNextCell(NORM_TETRA4,[3,2,5,0])
        m0.insertNextCell(NORM_TETRA4,[7,6,3,2])
        mm2=MEDFileUMesh()
        mm2[0]=m0 ; mm2[-1]=m
        #
        ff=MEDFileField1TS() ; ff.setFieldProfile(f,mm2,-1,pfl)
        #
        mm=MEDFileMesh.New(fname) ; f1ts=MEDFileField1TS(fname,fieldName,6,7)
        self.assertTrue(f.isEqual(f1ts.field(mm),1e-12,1e-12))
        pass

    def testExtractPart1(self):
        coo=DataArrayDouble([(0,0),(1,0),(2,0),(3,0),(4,0),(0,1),(1,1),(2,1),(3,1),(4,1),(0,2),(1,2),(2,2),(3,2),(4,2)])
        meshName="mesh"
        m0=MEDCouplingUMesh(meshName,2) ; m0.setCoords(coo) ; m0.allocateCells()
        m0.insertNextCell(NORM_TRI3,[8,4,3])
        m0.insertNextCell(NORM_TRI3,[8,9,4])
        m0.insertNextCell(NORM_TRI3,[7,13,8])
        m0.insertNextCell(NORM_TRI3,[7,12,13])
        m0.insertNextCell(NORM_TRI3,[0,6,1])
        m0.insertNextCell(NORM_TRI3,[0,5,6])
        m0.insertNextCell(NORM_QUAD4,[1,6,7,2])
        m0.insertNextCell(NORM_QUAD4,[2,7,8,3])
        m0.insertNextCell(NORM_QUAD4,[8,13,14,9])
        m0.insertNextCell(NORM_QUAD4,[6,11,12,7])
        m0.insertNextCell(NORM_QUAD4,[5,10,11,6])
        #
        m1=MEDCouplingUMesh(meshName,1) ; m1.setCoords(coo) ; m1.allocateCells()
        m1.insertNextCell(NORM_SEG2,[10,5])
        m1.insertNextCell(NORM_SEG2,[5,0])
        m1.insertNextCell(NORM_SEG2,[0,1])
        m1.insertNextCell(NORM_SEG2,[1,2])
        m1.insertNextCell(NORM_SEG2,[2,3])
        m1.insertNextCell(NORM_SEG2,[3,4])
        m1.insertNextCell(NORM_SEG2,[4,9])
        m1.insertNextCell(NORM_SEG2,[9,14])
        m1.insertNextCell(NORM_SEG2,[14,13])
        m1.insertNextCell(NORM_SEG2,[13,12])
        m1.insertNextCell(NORM_SEG2,[12,11])
        m1.insertNextCell(NORM_SEG2,[11,10])
        mm=MEDFileUMesh()
        mm[0]=m0 ; mm[-1]=m1
        arr0=DataArrayInt([0,1,2,3,4,6,7,8,12,13])
        tab={} #
        tab[0]=DataArrayInt([0,2,3,4,6,7])
        tab[-1]=DataArrayInt([2,3,4,5,9])
        fs=MEDFileFields()
        self.assertTrue(mm.deduceNodeSubPartFromCellSubPart(tab).isEqual(arr0))
        tab[1]=arr0
        #
        fname0="Field0"
        fmts=MEDFileFieldMultiTS() ; fs.pushField(fmts)
        t0=(16.5,3,4)
        ic=["toto [m]"]
        arr0_0=DataArrayDouble([100,101,102,103,104,105,106,107,108,109,110]) ; arr0_0.setInfoOnComponents(ic)
        f0=MEDCouplingFieldDouble(ON_CELLS) ; f0.setTime(*t0) ; f0.setArray(arr0_0)
        f0.setMesh(m0) ; f0.setName(fname0)
        f1=MEDCouplingFieldDouble(ON_CELLS) ; f1.setTime(*t0) ; f1.setArray(DataArrayDouble([200,201,202,203,204,205,206,207,208,209,210,211]))
        f1.setMesh(m1) ; f1.setName(fname0) ; f1.getArray().setInfoOnComponents(ic)
        f2=MEDCouplingFieldDouble(ON_NODES) ; f2.setTime(*t0) ; f2.setArray(DataArrayDouble([300,301,302,303,304,305,306,307,308,309,310,311,312,313,314]))
        f2.setMesh(m0) ; f2.setName(fname0) ; f2.getArray().setInfoOnComponents(ic)
        f1ts=MEDFileField1TS() ; f1ts.setFieldNoProfileSBT(f0) ; f1ts.setFieldNoProfileSBT(f1) ; f1ts.setFieldNoProfileSBT(f2)
        fmts.pushBackTimeStep(f1ts)
        #
        mmOut=mm.extractPart(tab)
        #
        fsPart0=fs.extractPart(tab,mm)
        self.assertEqual(len(fsPart0),1)
        fmtsP=fsPart0[0]
        self.assertEqual(len(fmtsP),1)
        f1ts=fmtsP[0]
        self.assertRaises(InterpKernelException,f1ts.field,mmOut)
        #
        self.assertTrue(mmOut[0].computeCellCenterOfMass().isEqual(m0[tab[0]].computeCellCenterOfMass(),1e-12))
        self.assertTrue(mmOut[-1].computeCellCenterOfMass().isEqual(m1[tab[-1]].computeCellCenterOfMass(),1e-12))
        #
        m0Part=m0.deepCopy()[tab[0]] ; m0Part.renumberNodes(tab[1].invertArrayN2O2O2N(mm.getNumberOfNodes()),len(tab[1])) ; m0Part.setName(m0.getName())
        self.assertTrue(mmOut[0].isEqual(m0Part,1e-12))
        m1Part=m1.deepCopy()[tab[-1]] ; m1Part.renumberNodes(tab[1].invertArrayN2O2O2N(mm.getNumberOfNodes()),len(tab[1])) ; m1Part.setName(m0.getName())
        self.assertTrue(mmOut[0].isEqual(m0Part,1e-12))
        self.assertTrue(mmOut[-1].isEqual(m1Part,1e-12))
        #
        f0Part=f1ts.getFieldOnMeshAtLevel(ON_CELLS,0,mmOut) ; f0Part.checkConsistencyLight()
        self.assertEqual(f0Part.getTypeOfField(),ON_CELLS)
        self.assertTrue(f0Part.getMesh().isEqual(m0Part,1e-12))
        arr0Exp=DataArrayDouble([100,102,103,104,106,107]) ; arr0Exp.setInfoOnComponents(ic)
        self.assertTrue(f0Part.getArray().isEqual(arr0Exp,1e-12)) ; self.assertEqual(f0Part.getTime(),list(t0))
        f1Part=f1ts.getFieldOnMeshAtLevel(ON_CELLS,-1,mmOut) ; f1Part.checkConsistencyLight()
        self.assertEqual(f1Part.getTypeOfField(),ON_CELLS)
        self.assertTrue(f1Part.getMesh().isEqual(m1Part,1e-12))
        arr1Exp=DataArrayDouble([202,203,204,205,209]) ; arr1Exp.setInfoOnComponents(ic)
        self.assertTrue(f1Part.getArray().isEqual(arr1Exp,1e-12)) ; self.assertEqual(f1Part.getTime(),list(t0))
        #
        f2Part=f1ts.getFieldOnMeshAtLevel(ON_NODES,0,mmOut) ; f2Part.checkConsistencyLight()
        arr2Exp=DataArrayDouble([300,301,302,303,304,306,307,308,312,313]) ; arr2Exp.setInfoOnComponents(ic)
        self.assertTrue(f2Part.getArray().isEqual(arr2Exp,1e-12)) ; self.assertEqual(f2Part.getTime(),list(t0))
        # multisteps
        fs=MEDFileFields() ; fmts=MEDFileFieldMultiTS() ; fs.pushField(fmts)
        tss=[(16.5,3,4),(17.5,4,5),(18.5,5,6)]
        for i,tt in enumerate(tss):
            f0=MEDCouplingFieldDouble(ON_CELLS) ; f0.setTime(*tt)
            myarr=arr0_0+i*1000.
            f0.setArray(myarr)
            f0.setMesh(m0) ; f0.setName(fname0) ; f0.getArray().setInfoOnComponents(ic)
            f1ts=MEDFileField1TS() ; f1ts.setFieldNoProfileSBT(f0) ; fmts.pushBackTimeStep(f1ts)
            pass
        fsPart1=fs.extractPart(tab,mm)
        self.assertEqual(len(fsPart1),1)
        fmtsP=fsPart1[0]
        self.assertEqual(len(fmtsP),len(tss))
        for i,(f1tsP,tt) in enumerate(zip(fmtsP,tss)):
            fPart=f1tsP.field(mmOut) ; fPart.checkConsistencyLight()
            self.assertEqual(fPart.getTypeOfField(),ON_CELLS)
            arr0Exp=DataArrayDouble([100,102,103,104,106,107]) ; arr0Exp.setInfoOnComponents(ic) ; arr0Exp+=i*1000.
            self.assertTrue(fPart.getMesh().isEqual(m0Part,1e-12))
            self.assertTrue(fPart.getArray().isEqual(arr0Exp,1e-12))
            self.assertEqual(fPart.getTime(),list(tt))
            pass
        pass

    def testSymmetryPlusAggregationMFD1(self):
        """ Testing of MEDFileData::Aggregate and MEDFileUMesh::Aggregate and MEDFileUMesh::getAllDistributionOfType """
        fname1="Pyfile106_1.med"
        fname2="Pyfile106_2.med"
        fname3="Pyfile106_3.med"
        meshName="mesh"
        mm1=MEDFileUMesh()
        da1=DataArrayDouble([1,2,10,3,4,11,5,6,12,7,8,13],4,3) ; da1.setInfoOnComponents(["aa [m]","bbb [kg]","cccc [MW]"])
        mm1.setCoords(da1)
        mm1_0=MEDCouplingUMesh(meshName,3) ; mm1_0.allocateCells()
        mm1_0.setCoords(da1)
        mm1_0.insertNextCell(NORM_TETRA4,[0,1,2,3])
        mm1_0.insertNextCell(NORM_TETRA4,[4,5,6,7])
        mm1_0.insertNextCell(NORM_PENTA6,[8,9,10,11,12,13])
        mm1_0.insertNextCell(NORM_PENTA6,[14,15,16,17,18,19])
        mm1_0.insertNextCell(NORM_PENTA6,[20,21,22,23,24,25])
        mm1[0]=mm1_0
        mm1.setFamilyFieldArr(0,DataArrayInt([1,2,3,4,5]))
        mm1.setRenumFieldArr(0,DataArrayInt([11,12,13,14,15]))
        #
        mm1_1=MEDCouplingUMesh(meshName,2) ; mm1_1.allocateCells()
        mm1_1.setCoords(da1)
        mm1_1.insertNextCell(NORM_TRI3,[0,1,2])
        mm1_1.insertNextCell(NORM_TRI3,[3,4,5])
        mm1_1.insertNextCell(NORM_QUAD4,[6,7,8,9])
        mm1_1.insertNextCell(NORM_QUAD4,[10,11,12,13])
        mm1_1.insertNextCell(NORM_QUAD4,[14,15,16,17])
        mm1_1.insertNextCell(NORM_QUAD4,[18,19,20,21])
        mm1[-1]=mm1_1
        mm1.setFamilyFieldArr(-1,DataArrayInt([6,7,8,9,10,11]))
        mm1.setRenumFieldArr(-1,DataArrayInt([16,17,18,19,20,21]))
        for i in range(1,10):
            mm1.setFamilyId("F%d"%i,i)
        mm1.setFamilyId("FAMILLE_ZERO",0)
        mm1.setFamilyId("H1",100)
        mm1.setFamiliesOnGroup("myGRP",["F2","F6"])
        mm1.setFamiliesOnGroup("myGRP1",["F2","F6"])
        mm1.setFamilyFieldArr(1,DataArrayInt([12,13,14,15]))
        mm1.setRenumFieldArr(1,DataArrayInt([22,23,24,25]))
        ##############
        mm2=MEDFileUMesh()
        da1=DataArrayDouble([9,10,30,11,12,31,13,14,32,15,16,33,17,18,34],5,3) ; da1.setInfoOnComponents(["aa [m]","bbb [kg]","cccc [MW]"])
        mm2.setCoords(da1)
        mm2_0=MEDCouplingUMesh(meshName,3) ; mm2_0.allocateCells()
        mm2_0.setCoords(da1)
        mm2_0.insertNextCell(NORM_TETRA4,[100,101,102,103])
        mm2_0.insertNextCell(NORM_TETRA4,[104,105,106,107])
        mm2_0.insertNextCell(NORM_TETRA4,[108,109,110,111])
        mm2_0.insertNextCell(NORM_PENTA6,[112,113,114,115,116,117])
        mm2[0]=mm2_0
        mm2.setFamilyFieldArr(0,DataArrayInt([40,41,42,43]))
        mm2.setRenumFieldArr(0,DataArrayInt([50,51,52,53]))
        #
        mm2_1=MEDCouplingUMesh(meshName,2) ; mm2_1.allocateCells()
        mm2_1.setCoords(da1)
        mm2_1.insertNextCell(NORM_TRI3,[100,101,102])
        mm2_1.insertNextCell(NORM_TRI3,[103,104,105])
        mm2_1.insertNextCell(NORM_TRI3,[106,107,108])
        mm2_1.insertNextCell(NORM_QUAD4,[109,110,111,112])
        mm2_1.insertNextCell(NORM_QUAD4,[113,114,115,116])
        mm2_1.insertNextCell(NORM_QUAD4,[117,118,119,120])
        mm2_1.insertNextCell(NORM_QUAD4,[121,122,123,124])
        mm2_1.insertNextCell(NORM_QUAD4,[125,126,127,128])
        mm2[-1]=mm2_1
        mm2.setFamilyFieldArr(-1,DataArrayInt([200,201,202,203,204,205,206,207]))
        mm2.setRenumFieldArr(-1,DataArrayInt([300,301,302,303,304,305,306,307]))
        for i in range(1,12):
            mm2.setFamilyId("G%d"%i,i+30)
        mm2.setFamilyId("H1",100)
        mm2.setFamilyId("FAMILLE_ZERO",0)
        mm2.setFamiliesOnGroup("myGRP",["G2","G6"])
        mm2.setFamiliesOnGroup("myGRP2",["G4","G7"])
        mm2.setFamilyFieldArr(1,DataArrayInt([112,113,114,115,116]))
        mm2.setRenumFieldArr(1,DataArrayInt([122,123,124,125,126]))
        #
        mm=MEDFileUMesh.Aggregate([mm1,mm2])
        #######
        def CheckMesh(tester,mm):
            cooExp=DataArrayDouble([(1,2,10),(3,4,11),(5,6,12),(7,8,13),(9,10,30),(11,12,31),(13,14,32),(15,16,33),(17,18,34)]) ; cooExp.setInfoOnComponents(["aa [m]","bbb [kg]","cccc [MW]"])
            tester.assertTrue(mm.getCoords().isEqual(cooExp,1e-12))
            tester.assertTrue(mm[0].getNodalConnectivity().isEqual(DataArrayInt([14,0,1,2,3,14,4,5,6,7,14,104,105,106,107,14,108,109,110,111,14,112,113,114,115,16,8,9,10,11,12,13,16,14,15,16,17,18,19,16,20,21,22,23,24,25,16,116,117,118,119,120,121])))
            tester.assertTrue(mm[0].getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25,32,39,46,53])))
            tester.assertTrue(mm[-1].getNodalConnectivity().isEqual(DataArrayInt([3,0,1,2,3,3,4,5,3,104,105,106,3,107,108,109,3,110,111,112,4,6,7,8,9,4,10,11,12,13,4,14,15,16,17,4,18,19,20,21,4,113,114,115,116,4,117,118,119,120,4,121,122,123,124,4,125,126,127,128,4,129,130,131,132])))
            tester.assertTrue(mm[-1].getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8,12,16,20,25,30,35,40,45,50,55,60,65])))
            tester.assertTrue(mm.getFamilyFieldAtLevel(0).isEqual(DataArrayInt([1,2,40,41,42,3,4,5,43])))
            tester.assertTrue(mm.getNumberFieldAtLevel(0).isEqual(DataArrayInt([11,12,50,51,52,13,14,15,53])))
            tester.assertTrue(mm.getFamilyFieldAtLevel(-1).isEqual(DataArrayInt([6,7,200,201,202,8,9,10,11,203,204,205,206,207])))
            tester.assertTrue(mm.getNumberFieldAtLevel(-1).isEqual(DataArrayInt([16,17,300,301,302,18,19,20,21,303,304,305,306,307])))
            refFamIds=[("FAMILLE_ZERO",0),('F1',1),('F2',2),('F3',3),('F4',4),('F5',5),('F6',6),('F7',7),('F8',8),('F9',9),('G1',31),('G10',40),('G11',41),('G2',32),('G3',33),('G4',34),('G5',35),('G6',36),('G7',37),('G8',38),('G9',39),("H1",100)]
            tester.assertEqual(set(mm.getFamiliesNames()),set([elt[0] for elt in refFamIds]))
            tester.assertEqual(set([mm.getFamilyId(elt) for elt in mm.getFamiliesNames()]),set([elt[1] for elt in refFamIds]))
            tester.assertEqual(mm.getGroupsNames(),('myGRP','myGRP1','myGRP2'))
            tester.assertEqual(mm.getAllDistributionOfTypes(),[(NORM_TRI3,5),(NORM_QUAD4,9),(NORM_TETRA4,5),(NORM_PENTA6,4),(NORM_ERROR,9)])
            pass
        CheckMesh(self,mm)
        ##
        fieldName="zeField"
        t1=(2.3,3,5)
        t2=(5.6,7,12)
        infoc=["dd [W]","eee [kA]"]
        ##
        fmts1=MEDFileFieldMultiTS()
        f1ts1=MEDFileField1TS()
        f1_1=MEDCouplingFieldDouble(ON_CELLS) ; f1_1.setMesh(mm1[0]) ; f1_1.setName(fieldName)
        arr1=DataArrayDouble([(10,110),(11,111),(12,112),(13,113),(14,114)])
        arr1.setInfoOnComponents(infoc)
        f1_1.setArray(arr1) ; f1_1.setTime(*t1) ; f1_1.setTimeUnit("ms")
        f1_1.checkConsistencyLight()
        f1ts1.setFieldNoProfileSBT(f1_1)
        #
        f1_2=MEDCouplingFieldDouble(ON_CELLS) ; f1_2.setMesh(mm1[-1]) ; f1_2.setName(fieldName)
        arr2=DataArrayDouble([(15,115),(16,116),(17,117),(18,118),(19,119),(20,120)])
        arr2.setInfoOnComponents(infoc)
        f1_2.setArray(arr2) ; f1_2.setTime(*t1) ; f1_2.setTimeUnit("ms")
        f1_2.checkConsistencyLight()
        f1ts1.setFieldNoProfileSBT(f1_2)
        f1_3=MEDCouplingFieldDouble(ON_NODES) ; f1_3.setMesh(mm1[0]) ; f1_3.setName(fieldName)
        arr3=DataArrayDouble([(21,121),(22,122),(23,123),(24,124)])
        arr3.setInfoOnComponents(infoc)
        f1_3.setArray(arr3) ; f1_3.setTime(*t1) ; f1_3.setTimeUnit("ms")
        f1_3.checkConsistencyLight()
        f1ts1.setFieldNoProfileSBT(f1_3)
        fmts1.pushBackTimeStep(f1ts1)
        #
        f1ts2=f1ts1.deepCopy()
        f1ts2.setTime(t2[1],t2[2],t2[0])
        f1ts2.getUndergroundDataArray()[:]+=2000
        fmts1.pushBackTimeStep(f1ts2)
        ### fmts2
        fmts2=MEDFileFieldMultiTS()
        f1ts3=MEDFileField1TS()
        f2_1=MEDCouplingFieldDouble(ON_CELLS) ; f2_1.setMesh(mm2[0]) ; f2_1.setName(fieldName)
        arr4=DataArrayDouble([(50,150),(51,151),(52,152),(53,153)])
        arr4.setInfoOnComponents(infoc)
        f2_1.setArray(arr4) ; f2_1.setTime(*t1) ; f2_1.setTimeUnit("ms")
        f2_1.checkConsistencyLight()
        f1ts3.setFieldNoProfileSBT(f2_1)
        f2_2=MEDCouplingFieldDouble(ON_CELLS) ; f2_2.setMesh(mm2[-1]) ; f2_2.setName(fieldName)
        arr5=DataArrayDouble([(54,154),(55,155),(56,156),(57,157),(158,158),(59,159),(60,160),(61,161)])
        arr5.setInfoOnComponents(infoc)
        f2_2.setArray(arr5) ; f2_2.setTime(*t1) ; f2_2.setTimeUnit("ms")
        f2_2.checkConsistencyLight()
        f1ts3.setFieldNoProfileSBT(f2_2)
        f2_3=MEDCouplingFieldDouble(ON_NODES) ; f2_3.setMesh(mm2[0]) ; f2_3.setName(fieldName)
        arr6=DataArrayDouble([(62,162),(63,163),(64,164),(65,165),(66,166)])
        arr6.setInfoOnComponents(infoc)
        f2_3.setArray(arr6) ; f2_3.setTime(*t1) ; f2_3.setTimeUnit("ms")
        f2_3.checkConsistencyLight()
        f1ts3.setFieldNoProfileSBT(f2_3)
        fmts2.pushBackTimeStep(f1ts3)
        #
        f1ts4=f1ts3.deepCopy()
        f1ts4.setTime(t2[1],t2[2],t2[0])
        f1ts4.getUndergroundDataArray()[:]+=2000
        fmts2.pushBackTimeStep(f1ts4)
        #
        mfd1=MEDFileData()
        mfd1.setMeshes(MEDFileMeshes())
        mfd1.getMeshes().pushMesh(mm1)
        mfd1.setFields(MEDFileFields())
        mfd1.getFields().pushField(fmts1)
        #
        mfd2=MEDFileData()
        mfd2.setMeshes(MEDFileMeshes())
        mfd2.getMeshes().pushMesh(mm2)
        mfd2.setFields(MEDFileFields())
        mfd2.getFields().pushField(fmts2)
        # ze Call !
        mfd=MEDFileData.Aggregate([mfd1,mfd2])
        def CheckMFD(tester,mfd):
            tester.assertEqual(len(mfd.getMeshes()),1)
            tester.assertEqual(len(mfd.getFields()),1)
            CheckMesh(self,mfd.getMeshes()[0])
            tester.assertEqual(len(mfd.getFields()[0]),2)
            zeF1=mfd.getFields()[0][0]
            zeF1_1=zeF1.getFieldOnMeshAtLevel(ON_CELLS,0,mfd.getMeshes()[0])
            ref=MEDCouplingFieldDouble.MergeFields([f1_1,f2_1])
            o2n=ref.getMesh().deepCopy().sortCellsInMEDFileFrmt()
            ref.renumberCells(o2n)
            tester.assertTrue(ref.isEqual(zeF1_1,1e-12,1e-12))
            zeF1_2=zeF1.getFieldOnMeshAtLevel(ON_CELLS,-1,mfd.getMeshes()[0])
            ref=MEDCouplingFieldDouble.MergeFields([f1_2,f2_2])
            o2n=ref.getMesh().deepCopy().sortCellsInMEDFileFrmt()
            ref.renumberCells(o2n)
            tester.assertTrue(ref.isEqual(zeF1_2,1e-12,1e-12))
            zeF1_3=zeF1.getFieldOnMeshAtLevel(ON_NODES,0,mfd.getMeshes()[0])
            ref=MEDCouplingFieldDouble.MergeFields([f1_3,f2_3])
            o2n=ref.getMesh().deepCopy().sortCellsInMEDFileFrmt()
            ref.renumberCells(o2n)
            tester.assertTrue(ref.isEqual(zeF1_3,1e-12,1e-12))
            #
            zeF2=mfd.getFields()[0][1]
            zeF2_1=zeF2.getFieldOnMeshAtLevel(ON_CELLS,0,mfd.getMeshes()[0])
            ref=MEDCouplingFieldDouble.MergeFields([f1_1,f2_1])
            o2n=ref.getMesh().deepCopy().sortCellsInMEDFileFrmt()
            ref.renumberCells(o2n)
            ref.setTime(*t2) ; ref.getArray()[:]+=2000
            tester.assertTrue(ref.isEqual(zeF2_1,1e-12,1e-12))
            zeF2_2=zeF2.getFieldOnMeshAtLevel(ON_CELLS,-1,mfd.getMeshes()[0])
            ref=MEDCouplingFieldDouble.MergeFields([f1_2,f2_2])
            o2n=ref.getMesh().deepCopy().sortCellsInMEDFileFrmt()
            ref.renumberCells(o2n)
            ref.setTime(*t2) ; ref.getArray()[:]+=2000
            tester.assertTrue(ref.isEqual(zeF2_2,1e-12,1e-12))
            zeF2_3=zeF2.getFieldOnMeshAtLevel(ON_NODES,0,mfd.getMeshes()[0])
            ref=MEDCouplingFieldDouble.MergeFields([f1_3,f2_3])
            o2n=ref.getMesh().deepCopy().sortCellsInMEDFileFrmt()
            ref.renumberCells(o2n)
            ref.setTime(*t2) ; ref.getArray()[:]+=2000
            tester.assertTrue(ref.isEqual(zeF2_3,1e-12,1e-12))
        CheckMFD(self,mfd)
        mfd1.write(fname1,2) ; mfd2.write(fname2,2)
        mfd=MEDFileData.Aggregate([MEDFileData(fname1),MEDFileData(fname2)])
        CheckMFD(self,mfd)
        pass

    def testExtrudedMesh1(self):
        fname="Pyfile107.med"
        arrX=DataArrayDouble([0,1,2,3]) ; arrY=DataArrayDouble([0,1,2,3,4]) ; arrZ=DataArrayDouble([0,1,2,3,4,5])
        mesh3D=MEDCouplingCMesh() ; mesh3D.setCoords(arrX,arrY,arrZ) ; mesh3D.setName("mesh")
        ex=MEDCouplingMappedExtrudedMesh(mesh3D)
        mm=MEDFileUMesh(ex)
        mm.write(fname,2)
        ex2=mm.convertToExtrudedMesh()
        mm2=MEDFileMesh.New(fname)
        ex3=mm2.convertToExtrudedMesh()
        self.assertTrue(ex.isEqual(ex2,1e-12))
        self.assertTrue(ex.isEqual(ex3,1e-12))
        pass
    
    @unittest.skipUnless(LooseVersion(MEDFileVersionStr())>=LooseVersion('3.2.1'),"This test requires at least MEDFile version 3.2.1")
    def testWriteInto30(self):
        fname="Pyfile108.med"
        fname2="Pyfile109.med"
        m=MEDCouplingUMesh("mesh",1) ; m.setCoords(DataArrayDouble([0,0,1,1],2,2)) ; m.allocateCells() ; m.insertNextCell(NORM_SEG2,[1,0])
        mm=MEDFileUMesh() ; mm[0]=m
        mm.setFamilyId("FAMILLE_ZERO",0)
        #
        mm.write30(fname,2)
        assert(LooseVersion(MEDFileVersionOfFileStr(fname)).version[:2]==[3,0]) # checks that just written MED file has a version == 3.0.x
        mm2=MEDFileUMesh(fname)
        self.assertTrue(mm.isEqual(mm2,1e-12))
        #
        mm.write(fname2,2)
        assert(LooseVersion(MEDFileVersionOfFileStr(fname2)).version[:2]==list(MEDFileVersion()[:2])) # checks that MED file version of written mesh is thoose of the current MED file lib
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def testPickelizationOfMEDFileObjects1(self):
        fname="Pyfile110.med"
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
        for i in range(9):
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
        fieldNode1=MEDCouplingFieldDouble(ON_NODES,ONE_TIME)
        fieldNode1.setName("fieldNode1")
        fieldNode1.setTime(tim,dt,it)
        pfl1=DataArrayInt([0,1,2,3,4,5,6]) ; pfl1.setName("PflIdentity1")
        arr1=DataArrayDouble([20,21,22,23,24,25,26])
        fieldNode1.setArray(arr1)
        f1=MEDFileField1TS()
        f1.setFieldProfile(fieldNode1,m,-1,pfl1)
        mfd=MEDFileData()
        mfd.setMeshes(MEDFileMeshes()) ; mfd.setFields(MEDFileFields())
        mfd.getMeshes().pushMesh(m)
        fmts=MEDFileFieldMultiTS() ; fmts.pushBackTimeStep(f0)
        mfd.getFields().pushField(fmts)
        # first start gently
        d=mfd.serialize()
        mfd2=MEDFileData(d)
        self.assertEqual(len(mfd2.getMeshes()),1)
        self.assertEqual(len(mfd2.getFields()),1)
        self.assertEqual(len(mfd2.getFields()[0]),1)
        self.assertTrue(mfd2.getMeshes()[0].isEqual(mfd.getMeshes()[0],1e-12))
        ff2=mfd2.getFields()[0][0].field(mfd2.getMeshes()[0])
        ff =mfd.getFields()[0][0].field(mfd.getMeshes()[0])
        self.assertTrue(ff2.isEqual(ff,1e-12,1e-12))
        # OK now end of joke -> serialization of MEDFileData
        st=pickle.dumps(mfd,pickle.HIGHEST_PROTOCOL)
        mfd3=pickle.loads(st)
        # check of object
        self.assertEqual(len(mfd3.getMeshes()),1)
        self.assertEqual(len(mfd3.getFields()),1)
        self.assertEqual(len(mfd3.getFields()[0]),1)
        self.assertTrue(mfd3.getMeshes()[0].isEqual(mfd.getMeshes()[0],1e-12))
        ff3=mfd3.getFields()[0][0].field(mfd3.getMeshes()[0])
        self.assertTrue(ff3.isEqual(ff,1e-12,1e-12))
        # serialization of MEDFileFields
        st=pickle.dumps(mfd.getFields(),pickle.HIGHEST_PROTOCOL)
        fs4=pickle.loads(st)
        ff4=fs4[0][0].field(mfd3.getMeshes()[0])
        self.assertTrue(ff4.isEqual(ff,1e-12,1e-12))
        # serialization of MEDFileFieldMulitTS
        st=pickle.dumps(mfd.getFields()[0],pickle.HIGHEST_PROTOCOL)
        fmts5=pickle.loads(st)
        ff5=fmts5[0].field(mfd3.getMeshes()[0])
        self.assertTrue(ff5.isEqual(ff,1e-12,1e-12))
        # serialization of MEDFileField1TS
        st=pickle.dumps(mfd.getFields()[0][0],pickle.HIGHEST_PROTOCOL)
        f1ts6=pickle.loads(st)
        ff6=f1ts6.field(mfd3.getMeshes()[0])
        self.assertTrue(ff6.isEqual(ff,1e-12,1e-12))
        # serialization of MEDFileMeshes
        st=pickle.dumps(mfd.getMeshes(),pickle.HIGHEST_PROTOCOL)
        ms7=pickle.loads(st)
        self.assertEqual(len(ms7),1)
        self.assertTrue(ms7[0].isEqual(mfd.getMeshes()[0],1e-12))
        pass

    @unittest.skipUnless(MEDCouplingHasNumPyBindings(),"requires numpy")
    def testPickelizationOfMEDFileObjects2(self):
        # CMesh
        self.testMEDMesh6() # generates MEDFileMesh5.med file
        mm=MEDFileMesh.New("MEDFileMesh5.med")
        self.assertTrue(isinstance(mm,MEDFileCMesh))
        st=pickle.dumps(mm,pickle.HIGHEST_PROTOCOL)
        mm2=pickle.loads(st)
        self.assertTrue(isinstance(mm2,MEDFileCMesh))
        self.assertTrue(mm.getMesh().isEqual(mm2.getMesh(),1e-12))
        # CurveLinear
        self.testCurveLinearMesh1() # generates Pyfile55.med
        mm=MEDFileMesh.New("Pyfile55.med")
        self.assertTrue(isinstance(mm,MEDFileCurveLinearMesh))
        st=pickle.dumps(mm,pickle.HIGHEST_PROTOCOL)
        mm3=pickle.loads(st)
        self.assertTrue(isinstance(mm3,MEDFileCurveLinearMesh))
        self.assertTrue(mm.getMesh().isEqual(mm3.getMesh(),1e-12))
        self.testInt32InMEDFileFieldStar1()# generates Pyfile63.med
        # MEDFileIntFieldMultiTS
        fs4=MEDFileFields("Pyfile63.med")
        ms4=MEDFileMeshes("Pyfile63.med")
        self.assertTrue(isinstance(fs4[0],MEDFileIntFieldMultiTS))
        st=pickle.dumps(fs4[0],pickle.HIGHEST_PROTOCOL)
        fmts5=pickle.loads(st)
        self.assertEqual(len(fs4[0]),len(fmts5))
        self.assertTrue(isinstance(fmts5,MEDFileIntFieldMultiTS))
        self.assertTrue(fmts5[0].field(ms4[0]).isEqual((fs4[0][0]).field(ms4[0]),1e-12,0))
        # MEDFileIntField1TS
        st=pickle.dumps(fs4[0][0],pickle.HIGHEST_PROTOCOL)
        f1ts6=pickle.loads(st)
        self.assertTrue(isinstance(f1ts6,MEDFileIntField1TS))
        self.assertTrue(f1ts6.field(ms4[0]).isEqual((fs4[0][0]).field(ms4[0]),1e-12,0))
        # MEDFileParameters
        self.testParameters1()# generates Pyfile56.med
        params=MEDFileParameters("Pyfile56.med")
        st=pickle.dumps(params,pickle.HIGHEST_PROTOCOL)
        params7=pickle.loads(st)
        self.assertEqual(len(params),len(params7))
        for i in range(len(params)):
            self.assertTrue(params[i].isEqual(params7[i],1e-12)[0])
            pass
        pass

    def testGlobalNumOnNodes1(self):
        """Test global number on nodes here. Used by partitionners."""
        fname="Pyfile112.med"
        arr=DataArrayDouble(5) ; arr.iota()
        m=MEDCouplingUMesh.Build1DMeshFromCoords(arr)
        m.setName("mesh")
        mm=MEDFileUMesh()
        mm[0]=m
        self.assertTrue(not mm.getGlobalNumFieldAtLevel(1))
        d=DataArrayInt([7,8,9,2,0])
        dRef=d.deepCopy()
        mm.setGlobalNumFieldAtLevel(1,d)
        mm.checkConsistency()
        self.assertRaises(InterpKernelException,mm.setGlobalNumFieldAtLevel,1,d[::2])
        mm.checkConsistency()
        self.assertEqual(d.getHiddenCppPointer(),mm.getGlobalNumFieldAtLevel(1).getHiddenCppPointer())
        self.assertTrue(mm.getGlobalNumFieldAtLevel(1).isEqual(dRef))
        mm.write(fname,2)
        mm2=MEDFileMesh.New(fname)
        self.assertTrue(mm.isEqual(mm2,1e-12)[0])
        self.assertTrue(mm2.getGlobalNumFieldAtLevel(1).isEqual(dRef))
        mm2.getGlobalNumFieldAtLevel(1).setIJ(0,0,10)
        self.assertTrue(not mm.isEqual(mm2,1e-12)[0])
        mm2.getGlobalNumFieldAtLevel(1).setIJ(0,0,7)
        self.assertTrue(mm.isEqual(mm2,1e-12)[0])
        pass

    def testPartialReadOfEntities1(self):
        """Test for advanced API on read to speed up read phase for users with "huge" number of time steps (more than 10 000)."""
        fname="Pyfile113.med"
        arr=DataArrayDouble(5) ; arr.iota()
        m=MEDCouplingUMesh.Build1DMeshFromCoords(arr)
        m.setName("mesh")
        mm=MEDFileUMesh()
        mm[0]=m
        #
        fieldName="Field"
        ts1=(5.,1,2)
        f1=MEDCouplingFieldDouble(ON_NODES) ; f1.setMesh(m) ; f1.setName(fieldName)
        f1.setArray(DataArrayDouble([0.,0.1,0.2,0.3,0.4]))
        f1.setTime(*ts1)
        f2=MEDCouplingFieldDouble(ON_CELLS) ; f2.setMesh(m) ; f2.setName(fieldName)
        f2.setArray(DataArrayDouble([1.,1.1,1.2,1.3]))
        f2.setTime(*ts1)
        f1ts=MEDFileField1TS()
        f1ts.setFieldNoProfileSBT(f1)
        f1ts.setFieldNoProfileSBT(f2)
        self.assertEqual(set(f1ts.getTypesOfFieldAvailable()),set([ON_NODES,ON_CELLS]))
        f1ts_2=f1ts.deepCopy()
        f1ts_2.getUndergroundDataArray()[:]+=2
        f1ts_2.setTime(3,4,6.)
        fmts=MEDFileFieldMultiTS()
        fmts.pushBackTimeStep(f1ts)
        fmts.pushBackTimeStep(f1ts_2)
        #
        mm.write(fname,2)
        fmts.write(fname,0)
        #
        ent=MEDFileEntities.BuildFrom([(ON_NODES,NORM_ERROR)])
        mm=MEDFileMesh.New(fname)
        fs=MEDFileFields(fname,False,ent) # the important line is here - We specify to MEDFileFields to read only nodes part to speed up read phase (by avoiding to scan all entities time geo types)
        fs.loadArrays()
        self.assertEqual(len(fs),1)
        fmts=fs[0]
        self.assertEqual(len(fmts),2)
        ff0=fmts[0] ; ff1=fmts[1]
        self.assertEqual(ff0.getTypesOfFieldAvailable(),[ON_NODES]) # only NODES have been loaded
        self.assertTrue(ff0.field(mm).isEqual(f1,1e-12,1e-12))
        f3=f1.deepCopy() ; f3+=2. ; f3.setTime(6.,3,4)
        self.assertTrue(ff1.field(mm).isEqual(f3,1e-12,1e-12))
        pass
    
    def testFloat32InMEDFileFieldStar1(self):
        """Like testInt32InMEDFileFieldStar1 but with float32 :)"""
        fname="Pyfile114.med"
        f1=MEDLoaderDataForTest.buildVecFieldOnCells_1();
        f1=f1.convertToFloatField()
        m1=f1.getMesh()
        mm1=MEDFileUMesh.New()
        mm1.setCoords(m1.getCoords())
        mm1.setMeshAtLevel(0,m1)
        mm1.setName(m1.getName())
        mm1.write(fname,2)
        ff1=MEDFileFloatField1TS()
        ff1.setFieldNoProfileSBT(f1)
        a=ff1.getFieldOnMeshAtLevel(0,ON_CELLS,mm1)
        self.assertEqual(a.getArray().getInfoOnComponents(),['power [MW/m^3]','density [g/cm^3]','temperature [K]'])
        self.assertTrue(a.isEqual(f1,1e-12,1e-12))
        ff1.write(fname,0)
        ff2=MEDFileAnyTypeField1TS.New(fname)
        self.assertEqual(ff2.getName(),"VectorFieldOnCells")
        self.assertEqual(ff2.getTime(),[0,1,2.0])
        self.assertTrue(isinstance(ff2,MEDFileFloatField1TS))
        a=ff1.getFieldOnMeshAtLevel(ON_CELLS,0,mm1)
        self.assertEqual(a.getArray().getInfoOnComponents(),['power [MW/m^3]','density [g/cm^3]','temperature [K]'])
        self.assertTrue(a.isEqual(f1,1e-12,1e-12))
        ff2.setTime(1,2,3.)
        c=ff2.getUndergroundDataArray() ; c*=2
        ff2.write(fname,0) # 2 time steps in 
        ffs1=MEDFileAnyTypeFieldMultiTS.New(fname,"VectorFieldOnCells")
        self.assertEqual(ffs1.getTimeSteps(),[(0, 1, 2.0), (1, 2, 3.0)])
        self.assertEqual(len(ffs1),2)
        self.assertTrue(isinstance(ffs1,MEDFileFloatFieldMultiTS))
        a=ffs1[2.].getFieldOnMeshAtLevel(ON_CELLS,0,mm1)
        self.assertTrue(a.isEqual(f1,1e-12,1e-12))
        a=ffs1.getFieldOnMeshAtLevel(ON_CELLS,0,1,0,mm1)
        self.assertTrue(a.isEqual(f1,1e-12,1e-12))
        it=ffs1.__iter__() ; it.next() ; ff2bis=it.next()
        a=ff2bis.getFieldOnMeshAtLevel(0,ON_CELLS,mm1)
        self.assertTrue(a.getArray().isEqual(2*f1.getArray(),1e-7))
        f1.setTime(3.,1,2) ; f1.getArray()[:]*=2
        self.assertTrue(a.isEqual(f1,1e-12,1e-12)) ; f1.getArray()[:]/=2
        bc=DataArrayFloat(6,3) ; bc[:]=0 ; bc.setInfoOnComponents(['power [MW/m^3]','density [g/cm^3]','temperature [K]'])
        for it in ffs1:
            a=it.getFieldOnMeshAtLevel(ON_CELLS,0,mm1)
            bc+=a.getArray()
            pass
        self.assertTrue(bc.isEqual(3*f1.getArray(),1e-7))
        nf1=MEDCouplingFieldFloat(ON_NODES)
        nf1.setTime(9.,10,-1)
        nf1.setMesh(f1.getMesh())
        narr=DataArrayFloat(12,2) ; narr.setInfoOnComponents(["aa [u1]","bbbvv [ppp]"]) ; narr[:,0]=list(range(12)) ; narr[:,1]=2*narr[:,0]
        nf1.setName("VectorFieldOnNodes") ; nf1.setArray(narr)
        nff1=MEDFileFloatField1TS.New()
        nff1.setFieldNoProfileSBT(nf1)
        self.assertEqual(nff1.getInfo(),('aa [u1]','bbbvv [ppp]'))
        self.assertEqual(nff1.getTime(),[10,-1,9.0])
        nff1.write(fname,0)
        #
        nf2=MEDCouplingFieldFloat(ON_NODES)
        nf2.setTime(19.,20,-11)
        nf2.setMesh(f1.getMesh())
        narr2=DataArrayFloat(8,2) ; narr.setInfoOnComponents(["aapfl [u1]","bbbvvpfl [ppp]"]) ; narr2[:,0]=list(range(8)) ; narr2[:,0]+=10  ; narr2[:,1]=3*narr2[:,0]
        nf2.setName("VectorFieldOnNodesPfl") ; narr2.setName(nf2.getName()) ; nf2.setArray(narr2)
        nff2=MEDFileFloatField1TS.New()
        npfl=DataArrayInt([1,2,4,5,6,7,10,11]) ; npfl.setName("npfl")
        nff2.setFieldProfile(nf2,mm1,0,npfl)
        nff2.getFieldWithProfile(ON_NODES,0,mm1)
        a,b=nff2.getFieldWithProfile(ON_NODES,0,mm1) ; b.setName(npfl.getName())
        self.assertTrue(b.isEqual(npfl))
        self.assertTrue(a.isEqual(narr2,1e-7))
        nff2.write(fname,0)
        nff2bis=MEDFileFloatField1TS(fname,"VectorFieldOnNodesPfl")
        a,b=nff2bis.getFieldWithProfile(ON_NODES,0,mm1) ; b.setName(npfl.getName())
        self.assertTrue(b.isEqual(npfl))
        self.assertTrue(a.isEqual(narr2,1e-7))
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
        self.assertTrue(isinstance(ffs[0],MEDFileFloatFieldMultiTS))
        self.assertTrue(isinstance(ffs[1],MEDFileFloatFieldMultiTS))
        self.assertTrue(isinstance(ffs[2],MEDFileFieldMultiTS))
        self.assertTrue(isinstance(ffs[3],MEDFileFloatFieldMultiTS))
        #
        self.assertTrue(fs["VectorFieldOnCells"][0].getUndergroundDataArray().isEqualWithoutConsideringStr(f1.getArray(),1e-7))
        self.assertTrue(fs["VectorFieldOnCells"][1,2].getUndergroundDataArray().isEqualWithoutConsideringStr(2*f1.getArray(),1e-7))
        self.assertTrue(fs["VectorFieldOnNodesPfl"][0].getUndergroundDataArray().isEqualWithoutConsideringStr(narr2,1e-7))
        self.assertTrue(fs["VectorFieldOnNodes"][9.].getUndergroundDataArray().isEqualWithoutConsideringStr(narr,1e-7))
        self.assertTrue(fs["VectorFieldOnNodesDouble"][29.].getUndergroundDataArray().isEqualWithoutConsideringStr(f1.getMesh().getCoords(),1e-12))
        #
        nf3_read=MEDFileFieldMultiTS(fname,"VectorFieldOnNodesDouble")
        self.assertTrue(nf3_read[29.].getUndergroundDataArray().isEqualWithoutConsideringStr(f1.getMesh().getCoords(),1e-12))
        self.assertRaises(InterpKernelException,MEDFileFloatFieldMultiTS.New,fname,"VectorFieldOnNodesDouble")# exception because trying to read a double field with int instance
        self.assertRaises(InterpKernelException,MEDFileFieldMultiTS.New,fname,"VectorFieldOnNodes")# exception because trying to read a int field with double instance
        MEDFileField1TS.New(fname,"VectorFieldOnNodesDouble",30,-21)
        self.assertRaises(InterpKernelException,MEDFileFloatField1TS.New,fname,"VectorFieldOnNodesDouble",30,-21)# exception because trying to read a double field with int instance
        MEDFileFloatField1TS.New(fname,"VectorFieldOnNodes",10,-1)
        self.assertRaises(InterpKernelException,MEDFileField1TS.New,fname,"VectorFieldOnNodes",10,-1)# exception because trying to read a double field with int instance
        #
        self.assertEqual(fs.getMeshesNames(),('3DSurfMesh_1','3DSurfMesh_1','3DSurfMesh_1','3DSurfMesh_1'))
        self.assertTrue(fs.changeMeshNames([('3DSurfMesh_1','3DSurfMesh')]))
        self.assertEqual(fs.getMeshesNames(),('3DSurfMesh','3DSurfMesh','3DSurfMesh','3DSurfMesh'))
        self.assertTrue(not fs.changeMeshNames([('3DSurfMesh_1','3DSurfMesh')]))
        pass

    def testPenta18_1(self):
        """EDF8478 : Test of read/write of penta18"""
        fname="Pyfile115.med"
        arr=DataArrayDouble([
            (0.,1.,1.),(0.,0.,1.),(1.,0.,1.),
            (0.,1.,0.),(0.,0.,0.),(1.,0.,0.),
            (0.,0.5,1.),(0.5,0.,1.),(0.5,0.5,1.),
            (0.,0.5,0.),(0.5,0.,0.),(0.5,0.5,0.),
            (0.,1.,0.5),(0.,0.,0.5),(1.,0.,0.5),
            (0.,0.5,0.5),(0.5,0.,0.5),(0.5,0.5,0.5)])
        m=MEDCouplingUMesh("mesh",3)
        m.setCoords(arr)
        m.allocateCells(1)
        m.insertNextCell(NORM_PENTA18,list(range(18)))
        m.checkConsistencyLight()
        #
        f=MEDCouplingFieldDouble(ON_NODES)
        f.setMesh(m)
        f.setName("FieldOnPenta18")
        f.setArray(DataArrayDouble(list(range(18))))
        f.checkConsistencyLight()
        #
        m2,d,di,rd,rdi=m.buildDescendingConnectivity()
        #
        f2=MEDCouplingFieldDouble(ON_NODES)
        f2.setMesh(m)
        f2.setName("FieldOnPenta18Sub")
        f2.setArray(DataArrayDouble(list(range(18))))
        f2.checkConsistencyLight()
        WriteField(fname,f2,True)
        f3=ReadField(fname)
        self.assertTrue(f2.isEqual(f3,1e-12,1e-12))
        self.assertEqual(f3.getMesh().getNumberOfCells(),1)
        self.assertEqual(f3.getMesh().getTypeOfCell(0),NORM_PENTA18)
        pass
    
    pass

if __name__ == "__main__":
    unittest.main()
