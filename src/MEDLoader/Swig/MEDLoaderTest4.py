#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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
# Author : Anthony Geay (CEA/DEN)

from MEDLoader import *
import unittest
from math import pi,e,sqrt

class MEDLoaderTest4(unittest.TestCase):
    """
    Test series to emulate the future MEDReader plugin for PARAVIS.
    """
    def test1(self):
        """
        This test is the most simple one. One time serie of one field with only cell fields with no profiles.
        The only "difficulty" is that the cell field is lying on different levels (2D and 1D) to maximize the compatibility with ParaVIS.
        """
        fname="ForMEDReader1.med"
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
        fam=DataArrayInt(9) ; fam.iota(0) ; mm.setFamilyFieldArr(0,fam)
        fam=DataArrayInt(32) ; fam.iota(20) ; mm.setFamilyFieldArr(-1,fam) ; del fam
        num=DataArrayInt(9) ; num.iota(100) ; mm.setRenumFieldArr(0,num)
        num=DataArrayInt(32) ; num.iota(120) ; mm.setRenumFieldArr(-1,num) ; del num
        #
        fieldName="zeField"
        fs=MEDFileFieldMultiTS()
        ##### Time step 0
        i=0
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        fs.pushBackTimeStep(f)
        ##### Time step 1
        i=1
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(1100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(1200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        fs.pushBackTimeStep(f)
        ##### Time step 2
        i=2
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(2100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(2200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        fs.pushBackTimeStep(f)
        ##### Time step 3
        i=3
        f=MEDFileField1TS()
        #
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(3100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(3200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        #
        fs.pushBackTimeStep(f)
        ##### Time step 4
        i=4
        f=MEDFileField1TS()
        #
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(4100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(4200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        fs.pushBackTimeStep(f)
        mm.write(fname,2)
        fs.write(fname,0)
        a0Exp=mm.getCoords().deepCpy()
        del m,m1,mm,fs,f,fCell0,fCell1
        ########## GO for reading in MEDReader, by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False) # False is important to not read the values
        refMem=fields.getHeapMemorySize()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),1)
        for fmts in allFMTSLeavesToDisplay[0]:
            self.assertEqual(fmts.getTimeSteps(),[(0,0,0.),(1,0,1.),(2,0,2.),(3,0,3.),(4,0,4.)]) # All discretizations have the same time series
            pass
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),1)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        #
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        for i in xrange(1,5):
            self.assertTrue(fcscp.isDataSetSupportEqualToThePreviousOne(i,fields))
            pass
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,5,5,5,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,100,104,108,112,117,122,127,132])))
        self.assertTrue(a3.isEqual(DataArrayInt([2,0,1,2,1,2,2,2,0,2,3,4,2,4,5,2,5,3,2,6,7,2,7,8,2,8,6,2,9,10,2,10,11,2,11,9,2,12,13,2,13,14,2,14,15,2,15,12,2,16,17,2,17,18,2,18,19,2,19,16,2,20,21,2,21,22,2,22,23,2,23,20,2,24,25,2,25,26,2,26,27,2,27,24,2,28,29,2,29,30,2,30,31,2,31,28,3,0,1,2,3,3,4,5,3,6,7,8,3,9,10,11,4,12,13,14,15,4,16,17,18,19,4,20,21,22,23,4,24,25,26,27,4,28,29,30,31])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,0,1,2,3,4,5,6,7,8])))
        self.assertTrue(not a7)
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,100,101,102,103,104,105,106,107,108])))
        self.assertTrue(not a9)
        for i in xrange(5):
            fsst=MEDFileField1TSStructItem.BuildItemFrom(fields[0][i],mst)
            fields[0][i].loadArraysIfNecessary()
            tmpMem=fields.getHeapMemorySize()
            self.assertTrue(tmpMem-refMem>=41*2*8)
            refMem=tmpMem
            v=mml.buildDataArray(fsst,fields,fields[0][i].getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),fields[0][i].getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([200.,201.,202.,203.,204.,205.,206.,207.,208.,209.,210.,211.,212.,213.,214.,215.,216.,217.,218.,219.,220.,221.,222.,223.,224.,225.,226.,227.,228.,229.,230.,231.,232.,233.,234.,235.,236.,237.,238.,239.,240.,241.,242.,243.,244.,245.,246.,247.,248.,249.,250.,251.,252.,253.,254.,255.,256.,257.,258.,259.,260.,261.,262.,263.,100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,110.,111.,112.,113.,114.,115.,116.,117.],41,2) ; vExp.setInfoOnComponents(['Comp1 [m]','Com2 [s^2]']) ; vExp+=i*1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass

    def test2(self):
        """
        One time serie of one field with cell and node discretization in the same field with no profiles.
        Here as there is no profile only one VTK support is requested.
        """
        fname="ForMEDReader2.med"
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
        fieldName="zeField"
        fs=MEDFileFieldMultiTS()
        ##### Time step 0
        i=0
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        #
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs.pushBackTimeStep(f)
        ##### Time step 1
        i=1
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(1100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(1200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        #
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(1300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs.pushBackTimeStep(f)
        ##### Time step 2
        i=2
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(2300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        #
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(2100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(2200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        fs.pushBackTimeStep(f)
        ##### Time step 3
        i=3
        f=MEDFileField1TS()
        #
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(3100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(3200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        #
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(3300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        #
        fs.pushBackTimeStep(f)
        ##### Time step 4
        i=4
        f=MEDFileField1TS()
        #
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(4100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(4200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        #
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(4300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        #
        fs.pushBackTimeStep(f)
        mm.write(fname,2)
        fs.write(fname,0)
        a0Exp=mm.getCoords().deepCpy()
        del m,m1,mm,fs,f,fCell0,fCell1
        ########## GO for reading in MEDReader, by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),2)
        for fmts in allFMTSLeavesToDisplay[0]:
            self.assertEqual(fmts.getTimeSteps(),[(0,0,0.),(1,0,1.),(2,0,2.),(3,0,3.),(4,0,4.)]) # All discretizations have the same time series
            pass
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),2)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),2)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        assert isinstance(mml2,MEDUMeshMultiLev)
        for i in xrange(1,5):
            self.assertTrue(fcscp.isDataSetSupportEqualToThePreviousOne(i,fields))
            pass
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,5,5,5,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,100,104,108,112,117,122,127,132])))
        self.assertTrue(a3.isEqual(DataArrayInt([2,0,1,2,1,2,2,2,0,2,3,4,2,4,5,2,5,3,2,6,7,2,7,8,2,8,6,2,9,10,2,10,11,2,11,9,2,12,13,2,13,14,2,14,15,2,15,12,2,16,17,2,17,18,2,18,19,2,19,16,2,20,21,2,21,22,2,22,23,2,23,20,2,24,25,2,25,26,2,26,27,2,27,24,2,28,29,2,29,30,2,30,31,2,31,28,3,0,1,2,3,3,4,5,3,6,7,8,3,9,10,11,4,12,13,14,15,4,16,17,18,19,4,20,21,22,23,4,24,25,26,27,4,28,29,30,31])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        # for cells
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)# Second 0 is for cells
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([200.,201.,202.,203.,204.,205.,206.,207.,208.,209.,210.,211.,212.,213.,214.,215.,216.,217.,218.,219.,220.,221.,222.,223.,224.,225.,226.,227.,228.,229.,230.,231.,232.,233.,234.,235.,236.,237.,238.,239.,240.,241.,242.,243.,244.,245.,246.,247.,248.,249.,250.,251.,252.,253.,254.,255.,256.,257.,258.,259.,260.,261.,262.,263.,100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,110.,111.,112.,113.,114.,115.,116.,117.],41,2) ; vExp.setInfoOnComponents(['Comp1 [m]','Com2 [s^2]']) ; vExp+=i*1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)# Second 0 is for cells
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([300.,301.,302.,303.,304.,305.,306.,307.,308.,309.,310.,311.,312.,313.,314.,315.,316.,317.,318.,319.,320.,321.,322.,323.,324.,325.,326.,327.,328.,329.,330.,331.,332.,333.,334.,335.,336.,337.,338.,339.,340.,341.,342.,343.,344.,345.,346.,347.,348.,349.,350.,351.,352.,353.,354.,355.,356.,357.,358.,359.,360.,361.,362.,363.],32,2) ; vExp.setInfoOnComponents(['Comp1 [m]','Com2 [s^2]']) ; vExp.setInfoOnComponents(['Comp1 [m]','Com2 [s^2]']) ; vExp+=i*1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass

    def test3(self):
        """ This test is more advanced a same field is defined on CELLS for time steps 0, 2 and 4, and on NODES for time steps 1 and 3.
        So two time step series on the same field. No profile here neither on cells nor on nodes.
        """
        fname="ForMEDReader3.med"
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
        fieldName="zeField"
        fs=MEDFileFieldMultiTS()
        ##### Time step 0 on cells
        i=0
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        fs.pushBackTimeStep(f)
        ##### Time step 1 on nodes
        i=1
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(1300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs.pushBackTimeStep(f)
        ##### Time step 2 on cells
        i=2
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(2100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(2200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        fs.pushBackTimeStep(f)
        ##### Time step 3 on nodes
        i=3
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(3300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs.pushBackTimeStep(f)
        ##### Time step 4
        i=4
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(4100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkCoherency()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(4200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkCoherency()
        f.setFieldNoProfileSBT(fCell1)
        #
        fs.pushBackTimeStep(f)
        mm.write(fname,2)
        fs.write(fname,0)
        a0Exp=mm.getCoords().deepCpy()
        del m,m1,mm,fs,f,fCell0,fCell1
        ########## GO for reading in MEDReader, by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),2)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),2) # two time series here : one for the cells, the second one for the nodes
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[1]),1)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=allFMTSLeavesPerCommonSupport[0][1] # start with the cells
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        for i in xrange(1,3):
            self.assertTrue(fcscp.isDataSetSupportEqualToThePreviousOne(i,fields))
            pass
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,5,5,5,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,100,104,108,112,117,122,127,132])))
        self.assertTrue(a3.isEqual(DataArrayInt([2,0,1,2,1,2,2,2,0,2,3,4,2,4,5,2,5,3,2,6,7,2,7,8,2,8,6,2,9,10,2,10,11,2,11,9,2,12,13,2,13,14,2,14,15,2,15,12,2,16,17,2,17,18,2,18,19,2,19,16,2,20,21,2,21,22,2,22,23,2,23,20,2,24,25,2,25,26,2,26,27,2,27,24,2,28,29,2,29,30,2,30,31,2,31,28,3,0,1,2,3,3,4,5,3,6,7,8,3,9,10,11,4,12,13,14,15,4,16,17,18,19,4,20,21,22,23,4,24,25,26,27,4,28,29,30,31])))
        assert a4 is None
        assert a5 is None
        # for cells
        for i in xrange(3):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)# Second 0 is for cells
            f.loadArraysIfNecessary()
            self.assertEqual(f.getName(),"zeField")
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([200.,201.,202.,203.,204.,205.,206.,207.,208.,209.,210.,211.,212.,213.,214.,215.,216.,217.,218.,219.,220.,221.,222.,223.,224.,225.,226.,227.,228.,229.,230.,231.,232.,233.,234.,235.,236.,237.,238.,239.,240.,241.,242.,243.,244.,245.,246.,247.,248.,249.,250.,251.,252.,253.,254.,255.,256.,257.,258.,259.,260.,261.,262.,263.,100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,110.,111.,112.,113.,114.,115.,116.,117.],41,2) ; vExp.setInfoOnComponents(['Comp1 [m]','Com2 [s^2]']) ; vExp+=i*2000
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        # for nodes
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[1],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),1)
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        for i in xrange(1,2):
            self.assertTrue(fcscp.isDataSetSupportEqualToThePreviousOne(i,fields))
            pass
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,5,5,5,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,8,12,16,21,26,31,36])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,0,1,2,3,3,4,5,3,6,7,8,3,9,10,11,4,12,13,14,15,4,16,17,18,19,4,20,21,22,23,4,24,25,26,27,4,28,29,30,31])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        for i in xrange(2):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)# Second 0 is for cells
            f.loadArraysIfNecessary()
            self.assertEqual(f.getName(),"zeField")
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([300.,301.,302.,303.,304.,305.,306.,307.,308.,309.,310.,311.,312.,313.,314.,315.,316.,317.,318.,319.,320.,321.,322.,323.,324.,325.,326.,327.,328.,329.,330.,331.,332.,333.,334.,335.,336.,337.,338.,339.,340.,341.,342.,343.,344.,345.,346.,347.,348.,349.,350.,351.,352.,353.,354.,355.,356.,357.,358.,359.,360.,361.,362.,363.],32,2) ; vExp.setInfoOnComponents(['Comp1 [m]','Com2 [s^2]']) ; vExp.setInfoOnComponents(['Comp1 [m]','Com2 [s^2]']) ; vExp+=i*2000+1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass

    def test4(self):
        """ This test defines 3 fields on nodes on the same mesh. All of these fields have no profile.
        """
        fname="ForMEDReader4.med"
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
        mm=MEDFileUMesh() ; mm.setMeshes([m])
        #
        fieldName1="zeField1"
        fieldName2="zeField2"
        fieldName3="zeField3"
        fs1=MEDFileFieldMultiTS() ; fs2=MEDFileFieldMultiTS() ; fs3=MEDFileFieldMultiTS()
        ##### Time step 0
        i=0
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName1) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(0+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs1.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName2) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(100+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs2.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName3) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(200+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs3.pushBackTimeStep(f)
        ##### Time step 1
        i=1
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName1) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(0+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs1.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName2) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(100+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs2.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName3) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(200+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs3.pushBackTimeStep(f)
        ##### Time step 2
        i=2
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName1) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(0+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs1.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName2) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(100+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs2.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName3) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(200+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs3.pushBackTimeStep(f)
        ##### Time step 3
        i=3
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName1) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(0+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs1.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName2) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(100+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs2.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName3) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(200+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs3.pushBackTimeStep(f)
        ##### Time step 4
        i=4
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName1) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(0+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs1.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName2) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(100+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs2.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName3) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(200+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
        fNode.checkCoherency()
        f.setFieldNoProfileSBT(fNode)
        fs3.pushBackTimeStep(f)
        #
        mm.write(fname,2)
        fs1.write(fname,0) ; fs2.write(fname,0) ; fs3.write(fname,0)
        a0Exp=mm.getCoords().deepCpy()
        del m,mm,fs1,fs2,fs3,f,fNode
        ########## GO for reading in MEDReader, by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),3)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 3 fields are defined on the same time steps
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),3)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),3)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),3)
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        for i in xrange(1,5):
            self.assertTrue(fcscp.isDataSetSupportEqualToThePreviousOne(i,fields))
            pass
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,5,5,5,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,8,12,16,21,26,31,36])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,0,1,2,3,3,4,5,3,6,7,8,3,9,10,11,4,12,13,14,15,4,16,17,18,19,4,20,21,22,23,4,24,25,26,27,4,28,29,30,31])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        # test all the time steps of the 1/1 time step serie, on field 1
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName1)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,46.,47.,48.,49.,50.,51.,52.,53.,54.,55.,56.,57.,58.,59.,60.,61.,62.,63.],32,2) ; vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]']) ; vExp+=i*1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        # test all the time steps of the 1/1 time step serie, on field 2
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName2)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,46.,47.,48.,49.,50.,51.,52.,53.,54.,55.,56.,57.,58.,59.,60.,61.,62.,63.],32,2) ; vExp.setInfoOnComponents(['Comp1_1 [m]','Com2_1 [s^2]']) ; vExp+=i*1000+100
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        # test all the time steps of the 1/1 time step serie, on field 3
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][2][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName3)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35.,36.,37.,38.,39.,40.,41.,42.,43.,44.,45.,46.,47.,48.,49.,50.,51.,52.,53.,54.,55.,56.,57.,58.,59.,60.,61.,62.,63.],32,2) ; vExp.setInfoOnComponents(['Comp1_2 [m]','Com2_2 [s^2]']) ; vExp+=i*1000+200
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass
    
    def test5(self):
        """ This test plays with profiles both cell profiles and node profiles. Two first fields (resp on cells and on nodes) lie on the same mesh support whereas the third
        mesh lies on a different mesh.
        """
        fname="ForMEDReader5.med"
        # building a mesh containing 6 tri3 + 5 quad4
        m=MEDCouplingUMesh("mesh",2)
        coords=DataArrayDouble([(0,0),(1,0),(2,0),(3,0),(4,0),(0,1),(1,1),(2,1),(3,1),(4,1),(0,2),(1,2),(2,2),(3,2),(4,2)]) ; coords.setInfoOnComponents(["XX [m]","YYY [km]"])
        m.setCoords(coords)
        m.allocateCells()
        m.insertNextCell(NORM_TRI3,[2,7,3]) ; m.insertNextCell(NORM_TRI3,[7,8,3]) ; m.insertNextCell(NORM_TRI3,[3,8,4]) ; m.insertNextCell(NORM_TRI3,[8,9,4])
        m.insertNextCell(NORM_TRI3,[13,9,8]) ; m.insertNextCell(NORM_TRI3,[13,14,9])
        m.insertNextCell(NORM_QUAD4,[0,5,6,1]) ; m.insertNextCell(NORM_QUAD4,[1,6,7,2]) ; m.insertNextCell(NORM_QUAD4,[5,10,11,6]) ; m.insertNextCell(NORM_QUAD4,[6,11,12,7])
        m.insertNextCell(NORM_QUAD4,[12,13,8,7])
        mm=MEDFileUMesh() ; mm.setMeshes([m])
        fam=DataArrayInt(11) ; fam.iota(0) ; mm.setFamilyFieldArr(0,fam) ; del fam
        num=DataArrayInt(11) ; num.iota(100) ; mm.setRenumFieldArr(0,num) ; del num
        #
        fieldName1="zeField1" ; pfl1=DataArrayInt([0,1,2,3,4,5]) ; pfl1.setName("pfl1") # on cells
        fieldName2="zeField2" ; pfl2=DataArrayInt([2,3,4,7,8,9,13,14]) ; pfl2.setName("pfl2") # on nodes
        fieldName3="zeField3" ; pfl3=DataArrayInt([0,1,2,3,4,5,9,10]) ; pfl3.setName("pfl3") # on cells but different support
        fs1=MEDFileFieldMultiTS() ; fs2=MEDFileFieldMultiTS() ; fs3=MEDFileFieldMultiTS()
        #
        for i in xrange(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName1)
            arr=DataArrayDouble(2*6) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl1)
            fs1.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName2)
            arr=DataArrayDouble(2*8) ; arr.iota(100+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl2)
            fs2.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName3)
            arr=DataArrayDouble(2*8) ; arr.iota(200+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl3)
            fs3.pushBackTimeStep(f)
            pass
        mm.write(fname,2)
        fs1.write(fname,0) ; fs2.write(fname,0) ; fs3.write(fname,0)
        a0Exp=mm.getCoords().deepCpy()
        del m,mm,fs1,fs2,fs3,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),3)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 3 fields are defined on the same time steps
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),3)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        ms[0].getDirectUndergroundSingleGeoTypeMeshes(0)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),2) # 2 support here
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[1][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        for i in xrange(1,5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp[pfl2].changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,5,5,5,5,5])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,8,12,16,20])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,0,3,1,3,3,4,1,3,1,4,2,3,4,5,2,3,6,5,4,3,6,7,5])))
        assert a4 is None
        assert a5 is None
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,1,2,3,4,5])))
        self.assertTrue(not a7)
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([100,101,102,103,104,105])))
        self.assertTrue(not a9)
        for i in xrange(5):
            nbOfT=[6,8]
            fieldNames=[fieldName1,fieldName2]
            for j in xrange(2):
                m={"i":j}
                f=allFMTSLeavesPerCommonSupport[0][0][j][i]
                fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
                f.loadArraysIfNecessary()
                v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
                self.assertEqual(f.getName(),fieldNames[j])
                self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
                vExp=DataArrayDouble(nbOfT[j]*2) ; vExp.iota(j*100+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_%(i)i [m]'%m,'Com2_%(i)i [s^2]'%m])
                self.assertTrue(v.isEqual(vExp,1e-12))
                pass
            pass
        # Let's go for the 2nd support
        fcscp=allFMTSLeavesPerCommonSupport[1][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        for i in xrange(1,5):
            self.assertTrue(fcscp.isDataSetSupportEqualToThePreviousOne(i,fields))
            pass
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,5,5,5,5,5,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,8,12,16,20,24,29])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,2,7,3,3,7,8,3,3,3,8,4,3,8,9,4,3,13,9,8,3,13,14,9,4,6,11,12,7,4,12,13,8,7])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,1,2,3,4,5,9,10])))
        self.assertTrue(not a7)
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([100,101,102,103,104,105,109,110])))
        self.assertTrue(not a9)
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[1][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),"zeField3")
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(8*2) ; vExp.iota(200+1000*i) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_2 [m]'%m,'Com2_2 [s^2]'%m])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass
    
    def test6(self):
        """ This test plays with cartesian mesh and profiles. When a sub cartesian mesh can also be considered as a cartesian mesh it is done.
        """
        fname="ForMEDReader6.med"
        m=MEDCouplingCMesh("mesh")
        coordsX=DataArrayDouble([0,1.1,2.2,3.3,4.4]) ; coordsX.setInfoOnComponents(["XX [m]"])
        coordsY=DataArrayDouble([0,1.7,3.4]) ; coordsY.setInfoOnComponents(["YYY [km]"])
        m.setCoords(coordsX,coordsY)
        mm=MEDFileCMesh() ; mm.setMesh(m)
        fam=DataArrayInt(8) ; fam.iota(0) ; mm.setFamilyFieldArr(0,fam) ; del fam
        num=DataArrayInt(8) ; num.iota(100) ; mm.setRenumFieldArr(0,num) ; del num
        #
        fieldName0="zeField0" ; # on cells
        fieldName1="zeField1" ; pfl1=DataArrayInt([2,3,6,7]) ; pfl1.setName("pfl1") # on cells
        fieldName2="zeField2" ; pfl2=DataArrayInt([2,3,4,7,8,9,12,13,14]) ; pfl2.setName("pfl2") # on nodes
        fieldName3="zeField3" ; pfl3=DataArrayInt([2,3,5,7]) ; pfl3.setName("pfl3") # on cells but different support
        fieldName4="zeField4" ;# on nodes
        fs0=MEDFileFieldMultiTS() ; fs1=MEDFileFieldMultiTS() ; fs2=MEDFileFieldMultiTS() ; fs3=MEDFileFieldMultiTS() ; fs4=MEDFileFieldMultiTS()
        #
        for i in xrange(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*8) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkCoherency()
            f.setFieldNoProfileSBT(fNode)
            fs0.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName1)
            arr=DataArrayDouble(2*4) ; arr.iota(100+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl1)
            self.assertEqual(pfl1.getName(),"pfl1")
            fs1.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName2)
            arr=DataArrayDouble(2*9) ; arr.iota(200+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl2)
            self.assertEqual(pfl2.getName(),"pfl2")
            fs2.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName3)
            arr=DataArrayDouble(2*4) ; arr.iota(300+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_3 [m]","Com2_3 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl3)
            self.assertEqual(pfl3.getName(),"pfl3")
            fs3.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName4) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*15) ; arr.iota(400+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_4 [m]","Com2_4 [s^2]"]) ; fNode.checkCoherency()
            f.setFieldNoProfileSBT(fNode)
            fs4.pushBackTimeStep(f)
            pass
        mm.write(fname,2)
        fs0.write(fname,0) ; fs1.write(fname,0) ; fs2.write(fname,0) ; fs3.write(fname,0) ; fs4.write(fname,0)
        del m,mm,fs1,fs2,fs3,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),5)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 5 fields are defined on the same time steps
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),5)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),3)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[1][0]),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[2][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        a,b=mml2.buildVTUArrays()
        self.assertTrue(a.isEqual(coordsX,1e-12))
        self.assertTrue(b.isEqual(coordsY,1e-12))
        self.assertTrue(isinstance(mml2,MEDCMeshMultiLev))
        for i in xrange(1,5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,1,2,3,4,5,6,7])))
        self.assertTrue(a7) # True because no copy
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([100,101,102,103,104,105,106,107])))
        self.assertTrue(a9) # True because no copy
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(8*2) ; vExp.iota(0+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName4)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(15*2) ; vExp.iota(400+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_4 [m]','Com2_4 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        
        fcscp=allFMTSLeavesPerCommonSupport[1][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDCMeshMultiLev)) # here the 2nd support is a part of CMesh that is also a CMesh -> CMesh not a UMesh
        a,b=mml2.buildVTUArrays()
        self.assertTrue(a.isEqual(coordsX[[2,3,4]],1e-12))
        self.assertTrue(b.isEqual(coordsY,1e-12))
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([2,3,6,7])))
        self.assertTrue(not a7) # False because copy
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([102,103,106,107])))
        self.assertTrue(not a9) # False because copy
        for i in xrange(1,5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[1][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName1)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(4*2) ; vExp.iota(100+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_1 [m]','Com2_1 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[1][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName2)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(9*2) ; vExp.iota(200+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_2 [m]','Com2_2 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        #
        fcscp=allFMTSLeavesPerCommonSupport[2][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev)) # here the 3rd support is a part of CMesh but impossible to simplify more than a UMesh
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        a0Exp=DataArrayDouble([0.,0.,1.1,0.,2.2,0.,3.3,0.,4.4,0.,0.,1.7,1.1,1.7,2.2,1.7,3.3,1.7,4.4,1.7,0.,3.4,1.1,3.4,2.2,3.4,3.3,3.4,4.4,3.4],15,2)
        a0Exp.setInfoOnComponents(["XX [m]","YYY [km]"])
        self.assertTrue(a0.isEqual(a0Exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5,10,15])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,3,2,7,8,4,4,3,8,9,4,7,6,11,12,4,9,8,13,14])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([2,3,5,7])))
        self.assertTrue(not a7) # False because copy
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([102,103,105,107])))
        self.assertTrue(not a9) # False because copy
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[2][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName3)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(4*2) ; vExp.iota(300+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_3 [m]','Com2_3 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass

    def test7(self):
        """ This test plays with curvilinear mesh and profiles. When a sub curvilinear mesh can also be considered as a cartesian mesh it is done.
        This test is very similar to the test6.
        """
        fname="ForMEDReader7.med"
        m=MEDCouplingCurveLinearMesh("mesh") ; m.setNodeGridStructure([5,3])
        a0Exp=DataArrayDouble([0.,0.,1.1,0.,2.2,0.,3.3,0.,4.4,0.,0.,1.7,1.1,1.7,2.2,1.7,3.3,1.7,4.4,1.7,0.,3.4,1.1,3.4,2.2,3.4,3.3,3.4,4.4,3.4],15,2)
        a0Exp.setInfoOnComponents(["XX [m]","YYY [km]"])
        m.setCoords(a0Exp)
        mm=MEDFileCurveLinearMesh() ; mm.setMesh(m)
        fam=DataArrayInt(8) ; fam.iota(0) ; mm.setFamilyFieldArr(0,fam) ; del fam
        num=DataArrayInt(8) ; num.iota(100) ; mm.setRenumFieldArr(0,num) ; del num
        #
        fieldName0="zeField0" ; # on cells
        fieldName1="zeField1" ; pfl1=DataArrayInt([2,3,6,7]) ; pfl1.setName("pfl1") # on cells
        fieldName2="zeField2" ; pfl2=DataArrayInt([2,3,4,7,8,9,12,13,14]) ; pfl2.setName("pfl2") # on nodes
        fieldName3="zeField3" ; pfl3=DataArrayInt([2,3,5,7]) ; pfl3.setName("pfl3") # on cells but different support
        fieldName4="zeField4" ;# on nodes
        fs0=MEDFileFieldMultiTS() ; fs1=MEDFileFieldMultiTS() ; fs2=MEDFileFieldMultiTS() ; fs3=MEDFileFieldMultiTS() ; fs4=MEDFileFieldMultiTS()
        #
        for i in xrange(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*8) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkCoherency()
            f.setFieldNoProfileSBT(fNode)
            fs0.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName1)
            arr=DataArrayDouble(2*4) ; arr.iota(100+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl1)
            self.assertEqual(pfl1.getName(),"pfl1")
            fs1.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName2)
            arr=DataArrayDouble(2*9) ; arr.iota(200+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl2)
            self.assertEqual(pfl2.getName(),"pfl2")
            fs2.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName3)
            arr=DataArrayDouble(2*4) ; arr.iota(300+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_3 [m]","Com2_3 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl3)
            self.assertEqual(pfl3.getName(),"pfl3")
            fs3.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName4) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*15) ; arr.iota(400+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_4 [m]","Com2_4 [s^2]"]) ; fNode.checkCoherency()
            f.setFieldNoProfileSBT(fNode)
            fs4.pushBackTimeStep(f)
            pass
        mm.write(fname,2)
        fs0.write(fname,0) ; fs1.write(fname,0) ; fs2.write(fname,0) ; fs3.write(fname,0) ; fs4.write(fname,0)
        del m,mm,fs1,fs2,fs3,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),5)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 5 fields are defined on the same time steps
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),5)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),3)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[1][0]),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[2][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDCurveLinearMeshMultiLev))
        a,b=mml2.buildVTUArrays()
        self.assertTrue(a.isEqual(a0Exp,1e-12))
        self.assertEqual(b,[5,3])
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,1,2,3,4,5,6,7])))
        self.assertTrue(a7) # True because no copy
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([100,101,102,103,104,105,106,107])))
        self.assertTrue(a9) # True because no copy
        for i in xrange(1,5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(8*2) ; vExp.iota(0+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName4)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(15*2) ; vExp.iota(400+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_4 [m]','Com2_4 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        #
        fcscp=allFMTSLeavesPerCommonSupport[1][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDCurveLinearMeshMultiLev)) # here the 2nd support is a part of CMesh that is also a CMesh -> CMesh not a UMesh
        a,b=mml2.buildVTUArrays()
        self.assertTrue(a.isEqual(a0Exp[pfl2],1e-12))
        self.assertEqual(b,[3,3])
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([2,3,6,7])))
        self.assertTrue(not a7) # False because copy
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([102,103,106,107])))
        self.assertTrue(not a9) # False because copy
        for i in xrange(1,5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[1][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName1)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(4*2) ; vExp.iota(100+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_1 [m]','Com2_1 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[1][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName2)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(9*2) ; vExp.iota(200+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_2 [m]','Com2_2 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        #
        fcscp=allFMTSLeavesPerCommonSupport[2][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev)) # here the 3rd support is a part of CMesh but impossible to simplify more than a UMesh
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        a0Exp=DataArrayDouble([0.,0.,1.1,0.,2.2,0.,3.3,0.,4.4,0.,0.,1.7,1.1,1.7,2.2,1.7,3.3,1.7,4.4,1.7,0.,3.4,1.1,3.4,2.2,3.4,3.3,3.4,4.4,3.4],15,2)
        a0Exp.setInfoOnComponents(["XX [m]","YYY [km]"])
        self.assertTrue(a0.isEqual(a0Exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5,10,15])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,3,2,7,8,4,4,3,8,9,4,7,6,11,12,4,9,8,13,14])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([2,3,5,7])))
        self.assertTrue(not a7) # False because copy
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([102,103,105,107])))
        self.assertTrue(not a9) # False because copy
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[2][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName3)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(4*2) ; vExp.iota(300+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_3 [m]','Com2_3 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass

    def test8(self):
        """ This test plays with with gauss fields with no profiles.
        """
        fname="ForMEDReader8.med"
        # building a mesh containing 6 tri3 + 5 quad4
        m=MEDCouplingUMesh("mesh",2)
        coords=DataArrayDouble([(0,0),(1,0),(2,0),(3,0),(4,0),(0,1),(1,1),(2,1),(3,1),(4,1),(0,2),(1,2),(2,2),(3,2),(4,2)]) ; coords.setInfoOnComponents(["XX [m]","YYY [km]"])
        m.setCoords(coords)
        m.allocateCells()
        m.insertNextCell(NORM_TRI3,[2,7,3]) ; m.insertNextCell(NORM_TRI3,[7,8,3]) ; m.insertNextCell(NORM_TRI3,[3,8,4]) ; m.insertNextCell(NORM_TRI3,[8,9,4])
        m.insertNextCell(NORM_TRI3,[13,9,8]) ; m.insertNextCell(NORM_TRI3,[13,14,9])
        m.insertNextCell(NORM_QUAD4,[0,5,6,1]) ; m.insertNextCell(NORM_QUAD4,[1,6,7,2]) ; m.insertNextCell(NORM_QUAD4,[5,10,11,6]) ; m.insertNextCell(NORM_QUAD4,[6,11,12,7])
        m.insertNextCell(NORM_QUAD4,[12,13,8,7])
        mm=MEDFileUMesh() ; mm.setMeshes([m])
        #
        fieldName0="zeField0"
        fieldName1="zeField1"
        fieldName2="zeField2"
        fieldName3="zeField3"
        fs0=MEDFileFieldMultiTS() ; fs1=MEDFileFieldMultiTS() ; fs2=MEDFileFieldMultiTS() ; fs3=MEDFileFieldMultiTS()
        for i in xrange(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_GAUSS_NE) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*38) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkCoherency()
            f.setFieldNoProfileSBT(fNode)
            fs0.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName1) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*11) ; arr.iota(100+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"]) ; fNode.checkCoherency()
            f.setFieldNoProfileSBT(fNode)
            fs1.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName2) ; fNode.setMesh(m)
            fNode.setGaussLocalizationOnCells([0,1,2,3],[0.,0.,1.,0.,0.,1.],[0.5,0.5,0.7,0.7],[0.8,0.2])
            fNode.setGaussLocalizationOnCells([4,5],[0.,0.,1.,0.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3],[0.8,0.05,0.1,0.04,0.01])
            fNode.setGaussLocalizationOnCells([6,7,8],[0.,0.,1.,0.,1.,1.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2],[0.8,0.05,0.1,0.04])
            fNode.setGaussLocalizationOnCells([9,10],[0.,0.,1.,0.,1.,1.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.8,0.8],[0.8,0.05,0.1,0.01,0.02,0.005,0.005])
            arr=DataArrayDouble(2*(4*2+2*5+3*4+2*7)) ; arr.iota(300+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"]) ; fNode.checkCoherency()
            f.setFieldNoProfileSBT(fNode)
            fs2.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName3) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*15) ; arr.iota(400+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_3 [m]","Com2_3 [s^2]"]) ; fNode.checkCoherency()
            f.setFieldNoProfileSBT(fNode)
            fs3.pushBackTimeStep(f)
            #
            pass
        #
        mm.write(fname,2)
        fs0.write(fname,0) ; fs1.write(fname,0) ; fs2.write(fname,0) ; fs3.write(fname,0)
        a0Exp=mm.getCoords().deepCpy()
        del m,mm,fs1,fs2,fs3,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),4)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 4 fields are defined on the same time steps
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),4)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),4)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,5,5,5,5,5,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,8,12,16,20,24,29,34,39,44])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,2,7,3,3,7,8,3,3,3,8,4,3,8,9,4,3,13,9,8,3,13,14,9,4,0,5,6,1,4,1,6,7,2,4,5,10,11,6,4,6,11,12,7,4,12,13,8,7])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        for i in xrange(1,5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(38*2) ; vExp.iota(0+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName1)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(11*2) ; vExp.iota(100+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_1 [m]','Com2_1 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[0][0][2][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName2)
            #self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer()) # not a bug
            vExp=DataArrayDouble(44*2) ; vExp.iota(300+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_2 [m]','Com2_2 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[0][0][3][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName3)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(15*2) ; vExp.iota(400+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_3 [m]','Com2_3 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        #
        pass

    def test9(self):
        """ This test plays with with gauss fields with profiles.
        """
        fname="ForMEDReader9.med"
        # building a mesh containing 6 tri3 + 5 quad4
        m=MEDCouplingUMesh("mesh",2)
        coords=DataArrayDouble([(0,0),(1,0),(2,0),(3,0),(4,0),(0,1),(1,1),(2,1),(3,1),(4,1),(0,2),(1,2),(2,2),(3,2),(4,2)]) ; coords.setInfoOnComponents(["XX [m]","YYY [km]"])
        m.setCoords(coords)
        m.allocateCells()
        m.insertNextCell(NORM_TRI3,[2,7,3]) ; m.insertNextCell(NORM_TRI3,[7,8,3]) ; m.insertNextCell(NORM_TRI3,[3,8,4]) ; m.insertNextCell(NORM_TRI3,[8,9,4])
        m.insertNextCell(NORM_TRI3,[13,9,8]) ; m.insertNextCell(NORM_TRI3,[13,14,9])
        m.insertNextCell(NORM_QUAD4,[0,5,6,1]) ; m.insertNextCell(NORM_QUAD4,[1,6,7,2]) ; m.insertNextCell(NORM_QUAD4,[5,10,11,6]) ; m.insertNextCell(NORM_QUAD4,[6,11,12,7])
        m.insertNextCell(NORM_QUAD4,[12,13,8,7])
        mm=MEDFileUMesh() ; mm.setMeshes([m])
        #
        fieldName0="zeField0"
        fieldName1="zeField1"
        fieldName2="zeField2"
        fieldName3="zeField3"
        pfl1=DataArrayInt([0,1,7,9,10]) ; pfl1.setName("pfl1") # on cells
        pfl2=DataArrayInt([1,2,3,6,7,8,11,12,13]) ; pfl2.setName("pfl2") # on nodes
        fs0=MEDFileFieldMultiTS() ; fs1=MEDFileFieldMultiTS() ; fs2=MEDFileFieldMultiTS() ; fs3=MEDFileFieldMultiTS()
        for i in xrange(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_GAUSS_NE) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0)
            arr=DataArrayDouble(2*18) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl1)
            fs0.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName1)
            arr=DataArrayDouble(2*5) ; arr.iota(100+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl1)
            fs1.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName2) ; fNode.setMesh(m[pfl1])
            fNode.setGaussLocalizationOnCells([0],[0.,0.,1.,0.,0.,1.],[0.5,0.5,0.7,0.7],[0.8,0.2])
            fNode.setGaussLocalizationOnCells([1],[0.,0.,1.,0.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3],[0.8,0.05,0.1,0.04,0.01])
            fNode.setGaussLocalizationOnCells([2,3],[0.,0.,1.,0.,1.,1.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2],[0.8,0.05,0.1,0.04])
            fNode.setGaussLocalizationOnCells([4],[0.,0.,1.,0.,1.,1.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.8,0.8],[0.8,0.05,0.1,0.01,0.02,0.005,0.005])
            arr=DataArrayDouble(2*(2*1+5*1+4*2+7*1)) ; arr.iota(300+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"]) ; fNode.checkCoherency()
            f.setFieldProfile(fNode,mm,0,pfl1)
            fs2.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName3)
            arr=DataArrayDouble(2*9) ; arr.iota(400+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_3 [m]","Com2_3 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl2)
            fs3.pushBackTimeStep(f)
            #
            pass
        #
        mm.write(fname,2)
        fs0.write(fname,0) ; fs1.write(fname,0) ; fs2.write(fname,0) ; fs3.write(fname,0)
        a0Exp=mm.getCoords().deepCpy()
        del m,mm,fs1,fs2,fs3,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),4)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 4 fields are defined on the same time steps
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),4)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),4)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp[pfl2].changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,5,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,8,13,18])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,1,4,2,3,4,5,2,4,0,3,4,1,4,3,6,7,4,4,7,8,5,4])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        for i in xrange(1,5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(18*2) ; vExp.iota(0+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName1)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(5*2) ; vExp.iota(100+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_1 [m]','Com2_1 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[0][0][2][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName2)
            #self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer()) # not a bug
            vExp=DataArrayDouble(22*2) ; vExp.iota(300+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_2 [m]','Com2_2 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[0][0][3][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName3)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(9*2) ; vExp.iota(400+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_3 [m]','Com2_3 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass
    
    def test10(self):
        """ This test plays with fields only on nodes containing profiles.
        """
        fname="ForMEDReader10.med"
        # building a mesh containing 6 tri3 + 5 quad4
        m=MEDCouplingUMesh("mesh",2)
        coords=DataArrayDouble([(0,0),(1,0),(2,0),(3,0),(4,0),(0,1),(1,1),(2,1),(3,1),(4,1),(0,2),(1,2),(2,2),(3,2),(4,2)]) ; coords.setInfoOnComponents(["XX [m]","YYY [km]"])
        m.setCoords(coords)
        m.allocateCells()
        m.insertNextCell(NORM_TRI3,[2,7,3]) ; m.insertNextCell(NORM_TRI3,[7,8,3]) ; m.insertNextCell(NORM_TRI3,[3,8,4]) ; m.insertNextCell(NORM_TRI3,[8,9,4])
        m.insertNextCell(NORM_TRI3,[13,9,8]) ; m.insertNextCell(NORM_TRI3,[13,14,9])
        m.insertNextCell(NORM_QUAD4,[0,5,6,1]) ; m.insertNextCell(NORM_QUAD4,[1,6,7,2]) ; m.insertNextCell(NORM_QUAD4,[5,10,11,6]) ; m.insertNextCell(NORM_QUAD4,[6,11,12,7])
        m.insertNextCell(NORM_QUAD4,[12,13,8,7])
        mm=MEDFileUMesh() ; mm.setMeshes([m])
        #
        fieldName0="zeField0"
        fieldName1="zeField1"
        fieldName2="zeField2"
        pfl1=DataArrayInt([1,2,3,6,7,8,11,12,13]) ; pfl1.setName("pfl1") # on nodes
        fs0=MEDFileFieldMultiTS() ; fs1=MEDFileFieldMultiTS() ; fs2=MEDFileFieldMultiTS()
        for i in xrange(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0)
            arr=DataArrayDouble(2*9) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl1)
            fs0.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName1)
            arr=DataArrayDouble(2*9) ; arr.iota(100+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl1)
            fs1.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName2)
            arr=DataArrayDouble(2*9) ; arr.iota(200+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
            f.setFieldProfile(fNode,mm,0,pfl1)
            fs2.pushBackTimeStep(f)
            #
            pass
        #
        mm.write(fname,2)
        fs0.write(fname,0) ; fs1.write(fname,0) ; fs2.write(fname,0)
        a0Exp=mm.getCoords().deepCpy()
        del m,mm,fs1,fs2,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),3)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 4 fields are defined on the same time steps
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),3)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),3)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp[pfl1].changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,5,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,8,13,18])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,1,4,2,3,4,5,2,4,0,3,4,1,4,3,6,7,4,4,7,8,5,4])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        for i in xrange(1,5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            #self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer()) # not a bug
            vExp=DataArrayDouble(9*2) ; vExp.iota(0+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName1)
            #self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer()) # not a bug
            vExp=DataArrayDouble(9*2) ; vExp.iota(100+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_1 [m]','Com2_1 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[0][0][2][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName2)
            #self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer()) # not a bug
            vExp=DataArrayDouble(9*2) ; vExp.iota(200+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_2 [m]','Com2_2 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass
    
    def test11(self):
        """ This test is the ultimate test for the profiles with gauss points. It tests that even if there is non contiguous parts in definition of gauss points, it works !
        WARNING here, as no other discretizations exists, the priority is given to the field -> the mesh is renumbered to accelerate the build of array of field.
        """
        fname="ForMEDReader11.med"
        m=MEDCouplingCMesh("mesh")
        arr=DataArrayDouble(5) ; arr.iota()
        m.setCoords(arr,arr)
        m=m.buildUnstructured() ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        mm=MEDFileUMesh() ; mm.setMeshes([m])
        #
        fieldName0="zeField0"
        fs0=MEDFileFieldMultiTS()
        for i in xrange(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0) ; fNode.setMesh(m)
            fNode.setGaussLocalizationOnCells([0,2,3,4,7,15],[0.,0.,1.,0.,1.,1.,0.,1.],[0.5,0.5,0.7,0.7],[0.8,0.2])
            fNode.setGaussLocalizationOnCells([1,5,8,9],[0.,0.,1.,0.,1.,1.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3],[0.8,0.05,0.1,0.04,0.01])
            fNode.setGaussLocalizationOnCells([6,10,13],[0.,0.,1.,0.,1.,1.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2],[0.8,0.05,0.1,0.04])
            fNode.setGaussLocalizationOnCells([11,12,14],[0.,0.,1.,0.,1.,1.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.8,0.8],[0.8,0.05,0.1,0.01,0.02,0.005,0.005])
            arr=DataArrayDouble(2*(2*6+5*4+4*3+7*3)) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkCoherency()
            f.setFieldNoProfileSBT(fNode)
            fs0.pushBackTimeStep(f)
            pass
        mm.write(fname,2)
        fs0.write(fname,0)
        a0Exp=mm.getCoords().deepCpy()
        del m,mm,fs0,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),1)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 1 field is defined on the same time steps
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),1)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,1,0,5,6,4,3,2,7,8,4,4,3,8,9,4,6,5,10,11,4,9,8,13,14,4,19,18,23,24,4,2,1,6,7,4,7,6,11,12,4,11,10,15,16,4,12,11,16,17,4,8,7,12,13,4,13,12,17,18,4,17,16,21,22,4,14,13,18,19,4,16,15,20,21,4,18,17,22,23]))) # <- here the mesh is renumbered : the mesh is equal to m[[0,2,3,4,7,15, 1,5,8,9, 6,10,13, 11,12,14]]
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        for i in xrange(1,5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([0.,1.,2.,3.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,44.,45.,46.,47.,126.,127.,128.,129.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,26.,27.,28.,29.,30.,31.,32.,33.,34.,35.,48.,49.,50.,51.,52.,53.,54.,55.,56.,57.,58.,59.,60.,61.,62.,63.,64.,65.,66.,67.,36.,37.,38.,39.,40.,41.,42.,43.,68.,69.,70.,71.,72.,73.,74.,75.,104.,105.,106.,107.,108.,109.,110.,111.,76.,77.,78.,79.,80.,81.,82.,83.,84.,85.,86.,87.,88.,89.,90.,91.,92.,93.,94.,95.,96.,97.,98.,99.,100.,101.,102.,103.,112.,113.,114.,115.,116.,117.,118.,119.,120.,121.,122.,123.,124.,125.],65,2) ; vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
            vExp+=i*1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass

    def test12(self):
        """ This test is the second ultimate test for the profiles with gauss points.
        This test is close to test11 but here a 2nd field on cells without profile. So here the mesh is expected to be the same than m.
        """
        fname="ForMEDReader12.med"
        m=MEDCouplingCMesh("mesh")
        arr=DataArrayDouble(5) ; arr.iota()
        m.setCoords(arr,arr)
        m=m.buildUnstructured() ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        mm=MEDFileUMesh() ; mm.setMeshes([m])
        #
        fieldName0="zeField0"
        fieldName1="zeField1"
        fs0=MEDFileFieldMultiTS() ; fs1=MEDFileFieldMultiTS()
        for i in xrange(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0) ; fNode.setMesh(m)
            fNode.setGaussLocalizationOnCells([0,2,3,4,7,15],[0.,0.,1.,0.,1.,1.,0.,1.],[0.5,0.5,0.7,0.7],[0.8,0.2])
            fNode.setGaussLocalizationOnCells([1,5,8,9],[0.,0.,1.,0.,1.,1.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3],[0.8,0.05,0.1,0.04,0.01])
            fNode.setGaussLocalizationOnCells([6,10,13],[0.,0.,1.,0.,1.,1.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2],[0.8,0.05,0.1,0.04])
            fNode.setGaussLocalizationOnCells([11,12,14],[0.,0.,1.,0.,1.,1.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.8,0.8],[0.8,0.05,0.1,0.01,0.02,0.005,0.005])
            arr=DataArrayDouble(2*(2*6+5*4+4*3+7*3)) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkCoherency()
            f.setFieldNoProfileSBT(fNode)
            fs0.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName1) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*16) ; arr.iota(300+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"]) ; fNode.checkCoherency()
            f.setFieldNoProfileSBT(fNode)
            fs1.pushBackTimeStep(f)
            pass
        mm.write(fname,2)
        fs0.write(fname,0) ; fs1.write(fname,0)
        a0Exp=mm.getCoords().deepCpy()
        del m,mm,fs0,fs1,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),2)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 2 fields are defined on the same time steps
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),2)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),2)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,1,0,5,6,4,2,1,6,7,4,3,2,7,8,4,4,3,8,9,4,6,5,10,11,4,7,6,11,12,4,8,7,12,13,4,9,8,13,14,4,11,10,15,16,4,12,11,16,17,4,13,12,17,18,4,14,13,18,19,4,16,15,20,21,4,17,16,21,22,4,18,17,22,23,4,19,18,23,24]))) # <- here the mesh is NOT renumbered : the mesh is equal to m
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        for i in xrange(1,5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            #self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer()) # not a bug : huge reordering performed !
            vExp=DataArrayDouble(65*2) ; vExp.iota(0+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName1)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer()) # not a bug : huge reordering performed !
            vExp=DataArrayDouble(16*2) ; vExp.iota(300+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_1 [m]','Com2_1 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass

    def test13(self):
            """ Testing polyhedrons mixed with hexa8"""
            fname="ForMEDReader13.med"
            m=MEDCouplingUMesh("mesh",3)
            m.allocateCells()
            m.insertNextCell(NORM_HEXA8,[1,0,6,7,13,12,18,19]) ; m.insertNextCell(NORM_HEXA8,[2,1,7,8,14,13,19,20])
            m.insertNextCell(NORM_POLYHED,[3,2,8,9,-1,15,21,20,14,-1,3,15,14,2,-1,2,14,20,8,-1,8,20,21,9,-1,9,21,15,3])
            m.insertNextCell(NORM_POLYHED,[4,3,9,10,-1,16,22,21,15,-1,4,16,15,3,-1,3,15,21,9,-1,9,21,22,10,-1,10,22,16,4])
            m.insertNextCell(NORM_POLYHED,[5,4,10,11,-1,17,23,22,16,-1,5,17,16,4,-1,4,16,22,10,-1,10,22,23,11,-1,11,23,17,5])
            coords=DataArrayDouble([0.,0.,0.,1.,0.,0.,2.,0.,0.,3.,0.,0.,4.,0.,0.,5.,0.,0.,0.,1.,0.,1.,1.,0.,2.,1.,0.,3.,1.,0.,4.,1.,0.,5.,1.,0.,0.,0.,1.,1.,0.,1.,2.,0.,1.,3.,0.,1.,4.,0.,1.,5.,0.,1.,0.,1.,1.,1.,1.,1.,2.,1.,1.,3.,1.,1.,4.,1.,1.,5.,1.,1.],24,3) ; coords.setInfoOnComponents(["XX [m]","YYY [km]","ZZZZ [Mm]"])
            m.setCoords(coords)
            mm=MEDFileUMesh() ; mm.setMeshes([m])
            fs0=MEDFileFieldMultiTS() ; fs1=MEDFileFieldMultiTS() ; fs2=MEDFileFieldMultiTS() ; fs3=MEDFileFieldMultiTS()
            fieldName0="zeField0"
            fieldName1="zeField1"
            fieldName2="zeField2" ; pfl1=DataArrayInt([2,3]) ; pfl1.setName("pfl1")
            fieldName3="zefield3" ; pfl2=DataArrayInt([2,3,4]) ; pfl2.setName("pfl2")
            for i in xrange(5):
                f=MEDFileField1TS()
                fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
                fNode.setName(fieldName0) ; fNode.setMesh(m)
                arr=DataArrayDouble(2*5) ; arr.iota(0+1000*i) ; arr.rearrange(2)
                fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkCoherency()
                f.setFieldNoProfileSBT(fNode)
                fs0.pushBackTimeStep(f)
                #
                f=MEDFileField1TS()
                fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
                fNode.setName(fieldName1) ; fNode.setMesh(m)
                arr=DataArrayDouble(2*5) ; arr.iota(100+1000*i) ; arr.rearrange(2)
                fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"]) ; fNode.checkCoherency()
                f.setFieldNoProfileSBT(fNode)
                fs1.pushBackTimeStep(f)
                #
                f=MEDFileField1TS()
                fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
                fNode.setName(fieldName2) ; fNode.setMesh(m[pfl1])
                arr=DataArrayDouble(2*2) ; arr.iota(200+1000*i) ; arr.rearrange(2)
                fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"]) ; fNode.checkCoherency()
                f.setFieldProfile(fNode,mm,0,pfl1)
                fs2.pushBackTimeStep(f)
                #
                f=MEDFileField1TS()
                fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
                fNode.setName(fieldName3) ; fNode.setMesh(m[pfl2])
                arr=DataArrayDouble(2*3) ; arr.iota(300+1000*i) ; arr.rearrange(2)
                fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_3 [m]","Com2_3 [s^2]"]) ; fNode.checkCoherency()
                f.setFieldProfile(fNode,mm,0,pfl2)
                fs3.pushBackTimeStep(f)
                pass
            mm.write(fname,2)
            fs0.write(fname,0) ; fs1.write(fname,0) ; fs2.write(fname,0) ; fs3.write(fname,0)
            a0Exp=mm.getCoords().deepCpy()
            del m,mm,fs0
            ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
            ms=MEDFileMeshes(fname)
            fields=MEDFileFields(fname,False)
            fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
            allFMTSLeavesToDisplay=[]
            for fields in fields_per_mesh:
                allFMTSLeavesToDisplay2=[]
                for fmts in fields:
                    allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                    pass
                allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
                pass
            self.assertEqual(len(allFMTSLeavesToDisplay),1)
            self.assertEqual(len(allFMTSLeavesToDisplay[0]),4)
            allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
            self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 2 fields are defined on the same time steps
            self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),4)
            allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
            self.assertEqual(len(allFMTSLeavesPerCommonSupport),3)
            self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),2)
            self.assertEqual(len(allFMTSLeavesPerCommonSupport[1][0]),1)
            self.assertEqual(len(allFMTSLeavesPerCommonSupport[2][0]),1)
            #
            mst=MEDFileMeshStruct.New(ms[0])
            #
            fcscp=allFMTSLeavesPerCommonSupport[0][1]
            mml=fcscp.buildFromScratchDataSetSupport(0,fields)
            mml2=mml.prepare()
            self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
            ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
            self.assertTrue(ncc)
            self.assertTrue(a0.isEqual(a0Exp,1e-12))
            self.assertTrue(a1.isEqual(DataArrayByte([12,12,42,42,42])))
            self.assertTrue(a2.isEqual(DataArrayInt([0,9,18,27,36])))
            self.assertTrue(a3.isEqual(DataArrayInt([8,1,0,6,7,13,12,18,19,8,2,1,7,8,14,13,19,20,8,2,3,8,9,14,15,20,21,8,3,4,9,10,15,16,21,22,8,4,5,10,11,16,17,22,23])))
            self.assertTrue(a4.isEqual(DataArrayInt([-1,-1,0,31,62])))
            self.assertTrue(a5.isEqual(DataArrayInt([6,4,3,2,8,9,4,15,21,20,14,4,3,15,14,2,4,2,14,20,8,4,8,20,21,9,4,9,21,15,3,6,4,4,3,9,10,4,16,22,21,15,4,4,16,15,3,4,3,15,21,9,4,9,21,22,10,4,10,22,16,4,6,4,5,4,10,11,4,17,23,22,16,4,5,17,16,4,4,4,16,22,10,4,10,22,23,11,4,11,23,17,5])))
            for i in xrange(1,5):
                self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
                pass
            pass
            for i in xrange(5):
                f=allFMTSLeavesPerCommonSupport[0][0][0][i]
                fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
                f.loadArraysIfNecessary()
                v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
                self.assertEqual(f.getName(),fieldName0)
                self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
                vExp=DataArrayDouble(5*2) ; vExp.iota(0+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
                self.assertTrue(v.isEqual(vExp,1e-12))
                #
                f=allFMTSLeavesPerCommonSupport[0][0][1][i]
                fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
                f.loadArraysIfNecessary()
                v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
                self.assertEqual(f.getName(),fieldName1)
                self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
                vExp=DataArrayDouble(5*2) ; vExp.iota(100+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_1 [m]','Com2_1 [s^2]'])
                self.assertTrue(v.isEqual(vExp,1e-12))
                pass
            #
            fcscp=allFMTSLeavesPerCommonSupport[1][1]
            mml=fcscp.buildFromScratchDataSetSupport(0,fields)
            mml2=mml.prepare()
            self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
            ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
            self.assertTrue(ncc)
            self.assertTrue(a0.isEqual(a0Exp,1e-12))
            self.assertTrue(a1.isEqual(DataArrayByte([42,42])))
            self.assertTrue(a2.isEqual(DataArrayInt([0,9])))
            self.assertTrue(a3.isEqual(DataArrayInt([8,2,3,8,9,14,15,20,21,8,3,4,9,10,15,16,21,22])))
            self.assertTrue(a4.isEqual(DataArrayInt([0,31])))
            self.assertTrue(a5.isEqual(DataArrayInt([6,4,3,2,8,9,4,15,21,20,14,4,3,15,14,2,4,2,14,20,8,4,8,20,21,9,4,9,21,15,3,6,4,4,3,9,10,4,16,22,21,15,4,4,16,15,3,4,3,15,21,9,4,9,21,22,10,4,10,22,16,4])))
            for i in xrange(1,5):
                self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
                pass
            pass
            for i in xrange(5):
                f=allFMTSLeavesPerCommonSupport[1][0][0][i]
                fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
                f.loadArraysIfNecessary()
                v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
                self.assertEqual(f.getName(),fieldName2)
                self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
                vExp=DataArrayDouble(2*2) ; vExp.iota(200+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_2 [m]','Com2_2 [s^2]'])
                self.assertTrue(v.isEqual(vExp,1e-12))
                pass
            #
            fcscp=allFMTSLeavesPerCommonSupport[2][1]
            mml=fcscp.buildFromScratchDataSetSupport(0,fields)
            mml2=mml.prepare()
            self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
            ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
            self.assertTrue(ncc)
            self.assertTrue(a0.isEqual(a0Exp,1e-12))
            self.assertTrue(a1.isEqual(DataArrayByte([42,42,42])))
            self.assertTrue(a2.isEqual(DataArrayInt([0,9,18])))
            self.assertTrue(a3.isEqual(DataArrayInt([8,2,3,8,9,14,15,20,21,8,3,4,9,10,15,16,21,22,8,4,5,10,11,16,17,22,23])))
            self.assertTrue(a4.isEqual(DataArrayInt([0,31,62])))
            self.assertTrue(a5.isEqual(DataArrayInt([6,4,3,2,8,9,4,15,21,20,14,4,3,15,14,2,4,2,14,20,8,4,8,20,21,9,4,9,21,15,3,6,4,4,3,9,10,4,16,22,21,15,4,4,16,15,3,4,3,15,21,9,4,9,21,22,10,4,10,22,16,4,6,4,5,4,10,11,4,17,23,22,16,4,5,17,16,4,4,4,16,22,10,4,10,22,23,11,4,11,23,17,5])))
            for i in xrange(1,5):
                self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
                pass
            pass
            for i in xrange(5):
                f=allFMTSLeavesPerCommonSupport[2][0][0][i]
                fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
                f.loadArraysIfNecessary()
                v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
                self.assertEqual(f.getName(),fieldName3)
                self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
                vExp=DataArrayDouble(3*2) ; vExp.iota(300+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_3 [m]','Com2_3 [s^2]'])
                self.assertTrue(v.isEqual(vExp,1e-12))
                pass
            pass

    def test14(self):
            """ Testing only polyhedrons"""
            fname="ForMEDReader14.med"
            m=MEDCouplingUMesh("mesh",3)
            m.allocateCells()
            m.insertNextCell(NORM_POLYHED,[3,2,8,9,-1,15,21,20,14,-1,3,15,14,2,-1,2,14,20,8,-1,8,20,21,9,-1,9,21,15,3])
            m.insertNextCell(NORM_POLYHED,[4,3,9,10,-1,16,22,21,15,-1,4,16,15,3,-1,3,15,21,9,-1,9,21,22,10,-1,10,22,16,4])
            m.insertNextCell(NORM_POLYHED,[5,4,10,11,-1,17,23,22,16,-1,5,17,16,4,-1,4,16,22,10,-1,10,22,23,11,-1,11,23,17,5])
            coords=DataArrayDouble([0.,0.,0.,1.,0.,0.,2.,0.,0.,3.,0.,0.,4.,0.,0.,5.,0.,0.,0.,1.,0.,1.,1.,0.,2.,1.,0.,3.,1.,0.,4.,1.,0.,5.,1.,0.,0.,0.,1.,1.,0.,1.,2.,0.,1.,3.,0.,1.,4.,0.,1.,5.,0.,1.,0.,1.,1.,1.,1.,1.,2.,1.,1.,3.,1.,1.,4.,1.,1.,5.,1.,1.],24,3) ; coords.setInfoOnComponents(["XX [m]","YYY [km]","ZZZZ [Mm]"])
            m.setCoords(coords)
            mm=MEDFileUMesh() ; mm.setMeshes([m])
            fs0=MEDFileFieldMultiTS() ; fs1=MEDFileFieldMultiTS()
            fieldName0="zeField0"
            fieldName1="zeField1"
            for i in xrange(5):
                f=MEDFileField1TS()
                fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
                fNode.setName(fieldName0) ; fNode.setMesh(m)
                arr=DataArrayDouble(2*3) ; arr.iota(0+1000*i) ; arr.rearrange(2)
                fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkCoherency()
                f.setFieldNoProfileSBT(fNode)
                fs0.pushBackTimeStep(f)
                #
                f=MEDFileField1TS()
                fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
                fNode.setName(fieldName1) ; fNode.setMesh(m)
                arr=DataArrayDouble(2*3) ; arr.iota(100+1000*i) ; arr.rearrange(2)
                fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"]) ; fNode.checkCoherency()
                f.setFieldNoProfileSBT(fNode)
                fs1.pushBackTimeStep(f)
                pass
            mm.write(fname,2)
            fs0.write(fname,0) ; fs1.write(fname,0)
            a0Exp=mm.getCoords().deepCpy()
            del m,mm,fs0
            ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
            ms=MEDFileMeshes(fname)
            fields=MEDFileFields(fname,False)
            fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
            allFMTSLeavesToDisplay=[]
            for fields in fields_per_mesh:
                allFMTSLeavesToDisplay2=[]
                for fmts in fields:
                    allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                    pass
                allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
                pass
            self.assertEqual(len(allFMTSLeavesToDisplay),1)
            self.assertEqual(len(allFMTSLeavesToDisplay[0]),2)
            allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
            self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 2 fields are defined on the same time steps
            self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),2)
            allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
            self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
            self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),2)
            #
            mst=MEDFileMeshStruct.New(ms[0])
            #
            fcscp=allFMTSLeavesPerCommonSupport[0][1]
            mml=fcscp.buildFromScratchDataSetSupport(0,fields)
            mml2=mml.prepare()
            self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
            ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
            self.assertTrue(ncc)
            self.assertTrue(a0.isEqual(a0Exp,1e-12))
            self.assertTrue(a1.isEqual(DataArrayByte([42,42,42])))
            self.assertTrue(a2.isEqual(DataArrayInt([0,9,18])))
            self.assertTrue(a3.isEqual(DataArrayInt([8,2,3,8,9,14,15,20,21,8,3,4,9,10,15,16,21,22,8,4,5,10,11,16,17,22,23])))
            self.assertTrue(a4.isEqual(DataArrayInt([0,31,62])))
            self.assertTrue(a5.isEqual(DataArrayInt([6,4,3,2,8,9,4,15,21,20,14,4,3,15,14,2,4,2,14,20,8,4,8,20,21,9,4,9,21,15,3,6,4,4,3,9,10,4,16,22,21,15,4,4,16,15,3,4,3,15,21,9,4,9,21,22,10,4,10,22,16,4,6,4,5,4,10,11,4,17,23,22,16,4,5,17,16,4,4,4,16,22,10,4,10,22,23,11,4,11,23,17,5])))
            for i in xrange(1,5):
                self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
                pass
            a6,a7=mml2.retrieveFamilyIdsOnCells()
            self.assertTrue(a6.isEqual(DataArrayInt([0,0,0])))
            self.assertTrue(a7)
            a8,a9=mml2.retrieveNumberIdsOnCells()
            self.assertTrue(a8 is None)
            self.assertTrue(a9)
            for i in xrange(5):
                f=allFMTSLeavesPerCommonSupport[0][0][0][i]
                fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
                f.loadArraysIfNecessary()
                v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
                self.assertEqual(f.getName(),fieldName0)
                self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
                vExp=DataArrayDouble(3*2) ; vExp.iota(0+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
                self.assertTrue(v.isEqual(vExp,1e-12))
                #
                f=allFMTSLeavesPerCommonSupport[0][0][1][i]
                fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
                f.loadArraysIfNecessary()
                v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
                self.assertEqual(f.getName(),fieldName1)
                self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
                vExp=DataArrayDouble(3*2) ; vExp.iota(100+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_1 [m]','Com2_1 [s^2]'])
                self.assertTrue(v.isEqual(vExp,1e-12))
                pass
            pass

    def test15(self):
        """
        "ForMEDReader15.med" file has a spaceDim 3 mesh "mesh" (it is important !)
        and a field "zeField" lying on a single geometric type for Cell discr and node part.
        Test that can appear the most simple but it hides a big issue of MEDReader
        that copies are reduced at most. So it can leads to SIGSEGV if the memory management is not OK for int* and double * similar between VTK and MEDCoupling.
        """
        fname="ForMEDReader15.med"
        m0=MEDCouplingCMesh()
        arr=DataArrayDouble(3) ; arr.iota(0)
        m0.setCoords(arr,arr,arr)
        m0.setName("mesh")
        m0=m0.buildUnstructured()
        #
        fieldName="zeField"
        fCell=MEDCouplingFieldDouble(ON_CELLS)
        fCell.setName(fieldName)
        fCell.setMesh(m0)
        #
        fNode=MEDCouplingFieldDouble(ON_NODES)
        fNode.setName(fieldName)
        fNode.setMesh(m0)
        #
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m0)
        fam=DataArrayInt(8) ; fam.iota(0) ; mm.setFamilyFieldArr(0,fam) ; del fam
        num=DataArrayInt(8) ; num.iota(100) ; mm.setRenumFieldArr(0,num) ; del num
        #
        ffs=MEDFileFieldMultiTS()
        # TimeStep 0
        t=(1.,0,0) ; off=0.
        f1ts=MEDFileField1TS()
        a=DataArrayDouble(m0.getNumberOfCells()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        fCell.setArray(a)
        fCell.setTime(*t)
        fCell.checkCoherency()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode.setArray(a)
        fNode.setTime(*t)
        fNode.checkCoherency()
        f1ts.setFieldNoProfileSBT(fCell)
        f1ts.setFieldNoProfileSBT(fNode)
        ffs.pushBackTimeStep(f1ts)
        # TimeStep 1
        t=(2.1,1,0) ; off=100.
        f1ts=MEDFileField1TS()
        a=DataArrayDouble(m0.getNumberOfCells()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        fCell.setArray(a)
        fCell.setTime(*t)
        fCell.checkCoherency()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode.setArray(a)
        fNode.setTime(*t)
        fNode.checkCoherency()
        f1ts.setFieldNoProfileSBT(fCell)
        f1ts.setFieldNoProfileSBT(fNode)
        ffs.pushBackTimeStep(f1ts)
        # TimeStep 2
        t=(3.2,2,0) ; off=200.
        f1ts=MEDFileField1TS()
        a=DataArrayDouble(m0.getNumberOfCells()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        fCell.setArray(a)
        fCell.setTime(*t)
        fCell.checkCoherency()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode.setArray(a)
        fNode.setTime(*t)
        fNode.checkCoherency()
        f1ts.setFieldNoProfileSBT(fCell)
        f1ts.setFieldNoProfileSBT(fNode)
        ffs.pushBackTimeStep(f1ts)
        # TimeStep 3
        t=(4.3,3,1) ; off=300.
        f1ts=MEDFileField1TS()
        a=DataArrayDouble(m0.getNumberOfCells()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        fCell.setArray(a)
        fCell.setTime(*t)
        fCell.checkCoherency()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode.setArray(a)
        fNode.setTime(*t)
        fNode.checkCoherency()
        f1ts.setFieldNoProfileSBT(fCell)
        f1ts.setFieldNoProfileSBT(fNode)
        ffs.pushBackTimeStep(f1ts)
        #
        mm.write(fname,2)
        ffs.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),2)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 2 fields are defined on the same time steps
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),2)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),2)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(ncc)
        self.assertTrue(a0.isEqual(DataArrayDouble([0.,0.,0.,1.,0.,0.,2.,0.,0.,0.,1.,0.,1.,1.,0.,2.,1.,0.,0.,2.,0.,1.,2.,0.,2.,2.,0.,0.,0.,1.,1.,0.,1.,2.,0.,1.,0.,1.,1.,1.,1.,1.,2.,1.,1.,0.,2.,1.,1.,2.,1.,2.,2.,1.,0.,0.,2.,1.,0.,2.,2.,0.,2.,0.,1.,2.,1.,1.,2.,2.,1.,2.,0.,2.,2.,1.,2.,2.,2.,2.,2.0],27,3),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([12,12,12,12,12,12,12,12])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,9,18,27,36,45,54,63])))
        self.assertTrue(a3.isEqual(DataArrayInt([8,1,0,3,4,10,9,12,13,8,2,1,4,5,11,10,13,14,8,4,3,6,7,13,12,15,16,8,5,4,7,8,14,13,16,17,8,10,9,12,13,19,18,21,22,8,11,10,13,14,20,19,22,23,8,13,12,15,16,22,21,24,25,8,14,13,16,17,23,22,25,26])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,1,2,3,4,5,6,7])))
        self.assertTrue(a7) # no copy here
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([100,101,102,103,104,105,106,107])))
        self.assertTrue(a9) # no copy here
        pass

    def test16(self):
        """ Here 2 meshes "mesh1" and "mesh2" and 4 fields (no profiles here) :
        - "zeField1_0" (CELLS) and "zeField2_0" (NODES) on "mesh1"
        - "zeField3_1" (CELLS) and "zeField4_1" (NODES) on "mesh2"
        time steps series are the same for the whole 4 fields
        """
        fname="ForMEDReader16.med"
        m0=MEDCouplingCMesh()
        arr=DataArrayDouble(3) ; arr.iota(0)
        m0.setCoords(arr,arr,arr)
        m0.setName("mesh1")
        m0=m0.buildUnstructured()
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS)
        fCell1.setName("zeField1_0")
        fCell1.setMesh(m0)
        #
        fNode1=MEDCouplingFieldDouble(ON_NODES)
        fNode1.setName("zeField2_0")
        fNode1.setMesh(m0)
        #
        mms=MEDFileMeshes()
        mm1=MEDFileUMesh()
        mm1.setMeshAtLevel(0,m0)
        fam=DataArrayInt([0,1,0,1,2,3,2,3]); mm1.setFamilyFieldArr(0,fam) ; del fam
        num=DataArrayInt(8) ; num.iota(100) ; mm1.setRenumFieldArr(0,num) ; del num
        mm1.setFamilyId("FAMILLE_ZERO",0) ; mm1.setFamilyId("Family1_1",1) ; mm1.setFamilyId("Family1_2",2) ; mm1.setFamilyId("Family1_3",3) ; mm1.setFamilyId("Family1_4",4)
        mm1.setFamiliesIdsOnGroup("Grp1_1",[0,1]) ; mm1.setFamiliesIdsOnGroup("Grp1_2",[2,3])
        mms.pushMesh(mm1) ; del mm1
        #
        m1=m0.deepCpy() ; m1.translate([2.5,0.,0.]) ; m1.setName("mesh2")
        #
        fCell2=MEDCouplingFieldDouble(ON_CELLS)
        fCell2.setName("zeField3_1")
        fCell2.setMesh(m1)
        #
        fNode2=MEDCouplingFieldDouble(ON_NODES)
        fNode2.setName("zeField4_1")
        fNode2.setMesh(m1)
        #
        mm2=MEDFileUMesh()
        mm2.setMeshAtLevel(0,m1)
        fam=DataArrayInt([0,1,0,1,2,3,2,3]); mm2.setFamilyFieldArr(0,fam) ; del fam
        num=DataArrayInt(8) ; num.iota(200) ; mm2.setRenumFieldArr(0,num) ; del num
        mm2.setFamilyId("FAMILLE_ZERO",0) ; mm2.setFamilyId("Family2_1",1) ; mm2.setFamilyId("Family2_2",2) ; mm2.setFamilyId("Family2_3",3) ; mm2.setFamilyId("Family2_4",4)
        mm2.setFamiliesIdsOnGroup("Grp2_1",[0,1]) ; mm2.setFamiliesIdsOnGroup("Grp2_2",[2,3]) ; mm2.setFamiliesIdsOnGroup("Grp2_3",[1,2,3])
        mms.pushMesh(mm2) ; del mm2
        ffs1_1=MEDFileFieldMultiTS()
        ffs1_2=MEDFileFieldMultiTS()
        ffs2_1=MEDFileFieldMultiTS()
        ffs2_2=MEDFileFieldMultiTS()
        mts=MEDFileFields()
        for elt in ffs1_1,ffs1_2,ffs2_1,ffs2_2:
            mts.pushField(elt)
            pass
        # TimeStep 0
        t=(1.,0,0) ; off=0.
        f1ts1=MEDFileField1TS()
        f1ts2=MEDFileField1TS()
        a=DataArrayDouble(m0.getNumberOfCells()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        fCell1.setArray(a)
        fCell1.setTime(*t)
        fCell1.checkCoherency()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode1.setArray(a)
        fNode1.setTime(*t)
        fNode1.checkCoherency()
        f1ts1.setFieldNoProfileSBT(fCell1) ; ffs1_1.pushBackTimeStep(f1ts1)
        f1ts2.setFieldNoProfileSBT(fNode1) ; ffs1_2.pushBackTimeStep(f1ts2)
        #
        f1ts1=MEDFileField1TS()
        f1ts2=MEDFileField1TS()
        a=DataArrayDouble(m1.getNumberOfCells()) ; a.iota(1000.+off) ; a.setInfoOnComponents(["xx [m]"])
        fCell2.setArray(a)
        fCell2.setTime(*t)
        fCell2.checkCoherency()
        a=DataArrayDouble(m1.getNumberOfNodes()) ; a.iota(1000+off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode2.setArray(a)
        fNode2.setTime(*t)
        fNode2.checkCoherency()
        f1ts1.setFieldNoProfileSBT(fCell2) ; ffs2_1.pushBackTimeStep(f1ts1)
        f1ts2.setFieldNoProfileSBT(fNode2) ; ffs2_2.pushBackTimeStep(f1ts2)
        # TimeStep 1
        t=(2.1,1,0) ; off=100.
        f1ts1=MEDFileField1TS()
        f1ts2=MEDFileField1TS()
        a=DataArrayDouble(m0.getNumberOfCells()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        fCell1.setArray(a)
        fCell1.setTime(*t)
        fCell1.checkCoherency()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode1.setArray(a)
        fNode1.setTime(*t)
        fNode1.checkCoherency()
        f1ts1.setFieldNoProfileSBT(fCell1) ; ffs1_1.pushBackTimeStep(f1ts1)
        f1ts2.setFieldNoProfileSBT(fNode1) ; ffs1_2.pushBackTimeStep(f1ts2)
        #
        f1ts1=MEDFileField1TS()
        f1ts2=MEDFileField1TS()
        a=DataArrayDouble(m1.getNumberOfCells()) ; a.iota(1000.+off) ; a.setInfoOnComponents(["xx [m]"])
        fCell2.setArray(a)
        fCell2.setTime(*t)
        fCell2.checkCoherency()
        a=DataArrayDouble(m1.getNumberOfNodes()) ; a.iota(1000+off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode2.setArray(a)
        fNode2.setTime(*t)
        fNode2.checkCoherency()
        f1ts1.setFieldNoProfileSBT(fCell2) ; ffs2_1.pushBackTimeStep(f1ts1)
        f1ts2.setFieldNoProfileSBT(fNode2) ; ffs2_2.pushBackTimeStep(f1ts2)
        # TimeStep 2
        t=(3.1,2,0) ; off=200.
        f1ts1=MEDFileField1TS()
        f1ts2=MEDFileField1TS()
        a=DataArrayDouble(m0.getNumberOfCells()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        fCell1.setArray(a)
        fCell1.setTime(*t)
        fCell1.checkCoherency()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode1.setArray(a)
        fNode1.setTime(*t)
        fNode1.checkCoherency()
        f1ts1.setFieldNoProfileSBT(fCell1) ; ffs1_1.pushBackTimeStep(f1ts1)
        f1ts2.setFieldNoProfileSBT(fNode1) ; ffs1_2.pushBackTimeStep(f1ts2)
        #
        f1ts1=MEDFileField1TS()
        f1ts2=MEDFileField1TS()
        a=DataArrayDouble(m1.getNumberOfCells()) ; a.iota(1000.+off) ; a.setInfoOnComponents(["xx [m]"])
        fCell2.setArray(a)
        fCell2.setTime(*t)
        fCell2.checkCoherency()
        a=DataArrayDouble(m1.getNumberOfNodes()) ; a.iota(1000+off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode2.setArray(a)
        fNode2.setTime(*t)
        fNode2.checkCoherency()
        f1ts1.setFieldNoProfileSBT(fCell2) ; ffs2_1.pushBackTimeStep(f1ts1)
        f1ts2.setFieldNoProfileSBT(fNode2) ; ffs2_2.pushBackTimeStep(f1ts2)
        #
        mms.write(fname,2) ; mts.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),2)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),2)
        self.assertEqual(len(allFMTSLeavesToDisplay[1]),2)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 2 fields are defined on the same time steps
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),4)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),2)
        allFMTSLeavesPerCommonSupport2=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport2),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport2[0][0]),2)
        pass

    def test17(self):
        """ First test on GAUSS_NE (Elno). Here no Profiles.
        2 times steps.
        """
        fname="ForMEDReader17.med"
        fieldName1="MyFirstElno"
        fieldName2="ACellField"
        fieldName3="ANodeField"
        coo=DataArrayDouble([0.,0.,1.,0.,2.,0.,0.,1.,1.,1.,2.,1.],6,2)
        m=MEDCouplingUMesh("mesh",2)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_QUAD4,[0,3,4,1])
        m.insertNextCell(NORM_QUAD4,[1,4,5,2])
        m.checkCoherency2()
        #
        t=(1.1,0,-1)
        f=MEDCouplingFieldDouble(ON_GAUSS_NE) ; f.setTime(*t) ; f.setMesh(m)
        f.setArray(DataArrayDouble([3.,5.,7.,6.,2.,3.,11.,8.]))
        f.setName(fieldName1)
        f.checkCoherency()
        MEDLoader.WriteField(fname,f,True)
        f2=MEDCouplingFieldDouble(ON_CELLS) ; f2.setTime(*t) ; f2.setMesh(m)
        f2.setArray(DataArrayDouble([7.,11.],2,1))
        f2.setName(fieldName2)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f2)
        f3=MEDCouplingFieldDouble(ON_NODES) ; f3.setTime(*t) ; f3.setMesh(m)
        f3.setArray(DataArrayDouble([1.,2.,4.,1.,2.,4.],6,1))
        f3.setName(fieldName3)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f3)
        #
        t=(2.1,1,-1)
        f=MEDCouplingFieldDouble(ON_GAUSS_NE) ; f.setTime(*t) ; f.setMesh(m)
        f.setArray(DataArrayDouble([7.,6.,3.,5.,11.,8.,2.,3.]))
        f.setName(fieldName1)
        f.checkCoherency()
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f)
        f2=MEDCouplingFieldDouble(ON_CELLS) ; f2.setTime(*t) ; f2.setMesh(m)
        f2.setArray(DataArrayDouble([11.,7.],2,1))
        f2.setName(fieldName2)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f2)
        f3=MEDCouplingFieldDouble(ON_NODES) ; f3.setTime(*t) ; f3.setMesh(m)
        f3.setArray(DataArrayDouble([4.,2.,1.,4.,2.,1.],6,1))
        f3.setName(fieldName3)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f3)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),3)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),3)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),3)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc) # spaceDim 2 -> VTK wants 3D
        self.assertTrue(a0.isEqual(DataArrayDouble([0.,0.,0.,1.,0.,0.,2.,0.,0.,0.,1.,0.,1.,1.,0.,2.,1.,0.],6,3),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,0,3,4,1,4,1,4,5,2])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,0])))
        self.assertTrue(a7) # no copy here
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([0,1])))
        self.assertTrue(a9) # no copy here
        for i in xrange(1,2):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        vExp0=[DataArrayDouble([7.,11.]),DataArrayDouble([11.,7.])]
        vExp1=[DataArrayDouble([3.,5.,7.,6.,2.,3.,11.,8.]),DataArrayDouble([7.,6.,3.,5.,11.,8.,2.,3.])]
        for i in xrange(2):
            f=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName2)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            self.assertTrue(v.isEqual(vExp0[i],1e-12))
            #
            f=allFMTSLeavesPerCommonSupport1[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName1)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            self.assertTrue(v.isEqual(vExp1[i],1e-12))
            pass
        pass
    
    def test18(self):
        """ First test on GAUSS_PT. Here no Profiles. 2 times steps.
        """
        fname="ForMEDReader18.med"
        fieldName1="MyFirstGauss"
        fieldName2="ACellField"
        fieldName3="ANodeField"
        coo=DataArrayDouble([0.,0.,1.,0.,2.,0.,0.,1.,1.,1.,2.,1.],6,2)
        m=MEDCouplingUMesh("mesh",2)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_QUAD4,[0,3,4,1])
        m.insertNextCell(NORM_QUAD4,[1,4,5,2])
        m.checkCoherency2()
        #
        t=(1.1,0,-1)
        f=MEDCouplingFieldDouble(ON_GAUSS_PT) ; f.setTime(*t) ; f.setMesh(m)
        f.setGaussLocalizationOnType(NORM_QUAD4,[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.2,0.2,0.8,0.8],[0.7,0.3])
        f.setArray(DataArrayDouble([3.,5.,4.,6.])) ; f.getArray().setInfoOnComponents(["Smth"])
        f.setName(fieldName1)
        f.checkCoherency()
        MEDLoader.WriteField(fname,f,True)
        f2=MEDCouplingFieldDouble(ON_CELLS) ; f2.setTime(*t) ; f2.setMesh(m)
        f2.setArray(DataArrayDouble([7.,11.],2,1))
        f2.setName(fieldName2)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f2)
        f3=MEDCouplingFieldDouble(ON_NODES) ; f3.setTime(*t) ; f3.setMesh(m)
        f3.setArray(DataArrayDouble([1.,2.,4.,1.,2.,4.],6,1))
        f3.setName(fieldName3)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f3)
        #
        t=(2.1,1,-1)
        f=MEDCouplingFieldDouble(ON_GAUSS_PT) ; f.setTime(*t) ; f.setMesh(m)
        f.setGaussLocalizationOnType(NORM_QUAD4,[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.2,0.2,0.8,0.8],[0.7,0.3])
        f.setArray(DataArrayDouble([5.,3.,6.,4.])) ; f.getArray().setInfoOnComponents(["Smth"])
        f.setName(fieldName1)
        f.checkCoherency()
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f)
        f2=MEDCouplingFieldDouble(ON_CELLS) ; f2.setTime(*t) ; f2.setMesh(m)
        f2.setArray(DataArrayDouble([11.,7.],2,1))
        f2.setName(fieldName2)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f2)
        f3=MEDCouplingFieldDouble(ON_NODES) ; f3.setTime(*t) ; f3.setMesh(m)
        f3.setArray(DataArrayDouble([4.,2.,1.,4.,2.,1.],6,1))
        f3.setName(fieldName3)
        MEDLoader.WriteFieldUsingAlreadyWrittenMesh(fname,f3)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),3)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),3)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),3)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc) # spaceDim 2 -> VTK wants 3D
        self.assertTrue(a0.isEqual(DataArrayDouble([0.,0.,0.,1.,0.,0.,2.,0.,0.,0.,1.,0.,1.,1.,0.,2.,1.,0.],6,3),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,0,3,4,1,4,1,4,5,2])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,0])))
        self.assertTrue(a7) # no copy here
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([0,1])))
        self.assertTrue(a9) # no copy here
        for i in xrange(1,2):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        vExp0=[DataArrayDouble([7.,11.]),DataArrayDouble([11.,7.])]
        vExp1=[DataArrayDouble([3.,5.,4.,6.]),DataArrayDouble([5.,3.,6.,4.])]
        for i in xrange(2):
            f=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName2)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            self.assertTrue(v.isEqual(vExp0[i],1e-12))
            #
            f=allFMTSLeavesPerCommonSupport1[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName1)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp1[i].setInfoOnComponents(["Smth"])
            self.assertTrue(v.isEqual(vExp1[i],1e-12))
            pass
        ## Now same exercise but with a different load strategy. All is load directly.
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname) # here all is read, the SauvReader (or other Reader) is emulated
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),3)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),3)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),3)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc) # spaceDim 2 -> VTK wants 3D
        self.assertTrue(a0.isEqual(DataArrayDouble([0.,0.,0.,1.,0.,0.,2.,0.,0.,0.,1.,0.,1.,1.,0.,2.,1.,0.],6,3),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,0,3,4,1,4,1,4,5,2])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,0])))
        self.assertTrue(a7) # no copy here
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([0,1])))
        self.assertTrue(a9) # no copy here
        for i in xrange(1,2):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        vExp0=[DataArrayDouble([7.,11.]),DataArrayDouble([11.,7.])]
        vExp1=[DataArrayDouble([3.,5.,4.,6.]),DataArrayDouble([5.,3.,6.,4.])]
        for i in xrange(2):
            f=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst) # no load needed here
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName2)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            self.assertTrue(v.isEqual(vExp0[i],1e-12))
            #
            f=allFMTSLeavesPerCommonSupport1[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst) # no load needed here
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName1)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp1[i].setInfoOnComponents(["Smth"])
            self.assertTrue(v.isEqual(vExp1[i],1e-12))
            pass
        pass
    
    def test19(self):
        """
        This test is a simple non profile CELL field but lying on cells of dimension -1 (not 0 as "usual").
        """
        fname="ForMEDReader19.med"
        fieldName="ACellFieldOnDimM1"
        coo=DataArrayDouble(3) ; coo.iota()
        m=MEDCouplingCMesh() ; m.setCoords(coo,coo,coo) ; m.setName("mesh")
        m0=m.buildUnstructured() ; del m
        m1=m0.computeSkin()
        #
        mm=MEDFileUMesh()                                
        mm.setMeshAtLevel(0,m0)
        mm.setMeshAtLevel(-1,m1)
        ff=MEDFileFieldMultiTS()
        # time 0
        t=(1.1,1,-1)
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setTime(*t) ; f.setMesh(m1)
        f.setName(fieldName)
        arr=DataArrayDouble(24) ; arr.iota() ; arr.setInfoOnComponents(["AStr"])
        f.setArray(arr)
        f.checkCoherency()
        f1ts=MEDFileField1TS() ; f1ts.setFieldNoProfileSBT(f)
        ff.pushBackTimeStep(f1ts)
        # time 1
        t=(2.1,2,-2)
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setTime(*t) ; f.setMesh(m1)
        f.setName(fieldName)
        arr=DataArrayDouble(24) ; arr.iota() ; arr.reverse() ; arr.setInfoOnComponents(["AStr"])
        f.setArray(arr)
        f.checkCoherency()
        f1ts=MEDFileField1TS() ; f1ts.setFieldNoProfileSBT(f)
        ff.pushBackTimeStep(f1ts)
        #
        mm.write(fname,2)
        ff.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                allFMTSLeavesToDisplay2+=fmts.splitDiscretizations()
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),1)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),1)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(ncc)
        self.assertTrue(a0.isEqual(DataArrayDouble([0.,0.,0.,1.,0.,0.,2.,0.,0.,0.,1.,0.,1.,1.,0.,2.,1.,0.,0.,2.,0.,1.,2.,0.,2.,2.,0.,0.,0.,1.,1.,0.,1.,2.,0.,1.,0.,1.,1.,1.,1.,1.,2.,1.,1.,0.,2.,1.,1.,2.,1.,2.,2.,1.,0.,0.,2.,1.,0.,2.,2.,0.,2.,0.,1.,2.,1.,1.,2.,2.,1.,2.,0.,2.,2.,1.,2.,2.,2.,2.,2.],27,3),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,1,0,3,4,4,1,10,9,0,4,0,9,12,3,4,2,1,4,5,4,2,11,10,1,4,5,14,11,2,4,4,3,6,7,4,3,12,15,6,4,6,15,16,7,4,5,4,7,8,4,7,16,17,8,4,8,17,14,5,4,19,22,21,18,4,10,19,18,9,4,9,18,21,12,4,20,23,22,19,4,11,20,19,10,4,14,23,20,11,4,22,25,24,21,4,12,21,24,15,4,15,24,25,16,4,23,26,25,22,4,16,25,26,17,4,17,26,23,14])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])))
        self.assertTrue(a7) # no copy here
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8 is None)
        self.assertTrue(a9) # no copy here
        for i in xrange(1,2):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        vExp0=[DataArrayDouble([7.,11.]),DataArrayDouble([11.,7.])]
        for i in xrange(2):
            f=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(24) ; vExp.iota()
            if i==1: vExp.reverse()
            vExp.setInfoOnComponents(["AStr"])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass

    pass

unittest.main()
