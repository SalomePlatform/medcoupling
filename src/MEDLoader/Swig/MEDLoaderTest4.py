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
        fcscp=MEDFileFastCellSupportComparator.New(mst,allFMTSLeavesPerCommonSupport[0][0])
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        self.assertTrue(isinstance(mml,MEDUMeshMultiLev))
        for i in xrange(1,5):
            self.assertTrue(fcscp.isDataSetSupportEqualToThePreviousOne(i,fields))
            pass
        a0,a1,a2,a3,a4,a5=mml.buildVTUArrays()
        self.assertTrue(a0.isEqual(a0Exp,1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,5,5,5,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,67,70,73,76,80,84,88,92,96])))
        self.assertTrue(a3.isEqual(DataArrayInt([2,0,1,2,1,2,2,2,0,2,3,4,2,4,5,2,5,3,2,6,7,2,7,8,2,8,6,2,9,10,2,10,11,2,11,9,2,12,13,2,13,14,2,14,15,2,15,12,2,16,17,2,17,18,2,18,19,2,19,16,2,20,21,2,21,22,2,22,23,2,23,20,2,24,25,2,25,26,2,26,27,2,27,24,2,28,29,2,29,30,2,30,31,2,31,28,3,0,1,2,3,3,4,5,3,6,7,8,3,9,10,11,4,12,13,14,15,4,16,17,18,19,4,20,21,22,23,4,24,25,26,27,4,28,29,30,31])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
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
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0]),2)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=MEDFileFastCellSupportComparator.New(mst,allFMTSLeavesPerCommonSupport[0][0])
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        assert isinstance(mml,MEDUMeshMultiLev)
        for i in xrange(1,5):
            self.assertTrue(fcscp.isDataSetSupportEqualToThePreviousOne(i,fields))
            pass
        a0,a1,a2,a3,a4,a5=mml.buildVTUArrays()
        self.assertTrue(a0.isEqual(a0Exp,1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5,5,5,5,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,67,70,73,76,80,84,88,92,96])))
        self.assertTrue(a3.isEqual(DataArrayInt([2,0,1,2,1,2,2,2,0,2,3,4,2,4,5,2,5,3,2,6,7,2,7,8,2,8,6,2,9,10,2,10,11,2,11,9,2,12,13,2,13,14,2,14,15,2,15,12,2,16,17,2,17,18,2,18,19,2,19,16,2,20,21,2,21,22,2,22,23,2,23,20,2,24,25,2,25,26,2,26,27,2,27,24,2,28,29,2,29,30,2,30,31,2,31,28,3,0,1,2,3,3,4,5,3,6,7,8,3,9,10,11,4,12,13,14,15,4,16,17,18,19,4,20,21,22,23,4,24,25,26,27,4,28,29,30,31])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        # for cells
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)# Second 0 is for cells
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([200.,201.,202.,203.,204.,205.,206.,207.,208.,209.,210.,211.,212.,213.,214.,215.,216.,217.,218.,219.,220.,221.,222.,223.,224.,225.,226.,227.,228.,229.,230.,231.,232.,233.,234.,235.,236.,237.,238.,239.,240.,241.,242.,243.,244.,245.,246.,247.,248.,249.,250.,251.,252.,253.,254.,255.,256.,257.,258.,259.,260.,261.,262.,263.,100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,110.,111.,112.,113.,114.,115.,116.,117.],41,2) ; vExp.setInfoOnComponents(['Comp1 [m]','Com2 [s^2]']) ; vExp+=i*1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        for i in xrange(5):
            f=allFMTSLeavesPerCommonSupport[0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)# Second 0 is for cells
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([300.,301.,302.,303.,304.,305.,306.,307.,308.,309.,310.,311.,312.,313.,314.,315.,316.,317.,318.,319.,320.,321.,322.,323.,324.,325.,326.,327.,328.,329.,330.,331.,332.,333.,334.,335.,336.,337.,338.,339.,340.,341.,342.,343.,344.,345.,346.,347.,348.,349.,350.,351.,352.,353.,354.,355.,356.,357.,358.,359.,360.,361.,362.,363.],32,2) ; vExp.setInfoOnComponents(['Comp1 [m]','Com2 [s^2]']) ; vExp.setInfoOnComponents(['Comp1 [m]','Com2 [s^2]']) ; vExp+=i*1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass
    pass

unittest.main()
