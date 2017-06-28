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
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell1)
        fs.pushBackTimeStep(f)
        ##### Time step 1
        i=1
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(1100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(1200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell1)
        fs.pushBackTimeStep(f)
        ##### Time step 2
        i=2
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(2100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(2200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
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
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(3200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
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
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(4200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell1)
        fs.pushBackTimeStep(f)
        mm.write(fname,2)
        fs.write(fname,0)
        a0Exp=mm.getCoords().deepCopy()
        del m,m1,mm,fs,f,fCell0,fCell1
        ########## GO for reading in MEDReader, by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False) # False is important to not read the values
        fields.removeFieldsWithoutAnyTimeStep()
        refMem=fields.getHeapMemorySize()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        for i in range(1, 5):
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(5):
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
        fieldName="zeField"
        fs=MEDFileFieldMultiTS()
        ##### Time step 0
        i=0
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell1)
        #
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs.pushBackTimeStep(f)
        ##### Time step 1
        i=1
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(1100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(1200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell1)
        #
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(1300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs.pushBackTimeStep(f)
        ##### Time step 2
        i=2
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(2300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        #
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(2100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(2200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
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
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(3200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell1)
        #
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(3300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkConsistencyLight()
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
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(4200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell1)
        #
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(4300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        #
        fs.pushBackTimeStep(f)
        mm.write(fname,2)
        fs.write(fname,0)
        a0Exp=mm.getCoords().deepCopy()
        del m,m1,mm,fs,f,fCell0,fCell1
        ########## GO for reading in MEDReader, by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        for i in range(1, 5):
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        # for cells
        for i in range(5):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)# Second 0 is for cells
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([200.,201.,202.,203.,204.,205.,206.,207.,208.,209.,210.,211.,212.,213.,214.,215.,216.,217.,218.,219.,220.,221.,222.,223.,224.,225.,226.,227.,228.,229.,230.,231.,232.,233.,234.,235.,236.,237.,238.,239.,240.,241.,242.,243.,244.,245.,246.,247.,248.,249.,250.,251.,252.,253.,254.,255.,256.,257.,258.,259.,260.,261.,262.,263.,100.,101.,102.,103.,104.,105.,106.,107.,108.,109.,110.,111.,112.,113.,114.,115.,116.,117.],41,2) ; vExp.setInfoOnComponents(['Comp1 [m]','Com2 [s^2]']) ; vExp+=i*1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        for i in range(5):
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
        fieldName="zeField"
        fs=MEDFileFieldMultiTS()
        ##### Time step 0 on cells
        i=0
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell1)
        fs.pushBackTimeStep(f)
        ##### Time step 1 on nodes
        i=1
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(1300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs.pushBackTimeStep(f)
        ##### Time step 2 on cells
        i=2
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(2100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(2200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell1)
        fs.pushBackTimeStep(f)
        ##### Time step 3 on nodes
        i=3
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName) ; fNode.setMesh(m1)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(3300) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs.pushBackTimeStep(f)
        ##### Time step 4
        i=4
        f=MEDFileField1TS()
        fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i),i,0)
        fCell0.setName(fieldName) ; fCell0.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfCells()) ; arr.iota(4100) ; arr.rearrange(2)
        fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell0.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell0)
        #
        fCell1=MEDCouplingFieldDouble(ON_CELLS) ; fCell1.setTime(float(i),i,0)
        fCell1.setName(fieldName) ; fCell1.setMesh(m1)
        arr=DataArrayDouble(2*m1.getNumberOfCells()) ; arr.iota(4200) ; arr.rearrange(2)
        fCell1.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]","Com2 [s^2]"])
        fCell1.checkConsistencyLight()
        f.setFieldNoProfileSBT(fCell1)
        #
        fs.pushBackTimeStep(f)
        mm.write(fname,2)
        fs.write(fname,0)
        a0Exp=mm.getCoords().deepCopy()
        del m,m1,mm,fs,f,fCell0,fCell1
        ########## GO for reading in MEDReader, by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        for i in range(1, 3):
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        # for cells
        for i in range(3):
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
        for i in range(1, 2):
            self.assertTrue(fcscp.isDataSetSupportEqualToThePreviousOne(i,fields))
            pass
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0Exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,5,5,5,9,9,9,9,9,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,8,12,16,21,26,31,36,41,44,47,50,53,56,59,62,65,68,71,74,77,80,83,86,89,92,95,98,101,104,107,110,113,116,119,122,125,128,131,134])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,0,1,2,3,3,4,5,3,6,7,8,3,9,10,11,4,12,13,14,15,4,16,17,18,19,4,20,21,22,23,4,24,25,26,27,4,28,29,30,31,2,0,1,2,1,2,2,2,0,2,3,4,2,4,5,2,5,3,2,6,7,2,7,8,2,8,6,2,9,10,2,10,11,2,11,9,2,12,13,2,13,14,2,14,15,2,15,12,2,16,17,2,17,18,2,18,19,2,19,16,2,20,21,2,21,22,2,22,23,2,23,20,2,24,25,2,25,26,2,26,27,2,27,24,2,28,29,2,29,30,2,30,31,2,31,28])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(2):
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
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs1.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName2) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(100+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs2.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName3) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(200+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs3.pushBackTimeStep(f)
        ##### Time step 1
        i=1
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName1) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(0+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs1.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName2) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(100+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs2.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName3) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(200+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs3.pushBackTimeStep(f)
        ##### Time step 2
        i=2
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName1) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(0+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs1.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName2) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(100+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs2.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName3) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(200+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs3.pushBackTimeStep(f)
        ##### Time step 3
        i=3
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName1) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(0+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs1.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName2) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(100+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs2.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName3) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(200+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs3.pushBackTimeStep(f)
        ##### Time step 4
        i=4
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName1) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(0+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs1.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName2) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(100+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs2.pushBackTimeStep(f)
        #
        f=MEDFileField1TS()
        fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
        fNode.setName(fieldName3) ; fNode.setMesh(m)
        arr=DataArrayDouble(2*m.getNumberOfNodes()) ; arr.iota(200+1000*i) ; arr.rearrange(2)
        fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"])
        fNode.checkConsistencyLight()
        f.setFieldNoProfileSBT(fNode)
        fs3.pushBackTimeStep(f)
        #
        mm.write(fname,2)
        fs1.write(fname,0) ; fs2.write(fname,0) ; fs3.write(fname,0)
        a0Exp=mm.getCoords().deepCopy()
        del m,mm,fs1,fs2,fs3,f,fNode
        ########## GO for reading in MEDReader, by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        for i in range(1, 5):
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        # test all the time steps of the 1/1 time step serie, on field 1
        for i in range(5):
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
        for i in range(5):
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
        for i in range(5):
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
        for i in range(5):
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
        a0Exp=mm.getCoords().deepCopy()
        del m,mm,fs1,fs2,fs3,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        for i in range(1, 5):
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(5):
            nbOfT=[6,8]
            fieldNames=[fieldName1,fieldName2]
            for j in range(2):
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
        for i in range(1, 5):
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(5):
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
        num=DataArrayInt(15) ; num.iota(200) ; mm.setRenumFieldArr(1,num) ; del num
        #
        fieldName0="zeField0" ; # on cells
        fieldName1="zeField1" ; pfl1=DataArrayInt([2,3,6,7]) ; pfl1.setName("pfl1") # on cells
        fieldName2="zeField2" ; pfl2=DataArrayInt([2,3,4,7,8,9,12,13,14]) ; pfl2.setName("pfl2") # on nodes
        fieldName3="zeField3" ; pfl3=DataArrayInt([2,3,5,7]) ; pfl3.setName("pfl3") # on cells but different support
        fieldName4="zeField4" ;# on nodes
        fs0=MEDFileFieldMultiTS() ; fs1=MEDFileFieldMultiTS() ; fs2=MEDFileFieldMultiTS() ; fs3=MEDFileFieldMultiTS() ; fs4=MEDFileFieldMultiTS()
        #
        for i in range(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*8) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkConsistencyLight()
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
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_4 [m]","Com2_4 [s^2]"]) ; fNode.checkConsistencyLight()
            f.setFieldNoProfileSBT(fNode)
            fs4.pushBackTimeStep(f)
            pass
        mm.write(fname,2)
        fs0.write(fname,0) ; fs1.write(fname,0) ; fs2.write(fname,0) ; fs3.write(fname,0) ; fs4.write(fname,0)
        del m,mm,fs1,fs2,fs3,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        (a,b),c=mml2.buildVTUArrays()
        self.assertTrue(c)# c is True here because the returned array is directly those coming from internal structure
        self.assertTrue(a.isEqual(coordsX,1e-12))
        self.assertTrue(b.isEqual(coordsY,1e-12))
        self.assertTrue(isinstance(mml2,MEDCMeshMultiLev))
        for i in range(1, 5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,1,2,3,4,5,6,7])))
        self.assertTrue(a7) # True because no copy
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([100,101,102,103,104,105,106,107])))
        self.assertTrue(a9) # True because no copy
        a10,a11=mml2.retrieveNumberIdsOnNodes()
        self.assertTrue(a10.isEqual(DataArrayInt([200,201,202,203,204,205,206,207,208,209,210,211,212,213,214])))
        self.assertTrue(a11) # True because no copy
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(5):
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
        (a,b),c=mml2.buildVTUArrays()
        self.assertTrue(not c)# c is False because this a sub support specialy built for buildVTUArrays
        self.assertTrue(a.isEqual(coordsX[[2,3,4]],1e-12))
        self.assertTrue(b.isEqual(coordsY,1e-12))
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([2,3,6,7])))
        self.assertTrue(not a7) # False because copy
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([102,103,106,107])))
        self.assertTrue(not a9) # False because copy
        a10,a11=mml2.retrieveNumberIdsOnNodes()
        self.assertTrue(a10.isEqual(DataArrayInt([202,203,204,207,208,209,212,213,214])))
        self.assertTrue(not a11) # False because copy
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in range(5):
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(5):
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
        for i in range(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*8) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkConsistencyLight()
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
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_4 [m]","Com2_4 [s^2]"]) ; fNode.checkConsistencyLight()
            f.setFieldNoProfileSBT(fNode)
            fs4.pushBackTimeStep(f)
            pass
        mm.write(fname,2)
        fs0.write(fname,0) ; fs1.write(fname,0) ; fs2.write(fname,0) ; fs3.write(fname,0) ; fs4.write(fname,0)
        del m,mm,fs1,fs2,fs3,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        a,b,c=mml2.buildVTUArrays()
        self.assertTrue(c)#True here because a is directly coming from internal data without copy
        self.assertTrue(a.isEqual(a0Exp,1e-12))
        self.assertEqual(b,[5,3])
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,1,2,3,4,5,6,7])))
        self.assertTrue(a7) # True because no copy
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([100,101,102,103,104,105,106,107])))
        self.assertTrue(a9) # True because no copy
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in range(5):
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
        a,b,c=mml2.buildVTUArrays()
        self.assertTrue(not c)#False here because a is the result of a computation not the internal strucutre
        self.assertTrue(a.isEqual(a0Exp[pfl2],1e-12))
        self.assertEqual(b,[3,3])
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([2,3,6,7])))
        self.assertTrue(not a7) # False because copy
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([102,103,106,107])))
        self.assertTrue(not a9) # False because copy
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in range(5):
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(5):
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
        for i in range(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_GAUSS_NE) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*38) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkConsistencyLight()
            f.setFieldNoProfileSBT(fNode)
            fs0.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName1) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*11) ; arr.iota(100+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"]) ; fNode.checkConsistencyLight()
            f.setFieldNoProfileSBT(fNode)
            fs1.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName2) ; fNode.setMesh(m)
            fNode.setGaussLocalizationOnCells([0,1,2,3],[0.,0.,1.,0.,0.,1.],[0.5,0.5,0.7,0.7],[0.8,0.2])
            fNode.setGaussLocalizationOnCells([4,5],[0.,0.,1.,0.,0.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3],[0.8,0.05,0.1,0.04,0.01])
            fNode.setGaussLocalizationOnCells([6,7,8],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2],[0.8,0.05,0.1,0.04])
            fNode.setGaussLocalizationOnCells([9,10],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.8,0.8],[0.8,0.05,0.1,0.01,0.02,0.005,0.005])
            arr=DataArrayDouble(2*(4*2+2*5+3*4+2*7)) ; arr.iota(300+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"]) ; fNode.checkConsistencyLight()
            f.setFieldNoProfileSBT(fNode)
            fs2.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_NODES) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName3) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*15) ; arr.iota(400+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_3 [m]","Com2_3 [s^2]"]) ; fNode.checkConsistencyLight()
            f.setFieldNoProfileSBT(fNode)
            fs3.pushBackTimeStep(f)
            #
            pass
        #
        mm.write(fname,2)
        fs0.write(fname,0) ; fs1.write(fname,0) ; fs2.write(fname,0) ; fs3.write(fname,0)
        a0Exp=mm.getCoords().deepCopy()
        del m,mm,fs1,fs2,fs3,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                #for itmp in tmp:
                #    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                #    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in range(5):
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
        for i in range(5):
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
            fNode.setGaussLocalizationOnCells([2,3],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2],[0.8,0.05,0.1,0.04])
            fNode.setGaussLocalizationOnCells([4],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.8,0.8],[0.8,0.05,0.1,0.01,0.02,0.005,0.005])
            arr=DataArrayDouble(2*(2*1+5*1+4*2+7*1)) ; arr.iota(300+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"]) ; fNode.checkConsistencyLight()
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
        a0Exp=mm.getCoords().deepCopy()
        del m,mm,fs1,fs2,fs3,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                #for itmp in tmp:
                #    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                #    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in range(5):
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
        for i in range(5):
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
        a0Exp=mm.getCoords().deepCopy()
        del m,mm,fs1,fs2,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in range(5):
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
        for i in range(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0) ; fNode.setMesh(m)
            fNode.setGaussLocalizationOnCells([0,2,3,4,7,15],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7],[0.8,0.2])
            fNode.setGaussLocalizationOnCells([1,5,8,9],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3],[0.8,0.05,0.1,0.04,0.01])
            fNode.setGaussLocalizationOnCells([6,10,13],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2],[0.8,0.05,0.1,0.04])
            fNode.setGaussLocalizationOnCells([11,12,14],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.8,0.8],[0.8,0.05,0.1,0.01,0.02,0.005,0.005])
            arr=DataArrayDouble(2*(2*6+5*4+4*3+7*3)) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkConsistencyLight()
            f.setFieldNoProfileSBT(fNode)
            fs0.pushBackTimeStep(f)
            pass
        mm.write(fname,2)
        fs0.write(fname,0)
        a0Exp=mm.getCoords().deepCopy()
        del m,mm,fs0,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                #for itmp in tmp:
                #    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                #    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in range(5):
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
        for i in range(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0) ; fNode.setMesh(m)
            fNode.setGaussLocalizationOnCells([0,2,3,4,7,15],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7],[0.8,0.2])
            fNode.setGaussLocalizationOnCells([1,5,8,9],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3],[0.8,0.05,0.1,0.04,0.01])
            fNode.setGaussLocalizationOnCells([6,10,13],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2],[0.8,0.05,0.1,0.04])
            fNode.setGaussLocalizationOnCells([11,12,14],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.8,0.8],[0.8,0.05,0.1,0.01,0.02,0.005,0.005])
            arr=DataArrayDouble(2*(2*6+5*4+4*3+7*3)) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkConsistencyLight()
            f.setFieldNoProfileSBT(fNode)
            fs0.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName1) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*16) ; arr.iota(300+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"]) ; fNode.checkConsistencyLight()
            f.setFieldNoProfileSBT(fNode)
            fs1.pushBackTimeStep(f)
            pass
        mm.write(fname,2)
        fs0.write(fname,0) ; fs1.write(fname,0)
        a0Exp=mm.getCoords().deepCopy()
        del m,mm,fs0,fs1,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                #for itmp in tmp:
                #    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                #    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in range(5):
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
            for i in range(5):
                f=MEDFileField1TS()
                fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
                fNode.setName(fieldName0) ; fNode.setMesh(m)
                arr=DataArrayDouble(2*5) ; arr.iota(0+1000*i) ; arr.rearrange(2)
                fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkConsistencyLight()
                f.setFieldNoProfileSBT(fNode)
                fs0.pushBackTimeStep(f)
                #
                f=MEDFileField1TS()
                fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
                fNode.setName(fieldName1) ; fNode.setMesh(m)
                arr=DataArrayDouble(2*5) ; arr.iota(100+1000*i) ; arr.rearrange(2)
                fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"]) ; fNode.checkConsistencyLight()
                f.setFieldNoProfileSBT(fNode)
                fs1.pushBackTimeStep(f)
                #
                f=MEDFileField1TS()
                fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
                fNode.setName(fieldName2) ; fNode.setMesh(m[pfl1])
                arr=DataArrayDouble(2*2) ; arr.iota(200+1000*i) ; arr.rearrange(2)
                fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_2 [m]","Com2_2 [s^2]"]) ; fNode.checkConsistencyLight()
                f.setFieldProfile(fNode,mm,0,pfl1)
                fs2.pushBackTimeStep(f)
                #
                f=MEDFileField1TS()
                fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
                fNode.setName(fieldName3) ; fNode.setMesh(m[pfl2])
                arr=DataArrayDouble(2*3) ; arr.iota(300+1000*i) ; arr.rearrange(2)
                fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_3 [m]","Com2_3 [s^2]"]) ; fNode.checkConsistencyLight()
                f.setFieldProfile(fNode,mm,0,pfl2)
                fs3.pushBackTimeStep(f)
                pass
            mm.write(fname,2)
            fs0.write(fname,0) ; fs1.write(fname,0) ; fs2.write(fname,0) ; fs3.write(fname,0)
            a0Exp=mm.getCoords().deepCopy()
            del m,mm,fs0
            ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
            ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
            fields=MEDFileFields(fname,False)
            fields.removeFieldsWithoutAnyTimeStep()
            fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
            allFMTSLeavesToDisplay=[]
            for fields in fields_per_mesh:
                allFMTSLeavesToDisplay2=[]
                for fmts in fields:
                    tmp=fmts.splitDiscretizations()
                    for itmp in tmp:
                        self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                        pass
                    allFMTSLeavesToDisplay2+=tmp
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
            self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
            for i in range(1, 5):
                self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
                pass
            pass
            for i in range(5):
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
            self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
            for i in range(1, 5):
                self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
                pass
            pass
            for i in range(5):
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
            for i in range(1, 5):
                self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
                pass
            pass
            for i in range(5):
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
            for i in range(5):
                f=MEDFileField1TS()
                fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
                fNode.setName(fieldName0) ; fNode.setMesh(m)
                arr=DataArrayDouble(2*3) ; arr.iota(0+1000*i) ; arr.rearrange(2)
                fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkConsistencyLight()
                f.setFieldNoProfileSBT(fNode)
                fs0.pushBackTimeStep(f)
                #
                f=MEDFileField1TS()
                fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
                fNode.setName(fieldName1) ; fNode.setMesh(m)
                arr=DataArrayDouble(2*3) ; arr.iota(100+1000*i) ; arr.rearrange(2)
                fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"]) ; fNode.checkConsistencyLight()
                f.setFieldNoProfileSBT(fNode)
                fs1.pushBackTimeStep(f)
                pass
            mm.write(fname,2)
            fs0.write(fname,0) ; fs1.write(fname,0)
            a0Exp=mm.getCoords().deepCopy()
            del m,mm,fs0
            ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
            ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
            fields=MEDFileFields(fname,False)
            fields.removeFieldsWithoutAnyTimeStep()
            fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
            allFMTSLeavesToDisplay=[]
            for fields in fields_per_mesh:
                allFMTSLeavesToDisplay2=[]
                for fmts in fields:
                    tmp=fmts.splitDiscretizations()
                    for itmp in tmp:
                        self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                        pass
                    allFMTSLeavesToDisplay2+=tmp
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
            self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
            for i in range(1, 5):
                self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
                pass
            a6,a7=mml2.retrieveFamilyIdsOnCells()
            self.assertTrue(a6.isEqual(DataArrayInt([0,0,0])))
            self.assertTrue(a7)
            a8,a9=mml2.retrieveNumberIdsOnCells()
            self.assertTrue(a8 is None)
            self.assertTrue(a9)
            for i in range(5):
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
        fCell.checkConsistencyLight()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode.setArray(a)
        fNode.setTime(*t)
        fNode.checkConsistencyLight()
        f1ts.setFieldNoProfileSBT(fCell)
        f1ts.setFieldNoProfileSBT(fNode)
        ffs.pushBackTimeStep(f1ts)
        # TimeStep 1
        t=(2.1,1,0) ; off=100.
        f1ts=MEDFileField1TS()
        a=DataArrayDouble(m0.getNumberOfCells()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        fCell.setArray(a)
        fCell.setTime(*t)
        fCell.checkConsistencyLight()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode.setArray(a)
        fNode.setTime(*t)
        fNode.checkConsistencyLight()
        f1ts.setFieldNoProfileSBT(fCell)
        f1ts.setFieldNoProfileSBT(fNode)
        ffs.pushBackTimeStep(f1ts)
        # TimeStep 2
        t=(3.2,2,0) ; off=200.
        f1ts=MEDFileField1TS()
        a=DataArrayDouble(m0.getNumberOfCells()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        fCell.setArray(a)
        fCell.setTime(*t)
        fCell.checkConsistencyLight()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode.setArray(a)
        fNode.setTime(*t)
        fNode.checkConsistencyLight()
        f1ts.setFieldNoProfileSBT(fCell)
        f1ts.setFieldNoProfileSBT(fNode)
        ffs.pushBackTimeStep(f1ts)
        # TimeStep 3
        t=(4.3,3,1) ; off=300.
        f1ts=MEDFileField1TS()
        a=DataArrayDouble(m0.getNumberOfCells()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        fCell.setArray(a)
        fCell.setTime(*t)
        fCell.checkConsistencyLight()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode.setArray(a)
        fNode.setTime(*t)
        fNode.checkConsistencyLight()
        f1ts.setFieldNoProfileSBT(fCell)
        f1ts.setFieldNoProfileSBT(fNode)
        ffs.pushBackTimeStep(f1ts)
        #
        mm.write(fname,2)
        ffs.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
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
        m1=m0.deepCopy() ; m1.translate([2.5,0.,0.]) ; m1.setName("mesh2")
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
        fCell1.checkConsistencyLight()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode1.setArray(a)
        fNode1.setTime(*t)
        fNode1.checkConsistencyLight()
        f1ts1.setFieldNoProfileSBT(fCell1) ; ffs1_1.pushBackTimeStep(f1ts1)
        f1ts2.setFieldNoProfileSBT(fNode1) ; ffs1_2.pushBackTimeStep(f1ts2)
        #
        f1ts1=MEDFileField1TS()
        f1ts2=MEDFileField1TS()
        a=DataArrayDouble(m1.getNumberOfCells()) ; a.iota(1000.+off) ; a.setInfoOnComponents(["xx [m]"])
        fCell2.setArray(a)
        fCell2.setTime(*t)
        fCell2.checkConsistencyLight()
        a=DataArrayDouble(m1.getNumberOfNodes()) ; a.iota(1000+off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode2.setArray(a)
        fNode2.setTime(*t)
        fNode2.checkConsistencyLight()
        f1ts1.setFieldNoProfileSBT(fCell2) ; ffs2_1.pushBackTimeStep(f1ts1)
        f1ts2.setFieldNoProfileSBT(fNode2) ; ffs2_2.pushBackTimeStep(f1ts2)
        # TimeStep 1
        t=(2.1,1,0) ; off=100.
        f1ts1=MEDFileField1TS()
        f1ts2=MEDFileField1TS()
        a=DataArrayDouble(m0.getNumberOfCells()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        fCell1.setArray(a)
        fCell1.setTime(*t)
        fCell1.checkConsistencyLight()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode1.setArray(a)
        fNode1.setTime(*t)
        fNode1.checkConsistencyLight()
        f1ts1.setFieldNoProfileSBT(fCell1) ; ffs1_1.pushBackTimeStep(f1ts1)
        f1ts2.setFieldNoProfileSBT(fNode1) ; ffs1_2.pushBackTimeStep(f1ts2)
        #
        f1ts1=MEDFileField1TS()
        f1ts2=MEDFileField1TS()
        a=DataArrayDouble(m1.getNumberOfCells()) ; a.iota(1000.+off) ; a.setInfoOnComponents(["xx [m]"])
        fCell2.setArray(a)
        fCell2.setTime(*t)
        fCell2.checkConsistencyLight()
        a=DataArrayDouble(m1.getNumberOfNodes()) ; a.iota(1000+off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode2.setArray(a)
        fNode2.setTime(*t)
        fNode2.checkConsistencyLight()
        f1ts1.setFieldNoProfileSBT(fCell2) ; ffs2_1.pushBackTimeStep(f1ts1)
        f1ts2.setFieldNoProfileSBT(fNode2) ; ffs2_2.pushBackTimeStep(f1ts2)
        # TimeStep 2
        t=(3.1,2,0) ; off=200.
        f1ts1=MEDFileField1TS()
        f1ts2=MEDFileField1TS()
        a=DataArrayDouble(m0.getNumberOfCells()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        fCell1.setArray(a)
        fCell1.setTime(*t)
        fCell1.checkConsistencyLight()
        a=DataArrayDouble(m0.getNumberOfNodes()) ; a.iota(off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode1.setArray(a)
        fNode1.setTime(*t)
        fNode1.checkConsistencyLight()
        f1ts1.setFieldNoProfileSBT(fCell1) ; ffs1_1.pushBackTimeStep(f1ts1)
        f1ts2.setFieldNoProfileSBT(fNode1) ; ffs1_2.pushBackTimeStep(f1ts2)
        #
        f1ts1=MEDFileField1TS()
        f1ts2=MEDFileField1TS()
        a=DataArrayDouble(m1.getNumberOfCells()) ; a.iota(1000.+off) ; a.setInfoOnComponents(["xx [m]"])
        fCell2.setArray(a)
        fCell2.setTime(*t)
        fCell2.checkConsistencyLight()
        a=DataArrayDouble(m1.getNumberOfNodes()) ; a.iota(1000+off) ; a.setInfoOnComponents(["xx [m]"])
        a=a.negate()
        fNode2.setArray(a)
        fNode2.setTime(*t)
        fNode2.checkConsistencyLight()
        f1ts1.setFieldNoProfileSBT(fCell2) ; ffs2_1.pushBackTimeStep(f1ts1)
        f1ts2.setFieldNoProfileSBT(fNode2) ; ffs2_2.pushBackTimeStep(f1ts2)
        #
        mms.write(fname,2) ; mts.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        m.checkConsistency()
        #
        t=(1.1,0,-1)
        f=MEDCouplingFieldDouble(ON_GAUSS_NE) ; f.setTime(*t) ; f.setMesh(m)
        f.setArray(DataArrayDouble([3.,5.,7.,6.,2.,3.,11.,8.]))
        f.setName(fieldName1)
        f.checkConsistencyLight()
        WriteField(fname,f,True)
        f2=MEDCouplingFieldDouble(ON_CELLS) ; f2.setTime(*t) ; f2.setMesh(m)
        f2.setArray(DataArrayDouble([7.,11.],2,1))
        f2.setName(fieldName2)
        WriteFieldUsingAlreadyWrittenMesh(fname,f2)
        f3=MEDCouplingFieldDouble(ON_NODES) ; f3.setTime(*t) ; f3.setMesh(m)
        f3.setArray(DataArrayDouble([1.,2.,4.,1.,2.,4.],6,1))
        f3.setName(fieldName3)
        WriteFieldUsingAlreadyWrittenMesh(fname,f3)
        #
        t=(2.1,1,-1)
        f=MEDCouplingFieldDouble(ON_GAUSS_NE) ; f.setTime(*t) ; f.setMesh(m)
        f.setArray(DataArrayDouble([7.,6.,3.,5.,11.,8.,2.,3.]))
        f.setName(fieldName1)
        f.checkConsistencyLight()
        WriteFieldUsingAlreadyWrittenMesh(fname,f)
        f2=MEDCouplingFieldDouble(ON_CELLS) ; f2.setTime(*t) ; f2.setMesh(m)
        f2.setArray(DataArrayDouble([11.,7.],2,1))
        f2.setName(fieldName2)
        WriteFieldUsingAlreadyWrittenMesh(fname,f2)
        f3=MEDCouplingFieldDouble(ON_NODES) ; f3.setTime(*t) ; f3.setMesh(m)
        f3.setArray(DataArrayDouble([4.,2.,1.,4.,2.,1.],6,1))
        f3.setName(fieldName3)
        WriteFieldUsingAlreadyWrittenMesh(fname,f3)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 2):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        vExp0=[DataArrayDouble([7.,11.]),DataArrayDouble([11.,7.])]
        vExp1=[DataArrayDouble([3.,5.,7.,6.,2.,3.,11.,8.]),DataArrayDouble([7.,6.,3.,5.,11.,8.,2.,3.])]
        for i in range(2):
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
        m.checkConsistency()
        #
        t=(1.1,0,-1)
        f=MEDCouplingFieldDouble(ON_GAUSS_PT) ; f.setTime(*t) ; f.setMesh(m)
        f.setGaussLocalizationOnType(NORM_QUAD4,[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.2,0.2,0.8,0.8],[0.7,0.3])
        f.setArray(DataArrayDouble([3.,5.,4.,6.])) ; f.getArray().setInfoOnComponents(["Smth"])
        f.setName(fieldName1)
        f.checkConsistencyLight()
        WriteField(fname,f,True)
        f2=MEDCouplingFieldDouble(ON_CELLS) ; f2.setTime(*t) ; f2.setMesh(m)
        f2.setArray(DataArrayDouble([7.,11.],2,1))
        f2.setName(fieldName2)
        WriteFieldUsingAlreadyWrittenMesh(fname,f2)
        f3=MEDCouplingFieldDouble(ON_NODES) ; f3.setTime(*t) ; f3.setMesh(m)
        f3.setArray(DataArrayDouble([1.,2.,4.,1.,2.,4.],6,1))
        f3.setName(fieldName3)
        WriteFieldUsingAlreadyWrittenMesh(fname,f3)
        #
        t=(2.1,1,-1)
        f=MEDCouplingFieldDouble(ON_GAUSS_PT) ; f.setTime(*t) ; f.setMesh(m)
        f.setGaussLocalizationOnType(NORM_QUAD4,[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.2,0.2,0.8,0.8],[0.7,0.3])
        f.setArray(DataArrayDouble([5.,3.,6.,4.])) ; f.getArray().setInfoOnComponents(["Smth"])
        f.setName(fieldName1)
        f.checkConsistencyLight()
        WriteFieldUsingAlreadyWrittenMesh(fname,f)
        f2=MEDCouplingFieldDouble(ON_CELLS) ; f2.setTime(*t) ; f2.setMesh(m)
        f2.setArray(DataArrayDouble([11.,7.],2,1))
        f2.setName(fieldName2)
        WriteFieldUsingAlreadyWrittenMesh(fname,f2)
        f3=MEDCouplingFieldDouble(ON_NODES) ; f3.setTime(*t) ; f3.setMesh(m)
        f3.setArray(DataArrayDouble([4.,2.,1.,4.,2.,1.],6,1))
        f3.setName(fieldName3)
        WriteFieldUsingAlreadyWrittenMesh(fname,f3)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertEqual([NORM_QUAD4],fcscp.getGeoTypesAt(0,ms[0]))
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 2):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        vExp0=[DataArrayDouble([7.,11.]),DataArrayDouble([11.,7.])]
        vExp1=[DataArrayDouble([3.,5.,4.,6.]),DataArrayDouble([5.,3.,6.,4.])]
        for i in range(2):
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
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname) # here all is read, the SauvReader (or other Reader) is emulated
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 2):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        vExp0=[DataArrayDouble([7.,11.]),DataArrayDouble([11.,7.])]
        vExp1=[DataArrayDouble([3.,5.,4.,6.]),DataArrayDouble([5.,3.,6.,4.])]
        for i in range(2):
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
        f.checkConsistencyLight()
        f1ts=MEDFileField1TS() ; f1ts.setFieldNoProfileSBT(f)
        ff.pushBackTimeStep(f1ts)
        # time 1
        t=(2.1,2,-2)
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setTime(*t) ; f.setMesh(m1)
        f.setName(fieldName)
        arr=DataArrayDouble(24) ; arr.iota() ; arr.reverse() ; arr.setInfoOnComponents(["AStr"])
        f.setArray(arr)
        f.checkConsistencyLight()
        f1ts=MEDFileField1TS() ; f1ts.setFieldNoProfileSBT(f)
        ff.pushBackTimeStep(f1ts)
        #
        mm.write(fname,2)
        ff.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        a10,a11=mml2.retrieveFamilyIdsOnNodes()
        self.assertTrue(not a10)
        self.assertTrue(a11) # no copy here
        a12,a13=mml2.retrieveNumberIdsOnNodes()
        self.assertTrue(not a12)
        self.assertTrue(a13) # no copy here
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 2):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in range(2):
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

    def test20(self):
        """ This test works with groups/families on cells AND on nodes. Here 4 fields each on same time steps (2).
        1 field on CELLS without profile, 1 field on CELLS with profile, 1 field on NODES without profile, 1 field on NODES with profile.
        All of these 4 fields lies on a single mesh "mesh". The 2 fields on profile lies on a same support.
        One drawback of this test : no multi geom type. Coming soon !
        """
        fname="ForMEDReader20.med"
        fieldName0="ANodeField"
        fieldName1="ACellField"
        fieldName2="ANodeFieldPfl"
        fieldName3="ACellFieldPfl"
        pfl2=DataArrayInt([5,6,7,10,11,12,15,16,17,20,21,22]) ; pfl2.setName("pfl2")
        pfl3=DataArrayInt([4,5,8,9,12,13]) ; pfl3.setName("pfl3")
        #
        arr=DataArrayDouble(5) ; arr.iota()
        m=MEDCouplingCMesh("mesh") ; m.setCoords(arr,arr)
        m=m.buildUnstructured()
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        fs=MEDFileFields()
        fmts0=MEDFileFieldMultiTS() ; fs.pushField(fmts0)
        fmts0.setDtUnit("s")
        fmts1=MEDFileFieldMultiTS() ; fs.pushField(fmts1)
        fmts1.setDtUnit("s")
        fmts2=MEDFileFieldMultiTS() ; fs.pushField(fmts2)
        fmts2.setDtUnit("s")
        fmts3=MEDFileFieldMultiTS() ; fs.pushField(fmts3)
        fmts3.setDtUnit("s")
        ####
        t=(1.1,0,-2)
        f0=MEDCouplingFieldDouble(ON_NODES) ; f0.setMesh(m)
        f0.setName(fieldName0) ; f0.setTime(*t)
        da=m.getCoords().magnitude() ; da.setInfoOnComponents(["zeInfo"])
        f0.setArray(da)
        f0.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldNoProfileSBT(f0)
        fmts0.pushBackTimeStep(f1ts)
        #
        f1=MEDCouplingFieldDouble(ON_CELLS) ; f1.setMesh(m)
        f1.setName(fieldName1) ; f1.setTime(*t)
        da=m.computeCellCenterOfMass().magnitude() ; da.setInfoOnComponents(["zeInfoCell"])
        f1.setArray(da)
        f1.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldNoProfileSBT(f1)
        fmts1.pushBackTimeStep(f1ts)
        #
        f2=MEDCouplingFieldDouble(ON_NODES) ; mTmp=m[pfl3] ; mTmp.zipCoords() ; mTmp.setName(m.getName()) ; f2.setMesh(mTmp)
        f2.setName(fieldName2) ; f2.setTime(*t)
        da=m.getCoords().magnitude()[pfl2] ; da.setInfoOnComponents(["zzzz"])
        f2.setArray(da)
        f2.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldProfile(f2,mm,0,pfl2)
        fmts2.pushBackTimeStep(f1ts)
        #
        f3=MEDCouplingFieldDouble(ON_CELLS) ; mTmp=m[pfl3] ; mTmp.setName(m.getName()) ; f3.setMesh(mTmp)
        f3.setName(fieldName3) ; f3.setTime(*t)
        da=mTmp.computeCellCenterOfMass().magnitude() ; da.setInfoOnComponents(["abcdefg"])
        f3.setArray(da)
        f3.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldProfile(f3,mm,0,pfl3)
        fmts3.pushBackTimeStep(f1ts)
        ####
        t=(2.1,1,-3)
        f0=MEDCouplingFieldDouble(ON_NODES) ; f0.setMesh(m)
        f0.setName(fieldName0) ; f0.setTime(*t)
        da=m.getCoords().magnitude() ; da.reverse() ; da.setInfoOnComponents(["zeInfo"])
        f0.setArray(da)
        f0.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldNoProfileSBT(f0)
        fmts0.pushBackTimeStep(f1ts)
        #
        f1=MEDCouplingFieldDouble(ON_CELLS) ; f1.setMesh(m)
        f1.setName(fieldName1) ; f1.setTime(*t)
        da=m.computeCellCenterOfMass().magnitude() ; da.reverse() ; da.setInfoOnComponents(["zeInfoCell"])
        f1.setArray(da)
        f1.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldNoProfileSBT(f1)
        fmts1.pushBackTimeStep(f1ts)
        #
        f2=MEDCouplingFieldDouble(ON_NODES) ; mTmp=m[pfl3] ; mTmp.zipCoords() ; mTmp.setName(m.getName()) ; f2.setMesh(mTmp)
        f2.setName(fieldName2) ; f2.setTime(*t)
        da=m.getCoords().magnitude()[pfl2] ; da.reverse() ; da.setInfoOnComponents(["zzzz"])
        f2.setArray(da)
        f2.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldProfile(f2,mm,0,pfl2)
        fmts2.pushBackTimeStep(f1ts)
        #
        f3=MEDCouplingFieldDouble(ON_CELLS) ; mTmp=m[pfl3] ; mTmp.setName(m.getName()) ; f3.setMesh(mTmp)
        f3.setName(fieldName3) ; f3.setTime(*t)
        da=mTmp.computeCellCenterOfMass().magnitude() ; da.reverse() ; da.setInfoOnComponents(["abcdefg"])
        f3.setArray(da)
        f3.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldProfile(f3,mm,0,pfl3)
        fmts3.pushBackTimeStep(f1ts)
        ####
        grp1=DataArrayInt([6,7,8,11,12,13,16,17,18]) ; grp1.setName("grp1")
        grp2=DataArrayInt([10,11,15,16,20,21]) ; grp2.setName("grp2")
        mm.setGroupsAtLevel(1,[grp1,grp2])
        grp3=DataArrayInt([4,5,6]) ; grp3.setName("grp3")
        grp4=DataArrayInt([8,9,10]) ; grp4.setName("grp4")
        mm.setGroupsAtLevel(0,[grp3,grp4])
        d=DataArrayInt(25) ; d.iota() ; d*=10 ;  mm.setRenumFieldArr(1,d)
        d=DataArrayInt(16) ; d.iota() ; d*=11 ;  mm.setRenumFieldArr(0,d)
        mm.write(fname,2)
        fs.appendGlobs(fmts2,1e-12)
        fs.appendGlobs(fmts3,1e-12)
        fs.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),4)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),4)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[1][0]),2)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(DataArrayDouble([0.,0.,0.,1.,0.,0.,2.,0.,0.,3.,0.,0.,4.,0.,0.,0.,1.,0.,1.,1.,0.,2.,1.,0.,3.,1.,0.,4.,1.,0.,0.,2.,0.,1.,2.,0.,2.,2.,0.,3.,2.,0.,4.,2.,0.,0.,3.,0.,1.,3.,0.,2.,3.,0.,3.,3.,0.,4.,3.,0.,0.,4.,0.,1.,4.,0.,2.,4.,0.,3.,4.,0.,4.,4.,0.],25,3),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,1,0,5,6,4,2,1,6,7,4,3,2,7,8,4,4,3,8,9,4,6,5,10,11,4,7,6,11,12,4,8,7,12,13,4,9,8,13,14,4,11,10,15,16,4,12,11,16,17,4,13,12,17,18,4,14,13,18,19,4,16,15,20,21,4,17,16,21,22,4,18,17,22,23,4,19,18,23,24])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([-5,-5,-5,-5,-6,-6,-6,-5,-7,-7,-7,-5,-5,-5,-5,-5])))
        self.assertTrue(a7) # no copy here
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([0,11,22,33,44,55,66,77,88,99,110,121,132,143,154,165])))
        self.assertTrue(a9) # no copy here
        a10,a11=mml2.retrieveFamilyIdsOnNodes()
        self.assertTrue(a10.isEqual(DataArrayInt([1,1,1,1,1,1,2,2,2,1,3,4,2,2,1,3,4,2,2,1,3,3,1,1,1])))
        self.assertTrue(a11) # no copy here
        a12,a13=mml2.retrieveNumberIdsOnNodes()
        self.assertTrue(a12.isEqual(DataArrayInt([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240])))
        self.assertTrue(a13) # no copy here
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 2):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in range(2):
            f=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName1)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([0.7071067811865476,1.5811388300841898,2.5495097567963922,3.5355339059327378,1.5811388300841898,2.1213203435596424,2.9154759474226504,3.8078865529319543,2.5495097567963922,2.9154759474226504,3.5355339059327378,4.301162633521313,3.5355339059327378,3.8078865529319543,4.301162633521313,4.949747468305833])
            if i==1: vExp.reverse()
            vExp.setInfoOnComponents(["zeInfoCell"])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport1[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([0.,1.,2.,3.,4.,1.,1.4142135623730951,2.23606797749979,3.1622776601683795,4.123105625617661,2.,2.23606797749979,2.8284271247461903,3.605551275463989,4.47213595499958,3.,3.1622776601683795,3.605551275463989,4.242640687119285,5.,4.,4.123105625617661,4.47213595499958,5.,5.656854249492381])
            if i==1: vExp.reverse()
            vExp.setInfoOnComponents(["zeInfo"])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        ### Testing the 2nd support
        fcscp=allFMTSLeavesPerCommonSupport1[1][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(DataArrayDouble([0.,1.,0.,1.,1.,0.,2.,1.,0.,0.,2.,0.,1.,2.,0.,2.,2.,0.,0.,3.,0.,1.,3.,0.,2.,3.,0.,0.,4.,0.,1.,4.,0.,2.,4.,0.],12,3),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5,10,15,20,25])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,1,0,3,4,4,2,1,4,5,4,4,3,6,7,4,5,4,7,8,4,7,6,9,10,4,8,7,10,11])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([-6,-6,-7,-7,-5,-5])))
        self.assertTrue(not a7) # copy here
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([44,55,88,99,132,143])))
        self.assertTrue(not a9) # copy here
        a10,a11=mml2.retrieveFamilyIdsOnNodes()
        self.assertTrue(a10.isEqual(DataArrayInt([1,2,2,3,4,2,3,4,2,3,3,1])))
        self.assertTrue(not a11) # copy here
        a12,a13=mml2.retrieveNumberIdsOnNodes()
        self.assertTrue(a12.isEqual(DataArrayInt([50,60,70,100,110,120,150,160,170,200,210,220])))
        self.assertTrue(not a13) # copy here
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(2):
            f=allFMTSLeavesPerCommonSupport1[1][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName3)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([1.5811388300842,2.1213203435596,2.5495097567964,2.9154759474227,3.5355339059327,3.807886552932])
            if i==1: vExp.reverse()
            vExp.setInfoOnComponents(["abcdefg"])
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport1[1][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName2)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble([1.,1.4142135623731,2.2360679774998,2.,2.2360679774998,2.8284271247462,3.,3.1622776601684,3.605551275464,4.,4.1231056256177,4.4721359549996])
            if i==1: vExp.reverse()
            vExp.setInfoOnComponents(["zzzz"])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass

    def test21(self):
        """ Here the created MED file contains only a mesh. The aim here is to test capability of MEDReader to support no fields.
        This test checks nothing but write a MED file to be used by MEDReader tests.
        """
        fname="ForMEDReader21.med"
        mm=MEDFileUMesh()
        #
        m0=MEDCouplingCMesh("mesh") ; arr=DataArrayDouble(5) ; arr.iota() ; m0.setCoords(arr,arr) ; m0=m0.buildUnstructured()
        mm.setMeshAtLevel(0,m0)
        grp0=DataArrayInt([5,6,9,10]) ; grp0.setName("Inside2D")
        grp1=DataArrayInt([0,1,2,3,4,7,8,11,12,13,14,15]) ; grp1.setName("Border2D")
        grp2=DataArrayInt([2,3,6,7]) ; grp2.setName("LowerRight2D")
        mm.setGroupsAtLevel(0,[grp0,grp1,grp2])
        #
        m1=MEDCouplingUMesh(m0.getName(),1) ; m1.setCoords(m0.getCoords()) ; m1.allocateCells()
        for elt in [[0,1],[1,2],[2,3],[3,4],[4,9],[9,14],[14,19],[19,24],[24,23],[23,22],[22,21],[21,20],[20,15],[15,10],[10,5],[5,0],[2,7],[7,12],[12,17],[17,22],
                    [10,11],[11,12],[12,13],[13,14]]:
            m1.insertNextCell(NORM_SEG2,elt)
            pass
        mm.setMeshAtLevel(-1,m1)
        grp4=DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]) ; grp4.setName("Border1D")
        grp5=DataArrayInt([16,17,18,19,20,21,22,23]) ; grp5.setName("Inside1D")
        grp6=DataArrayInt([18,19,22,23]) ; grp6.setName("UpperRight1D")
        mm.setGroupsAtLevel(-1,[grp4,grp5,grp6])
        #
        grp7=DataArrayInt([1,2,3,6,7,8,11,12,13,16,17,18,21,22,23]) ; grp7.setName("InsideYNode")
        grp8=DataArrayInt([5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]) ; grp8.setName("InsideXNode")
        mm.setGroupsAtLevel(1,[grp7,grp8])
        #
        mm.write(fname,2)
        pass

    def test22(self):
        """ Use case where a field on nodes (ANodeField) on a mesh defined both in meshdim 2 and meshdim 1.
        The only possible geometrical support that suits the field is those with meshdim equal to 1 (-1 in relative).
        """
        fname="ForMEDReader22.med"
        fieldName0="ANodeField"
        mm=MEDFileUMesh()
        coo=DataArrayDouble([(4.,3.),(7.,3.),(2.,5.),(6.,5.),(9.,5.),(4.,7.),(8.,7.),(3.,8.),(9.,8.)])
        m0=MEDCouplingUMesh("mesh",2) ; m0.setCoords(coo) ; m0.allocateCells() ; m0.insertNextCell(NORM_TRI3,[2,3,0]) ; m0.insertNextCell(NORM_TRI3,[3,1,0]) ; m0.insertNextCell(NORM_TRI3,[3,4,1])
        mm.setMeshAtLevel(0,m0)
        m1=MEDCouplingUMesh("mesh",1) ; m1.setCoords(coo) ; m1.allocateCells() ; m1.insertNextCell(NORM_SEG2,[2,0]) ;  m1.insertNextCell(NORM_SEG2,[0,1]) ; m1.insertNextCell(NORM_SEG2,[1,4])
        m1.insertNextCell(NORM_SEG2,[3,5]) ; m1.insertNextCell(NORM_SEG2,[5,7]) ; m1.insertNextCell(NORM_SEG2,[3,6]) ; m1.insertNextCell(NORM_SEG2,[6,8])
        mm.setMeshAtLevel(-1,m1)
        fs=MEDFileFields()
        fmts0=MEDFileFieldMultiTS() ; fs.pushField(fmts0)
        fmts0.setDtUnit("s")
        #
        t=(1.1,0,-2)
        f0=MEDCouplingFieldDouble(ON_NODES) ; f0.setMesh(m1)
        f0.setName(fieldName0) ; f0.setTime(*t)
        da=DataArrayDouble(9) ; da.iota() ; da.setInfoOnComponents(["zeInfo"])
        f0.setArray(da)
        f0.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldNoProfileSBT(f0)
        fmts0.pushBackTimeStep(f1ts)
        #
        t=(2.1,1,-3)
        f0=MEDCouplingFieldDouble(ON_NODES) ; f0.setMesh(m1)
        f0.setName(fieldName0) ; f0.setTime(*t)
        da=DataArrayDouble(9) ; da.iota() ; da.reverse() ; da.setInfoOnComponents(["zeInfo"])
        f0.setArray(da)
        f0.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldNoProfileSBT(f0)
        fmts0.pushBackTimeStep(f1ts)
        #
        mm.write(fname,2)
        fs.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertEqual([NORM_TRI3,NORM_SEG2],fcscp.getGeoTypesAt(0,ms[0]))#contains all cell types of underlying mesh because only nodes with no profiles
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(DataArrayDouble([(4.,3.,0.),(7.,3.,0.),(2.,5.,0.),(6.,5.,0.),(9.,5.,0.),(4.,7.,0.),(8.,7.,0.),(3.,8.,0.),(9.,8.,0.)]),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,5,5,3,3,3,3,3,3,3])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,8,12,15,18,21,24,27,30])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,2,3,0,3,3,1,0,3,3,4,1,2,2,0,2,0,1,2,1,4,2,3,5,2,5,7,2,3,6,2,6,8])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,0,0,0,0,0,0,0,0,0])))
        self.assertTrue(not a7) # copy here
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(not a8)
        self.assertTrue(a9) # nocopy here
        a10,a11=mml2.retrieveFamilyIdsOnNodes()
        self.assertTrue(not a10)
        self.assertTrue(a11) # no copy here
        a12,a13=mml2.retrieveNumberIdsOnNodes()
        self.assertTrue(not a12)
        self.assertTrue(a13) # no copy here
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        #
        f=allFMTSLeavesPerCommonSupport1[0][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
        f.loadArraysIfNecessary()
        v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
        self.assertEqual(f.getName(),fieldName0)
        self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
        vExp=DataArrayDouble(9) ; vExp.iota() ; vExp.setInfoOnComponents(["zeInfo"])
        self.assertTrue(v.isEqual(vExp,1e-12))
        #
        f=allFMTSLeavesPerCommonSupport1[0][0][0][1]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
        f.loadArraysIfNecessary()
        v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
        self.assertEqual(f.getName(),fieldName0)
        self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
        vExp=DataArrayDouble(9) ; vExp.iota() ; vExp.setInfoOnComponents(["zeInfo"]) ; vExp.reverse()
        self.assertTrue(v.isEqual(vExp,1e-12))
        pass 
    
    def test23(self):
        """ Non regression test 2219 of modes. Idem than test22 except that here the node field is on profile.
        """
        fname="ForMEDReader23.med"
        fieldName0="ANodeField"
        mm=MEDFileUMesh()
        coo=DataArrayDouble([(4.,3.),(7.,3.),(2.,5.),(6.,5.),(9.,5.),(4.,7.),(8.,7.),(3.,8.),(9.,8.)])
        m0=MEDCouplingUMesh("mesh",2) ; m0.setCoords(coo) ; m0.allocateCells() ; m0.insertNextCell(NORM_TRI3,[2,3,0]) ; m0.insertNextCell(NORM_TRI3,[3,1,0]) ; m0.insertNextCell(NORM_TRI3,[3,4,1])
        mm.setMeshAtLevel(0,m0)
        m1=MEDCouplingUMesh("mesh",1) ; m1.setCoords(coo) ; m1.allocateCells() ; m1.insertNextCell(NORM_SEG2,[2,0]) ;  m1.insertNextCell(NORM_SEG2,[0,1]) ; m1.insertNextCell(NORM_SEG2,[1,4])
        m1.insertNextCell(NORM_SEG2,[3,5]) ; m1.insertNextCell(NORM_SEG2,[5,7]) ; m1.insertNextCell(NORM_SEG2,[3,6]) ; m1.insertNextCell(NORM_SEG2,[6,8])
        mm.setMeshAtLevel(-1,m1)
        fmts0=MEDFileFieldMultiTS()
        fmts0.setDtUnit("s")
        #
        pfl=DataArrayInt([0,1,2,4]) ; pfl.setName("pfl")
        pflCell=DataArrayInt([0,1,2]) ; m1Part=m1[pflCell] ; m1Part.zipCoords()
        #
        t=(1.1,0,-2)
        f0=MEDCouplingFieldDouble(ON_NODES) ; f0.setMesh(m1Part)
        f0.setName(fieldName0) ; f0.setTime(*t)
        da=DataArrayDouble(4) ; da.iota() ; da.setInfoOnComponents(["zeInfo"])
        f0.setArray(da)
        f0.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldProfile(f0,mm,-1,pfl)
        fmts0.pushBackTimeStep(f1ts)
        #
        t=(2.1,1,-3)
        f0=MEDCouplingFieldDouble(ON_NODES) ; f0.setMesh(m1Part)
        f0.setName(fieldName0) ; f0.setTime(*t)
        da=DataArrayDouble(4) ; da.iota() ; da.reverse() ; da.setInfoOnComponents(["zeInfo"])
        f0.setArray(da)
        f0.checkConsistencyLight()
        f1ts=MEDFileField1TS()
        f1ts.setFieldProfile(f0,mm,-1,pfl)
        fmts0.pushBackTimeStep(f1ts)
        mm.write(fname,2)
        fmts0.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(DataArrayDouble([(4.,3.,0.),(7.,3.,0.),(2.,5.,0.),(9.,5.,0.)]),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([3,3,3])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,3,6])))
        self.assertTrue(a3.isEqual(DataArrayInt([2,2,0,2,0,1,2,1,3])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([0,0,0])))
        self.assertTrue(not a7) # copy here
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(not a8)
        self.assertTrue(a9) # nocopy here
        a10,a11=mml2.retrieveFamilyIdsOnNodes()
        self.assertTrue(not a10)
        self.assertTrue(a11) # no copy here
        a12,a13=mml2.retrieveNumberIdsOnNodes()
        self.assertTrue(not a12)
        self.assertTrue(a13) # no copy here
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        #
        f=allFMTSLeavesPerCommonSupport1[0][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
        f.loadArraysIfNecessary()
        v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
        self.assertEqual(f.getName(),fieldName0)
        vExp=DataArrayDouble(4) ; vExp.iota() ; vExp.setInfoOnComponents(["zeInfo"])
        self.assertTrue(v.isEqual(vExp,1e-12))
        #
        f=allFMTSLeavesPerCommonSupport1[0][0][0][1]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
        f.loadArraysIfNecessary()
        v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
        self.assertEqual(f.getName(),fieldName0)
        vExp=DataArrayDouble(4) ; vExp.iota() ; vExp.setInfoOnComponents(["zeInfo"]) ; vExp.reverse()
        self.assertTrue(v.isEqual(vExp,1e-12))
        pass

    def test24(self):
        """ Non regression test for cartesian mesh whose the 3rd direction has only one node. It a false 3D mesh.
        """
        fname="ForMEDReader24.med"
        fieldName0="zeFieldNode"
        cmesh=MEDCouplingCMesh("mesh")
        arr0=DataArrayDouble([0.,1.1,2.2,3.3,4.4])
        arr1=DataArrayDouble([0.,1.4,2.3])
        arr2=DataArrayDouble([5.])
        cmesh.setCoords(arr0,arr1,arr2)
        fmts0=MEDFileFieldMultiTS()
        fmts0.setDtUnit("s")
        #
        t=(1.1,2,3)
        f=MEDCouplingFieldDouble(ON_NODES) ; f.setName(fieldName0)
        f.setMesh(cmesh)
        arr=DataArrayDouble(15) ; arr.setInfoOnComponents(["tutu"]) ; arr.iota()
        f.setArray(arr)
        f.setTime(*t)
        f1ts=MEDFileField1TS()
        f1ts.setFieldNoProfileSBT(f)
        fmts0.pushBackTimeStep(f1ts)
        #
        t=(3.3,4,5)
        arr=DataArrayDouble(15) ; arr.setInfoOnComponents(["tutu"]) ; arr.iota()
        arr.reverse()
        f.setArray(arr)
        f.setTime(*t)
        f1ts=MEDFileField1TS()
        f1ts.setFieldNoProfileSBT(f)
        fmts0.pushBackTimeStep(f1ts)
        #
        mm=MEDFileCMesh() ; mm.setMesh(cmesh)
        mm.write(fname,2)
        fmts0.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),1)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),1)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport[0][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDCMeshMultiLev))
        (a,b,c),d=mml2.buildVTUArrays()
        self.assertTrue(d)#d is True because the a,b and c are directly those in the internal data structure
        self.assertTrue(a.isEqual(arr0,1e-12))
        self.assertTrue(b.isEqual(arr1,1e-12))
        self.assertTrue(c.isEqual(arr2,1e-12))
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(2):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(15) ; vExp.iota(0) ; vExp.setInfoOnComponents(["tutu"])
            if i==1:
                vExp.reverse()
                pass
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass

    def test25(self):
        """ A tricky test that reproduces an invalid behaviour
        Here a same field is defined both on CELLS and GAUSS_PT, with a profile for each.
        The problem appears on array computation when performing CELLS then GAUSS_PT and CELLS again.
        """
        fname="ForMEDReader25.med"
        m=MEDFileUMesh()
        coords=DataArrayDouble([0.,0.,1.,0.,2.,0.,0.,1.,1.,1.,2.,1.,0.,2.,1.,2.,2.,2.,0.,3.,1.,3.,2.,3.,1.,4.,1.,5.,1.,6.],15,2)
        m0=MEDCouplingUMesh("mesh",2) ; m0.setCoords(coords)
        m0.allocateCells()
        m0.insertNextCell(NORM_QUAD4,[0,3,4,1])
        m0.insertNextCell(NORM_QUAD4,[1,4,5,2])
        m0.insertNextCell(NORM_QUAD4,[3,6,7,4])
        m0.insertNextCell(NORM_QUAD4,[4,7,8,5])
        m0.insertNextCell(NORM_QUAD4,[6,9,10,7])
        m0.insertNextCell(NORM_QUAD4,[7,10,11,8])
        m.setMeshAtLevel(0,m0)
        m1=MEDCouplingUMesh("mesh",1) ; m1.setCoords(coords)
        m1.allocateCells()
        m1.insertNextCell(NORM_SEG2,[10,12])
        m1.insertNextCell(NORM_SEG2,[12,13])
        m1.insertNextCell(NORM_SEG2,[13,14])
        m.setMeshAtLevel(-1,m1)
        m.setFamilyFieldArr(0,DataArrayInt([-1,-2,-3,-4,-5,-6]))
        m.setFamilyFieldArr(-1,DataArrayInt([-7,-8,-9]))
        m.setFamilyFieldArr(1,DataArrayInt([3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]))
        m.setRenumFieldArr(0,DataArrayInt([101,102,103,104,105,106]))
        m.setRenumFieldArr(-1,DataArrayInt([107,108,109]))
        m.setRenumFieldArr(1,DataArrayInt([203,204,205,206,207,208,209,210,211,212,213,214,215,216,217]))
        #
        fmts=MEDFileFieldMultiTS()
        info0=["aa","bbb"]
        name0="zeField"
        pflName0="pfl"
        pflName1="pfl2"
        #
        f1ts=MEDFileField1TS()
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setName(name0)
        arr=DataArrayDouble([(-1,-11),(-2,-22)]) ; arr.setInfoOnComponents(info0)
        f.setArray(arr)
        pfl0=DataArrayInt([0,1]) ; pfl0.setName(pflName0)
        f1ts.setFieldProfile(f,m,-1,pfl0)
        del f
        f2=MEDCouplingFieldDouble(ON_GAUSS_PT) ; f2.setName(name0)
        arr=DataArrayDouble(15) ; arr.iota(1)
        arr=DataArrayDouble.Meld(arr,arr+10) ; arr.setInfoOnComponents(info0)
        f2.setArray(arr)
        pfl1=DataArrayInt([1,3,5]) ; pfl1.setName(pflName1)
        tmp=m0[pfl1] ; f2.setMesh(tmp)
        f2.setGaussLocalizationOnType(NORM_QUAD4,[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[-0.5,-0.5,0.5,-0.5,0.5,0.5,-0.5,0.5,0.,0.],[0.1,0.1,0.1,0.1,0.6])
        f2.checkConsistencyLight()
        f1ts.setFieldProfile(f2,m,0,pfl1)
        fmts.pushBackTimeStep(f1ts)
        #
        m.write(fname,2)
        fmts.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False) # false is absolutely necessary for the test
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),2)
        ### here the test is important !!! Pointers must be different !
        self.assertTrue(allFMTSLeavesToDisplay[0][0][0].getUndergroundDataArray().getHiddenCppPointer()!=allFMTSLeavesToDisplay[0][1][0].getUndergroundDataArray().getHiddenCppPointer())
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),2)
        ### here the test is important !!! Pointers must be different !
        self.assertTrue(allFMTSLeavesToDisplay[0][0][0].getUndergroundDataArray().getHiddenCppPointer()!=allFMTSLeavesToDisplay[0][1][0].getUndergroundDataArray().getHiddenCppPointer())
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertTrue(allFMTSLeavesToDisplay[0][0][0].getUndergroundDataArray().getHiddenCppPointer()!=allFMTSLeavesToDisplay[0][1][0].getUndergroundDataArray().getHiddenCppPointer())
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[1][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        # emulate first click
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        self.assertEqual([NORM_SEG2],fcscp.getGeoTypesAt(0,ms[0]))
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc) # copy here because 2D -> 3D
        expCoords=coords.changeNbOfComponents(3,0.)
        self.assertTrue(a0.isEqual(expCoords,1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([3,3])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,3])))
        self.assertTrue(a3.isEqual(DataArrayInt([2,10,12,2,12,13])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([-7,-8])))
        self.assertTrue(not a7) # copy here because profile on cells
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([107,108])))
        self.assertTrue(not a9) # copy here because profile on cells
        a10,a11=mml2.retrieveFamilyIdsOnNodes()
        self.assertTrue(a10.isEqual(DataArrayInt([3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])))
        self.assertTrue(a11) # no copy here
        a12,a13=mml2.retrieveNumberIdsOnNodes()
        self.assertTrue(a12.isEqual(DataArrayInt([203,204,205,206,207,208,209,210,211,212,213,214,215,216,217])))
        self.assertTrue(a13) # no copy here
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        fff0=allFMTSLeavesPerCommonSupport1[0][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(fff0,mst)
        fff0.loadArraysIfNecessary()
        self.assertEqual([ON_CELLS],fff0.getTypesOfFieldAvailable())
        v=mml.buildDataArray(fsst,fields,fff0.getUndergroundDataArray())
        self.assertEqual(fff0.getName(),name0)
        self.assertEqual(v.getHiddenCppPointer(),fff0.getUndergroundDataArray().getHiddenCppPointer())
        vExp=DataArrayDouble([(-1,-11),(-2,-22)]) ; vExp.setInfoOnComponents(info0)
        self.assertTrue(v.isEqual(vExp,1e-12))
        del fff0
        # emulate second click
        fcscp=allFMTSLeavesPerCommonSupport1[1][1]
        self.assertEqual([NORM_QUAD4],fcscp.getGeoTypesAt(0,ms[0]))
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc) # copy here because 2D -> 3D
        expCoords=coords.changeNbOfComponents(3,0.)
        self.assertTrue(a0.isEqual(expCoords,1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5,10])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,1,4,5,2,4,4,7,8,5,4,7,10,11,8])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([-2,-4,-6])))
        self.assertTrue(not a7) # copy here because profile on cells
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([102,104,106])))
        self.assertTrue(not a9) # copy here because profile on cells
        a10,a11=mml2.retrieveFamilyIdsOnNodes()
        self.assertTrue(a10.isEqual(DataArrayInt([3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])))
        self.assertTrue(a11) # no copy here
        a12,a13=mml2.retrieveNumberIdsOnNodes()
        self.assertTrue(a12.isEqual(DataArrayInt([203,204,205,206,207,208,209,210,211,212,213,214,215,216,217])))
        self.assertTrue(a13) # no copy here
        fff1=allFMTSLeavesPerCommonSupport1[1][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(fff1,mst)
        fff1.loadArraysIfNecessary()
        self.assertEqual([ON_GAUSS_PT],fff1.getTypesOfFieldAvailable())
        v=mml.buildDataArray(fsst,fields,fff1.getUndergroundDataArray())
        self.assertEqual(fff1.getName(),name0)
        self.assertEqual(v.getHiddenCppPointer(),fff1.getUndergroundDataArray().getHiddenCppPointer())
        vExp=DataArrayDouble([1.,11.,2.,12.,3.,13.,4.,14.,5.,15.,6.,16.,7.,17.,8.,18.,9.,19.,10.,20.,11.,21.,12.,22.,13.,23.,14.,24.,15.,25.],15,2) ; vExp.setInfoOnComponents(info0)
        self.assertTrue(v.isEqual(vExp,1e-12))
        # emulate third click
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc) # copy here because 2D -> 3D
        expCoords=coords.changeNbOfComponents(3,0.)
        self.assertTrue(a0.isEqual(expCoords,1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([3,3])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,3])))
        self.assertTrue(a3.isEqual(DataArrayInt([2,10,12,2,12,13])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([-7,-8])))
        self.assertTrue(not a7) # copy here because profile on cells
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8.isEqual(DataArrayInt([107,108])))
        self.assertTrue(not a9) # copy here because profile on cells
        a10,a11=mml2.retrieveFamilyIdsOnNodes()
        self.assertTrue(a10.isEqual(DataArrayInt([3,4,5,6,7,8,9,10,11,12,13,14,15,16,17])))
        self.assertTrue(a11) # no copy here
        a12,a13=mml2.retrieveNumberIdsOnNodes()
        self.assertTrue(a12.isEqual(DataArrayInt([203,204,205,206,207,208,209,210,211,212,213,214,215,216,217])))
        self.assertTrue(a13) # no copy here
        fff0=allFMTSLeavesPerCommonSupport1[0][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(fff0,mst)
        fff0.loadArraysIfNecessary()
        self.assertEqual([ON_CELLS],fff0.getTypesOfFieldAvailable())
        v=mml.buildDataArray(fsst,fields,fff0.getUndergroundDataArray())
        self.assertEqual(fff0.getName(),name0)
        self.assertEqual(v.getHiddenCppPointer(),fff0.getUndergroundDataArray().getHiddenCppPointer())
        vExp=DataArrayDouble([(-1,-11),(-2,-22)]) ; vExp.setInfoOnComponents(info0)
        self.assertTrue(v.isEqual(vExp,1e-12)) # <- THE test is here !!!
        del fff0
        pass

    def test26(self):
        """ Test focused on field on nodes (here f0Node and f1Node) lying on a profile of nodes that do not match perfectly a sub set of cells of its underlying mesh. See bug EDF 2405 and 2177.
        For this type of fields the support will contain only vertices.
        """
        fname="ForMEDReader26.med"
        coords=DataArrayDouble([(0.,0.,0.),(1.,0.,0.),(2.,0.,0.),(3.,0.,0.),(0.,1.,0.),(1.,1.,0.),(2.,1.,0.),(3.,1.,0.),(0.,2.,0.),(1.,2.,0.),(2.,2.,0.),(3.,2.,0.),(0.,3.,0.),(1.,3.,0.),(2.,3.,0.),(3.,3.,0.)])
        m0=MEDCouplingUMesh("mesh",2)
        m0.allocateCells()
        for elt in [[2,6,3],[6,7,3],[9,6,5],[9,10,6]]:
            m0.insertNextCell(NORM_TRI3,elt)
            pass
        for elt in [[0,4,5,1],[1,5,6,2],[4,8,9,5],[6,10,11,7],[8,12,13,9],[9,13,14,10],[10,14,15,11]]:
            m0.insertNextCell(NORM_QUAD4,elt)
            pass
        m0.setCoords(coords)
        ##
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m0)
        mm.setFamilyFieldArr(0,DataArrayInt([-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11]))
        mm.setFamilyFieldArr(1,DataArrayInt([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]))
        #
        f1ts0Node=MEDFileField1TS()
        f1ts1Node=MEDFileField1TS()
        f1ts2Cell=MEDFileField1TS()
        f1ts3Cell=MEDFileField1TS()
        f1ts4Cell=MEDFileField1TS()
        f1ts5Node=MEDFileField1TS()
        #
        pfl0=DataArrayInt([4,5,6,8,9,12]) ; pfl0.setName("pfl0")
        pfl1=DataArrayInt([0,1,4,5,7,10]) ; pfl1.setName("pfl1")
        pfl2=DataArrayInt([0,1,2,3,4,5,6,7,10,11,14,15]) ; pfl2.setName("pfl2")
        #
        f0Node=MEDCouplingFieldDouble(ON_NODES) ; f0Node.setName("f0Node")
        arr0=DataArrayDouble(6) ; arr0.iota()
        f0Node.setArray(arr0)
        f1ts0Node.setFieldProfile(f0Node,mm,0,pfl0)
        #
        f1Node=MEDCouplingFieldDouble(ON_NODES) ; f1Node.setName("f1Node")
        arr1=DataArrayDouble(6) ; arr1.iota() ; arr1.reverse()
        f1Node.setArray(arr1)
        f1ts1Node.setFieldProfile(f1Node,mm,0,pfl0)
        #
        f2Cell=MEDCouplingFieldDouble(ON_CELLS) ; f2Cell.setName("f2Cell")
        arr2=DataArrayDouble([2,3,0,1,4,5])
        f2Cell.setArray(arr2)
        f1ts2Cell.setFieldProfile(f2Cell,mm,0,pfl1)
        #
        f3Cell=MEDCouplingFieldDouble(ON_CELLS) ; f3Cell.setName("f3Cell")
        arr3=DataArrayDouble([5,4,3,2,1,0])
        f3Cell.setArray(arr3)
        f1ts3Cell.setFieldProfile(f3Cell,mm,0,pfl1)
        #
        f4Cell=MEDCouplingFieldDouble(ON_CELLS) ; f4Cell.setName("f4Cell")
        arr4=DataArrayDouble([2,2,0,1,1,0])
        f4Cell.setArray(arr4)
        f1ts4Cell.setFieldProfile(f4Cell,mm,0,pfl1)
        #
        f5Node=MEDCouplingFieldDouble(ON_NODES) ; f5Node.setName("f5Node")
        arr5=DataArrayDouble([0,1,2,3,10,11,13,2,11,1,10,0])
        f5Node.setArray(arr5)
        f1ts5Node.setFieldProfile(f5Node,mm,0,pfl2)
        #
        fs=MEDFileFields()
        for f in [f1ts0Node,f1ts1Node,f1ts2Cell,f1ts3Cell,f1ts4Cell,f1ts5Node]:
            fmts=MEDFileFieldMultiTS()
            fmts.pushBackTimeStep(f)
            fs.pushField(fmts)
            pass
        mm.write(fname,2)
        fs.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),6)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),6)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),4)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[1][0]),2)# <- the smart one is here
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport1[1][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        self.assertEqual([3,4,0],mml2.getGeoTypes())
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(DataArrayDouble([0.,1.,0.,1.,1.,0.,2.,1.,0.,0.,2.,0.,1.,2.,0.,0.,3.,0.],6,3),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,9,1])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,9])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,4,2,1,4,0,3,4,1,1,5])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([-3,-7,13])))
        self.assertTrue(not a7) # copy here because profile on cells
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8 is None)
        self.assertTrue(a9) # no copy here because no number field
        a10,a11=mml2.retrieveFamilyIdsOnNodes()
        self.assertTrue(a10.isEqual(DataArrayInt([5,6,7,9,10,13])))
        self.assertTrue(not a11) # copy here
        a12,a13=mml2.retrieveNumberIdsOnNodes()
        self.assertTrue(a12 is None)
        self.assertTrue(a13) # no copy here because no number field
        #
        fff0=allFMTSLeavesPerCommonSupport1[1][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(fff0,mst)
        fff0.loadArraysIfNecessary()
        v=mml2.buildDataArray(fsst,fields,fff0.getUndergroundDataArray())
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        self.assertEqual(fff0.getName(),"f0Node")
        self.assertEqual(v.getHiddenCppPointer(),fff0.getUndergroundDataArray().getHiddenCppPointer())
        vExp=DataArrayDouble([0.,1.,2.,3.,4.,5.])
        self.assertTrue(v.isEqual(vExp,1e-12)) # <- THE test is here !!!
        #
        fff1=allFMTSLeavesPerCommonSupport1[1][0][1][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(fff1,mst)
        fff1.loadArraysIfNecessary()
        v=mml2.buildDataArray(fsst,fields,fff1.getUndergroundDataArray())
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        self.assertEqual(fff1.getName(),"f1Node")
        self.assertEqual(v.getHiddenCppPointer(),fff1.getUndergroundDataArray().getHiddenCppPointer())
        vExp=DataArrayDouble([5.,4.,3.,2.,1.,0.])
        self.assertTrue(v.isEqual(vExp,1e-12)) # <- THE test is here !!!
        pass

    def test27(self):
        """ This test defines 2 fields f0 and f1 on nodes lying on an unstructured mesh with no cells.
        f0 is a field on all nodes. f1 is a partial field on nodes.
        """
        fname="ForMEDReader27.med"
        coords=DataArrayDouble([(0.,0.,0.),(1.,0.,0.),(2.,0.,0.),(3.,0.,0.),(0.,1.,0.),(1.,1.,0.),(2.,1.,0.),(3.,1.,0.),(0.,2.,0.),(1.,2.,0.),(2.,2.,0.),(3.,2.,0.),(0.,3.,0.),(1.,3.,0.),(2.,3.,0.),(3.,3.,0.)])
        m0=MEDCouplingUMesh("mesh",2)
        m0.allocateCells()
        m0.setCoords(coords)
        ##
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m0)
        mm.setFamilyFieldArr(1,DataArrayInt([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]))
        #
        f1ts0=MEDFileField1TS()
        f1ts1=MEDFileField1TS()
        #
        f0=MEDCouplingFieldDouble(ON_NODES) ; f0.setMesh(m0) ; f0.setName("f0NoPfl")
        arr0=DataArrayDouble([0.,1.,2.,3.,1.,1.5,2.2,3.1,2.,2.2,3.,3.1,3.,3.1,3.5,4.])
        f0.setArray(arr0)
        f0.checkConsistencyLight()
        f1ts0.setFieldNoProfileSBT(f0)
        self.assertEqual(f1ts0.getMeshName(),"mesh")
        #
        pfl1=DataArrayInt([0,1,2,3,4,5,6,8,9,12]) ; pfl1.setName("pfl1")
        f1=MEDCouplingFieldDouble(ON_NODES) ; f1.setName("f1Pfl")
        arr1=DataArrayDouble([3.,2.,1.,0.,2.,1.5,0.,1.,0.,0.2])
        f1.setArray(arr1)
        f1ts1.setFieldProfile(f1,mm,0,pfl1)
        self.assertEqual(f1ts1.getMeshName(),"mesh")
        #
        fs=MEDFileFields()
        fmts0=MEDFileFieldMultiTS()
        fmts0.pushBackTimeStep(f1ts0)
        fmts1=MEDFileFieldMultiTS()
        fmts1.pushBackTimeStep(f1ts1)
        fs.pushField(fmts0) ; fs.pushField(fmts1)
        self.assertEqual(fs[0].getMeshName(),"mesh")
        self.assertEqual(fs[1].getMeshName(),"mesh")
        mm.write(fname,2)
        fs.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        self.assertEqual(fields[0].getMeshName(),"mesh")
        self.assertEqual(fields[1].getMeshName(),"mesh")
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        self.assertEqual(fields_per_mesh[0][0].getMeshName(),"mesh")
        self.assertEqual(fields_per_mesh[0][1].getMeshName(),"mesh")
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),2)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),2)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[1][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(ncc)
        self.assertTrue(a0.isEqual(DataArrayDouble([(0.,0.,0.),(1.,0.,0.),(2.,0.,0.),(3.,0.,0.),(0.,1.,0.),(1.,1.,0.),(2.,1.,0.),(3.,1.,0.),(0.,2.,0.),(1.,2.,0.),(2.,2.,0.),(3.,2.,0.),(0.,3.,0.),(1.,3.,0.),(2.,3.,0.),(3.,3.,0.)]),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([])))
        self.assertTrue(a2.isEqual(DataArrayInt([])))
        self.assertTrue(a3.isEqual(DataArrayInt([])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        #
        fff0=allFMTSLeavesPerCommonSupport1[0][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(fff0,mst)
        fff0.loadArraysIfNecessary()
        v=mml2.buildDataArray(fsst,fields,fff0.getUndergroundDataArray())
        self.assertEqual(fff0.getName(),"f0NoPfl")
        self.assertEqual(v.getHiddenCppPointer(),fff0.getUndergroundDataArray().getHiddenCppPointer())
        vExp=DataArrayDouble([0.,1.,2.,3.,1.,1.5,2.2,3.1,2.,2.2,3.,3.1,3.,3.1,3.5,4])
        self.assertTrue(v.isEqual(vExp,1e-12))
        #
        fcscp=allFMTSLeavesPerCommonSupport1[1][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(DataArrayDouble([(0,0,0),(1,0,0),(2,0,0),(3,0,0),(0,1,0),(1,1,0),(2,1,0),(0,2,0),(1,2,0),(0,3,0)]),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([])))
        self.assertTrue(a2.isEqual(DataArrayInt([])))
        self.assertTrue(a3.isEqual(DataArrayInt([])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        fff1=allFMTSLeavesPerCommonSupport1[1][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(fff1,mst)
        fff1.loadArraysIfNecessary()
        v=mml2.buildDataArray(fsst,fields,fff1.getUndergroundDataArray())
        self.assertEqual(fff1.getName(),"f1Pfl")
        self.assertNotEqual(v.getHiddenCppPointer(),fff1.getUndergroundDataArray().getHiddenCppPointer()) # pointers are not equal because Profile
        vExp=DataArrayDouble([3.,2.,1.,0.,2.,1.5,0.,1.,0.,0.2])
        self.assertTrue(v.isEqual(vExp,1e-12))
        pass

    def test28(self):
        """ This test defines 2 fields f0,f1,f2,f3 lying on an unstructured mesh whith cells including NORM_POINT1.
        Both f0 and f1 are on NODES and f2 and f3 are on cells. f1 and f2 share the same support.
        f0 is on a nodal support that is not matchable with any cells (including NORM_POINT1)
        This test is a more aggressive version of test26.
        """
        fname="ForMEDReader28.med"
        coords=DataArrayDouble([(0.,0.,0.),(1.,0.,0.),(2.,0.,0.),(3.,0.,0.),(0.,1.,0.),(1.,1.,0.),(2.,1.,0.),(3.,1.,0.),(0.,2.,0.),(1.,2.,0.),(2.,2.,0.),(3.,2.,0.),(0.,3.,0.),(1.,3.,0.),(2.,3.,0.),(3.,3.,0.)])
        m0=MEDCouplingUMesh("mesh",2)
        m0.allocateCells()
        for elt in [[2,6,3],[6,7,3],[9,6,5],[9,10,6]]:
            m0.insertNextCell(NORM_TRI3,elt)
            pass
        for elt in [[0,4,5,1],[1,5,6,2],[4,8,9,5],[6,10,11,7],[8,12,13,9],[9,13,14,10],[10,14,15,11]]:
            m0.insertNextCell(NORM_QUAD4,elt)
            pass
        m0.setCoords(coords)
        m2=MEDCouplingUMesh("mesh",0) ; m2.setCoords(coords)
        m2.allocateCells()
        for elt in [[8],[13]]:
            m2.insertNextCell(NORM_POINT1,elt)
            pass
        ##
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m0)
        mm.setMeshAtLevel(-2,m2)
        mm.setFamilyFieldArr(0,DataArrayInt([-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11]))
        mm.setFamilyFieldArr(-2,DataArrayInt([-12,-13]))
        mm.setFamilyFieldArr(1,DataArrayInt([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]))
        #
        f1ts0Node=MEDFileField1TS()
        f1ts1Node=MEDFileField1TS()
        f1ts2Cell=MEDFileField1TS()
        f1ts3Cell=MEDFileField1TS()
        #
        pfl0=DataArrayInt([4,5,6,8,9,12]) ; pfl0.setName("pfl0")
        pfl1=DataArrayInt([0,1,4,5,7,10]) ; pfl1.setName("pfl1")
        pfl2=DataArrayInt([0,1,2,3,4,5,6,7,10,11,14,15]) ; pfl2.setName("pfl2")
        #
        f0Node=MEDCouplingFieldDouble(ON_NODES) ; f0Node.setName("f0Node")
        arr0=DataArrayDouble(6) ; arr0.iota()
        f0Node.setArray(arr0)
        f1ts0Node.setFieldProfile(f0Node,mm,0,pfl0)
        #
        f1Node=MEDCouplingFieldDouble(ON_NODES) ; f1Node.setName("f1Node")
        arr1=DataArrayDouble(12) ; arr1.iota() ; arr1.reverse()
        f1Node.setArray(arr1)
        f1ts1Node.setFieldProfile(f1Node,mm,0,pfl2)
        #
        f2Cell=MEDCouplingFieldDouble(ON_CELLS) ; f2Cell.setName("f2Cell")
        arr2=DataArrayDouble([2,3,0,1,4,5])
        f2Cell.setArray(arr2)
        f1ts2Cell.setFieldProfile(f2Cell,mm,0,pfl1)
        #
        f3Cell=MEDCouplingFieldDouble(ON_CELLS) ; f3Cell.setName("f3Cell")
        arr3=DataArrayDouble([5,4,3,2,1,0]) ; f3Cell.setArray(arr3)
        f1ts3Cell.setFieldProfile(f3Cell,mm,0,pfl1)
        f3Cell.setMesh(m2)
        arr3=DataArrayDouble([-1.1,-3.1]) ; f3Cell.setArray(arr3)
        f1ts3Cell.setFieldNoProfileSBT(f3Cell)
        #
        fs=MEDFileFields()
        for f in [f1ts0Node,f1ts1Node,f1ts2Cell,f1ts3Cell]:
            fmts=MEDFileFieldMultiTS()
            fmts.pushBackTimeStep(f)
            fs.pushField(fmts)
            pass
        mm.write(fname,2)
        fs.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),4)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),4)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),3)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[1][0]),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[2][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport1[2][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(DataArrayDouble([0.,1.,0.,1.,1.,0.,2.,1.,0.,0.,2.,0.,1.,2.,0.,0.,3.,0.],6,3),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,9,1,1])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,9,11])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,4,2,1,4,0,3,4,1,1,3,1,5])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6.isEqual(DataArrayInt([-3,-7,-12,13])))
        self.assertTrue(not a7) # copy here because profile on cells
        a8,a9=mml2.retrieveNumberIdsOnCells()
        self.assertTrue(a8 is None)
        self.assertTrue(a9) # no copy here because no number field
        a10,a11=mml2.retrieveFamilyIdsOnNodes()
        self.assertTrue(a10.isEqual(DataArrayInt([5,6,7,9,10,13])))
        self.assertTrue(not a11) # copy here
        a12,a13=mml2.retrieveNumberIdsOnNodes()
        self.assertTrue(a12 is None)
        self.assertTrue(a13) # no copy here because no number field
        #
        fff0=allFMTSLeavesPerCommonSupport1[2][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(fff0,mst)
        fff0.loadArraysIfNecessary()
        v=mml2.buildDataArray(fsst,fields,fff0.getUndergroundDataArray())
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        self.assertEqual(fff0.getName(),"f0Node")
        self.assertEqual(v.getHiddenCppPointer(),fff0.getUndergroundDataArray().getHiddenCppPointer())
        vExp=DataArrayDouble([0.,1.,2.,3.,4.,5.])
        self.assertTrue(v.isEqual(vExp,1e-12)) # <- THE test is here !!!
        ###
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(DataArrayDouble([(0,0,0),(1,0,0),(2,0,0),(3,0,0),(0,1,0),(1,1,0),(2,1,0),(3,1,0),(2,2,0),(3,2,0),(2,3,0),(3,3,0)]),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,5,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,8,13,18,23])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,2,6,3,3,6,7,3,4,0,4,5,1,4,1,5,6,2,4,6,8,9,7,4,8,10,11,9])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        fff1=allFMTSLeavesPerCommonSupport1[0][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(fff1,mst)
        fff1.loadArraysIfNecessary()
        v=mml2.buildDataArray(fsst,fields,fff1.getUndergroundDataArray())
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        self.assertEqual(fff1.getName(),"f2Cell")
        self.assertNotEqual(v.getHiddenCppPointer(),fff0.getUndergroundDataArray().getHiddenCppPointer())
        vExp=DataArrayDouble([2,3,0,1,4,5])
        self.assertTrue(v.isEqual(vExp,1e-12))
        fff2=allFMTSLeavesPerCommonSupport1[0][0][1][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(fff2,mst)
        fff2.loadArraysIfNecessary()
        v=mml2.buildDataArray(fsst,fields,fff2.getUndergroundDataArray())
        self.assertEqual(fff2.getName(),"f1Node")
        self.assertNotEqual(v.getHiddenCppPointer(),fff0.getUndergroundDataArray().getHiddenCppPointer())
        vExp=DataArrayDouble([11,10,9,8,7,6,5,4,3,2,1,0])
        self.assertTrue(v.isEqual(vExp,1e-12))
        ###
        fcscp=allFMTSLeavesPerCommonSupport1[1][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(ncc)# here all the 16 nodes are taken
        self.assertTrue(a0.isEqual(DataArrayDouble([(0.,0.,0.),(1.,0.,0.),(2.,0.,0.),(3.,0.,0.),(0.,1.,0.),(1.,1.,0.),(2.,1.,0.),(3.,1.,0.),(0.,2.,0.),(1.,2.,0.),(2.,2.,0.),(3.,2.,0.),(0.,3.,0.),(1.,3.,0.),(2.,3.,0.),(3.,3.,0.)]),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([1,1,5,5,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,2,4,8,12,17,22,27])))
        self.assertTrue(a3.isEqual(DataArrayInt([1,8,1,13,3,2,6,3,3,6,7,3,4,0,4,5,1,4,1,5,6,2,4,6,10,11,7,4,10,14,15,11])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        fff3=allFMTSLeavesPerCommonSupport1[1][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(fff3,mst)
        fff3.loadArraysIfNecessary()
        v=mml2.buildDataArray(fsst,fields,fff3.getUndergroundDataArray())
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        self.assertEqual(fff3.getName(),"f3Cell")
        self.assertNotEqual(v.getHiddenCppPointer(),fff0.getUndergroundDataArray().getHiddenCppPointer())
        vExp=DataArrayDouble([-1.1,-3.1,5,4,3,2,1,0])
        self.assertTrue(v.isEqual(vExp,1e-12))
        pass

    def test29(self):
        """ This test focused on HEXA27 cell for which the MED numbering is not equal to the VTK numbering. So here the HEXA27 cell is those in MED file documentation (reference element).
        """
        fname="ForMEDReader29.med"
        coo=DataArrayDouble([[0.,2.,2.],[0.,0.,2.],[2.,0.,2.],[2.,2.,2.],[0.,2.,0.],[0.,0.,0.],[2.,0.,0.],[2.,2.,0.], [0.,1.,2.],[1.,0.,2.],[2.,1.,2.],[1.,2.,2.], [0.,1.,0.],[1.,0.,0.],[2.,1.,0.],[1.,2.,0.], [0.,2.,1.],[0.,0.,1.],[2.,0.,1.],[2.,2.,1.], [1.,1.,2.], [0.,1.,1.],[1.,0.,1.],[2.,1.,1.],[1.,2.,1.], [1.,1.,0.], [1.,1.,1.]])
        m=MEDCouplingUMesh("mesh",3) ; m.setCoords(coo)
        m.allocateCells()
        # MED = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26]
        # VTK = [0,1,2,3,4,5,6,7, 8,9,10,11,12,13,14,15,16,17,18,19,24,22,21,23,20,25,26]
        m.insertNextCell(NORM_HEXA27,[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26])
        fCell=MEDCouplingFieldDouble(ON_CELLS) ; fCell.setName("fCell")
        arrCell=DataArrayDouble([7.]) ; arrCell.setInfoOnComponent(0,"smth") ; fCell.setArray(arrCell)
        fCell.setMesh(m)
        WriteField(fname,fCell,True)
        refCoo=[-1.,-1.,-1.,-1.,1.,-1.,1.,1.,-1.,1.,-1.,-1.,-1.,-1.,1.,-1.,1.,1.,1.,1.,1.,1.,-1.,1.,-1.,0.,-1.,0.,1.,-1.,1.,0.,-1.,0.,-1.,-1.,-1.,0.,1.,0.,1.,1.,1.,0.,1.,0.,-1.,1.,-1.,-1.,0.,-1.,1.,0.,1.,1.,0.,1.,-1.,0.,0.,0.,-1.,-1.,0.,0.,0.,1.,0.,1.,0.,0.,0.,-1.,0.,0.,0.,1.,0.,0.,0.]
        weights=[0.1714677640603571,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.43895747599451346,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.43895747599451346,0.27434842249657115,0.43895747599451346,0.7023319615912209,0.43895747599451346,0.27434842249657115,0.43895747599451346,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.43895747599451346,0.27434842249657115,0.1714677640603571,0.27434842249657115,0.1714677640603571]
        gCoords=[-0.774596669241483,-0.774596669241483,-0.774596669241483,-0.774596669241483,-0.774596669241483,0.0,-0.774596669241483,-0.774596669241483,0.774596669241483,-0.774596669241483,0.0,-0.774596669241483,-0.774596669241483,0.0,0.0,-0.774596669241483,0.0,0.774596669241483,-0.774596669241483,0.774596669241483,-0.774596669241483,-0.774596669241483,0.774596669241483,0.0,-0.774596669241483,0.774596669241483,0.774596669241483,0.0,-0.774596669241483,-0.774596669241483,0.0,-0.774596669241483,0.0,0.0,-0.774596669241483,0.774596669241483,0.0,0.0,-0.774596669241483,0.0,0.0,0.0,0.0,0.0,0.774596669241483,0.0,0.774596669241483,-0.774596669241483,0.0,0.774596669241483,0.0,0.0,0.774596669241483,0.774596669241483,0.774596669241483,-0.774596669241483,-0.774596669241483,0.774596669241483,-0.774596669241483,0.0,0.774596669241483,-0.774596669241483,0.774596669241483,0.774596669241483,0.0,-0.774596669241483,0.774596669241483,0.0,0.0,0.774596669241483,0.0,0.774596669241483,0.774596669241483,0.774596669241483,-0.774596669241483,0.774596669241483,0.774596669241483,0.0,0.774596669241483,0.774596669241483,0.774596669241483]
        fGauss=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fGauss.setName("fGauss")
        fGauss.setMesh(m)
        fGauss.setGaussLocalizationOnType(NORM_HEXA27,refCoo,gCoords,weights)
        arrGauss=DataArrayDouble(fGauss.getNumberOfTuplesExpected()) ; arrGauss.setInfoOnComponent(0,"gaussc") ; arrGauss.iota()
        fGauss.setArray(arrGauss)
        WriteFieldUsingAlreadyWrittenMesh(fname,fGauss)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),2)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),2)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),2)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        self.assertTrue(a0.isEqual(coo,1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([29])))
        self.assertTrue(a2.isEqual(DataArrayInt([0])))
        # the connectivity must be not a iota as declared in m.insertNextCell
        self.assertTrue(a3.isEqual(DataArrayInt([27,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,24,22,21,23,20,25,26])))# the test is on this line to check that connectivity has been processed for HEXA27
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        ffCell=allFMTSLeavesPerCommonSupport1[0][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
        ffCell.loadArraysIfNecessary()
        v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
        self.assertEqual(ffCell.getName(),"fCell")
        self.assertTrue(v.isEqual(arrCell,1e-12)) ; self.assertTrue(v.isEqualWithoutConsideringStr(DataArrayDouble([7.]),1e-12)) ; self.assertEqual(v.getInfoOnComponents(),["smth"])
        del ffCell
        #
        ffGauss=allFMTSLeavesPerCommonSupport1[0][0][1][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(ffGauss,mst)
        ffGauss.loadArraysIfNecessary()
        v=mml2.buildDataArray(fsst,fields,ffGauss.getUndergroundDataArray())
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        self.assertEqual(v.getHiddenCppPointer(),ffGauss.getUndergroundDataArray().getHiddenCppPointer())
        self.assertEqual(ffGauss.getName(),"fGauss")
        self.assertTrue(v.isEqual(arrGauss, 1e-12)) ; self.assertTrue(v.isEqualWithoutConsideringStr(DataArrayDouble(list(range(27))), 1e-12)) ; self.assertEqual(v.getInfoOnComponents(), ["gaussc"])
        ffGauss=allFMTSLeavesPerCommonSupport1[0][0][1][0]
        pass

    def test30(self):
        """ This test is focused on cartesian meshes. Here the cartesian mesh "CartMesh" has a field on HEXA8 (FieldOnCells) and a field on QUAD4 (FieldOnFaces).
        So the first one (FieldOnCells) lies on a cartesian mesh whereas the second one lies on unstructured one.
        """
        fname="ForMEDReader30.med"
        c=MEDCouplingCMesh()
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(4) ; arrY.iota()
        arrZ=DataArrayDouble(5) ; arrZ.iota()
        c.setCoords(arrX,arrY,arrZ)
        c.setName("CartMesh")
        cc=MEDFileCMesh()
        cc.setMesh(c)
        tmpFacesMesh=c.build1SGTSubLevelMesh()
        famIdFaces=DataArrayInt(98) ; famIdFaces[:36]=-1 ; famIdFaces[36:68]=-2 ; famIdFaces[68:]=-3
        famIdCells=DataArrayInt(24) ; famIdCells[:]=0
        #cc.setFamilyFieldArr(0,famIdCells)
        #cc.setFamilyFieldArr(-1,famIdFaces)
        cc.addFamily("FacesX",-1) ; cc.addFamily("FacesY",-2) ; cc.addFamily("FacesZ",-3)
        cc.setFamiliesOnGroup("FacesX1",["FacesX"])
        cc.setFamiliesOnGroup("FacesY1",["FacesY"])
        cc.setFamiliesOnGroup("FacesZ1",["FacesZ"])
        #
        fmts0=MEDFileFieldMultiTS()
        fmts1=MEDFileFieldMultiTS()
        for i in range(30):
            f1ts=MEDFileField1TS()
            fFaces=MEDCouplingFieldDouble(ON_CELLS) ; fFaces.setName("FieldOnFaces")
            arr=DataArrayDouble(98) ; arr.iota() ; arr[i]=100.
            fFaces.setArray(arr)
            fFaces.setTime(float(i)+0.1,i,-1)
            fFaces.setMesh(tmpFacesMesh)
            f1ts.setFieldNoProfileSBT(fFaces)
            fmts0.pushBackTimeStep(f1ts)
            #
            f1ts=MEDFileField1TS()
            fCells=MEDCouplingFieldDouble(ON_CELLS) ; fCells.setName("FieldOnCells")
            arr=DataArrayDouble(24) ; arr.iota() ; arr[i%24]=30.
            fCells.setArray(arr)
            fCells.setTime(float(i)+0.1,i,-1)
            fCells.setMesh(c)
            f1ts.setFieldNoProfileSBT(fCells)
            fmts1.pushBackTimeStep(f1ts)
            pass
        fs=MEDFileFields()
        fs.pushField(fmts0)
        fs.pushField(fmts1)
        cc.write(fname,2)
        fs.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),2)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),2)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[1][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDCMeshMultiLev)) # here CMesh is important
        (a,b,c),d=mml2.buildVTUArrays()
        self.assertTrue(d)#d is True because the a,b and c are directly those in the internal data structure
        self.assertTrue(a.isEqual(arrX,1e-12))
        self.assertTrue(b.isEqual(arrY,1e-12))
        self.assertTrue(c.isEqual(arrZ,1e-12))
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(30):
            ffCell=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            myarr=DataArrayDouble(24) ; myarr.iota() ; myarr[i%24]=30.
            self.assertEqual(ffCell.getName(),"FieldOnCells")
            self.assertTrue(v.isEqual(myarr,1e-12))
            pass
        #
        fcscp=allFMTSLeavesPerCommonSupport1[1][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev)) # here UMesh is important
        ref=ms[0].getImplicitFaceMesh().getCoords().getHiddenCppPointer()
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertEqual(ref,a0.getHiddenCppPointer())
        self.assertTrue(ncc)
        self.assertTrue(a1.isEqual(DataArrayByte([9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310,315,320,325,330,335,340,345,350,355,360,365,370,375,380,385,390,395,400,405,410,415,420,425,430,435,440,445,450,455,460,465,470,475,480,485])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,0,12,15,3,4,12,24,27,15,4,24,36,39,27,4,36,48,51,39,4,3,15,18,6,4,15,27,30,18,4,27,39,42,30,4,39,51,54,42,4,6,18,21,9,4,18,30,33,21,4,30,42,45,33,4,42,54,57,45,4,1,13,16,4,4,13,25,28,16,4,25,37,40,28,4,37,49,52,40,4,4,16,19,7,4,16,28,31,19,4,28,40,43,31,4,40,52,55,43,4,7,19,22,10,4,19,31,34,22,4,31,43,46,34,4,43,55,58,46,4,2,14,17,5,4,14,26,29,17,4,26,38,41,29,4,38,50,53,41,4,5,17,20,8,4,17,29,32,20,4,29,41,44,32,4,41,53,56,44,4,8,20,23,11,4,20,32,35,23,4,32,44,47,35,4,44,56,59,47,4,0,12,13,1,4,12,24,25,13,4,24,36,37,25,4,36,48,49,37,4,1,13,14,2,4,13,25,26,14,4,25,37,38,26,4,37,49,50,38,4,3,15,16,4,4,15,27,28,16,4,27,39,40,28,4,39,51,52,40,4,4,16,17,5,4,16,28,29,17,4,28,40,41,29,4,40,52,53,41,4,6,18,19,7,4,18,30,31,19,4,30,42,43,31,4,42,54,55,43,4,7,19,20,8,4,19,31,32,20,4,31,43,44,32,4,43,55,56,44,4,9,21,22,10,4,21,33,34,22,4,33,45,46,34,4,45,57,58,46,4,10,22,23,11,4,22,34,35,23,4,34,46,47,35,4,46,58,59,47,4,0,1,4,3,4,3,4,7,6,4,6,7,10,9,4,1,2,5,4,4,4,5,8,7,4,7,8,11,10,4,12,13,16,15,4,15,16,19,18,4,18,19,22,21,4,13,14,17,16,4,16,17,20,19,4,19,20,23,22,4,24,25,28,27,4,27,28,31,30,4,30,31,34,33,4,25,26,29,28,4,28,29,32,31,4,31,32,35,34,4,36,37,40,39,4,39,40,43,42,4,42,43,46,45,4,37,38,41,40,4,40,41,44,43,4,43,44,47,46,4,48,49,52,51,4,51,52,55,54,4,54,55,58,57,4,49,50,53,52,4,52,53,56,55,4,55,56,59,58])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        for i in range(30):
            ffCell=allFMTSLeavesPerCommonSupport1[1][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            myarr=DataArrayDouble(98) ; myarr.iota() ; myarr[i]=100.
            self.assertEqual(ffCell.getName(),"FieldOnFaces")
            self.assertTrue(v.isEqual(myarr,1e-12))
            pass
        pass

    def test31(self):
        """non regression test of EDF 7972"""
        fname="ForMEDReader31.med"
        c=MEDCouplingCMesh()
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(4) ; arrY.iota()
        arrZ=DataArrayDouble(5) ; arrZ.iota()
        c.setCoords(arrX,arrY,arrZ)
        c.setName("CartMesh")
        cc=MEDFileCMesh()
        cc.setMesh(c)
        famIdCells=DataArrayInt(24) ; famIdCells[:]=0
        cc.setFamilyFieldArr(0,famIdCells)
        #cc.setFamilyFieldArr(-1,famIdFaces)
        cc.addFamily("FacesX",-1) ; cc.addFamily("FacesY",-2) ; cc.addFamily("FacesZ",-3)
        cc.setFamiliesOnGroup("FacesX1",["FacesX"])
        cc.setFamiliesOnGroup("FacesY1",["FacesY"])
        cc.setFamiliesOnGroup("FacesZ1",["FacesZ"])
        fmts0=MEDFileFieldMultiTS()
        fmts1=MEDFileFieldMultiTS()
        pfl=DataArrayInt(11) ; pfl.iota() ; pfl.setName("PflOnHECA8")
        for i in range(30):
            f1ts=MEDFileField1TS()
            fFaces=MEDCouplingFieldDouble(ON_CELLS) ; fFaces.setName("FieldOnCells")
            arr=DataArrayDouble(11) ; arr.iota() ; arr[i%11]=100.
            fFaces.setArray(arr)
            fFaces.setTime(float(i)+0.1,i,-1)
            fFaces.setMesh(c.buildUnstructured()[:11])
            f1ts.setFieldProfile(fFaces,cc,0,pfl)# here, a test is done to check that "NORM_HEXA8" string is not 30 times appended at the end of pfl name.
            self.assertEqual("PflOnHECA8",pfl.getName())
            fmts0.pushBackTimeStep(f1ts)
            pass
        fs=MEDFileFields()
        fs.pushField(fmts0)
        cc.write(fname,2)
        fs.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev)) # here UMesh is important
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        self.assertTrue(not ncc)# here ncc=False because the coordinates are not in ms neither in children. This is the most important line in the test.
        self.assertTrue(a0.isEqual(DataArrayDouble([0.,0.,0.,1.,0.,0.,2.,0.,0.,0.,1.,0.,1.,1.,0.,2.,1.,0.,0.,2.,0.,1.,2.,0.,2.,2.,0.,0.,3.,0.,1.,3.,0.,2.,3.,0.,0.,0.,1.,1.,0.,1.,2.,0.,1.,0.,1.,1.,1.,1.,1.,2.,1.,1.,0.,2.,1.,1.,2.,1.,2.,2.,1.,0.,3.,1.,1.,3.,1.,2.,3.,1.,0.,0.,2.,1.,0.,2.,2.,0.,2.,0.,1.,2.,1.,1.,2.,2.,1.,2.,0.,2.,2.,1.,2.,2.,2.,2.,2.,0.,3.,2.,1.,3.,2.,2.,3.,2.,0.,0.,3.,1.,0.,3.,2.,0.,3.,0.,1.,3.,1.,1.,3.,2.,1.,3.,0.,2.,3.,1.,2.,3.,2.,2.,3.,0.,3.,3.,1.,3.,3.,2.,3.,3.,0.,0.,4.,1.,0.,4.,2.,0.,4.,0.,1.,4.,1.,1.,4.,2.,1.,4.,0.,2.,4.,1.,2.,4.,2.,2.,4.,0.,3.,4.,1.,3.,4.,2.,3.,4.],60,3),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([12,12,12,12,12,12,12,12,12,12,12])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,9,18,27,36,45,54,63,72,81,90])))
        self.assertTrue(a3.isEqual(DataArrayInt([8,1,0,3,4,13,12,15,16,8,2,1,4,5,14,13,16,17,8,4,3,6,7,16,15,18,19,8,5,4,7,8,17,16,19,20,8,7,6,9,10,19,18,21,22,8,8,7,10,11,20,19,22,23,8,13,12,15,16,25,24,27,28,8,14,13,16,17,26,25,28,29,8,16,15,18,19,28,27,30,31,8,17,16,19,20,29,28,31,32,8,19,18,21,22,31,30,33,34])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        for i in range(30):
            ffCell=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            # self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer()) # to be improved... maybe this line could be true
            myarr=DataArrayDouble(11) ; myarr.iota() ; myarr[i%11]=100.
            self.assertEqual(ffCell.getName(),"FieldOnCells")
            self.assertTrue(v.isEqual(myarr,1e-12))
            pass
        pass

    def test32(self):
        """ This test is close to test30 except that here the profiles on dim-1 of structured mesh is considered here."""
        fname="ForMEDReader32.med"
        c=MEDCouplingCMesh()
        arrX=DataArrayDouble(3) ; arrX.iota()
        arrY=DataArrayDouble(4) ; arrY.iota()
        arrZ=DataArrayDouble(5) ; arrZ.iota()
        c.setCoords(arrX,arrY,arrZ)
        c.setName("CartMesh")
        cc=MEDFileCMesh()
        cc.setMesh(c)
        tmpFacesMesh=c.build1SGTSubLevelMesh()
        famIdFaces=DataArrayInt(98) ; famIdFaces[:36]=-1 ; famIdFaces[36:68]=-2 ; famIdFaces[68:]=-3
        famIdCells=DataArrayInt(24) ; famIdCells[:]=0
        cc.setFamilyFieldArr(0,famIdCells)
        #cc.setFamilyFieldArr(-1,famIdFaces)
        cc.addFamily("FacesX",-1) ; cc.addFamily("FacesY",-2) ; cc.addFamily("FacesZ",-3)
        cc.setFamiliesOnGroup("FacesX1",["FacesX"])
        cc.setFamiliesOnGroup("FacesY1",["FacesY"])
        cc.setFamiliesOnGroup("FacesZ1",["FacesZ"])
        fmts0=MEDFileFieldMultiTS()
        fmts1=MEDFileFieldMultiTS()
        pfl=DataArrayInt(31) ; pfl.iota() ; pfl.setName("PflOnQUAD4")
        for i in range(30):
            f1ts=MEDFileField1TS()
            fFaces=MEDCouplingFieldDouble(ON_CELLS) ; fFaces.setName("FieldOnFaces")
            arr=DataArrayDouble(31) ; arr.iota() ; arr[i]=100.
            fFaces.setArray(arr)
            fFaces.setTime(float(i)+0.1,i,-1)
            fFaces.setMesh(tmpFacesMesh[:31])
            f1ts.setFieldProfile(fFaces,cc,-1,pfl)# here, a test is done to check that "NORM_QUAD4" string is not 30 times appended at the end of pfl name.
            self.assertEqual("PflOnQUAD4",pfl.getName())
            fmts0.pushBackTimeStep(f1ts)
            #
            f1ts=MEDFileField1TS()
            fCells=MEDCouplingFieldDouble(ON_CELLS) ; fCells.setName("FieldOnCells")
            arr=DataArrayDouble(24) ; arr.iota() ; arr[i%24]=30.
            fCells.setArray(arr)
            fCells.setTime(float(i)+0.1,i,-1)
            fCells.setMesh(c)
            f1ts.setFieldNoProfileSBT(fCells)
            fmts1.pushBackTimeStep(f1ts)
            pass
        fs=MEDFileFields()
        fs.pushField(fmts0)
        fs.pushField(fmts1)
        cc.write(fname,2)
        fs.write(fname,0)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),2)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),2)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),1)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[1][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDCMeshMultiLev)) # here CMesh is important
        (a,b,c),d=mml2.buildVTUArrays()
        self.assertTrue(d)#d is True because the a,b and c are directly those in the internal data structure
        self.assertTrue(a.isEqual(arrX,1e-12))
        self.assertTrue(b.isEqual(arrY,1e-12))
        self.assertTrue(c.isEqual(arrZ,1e-12))
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(30):
            ffCell=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            myarr=DataArrayDouble(24) ; myarr.iota() ; myarr[i%24]=30.
            self.assertEqual(ffCell.getName(),"FieldOnCells")
            self.assertTrue(v.isEqual(myarr,1e-12))
            pass
        #
        fcscp=allFMTSLeavesPerCommonSupport1[1][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev)) # here UMesh is important
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        self.assertTrue(ncc)# True because, the coords are computed by the implicit unstructured level -1 structured mesh
        self.assertTrue(a0.isEqual(DataArrayDouble([0.,0.,0.,1.,0.,0.,2.,0.,0.,0.,1.,0.,1.,1.,0.,2.,1.,0.,0.,2.,0.,1.,2.,0.,2.,2.,0.,0.,3.,0.,1.,3.,0.,2.,3.,0.,0.,0.,1.,1.,0.,1.,2.,0.,1.,0.,1.,1.,1.,1.,1.,2.,1.,1.,0.,2.,1.,1.,2.,1.,2.,2.,1.,0.,3.,1.,1.,3.,1.,2.,3.,1.,0.,0.,2.,1.,0.,2.,2.,0.,2.,0.,1.,2.,1.,1.,2.,2.,1.,2.,0.,2.,2.,1.,2.,2.,2.,2.,2.,0.,3.,2.,1.,3.,2.,2.,3.,2.,0.,0.,3.,1.,0.,3.,2.,0.,3.,0.,1.,3.,1.,1.,3.,2.,1.,3.,0.,2.,3.,1.,2.,3.,2.,2.,3.,0.,3.,3.,1.,3.,3.,2.,3.,3.,0.,0.,4.,1.,0.,4.,2.,0.,4.,0.,1.,4.,1.,1.,4.,2.,1.,4.,0.,2.,4.,1.,2.,4.,2.,2.,4.,0.,3.,4.,1.,3.,4.,2.,3.,4.],60,3),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,0,12,15,3,4,12,24,27,15,4,24,36,39,27,4,36,48,51,39,4,3,15,18,6,4,15,27,30,18,4,27,39,42,30,4,39,51,54,42,4,6,18,21,9,4,18,30,33,21,4,30,42,45,33,4,42,54,57,45,4,1,13,16,4,4,13,25,28,16,4,25,37,40,28,4,37,49,52,40,4,4,16,19,7,4,16,28,31,19,4,28,40,43,31,4,40,52,55,43,4,7,19,22,10,4,19,31,34,22,4,31,43,46,34,4,43,55,58,46,4,2,14,17,5,4,14,26,29,17,4,26,38,41,29,4,38,50,53,41,4,5,17,20,8,4,17,29,32,20,4,29,41,44,32])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        a6,a7=mml2.retrieveFamilyIdsOnCells()
        self.assertTrue(a6 is None)
        self.assertTrue(a7)
        for i in range(30):
            ffCell=allFMTSLeavesPerCommonSupport1[1][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            myarr=DataArrayDouble(31) ; myarr.iota() ; myarr[i]=100.
            self.assertEqual(ffCell.getName(),"FieldOnFaces")
            self.assertTrue(v.isEqual(myarr,1e-12))
            pass
        pass
    
    def test33(self):
        """Non regression test concerning polygons. Thanks Adrien. This bug can't be shown by simply reading an displaying a MED file containing only polygons. A filter must be applied on it to show it. The a2 array was responsible of that bug."""
        fname="ForMEDReader33.med"
        fieldName="ACellField"
        coo=DataArrayDouble([(5.5,0.5),(5.5,-0.5),(6.5,0.5),(6.5,-0.5),(6.5,1.5),(7.5,0.5),(7.5,-0.5),(7.5,1.5),(7.5,2.5),(8.5,0.5),(8.5,-0.5),(8.5,1.5),(8.5,2.5),(8.5,3.5),(8.55,0.5),(8.55,-0.5),(8.55,1.5),(8.55,2.5),(8.55,3.5)])
        m=MEDCouplingUMesh("mesh",2)
        m.setCoords(coo)
        m.allocateCells()
        for i,c in enumerate([(1,0,2,3),(3,2,5,6),(2,4,7,5),(6,5,9,10),(5,7,11,9),(7,8,12,11),(10,9,14,15),(9,11,16,14),(11,12,17,16),(12,13,18,17)]):
            if i<6:
                typ=NORM_QUAD4
                pass
            else:
                typ=NORM_POLYGON
                pass
            m.insertNextCell(typ,c)
            pass
        mm=MEDFileUMesh()
        mm.setMeshAtLevel(0,m)
        mm.write(fname,2)
        for i in range(15):
            fCell0=MEDCouplingFieldDouble(ON_CELLS) ; fCell0.setTime(float(i)+0.1,i,0)
            fCell0.setName(fieldName) ; fCell0.setMesh(m)
            arr=DataArrayDouble(m.getNumberOfCells()) ; arr.iota(0) ; arr[i%10]=100.
            fCell0.setArray(arr) ; arr.setInfoOnComponents(["Comp1 [m]"])
            fCell0.checkConsistencyLight()
            WriteFieldUsingAlreadyWrittenMesh(fname,fCell0)
            pass
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
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
        self.assertTrue(not ncc)# false beacause 2D in MED file
        self.assertTrue(a0.isEqual(DataArrayDouble([(5.5,0.5,0),(5.5,-0.5,0),(6.5,0.5,0),(6.5,-0.5,0),(6.5,1.5,0),(7.5,0.5,0),(7.5,-0.5,0),(7.5,1.5,0),(7.5,2.5,0),(8.5,0.5,0),(8.5,-0.5,0),(8.5,1.5,0),(8.5,2.5,0),(8.5,3.5,0),(8.55,0.5,0),(8.55,-0.5,0),(8.55,1.5,0),(8.55,2.5,0),(8.55,3.5,0)]),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([9,9,9,9,9,9,7,7,7,7])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5,10,15,20,25,30,35,40,45])))# the bug was here.
        self.assertTrue(a3.isEqual(DataArrayInt([4,1,0,2,3,4,3,2,5,6,4,2,4,7,5,4,6,5,9,10,4,5,7,11,9,4,7,8,12,11,4,10,9,14,15,4,9,11,16,14,4,11,12,17,16,4,12,13,18,17])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(15):
            ffCell=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            myarr=DataArrayDouble(10) ; myarr.iota() ; myarr[i%10]=100. ; myarr.setInfoOnComponent(0,"Comp1 [m]")
            self.assertEqual(ffCell.getName(),fieldName)
            self.assertTrue(v.isEqual(myarr,1e-12))
            pass
        pass

    def test34(self):
        """ This test is the thirs ultimate test (base on test12) for the profiles with gauss points.
        This test highlight the hidden imp linked to bug #8655.
        This test is close to test11 but here a 2nd field on cells without profile. So here the mesh is expected to be the same than m.
        """
        fname="ForMEDReader34.med"
        m=MEDCouplingCMesh("mesh")
        arr=DataArrayDouble(5) ; arr.iota()
        m.setCoords(arr,arr)
        m=m.buildUnstructured() ; m.getCoords().setInfoOnComponents(["XX [m]","YYY [km]"])
        mm=MEDFileUMesh() ; mm.setMeshes([m])
        #
        fieldName0="zeField0"
        fieldName1="zeField1"
        fs0=MEDFileFieldMultiTS() ; fs1=MEDFileFieldMultiTS()
        for i in range(5):
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName0) ; fNode.setMesh(m)
            fNode.setGaussLocalizationOnCells([0,2,3,4,7,15],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7],[0.8,0.2])
            fNode.setGaussLocalizationOnCells([1,5,8,9],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3],[0.8,0.05,0.1,0.04,0.01])
            fNode.setGaussLocalizationOnCells([6,10,13],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2],[0.8,0.05,0.1,0.04])
            fNode.setGaussLocalizationOnCells([11,12,14],[-1.,-1.,1.,-1.,1.,1.,-1.,1.],[0.5,0.5,0.7,0.7,0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.8,0.8],[0.8,0.05,0.1,0.01,0.02,0.005,0.005])
            arr=DataArrayDouble(2*(2*6+5*4+4*3+7*3)) ; arr.iota(0+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_0 [m]","Com2_0 [s^2]"]) ; fNode.checkConsistencyLight()
            f.setFieldNoProfileSBT(fNode)
            fs0.pushBackTimeStep(f)
            #
            f=MEDFileField1TS()
            fNode=MEDCouplingFieldDouble(ON_CELLS) ; fNode.setTime(float(i),i,0)
            fNode.setName(fieldName1) ; fNode.setMesh(m)
            arr=DataArrayDouble(2*16) ; arr.iota(300+1000*i) ; arr.rearrange(2)
            fNode.setArray(arr) ; arr.setInfoOnComponents(["Comp1_1 [m]","Com2_1 [s^2]"]) ; fNode.checkConsistencyLight()
            f.setFieldNoProfileSBT(fNode)
            fs1.pushBackTimeStep(f)
            pass
        mm.write(fname,2)
        fs0.write(fname,0) ; fs1.write(fname,0)
        a0Exp=mm.getCoords().deepCopy()
        del m,mm,fs0,fs1,f,fNode
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    if itmp.presenceOfMultiDiscPerGeoType():
                        tmp2=itmp.splitMultiDiscrPerGeoTypes()
                        for iii,itmp2 in enumerate(tmp2):
                            name="%s_%i"%(itmp2.getName(),iii)
                            itmp2.setName(name)
                            allFMTSLeavesToDisplay2.append(itmp2)
                            pass
                        pass
                    else:
                        allFMTSLeavesToDisplay2.append(itmp)
                        pass
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        # Here 2 MED fields in input and at the end 5 ! 1+4 ! 4 fields have been built from zeField0 due to subspliting per dis / per geo type
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),5)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 2 fields are defined on the same time steps
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),5)
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport),5)
        for i in range(5):
            self.assertEqual(len(allFMTSLeavesPerCommonSupport[i][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        #
        fcscp=allFMTSLeavesPerCommonSupport[4][1]
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
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1, 5):
            self.assertTrue((fcscp.isDataSetSupportEqualToThePreviousOne(i,fields)))
            pass
        for i in range(5):
            f=allFMTSLeavesPerCommonSupport[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            vExp=DataArrayDouble([(0.,1.),(2.,3.),(14.,15.),(16.,17.),(18.,19.),(20.,21.),(22.,23.),(24.,25.),(44.,45.),(46.,47.),(126.,127.),(128.,129.)])
            vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
            vExp+=i*1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[1][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            vExp=DataArrayDouble([(4.,5.),(6.,7.),(8.,9.),(10.,11.),(12.,13.),(26.,27.),(28.,29.),(30.,31.),(32.,33.),(34.,35.),(48.,49.),(50.,51.),(52.,53.),(54.,55.),(56.,57.),(58.,59.),(60.,61.),(62.,63.),(64.,65.),(66.,67.)])
            vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
            vExp+=i*1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[2][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            vExp=DataArrayDouble([(36.,37.),(38.,39.),(40.,41.),(42.,43.),(68.,69.),(70.,71.),(72.,73.),(74.,75.),(104.,105.),(106.,107.),(108.,109.),(110.,111.)])
            vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
            vExp+=i*1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[3][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName0)
            vExp=DataArrayDouble([(76,77),(78,79),(80,81),(82,83),(84,85),(86,87),(88,89),(90,91),(92,93),(94,95),(96,97),(98,99),(100,101),(102,103),(112,113),(114,115),(116,117),(118,119),(120,121),(122,123),(124,125)])
            vExp.setInfoOnComponents(['Comp1_0 [m]','Com2_0 [s^2]'])
            vExp+=i*1000
            self.assertTrue(v.isEqual(vExp,1e-12))
            #
            f=allFMTSLeavesPerCommonSupport[4][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(f,mst)
            f.loadArraysIfNecessary()
            v=mml.buildDataArray(fsst,fields,f.getUndergroundDataArray())
            self.assertEqual(f.getName(),fieldName1)
            self.assertEqual(v.getHiddenCppPointer(),f.getUndergroundDataArray().getHiddenCppPointer())
            vExp=DataArrayDouble(16*2) ; vExp.iota(300+i*1000) ; vExp.rearrange(2) ; vExp.setInfoOnComponents(['Comp1_1 [m]','Com2_1 [s^2]'])
            self.assertTrue(v.isEqual(vExp,1e-12))
            pass
        pass

    def test35(self):
        """ Emulate MEDReader in // mode context. Here a Simple mesh having more nodes than really needed. This test focuses on that point particulary."""
        fname="ForMEDReader35.med"
        arrX=DataArrayDouble(7) ; arrX.iota()
        arrY=DataArrayDouble([0.,1.])
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY) ; m=m.buildUnstructured() ; m=m[[0,5,1,4,2,3]] ; m.changeSpaceDimension(3,0.) ; m.setName("Mesh")
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(m) ; f.setName("Field") ; f.setArray(DataArrayDouble([(0.1,1.1),(2.1,3.1),(4.1,5.1),(6.1,7.1),(8.1,9.1),(10.1,11.1)])) ; f.getArray().setInfoOnComponents(["aa","bbb"])
        WriteUMesh(fname,m,True)
        WriteFieldUsingAlreadyWrittenMesh(fname,f)
        ########## GO for reading in MEDReader,by not loading all. Mesh is fully loaded but not fields values
        ms=MEDFileMeshes() # here we reproduce what is done by ParaMEDFileMeshes.ParaNew
        ms.pushMesh(MEDFileUMesh.LoadPartOf(fname,"Mesh",[NORM_QUAD4],[0,2,1],-1,-1));
        ms[0].zipCoords()
        ms.cartesianizeMe()
        #
        fields=MEDFileFields.LoadPartOf(fname,False,ms);
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    if itmp.presenceOfMultiDiscPerGeoType():
                        tmp2=itmp.splitMultiDiscrPerGeoTypes()
                        for iii,itmp2 in enumerate(tmp2):
                            name="%s_%i"%(itmp2.getName(),iii)
                            itmp2.setName(name)
                            allFMTSLeavesToDisplay2.append(itmp2)
                            pass
                        pass
                    else:
                        allFMTSLeavesToDisplay2.append(itmp)
                        pass
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        #
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),1)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1) # one time serie here : because the 2 fields are defined on the same time steps
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
        self.assertTrue(ncc)
        self.assertTrue(a0.isEqual(m.getCoords()[[0,1,5,6,7,8,12,13]],1e-12))# <- the aim of the test
        self.assertTrue(a1.isEqual(DataArrayByte([9,9])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,5])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,1,0,4,5,4,3,2,6,7])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny().isEqual(DataArrayInt([0,1,5,6,7,8,12,13])))
        f2=allFMTSLeavesPerCommonSupport[0][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(f2,mst)
        f2.loadArraysIfNecessary()
        v=mml.buildDataArray(fsst,fields,f2.getUndergroundDataArray())
        self.assertEqual(f2.getName(),f.getName())
        vExp=DataArrayDouble([(0.1,1.1),(2.1,3.1)])
        vExp.setInfoOnComponents(['aa','bbb'])
        self.assertTrue(v.isEqual(vExp,1e-12))
        pass

    def test36(self):
        """Bug EDF11027. Here mesh at level 0 (TRI3) does not fetch all the nodes. Level -1 (SEG2) does not fetch all the nodes neither. But all TRI3 + all SEG2 fetch all nodes.
        aaa field on GAUSSPoints lying only on TRI3 share the same support than profile node field ccc.
        But bbb field on all nodes is not on the same support. Past optimization that make the assumtion a support on all lev0 cells lies on all nodes is now over."""
        meshName="mesh"
        fname="ForMEDReader36.med"
        c=DataArrayDouble([(0,0),(1,0),(1,1),(0,1),(2,0),(-1,0),(1,2)])
        m0=MEDCoupling1SGTUMesh(meshName,NORM_TRI3)
        m0.setCoords(c)
        m0.setNodalConnectivity(DataArrayInt([0,2,1,3,2,0,2,4,1]))
        mm=MEDFileUMesh()
        mm[0]=m0
        m1=MEDCoupling1SGTUMesh(meshName,NORM_SEG2)
        m1.setCoords(c)
        m1.setNodalConnectivity(DataArrayInt([5,0,0,3,3,2,2,6]))
        mm[-1]=m1
        #
        zeTime=(1.1,2,3)
        ff1=MEDFileField1TS()
        f1=MEDCouplingFieldDouble(ON_NODES) ; f1.setMesh(m0)
        arr=DataArrayDouble(7) ; arr.iota(2000)
        f1.setArray(arr)
        f1.setName("bbb")
        f1.checkConsistencyLight()
        f1.setTime(*zeTime)
        ff1.setFieldNoProfileSBT(f1)
        #
        ff2=MEDFileField1TS()
        f2=MEDCouplingFieldDouble(ON_GAUSS_NE) ; f2.setMesh(m0)
        arr=DataArrayDouble(9) ; arr.iota(4000)
        f2.setArray(arr)
        f2.setName("ddd")
        f2.checkConsistencyLight()
        f2.setTime(*zeTime)
        ff2.setFieldNoProfileSBT(f2)
        #
        ff3=MEDFileField1TS()
        f3=MEDCouplingFieldDouble(ON_GAUSS_PT) ; f3.setMesh(m0)
        f3.setGaussLocalizationOnType(NORM_TRI3,[0,0,1,0,0,1],[0.333333,0.333333],[0.5])
        arr=DataArrayDouble(3) ; arr.iota(1000)
        f3.setArray(arr)
        f3.checkConsistencyLight()
        f3.setTime(*zeTime)
        f3.setName("aaa")
        ff3.setFieldNoProfileSBT(f3)
        #
        ff4=MEDFileField1TS()
        m0d=m0.deepCopy() ; m0d.zipCoords()
        f4=MEDCouplingFieldDouble(ON_NODES) ; f4.setMesh(m0d)
        arr=DataArrayDouble(5) ; arr.iota(3000)
        f4.setArray(arr)
        f4.setName("ccc")
        f4.checkConsistencyLight()
        f4.setTime(*zeTime)
        pfl=DataArrayInt([0,1,2,3,4]) ; pfl.setName("PFL")
        ff4.setFieldProfile(f4,mm,0,pfl)
        #
        mm.write(fname,2)
        ff3.write(fname,0)
        ff1.write(fname,0)
        ff4.write(fname,0)
        ###
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),3)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        self.assertEqual(len(allFMTSLeavesPerTimeSeries[0]),3)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[0][0]),2)
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1[1][0]),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)# here ncc=False because the coordinates are not in ms neither in children.
        self.assertTrue(a0.isEqual(DataArrayDouble([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(2,0,0)]),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,5,5])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,8])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,0,2,1,3,3,2,0,3,2,4,1])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1):
            ffCell=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            v.isEqual(DataArrayDouble([1000,1001,1002]),1e-12)
            #
            ffCell=allFMTSLeavesPerCommonSupport1[0][0][1][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            v.isEqual(DataArrayDouble([3000,3001,3002,3003,3004]),1e-12)
            pass
        fcscp=allFMTSLeavesPerCommonSupport1[1][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)# here ncc=False because the coordinates are not in ms neither in children.
        self.assertTrue(a0.isEqual(DataArrayDouble([(0,0,0),(1,0,0),(1,1,0),(0,1,0),(2,0,0),(-1,0,0),(1,2,0)]),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([5,5,5,3,3,3,3])))
        self.assertTrue(a2.isEqual(DataArrayInt([0,4,8,12,15,18,21])))
        self.assertTrue(a3.isEqual(DataArrayInt([3,0,2,1,3,3,2,0,3,2,4,1,2,5,0,2,0,3,2,3,2,2,2,6])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1):
            ffCell=allFMTSLeavesPerCommonSupport1[1][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            v.isEqual(DataArrayDouble([2000,2001,2002,2003,2004,2005,2006]),1e-12)
            pass
        pass

    def test37(self):
        """ Introduction of non cartesian meshes management. Here cylindrical."""
        fname="ForMEDReader37.med"
        meshName="mesh"
        description="Cylindrical grid"
        comps=["X [cm]","Y [cm]","Z [cm]"]
        arrX=DataArrayDouble(3) ; arrX.iota() ; arrX*=0.8 ; arrX.setInfoOnComponent(0,comps[0])
        arrY=DataArrayDouble(4) ; arrY.iota() ; arrY*=pi/(len(arrY)-1) ; arrY.setInfoOnComponent(0,comps[1])
        arrZ=DataArrayDouble(5) ; arrZ.iota() ; arrZ*=1.6 ; arrZ-=8. ; arrZ.setInfoOnComponent(0,comps[2])
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY,arrZ) ; m.setName(meshName)
        mm=MEDFileCMesh() ; mm.setMesh(m) ; mm.setDescription(description)
        mm.setAxisType(AX_CYL) # the test is here !
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(m) ; f.setName("Field")
        arr=DataArrayDouble(m.getNumberOfCells()) ; arr.iota() ; arr*=0.1 ; f.setArray(arr) ; f.checkConsistencyLight()
        ff=MEDFileField1TS() ; ff.setFieldNoProfileSBT(f)
        fmts=MEDFileFieldMultiTS() ; fmts.pushBackTimeStep(ff)
        #
        ms=MEDFileMeshes() ; ms.pushMesh(mm)
        fields=MEDFileFields() ; fields.pushField(fmts)
        ms.write(fname,2) ; fields.write(fname,0)
        #
        del mm,fmts,fields,ms
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        ms.cartesianizeMe()
        #
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        #
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),1)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDCurveLinearMeshMultiLev))# <- hehe it is a CurveLinear no more a CMesh !
        a,b,c=mml2.buildVTUArrays()
        self.assertTrue(c)# the array is thoose in structure
        ref_a=DataArrayDouble([0.,0.,-8.,0.8,0.,-8.,1.6,0.,-8.,0.,0.,-8.,0.4,0.6928203230275509,-8.,0.8,1.3856406460551018,-8.,-0.,0.,-8.,-0.4,0.692820323027551,-8.,-0.8,1.385640646055102,-8.,-0.,0.,-8.,-0.8,0.,-8.,-1.6,0.,-8.,0.,0.,-6.4,0.8,0.,-6.4,1.6,0.,-6.4,0.,0.,-6.4,0.4,0.6928203230275509,-6.4,0.8,1.3856406460551018,-6.4,-0.,0.,-6.4,-0.4,0.692820323027551,-6.4,-0.8,1.385640646055102,-6.4,-0.,0.,-6.4,-0.8,0.,-6.4,-1.6,0.,-6.4,0.,0.,-4.8,0.8,0.,-4.8,1.6,0.,-4.8,0.,0.,-4.8,0.4,0.6928203230275509,-4.8,0.8,1.3856406460551018,-4.8,-0.,0.,-4.8,-0.4,0.692820323027551,-4.8,-0.8,1.385640646055102,-4.8,-0.,0.,-4.8,-0.8,0.,-4.8,-1.6,0.,-4.8,0.,0.,-3.2,0.8,0.,-3.2,1.6,0.,-3.2,0.,0.,-3.2,0.4,0.6928203230275509,-3.2,0.8,1.3856406460551018,-3.2,-0.,0.,-3.2,-0.4,0.692820323027551,-3.2,-0.8,1.385640646055102,-3.2,-0.,0.,-3.2,-0.8,0.,-3.2,-1.6,0.,-3.2,0.,0.,-1.6,0.8,0.,-1.6,1.6,0.,-1.6,0.,0.,-1.6,0.4,0.6928203230275509,-1.6,0.8,1.3856406460551018,-1.6,-0.,0.,-1.6,-0.4,0.692820323027551,-1.6,-0.8,1.385640646055102,-1.6,-0.,0.,-1.6,-0.8,0.,-1.6,-1.6,0.,-1.6],60,3)
        ref_a.setInfoOnComponents(comps)
        self.assertTrue(a.isEqual(ref_a,1e-14))
        self.assertEqual(b,[3,4,5])
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1):
            ffCell=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            self.assertTrue(v.isEqual(DataArrayDouble([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3]),1e-14))
            pass
        pass

    def test38(self):
        """ Introduction of non cartesian meshes management. Here spherical."""
        fname="ForMEDReader38.med"
        meshName="mesh"
        description="Spherical grid"
        comps=["X [cm]","Y [cm]","Z [cm]"]
        arrX=DataArrayDouble(3) ; arrX.iota() ; arrX*=0.8 ; arrX.setInfoOnComponent(0,comps[0])
        arrY=DataArrayDouble(4) ; arrY.iota() ; arrY*=pi/(len(arrY)-1) ; arrY.setInfoOnComponent(0,comps[1])
        arrZ=DataArrayDouble(5) ; arrZ.iota() ; arrZ*=2*pi/(len(arrZ)-1) ; arrZ.setInfoOnComponent(0,comps[2])
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY,arrZ) ; m.setName(meshName)
        mm=MEDFileCMesh() ; mm.setMesh(m) ; mm.setDescription(description)
        mm.setAxisType(AX_SPHER) # the test is here !
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(m) ; f.setName("Field")
        arr=DataArrayDouble(m.getNumberOfCells()) ; arr.iota() ; arr*=0.1 ; f.setArray(arr) ; f.checkConsistencyLight()
        ff=MEDFileField1TS() ; ff.setFieldNoProfileSBT(f)
        fmts=MEDFileFieldMultiTS() ; fmts.pushBackTimeStep(ff)
        #
        ms=MEDFileMeshes() ; ms.pushMesh(mm)
        fields=MEDFileFields() ; fields.pushField(fmts)
        ms.write(fname,2) ; fields.write(fname,0)
        #
        del mm,fmts,fields,ms
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        ms.cartesianizeMe()
        #
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        #
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),1)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDCurveLinearMeshMultiLev))
        a,b,c=mml2.buildVTUArrays()
        self.assertTrue(c)# the array is thoose in structure
        ref_a=DataArrayDouble([0.,0.,0.,0.,0.,0.8,0.,0.,1.6,0.,0.,0.,0.6928203230275509,0.,0.4,1.3856406460551018,0.,0.8,0.,0.,-0.,0.692820323027551,0.,-0.4,1.385640646055102,0.,-0.8,0.,0.,-0.,0.,0.,-0.8,0.,0.,-1.6,0.,0.,0.,0.,0.,0.8,0.,0.,1.6,0.,0.,0.,0.,0.6928203230275509,0.4,0.,1.3856406460551018,0.8,0.,0.,-0.,0.,0.692820323027551,-0.4,0.,1.385640646055102,-0.8,0.,0.,-0.,0.,0.,-0.8,0.,0.,-1.6,-0.,0.,0.,-0.,0.,0.8,-0.,0.,1.6,-0.,0.,0.,-0.6928203230275509,0.,0.4,-1.3856406460551018,0.,0.8,-0.,0.,-0.,-0.692820323027551,0.,-0.4,-1.385640646055102,0.,-0.8,-0.,0.,-0.,0.,0.,-0.8,0.,0.,-1.6,-0.,-0.,0.,-0.,-0.,0.8,-0.,-0.,1.6,-0.,-0.,0.,0.,-0.6928203230275509,0.4,0.,-1.3856406460551018,0.8,-0.,-0.,-0.,0.,-0.692820323027551,-0.4,0.,-1.385640646055102,-0.8,-0.,-0.,-0.,0.,0.,-0.8,0.,0.,-1.6,0.,-0.,0.,0.,-0.,0.8,0.,-0.,1.6,0.,-0.,0.,0.6928203230275509,0.,0.4,1.3856406460551018,0.,0.8,0.,-0.,-0.,0.692820323027551,0.,-0.4,1.385640646055102,0.,-0.8,0.,-0.,-0.,0.,0.,-0.8,0.,0.,-1.6],60,3)
        ref_a.setInfoOnComponents(comps)
        self.assertTrue(a.isEqual(ref_a,1e-14))
        self.assertEqual(b,[3,4,5])
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1):
            ffCell=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            self.assertTrue(v.isEqual(DataArrayDouble([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3]),1e-14))
            pass
        pass

    def test39(self):
        """Idem test37, test38, test39, test40 except that here it is an unstructured mesh."""
        fname="ForMEDReader39.med"
        meshName="mesh"
        description="Spherical grid"
        comps=["X [cm]","Y [cm]","Z [cm]"]
        arrX=DataArrayDouble(3) ; arrX.iota() ; arrX*=0.8 ; arrX.setInfoOnComponent(0,comps[0])
        arrY=DataArrayDouble(4) ; arrY.iota() ; arrY*=pi/(len(arrY)-1) ; arrY.setInfoOnComponent(0,comps[1])
        arrZ=DataArrayDouble(5) ; arrZ.iota() ; arrZ*=2*pi/(len(arrZ)-1) ; arrZ.setInfoOnComponent(0,comps[2])
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY,arrZ) ; m.setName(meshName) ; m=m.buildUnstructured()
        mm=MEDFileUMesh() ; mm[0]=m ; mm.setDescription(description) # the test is here : UMesh !
        mm.setAxisType(AX_SPHER) # the test is here !
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(m) ; f.setName("Field")
        arr=DataArrayDouble(m.getNumberOfCells()) ; arr.iota() ; arr*=0.1 ; f.setArray(arr) ; f.checkConsistencyLight()
        ff=MEDFileField1TS() ; ff.setFieldNoProfileSBT(f)
        fmts=MEDFileFieldMultiTS() ; fmts.pushBackTimeStep(ff)
        #
        ms=MEDFileMeshes() ; ms.pushMesh(mm)
        fields=MEDFileFields() ; fields.pushField(fmts)
        ms.write(fname,2) ; fields.write(fname,0)
        #
        del mm,fmts,fields,ms
        ms=MEDFileMeshes(fname)
        fields=MEDFileFields(fname,False)
        #
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        #
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),1)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(ncc)
        ref_a=DataArrayDouble([0.,0.,0.,0.,0.,0.8,0.,0.,1.6,0.,0.,0.,0.6928203230275509,0.,0.4,1.3856406460551018,0.,0.8,0.,0.,-0.,0.692820323027551,0.,-0.4,1.385640646055102,0.,-0.8,0.,0.,-0.,0.,0.,-0.8,0.,0.,-1.6,0.,0.,0.,0.,0.,0.8,0.,0.,1.6,0.,0.,0.,0.,0.6928203230275509,0.4,0.,1.3856406460551018,0.8,0.,0.,-0.,0.,0.692820323027551,-0.4,0.,1.385640646055102,-0.8,0.,0.,-0.,0.,0.,-0.8,0.,0.,-1.6,-0.,0.,0.,-0.,0.,0.8,-0.,0.,1.6,-0.,0.,0.,-0.6928203230275509,0.,0.4,-1.3856406460551018,0.,0.8,-0.,0.,-0.,-0.692820323027551,0.,-0.4,-1.385640646055102,0.,-0.8,-0.,0.,-0.,0.,0.,-0.8,0.,0.,-1.6,-0.,-0.,0.,-0.,-0.,0.8,-0.,-0.,1.6,-0.,-0.,0.,0.,-0.6928203230275509,0.4,0.,-1.3856406460551018,0.8,-0.,-0.,-0.,0.,-0.692820323027551,-0.4,0.,-1.385640646055102,-0.8,-0.,-0.,-0.,0.,0.,-0.8,0.,0.,-1.6,0.,-0.,0.,0.,-0.,0.8,0.,-0.,1.6,0.,-0.,0.,0.6928203230275509,0.,0.4,1.3856406460551018,0.,0.8,0.,-0.,-0.,0.692820323027551,0.,-0.4,1.385640646055102,0.,-0.8,0.,-0.,-0.,0.,0.,-0.8,0.,0.,-1.6],60,3)
        ref_a.setInfoOnComponents(comps)
        self.assertTrue(a0.isEqual(ref_a,1e-14))#<- Test is here
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1):
            ffCell=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            self.assertTrue(v.isEqual(DataArrayDouble([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3]),1e-14))
        pass

    def test40(self):
        """Idem test37, test38, test39, test40 except that here it is a CL mesh."""
        fname="ForMEDReader40.med"
        meshName="mesh"
        description="Spherical grid"
        comps=["X [cm]","Y [cm]","Z [cm]"]
        arrX=DataArrayDouble(3) ; arrX.iota() ; arrX*=0.8 ; arrX.setInfoOnComponent(0,comps[0])
        arrY=DataArrayDouble(4) ; arrY.iota() ; arrY*=pi/(len(arrY)-1) ; arrY.setInfoOnComponent(0,comps[1])
        arrZ=DataArrayDouble(5) ; arrZ.iota() ; arrZ*=2*pi/(len(arrZ)-1) ; arrZ.setInfoOnComponent(0,comps[2])
        m=MEDCouplingCMesh() ; m.setCoords(arrX,arrY,arrZ) ; m.setName(meshName) ; m=m.buildCurveLinear()
        mm=MEDFileCurveLinearMesh() ; mm.setMesh(m) ; mm.setDescription(description) # the test is here CLMesh!
        mm.setAxisType(AX_SPHER) # the test is here !
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(m) ; f.setName("Field")
        arr=DataArrayDouble(m.getNumberOfCells()) ; arr.iota() ; arr*=0.1 ; f.setArray(arr) ; f.checkConsistencyLight()
        ff=MEDFileField1TS() ; ff.setFieldNoProfileSBT(f)
        fmts=MEDFileFieldMultiTS() ; fmts.pushBackTimeStep(ff)
        #
        ms=MEDFileMeshes() ; ms.pushMesh(mm)
        fields=MEDFileFields() ; fields.pushField(fmts)
        ms.write(fname,2) ; fields.write(fname,0)
        #
        ms=MEDFileMeshes(fname) ; ms.cartesianizeMe()
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        #
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),1)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDCurveLinearMeshMultiLev))
        a,b,c=mml2.buildVTUArrays()
        self.assertTrue(c)
        ref_a=DataArrayDouble([0.,0.,0.,0.,0.,0.8,0.,0.,1.6,0.,0.,0.,0.6928203230275509,0.,0.4,1.3856406460551018,0.,0.8,0.,0.,-0.,0.692820323027551,0.,-0.4,1.385640646055102,0.,-0.8,0.,0.,-0.,0.,0.,-0.8,0.,0.,-1.6,0.,0.,0.,0.,0.,0.8,0.,0.,1.6,0.,0.,0.,0.,0.6928203230275509,0.4,0.,1.3856406460551018,0.8,0.,0.,-0.,0.,0.692820323027551,-0.4,0.,1.385640646055102,-0.8,0.,0.,-0.,0.,0.,-0.8,0.,0.,-1.6,-0.,0.,0.,-0.,0.,0.8,-0.,0.,1.6,-0.,0.,0.,-0.6928203230275509,0.,0.4,-1.3856406460551018,0.,0.8,-0.,0.,-0.,-0.692820323027551,0.,-0.4,-1.385640646055102,0.,-0.8,-0.,0.,-0.,0.,0.,-0.8,0.,0.,-1.6,-0.,-0.,0.,-0.,-0.,0.8,-0.,-0.,1.6,-0.,-0.,0.,0.,-0.6928203230275509,0.4,0.,-1.3856406460551018,0.8,-0.,-0.,-0.,0.,-0.692820323027551,-0.4,0.,-1.385640646055102,-0.8,-0.,-0.,-0.,0.,0.,-0.8,0.,0.,-1.6,0.,-0.,0.,0.,-0.,0.8,0.,-0.,1.6,0.,-0.,0.,0.6928203230275509,0.,0.4,1.3856406460551018,0.,0.8,0.,-0.,-0.,0.692820323027551,0.,-0.4,1.385640646055102,0.,-0.8,0.,-0.,-0.,0.,0.,-0.8,0.,0.,-1.6],60,3)
        ref_a.setInfoOnComponents(comps)
        self.assertTrue(a.isEqual(ref_a,1e-14))#<- Test is here
        self.assertEqual(b,[3,4,5])
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1):
            ffCell=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            self.assertTrue(v.isEqual(DataArrayDouble([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3]),1e-14))
        pass

    def test41(self):
        """This test focused on bug revealed with // load of multi nodes field with no profile. The error was the first node field (dataarray partdef) change the partdef for the others ! """
        fname="ForMEDReader41.med"
        meshName="mesh"
        nx=5
        arr=DataArrayDouble(nx) ; arr.iota()
        m=MEDCouplingCMesh() ; m.setCoords(arr,arr) ; m=m.buildUnstructured() ; m.setName(meshName)
        renum=DataArrayInt.Aggregate([DataArrayInt.Range(0,m.getNumberOfCells(),2),DataArrayInt.Range(1,m.getNumberOfCells(),2)])
        m=m[renum] # important think to renum if not we see nothing. The goal if to have dynamic_cast<DataPartDefinition>!=NULL
        mm=MEDFileUMesh() ; mm[0]=m
        mm.write(fname,2)
        f0=MEDCouplingFieldDouble(ON_NODES) ; f0.setMesh(m) ; f0.setName("aaa")
        arr0=DataArrayDouble(nx*nx) ; arr0.iota() ; f0.setArray(arr0)
        ff0=MEDFileField1TS() ; ff0.setFieldNoProfileSBT(f0)
        f1=MEDCouplingFieldDouble(ON_NODES) ; f1.setMesh(m) ; f1.setName("bbb")
        arr1=DataArrayDouble(nx*nx) ; arr1.iota() ; arr1+=100 ; f1.setArray(arr1)
        ff1=MEDFileField1TS() ; ff1.setFieldNoProfileSBT(f1)
        ff0.write(fname,0) ; ff1.write(fname,0)
        # 
        a=8 ; b=16
        ms=MEDFileMeshes()
        mm=MEDFileUMesh.LoadPartOf(fname,meshName,[NORM_QUAD4],[a,b,1],-1,-1)
        ms.pushMesh(mm)
        ms[0].zipCoords()
        ms.cartesianizeMe()
        fields=MEDFileFields.LoadPartOf(fname,False,ms);
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    if itmp.presenceOfMultiDiscPerGeoType():
                        tmp2=itmp.splitMultiDiscrPerGeoTypes()
                        for iii,itmp2 in enumerate(tmp2):
                            name="%s_%i"%(itmp2.getName(),iii)
                            itmp2.setName(name)
                            allFMTSLeavesToDisplay2.append(itmp2)
                            pass
                        pass
                    else:
                        allFMTSLeavesToDisplay2.append(itmp)
                        pass
                    pass
                allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
                pass
        # GO for reading in MEDReader, by not loading all. Mesh is fully loaded but not fields values
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        allFMTSLeavesPerCommonSupport=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesPerTimeSeries[0],ms[ms.getMeshesNames()[0]])
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=allFMTSLeavesPerCommonSupport[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        #
        f2=allFMTSLeavesPerCommonSupport[0][0][0][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(f2,mst)
        f2.loadArraysIfNecessary()
        v0=mml.buildDataArray(fsst,fields,f2.getUndergroundDataArray())
        assert(v0.isEqual(DataArrayDouble([1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24]),1e-12))
        #
        f2=allFMTSLeavesPerCommonSupport[0][0][1][0]
        fsst=MEDFileField1TSStructItem.BuildItemFrom(f2,mst)
        f2.loadArraysIfNecessary()
        v1=mml.buildDataArray(fsst,fields,f2.getUndergroundDataArray())
        assert(v1.isEqual(DataArrayDouble([101,102,103,104,106,107,108,109,111,112,113,114,116,117,118,119,121,122,123,124]),1e-12))
        pass

    def test42(self):
        """ EDF14869 - SEG4 """
        fname="ForMEDReader42.med"
        meshName="mesh"
        #
        a0exp=DataArrayDouble([0.,1.,0.3,0.7])
        m=MEDCouplingUMesh("mesh",1)
        m.setCoords(a0exp)
        m.allocateCells()
        m.insertNextCell(NORM_SEG4,[0,1,2,3])
        mm=MEDFileUMesh() ; mm[0]=m
        #
        f=MEDCouplingFieldDouble(ON_CELLS) ; f.setMesh(m) ; f.setName("Field")
        arr=DataArrayDouble(m.getNumberOfCells()) ; arr.iota() ; arr*=0.1 ; f.setArray(arr) ; f.checkConsistencyLight()
        ff=MEDFileField1TS() ; ff.setFieldNoProfileSBT(f)
        fmts=MEDFileFieldMultiTS() ; fmts.pushBackTimeStep(ff)
        #
        ms=MEDFileMeshes() ; ms.pushMesh(mm)
        fields=MEDFileFields() ; fields.pushField(fmts)
        ms.write(fname,2) ; fields.write(fname,0)
        #
        ms=MEDFileMeshes(fname) 
        fields=MEDFileFields(fname,False)
        fields.removeFieldsWithoutAnyTimeStep()
        fields_per_mesh=[fields.partOfThisLyingOnSpecifiedMeshName(meshName) for meshName in ms.getMeshesNames()]
        allFMTSLeavesToDisplay=[]
        for fields in fields_per_mesh:
            allFMTSLeavesToDisplay2=[]
            for fmts in fields:
                tmp=fmts.splitDiscretizations()
                for itmp in tmp:
                    self.assertTrue(not itmp.presenceOfMultiDiscPerGeoType())
                    pass
                allFMTSLeavesToDisplay2+=tmp
                pass
            allFMTSLeavesToDisplay.append(allFMTSLeavesToDisplay2)
            pass
        #
        self.assertEqual(len(allFMTSLeavesToDisplay),1)
        self.assertEqual(len(allFMTSLeavesToDisplay[0]),1)
        allFMTSLeavesPerTimeSeries=MEDFileAnyTypeFieldMultiTS.SplitIntoCommonTimeSeries(sum(allFMTSLeavesToDisplay,[]))
        self.assertEqual(len(allFMTSLeavesPerTimeSeries),1)
        allFMTSLeavesPerCommonSupport1=MEDFileAnyTypeFieldMultiTS.SplitPerCommonSupport(allFMTSLeavesToDisplay[0],ms[ms.getMeshesNames()[0]])
        self.assertEqual(len(allFMTSLeavesPerCommonSupport1),1)
        #
        mst=MEDFileMeshStruct.New(ms[0])
        fcscp=allFMTSLeavesPerCommonSupport1[0][1]
        mml=fcscp.buildFromScratchDataSetSupport(0,fields)
        mml2=mml.prepare()
        self.assertTrue(isinstance(mml2,MEDUMeshMultiLev))
        ncc,a0,a1,a2,a3,a4,a5=mml2.buildVTUArrays()
        self.assertTrue(not ncc)
        self.assertTrue(a0.isEqual(a0exp.changeNbOfComponents(3,0.),1e-12))
        self.assertTrue(a1.isEqual(DataArrayByte([35])))# VTK_CUBIC_LINE
        self.assertTrue(a2.isEqual(DataArrayInt([0])))
        self.assertTrue(a3.isEqual(DataArrayInt([4,0,1,2,3])))
        self.assertTrue(a4 is None)
        self.assertTrue(a5 is None)
        self.assertTrue(mml2.retrieveGlobalNodeIdsIfAny() is None)
        for i in range(1):
            ffCell=allFMTSLeavesPerCommonSupport1[0][0][0][i]
            fsst=MEDFileField1TSStructItem.BuildItemFrom(ffCell,mst)
            ffCell.loadArraysIfNecessary()
            v=mml2.buildDataArray(fsst,fields,ffCell.getUndergroundDataArray())
            self.assertEqual(v.getHiddenCppPointer(),ffCell.getUndergroundDataArray().getHiddenCppPointer())
            self.assertTrue(v.isEqual(DataArrayDouble([0.0]),1e-14))
        pass
    
    pass

if __name__ == "__main__":
  unittest.main()
