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

from MEDCoupling import *
import unittest

class MEDCouplingBasicsTest(unittest.TestCase):
    def testExampleFieldDoubleBuildSubPart1(self):
        from MEDCouplingDataForTest import MEDCouplingDataForTest
#! [PySnippetFieldDoubleBuildSubPart1_1]
        mesh1=MEDCouplingDataForTest.build2DTargetMesh_1()
        f1=MEDCouplingFieldDouble.New(ON_CELLS,ONE_TIME)
        f1.setTime(2.3,5,6)
        f1.setMesh(mesh1)
        array=DataArrayDouble.New()
        arr1=[3.,103.,4.,104.,5.,105.,6.,106.,7.,107.]
        array.setValues(arr1,mesh1.getNumberOfCells(),2)
        f1.setArray(array)
# ! [PySnippetFieldDoubleBuildSubPart1_1]
# ! [PySnippetFieldDoubleBuildSubPart1_2]
        part1=[2,1,4]
        f2=f1.buildSubPart(part1)
# ! [PySnippetFieldDoubleBuildSubPart1_2]
        f2.zipCoords()
        self.failUnlessEqual(3,f2.getNumberOfTuples())
        self.failUnlessEqual(2,f2.getNumberOfComponents())
        expected1=[5.,105.,4.,104.,7.,107.]
        for i in xrange(6):
            self.assertAlmostEqual(f2.getIJ(0,i),expected1[i],12)
            pass
        self.failUnlessEqual(3,f2.getMesh().getNumberOfCells())
        self.failUnlessEqual(6,f2.getMesh().getNumberOfNodes())
        self.failUnlessEqual(2,f2.getMesh().getSpaceDimension())
        self.failUnlessEqual(2,f2.getMesh().getMeshDimension())
        m2C=f2.getMesh()
        self.failUnlessEqual(13,m2C.getMeshLength())
        expected2=[0.2, -0.3, 0.7, -0.3, 0.2, 0.2, 0.7, 0.2, 0.2, 0.7, 0.7, 0.7]
        for i in xrange(12):
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12)
            pass
        expected3=[3,2,3,1,3,0,2,1,4,4,5,3,2]
        self.failUnlessEqual(expected3,list(m2C.getNodalConnectivity().getValues()))
        expected4=[0,4,8,13]
        self.failUnlessEqual(expected4,list(m2C.getNodalConnectivityIndex().getValues()))
        # Test with field on nodes.
# ! [PySnippetFieldDoubleBuildSubPart1_3]
        f1=MEDCouplingFieldDouble.New(ON_NODES,ONE_TIME)
        f1.setTime(2.3,5,6)
        f1.setMesh(mesh1)
        array=DataArrayDouble.New()
        arr2=[3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.]
        array.setValues(arr2,mesh1.getNumberOfNodes(),2)
        f1.setArray(array)
# ! [PySnippetFieldDoubleBuildSubPart1_3]
# ! [PySnippetFieldDoubleBuildSubPart1_4]
        part2=[1,2]
        f2=f1.buildSubPart(part2)
# ! [PySnippetFieldDoubleBuildSubPart1_4]
        self.failUnlessEqual(4,f2.getNumberOfTuples())
        self.failUnlessEqual(2,f2.getNumberOfComponents())
        expected5=[4.,104.,5.,105.,7.,107.,8.,108.]
        for i in xrange(8):
            self.assertAlmostEqual(f2.getIJ(0,i),expected5[i],12)
            pass
        self.failUnlessEqual(2,f2.getMesh().getNumberOfCells())
        self.failUnlessEqual(4,f2.getMesh().getNumberOfNodes())
        self.failUnlessEqual(2,f2.getMesh().getSpaceDimension())
        self.failUnlessEqual(2,f2.getMesh().getMeshDimension())
        m2C=f2.getMesh()
        self.failUnlessEqual(8,m2C.getMeshLength())
        for i in xrange(8):#8 is not an error
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12)
            pass
        self.failUnlessEqual(expected3[:4],list(m2C.getNodalConnectivity().getValues())[4:])
        self.failUnlessEqual(expected3[4:8],list(m2C.getNodalConnectivity().getValues())[:4])
        self.failUnlessEqual(expected4[:3],list(m2C.getNodalConnectivityIndex().getValues()))
        #idem previous because nodes of cell#4 are not fully present in part3
        part3=[1,2]
        arrr=DataArrayInt.New()
        arrr.setValues(part3,2,1)
        f2=f1.buildSubPart(arrr)
        self.failUnlessEqual(4,f2.getNumberOfTuples())
        self.failUnlessEqual(2,f2.getNumberOfComponents())
        for i in xrange(8):
            self.assertAlmostEqual(f2.getIJ(0,i),expected5[i],12)
            pass
        self.failUnlessEqual(2,f2.getMesh().getNumberOfCells())
        self.failUnlessEqual(4,f2.getMesh().getNumberOfNodes())
        self.failUnlessEqual(2,f2.getMesh().getSpaceDimension())
        self.failUnlessEqual(2,f2.getMesh().getMeshDimension())
        m2C=f2.getMesh()
        self.failUnlessEqual(8,m2C.getMeshLength())
        for i in xrange(8):#8 is not an error
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12)
            pass
        self.failUnlessEqual(expected3[:4],list(m2C.getNodalConnectivity().getValues())[4:8])
        self.failUnlessEqual(expected3[4:8],list(m2C.getNodalConnectivity().getValues())[:4])
        self.failUnlessEqual(expected4[:3],list(m2C.getNodalConnectivityIndex().getValues()))
        part4=[1,2,4]
        f2=f1.buildSubPart(part4)
        self.failUnlessEqual(6,f2.getNumberOfTuples())
        self.failUnlessEqual(2,f2.getNumberOfComponents())
        expected6=[4.,104.,5.,105.,7.,107.,8.,108.,10.,110.,11.,111.]
        for i in xrange(12):
            self.assertAlmostEqual(f2.getIJ(0,i),expected6[i],12)
            pass
        self.failUnlessEqual(3,f2.getMesh().getNumberOfCells())
        self.failUnlessEqual(6,f2.getMesh().getNumberOfNodes())
        self.failUnlessEqual(2,f2.getMesh().getSpaceDimension())
        self.failUnlessEqual(2,f2.getMesh().getMeshDimension())
        m2C=f2.getMesh()
        self.failUnlessEqual(13,m2C.getMeshLength())
        for i in xrange(12):
            self.assertAlmostEqual(expected2[i],m2C.getCoords().getIJ(0,i),12)
            pass
        self.failUnlessEqual(expected3[0:4],list(m2C.getNodalConnectivity().getValues())[4:8])
        self.failUnlessEqual(expected3[4:8],list(m2C.getNodalConnectivity().getValues())[0:4])
        self.failUnlessEqual(expected3[8:13],list(m2C.getNodalConnectivity().getValues())[8:13])
        self.failUnlessEqual(expected4,list(m2C.getNodalConnectivityIndex().getValues()))
        pass

    def testExampleUMeshStdBuild1(self):
# ! [PySnippetUMeshStdBuild1_1]
        coords=[-0.3,-0.3,0.,   0.2,-0.3,0.,   0.7,-0.3,0.,   -0.3,0.2,0.,   0.2,0.2,0., 
                 0.7,0.2,0.,    -0.3,0.7,0.,    0.2,0.7,0.,     0.7,0.7,0. ]
        nodalConnPerCell=[0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4]
# ! [PySnippetUMeshStdBuild1_1]
# ! [PySnippetUMeshStdBuild1_2]
        mesh=MEDCouplingUMesh.New("My2DMesh",2)
# ! [PySnippetUMeshStdBuild1_2]
# ! [PySnippetUMeshStdBuild1_3]
        mesh.allocateCells(5)#You can put more than 5 if you want but not less.
        mesh.insertNextCell(NORM_QUAD4,nodalConnPerCell[:4])
        mesh.insertNextCell(NORM_TRI3,nodalConnPerCell[4:7])
        mesh.insertNextCell(NORM_TRI3,nodalConnPerCell[7:10])
        mesh.insertNextCell(NORM_QUAD4,nodalConnPerCell[10:14])
        mesh.insertNextCell(NORM_QUAD4,nodalConnPerCell[14:])
        mesh.finishInsertingCells()
# ! [PySnippetUMeshStdBuild1_3]
# ! [PySnippetUMeshStdBuild1_4]
        myCoords=DataArrayDouble.New(coords,9,3)#here myCoords are declared to have 3 components, mesh will deduce that its spaceDim==3. 
        mesh.setCoords(myCoords)#myCorrds contains 9 tuples, that is to say mesh contains 9 nodes.
# ! [PySnippetUMeshStdBuild1_4]
# ! [PySnippetUMeshStdBuild1_5]
# ! [PySnippetUMeshStdBuild1_5]
        mesh.checkCoherency()
        pass

    def testExampleUMeshAdvBuild1(self):
# ! [PySnippetUMeshAdvBuild1_1]
        coords=[-0.3,-0.3,0.,   0.2,-0.3,0.,   0.7,-0.3,0.,   -0.3,0.2,0.,   0.2,0.2,0., 
                 0.7,0.2,0.,    -0.3,0.7,0.,    0.2,0.7,0.,     0.7,0.7,0. ]
        nodalConnPerCell=[4,0,3,4,1, 3,1,4,2, 3,4,5,2, 4,6,7,4,3, 4,7,8,5,4]
        nodalConnPerCellIndex=[0,5,9,13,18,23]
# ! [PySnippetUMeshAdvBuild1_1]
# ! [PySnippetUMeshAdvBuild1_2]
        mesh=MEDCouplingUMesh.New("My2DMesh",2)
# ! [PySnippetUMeshAdvBuild1_2]
# ! [PySnippetUMeshAdvBuild1_3]
        nodalConn=DataArrayInt.New(nodalConnPerCell,23,1)
        nodalConnI=DataArrayInt.New(nodalConnPerCellIndex,6,1)
        mesh.setConnectivity(nodalConn,nodalConnI,True)
# ! [PySnippetUMeshAdvBuild1_3]
# ! [PySnippetUMeshAdvBuild1_4]
        myCoords=DataArrayDouble.New(coords,9,3)#here myCoords are declared to have 3 components, mesh will deduce that its spaceDim==3. 
        mesh.setCoords(myCoords)#myCorrds contains 9 tuples, that is to say mesh contains 9 nodes.
# ! [PySnippetUMeshAdvBuild1_4]
# ! [PySnippetUMeshAdvBuild1_5]
# ! [PySnippetUMeshAdvBuild1_5]
        mesh.checkCoherency()
        pass
    pass

unittest.main()
