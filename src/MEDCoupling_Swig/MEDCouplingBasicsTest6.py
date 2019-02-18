#  -*- coding: utf-8 -*-
# Copyright (C) 2017-2019  CEA/DEN, EDF R&D
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


import sys
if sys.platform == "win32":
    from MEDCouplingCompat import *
else:
    from medcoupling import *
import unittest
from math import pi,e,sqrt,cos,sin
from datetime import datetime
import rlcompleter,readline # this line has to be here, to ensure a usability of MEDCoupling/MEDLoader. B4 removing it please notify to anthony.geay@edf.fr

class MEDCouplingBasicsTest6(unittest.TestCase):
    def testPointSetInvertOrientationOfAllCells(self):
        """ Test of inversion of orientation of cells on a different geo types"""
        def level1(self):
            m=MEDCouplingUMesh("",1)
            m.allocateCells()
            m.insertNextCell(NORM_SEG2,[3,4])
            m.insertNextCell(NORM_SEG2,[13,14])
            m.insertNextCell(NORM_SEG3,[5,6,7])
            m.insertNextCell(NORM_SEG2,[23,24])
            m.insertNextCell(NORM_SEG3,[8,9,10])
            ref0=DataArrayInt([0,3,6,10,13,17])
            self.assertTrue(m.getNodalConnectivityIndex().isEqual(ref0))
            m.invertOrientationOfAllCells()
            self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([1,4,3, 1,14,13, 2,7,6,5, 1,24,23, 2,10,9,8])))
            self.assertTrue(m.getNodalConnectivityIndex().isEqual(ref0))
            pass

        def level2(self):
            m=MEDCouplingUMesh("",2)
            m.allocateCells()
            m.insertNextCell(NORM_TRI3,[1,2,3])
            m.insertNextCell(NORM_QUAD4,[4,5,6,7])
            m.insertNextCell(NORM_POLYGON,[8,9,10,11,12,13])
            m.insertNextCell(NORM_TRI6,[14,15,16,17,18,19])
            m.insertNextCell(NORM_QUAD8,[20,21,22,23,24,25,26,27])
            m.insertNextCell(NORM_QPOLYG,[30,31,32,33,34,35, 36,37,38,39,40,41])
            ref0=DataArrayInt([0,4,9,16,23,32,45])
            self.assertTrue(m.getNodalConnectivityIndex().isEqual(ref0))
            m.invertOrientationOfAllCells()
            self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([3,1,3,2, 4,4,7,6,5, 5,8,13,12,11,10,9, 6,14,16,15,19,18,17, 8,20,23,22,21,27,26,25,24, 32,30,35,34,33,32,31,41,40,39,38,37,36])))
            self.assertTrue(m.getNodalConnectivityIndex().isEqual(ref0))
            pass

        def level3(self):
            m=MEDCouplingUMesh("",3)
            m.allocateCells()
            m.insertNextCell(NORM_TETRA4,[1,2,3,4])
            m.insertNextCell(NORM_PYRA5,[5,6,7,8,9])
            m.insertNextCell(NORM_PENTA6,[10,11,12,13,14,15])
            m.insertNextCell(NORM_HEXA8,[20,21,22,23,24,25,26,27])
            m.insertNextCell(NORM_TETRA10,[30,31,32,33,34,35,36,37,38,39])
            m.insertNextCell(NORM_PYRA13,[40,41,42,43,44,45,46,47,48,49,50,51,52])
            m.insertNextCell(NORM_HEXA20,[60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79])
            m.insertNextCell(NORM_PENTA15,[80,81,82,83,84,85,86,87,88,89,90,91,92,93,94])
            ref0=DataArrayInt([0,5,11,18,27,38,52,73,89])
            self.assertTrue(m.getNodalConnectivityIndex().isEqual(ref0))
            m.invertOrientationOfAllCells()
            self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([14,1,3,2,4, 15,5,8,7,6,9, 16,10,12,11,13,15,14, 18,20,23,22,21,24,27,26,25, 20,30,32,31,33,36,35,34,37,39,38, 23,40,43,42,41,44,48,47,46,45,49,52,51,50, 30,60,63,62,61,64,67,66,65,71,70,69,68,75,74,73,72,76,79,78,77, 25,80,82,81,83,85,84,88,87,86,91,90,89,92,94,93])))
            self.assertTrue(m.getNodalConnectivityIndex().isEqual(ref0))
            pass

        def gtumesh(self):
            m=MEDCoupling1SGTUMesh("",NORM_SEG2)
            m.setNodalConnectivity(DataArrayInt([1,2,3,4,5,6,7,8]))
            self.assertEqual(m.getNumberOfCells(),4)
            m2=m.deepCopy()
            self.assertTrue(m2.isEqual(m,0.))
            m.invertOrientationOfAllCells()
            self.assertTrue(not m2.isEqual(m,0.))
            m.getNodalConnectivity().isEqual(DataArrayInt([2,1,4,3,6,5,8,7]))
            m.invertOrientationOfAllCells()
            self.assertTrue(m2.isEqual(m,0.))
            #
            p=MEDCoupling1DGTUMesh("",NORM_POLYGON)
            ref0=DataArrayInt([0,3,7,12])
            p.setNodalConnectivity(DataArrayInt([1,2,3, 10,11,12,13, 20,21,22,23,24]),ref0)
            p2=p.deepCopy()
            self.assertTrue(p2.isEqual(p,0.))
            self.assertEqual(p.getNumberOfCells(),3)
            p.invertOrientationOfAllCells()
            self.assertTrue(not p2.isEqual(p,0.))
            self.assertTrue(p.getNodalConnectivityIndex().isEqual(ref0))
            self.assertTrue(p.getNodalConnectivity().isEqual(DataArrayInt([1,3,2, 10,13,12,11, 20,24,23,22,21])))
            p.invertOrientationOfAllCells()
            self.assertTrue(p2.isEqual(p,0.))
            pass
        level1(self)
        level2(self)
        level3(self)
        gtumesh(self)
        pass

    def testPenta18_1(self):
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
        self.assertTrue(m.getMeasureField(True).getArray().isEqual(DataArrayDouble([0.5]),1e-12))
        #
        f=MEDCouplingFieldDouble(ON_NODES)
        f.setMesh(m)
        f.setName("FieldOnPenta18")
        f.setArray(DataArrayDouble(list(range(18))))
        f.checkConsistencyLight()
        #
        m2,d,di,rd,rdi=m.buildDescendingConnectivity()
        self.assertTrue(m2.getNodalConnectivity().isEqual(DataArrayInt([6,0,1,2,6,7,8,6,3,5,4,11,10,9,9,0,3,4,1,12,9,13,6,15,9,1,4,5,2,13,10,14,7,16,9,2,4,5,0,14,11,12,8,17])))
        self.assertTrue(m2.getNodalConnectivityIndex().isEqual(DataArrayInt([0,7,14,24,34,44])))
        self.assertTrue(d.isEqual(DataArrayInt([0,1,2,3,4])))
        self.assertTrue(di.isEqual(DataArrayInt([0,5])))
        self.assertTrue(rd.isEqual(DataArrayInt([0,0,0,0,0])))
        self.assertTrue(rdi.isEqual(DataArrayInt([0,1,2,3,4,5])))
        #
        f2=MEDCouplingFieldDouble(ON_NODES)
        f2.setMesh(m)
        f2.setName("FieldOnPenta18Sub")
        f2.setArray(DataArrayDouble(list(range(18))))
        f2.checkConsistencyLight()
        pass

    def testSKLAReplaceDeletePacks(self):
        index=DataArrayInt([0,3,5,6,6])
        value=DataArrayInt([1,2,3, 2,3, 3  ])
        sla=MEDCouplingSkyLineArray(index,value)
        idx=DataArrayInt([0,3])
        packs=[DataArrayInt([4,5]),DataArrayInt([6,7,8])]
        sla.replaceSimplePacks(idx,packs)
        self.assertTrue(sla.getIndexArray().isEqual(DataArrayInt([0,2,4,5,8])))
        self.assertTrue(sla.getValuesArray().isEqual(DataArrayInt([4,5, 2,3, 3, 6,7,8])))
        sla.deleteSimplePacks(idx)
        self.assertTrue(sla.getIndexArray().isEqual(DataArrayInt([0,2,3])))
        self.assertTrue(sla.getValuesArray().isEqual(DataArrayInt([2,3, 3])))
        sla.deleteSimplePack(1)
        self.assertTrue(sla.getIndexArray().isEqual(DataArrayInt([0,2])))
        self.assertTrue(sla.getValuesArray().isEqual(DataArrayInt([2,3])))
        pass

    def testDADAsArcOfCircle(self):
        d=DataArrayDouble([3.06915124862645,2.1464466094067824,2.85355345827285,2.3620444674400574,2.637955532559882,2.1464467447661937],3,2)
        center,radius,ang=d.asArcOfCircle()
        self.assertTrue((d-center).magnitude().isUniform(radius,1e-10))
        self.assertAlmostEqual(ang,-4.712389294301196,12)
        pass

    def testDAMaxAbsValue(self):
        d=DataArrayDouble([-2,3,1.2,-2.9])
        a,b=d.getMaxAbsValue()
        self.assertAlmostEqual(a,3.,13)
        self.assertEqual(b,1)
        a,b=(-d).getMaxAbsValue()
        self.assertAlmostEqual(a,-3.,13)
        self.assertEqual(b,1)
        self.assertAlmostEqual((-d).getMaxAbsValueInArray(),-3.,13)
        pass

    def testDAIFindIdForEach1(self):
        a1=DataArrayInt([17,27,2,10,-4,3,12,27,16])
        b1=DataArrayInt([3,16,-4,27,17])
        ret=a1.findIdForEach(b1)
        self.assertTrue(ret.isEqual(DataArrayInt([5,8,4,7,0])))
        self.assertTrue(a1[ret].isEqual(b1))
        b2=DataArrayInt([3,16,22,27,17])
        self.assertRaises(InterpKernelException,a1.findIdForEach,b2) # 22 not in a1 !
        a1.rearrange(3)
        self.assertRaises(InterpKernelException,a1.findIdForEach,b1) # a1 is not single component
        pass

    def testAttractSeg3MidPtsAroundNodes1(self):
        """ Test of MEDCouplingUMesh.attractSeg3MidPtsAroundNodes methods """
        ptsExpToBeModified=DataArrayInt([95,96,97,98,101,103,104,106,108,110])
        eps=1e-12
        a=2./3.
        b=1./3.
        coo=DataArrayDouble([10,0,0,10,10,0,10,0,3.+b,10,0,6.+a,10,10,3.+b,10,10,6.+a,10,3.+b,0,10,6.+a,0,3.+b,0,0,6.+a,0,0,3.+b,10,0,6.+a,10,0,10,3.+b,6.+a,10,6.+a,6.+a,10,3.+b,3.+b,10,6.+a,3.+b,3.+b,0,3.+b,3.+b,0,6.+a,6.+a,0,3.+b,6.+a,0,6.+a,6.+a,10,3.+b,6.+a,10,6.+a,3.+b,10,3.+b,3.+b,10,6.+a,3.+b,3.+b,0,6.+a,3.+b,0,3.+b,6.+a,0,6.+a,6.+a,0,3.+b,3.+b,6.+a,6.+a,3.+b,6.+a,3.+b,6.+a,6.+a,6.+a,6.+a,6.+a,3.+b,3.+b,3.+b,6.+a,3.+b,3.+b,3.+b,6.+a,3.+b,6.+a,6.+a,3.+b,10,0,1.+a,10,0,5.,10,10,1.+a,10,10,5.,10,1.+a,0,10,5.,0,10,8.+b,0,5.,0,0,8.+b,0,0,5.,10,0,8.+b,10,0,10,1.+a,6.+a,10,5.,6.+a,10,8.+b,6.+a,10,1.+a,3.+b,10,3.+b,5.,10,5.,3.+b,10,6.+a,5.,10,8.+b,3.+b,10,3.+b,1.+a,10,6.+a,1.+a,3.+b,0,1.+a,3.+b,0,5.,6.+a,0,1.+a,5.,0,3.+b,6.+a,0,5.,5.,0,6.+a,8.+b,0,3.+b,8.+b,0,6.+a,6.+a,10,1.+a,8.+b,10,3.+b,6.+a,10,5.,8.+b,10,6.+a,3.+b,10,1.+a,5.,10,3.+b,3.+b,10,5.,5.,10,6.+a,3.+b,1.+a,0,5.,3.+b,0,6.+a,1.+a,0,8.+b,3.+b,0,3.+b,5.,0,5.,6.+a,0,6.+a,5.,0,8.+b,6.+a,0,3.+b,8.+b,0,6.+a,8.+b,0,3.+b,1.+a,6.+a,6.+a,1.+a,6.+a,5.,3.+b,6.+a,8.+b,3.+b,6.+a,3.+b,5.,6.+a,6.+a,5.,6.+a,5.,6.+a,6.+a,8.+b,6.+a,6.+a,3.+b,8.+b,6.+a,6.+a,8.+b,6.+a,3.+b,3.+b,5,3.+b,1.+a,3.+b,6.+a,3.+b,5.,6.+a,1.+a,3.+b,5.,3.+b,3.+b,8.+b,3.+b,3.+b,3.+b,6.+a,5.,3.+b,5.,3.+b,6.+a,6.+a,5.,6.+a,5.,3.+b,5.,6.+a,3.+b,8.+b,6.+a,3.+b,3.+b,8.+b,3.+b,6.+a,8.+b,3.+b,3.+b,3.+b,1.+a,6.+a,3.+b,1.+a,3.+b,6.+a,1.+a,6.+a,6.+a,1.+a],111,3)
        conn=DataArrayInt([30,17,28,32,16,19,29,33,18,83,93,94,58,84,95,96,61,62,85,97,60,30,19,29,33,18,3,12,14,2,84,95,96,61,47,51,50,37,64,86,98,63,30,28,30,34,32,29,31,35,33,87,99,100,93,88,101,102,95,85,89,103,97,30,29,31,35,33,12,13,15,14,88,101,102,95,48,53,52,51,86,90,104,98,30,30,23,22,34,31,21,20,35,91,71,105,99,92,67,106,101,89,72,70,103,30,31,21,20,35,13,5,4,15,92,67,106,101,49,39,54,53,90,68,66,104,30,16,32,24,8,18,33,25,9,94,107,73,57,96,108,75,59,60,97,74,43,30,18,33,25,9,2,14,6,0,96,108,75,59,50,55,40,36,63,98,76,44,30,32,34,26,24,33,35,27,25,100,109,77,107,102,110,79,108,97,103,78,74,30,33,35,27,25,14,15,7,6,102,110,79,108,52,56,41,55,98,104,80,76,30,34,22,10,26,35,20,11,27,105,69,81,109,106,65,82,110,103,70,45,78,30,35,20,11,27,15,4,1,7,106,65,82,110,54,38,42,56,104,66,46,80])
        connI=DataArrayInt([0,21,42,63,84,105,126,147,168,189,210,231,252])
        m=MEDCouplingUMesh("mesh",3)
        m.setConnectivity(conn,connI,True)
        m.setCoords(coo.deepCopy())# deep copy coo because next line is going to modify it, if it works normaly
        m.attractSeg3MidPtsAroundNodes(0.1,DataArrayInt([33,35])) # ze call is here !
        self.assertTrue(not m.getCoords().isEqual(coo,eps)) # some points have had their position changed...
        ptsExpNotToBeModified=ptsExpToBeModified.buildComplement(len(coo))
        self.assertTrue(m.getCoords()[ptsExpNotToBeModified].isEqual(coo[ptsExpNotToBeModified],eps))
        self.assertTrue((m.getCoords()[ptsExpToBeModified]-coo[ptsExpToBeModified]).magnitude().isUniform(4./3.,1e-12))
        ptsPosExp=DataArrayDouble([6.+a,3.+b,3.+a,6.+a,3.,3.+b,6.+b,3.+b,3.+b,7.,3.+b,3.+b,6.+a,6.+a,3.+a,6.+b,6.+a,3.+b,7.,6.+a,3.+b,6.+a,7.,3.+b,6.+a,3.+b,3.,6.+a,6.+a,3.],10,3)
        self.assertTrue(m.getCoords()[ptsExpToBeModified].isEqual(ptsPosExp,1e-12))
        pass

    def testRenumberNodesInConnOpt(self):
        """ Test of MEDCouplingPointSet.renumberNodesInConn with map as input coming from DataArrayInt.invertArrayN2O2O2NOptimized
        """
        m=MEDCouplingUMesh("mesh",2)
        m.allocateCells()
        m.insertNextCell(NORM_QUAD4,[10000,10002,10001,10003])
        coo=DataArrayDouble([(0,0),(1,1),(1,0),(0,1)])
        m.setCoords(coo)
        m.checkConsistencyLight()
        #
        d=DataArrayInt([10000,10001,10002,10003])
        myMap=d.invertArrayN2O2O2NOptimized()
        myMap2=d.giveN2OOptimized()
        m.checkConsistencyLight()
        #
        m.renumberNodesInConn(myMap) # <- test is here for UMesh
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([4,0,2,1,3])))
        m.renumberNodesInConn(myMap2) # <- test is here for UMesh
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([4,10000,10002,10001,10003])))
        #
        m=MEDCoupling1SGTUMesh("mesh",NORM_QUAD4)
        m.setNodalConnectivity(DataArrayInt([10000,10002,10001,10003]))
        m.setCoords(coo)
        m.checkConsistencyLight()
        m.renumberNodesInConn(myMap) # <- test is here for 1SGTUMesh
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([0,2,1,3])))
        #
        m=MEDCoupling1DGTUMesh("mesh",NORM_POLYGON)
        m.setCoords(coo)
        m.setNodalConnectivity(DataArrayInt([10000,10002,10001,10003]),DataArrayInt([0,4]))
        m.checkConsistencyLight()
        m.renumberNodesInConn(myMap) # <- test is here for 1DGTUMesh
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([0,2,1,3])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4])))
        pass

    def testSeg2bGP(self):
        """Test of Gauss points on SEG2 using SEG2B style as ref coords
        """
        coo=DataArrayDouble([[0.,0.,0.],[1.,1.,1.]])
        m=MEDCouplingUMesh("mesh",1) ; m.setCoords(coo)
        m.allocateCells()
        # the cell description is exactly those described in the description of HEXA27 in MED file 3.0.7 documentation
        m.insertNextCell(NORM_SEG2,[0,1])
        refCoo=[0.,1.]
        weights=[0.8,0.1,0.1]
        gCoords=[0.2,0.5,0.9]
        fGauss=MEDCouplingFieldDouble(ON_GAUSS_PT) ; fGauss.setName("fGauss")
        fGauss.setMesh(m)
        fGauss.setGaussLocalizationOnType(NORM_SEG2,refCoo,gCoords,weights)
        arr=DataArrayDouble(fGauss.getNumberOfTuplesExpected()) ; arr.iota()
        fGauss.setArray(arr)
        fGauss.checkConsistencyLight()
        arrOfDisc=fGauss.getLocalizationOfDiscr()
        self.assertTrue(arrOfDisc.isEqual(DataArrayDouble([0.2,0.2,0.2,0.5,0.5,0.5,0.9,0.9,0.9],3,3),1e-12))
        pass

    def testUMeshGetCellsContainingPtOn2DNonDynQuadraticCells(self):
        """getCellsContainingPoint is now dealing curves of quadratic 2D elements.
This test is a mesh containing 2 QUAD8 cells. The input point is located at a special loc.
If true geometry (with curve as edges) is considered the result of getCellsContainingPoint is not the same as if only linear part of cells is considered."""
        coords=DataArrayDouble([-0.9428090415820631,0.9428090415820631,-1.06066017177982,1.06066017177982,-1.1785113019775801,1.1785113019775801,-1.2963624321753402,1.2963624321753402,-1.4142135623731,1.41421356237309,-0.7653668647301801,1.8477590650225701,-0.6378057206084831,1.53979922085214,-0.510244576486786,1.23183937668172,-0.701586292669331,1.6937791429373599,-0.574025148547635,1.38581929876693,-0.9259503883660041,1.38578268717091,-0.740760310692803,1.10862614973673,-1.1111404660392,1.66293922460509],13,2)
        m=MEDCouplingUMesh("mesh",2)
        m.setCoords(coords)
        m.allocateCells()
        m.insertNextCell(NORM_QUAD8,[4,2,6,5,3,10,8,12])
        m.insertNextCell(NORM_QUAD8,[2,0,7,6,1,11,9,10])
        #
        zePt=DataArrayDouble([-0.85863751450784975,1.4203162316045934],1,2)
        a,b=m.getCellsContainingPoints(zePt,1e-12)
        self.assertTrue(b.isEqual(DataArrayInt([0,1])))
        self.assertTrue(a.isEqual(DataArrayInt([1]))) # <- test is here. 0 if only linear parts are considered.
        #
        a,b=m.getCellsContainingPointsLinearPartOnlyOnNonDynType(zePt,1e-12)
        self.assertTrue(b.isEqual(DataArrayInt([0,1])))
        self.assertTrue(a.isEqual(DataArrayInt([0]))) # <- like before
        pass

    def testComputeIntegralOfSeg2IntoTri3_1(self):
        seg2 = [(90453.702115782813, 36372.66281307926), (90457.969790110554, 36373.365088601546)]
        tri3 = [(90466.90625, 36376.9375), (90446.5, 36404), (90453.1875, 36365.75)]
        a,b=DataArrayDouble.ComputeIntegralOfSeg2IntoTri3(seg2,tri3)
        self.assertEqual(len(a),3)
        self.assertAlmostEqual(a[0],0.2460689650955214,12)
        self.assertAlmostEqual(a[1],0.10875598777133343,12)
        self.assertAlmostEqual(a[2],0.6451750471331451,12)
        self.assertAlmostEqual(b,4.32507052854159,12)
        pass

    def testRemoveDegenerated1DCells1(self):
        m=MEDCoupling1SGTUMesh("mesh",NORM_SEG2)
        conn=DataArrayInt([1,2, 3,4, 5,5, 5,6, 6,6, 6,7, 19,19, 7,8])
        m.setNodalConnectivity(conn) # no coords set. It s not a bug. removeDegenerated1DCells doesn't care
        m=m.buildUnstructured()
        aa=m.getNodalConnectivity().getHiddenCppPointer()
        self.assertTrue(m.removeDegenerated1DCells()) # <- test is here
        bb=m.getNodalConnectivity().getHiddenCppPointer()
        self.assertNotEqual(aa,bb)
        expConn=DataArrayInt([1,1,2,1,3,4,1,5,6,1,6,7,1,7,8])
        expConnI=DataArrayInt.Range(0,16,3)
        self.assertTrue(m.getNodalConnectivity().isEqual(expConn))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(expConnI))
        self.assertTrue(not m.removeDegenerated1DCells())
        cc=m.getNodalConnectivity().getHiddenCppPointer()
        self.assertEqual(bb,cc)
        self.assertTrue(m.getNodalConnectivity().isEqual(expConn))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(expConnI))
        pass

    def testMergeFieldsOnGauss1(self):
        mName="mesh"
        fieldName="field"
        #
        _a=0.446948490915965;
        _b=0.091576213509771;
        _p1=0.11169079483905;
        _p2=0.0549758718227661;
        refCoo1=[ 0.,0., 1.,0., 0.,1. ]
        gsCoo1=[ 2*_b-1, 1-4*_b, 2*_b-1, 2.07*_b-1, 1-4*_b,
                 2*_b-1, 1-4*_a, 2*_a-1, 2*_a-1, 1-4*_a, 2*_a-1, 2*_a-1 ]
        wg1=[ 4*_p2, 4*_p2, 4*_p2, 4*_p1, 4*_p1, 4*_p1 ]
        #
        refCoo2=[ -1.,-1., 1.,-1., 1.,1., -1.,1. ]
        gsCoo2=[0.1,0.1, 0.2,0.2, 0.5,0.5, 0.6,0.6, 0.7,0.7]
        wg2=[0.1,0.2,0.3,0.4,0.5]
        #
        coo=DataArrayDouble([0,0,1,0,2,0,0,1,1,1,2,1,0,2,1,2,2,2],9,2)
        m1=MEDCouplingUMesh(mName,2)
        m1.allocateCells() ; m1.setCoords(coo)
        m1.insertNextCell(NORM_TRI3,[1,4,2])
        m1.insertNextCell(NORM_TRI3,[4,5,2])
        m1.insertNextCell(NORM_QUAD4,[4,7,8,5])
        f1=MEDCouplingFieldDouble(ON_GAUSS_PT) ; f1.setName(fieldName)
        f1.setMesh(m1)
        f1.setGaussLocalizationOnType(NORM_TRI3,refCoo1,gsCoo1,wg1)
        f1.setGaussLocalizationOnType(NORM_QUAD4,refCoo2,gsCoo2,wg2)
        arr=DataArrayDouble(f1.getNumberOfTuplesExpected())
        arr.iota()
        f1.setArray(arr)
        f1.checkConsistencyLight()
        #
        m2=MEDCouplingUMesh(mName,2)
        m2.allocateCells() ; m2.setCoords(coo)
        m2.insertNextCell(NORM_QUAD4,[0,3,4,1])
        m2.insertNextCell(NORM_QUAD4,[3,6,7,4])
        ###################
        self.assertTrue(f1.getMesh().getCoords().isEqual(m2.getCoords(),1e-12))
        f1.getMesh().setCoords(m2.getCoords())
        #
        f2=MEDCouplingFieldDouble(ON_GAUSS_PT)
        f2.setMesh(m2)
        for gt in m2.getAllGeoTypes(): # on recopie les localisation en utilisant f1
            glt=f1.getGaussLocalizationIdOfOneType(gt)
            gloc=f1.getGaussLocalization(glt)
            f2.setGaussLocalizationOnType(gt,gloc.getRefCoords(),gloc.getGaussCoords(),gloc.getWeights())
        arr2=DataArrayDouble(f2.getNumberOfTuplesExpected())
        arr2[:]=0
        f2.setArray(arr2)
        f2.checkConsistencyLight()
        #
        fout1=MEDCouplingFieldDouble.MergeFields([f1,f2])
        fout2=MEDCouplingFieldDouble.MergeFields(f1,f2)
        #
        fOut=MEDCouplingFieldDouble(ON_GAUSS_PT)
        mOut=MEDCouplingUMesh.MergeUMeshes([f1.getMesh(),m2])
        mOut.setName(f1.getMesh().getName())
        fOut.setMesh(mOut)
        for gt in f1.getMesh().getAllGeoTypes(): # on recopie les localisation en utilisant f1
            glt=f1.getGaussLocalizationIdOfOneType(gt)
            gloc=f1.getGaussLocalization(glt)
            fOut.setGaussLocalizationOnType(gt,gloc.getRefCoords(),gloc.getGaussCoords(),gloc.getWeights())
        fOut.setArray(DataArrayDouble.Aggregate([f1.getArray(),arr2]))
        fOut.checkConsistencyLight()
        fOut.setName(f1.getName())
        fOut.getMesh().setName(f1.getMesh().getName())
        #
        self.assertTrue(fout1.isEqual(fOut,1e-12,1e-12))
        self.assertTrue(fout2.isEqual(fOut,1e-12,1e-12))
        pass

    pass

if __name__ == '__main__':
    unittest.main()
