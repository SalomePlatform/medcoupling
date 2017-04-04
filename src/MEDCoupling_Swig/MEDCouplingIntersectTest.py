#  -*- coding: utf-8 -*-
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

from MEDCoupling import *
import unittest
from math import pi,e,sqrt,cos,sin
from datetime import datetime
from MEDCouplingDataForTest import MEDCouplingDataForTest
import rlcompleter,readline # this line has to be here, to ensure a usability of MEDCoupling/MEDLoader. B4 removing it please notify to anthony.geay@cea.fr

class MEDCouplingIntersectTest(unittest.TestCase):
    def testSwig2NonRegressionBugIntersectMeshes1(self):
        src=MEDCouplingUMesh("src",2)
        src.setCoords(DataArrayDouble([-2.5,-3,-2.5,3,2.5,3],3,2))
        src.allocateCells()
        src.insertNextCell(NORM_TRI3,[0,1,2])
        #
        trg=MEDCouplingUMesh("trg",2)
        trg.setCoords(DataArrayDouble([-2.5,-3.,0.,-3.,0.,-2.,-2.,0.,-2.25,0.,-2.5,0.,-2.5,-1.5,0.,-2.5,-1.25,-3.,-1.414213562373095,-1.414213562373095],10,2))
        trg.allocateCells()
        trg.insertNextCell(NORM_QPOLYG,[2,1,0,5,3,7,8,6,4,9])
        #
        a,b,c=MEDCouplingUMesh.Intersect2DMeshes(src,trg,1.0e-8)
        a.mergeNodes(1e-8)
        self.assertTrue(a.getCoords().isEqual(DataArrayDouble([-2.5,-3.,-2.5,3.,2.5,3.,0.,-3.,0.,-2.,-2.,0.,-2.25,0.,-2.5,0.,-2.5,-1.5,0.,-2.5,-1.25,-3.,-1.414213562373095,-1.414213562373095,-1.2803687993289596,-1.5364425591947515,-1.8901843996644798,-2.2682212795973755,-1.81117884244736,-0.8483107924994473,-2.5,1.5,0.,3.,0.6098156003355202,0.7317787204026243],18,2),1e-12))
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([32,12,0,7,5,13,8,6,14,32,7,1,2,12,5,15,16,17,14,6])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,9,20])))
        self.assertTrue(b.isEqual(DataArrayInt([0,0])))
        self.assertTrue(c.isEqual(DataArrayInt([0,-1])))
        pass

    def testIntersect2DMeshesTmp1(self):
        m1c=MEDCouplingCMesh.New();
        coordsX=DataArrayDouble.New();
        arrX=[ -1., 1., 2., 4. ]
        coordsX.setValues(arrX,4,1);
        m1c.setCoordsAt(0,coordsX);
        coordsY=DataArrayDouble.New();
        arrY=[ -2., 2., 4., 8. ]
        coordsY.setValues(arrY,4,1);
        m1c.setCoordsAt(1,coordsY);
        m1=m1c.buildUnstructured()
        m1bis=m1.buildPartOfMySelf([3,4,5],False)
        m2=m1.deepCopy()
        m2=m2.buildPartOfMySelf([0,1,2],False)
        m2.translate([0.5,0.5])
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1bis,m2,1e-10)
        expected1=[0,0,1,1,1,2,2,2]
        expected2=[0,-1,0,1,-1,1,2,-1]
        self.assertEqual(8,d1.getNumberOfTuples());
        self.assertEqual(8,d2.getNumberOfTuples());
        self.assertEqual(8,m3.getNumberOfCells());
        self.assertEqual(22,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[5,17,1,16,12,5,16,0,4,5,17,12,5,18,1,17,13,5,19,2,18,13,5,17,5,6,19,13,5,20,2,19,14,5,21,3,20,14,5,19,6,7,21,14]
        expected4=[0,5,12,17,22,28,33,38,44]
        expected5=[-1.0,2.0,1.0,2.0,2.0,2.0,4.0,2.0,-1.0,4.0,1.0,4.0,2.0,4.0,4.0,4.0,-0.5,-1.5,1.5,-1.5,2.5,-1.5,4.5,-1.5,-0.5,2.5,1.5,2.5,2.5,2.5,4.5,2.5,-0.5,2.0,1.0,2.5,1.5,2.0,2.0,2.5,2.5,2.0,4.0,2.5]
        self.assertEqual(44,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(9,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in xrange(44):
            self.assertAlmostEqual(expected5[i],m3.getCoords().getIJ(0,i),12);
            pass
        pass

    def testIntersect2DMeshesTmp2(self):
        m1c=MEDCouplingCMesh.New();
        coordsX1=DataArrayDouble.New();
        arrX1=[ 0., 1., 1.5, 2. ]
        coordsX1.setValues(arrX1,4,1);
        m1c.setCoordsAt(0,coordsX1);
        coordsY1=DataArrayDouble.New();
        arrY1=[ 0., 1.5, 3.]
        coordsY1.setValues(arrY1,3,1);
        m1c.setCoordsAt(1,coordsY1);
        m1=m1c.buildUnstructured();
        m2c=MEDCouplingCMesh.New();
        coordsX2=DataArrayDouble.New();
        arrX2=[ 0., 1., 2. ]
        coordsX2.setValues(arrX2,3,1);
        m2c.setCoordsAt(0,coordsX2);
        coordsY2=DataArrayDouble.New();
        arrY2=[ 0., 1., 3.]
        coordsY2.setValues(arrY2,3,1);
        m2c.setCoordsAt(1,coordsY2);
        m2=m2c.buildUnstructured();
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1,m2,1e-10)
        #
        expected1=[0,0,1,1,2,2,3,4,5]
        expected2=[0,2,1,3,1,3,2,3,3]
        self.assertEqual(9,d1.getNumberOfTuples());
        self.assertEqual(9,d2.getNumberOfTuples());
        self.assertEqual(9,m3.getNumberOfCells());
        self.assertEqual(22,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[5,16,13,12,15,5,15,4,5,16,5,21,2,13,16,5,16,5,6,21,5,17,14,2,21,5,21,6,7,17,5,4,18,19,5,5,5,19,10,6,5,6,10,20,7]
        expected4=[0,5,10,15,20,25,30,35,40,45]
        expected5=[0.0,0.0,1.0,0.0,1.5,0.0,2.0,0.0,0.0,1.5,1.0,1.5,1.5,1.5,2.0,1.5,0.0,3.0,1.0,3.0,1.5,3.0,2.0,3.0,0.0,0.0,1.0,0.0,2.0,0.0,0.0,1.0,1.0,1.0,2.0,1.0,0.0,3.0,1.0,3.0,2.0,3.0,1.5,1.0]
        self.assertEqual(45,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(10,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in xrange(44):
            self.assertAlmostEqual(expected5[i],m3.getCoords().getIJ(0,i),12);
            pass
        pass

    def testIntersect2DMeshesTmp3(self):
        m1Coords=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1,0.,-1.5,0.5,0.,1.25,0.,0.70710678118654757,0.70710678118654757,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.70710678118654757,0.70710678118654757,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.70710678118654757,-0.70710678118654757,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.70710678118654757,-0.70710678118654757,1.0606601717798214,-1.0606601717798214];
        m1Conn=[0,3,1,13,11,9, 3,4,2,1,14,12,10,11, 5,3,0,15,13,17, 6,4,3,5,16,14,15,18, 5,0,7,17,21,19, 6,5,7,8,18,19,22,20, 0,1,7,9,23,21, 1,2,8,7,10,24,22,23];
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(2);
        m1.allocateCells(8);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[0:6]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[6:14]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[14:20]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[20:28]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[28:34]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[34:42]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[42:48]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[48:56]);
        m1.finishInsertingCells();
        myCoords1=DataArrayDouble.New();
        myCoords1.setValues(m1Coords,25,2);
        m1.setCoords(myCoords1);
        #
        m2Coords=[0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1,-1.1,-1.,0.,-1.,1.1,-1,1.7,-1.]
        m2Conn=[0,3,2,1, 1,2,5,4, 7,6,3,0, 8,9,6,7, 7,0,12,11, 8,7,11,10, 0,1,13,12, 1,4,14,13]
        m2=MEDCouplingUMesh.New();
        m2.setMeshDimension(2);
        m2.allocateCells(8);
        for i in xrange(8):
            m2.insertNextCell(NORM_QUAD4,4,m2Conn[4*i:4*(i+1)])
            pass
        m2.finishInsertingCells();
        myCoords2=DataArrayDouble.New();
        myCoords2.setValues(m2Coords,15,2);
        m2.setCoords(myCoords2);
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1,m2,1e-10)
        m3.unPolyze()
        #
        expected1=[0,1,1,1,2,3,3,3,4,5,5,5,6,7,7,7]
        expected2=[0,0,1,-1,2,2,3,-1,4,4,5,-1,6,6,7,-1]
        self.assertEqual(16,d1.getNumberOfTuples());
        self.assertEqual(16,d2.getNumberOfTuples());
        self.assertEqual(16,m3.getNumberOfCells());
        self.assertEqual(104,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[6,28,1,25,44,45,46,8,26,1,28,27,47,48,49,50,8,40,2,26,27,51,52,53,54,8,28,4,40,27,55,56,57,58,6,28,25,5,59,60,61,8,28,5,32,31,62,63,64,65,8,32,6,41,31,66,67,68,69,8,41,4,28,31,70,71,72,73,6,25,37,5,74,75,76,8,32,5,37,36,77,78,79,80,8,42,6,32,36,81,82,83,84,8,37,8,42,36,85,86,87,88,6,1,37,25,89,90,91,8,37,1,26,38,92,93,94,95,8,26,2,43,38,96,97,98,99,8,43,8,37,38,100,101,102,103]
        expected4=[0,7,16,25,34,41,50,59,68,75,84,93,102,109,118,127,136]
        expected5=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1.,-1.1,-1.,0.,-1.,1.1,-1.,1.7,-1.,1.118033988749895,1.,-1.118033988749895,1.,-1.118033988749895,-1.,1.118033988749895,-1.,0.7071067811865477,0.7071067811865476,0.5,0.,0.,0.5,1.05,0.,0.7071067811865475,0.7071067811865477,0.55,1.,1.1,0.5,1.4012585384440737,0.535233134659635,1.3,0.,1.1,0.5,1.1090169943749475,1.,0.,1.25,0.6123724356957946,1.369306393762915,1.1090169943749475,1.,0.55,1.,0.,0.5,-0.5,0.,-0.7071067811865477,0.7071067811865476,-0.7071067811865475,0.7071067811865477,-1.05,0.,-1.1,0.5,-0.55,1.,-1.3,0.,-1.4012585384440737,0.5352331346596344,-1.1090169943749475,1.,-1.1,0.5,-0.6123724356957941,1.3693063937629155,0.,1.25,-0.55,1.,-1.1090169943749475,1.,0.,-0.5,-0.7071067811865475,-0.7071067811865477,-0.5,0.,-1.05,0.,-0.7071067811865478,-0.7071067811865475,-0.55,-1.,-1.1,-0.5,-1.4012585384440734,-0.5352331346596354,-1.3,0.,-1.1,-0.5,-1.1090169943749475,-1.,0.,-1.25,-0.6123724356957945,-1.369306393762915,-1.1090169943749475,-1.,-0.55,-1.,0.7071067811865475,-0.7071067811865477,0.,-0.5,0.5,0.,0.7071067811865477,-0.7071067811865475,1.05,0.,1.1,-0.5,0.55,-1.,1.3,0.,1.4012585384440737,-0.535233134659635,1.1090169943749475,-1.,1.1,-0.5,0.6123724356957946,-1.369306393762915,0.,-1.25,0.55,-1.,1.1090169943749475,-1.0]
        self.assertEqual(136,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(17,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in xrange(208):
            self.assertAlmostEqual(expected5[i],m3.getCoords().getIJ(0,i),12);
            pass
        pass

    def testIntersect2DMeshesTmp4(self):
        m1Coords=[0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1,0.,-1.5,0.5,0.,1.25,0.,0.70710678118654757,0.70710678118654757,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.70710678118654757,0.70710678118654757,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.70710678118654757,-0.70710678118654757,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.70710678118654757,-0.70710678118654757,1.0606601717798214,-1.0606601717798214];
        m1Conn=[0,3,1,13,11,9, 3,4,2,1,14,12,10,11, 5,3,0,15,13,17, 6,4,3,5,16,14,15,18, 5,0,7,17,21,19, 6,5,7,8,18,19,22,20, 0,1,7,9,23,21, 1,2,8,7,10,24,22,23];
        m1=MEDCouplingUMesh.New();
        m1.setMeshDimension(2);
        m1.allocateCells(8);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[0:6]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[6:14]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[14:20]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[20:28]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[28:34]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[34:42]);
        m1.insertNextCell(NORM_TRI6,6,m1Conn[42:48]);
        m1.insertNextCell(NORM_QUAD8,8,m1Conn[48:56]);
        m1.finishInsertingCells();
        myCoords1=DataArrayDouble.New();
        myCoords1.setValues(m1Coords,25,2);
        m1.setCoords(myCoords1);
        #
        m2Coords=[0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1,-1.1,-1.,0.,-1.,1.1,-1,1.7,-1.]
        m2Conn=[0,3,2,1, 1,2,5,4, 7,6,3,0, 8,9,6,7, 7,0,12,11, 8,7,11,10, 0,1,13,12, 1,4,14,13]
        m2=MEDCouplingUMesh.New();
        m2.setMeshDimension(2);
        m2.allocateCells(8);
        for i in xrange(8):
            m2.insertNextCell(NORM_QUAD4,4,m2Conn[4*i:4*(i+1)])
            pass
        m2.finishInsertingCells();
        myCoords2=DataArrayDouble.New();
        myCoords2.setValues(m2Coords,15,2);
        m2.setCoords(myCoords2);
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m2,m1,1e-10)
        m3.unPolyze()
        #
        expected1=[0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7]
        expected2=[0,1,1,-1,2,3,3,-1,4,5,5,-1,6,7,7,-1]
        self.assertEqual(16,d1.getNumberOfTuples());
        self.assertEqual(16,d2.getNumberOfTuples());
        self.assertEqual(16,m3.getNumberOfCells());
        self.assertEqual(104,m3.getNumberOfNodes());
        self.assertEqual(2,m3.getSpaceDimension());
        self.assertEqual(expected1,d1.getValues());
        self.assertEqual(expected2,d2.getValues());
        expected3=[6,16,15,18,44,45,46,8,18,2,1,16,47,48,49,50,8,17,1,2,40,51,52,53,54,8,40,5,4,17,55,56,57,58,6,18,15,20,59,60,61,8,20,7,6,18,62,63,64,65,8,41,6,7,21,66,67,68,69,8,21,8,9,41,70,71,72,73,6,20,15,22,74,75,76,8,22,11,7,20,77,78,79,80,8,21,7,11,42,81,82,83,84,8,42,10,8,21,85,86,87,88,6,22,15,16,89,90,91,8,16,1,13,22,92,93,94,95,8,43,13,1,17,96,97,98,99,8,17,4,14,43,100,101,102,103]
        expected4=[0,7,16,25,34,41,50,59,68,75,84,93,102,109,118,127,136]
        expected5=[0.,0.,1.1, 0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1.,-1.1,-1.,0.,-1.,1.1,-1.,1.7,-1.,0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,1.1180339887498951,1.,-1.1180339887498951,1.,-1.1180339887498951,-1.,1.1180339887498951,-1.,0.5,0.,0.,0.5,0.7071067811865477,0.7071067811865476,0.55,1.,1.1,0.5,1.05,0.,0.7071067811865477,0.7071067811865475,1.3,0.,1.1,0.5,1.1090169943749475,1.,1.4012585384440737,0.535233134659635,1.4090169943749475,1.,1.7,0.5,1.6,0.,1.4012585384440737,0.535233134659635,0.,0.5,-0.5,0.,-0.7071067811865477,0.7071067811865476,-1.05,0.,-1.1,0.5,-0.55,1.,-0.7071067811865478,0.7071067811865475,-1.1090169943749475,1.,-1.1,0.5,-1.3,0.,-1.4012585384440737,0.5352331346596344,-1.6,0.,-1.7,0.5,-1.4090169943749475,1.,-1.4012585384440737,0.5352331346596344,-0.5,0.,0.,-0.5,-0.7071067811865475,-0.7071067811865477,-0.55,-1.,-1.1,-0.5,-1.05,0.,-0.7071067811865475,-0.7071067811865477,-1.3,0.,-1.1,-0.5,-1.1090169943749475,-1.,-1.4012585384440734,-0.5352331346596354,-1.4090169943749475,-1.,-1.7,-0.5,-1.6,0.,-1.4012585384440732,-0.5352331346596354,0.,-0.5,0.5,0.,0.7071067811865475,-0.7071067811865477,1.05,0.,1.1,-0.5,0.55,-1.,0.7071067811865475,-0.7071067811865477,1.1090169943749475,-1.,1.1,-0.5,1.3,0.,1.4012585384440737,-0.535233134659635,1.6,0.,1.7,-0.5,1.4090169943749475,-1.,1.4012585384440737,-0.535233134659635]
        self.assertEqual(136,m3.getNodalConnectivity().getNumberOfTuples());
        self.assertEqual(17,m3.getNodalConnectivityIndex().getNumberOfTuples());
        self.assertEqual(expected3,m3.getNodalConnectivity().getValues());
        self.assertEqual(expected4,m3.getNodalConnectivityIndex().getValues());
        for i in xrange(208):
            self.assertAlmostEqual(expected5[i],m3.getCoords().getIJ(0,i),12);
            pass
        pass

    def testIntersect2DMeshesTmp5(self):
        coords=DataArrayDouble.New([41,0,42,0,0,42,0,41,41.5,0,29.698484809834998,29.698484809834994,0,41.5,28.991378028648452,28.991378028648445,-42,0,-41,0,-29.698484809834994,29.698484809834998,-41.5,0,-28.991378028648445,28.991378028648452,0,-42,0,-41,-29.698484809835001,-29.698484809834994,0,-41.5,-28.991378028648455,-28.991378028648445,29.698484809834987,-29.698484809835001,28.991378028648441,-28.991378028648455,43,0,0,43,42.5,0,30.405591591021544,30.40559159102154,0,42.5,-43,0,-30.40559159102154,30.405591591021544,-42.5,0,0,-43,-30.405591591021551,-30.40559159102154,0,-42.5,30.405591591021537,-30.405591591021551,44,0,0,44,43.5,0,31.112698372208094,31.112698372208087,0,43.5,-44,0,-31.112698372208087,31.112698372208094,-43.5,0,0,-44,-31.112698372208097,-31.112698372208087,0,-43.5,31.112698372208083,-31.112698372208097,45,0,0,45,44.5,0,31.81980515339464,31.819805153394636,0,44.5,-45,0,-31.819805153394636,31.81980515339464,-44.5,0,0,-45,-31.819805153394647,-31.819805153394636,0,-44.5,31.819805153394629,-31.819805153394647,47,0,0,47,46,0,33.234018715767739,33.234018715767732,0,46,-47,0,-33.234018715767732,33.234018715767739,-46,0,0,-47,-33.234018715767739,-33.234018715767732,0,-46,33.234018715767725,-33.234018715767739,49,0,0,49,48,0,34.648232278140831,34.648232278140824,0,48,-49,0,-34.648232278140824,34.648232278140831,-48,0,0,-49,-34.648232278140839,-34.648232278140824,0,-48,34.648232278140817,-34.648232278140839,51,0,0,51,50,0,36.062445840513924,36.062445840513924,0,50,-51,0,-36.062445840513924,36.062445840513924,-50,0,0,-51,-36.062445840513931,-36.062445840513924,0,-50,36.062445840513917,-36.062445840513931,53,0,0,53,52,0,37.476659402887023,37.476659402887016,0,52,-53,0,-37.476659402887016,37.476659402887023,-52,0,0,-53,-37.47665940288703,-37.476659402887016,0,-52,37.476659402887009,-37.47665940288703,55,0,0,55,54,0,38.890872965260115,38.890872965260108,0,54,-55,0,-38.890872965260108,38.890872965260115,-54,0,0,-55,-38.890872965260122,-38.890872965260108,0,-54,38.890872965260101,-38.890872965260122,59,0,0,59,57,0,41.719300090006307,41.7193000900063,0,57,-59,0,-41.7193000900063,41.719300090006307,-57,0,0,-59,-41.719300090006314,-41.7193000900063,0,-57,41.719300090006293,-41.719300090006314,63,0,0,63,61,0,44.547727214752499,44.547727214752491,0,61,-63,0,-44.547727214752491,44.547727214752499,-61,0,0,-63,-44.547727214752506,-44.547727214752491,0,-61,44.547727214752484,-44.547727214752506,67,0,0,67,65,0,47.37615433949869,47.376154339498683,0,65,-67,0,-47.376154339498683,47.37615433949869,-65,0,0,-67,-47.376154339498697,-47.376154339498683,0,-65,47.376154339498676,-47.376154339498697,71,0,0,71,69,0,50.204581464244875,50.204581464244868,0,69,-71,0,-50.204581464244868,50.204581464244875,-69,0,0,-71,-50.204581464244889,-50.204581464244868,0,-69,50.20458146424486,-50.204581464244889,75,0,0,75,73,0,53.033008588991066,53.033008588991059,0,73,-75,0,-53.033008588991059,53.033008588991066,-73,0,0,-75,-53.033008588991073,-53.033008588991059,0,-73,53.033008588991052,-53.033008588991073,80,0,0,80,77.5,0,56.568542494923804,56.568542494923797,0,77.5,-80,0,-56.568542494923797,56.568542494923804,-77.5,0,0,-80,-56.568542494923818,-56.568542494923797,0,-77.5,56.56854249492379,-56.568542494923818],188,2)
        conn=DataArrayInt.New([8,0,1,2,3,4,5,6,7,8,3,2,8,9,6,10,11,12,8,9,8,13,14,11,15,16,17,8,14,13,1,0,16,18,4,19,8,1,20,21,2,22,23,24,5,8,2,21,25,8,24,26,27,10,8,8,25,28,13,27,29,30,15,8,13,28,20,1,30,31,22,18,8,20,32,33,21,34,35,36,23,8,21,33,37,25,36,38,39,26,8,25,37,40,28,39,41,42,29,8,28,40,32,20,42,43,34,31,8,32,44,45,33,46,47,48,35,8,33,45,49,37,48,50,51,38,8,37,49,52,40,51,53,54,41,8,40,52,44,32,54,55,46,43,8,44,56,57,45,58,59,60,47,8,45,57,61,49,60,62,63,50,8,49,61,64,52,63,65,66,53,8,52,64,56,44,66,67,58,55,8,56,68,69,57,70,71,72,59,8,57,69,73,61,72,74,75,62,8,61,73,76,64,75,77,78,65,8,64,76,68,56,78,79,70,67,8,68,80,81,69,82,83,84,71,8,69,81,85,73,84,86,87,74,8,73,85,88,76,87,89,90,77,8,76,88,80,68,90,91,82,79,8,80,92,93,81,94,95,96,83,8,81,93,97,85,96,98,99,86,8,85,97,100,88,99,101,102,89,8,88,100,92,80,102,103,94,91,8,92,104,105,93,106,107,108,95,8,93,105,109,97,108,110,111,98,8,97,109,112,100,111,113,114,101,8,100,112,104,92,114,115,106,103,8,104,116,117,105,118,119,120,107,8,105,117,121,109,120,122,123,110,8,109,121,124,112,123,125,126,113,8,112,124,116,104,126,127,118,115,8,116,128,129,117,130,131,132,119,8,117,129,133,121,132,134,135,122,8,121,133,136,124,135,137,138,125,8,124,136,128,116,138,139,130,127,8,128,140,141,129,142,143,144,131,8,129,141,145,133,144,146,147,134,8,133,145,148,136,147,149,150,137,8,136,148,140,128,150,151,142,139,8,140,152,153,141,154,155,156,143,8,141,153,157,145,156,158,159,146,8,145,157,160,148,159,161,162,149,8,148,160,152,140,162,163,154,151,8,152,164,165,153,166,167,168,155,8,153,165,169,157,168,170,171,158,8,157,169,172,160,171,173,174,161,8,160,172,164,152,174,175,166,163,8,164,176,177,165,178,179,180,167,8,165,177,181,169,180,182,183,170,8,169,181,184,172,183,185,186,173,8,172,184,176,164,186,187,178,175],540)
        connI=DataArrayInt.New([0,9,18,27,36,45,54,63,72,81,90,99,108,117,126,135,144,153,162,171,180,189,198,207,216,225,234,243,252,261,270,279,288,297,306,315,324,333,342,351,360,369,378,387,396,405,414,423,432,441,450,459,468,477,486,495,504,513,522,531,540],61)
        #
        m1=MEDCouplingUMesh.New("Fix",2);
        m1.setCoords(coords);
        m1.setConnectivity(conn,connI,True);
        #
        coords=DataArrayDouble([46.5,-2.5,53.5,-2.5,53.5,2.5,46.5,2.5,50,-2.5,53.5,0,50,2.5,46.5,0,60.5,-2.5,60.5,2.5,57,-2.5,60.5,0,57,2.5,53.5,7.5,46.5,7.5,53.5,5,50,7.5,46.5,5,60.5,7.5,60.5,5,57,7.5,-2,47,2,47,2,53,-2,53,0,47,2,50,0,53,-2,50,6,47,6,53,4,47,6,50,4,53,2,59,-2,59,2,56,0,59,-2,56,6,59,6,56,4,59],42,2)
        # connectivity
        conn=DataArrayInt([8,0,1,2,3,4,5,6,7,8,1,8,9,2,10,11,12,5,8,3,2,13,14,6,15,16,17,8,2,9,18,13,12,19,20,15,8,21,22,23,24,25,26,27,28,8,22,29,30,23,31,32,33,26,8,24,23,34,35,27,36,37,38,8,23,30,39,34,33,40,41,36],72);
        conn.setName("");
        connI=DataArrayInt([0,9,18,27,36,45,54,63,72],9)
        m2=MEDCouplingUMesh.New("Mobile",2);
        m2.setCoords(coords);
        m2.setConnectivity(conn,connI,True);
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1,m2,1e-10);
        self.assertEqual(105,m3.getNumberOfCells());
        self.assertEqual(105,d1.getNumberOfTuples());
        self.assertEqual(105,d2.getNumberOfTuples());
        self.assertEqual(704,m3.getNumberOfNodes());
        #
        areaExpected=[-65.18804756198824,-65.18804756198824,-65.18804756198824,-65.18804756198824,-66.75884388878285,-66.75884388878285,-66.7588438887833,-66.75884388878308,-68.32964021557768,-68.32964021557768,-68.32964021557814,-68.32964021557791,-69.9004365423732,-69.9004365423732,-69.90043654237297,-69.90043654237297,-1.194568659706448,-1.0869994447159463,-142.2316939607081,-144.51326206513068,-144.5132620651309,-1.1945686597064424,-143.3186934054243,-5.002264310862817,-10.0261332846393,-3.9727823117092953,-7.290862524642649,-124.504404940456,-3.9727823117093237,-146.82366506060032,-150.79644737231024,-5.002264310862776,-145.79418306144626,-5.00208651738126,-10.054764051268958,-4.001067863263231,-8.027932154428669,-129.99378209314813,-4.001067863263216,-153.07856481622616,-157.0796326794898,-5.0020865173811915,-152.07754616210832,-5.001928880064381,-10.050590216368969,-4.00098721602491,-8.025810856794209,-136.28350081741684,-4.000987216024939,-159.36183077064402,-163.36281798667005,-5.0019288800643285,-158.36088910660442,-1.2991516319851801,-3.702636830195414,-3.7815130030068254,-6.265364371195623,-0.02516260900254963,-0.6553944641345026,-3.975752765070567,-7.368528340442765,-142.57249927881398,-0.02516260900254963,-3.9757527650706095,-165.64508791977525,-169.64600329384803,-1.299151631985167,-3.7026368301953885,-164.6442148316677,-10.00321285677458,-20.08414323176165,-8.001644468035863,-16.042954878437143,-304.0096070742277,-8.00164446803587,-350.1399180412005,-358.1415625092368,-10.003212856774468,-348.13834965246224,-3.794150313030109,-8.65049239704272,-0.02260276689354157,-0.5885167811200915,-370.2185414798688,-0.022602766893559393,-383.2517009710623,-383.2743037379555,-3.7941503130300576,-379.48015342492505,-408.40704496667513,-408.4070449666742,-408.4070449666742,-408.4070449666742,-433.53978619538975,-433.5397861953902,-433.5397861953911,-433.53978619539066,-458.67252742410983,-458.6725274241094,-458.67252742410983,-458.6725274241089,-608.6835766330232,-608.6835766330232,-608.6835766330232,-608.6835766330241]
        expected1=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,16,16,17,18,19,19,20,20,20,20,20,21,21,22,23,23,24,24,24,24,24,25,25,26,27,27,28,28,28,28,28,29,29,30,31,31,32,32,32,32,32,32,32,32,32,33,33,33,34,35,35,35,36,36,36,36,36,37,37,38,39,39,40,40,40,40,40,41,41,42,43,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59]
        expected2=[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,2,-1,-1,-1,0,-1,0,2,4,5,-1,4,-1,-1,0,-1,0,2,4,5,-1,4,-1,-1,0,-1,0,2,4,5,-1,4,-1,-1,0,-1,0,1,2,3,4,5,6,7,-1,4,6,-1,-1,0,1,-1,1,3,6,7,-1,6,-1,-1,1,-1,1,3,6,7,-1,6,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]
        f3=m3.getMeasureField(False).getArray().getValues();
        for i in xrange(105):
            self.assertAlmostEqual(areaExpected[i],f3[i],10)
            pass
        self.assertEqual(expected1,d1.getValues())
        self.assertEqual(expected2,d2.getValues())
        pass

    def testIntersect2DMeshesTmp6(self):
        # coordinates
        coords=DataArrayDouble.New([2.7554552980815448e-15,45,-45,5.5109105961630896e-15,-31.819805153394636,31.81980515339464,2.8779199779962799e-15,47,2.8166876380389124e-15,46,-47,5.7558399559925599e-15,-33.234018715767732,33.234018715767739,-46,5.6333752760778247e-15],8,2);
        # connectivity
        conn=DataArrayInt.New([8,0,3,5,1,4,6,7,2])
        connI=DataArrayInt.New([0,9]);
        m1=MEDCouplingUMesh.New("Fixe",2);
        m1.setCoords(coords);
        m1.setConnectivity(conn,connI,True);
        #
        coords=DataArrayDouble.New([-7.3800475508445391,41.854329503018846,-3.7041190667754655,42.338274668899189,-3.7041190667754655,45.338274668899189,-7.3800475508445382,44.854329503018839,-5.5473631693521845,42.136406608386956,-3.7041190667754655,43.838274668899189,-5.5420833088100014,45.09630208595901,-7.3800475508445382,43.354329503018839,-3.7041190667754651,52.338274668899189,-7.3800475508445382,51.854329503018839,-3.7041190667754655,48.838274668899189,-5.5420833088100014,52.09630208595901,-7.3800475508445382,48.354329503018839],13,2);
        # connectivity
        conn=DataArrayInt.New([8,0,1,2,3,4,5,6,7,8,3,2,8,9,6,10,11,12]);
        connI=DataArrayInt.New([0,9,18]);
        #
        m2=MEDCouplingUMesh.New("Mobile",2);
        m2.setCoords(coords);
        m2.setConnectivity(conn,connI,True);
        #
        m3,d1,d2=MEDCouplingUMesh.Intersect2DMeshes(m1,m2,1e-10);
        self.assertTrue(d1.isEqual(DataArrayInt([0,0,0,0])));
        self.assertTrue(d2.isEqual(DataArrayInt([0,1,-1,-1])));
        self.assertEqual(4,m3.getNumberOfCells());
        self.assertEqual(4,d1.getNumberOfTuples());
        self.assertEqual(4,d2.getNumberOfTuples());
        self.assertEqual(43,m3.getNumberOfNodes());
        dI,areMerged,newNbOfNodes=m3.mergeNodes(1e-12)
        self.assertEqual(35,m3.getNumberOfNodes());
        m3.zipCoords();
        self.assertEqual(23,m3.getNumberOfNodes());
        #
        f=m3.getMeasureField(True);
        valuesExpected=DataArrayDouble([1.6603638692585716,5.747555728471923,129.68907101754394,7.4162714498559694])
        self.assertTrue(f.getArray().isEqual(valuesExpected,1e-12))
        pass


    def testSwig2Intersect2DMeshesQuadra1(self):
        import cmath
        def createDiagCircle(lX, lY, R, cells=[0,1]):
            """ A circle in a square box, cut along the diagonal.
            """
            c = []
            for i in range(8):
              c.append(cmath.rect(R, i*pi/4))

            coords = [0.0,0.0,          c[3].real,c[3].imag,       -lX/2.0, lY/2.0,
                      0.0, lY/2.0,      lX/2.0,lY/2.0,             lX/2.0,0.0,
                      #   6                  7                              8
                      lX/2.0,-lY/2.0,   c[7].real,c[7].imag,       c[1].real,c[1].imag,
                      #   9                  10                            11
                      c[5].real,c[5].imag,   -lX/2.0,-lY/2.0,      0.0, -lY/2.0,
                      #   12                  13                            14
                      -lX/2.0,0.0,         0.0,0.0,                  0.0, 0.0]
            # Points 13 (reps. 14) are average of points (6,7) (resp (1,2))
            coords[13*2]   = 0.5*(coords[6*2]+coords[7*2])
            coords[13*2+1] = 0.5*(coords[6*2+1]+coords[7*2+1])
            coords[14*2]   = 0.5*(coords[1*2]+coords[2*2])
            coords[14*2+1] = 0.5*(coords[1*2+1]+coords[2*2+1])
            connec  = [1,7,8,0]      # half circle up right
            connec3 = [6,7,1,2,4,13,8,14,3,5]

            baseMesh = MEDCouplingUMesh.New("box_circle", 2)
            baseMesh.allocateCells(2)
            meshCoords = DataArrayDouble.New(coords, len(coords)/2, 2)
            meshCoords.setInfoOnComponents(["X [au]", "Y [au]"])
            baseMesh.setCoords(meshCoords)

            if 0 in cells:
              baseMesh.insertNextCell(NORM_QPOLYG, connec)
            if 1 in cells:
              baseMesh.insertNextCell(NORM_QPOLYG, connec3)
            baseMesh.finishInsertingCells()
            baseMesh.checkConsistencyLight()
            return baseMesh

        eps = 1.0e-7
        m1 = createDiagCircle(1.0, 1.0, 0.5*0.90, cells=[0,1])
        m2 = createDiagCircle(1.0, 1.0, 0.5*0.95, cells=[0])
        m3, _, _= MEDCouplingUMesh.Intersect2DMeshes(m1, m2, eps)
        m3.mergeNodes(eps)
        m3.convertDegeneratedCells()
        m3.zipCoords()
        m4 = m3.deepCopy()
        m5, _, _ = MEDCouplingUMesh.Intersect2DMeshes(m3, m4, eps)
        m5.mergeNodes(eps)
        # Check coordinates:
        self.assertTrue(m3.getCoords().isEqual(m5.getCoords(), eps))

    def testIntersect2DMeshesTmp7(self):
        eps = 1.0e-8
        coords = [-0.5,-0.5,   -0.5, 0.5, 0.5, 0.5,    0.5,-0.5]
        connec = range(4)
        m1 = MEDCouplingUMesh.New("box", 2)
        m1.allocateCells(1)
        meshCoords = DataArrayDouble.New(coords, len(coords)/2, 2)
        m1.setCoords(meshCoords)
        m1.insertNextCell(NORM_POLYGON, connec)
        m1.finishInsertingCells()

        m2 = MEDCouplingDataForTest.buildCircle(0.25, 0.2, 0.4)
        # Was looping indefinitly:
        m_intersec, resToM1, resToM2 = MEDCouplingUMesh.Intersect2DMeshes(m1, m2, eps)
        m_intersec.zipCoords()
        coo_tgt = DataArrayDouble([-0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, -0.5, -0.03284271247461901, 0.4828427124746191,
          -0.014575131106459124, 0.5000000000000001, 0.5, -0.11224989991991996, 0.24271243444677046, 0.5, 0.5, 0.19387505004004,
          -0.04799910280454185, -0.06682678787499614, -0.023843325638122054, 0.4915644577163915, 0.5, -0.30612494995996, 0.0, -0.5,
          -0.5, 0.0, -0.25728756555322957, 0.5, -0.023843325638122026, 0.49156445771639157, -0.04799910280454181, -0.06682678787499613], 17 ,2)
        conn_tgt = [32, 5, 2, 6, 4, 7, 8, 9, 10, 32, 6, 3, 0, 1, 5, 4, 11, 12, 13, 14, 15, 16]
        connI_tgt = [0, 9, 22]
        res1_tgt  = [0, 0]
        res2_tgt = [0, -1]
        self.assert_(coo_tgt.isEqualWithoutConsideringStr(m_intersec.getCoords(), 1e-12))
        self.assertEqual(conn_tgt, m_intersec.getNodalConnectivity().getValues())
        self.assertEqual(connI_tgt, m_intersec.getNodalConnectivityIndex().getValues())
        self.assertEqual(res1_tgt, resToM1.getValues())
        self.assertEqual(res2_tgt, resToM2.getValues())

    def testSwig2Intersect2DMeshWith1DLine1(self):
        """A basic test with no colinearity between m1 and m2."""
        i=MEDCouplingIMesh("mesh",2,[5,5],[0.,0.],[1.,1.])
        m1=i.buildUnstructured()
        m2=MEDCouplingUMesh("mesh",1) ; m2.setCoords(DataArrayDouble([0.75,3.5,3.75,1.75],2,2)) ; m2.allocateCells() ; m2.insertNextCell(NORM_SEG2,[0,1])
        a,b,c,d=MEDCouplingUMesh.Intersect2DMeshWith1DLine(m1,m2,1e-12)
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([4,1,0,5,6,4,2,1,6,7,4,3,2,7,8,4,4,3,8,9,4,6,5,10,11,4,7,6,11,12,4,8,7,12,13,4,11,10,15,16,4,18,17,22,23,4,19,18,23,24,5,16,15,20,21,31,5,21,22,17,28,31,5,16,31,28,5,17,29,28,5,12,11,16,28,29,5,17,18,30,29,5,13,12,29,30,5,18,19,14,27,30,5,13,30,27,5,9,8,13,27,14])))
        self.assertTrue(b.getNodalConnectivity().isEqual(DataArrayInt([1,25,31,1,31,28,1,28,29,1,29,30,1,30,27,1,27,26])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25,30,35,40,45,50,56,62,66,70,76,81,86,92,96,102])))
        self.assertTrue(b.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6,9,12,15,18])))
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:25].isEqual(m1.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[25:25+2].isEqualWithoutConsideringStr(m2.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[27:].isEqualWithoutConsideringStr(DataArrayDouble([(3.3214285714285716,2.),(1.6071428571428572,3.),(2.,2.7708333333333335),(3.,2.1875),(1.,3.354166666666667)]),1e-12))
        self.assertTrue(c.isEqual(DataArrayInt([0,1,2,3,4,5,6,8,14,15,12,13,13,9,9,10,10,11,11,7])))
        self.assertTrue(d.isEqual(DataArrayInt([(10,10),(11,12),(13,14),(15,16),(17,18),(19,19)])))
        pass

    def testSwig2Intersect2DMeshWith1DLine2(self):
        """A basic test with colinearity between m1 and m2 and the last cell of m2 outside m1."""
        i=MEDCouplingIMesh("mesh",2,[5,5],[0.,0.],[1.,1.])
        m1=i.buildUnstructured()
        m2=MEDCouplingUMesh("mesh",1) ; m2.setCoords(DataArrayDouble([0.5,2.,2.25,2.,2.5,2.,2.75,2.,3.,2.,4.,2.,5.,2.],7,2)) ; m2.allocateCells()
        for i in xrange(6):
            m2.insertNextCell(NORM_SEG2,[i,i+1])
            pass
        a,b,c,d=MEDCouplingUMesh.Intersect2DMeshWith1DLine(m1,m2,1e-12)
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([4,1,0,5,6,4,2,1,6,7,4,3,2,7,8,4,4,3,8,9,4,16,15,20,21,4,17,16,21,22,4,18,17,22,23,4,19,18,23,24,5,6,5,10,25,11,5,7,6,11,12,5,8,7,12,26,27,28,13,5,9,8,13,14,5,11,25,10,15,16,5,12,11,16,17,5,13,28,27,26,12,17,18,5,14,13,18,19])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25,30,35,40,46,51,59,64,70,75,83,88])))
        self.assertTrue(b.getNodalConnectivity().isEqual(DataArrayInt([1,25,11,1,11,12,1,12,26,1,26,27,1,27,28,1,28,13,1,13,14,1,14,31])))
        self.assertTrue(b.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6,9,12,15,18,21,24])))
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:25].isEqual(m1.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[25:].isEqualWithoutConsideringStr(m2.getCoords(),1e-12))
        self.assertTrue(c.isEqual(DataArrayInt([0,1,2,3,12,13,14,15,4,5,6,7,8,9,10,11])))
        self.assertTrue(d.isEqual(DataArrayInt([(12,8),(13,9),(14,10),(14,10),(14,10),(14,10),(15,11),(-1,-1)])))
        pass

    def testSwig2Intersect2DMeshWith1DLine3(self):
        """m2 fully included in cell #12. of m1"""
        i=MEDCouplingIMesh("mesh",2,[5,5],[0.,0.],[1.,1.])
        m1=i.buildUnstructured()
        m2=MEDCouplingUMesh("mesh",1) ; m2.setCoords(DataArrayDouble([(0.75,3.25),(0.5,3.5),(0.25,3.25)])) ; m2.allocateCells()
        for i in xrange(2):
            m2.insertNextCell(NORM_SEG2,[i,i+1])
            pass
        a,b,c,d=MEDCouplingUMesh.Intersect2DMeshWith1DLine(m1,m2,1e-12)
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([4,1,0,5,6,4,2,1,6,7,4,3,2,7,8,4,4,3,8,9,4,6,5,10,11,4,7,6,11,12,4,8,7,12,13,4,9,8,13,14,4,11,10,15,16,4,12,11,16,17,4,13,12,17,18,4,14,13,18,19,4,17,16,21,22,4,18,17,22,23,4,19,18,23,24,5,16,15,20,21])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80])))
        self.assertTrue(b.getNodalConnectivity().isEqual(DataArrayInt([1,25,26,1,26,27])))
        self.assertTrue(b.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6])))
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:25].isEqual(m1.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[25:].isEqualWithoutConsideringStr(m2.getCoords(),1e-12))
        self.assertTrue(c.isEqual(DataArrayInt([0,1,2,3,4,5,6,7,8,9,10,11,13,14,15,12])))
        self.assertTrue(d.isEqual(DataArrayInt([(15,15),(15,15)])))
        pass

    def testSwig2Intersect2DMeshWith1DLine4(self):
        """A special case where an edge is simultaneously a cut and colinear. This tests also checks negative values in descending edges of m1."""
        i=MEDCouplingIMesh("mesh",2,[5,5],[0.,0.],[1.,1.])
        m1=i.buildUnstructured()
        part=DataArrayInt([0,1,2,3,4,7,8,11,12,13,14,15])
        m1_1=m1[part]
        m1_2=m1[part.buildComplement(m1.getNumberOfCells())]
        m1=MEDCouplingUMesh.MergeUMeshesOnSameCoords(m1_1,m1_2.buildSpreadZonesWithPoly())
        m1.zipCoords()
        m2=MEDCouplingUMesh("mesh",1) ; m2.setCoords(DataArrayDouble([(3.5,2.),(0.5,2.)])) ; m2.allocateCells()
        m2.insertNextCell(NORM_SEG2,[0,1])
        a,b,c,d=MEDCouplingUMesh.Intersect2DMeshWith1DLine(m1,m2,1e-12)
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([4,1,0,5,6,4,2,1,6,7,4,3,2,7,8,4,4,3,8,9,4,15,14,19,20,4,16,15,20,21,4,17,16,21,22,4,18,17,22,23,5,6,5,10,25,11,5,9,8,12,24,13,5,11,25,10,14,15,5,13,24,12,17,18,5,8,7,6,11,12,5,15,16,17,12,11])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25,30,35,40,46,52,58,64,70,76])))
        self.assertTrue(b.getNodalConnectivity().isEqual(DataArrayInt([1,24,12,1,12,11,1,11,25])))
        self.assertTrue(b.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6,9])))
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:24].isEqual(m1.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[24:].isEqualWithoutConsideringStr(m2.getCoords(),1e-12))
        self.assertTrue(c.isEqual(DataArrayInt([0,1,2,3,8,9,10,11,4,5,6,7,12,12])))
        self.assertTrue(d.isEqual(DataArrayInt([(9,11),(12,13),(8,10)])))
        pass

    def testSwig2Intersect2DMeshWith1DLine5(self):
        """A test focusing on a special case for cut."""
        i=MEDCouplingIMesh("mesh",2,[5,5],[0.,0.],[1.,1.])
        m1=i.buildUnstructured()
        m2=MEDCouplingUMesh("mesh",1) ; m2.setCoords(DataArrayDouble([(1.,0.),(3.,2.),(1.,4.)])) ; m2.allocateCells()
        for i in xrange(2):
            m2.insertNextCell(NORM_SEG2,[i,i+1])
            pass
        a,b,c,d=MEDCouplingUMesh.Intersect2DMeshWith1DLine(m1,m2,1e-12)
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([4,1,0,5,6,4,3,2,7,8,4,4,3,8,9,4,6,5,10,11,4,7,6,11,12,4,9,8,13,14,4,11,10,15,16,4,12,11,16,17,4,14,13,18,19,4,16,15,20,21,4,18,17,22,23,4,19,18,23,24,5,6,7,1,5,2,1,7,5,12,13,7,5,8,7,13,5,12,17,13,5,18,13,17,5,16,21,17,5,22,17,21])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,10,15,20,25,30,35,40,45,50,55,60,64,68,72,76,80,84,88,92])))
        self.assertTrue(b.getNodalConnectivity().isEqual(DataArrayInt([1,1,7,1,7,13,1,13,17,1,17,21])))
        self.assertTrue(b.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6,9,12])))
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:25].isEqual(m1.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[25:].isEqualWithoutConsideringStr(m2.getCoords(),1e-12))
        self.assertTrue(c.isEqual(DataArrayInt([0,2,3,4,5,7,8,9,11,12,14,15,1,1,6,6,10,10,13,13])))
        self.assertTrue(d.isEqual(DataArrayInt([(12,13),(14,15),(16,17),(18,19)])))
        pass

    def testIntersect2DMeshWith1DLine6(self):
        """ Basic test for Intersect2DMeshWith1DLine: a vertical line intersecting a square. """
        m1c = MEDCouplingCMesh()
        coordX = DataArrayDouble([-1., 1., 2])
        m1c.setCoordsAt(0,coordX)
        coordY = DataArrayDouble([0., 2.])
        m1c.setCoordsAt(1,coordY);
        m1 = m1c.buildUnstructured()

        # A simple line:
        m2 = MEDCouplingUMesh("bla", 1)
        coord2 = DataArrayDouble([0.,-1.0,  0.,1.,  0.,3.,  0.5,2.2], 4, 2)
        conn2 = DataArrayInt([NORM_SEG2,0,1,NORM_SEG3,1,2,3])
        connI2 = DataArrayInt([0,3,7])
        m2.setCoords(coord2)
        m2.setConnectivity(conn2, connI2)

        # End of construction of input meshes m1bis and m2 -> start of specific part of the test
        a,b,c,d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(m1, m2, 1e-10)
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([4,2,1,4,5,32,0,3,11,7,10,14,15,16,17,18,32,4,1,10,7,11,19,20,21,22,23])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,16,27])))
        self.assertTrue(b.getNodalConnectivity().isEqual(DataArrayInt([1,6,10,1,10,7,2,7,11,12,2,11,8,13])))
        self.assertTrue(b.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6,10,14])))
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:6].isEqual(m1.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[6:10].isEqual(m2.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[10:].isEqual(DataArrayDouble([(0.,0.),(0.5164175471673584,2.),(0.3796918047064557,1.43726403104512),(0.3796918047064557,2.56273596895488),(-1.,1.),(-0.24179122641632078,2.),(0.3796918047064558,1.4372640310451201),(0.,0.5),(-0.5,0.),(1.,1.),(0.5,0.),(0.,0.5),(0.3796918047064558,1.4372640310451201),(0.7582087735836792,2.)]),1e-12))
        self.assertTrue(c.isEqual(DataArrayInt([1,0,0])))
        self.assertTrue(d.isEqual(DataArrayInt([(-1,-1),(1,2),(1,2),(-1,-1)])))
        pass

    def testSwig2Intersect2DMeshWith1DLine7(self):
        """ Star pattern (a triangle intersecting another one upside down) """
        coords1 = DataArrayDouble([-2.,1.,   2.,1.,  0.,-2.], 3,2)
        coords2 = DataArrayDouble([0.,2.,   2.,-1.,  -2.,-1.,  0.,3.], 4,2)
        m1 = MEDCouplingUMesh("triangle", 2)
        m2 = MEDCouplingUMesh("tri_line", 1)
        m1.setCoords(coords1)
        m2.setCoords(coords2)
        m1.setConnectivity(DataArrayInt([NORM_TRI3, 0,1,2]), DataArrayInt([0,4]))
        m2.setConnectivity(DataArrayInt([NORM_SEG2,0,1,NORM_SEG2,1,2,NORM_SEG2,2,3]), DataArrayInt([0,3,6,9]))
    # End of construction of input meshes m1bis and m2 -> start of specific part of the test
        a,b,c,d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(m1, m2, 1e-10)
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([5,1,9,7,5,2,11,10,5,0,8,12,5,7,9,10,11,12,8])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,8,12,19])))
        self.assertTrue(b.getNodalConnectivity().isEqual(DataArrayInt([1,3,7,1,7,9,1,9,4,1,4,10,1,10,11,1,11,5,1,5,12,1,12,8,1,8,6])))
        self.assertTrue(b.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6,9,12,15,18,21,24,27])))
        self.assertTrue(a.getCoords()[:3].isEqual(m1.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[3:7].isEqual(m2.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[7:].isEqual(DataArrayDouble([(0.6666666666666666,1.),(-1.,1.),(1.3333333333333333,1.1102230246251565e-16),(0.6666666666666665,-0.9999999999999996),(-0.6666666666666667,-1.),(-1.4285714285714284,0.14285714285714302)]),1e-12))
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(c.isEqual(DataArrayInt([0,0,0,0])))
        self.assertTrue(d.isEqual(DataArrayInt([(-1,-1),(0,3),(-1,-1),(-1,-1),(1,3),(-1,-1),(-1,-1),(2,3),(-1,-1)])))
        pass

    def testSwig2Intersect2DMeshWith1DLine8(self):
        """ Line pieces ending (or fully located) in the middle of a cell """
        m1c = MEDCouplingCMesh()
        m1c.setCoordsAt(0,DataArrayDouble([-1., 1.]))
        m1c.setCoordsAt(1,DataArrayDouble([-1., 1.]));
        m1 = m1c.buildUnstructured()
        coords2 = DataArrayDouble([0.,0.,  0.,1.5, -1.5,0.,  0.5,0.0,  0.0,-0.5, 1.1,-0.6], 6,2)
        m2 = MEDCouplingUMesh("piecewise_line", 1)
        m2.setCoords(coords2)
        c = DataArrayInt([NORM_SEG2,2,1, NORM_SEG2,1,4, NORM_SEG2,4,3,  NORM_SEG2,3,5])
        cI = DataArrayInt([0,3,6,9,12])
        m2.setConnectivity(c, cI)
        a,b,c,d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(m1, m2, 1e-10)
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([5,2,11,10,5,3,13,7,8,12,5,1,0,10,11,12,8,7,13])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,4,10,19])))
        self.assertTrue(b.getNodalConnectivity().isEqual(DataArrayInt([1,6,10,1,10,11,1,11,5,1,5,12,1,12,8,1,8,7,1,7,13,1,13,9])))
        self.assertTrue(b.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6,9,12,15,18,21,24])))
        self.assertTrue(a.getCoords()[:4].isEqual(m1.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[4:10].isEqual(m2.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[10:].isEqual(DataArrayDouble([(-1.,0.5),(-0.5,1.),(0.,1.),(1.,-0.5)]),1e-12))
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(c.isEqual(DataArrayInt([0,0,0])))
        self.assertTrue(d.isEqual(DataArrayInt([(-1,-1),(0,2),(-1,-1),(-1,-1),(1,2),(1,2),(1,2),(-1,-1)])))
        pass

    def testSwig2Intersect2DMeshWith1DLine9(self):
        """ Intersection with a line whose connectivity is not consecutive """
        m1c = MEDCouplingCMesh()
        coordX = DataArrayDouble([-1., 1., 2])
        m1c.setCoordsAt(0,coordX)
        coordY = DataArrayDouble([0., 2.])
        m1c.setCoordsAt(1,coordY);
        m1 = m1c.buildUnstructured()
        # A simple line:
        m2 = MEDCouplingUMesh("bla", 1)
        coord2 = DataArrayDouble([0.,1.5,  0.5,1.,  0.0,0.5,  0.0,3.0,  0.0,-1.0], 5, 2)
        conn2 = DataArrayInt([NORM_SEG2,3,0,NORM_SEG3,0,2,1,NORM_SEG2,2,4])
        connI2 = DataArrayInt([0,3,7,10])
        m2.setCoords(coord2)
        m2.setConnectivity(conn2, connI2)
        # End of construction of input meshes m1bis and m2 -> start of specific part of the test
        a,b,c,d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(m1, m2, 1e-10)
        self.assertTrue(a.getNodalConnectivity().isEqual(DataArrayInt([4,2,1,4,5,32,4,1,11,8,6,12,14,15,16,17,18,19,32,0,3,12,6,8,11,20,21,22,23,24,25])))
        self.assertTrue(a.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,18,31])))
        self.assertTrue(b.getNodalConnectivity().isEqual(DataArrayInt([1,9,12,1,12,6,2,6,8,13,1,8,11,1,11,10])))
        self.assertTrue(b.getNodalConnectivityIndex().isEqual(DataArrayInt([0,3,6,10,13,16])))
        self.assertTrue(a.getCoords()[:6].isEqual(m1.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[6:11].isEqual(m2.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[11:].isEqual(DataArrayDouble([(0.,0.),(0.,2.),(0.5,1.),(1.,1.),(0.5,0.),(0.,0.25),(0.5,1.),(0.,1.75),(0.5,2.),(-1.,1.),(-0.5,2.),(0.,1.75),(0.5,1.),(0.,0.25),(-0.5,0.)]),1e-12))
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(c.isEqual(DataArrayInt([1,0,0])))
        self.assertTrue(d.isEqual(DataArrayInt([(-1,-1),(1,2),(1,2),(1,2),(-1,-1)])))
        pass

    def testSwig2Intersect2DMeshWith1DLine10(self):
        """ Intersection between a circle and various lines """
        eps = 1.0e-8
        m_circ = MEDCouplingDataForTest.buildCircle2(0.0, 0.0, 2.0)
        coords = [0.0,3.0,0.0,-3.0]
        connec = [0,1]
        m_line = MEDCouplingUMesh("seg", 1)
        m_line.allocateCells(1)
        meshCoords = DataArrayDouble.New(coords, len(coords)/2, 2)
        m_line.setCoords(meshCoords)
        m_line.insertNextCell(NORM_SEG2, connec)
        a, b, c, d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(m_circ, m_line, eps)
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:m_circ.getNumberOfNodes()].isEqual(m_circ.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[m_circ.getNumberOfNodes():m_circ.getNumberOfNodes()+m_line.getNumberOfNodes()].isEqual(m_line.getCoords(),1e-12))
        self.assertTrue(a.getCoords().isEqual(DataArrayDouble([(2.,0.),(1.4142135623730951,1.414213562373095),(0.,2.),(-1.414213562373095,1.4142135623730951),(-2.,0.),(-1.4142135623730954,-1.414213562373095),(0.,-2.),(1.4142135623730947,-1.4142135623730954),(0.,3.),(0.,-3.),(0.,-2.),(0.,2.),(2.,0.),(0.7653668647301797,-1.8477590650225735),(0.,0.),(0.7653668647301797,1.8477590650225735),(-2,0.),(-0.7653668647301795,1.8477590650225735),(0.,0.),(-0.7653668647301795,-1.8477590650225735)]),1e-12))
        self.assertEqual([32,1,7,10,11,12,13,14,15,32,5,3,11,10,16,17,18,19],a.getNodalConnectivity().getValues())
        self.assertEqual([0,9,18],  a.getNodalConnectivityIndex().getValues())
        self.assertEqual([1,8,11,1,11,10,1,10,9],b.getNodalConnectivity().getValues())
        self.assertEqual([0,3,6,9],b.getNodalConnectivityIndex().getValues())
        self.assertTrue(a.getCoords()[:8].isEqual(m_circ.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[8:10].isEqual(m_line.getCoords(),1e-12))
        coo_tgt = DataArrayDouble([2.,0.,1.4142135623730951,1.414213562373095,1.2246467991473532e-16,2.,-1.414213562373095,1.4142135623730951,-2.,0.,-1.4142135623730954,-1.414213562373095,-3.6739403974420594e-16,-2.,1.4142135623730947,-1.4142135623730954,0.,3.,0.,-3.,0.,-2.,0.,2.,2.,0.,0.7653668647301797,-1.8477590650225735,0.,0.,0.7653668647301797,1.8477590650225735,-2.,0.,-0.7653668647301795,1.8477590650225735,0.,0.,-0.7653668647301795,-1.8477590650225735])
        self.assertTrue(a.getCoords().isEqualWithoutConsideringStr(coo_tgt,1.0e-12))
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertEqual([0,0],c.getValues())
        self.assertEqual([-1,-1,0,1,-1,-1],d.getValues())

    def testSwig2Intersect2DMeshWith1DLine11(self):
        """ Quad line re-entering a square cell """
        eps = 1.0e-8
        m = MEDCouplingUMesh("box", 2)
        m.setCoords(DataArrayDouble([-1., -1., -1., 1., 1., 1., 1., -1.0],4,2))
        c, cI = [NORM_POLYGON, 0, 1, 2, 3], [0, 5]
        m.setConnectivity(DataArrayInt(c), DataArrayInt(cI))
        m.checkConsistencyLight()
        coords2 = [0., 1.3, -1.3, 0., -0.6, 0.6, 0., -1.3, -0.5, -0.5]
        connec2, cI2 = [NORM_SEG3, 0, 1, 2, NORM_SEG3, 1, 3, 4], [0,4,8]
        m_line = MEDCouplingUMesh("seg", 1)
        m_line.setCoords(DataArrayDouble(coords2, len(coords2)/2, 2))
        m_line.setConnectivity(DataArrayInt(connec2), DataArrayInt(cI2))
        a, b, c, d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(m, m_line, eps)
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:m.getNumberOfNodes()].isEqual(m.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[m.getNumberOfNodes():m.getNumberOfNodes()+m_line.getNumberOfNodes()].isEqual(m_line.getCoords(),1e-12))
        self.assertTrue(a.getCoords().isEqual(DataArrayDouble([(-1.,-1.),(-1.,1.),(1.,1.),(1.,-1.),(0.,1.3),(-1.3,0.),(-0.6,0.6),(0.,-1.3),(-0.5,-0.5),(-1.,0.23453685964236054),(-1.,-0.13033276368660177),(-0.2345368596423598,1.),(-0.1303327636866019,-1.),(-0.11489196370692323,1.1481421036683868),(-0.6,0.6),(-1.1481421036683859,0.11489196370692323),(-1.147455889106615,-0.0593103465193594),(-0.5,-0.5),(-0.0593103465193594,-1.147455889106615),(1.,0.),(0.4348336181566991,-1.),(-0.5651663818433009,-1.),(-1.,-0.5651663818433009),(-1.,0.05210204797787939),(-0.6,0.6),(0.3827315701788201,1.),(-0.6172684298211799,1.),(-0.6,0.6),(-1.,0.6172684298211802),(-0.6,0.6),(0.3827315701788201,1.),(1.,0.),(0.4348336181566991,-1.),(-0.5,-0.5),(-1.,0.05210204797787939),(-1.,-0.5651663818433009),(-0.5,-0.5),(-0.5651663818433009,-1.)]),1e-12))
        self.assertEqual([32,9,11,2,3,12,10,29,30,31,32,33,34,32,0,10,12,35,36,37,32,1,11,9,26,27,28],a.getNodalConnectivity().getValues())
        self.assertEqual([0,13,20,27],a.getNodalConnectivityIndex().getValues())
        self.assertEqual([2,4,11,13,2,11,9,14,2,9,5,15,2,5,10,16,2,10,12,17,2,12,7,18],b.getNodalConnectivity().getValues())
        self.assertEqual([0,4,8,12,16,20,24],b.getNodalConnectivityIndex().getValues())
        self.assertTrue(a.getCoords()[:4].isEqual(m.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[4:9].isEqual(m_line.getCoords(),1e-12))
        self.assertTrue(DataArrayInt([0,0,0]).isEqual(c))
        self.assertTrue(DataArrayInt([(-1,-1),(0,2),(-1,-1),(-1,-1),(0,1),(-1,-1)]).isEqual(d))
        pass

    def testSwig2Intersect2DMeshWith1DLine12(self):
        """ Two squares one in the other intersected by an horizontal line """
        eps = 1.0e-8
        m = MEDCouplingUMesh("boxbox", 2)
        m.setCoords(DataArrayDouble([-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,-0.5,-0.25,-0.25,-0.25,0.25,0.25,0.25,0.25,-0.25],8,2))
        c = [NORM_POLYGON, 4, 5, 6, 7, NORM_POLYGON, 0, 1, 5, 4, NORM_POLYGON, 1, 2, 3, 0, 4, 7, 6, 5]
        cI = [0, 5, 10, 19]
        m.setConnectivity(DataArrayInt(c), DataArrayInt(cI))
        m.checkConsistencyLight()
        coords2 = [-1., 0.25, 1., 0.25]
        connec2, cI2 = [NORM_SEG2, 0, 1], [0,3]
        m_line = MEDCouplingUMesh.New("seg", 1)
        m_line.setCoords(DataArrayDouble(coords2, len(coords2)/2, 2))
        m_line.setConnectivity(DataArrayInt(connec2), DataArrayInt(cI2))
        m_line2 = m_line.deepCopy()
        m2 = m.deepCopy()
        a, b, c, d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(m, m_line, eps)
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:m.getNumberOfNodes()].isEqual(m.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[m.getNumberOfNodes():m.getNumberOfNodes()+m_line.getNumberOfNodes()].isEqual(m_line.getCoords(),1e-12))
        self.assertTrue(a.getCoords().isEqual(DataArrayDouble([(-0.5,-0.5),(-0.5,0.5),(0.5,0.5),(0.5,-0.5),(-0.25,-0.25),(-0.25,0.25),(0.25,0.25),(0.25,-0.25),(-1.,0.25),(1.,0.25),(-0.5,0.25),(0.5,0.25)]),1e-12))
        self.assertEqual([5,4,5,6,7,5,1,5,10,5,4,0,10,5,5,5,1,2,11,6,5,3,0,4,7,6,11],a.getNodalConnectivity().getValues())
        self.assertEqual([0,5,9,14,20,27],a.getNodalConnectivityIndex().getValues())
        self.assertEqual([1,8,10,1,10,5,1,5,6,1,6,11,1,11,9],b.getNodalConnectivity().getValues())
        self.assertEqual([0,3,6,9,12,15],b.getNodalConnectivityIndex().getValues())
        self.assertTrue(c.isEqual(DataArrayInt([0,1,1,2,2])))
        self.assertTrue(d.isEqual(DataArrayInt([(-1,-1),(1,2),(3,0),(3,4),(-1,-1)])))
        pass

    def testSwig2Intersect2DMeshWith1DLine13(self):
        """ A square (side length) in a circle intersected by a simple horizontal line """
        import math
        eps = 1.0e-8
        m = MEDCouplingUMesh("boxcircle", 2)
        sq2 = math.sqrt(2.0)
        soth = (sq2+1.0)/2.0
        coo = [2., 0., sq2, sq2, 0., 2., -sq2, sq2, -2., 0., -sq2, -sq2, 0., -2., sq2, -sq2, -1., -1., -1., 1., 1.,
         1., 1., -1., -1., 0., 0., 1., 1., 0., 0., -1., -soth, soth, soth,soth]
        coo = DataArrayDouble(coo); coo.rearrange(2)
        m.setCoords(coo)
        c = [NORM_QPOLYG, 8, 9, 10, 11, 12, 13, 14, 15, NORM_QPOLYG, 3, 1, 10, 9, 2, 17, 13, 16, NORM_QPOLYG, 1, 7, 5, 3, 9, 8, 11, 10, 0, 6, 4, 16, 12, 15, 14, 17]
        cI = [0, 9, 18, 35]
        m.setConnectivity(DataArrayInt(c), DataArrayInt(cI))
        m.checkConsistencyLight()
        coords2 = [-2., 1., 2., 1.0]
        connec2, cI2 = [NORM_SEG2, 0, 1], [0,3]
        m_line = MEDCouplingUMesh("seg", 1)
        m_line.setCoords(DataArrayDouble(coords2, len(coords2)/2, 2))
        m_line.setConnectivity(DataArrayInt(connec2), DataArrayInt(cI2))
        a, b, c, d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(m, m_line, eps)
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:m.getNumberOfNodes()].isEqual(m.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[m.getNumberOfNodes():m.getNumberOfNodes()+m_line.getNumberOfNodes()].isEqual(m_line.getCoords(),1e-12))
        self.assertTrue(a.getCoords().isEqual(DataArrayDouble([(2.,0.),(1.4142135623730951,1.4142135623730951),(0.,2.),(-1.4142135623730951,1.4142135623730951),(-2.,0.),(-1.4142135623730951,-1.4142135623730951),(0.,-2.),(1.4142135623730951,-1.4142135623730951),(-1.,-1.),(-1.,1.),(1.,1.),(1.,-1.),(-1.,0.),(0.,1.),(1.,0.),(0.,-1.),(-1.2071067811865475,1.2071067811865475),(1.2071067811865475,1.2071067811865475),(-2.,1.),(2.,1.),(1.7320508075688772,1.),(-1.7320508075688772,1.),(-1.2071067811865475,1.2071067811865475),(-1.3660254037844386,1.),(-1.58670668058247,1.2175228580174415),(0.,-1.),(1.,0.),(1.2071067811865475,1.2071067811865475),(1.5867066805824703,1.2175228580174413),(1.9828897227476205,-0.26105238444010315),(0.,-2.),(-1.9828897227476205,-0.2610523844401032),(-1.3660254037844386,1.),(-1.,0.),(1.5867066805824703,1.2175228580174413),(1.3660254037844386,1.),(1.2071067811865475,1.2071067811865475),(0.,-2.),(-1.9828897227476205,-0.2610523844401032),(-1.3660254037844386,1.),(-1.,0.),(0.,-1.),(1.,0.),(1.3660254037844386,1.),(1.9828897227476205,-0.26105238444010315)]),1e-12))
        self.assertEqual([32,8,9,10,11,12,13,14,15,32,3,1,10,9,2,17,13,16,32,3,9,21,22,23,24,32,1,20,10,34,35,36,32,7,5,21,9,8,11,10,20,37,38,39,40,41,42,43,44],a.getNodalConnectivity().getValues())
        self.assertEqual([0,9,18,25,32,49],a.getNodalConnectivityIndex().getValues())
        self.assertEqual([1,18,21,1,21,9,1,9,10,1,10,20,1,20,19],b.getNodalConnectivity().getValues())
        self.assertEqual([0,3,6,9,12,15],b.getNodalConnectivityIndex().getValues())
        self.assertTrue(c.isEqual(DataArrayInt([0,1,2,2,2])))
        self.assertTrue(d.isEqual(DataArrayInt([(-1,-1),(2,4),(1,0),(3,4),(-1,-1)])))
        pass

    def testSwig2Intersect2DMeshWith1DLine14(self):
        """ A circle in a circle intersected by a simple horizontal line, not tangent to the circles """
        eps = 1.0e-8
        m = MEDCouplingUMesh("boxcircle", 2)
        coo = [2.,0.,1.4142135623730951,1.414213562373095,0.,2.,-1.414213562373095,1.4142135623730951,-2.,0.,-1.4142135623730954,-1.414213562373095,0.,-2.,
               1.4142135623730947,-1.4142135623730954,1.,0.,0.7071067811865476,0.7071067811865475,0.,1.,-0.7071067811865475,0.7071067811865476,-1.,0.,-0.7071067811865477,-0.7071067811865475,
               0.,-1.,0.7071067811865474,-0.7071067811865477,1.060660171779821,-1.0606601717798214,-1.0606601717798214,-1.0606601717798212]
        coo = DataArrayDouble(coo); coo.rearrange(2)
        m.setCoords(coo)
        c = [NORM_QPOLYG, 15, 13, 11, 9, 14, 12, 10, 8, NORM_QPOLYG, 7, 5, 13, 15, 6, 17, 14, 16, NORM_QPOLYG, 5, 3, 1, 7, 15, 9, 11, 13, 4, 2, 0, 16, 8, 10, 12, 17]
        cI = [0, 9, 18, 35]
        m.setConnectivity(DataArrayInt(c), DataArrayInt(cI))
        m.checkConsistencyLight()
        coords2 = [-2., 0., 2., 0.]
        connec2, cI2 = [NORM_SEG2, 0, 1], [0,3]
        m_line = MEDCouplingUMesh.New("seg", 1)
        m_line.setCoords(DataArrayDouble(coords2, len(coords2)/2, 2))
        m_line.setConnectivity(DataArrayInt(connec2), DataArrayInt(cI2))
        a, b, c, d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(m, m_line, eps)
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:m.getNumberOfNodes()].isEqual(m.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[m.getNumberOfNodes():m.getNumberOfNodes()+m_line.getNumberOfNodes()].isEqual(m_line.getCoords(),1e-12))
        self.assertTrue(a.getCoords().isEqual(DataArrayDouble([(2.,0.),(1.4142135623730951,1.414213562373095),(0.,2.),(-1.414213562373095,1.4142135623730951),(-2.,0.),(-1.4142135623730954,-1.414213562373095),(0.,-2.),(1.4142135623730947,-1.4142135623730954),(1.,0.),(0.7071067811865476,0.7071067811865475),(0.,1.),(-0.7071067811865475,0.7071067811865476),(-1.,0.),(-0.7071067811865477,-0.7071067811865475),(0.,-1.),(0.7071067811865474,-0.7071067811865477),(1.060660171779821,-1.0606601717798214),(-1.0606601717798214,-1.0606601717798212),(-2.,0.),(2.,0.),(-1.,0.),(1.,0.),(0.,2.),(1.8477590650225735,0.7653668647301795),(1.8477590650225735,-0.7653668647301797),(1.060660171779821,-1.0606601717798214),(0.9238795325112867,-0.38268343236508984),(0.9238795325112867,0.3826834323650897),(0.,1.),(-0.9238795325112867,0.3826834323650896),(-1.5,0.),(-1.8477590650225735,0.7653668647301792),(-1.0606601717798214,-1.0606601717798212),(-1.8477590650225733,-0.7653668647301799),(-1.5,0.),(-0.9238795325112866,-0.38268343236508995),(0.,1.),(-0.9238795325112867,0.3826834323650896),(-1.5,0.),(-1.8477590650225735,0.7653668647301792),(0.,2.),(1.8477590650225735,0.7653668647301795),(1.5,0.),(0.9238795325112867,0.3826834323650897),(1.060660171779821,-1.0606601717798214),(0.9238795325112867,-0.38268343236508984),(1.5,0.),(1.8477590650225735,-0.7653668647301797),(0.,1.),(0.9238795325112867,0.3826834323650897),(0.,0.),(-0.9238795325112867,0.3826834323650896),(0.,-1.),(-0.9238795325112866,-0.38268343236508995),(0.,0.),(0.9238795325112867,-0.38268343236508984)]),1e-12))
        self.assertEqual([32,7,5,13,15,6,17,14,16,32,9,11,20,18,3,1,19,21,36,37,38,39,40,41,42,43,32,7,15,21,19,44,45,46,47,32,13,5,18,20,32,33,34,35,32,11,9,21,20,48,49,50,51,32,15,13,20,21,52,53,54,55],a.getNodalConnectivity().getValues())
        self.assertEqual([0,9,26,35,44,53,62],a.getNodalConnectivityIndex().getValues())
        self.assertEqual([1,18,20,1,20,21,1,21,19],b.getNodalConnectivity().getValues())
        self.assertEqual([0,3,6,9],b.getNodalConnectivityIndex().getValues())
        self.assertTrue(c.isEqual(DataArrayInt([1,2,2,2,0,0])))
        self.assertTrue(d.isEqual(DataArrayInt([(1,3),(4,5),(1,2)])))
        pass

    def testSwig2Intersect2DMeshWith1DLine15(self):
        """ Same as testSwig2Intersect2DMeshWith1DLine13 except that the line is colinear AND splits on of the common edge of 2D mesh."""
        import math
        eps = 1.0e-8
        m = MEDCouplingUMesh("boxcircle", 2)
        sq2 = math.sqrt(2.0)
        soth = (sq2+1.0)/2.0
        coo = [2., 0., sq2, sq2, 0., 2., -sq2, sq2, -2., 0., -sq2, -sq2, 0., -2., sq2, -sq2, -1., -1., -1., 1., 1.,
         1., 1., -1., -1., 0., 0., 1., 1., 0., 0., -1., -soth, soth, soth,soth]
        coo = DataArrayDouble(coo); coo.rearrange(2)
        m.setCoords(coo)
        c = [NORM_QPOLYG, 8, 9, 10, 11, 12, 13, 14, 15, NORM_QPOLYG, 3, 1, 10, 9, 2, 17, 13, 16, NORM_QPOLYG, 1, 7, 5, 3, 9, 8, 11, 10, 0, 6, 4, 16, 12, 15, 14, 17]
        cI = [0, 9, 18, 35]
        m.setConnectivity(DataArrayInt(c), DataArrayInt(cI))
        m.checkConsistencyLight()
        coords2 = [(-2., 1.),(2.,1.),(0.,1)]
        connec2, cI2 = [NORM_SEG2, 0, 2, NORM_SEG2, 2, 1], [0,3,6]
        m_line = MEDCouplingUMesh("seg", 1)
        m_line.setCoords(DataArrayDouble(coords2))
        m_line.setConnectivity(DataArrayInt(connec2), DataArrayInt(cI2))
        a, b, c, d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(m, m_line, eps)
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:m.getNumberOfNodes()].isEqual(m.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[m.getNumberOfNodes():m.getNumberOfNodes()+m_line.getNumberOfNodes()].isEqual(m_line.getCoords(),1e-12))
        self.assertTrue(a.getCoords().isEqual(DataArrayDouble([(2.,0.),(1.4142135623730951,1.4142135623730951),(0.,2.),(-1.4142135623730951,1.4142135623730951),(-2.,0.),(-1.4142135623730951,-1.4142135623730951),(0.,-2.),(1.4142135623730951,-1.4142135623730951),(-1.,-1.),(-1.,1.),(1.,1.),(1.,-1.),(-1.,0.),(0.,1.),(1.,0.),(0.,-1.),(-1.2071067811865475,1.2071067811865475),(1.2071067811865475,1.2071067811865475),(-2.,1.),(2.,1.),(0.,1.),(1.7320508075688776,1.),(-1.7320508075688776,1.),(-0.5,1.),(0.5,1.),(0.5,1.),(-0.5,1.),(-1.2071067811865475,1.2071067811865475),(-1.3660254037844388,1.),(-1.58670668058247,1.2175228580174415),(0.,-1.),(1.,0.),(1.2071067811865475,1.2071067811865475),(1.5867066805824703,1.2175228580174413),(1.9828897227476205,-0.26105238444010315),(0.,-2.),(-1.9828897227476205,-0.2610523844401032),(-1.3660254037844388,1.),(-1.,0.),(1.5867066805824703,1.2175228580174413),(1.3660254037844388,1.),(1.2071067811865475,1.2071067811865475),(0.,-2.),(-1.9828897227476205,-0.2610523844401032),(-1.3660254037844388,1.),(-1.,0.),(0.,-1.),(1.,0.),(1.3660254037844388,1.),(1.9828897227476205,-0.26105238444010315)]),1e-12))
        self.assertEqual([32,8,9,20,10,11,12,23,24,14,15,32,3,1,10,20,9,2,17,25,26,16,32,3,9,22,27,28,29,32,1,21,10,39,40,41,32,7,5,22,9,8,11,10,21,42,43,44,45,46,47,48,49],a.getNodalConnectivity().getValues())
        self.assertEqual([0,11,22,29,36,53],a.getNodalConnectivityIndex().getValues())
        self.assertEqual([1,18,22,1,22,9,1,9,20,1,20,10,1,10,21,1,21,19],b.getNodalConnectivity().getValues())
        self.assertEqual([0,3,6,9,12,15,18],b.getNodalConnectivityIndex().getValues())
        self.assertTrue(c.isEqual(DataArrayInt([0,1,2,2,2])))
        self.assertTrue(d.isEqual(DataArrayInt([(-1,-1),(2,4),(1,0),(1,0),(3,4),(-1,-1)])))
        pass

    def testSwig2Intersect2DMeshWith1DLine16(self):
        """ Same than testSwig2Intersect2DMeshWith1DLine13 except it is a vertical line. Non regression test."""
        import math
        eps = 1.0e-8
        m = MEDCouplingUMesh("boxcircle", 2)
        sq2 = math.sqrt(2.0)
        soth = (sq2+1.0)/2.0
        coo = [2., 0., sq2, sq2, 0., 2., -sq2, sq2, -2., 0., -sq2, -sq2, 0., -2., sq2, -sq2, -1., -1., -1., 1., 1.,
         1., 1., -1., -1., 0., 0., 1., 1., 0., 0., -1., -soth, soth, soth,soth]
        coo = DataArrayDouble(coo); coo.rearrange(2)
        m.setCoords(coo)
        c = [NORM_QPOLYG, 8, 9, 10, 11, 12, 13, 14, 15, NORM_QPOLYG, 3, 1, 10, 9, 2, 17, 13, 16, NORM_QPOLYG, 1, 7, 5, 3, 9, 8, 11, 10, 0, 6, 4, 16, 12, 15, 14, 17]
        cI = [0, 9, 18, 35]
        m.setConnectivity(DataArrayInt(c), DataArrayInt(cI))
        m.checkConsistencyLight()
        coords2 = [1., 2., 1., -2.]
        connec2, cI2 = [NORM_SEG2, 0, 1], [0,3]
        m_line = MEDCouplingUMesh("seg", 1)
        m_line.setCoords(DataArrayDouble(coords2, len(coords2)/2, 2))
        m_line.setConnectivity(DataArrayInt(connec2), DataArrayInt(cI2))
        a, b, c, d = MEDCouplingUMesh.Intersect2DMeshWith1DLine(m, m_line, eps)
        self.assertTrue(a.getCoords().getHiddenCppPointer()==b.getCoords().getHiddenCppPointer())
        self.assertTrue(a.getCoords()[:m.getNumberOfNodes()].isEqual(m.getCoords(),1e-12))
        self.assertTrue(a.getCoords()[m.getNumberOfNodes():m.getNumberOfNodes()+m_line.getNumberOfNodes()].isEqual(m_line.getCoords(),1e-12))
        self.assertTrue(a.getCoords().isEqual(DataArrayDouble([(2., 0.),(1.4142135623730951,1.4142135623730951),(0.,2.),(-1.4142135623730951,1.4142135623730951),(-2.,0.),(-1.4142135623730951,-1.4142135623730951),(0.,-2.),(1.4142135623730951,-1.4142135623730951),(-1.,-1.),(-1.,1.),(1.,1.),(1.,-1.),(-1.,0.),(0.,1.),(1.,0.),(0.,-1.),(-1.2071067811865475,1.2071067811865475),(1.2071067811865475,1.2071067811865475),(1.,2.),(1.,-2.),(1.,1.7320508075688772),(1.,-1.7320508075688772),(1.2071067811865475,1.2071067811865475),(1.,1.3660254037844386),(1.217522858017441,1.5867066805824703),(-1.2071067811865475,1.2071067811865475),(-0.2610523844401028,1.9828897227476208),(1.,1.3660254037844386),(0.,1.),(1.2071067811865475,1.2071067811865475),(2.,0.),(1.217522858017441,-1.5867066805824703),(1.,-1.3660254037844386),(1.,0.),(-2.,0.),(-1.2071067811865475,1.2071067811865475),(-1.,0.),(0.,-1.),(1.,-1.3660254037844386),(-0.2610523844401028,-1.9828897227476208)]),1e-12))
        self.assertEqual([32,8,9,10,11,12,13,14,15,32,1,10,20,22,23,24,32,9,3,20,10,25,26,27,28,32,10,1,7,21,11,29,30,31,32,33,32,5,3,9,8,11,21,34,35,36,37,38,39],a.getNodalConnectivity().getValues())
        self.assertEqual([0,9,16,25,36,49],a.getNodalConnectivityIndex().getValues())
        self.assertEqual([1,18,20,1,20,10,1,10,11,1,11,21,1,21,19],b.getNodalConnectivity().getValues())
        self.assertEqual([0,3,6,9,12,15],b.getNodalConnectivityIndex().getValues())
        self.assertTrue(c.isEqual(DataArrayInt([0,1,1,2,2])))
        self.assertTrue(d.isEqual(DataArrayInt([(-1,-1),(1,2),(3,0),(3,4),(-1,-1)])))
        pass

    def testSwig2Conformize2D1(self):
        eps = 1.0e-8
        coo = [0.,-0.5,0.,0.,0.5,0.,0.5,-0.5,0.25,
               -0.1,0.25,0.,0.5,-0.1,0.,0.5,0.5,0.5,0.25,0.4,0.25,0.5,0.5,0.4]
        conn = [5,5,2,6,4,5,6,3,0,1,5,4,5,10,8,11,9,5,11,2,1,7,10,9]
        connI = [0,5,12,17,24]
        m = MEDCouplingUMesh("box",2)
        cooArr = DataArrayDouble(coo,len(coo)/2,2)
        m.setCoords(cooArr)
        m.setConnectivity(DataArrayInt(conn),DataArrayInt(connI))
        m.mergeNodes(eps)
        m.checkConsistencyLight()
        self.assertTrue(m.conformize2D(eps).isEqual(DataArrayInt([3])))
        self.assertEqual(m.getCoords().getHiddenCppPointer(),cooArr.getHiddenCppPointer()) # check that coordinates remain the same here
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,5,2,6,4,5,6,3,0,1,5,4,5,10,8,11,9,5,11,2,5,1,7,10,9])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,5,12,17,25])))
        pass

    def testSwig2Conformize2D2(self):
        eps = 1.0e-8
        coo=DataArrayDouble([-10,-6,0,-6,0,0,7,0,-10,2,0,2,0,6,7,6,0,8,7,8,-10,12,-4,12,0,12,0,11,7,11,-4,16,0,16,7,16],18,2)
        conn=DataArrayInt([2,3,7,6, 13,16,17,14, 4,10,12,5, 9,14,13,8, 8,9,7,6, 5,4,0,1, 16,12,11,15])
        m=MEDCoupling1SGTUMesh("mesh",NORM_QUAD4)
        m.setCoords(coo)
        m.setNodalConnectivity(conn)
        m=m.buildUnstructured()
        self.assertTrue(m.conformize2D(eps).isEqual(DataArrayInt([0,1,2,5])))
        self.assertEqual(m.getCoords().getHiddenCppPointer(),coo.getHiddenCppPointer()) # check that coordinates remain the same here
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([5,2,3,7,6,5, 5,13,12,16,17,14, 5,4,10,11,12,13,8,6,5, 4,9,14,13,8, 4,8,9,7,6, 5,5,4,0,1,2, 4,16,12,11,15])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,6,12,21,26,31,37,42])))
        pass

    def testSwigSplit2DCells1(self):
        coo=DataArrayDouble([[0,0],[1,0],[1,1],[0,1],[0.5,0],[1,0.5],[0.5,1],[0.,0.5]])
        m=MEDCouplingUMesh("mesh",2)
        m.setCoords(coo)
        m.allocateCells()
        m.insertNextCell(NORM_QUAD8,[0,1,2,3,4,5,6,7])
        _,d,di,_,_=m.buildDescendingConnectivity()
        subb=DataArrayInt([5])
        subbi=DataArrayInt([0,0,1,1,1])
        mid=DataArrayInt([-1,-1])
        midi=DataArrayInt([0,0,2,2,2])
        self.assertEqual(2,m.split2DCells(d,di,subb,subbi,mid,midi))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([32,0,1,5,2,3,4,8,9,6,7])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,11])))
        self.assertTrue(m.getCoords().isEqual(DataArrayDouble([[0,0],[1,0],[1,1],[0,1],[0.5,0],[1,0.5],[0.5,1],[0.,0.5],[1.,0.25],[1.,0.75]]),1e-12))
        pass

    def testSwig2Conformize2D3(self):
        eps = 1.0e-8
        coo=DataArrayDouble([-10,-6,0,-6,0,0,7,0,-10,2,0,2,0,6.5,7,6.5,0,8,7,8,-10,12,-4,12,0,12,0,11,7,11,-4,16,0,16,7,16],18,2)
        conn=DataArrayInt([2,3,7,6, 13,16,17,14, 4,10,12,5, 9,14,13,8, 8,9,7,6, 5,4,0,1, 16,12,11,15])
        m=MEDCoupling1SGTUMesh("mesh",NORM_QUAD4)
        m.setCoords(coo)
        m.setNodalConnectivity(conn)
        m=m.buildUnstructured()
        m.convertLinearCellsToQuadratic()
        self.assertTrue(m.conformize2D(eps).isEqual(DataArrayInt([0,1,2,5])))
        self.assertTrue(m.getCoords().getHiddenCppPointer()!=coo.getHiddenCppPointer()) # coordinates are not the same here contrary to testSwig2Conformize2D2 ...
        self.assertTrue(m.getCoords()[:18].isEqual(coo,1e-12)) # but the 18 first nodes are the same
        pass

    def testSwig2Conformize2D4(self):
        eps = 1.0e-8
        coo=DataArrayDouble([-10,-6,0,-6,0,0,7,0,-10,2,0,2,0,6.5,7,6.5,0,8,7,8,-10,12,-4,12,0,12,0,11,7,11,-4,16,0,16,7,16],18,2)
        conn=DataArrayInt([2,3,7,6, 13,16,17,14, 4,10,12,5, 9,14,13,8, 8,9,7,6, 5,4,0,1, 16,12,11,15])
        m=MEDCoupling1SGTUMesh("mesh",NORM_QUAD4)
        m.setCoords(coo)
        m.setNodalConnectivity(conn)
        m=m.buildUnstructured()
        m.convertLinearCellsToQuadratic()
        self.assertEqual(42,m.getNumberOfNodes())
        oldCoo=m.getCoords().deepCopy()
        m.conformize2D(eps)
        self.assertTrue(m.getCoords()[:42].isEqual(oldCoo,1e-12))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([32,2,3,7,6,5,18,19,20,42,43,32,13,12,16,17,14,44,38,23,24,25,32,4,10,11,12,13,8,6,5,26,45,39,44,31,34,42,29,8,9,14,13,8,30,25,31,32,8,8,9,7,6,32,33,20,34,32,5,4,0,1,2,29,35,36,46,43,8,16,12,11,15,38,39,40,41])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,11,22,39,48,57,68,77])))
        self.assertTrue(m.getCoords().isEqual(DataArrayDouble([[-10.,-6.0],[0.,-6.0],[0.,0.0],[7.,0.0],[-10.,2.0],[0.,2.0],[0.,6.5],[7.,6.5],[0.,8.0],[7.,8.0],[-10.,12.0],[-4.,12.0],[0.,12.0],[0.,11.0],[7.,11.0],[-4.,16.0],[0.,16.0],[7.,16.0],[3.5, 0.0],[7.,3.25],[3.5, 6.5],[0.,3.25],[0.,13.5],[3.5, 16.0],[7.,13.5],[3.5, 11.0],[-10.,7.0],[-5.,12.0],[0.,7.0],[-5.,2.0],[7.,9.5],[0.,9.5],[3.5, 8.0],[7.,7.25],[0.,7.25],[-10.,-2.0],[-5.,-6.0],[0.,-2.0],[0.,14.0],[-2.,12.0],[-4.,14.0],[-2.,16.0],[0.,4.25],[0.,1.0],[0.,11.5],[-7.,12.0],[0.,-3.]]),1e-12))
        pass

    def testSwig2Conformize2D5(self):
        eps=1e-8
        coo=DataArrayDouble([[2,2],[2,-6],[10,-2],[-2,-2],[6,0],[6,-4],[2,7],[2,4.5],[-1.4641016151377544,0],[-1.950753362380551,-1.3742621398390762],[-7,-3],[-0.8284271247461898,-4.82842712474619],[0.26794919243112281,3.5],[0,1.4641016151377548],[-4.4753766811902755,-2.1871310699195381],[-3.9142135623730949,-3.9142135623730949],[-1.8042260651806146,-3.23606797749979]])
        m=MEDCouplingUMesh("mesh",2)
        m.allocateCells()
        m.setCoords(coo)
        m.insertNextCell(NORM_TRI6,[1,2,0,5,4,3])
        m.insertNextCell(NORM_TRI6,[8,6,0,12,7,13])
        m.insertNextCell(NORM_TRI6,[11,9,10,16,14,15])
        self.assertTrue(m.conformize2D(eps).isEqual(DataArrayInt([0])))
        self.assertTrue(m.getCoords().isEqual(DataArrayDouble([2.,2.,2.,-6.,10.,-2.,-2.,-2.,6.,0.,6.,-4.,2.,7.,2.,4.5,-1.4641016151377544,0.,-1.950753362380551,-1.3742621398390762,-7.,-3.,-0.8284271247461898,-4.82842712474619,0.2679491924311228,3.5,8.881784197001252e-16,1.4641016151377548,-4.4753766811902755,-2.187131069919538,-3.914213562373095,-3.914213562373095,-1.8042260651806146,-3.236067977499789,-1.7705659643687133,-0.6647725630649153,0.46926627053963865,-5.695518130045146],19,2),1e-12))
        self.assertTrue(m.getNodalConnectivity().isEqual(DataArrayInt([32,1,2,0,8,9,11,5,4,13,17,16,18,6,8,6,0,12,7,13,6,11,9,10,16,14,15])))
        self.assertTrue(m.getNodalConnectivityIndex().isEqual(DataArrayInt([0,13,20,27])))
        pass

    def testSwig2Conformize3D1(self):
        """ Simple test where no edge merge is required, only face merging (first part of the algo) """
        mesh = MEDCouplingUMesh('merge', 3)
        coo = DataArrayDouble([(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1),(1,0.5,0.5),
                               (2,0,0),(2,1,0),(2,0,1),(2,1,1),(1,0,2),(2,0,2),(1,1,2),(2,1,2)])
        mesh.setCoords(coo)
        c = DataArrayInt([31, 1, 0, 2, 3, -1, 5, 7, 6, 4, -1, 1, 5, 4, 0, -1, 0, 4, 6, 2, -1, 2, 6, 7, 3, -1,
                           3, 7, 8, -1, 7, 5, 8, -1, 5, 1, 8, -1, 1, 3, 8, 31, 9, 1, 3, 10, -1, 11, 12, 7, 5, -1,
                           9, 11, 5, 1, -1, 1, 5, 7, 3, -1, 3, 7, 12, 10, -1, 10, 12, 11, 9, 31, 11, 5, 7, 12, -1,
                           14, 16, 15, 13, -1, 11, 14, 13, 5, -1, 5, 13, 15, 7, -1, 7, 15, 16, 12, -1, 12, 16, 14, 11])
        cI = DataArrayInt([0, 41, 71, 101])
        mesh.setConnectivity(c, cI)
        ret = mesh.conformize3D(1.0e-8)

        mretDesc, _, _, _, _ = mesh.buildDescendingConnectivity()
        c, cI = mretDesc.getNodalConnectivity().getValues(), mretDesc.getNodalConnectivityIndex().getValues()
        cRef = [5, 1, 0, 2, 3, 5, 5, 7, 6, 4, 5, 1, 5, 4, 0, 5, 0, 4, 6, 2, 5, 2, 6, 7, 3, 5, 3, 7, 8, 5,
                7, 5, 8, 5, 5, 1, 8, 5, 1, 3, 8, 5, 9, 1, 3, 10, 5, 11, 12, 7, 5, 5, 9, 11, 5, 1, 5, 3,
                7, 12, 10, 5, 10, 12, 11, 9, 5, 14, 16, 15, 13, 5, 11, 14, 13, 5, 5, 5, 13, 15, 7, 5, 7,
                15, 16, 12, 5, 12, 16, 14, 11]
        cIRef = [0, 5, 10, 15, 20, 25, 29, 33, 37, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86, 91]
        self.assertEqual(19, mretDesc.getNumberOfCells())
        self.assertEqual(cRef, c)
        self.assertEqual(cIRef, cI)
        self.assertEqual([1], ret.getValues())
        pass

    def testSwig2Conformize3D2(self):
        """ More advanced test where edge merge is required. """
        mesh = MEDCouplingUMesh('merge', 3)
        coo = DataArrayDouble([(0,0,0),(1,0,0),(0,1,0),(1,1,0),(0,0,1),(1,0,1),(0,1,1),(1,1,1),(0,0,2),(1,0,2),(0,1,2),
                               (1,1,2),(1,0.6,0),(1,0.6,1),(1,0.3,1),(1,0.3,2),(2,0,0),(2,1,0),(2,0,1),(2,1,1),(2,0,2),(2,1,2)])
        mesh.setCoords(coo)
        c = DataArrayInt([31, 1, 0, 2, 3, 12, -1, 5, 13, 7, 6, 4, -1, 1, 5, 4, 0, -1, 0, 4, 6, 2, -1, 2, 6, 7, 3, -1, 3, 7, 13, 12, -1,
                          5, 1, 12, 13, 31, 5, 4, 6, 7, 14, -1, 9, 15, 11, 10, 8, -1, 5, 9, 8, 4, -1, 4, 8, 10, 6, -1, 6, 10, 11, 7, -1,
                          7, 11, 15, 14, -1, 9, 5, 14, 15, 31, 16, 1, 3, 17, -1, 18, 19, 7, 5, -1, 16, 18, 5, 1, -1, 1, 5, 7, 3, -1,
                          3, 7, 19, 17, -1, 17, 19, 18, 16, 31, 18, 5, 7, 19, -1, 20, 21, 11, 9, -1, 18, 20, 9, 5, -1, 5, 9, 11, 7, -1,
                          7, 11, 21, 19, -1, 19, 21, 20, 18])
        cI = DataArrayInt([0, 37, 74, 104, 134])
        mesh.setConnectivity(c, cI)

        ret = mesh.conformize3D(1.0e-8)

        mretDesc, _, _, _, _ = mesh.buildDescendingConnectivity()
        mretDesc2, _, _, _, _ = mretDesc.buildDescendingConnectivity()
        c0, cI0 = mesh.getNodalConnectivity().getValues(), mesh.getNodalConnectivityIndex().getValues()
        c, cI = mretDesc.getNodalConnectivity().getValues(), mretDesc.getNodalConnectivityIndex().getValues()
        c2, cI2 = mretDesc2.getNodalConnectivity().getValues(), mretDesc2.getNodalConnectivityIndex().getValues()
        cRef0 = [31, 1, 0, 2, 3, 12, -1, 5, 14, 13, 7, 6, 4, -1, 1, 5, 4, 0, -1, 0, 4, 6, 2, -1, 2, 6, 7, 3, -1, 3, 7, 13, 12, -1, 5, 1, 12, 13, 14, 31, 9, 15, 11, 10, 8, -1, 5, 9, 8, 4, -1, 4, 8, 10, 6, -1, 6, 10, 11, 7, -1, 7, 11, 15, 14, 13, -1, 9, 5, 14, 15, -1, 5, 14, 13, 7, 6, 4, 31, 16, 1, 12, 3, 17, -1, 18, 19, 7, 13, 14, 5, -1, 16, 18, 5, 1, -1, 3, 7, 19, 17, -1, 17, 19, 18, 16, -1, 5, 1, 12, 13, 14, -1, 3, 7, 13, 12, 31, 5, 14, 13, 7, 19, 18, -1, 20, 21, 11, 15, 9, -1, 18, 20, 9, 5, -1, 7, 11, 21, 19, -1, 19, 21, 20, 18, -1, 9, 5, 14, 15, -1, 7, 11, 15, 14, 13]
        cIRef0 = [0, 39, 78, 117, 156]
        cRef = [5, 1, 0, 2, 3, 12, 5, 5, 14, 13, 7, 6, 4, 5, 1, 5, 4, 0, 5, 0, 4, 6, 2, 5, 2, 6, 7, 3, 5, 3, 7, 13, 12, 5, 5, 1, 12, 13, 14, 5, 9, 15, 11, 10, 8, 5, 5, 9, 8, 4, 5, 4, 8, 10, 6, 5, 6, 10, 11, 7, 5, 7, 11, 15, 14, 13, 5, 9, 5, 14, 15, 5, 16, 1, 12, 3, 17, 5, 18, 19, 7, 13, 14, 5, 5, 16, 18, 5, 1, 5, 3, 7, 19, 17, 5, 17, 19, 18, 16, 5, 20, 21, 11, 15, 9, 5, 18, 20, 9, 5, 5, 7, 11, 21, 19, 5, 19, 21, 20, 18]
        cIRef = [0, 6, 13, 18, 23, 28, 33, 39, 45, 50, 55, 60, 66, 71, 77, 84, 89, 94, 99, 105, 110, 115, 120]
        cRef2 = [1, 1, 0, 1, 0, 2, 1, 2, 3, 1, 3, 12, 1, 12, 1, 1, 5, 14, 1, 14, 13, 1, 13, 7, 1, 7, 6, 1, 6, 4, 1, 4, 5, 1, 1, 5, 1, 4, 0, 1, 6, 2, 1, 7, 3, 1, 13, 12, 1, 9, 15, 1, 15, 11, 1, 11, 10, 1, 10, 8, 1, 8, 9, 1, 5, 9, 1, 8, 4, 1, 10, 6, 1, 11, 7, 1, 15, 14, 1, 16, 1, 1, 3, 17, 1, 17, 16, 1, 18, 19, 1, 19, 7, 1, 5, 18, 1, 16, 18, 1, 19, 17, 1, 20, 21, 1, 21, 11, 1, 9, 20, 1, 18, 20, 1, 21, 19]
        cIRef2 = [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 117]
        self.assertEqual(22, mretDesc.getNumberOfCells())
        self.assertEqual(39, mretDesc2.getNumberOfCells())
        self.assertEqual(cRef0, c0)
        self.assertEqual(cIRef0, cI0)
        self.assertEqual(cRef, c)
        self.assertEqual(cIRef, cI)
        self.assertEqual(cRef2, c2)
        self.assertEqual(cIRef2, cI2)
        self.assertEqual(set([0,1,2,3]), set(ret.getValues()))
        pass

    def testSwig2Conformize3D3(self):
        """ LMEC's case (hexagonal prism) """
        eps = 1.0e-7   # 1.0e-8 is too fine
        mesh = MEDCouplingUMesh('merge', 3)
        coo = DataArrayDouble([(0,0,0),(-0.0054,-0.00935307,0),(0.0054,-0.00935307,0),(0.0108,0,0),(0.0054,0.00935307,0),(-0.0054,0.00935307,0),(-0.0108,6.93889e-18,0),(-0.0054,-0.0153031,0),(0.0054,-0.0153031,0),(0.0105529,-0.0123281,0),(0.0159529,-0.002975,0),(0.0159529,0.002975,0),(0.0105529,0.0123281,0),(0.0054,0.0153031,0),(-0.0054,0.0153031,0),(-0.0105529,0.0123281,0),(-0.0159529,0.002975,0),(-0.0159529,-0.002975,0),(-0.0105529,-0.0123281,0),(-0.00883523,-0.0153031,0),(0.00883523,-0.0153031,0),(0.0176705,0,0),(0.00883523,0.0153031,0),(-0.00883523,0.0153031,0),(-0.0176705,2.16401e-18,0),(0,0,0.05),(-0.0054,-0.00935307,0.05),(0.0054,-0.00935307,0.05),(0.0108,0,0.05),(0.0054,0.00935307,0.05),(-0.0054,0.00935307,0.05),(-0.0108,6.93889e-18,0.05),(-0.0054,-0.0153031,0.05),(0.0054,-0.0153031,0.05),(0.0105529,-0.0123281,0.05),(0.0159529,-0.002975,0.05),(0.0159529,0.002975,0.05),(0.0105529,0.0123281,0.05),(0.0054,0.0153031,0.05),(-0.0054,0.0153031,0.05),(-0.0105529,0.0123281,0.05),(-0.0159529,0.002975,0.05),(-0.0159529,-0.002975,0.05),(-0.0105529,-0.0123281,0.05),(-0.00883523,-0.0153031,0.05),(0.00883523,-0.0153031,0.05),(0.0176705,0,0.05),(0.00883523,0.0153031,0.05),(-0.00883523,0.0153031,0.05),(-0.0176705,2.16401e-18,0.05),(0.0176705,0,0.05),(0.00883523,0.0153031,0.05),(-0.00883523,0.0153031,0.05),(-0.0176705,2.16401e-18,0.05),(-0.00883523,-0.0153031,0.05),(0.00883523,-0.0153031,0.05),(0.0176705,0,0.15),(0.00883523,0.0153031,0.15),(-0.00883523,0.0153031,0.15),(-0.0176705,2.16401e-18,0.15),(-0.00883523,-0.0153031,0.15),(0.00883523,-0.0153031,0.15)])
        mesh.setCoords(coo)
        c = DataArrayInt([31, 1, 0, 2, -1, 25, 26, 27, -1, 2, 27, 26, 1, -1, 0, 25, 27, 2, -1, 1, 26, 25, 0, 31, 2, 0, 3, -1, 25, 27, 28, -1, 3, 28, 27, 2, -1, 0, 25, 28, 3, -1, 2, 27, 25, 0, 31, 3, 0, 4, -1, 25, 28, 29, -1, 4, 29, 28, 3, -1, 0, 25, 29, 4, -1, 3, 28, 25, 0, 31, 4, 0, 5, -1, 25, 29, 30, -1, 5, 30, 29, 4, -1, 0, 25, 30, 5, -1, 4, 29, 25, 0, 31, 5, 0, 6, -1, 25, 30, 31, -1, 6, 31, 30, 5, -1, 0, 25, 31, 6, -1, 5, 30, 25, 0, 31, 6, 0, 1, -1, 25, 31, 26, -1, 1, 26, 31, 6, -1, 0, 25, 26, 1, -1, 6, 31, 25, 0, 31, 1, 2, 8, 7, -1, 27, 26, 32, 33, -1, 7, 32, 26, 1, -1, 8, 33, 32, 7, -1, 2, 27, 33, 8, -1, 1, 26, 27, 2, 31, 2, 3, 10, 9, -1, 28, 27, 34, 35, -1, 9, 34, 27, 2, -1, 10, 35, 34, 9, -1, 3, 28, 35, 10, -1, 2, 27, 28, 3, 31, 3, 4, 12, 11, -1, 29, 28, 36, 37, -1, 11, 36, 28, 3, -1, 12, 37, 36, 11, -1, 4, 29, 37, 12, -1, 3, 28, 29, 4, 31, 4, 5, 14, 13, -1, 30, 29, 38, 39, -1, 13, 38, 29, 4, -1, 14, 39, 38, 13, -1, 5, 30, 39, 14, -1, 4, 29, 30, 5, 31, 5, 6, 16, 15, -1, 31, 30, 40, 41, -1, 15, 40, 30, 5, -1, 16, 41, 40, 15, -1, 6, 31, 41, 16, -1, 5, 30, 31, 6, 31, 6, 1, 18, 17, -1, 26, 31, 42, 43, -1, 17, 42, 31, 6, -1, 18, 43, 42, 17, -1, 1, 26, 43, 18, -1, 6, 31, 26, 1, 31, 19, 18, 1, 7, -1, 43, 44, 32, 26, -1, 7, 32, 44, 19, -1, 1, 26, 32, 7, -1, 18, 43, 26, 1, -1, 19, 44, 43, 18, 31, 20, 8, 2, 9, -1, 33, 45, 34, 27, -1, 9, 34, 45, 20, -1, 2, 27, 34, 9, -1, 8, 33, 27, 2, -1, 20, 45, 33, 8, 31, 21, 10, 3, 11, -1, 35, 46, 36, 28, -1, 11, 36, 46, 21, -1, 3, 28, 36, 11, -1, 10, 35, 28, 3, -1, 21, 46, 35, 10, 31, 22, 12, 4, 13, -1, 37, 47, 38, 29, -1, 13, 38, 47, 22, -1, 4, 29, 38, 13, -1, 12, 37, 29, 4, -1, 22, 47, 37, 12, 31, 23, 14, 5, 15, -1, 39, 48, 40, 30, -1, 15, 40, 48, 23, -1, 5, 30, 40, 15, -1, 14, 39, 30, 5, -1, 23, 48, 39, 14, 31, 24, 16, 6, 17, -1, 41, 49, 42, 31, -1, 17, 42, 49, 24, -1, 6, 31, 42, 17, -1, 16, 41, 31, 6, -1, 24, 49, 41, 16, 31, 50, 51, 52, 53, 54, 55, -1, 56, 61, 60, 59, 58, 57, -1, 50, 56, 57, 51, -1, 51, 57, 58, 52, -1, 52, 58, 59, 53, -1, 53, 59, 60, 54, -1, 54, 60, 61, 55, -1, 55, 61, 56, 50])
        cI = DataArrayInt([0, 23, 46, 69, 92, 115, 138, 168, 198, 228, 258, 288, 318, 348, 378, 408, 438, 468, 498, 542])
        mesh.setConnectivity(c, cI)
        mesh.mergeNodes(eps) # the initial case has double nodes

        ret = mesh.conformize3D(eps)

        mretDesc, _, _, _, _ = mesh.buildDescendingConnectivity()
        mretDesc2, _, _, _, _ = mretDesc.buildDescendingConnectivity()
        c0, cI0 = mesh.getNodalConnectivity().getValues(), mesh.getNodalConnectivityIndex().getValues()
        c, cI = mretDesc.getNodalConnectivity().getValues(), mretDesc.getNodalConnectivityIndex().getValues()
        c2, cI2 = mretDesc2.getNodalConnectivity().getValues(), mretDesc2.getNodalConnectivityIndex().getValues()
        cRef0 = [31, 1, 0, 2, -1, 25, 26, 27, -1, 2, 27, 26, 1, -1, 0, 25, 27, 2, -1, 1, 26, 25, 0, 31, 2, 0, 3, -1, 25, 27, 28, -1, 3, 28, 27, 2, -1, 0, 25, 28, 3, -1, 2, 27, 25, 0, 31, 3, 0, 4, -1, 25, 28, 29, -1, 4, 29, 28, 3, -1, 0, 25, 29, 4, -1, 3, 28, 25, 0, 31, 4, 0, 5, -1, 25, 29, 30, -1, 5, 30, 29, 4, -1, 0, 25, 30, 5, -1, 4, 29, 25, 0, 31, 5, 0, 6, -1, 25, 30, 31, -1, 6, 31, 30, 5, -1, 0, 25, 31, 6, -1, 5, 30, 25, 0, 31, 6, 0, 1, -1, 25, 31, 26, -1, 1, 26, 31, 6, -1, 0, 25, 26, 1, -1, 6, 31, 25, 0, 31, 1, 2, 8, 7, -1, 27, 26, 32, 33, -1, 7, 32, 26, 1, -1, 8, 33, 32, 7, -1, 2, 27, 33, 8, -1, 1, 26, 27, 2, 31, 2, 3, 10, 9, -1, 28, 27, 34, 35, -1, 9, 34, 27, 2, -1, 10, 35, 34, 9, -1, 3, 28, 35, 10, -1, 2, 27, 28, 3, 31, 3, 4, 12, 11, -1, 29, 28, 36, 37, -1, 11, 36, 28, 3, -1, 12, 37, 36, 11, -1, 4, 29, 37, 12, -1, 3, 28, 29, 4, 31, 4, 5, 14, 13, -1, 30, 29, 38, 39, -1, 13, 38, 29, 4, -1, 14, 39, 38, 13, -1, 5, 30, 39, 14, -1, 4, 29, 30, 5, 31, 5, 6, 16, 15, -1, 31, 30, 40, 41, -1, 15, 40, 30, 5, -1, 16, 41, 40, 15, -1, 6, 31, 41, 16, -1, 5, 30, 31, 6, 31, 6, 1, 18, 17, -1, 26, 31, 42, 43, -1, 17, 42, 31, 6, -1, 18, 43, 42, 17, -1, 1, 26, 43, 18, -1, 6, 31, 26, 1, 31, 19, 18, 1, 7, -1, 43, 44, 32, 26, -1, 7, 32, 44, 19, -1, 1, 26, 32, 7, -1, 18, 43, 26, 1, -1, 19, 44, 43, 18, 31, 20, 8, 2, 9, -1, 33, 45, 34, 27, -1, 9, 34, 45, 20, -1, 2, 27, 34, 9, -1, 8, 33, 27, 2, -1, 20, 45, 33, 8, 31, 21, 10, 3, 11, -1, 35, 46, 36, 28, -1, 11, 36, 46, 21, -1, 3, 28, 36, 11, -1, 10, 35, 28, 3, -1, 21, 46, 35, 10, 31, 22, 12, 4, 13, -1, 37, 47, 38, 29, -1, 13, 38, 47, 22, -1, 4, 29, 38, 13, -1, 12, 37, 29, 4, -1, 22, 47, 37, 12, 31, 23, 14, 5, 15, -1, 39, 48, 40, 30, -1, 15, 40, 48, 23, -1, 5, 30, 40, 15, -1, 14, 39, 30, 5, -1, 23, 48, 39, 14, 31, 24, 16, 6, 17, -1, 41, 49, 42, 31, -1, 17, 42, 49, 24, -1, 6, 31, 42, 17, -1, 16, 41, 31, 6, -1, 24, 49, 41, 16, 31, 50, 55, 54, 53, 52, 51, -1, 46, 50, 51, 47, 37, 36, -1, 47, 51, 52, 48, 39, 38, -1, 48, 52, 53, 49, 41, 40, -1, 49, 53, 54, 44, 43, 42, -1, 44, 54, 55, 45, 33, 32, -1, 45, 55, 50, 46, 35, 34, -1, 25, 26, 27, -1, 25, 31, 26, -1, 27, 26, 32, 33, -1, 26, 31, 42, 43, -1, 43, 44, 32, 26, -1, 41, 49, 42, 31, -1, 25, 29, 30, -1, 25, 30, 31, -1, 30, 29, 38, 39, -1, 31, 30, 40, 41, -1, 39, 48, 40, 30, -1, 25, 27, 28, -1, 28, 27, 34, 35, -1, 33, 45, 34, 27, -1, 35, 46, 36, 28, -1, 25, 28, 29, -1, 29, 28, 36, 37, -1, 37, 47, 38, 29]
        cIRef0 = [0, 23, 46, 69, 92, 115, 138, 168, 198, 228, 258, 288, 318, 348, 378, 408, 438, 468, 498, 631]
        cRef = [5, 1, 0, 2, 5, 25, 26, 27, 5, 2, 27, 26, 1, 5, 0, 25, 27, 2, 5, 1, 26, 25, 0, 5, 2, 0, 3, 5, 25, 27, 28, 5, 3, 28, 27, 2, 5, 0, 25, 28, 3, 5, 3, 0, 4, 5, 25, 28, 29, 5, 4, 29, 28, 3, 5, 0, 25, 29, 4, 5, 4, 0, 5, 5, 25, 29, 30, 5, 5, 30, 29, 4, 5, 0, 25, 30, 5, 5, 5, 0, 6, 5, 25, 30, 31, 5, 6, 31, 30, 5, 5, 0, 25, 31, 6, 5, 6, 0, 1, 5, 25, 31, 26, 5, 1, 26, 31, 6, 5, 1, 2, 8, 7, 5, 27, 26, 32, 33, 5, 7, 32, 26, 1, 5, 8, 33, 32, 7, 5, 2, 27, 33, 8, 5, 2, 3, 10, 9, 5, 28, 27, 34, 35, 5, 9, 34, 27, 2, 5, 10, 35, 34, 9, 5, 3, 28, 35, 10, 5, 3, 4, 12, 11, 5, 29, 28, 36, 37, 5, 11, 36, 28, 3, 5, 12, 37, 36, 11, 5, 4, 29, 37, 12, 5, 4, 5, 14, 13, 5, 30, 29, 38, 39, 5, 13, 38, 29, 4, 5, 14, 39, 38, 13, 5, 5, 30, 39, 14, 5, 5, 6, 16, 15, 5, 31, 30, 40, 41, 5, 15, 40, 30, 5, 5, 16, 41, 40, 15, 5, 6, 31, 41, 16, 5, 6, 1, 18, 17, 5, 26, 31, 42, 43, 5, 17, 42, 31, 6, 5, 18, 43, 42, 17, 5, 1, 26, 43, 18, 5, 19, 18, 1, 7, 5, 43, 44, 32, 26, 5, 7, 32, 44, 19, 5, 19, 44, 43, 18, 5, 20, 8, 2, 9, 5, 33, 45, 34, 27, 5, 9, 34, 45, 20, 5, 20, 45, 33, 8, 5, 21, 10, 3, 11, 5, 35, 46, 36, 28, 5, 11, 36, 46, 21, 5, 21, 46, 35, 10, 5, 22, 12, 4, 13, 5, 37, 47, 38, 29, 5, 13, 38, 47, 22, 5, 22, 47, 37, 12, 5, 23, 14, 5, 15, 5, 39, 48, 40, 30, 5, 15, 40, 48, 23, 5, 23, 48, 39, 14, 5, 24, 16, 6, 17, 5, 41, 49, 42, 31, 5, 17, 42, 49, 24, 5, 24, 49, 41, 16, 5, 50, 55, 54, 53, 52, 51, 5, 46, 50, 51, 47, 37, 36, 5, 47, 51, 52, 48, 39, 38, 5, 48, 52, 53, 49, 41, 40, 5, 49, 53, 54, 44, 43, 42, 5, 44, 54, 55, 45, 33, 32, 5, 45, 55, 50, 46, 35, 34]
        cIRef = [0, 4, 8, 13, 18, 23, 27, 31, 36, 41, 45, 49, 54, 59, 63, 67, 72, 77, 81, 85, 90, 95, 99, 103, 108, 113, 118, 123, 128, 133, 138, 143, 148, 153, 158, 163, 168, 173, 178, 183, 188, 193, 198, 203, 208, 213, 218, 223, 228, 233, 238, 243, 248, 253, 258, 263, 268, 273, 278, 283, 288, 293, 298, 303, 308, 313, 318, 323, 328, 333, 338, 343, 348, 353, 358, 363, 368, 373, 378, 385, 392, 399, 406, 413, 420, 427]
        cRef2 = [1, 1, 0, 1, 0, 2, 1, 2, 1, 1, 25, 26, 1, 26, 27, 1, 27, 25, 1, 2, 27, 1, 26, 1, 1, 0, 25, 1, 0, 3, 1, 3, 2, 1, 27, 28, 1, 28, 25, 1, 3, 28, 1, 0, 4, 1, 4, 3, 1, 28, 29, 1, 29, 25, 1, 4, 29, 1, 0, 5, 1, 5, 4, 1, 29, 30, 1, 30, 25, 1, 5, 30, 1, 0, 6, 1, 6, 5, 1, 30, 31, 1, 31, 25, 1, 6, 31, 1, 1, 6, 1, 31, 26, 1, 2, 8, 1, 8, 7, 1, 7, 1, 1, 26, 32, 1, 32, 33, 1, 33, 27, 1, 7, 32, 1, 8, 33, 1, 3, 10, 1, 10, 9, 1, 9, 2, 1, 27, 34, 1, 34, 35, 1, 35, 28, 1, 9, 34, 1, 10, 35, 1, 4, 12, 1, 12, 11, 1, 11, 3, 1, 28, 36, 1, 36, 37, 1, 37, 29, 1, 11, 36, 1, 12, 37, 1, 5, 14, 1, 14, 13, 1, 13, 4, 1, 29, 38, 1, 38, 39, 1, 39, 30, 1, 13, 38, 1, 14, 39, 1, 6, 16, 1, 16, 15, 1, 15, 5, 1, 30, 40, 1, 40, 41, 1, 41, 31, 1, 15, 40, 1, 16, 41, 1, 1, 18, 1, 18, 17, 1, 17, 6, 1, 31, 42, 1, 42, 43, 1, 43, 26, 1, 17, 42, 1, 18, 43, 1, 19, 18, 1, 7, 19, 1, 43, 44, 1, 44, 32, 1, 44, 19, 1, 20, 8, 1, 9, 20, 1, 33, 45, 1, 45, 34, 1, 45, 20, 1, 21, 10, 1, 11, 21, 1, 35, 46, 1, 46, 36, 1, 46, 21, 1, 22, 12, 1, 13, 22, 1, 37, 47, 1, 47, 38, 1, 47, 22, 1, 23, 14, 1, 15, 23, 1, 39, 48, 1, 48, 40, 1, 48, 23, 1, 24, 16, 1, 17, 24, 1, 41, 49, 1, 49, 42, 1, 49, 24, 1, 50, 55, 1, 55, 54, 1, 54, 53, 1, 53, 52, 1, 52, 51, 1, 51, 50, 1, 46, 50, 1, 51, 47, 1, 52, 48, 1, 53, 49, 1, 54, 44, 1, 55, 45]
        cIRef2 = [0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72, 75, 78, 81, 84, 87, 90, 93, 96, 99, 102, 105, 108, 111, 114, 117, 120, 123, 126, 129, 132, 135, 138, 141, 144, 147, 150, 153, 156, 159, 162, 165, 168, 171, 174, 177, 180, 183, 186, 189, 192, 195, 198, 201, 204, 207, 210, 213, 216, 219, 222, 225, 228, 231, 234, 237, 240, 243, 246, 249, 252, 255, 258, 261, 264, 267, 270, 273, 276, 279, 282, 285, 288, 291, 294, 297, 300, 303, 306, 309, 312, 315, 318, 321, 324, 327, 330, 333, 336, 339, 342, 345, 348, 351, 354, 357, 360, 363]
        self.assertEqual(85, mretDesc.getNumberOfCells())
        self.assertEqual(121, mretDesc2.getNumberOfCells())
        self.assertEqual(cRef0, c0)
        self.assertEqual(cIRef0, cI0)
        self.assertEqual(cRef, c)
        self.assertEqual(cIRef, cI)
        self.assertEqual(cRef2, c2)
        self.assertEqual(cIRef2, cI2)
        self.assertEqual(set([18]), set(ret.getValues()))
        pass

if __name__ == '__main__':
    unittest.main()
