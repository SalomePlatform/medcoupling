#  -*- coding: utf-8 -*-
# Copyright (C) 2017  CEA/DEN, EDF R&D
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
    
    pass

if __name__ == '__main__':
    unittest.main()
